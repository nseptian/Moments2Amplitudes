"""
Centralized environment resolution for Moments2Amplitudes.
Ensures external ROOT 6.36 and BruFit installations take strict precedence
over any internal conda/Pixi ROOT libraries.
"""

from __future__ import annotations
import os
import sys
import shutil
from pathlib import Path
from typing import Dict, Optional, Any

DEFAULT_ROOTSYS = "/d/grid17/septian/root_v6-36-08_cuda_build"
DEFAULT_BRUFIT = "/d/grid17/septian/brufit_dev/brufit_negative"


def get_default_max_workers() -> int:
    """Return default maximum parallel workers/threads (half of available CPU cores)."""
    cpus = os.cpu_count() or 2
    return max(1, cpus // 2)


def get_repo_root() -> Path:
    """Find repository root directory."""
    if "MOMS2AMPS" in os.environ:
        p = Path(os.environ["MOMS2AMPS"]).resolve()
        if p.is_dir():
            return p
    return Path(__file__).resolve().parent.parent.parent


def load_dot_env() -> Dict[str, str]:
    """Parse .env file from repository root if present."""
    repo_root = get_repo_root()
    env_file = repo_root / ".env"
    env_vars: Dict[str, str] = {}
    if env_file.is_file():
        try:
            with open(env_file, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    if "=" in line:
                        key, val = line.split("=", 1)
                        env_vars[key.strip()] = val.strip().strip("'\"")
        except Exception:
            pass
    return env_vars


def is_pixi_or_conda_path(path_str: Optional[str]) -> bool:
    """Return True if path points to a pixi/conda environment rather than external system ROOT."""
    if not path_str:
        return False
    p = str(path_str)
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix and p.startswith(conda_prefix):
        return True
    if ".pixi" in p or "miniconda" in p or "anaconda" in p or "envs/default" in p:
        return True
    return False


def resolve_rootsys(config: Optional[Any] = None) -> str:
    """
    Resolve ROOTSYS with strict preference for system ROOT 6.36,
    avoiding Pixi's bundled ROOT 6.28.
    """
    # 1. Config override
    if config and getattr(config, "environment", None) and config.environment.root_dir:
        cand = config.environment.root_dir
        if Path(cand).is_dir():
            return str(cand)

    # 2. .env file
    dot_env = load_dot_env()
    if "ROOTSYS" in dot_env and Path(dot_env["ROOTSYS"]).is_dir():
        return dot_env["ROOTSYS"]

    # 3. CONDA_BACKUP_ROOTSYS (saved by conda activate before pixi overwrote it)
    if "CONDA_BACKUP_ROOTSYS" in os.environ and os.environ["CONDA_BACKUP_ROOTSYS"]:
        cand = os.environ["CONDA_BACKUP_ROOTSYS"]
        if Path(cand).is_dir() and not is_pixi_or_conda_path(cand):
            return cand

    # 4. ROOTSYS in environment, if NOT pixi/conda
    if "ROOTSYS" in os.environ:
        cand = os.environ["ROOTSYS"]
        if not is_pixi_or_conda_path(cand) and Path(cand).is_dir():
            return cand

    # 5. Default system ROOT 6.36
    if Path(DEFAULT_ROOTSYS).is_dir():
        return DEFAULT_ROOTSYS

    # Fallback to whatever is in ROOTSYS or default
    return os.environ.get("ROOTSYS", DEFAULT_ROOTSYS)


def resolve_brufit(config: Optional[Any] = None) -> str:
    """Resolve BRUFIT directory."""
    if config and getattr(config, "environment", None) and config.environment.brufit_dir:
        cand = config.environment.brufit_dir
        if Path(cand).is_dir():
            return str(cand)

    dot_env = load_dot_env()
    if "BRUFIT" in dot_env and Path(dot_env["BRUFIT"]).is_dir():
        return dot_env["BRUFIT"]

    if "BRUFIT" in os.environ and Path(os.environ["BRUFIT"]).is_dir():
        return os.environ["BRUFIT"]

    return DEFAULT_BRUFIT


def get_root_binary(config: Optional[Any] = None) -> str:
    """Get the path to the ROOT 6.36 binary."""
    rootsys = resolve_rootsys(config)
    cand = Path(rootsys) / "bin" / "root"
    if cand.is_file() and os.access(cand, os.X_OK):
        return str(cand)

    # If default rootsys has binary
    default_cand = Path(DEFAULT_ROOTSYS) / "bin" / "root"
    if default_cand.is_file() and os.access(default_cand, os.X_OK):
        return str(default_cand)

    # Search PATH (skip pixi/conda root if possible)
    which_root = shutil.which("root")
    if which_root and not is_pixi_or_conda_path(which_root):
        return which_root

    # Last resort fallback: check sys.prefix / pixi root or any which_root
    pixi_root = Path(sys.prefix) / "bin" / "root"
    if pixi_root.is_file() and os.access(pixi_root, os.X_OK):
        return str(pixi_root)
    if which_root and os.access(which_root, os.X_OK):
        return which_root

    return str(cand)


def filter_pixi_paths(path_str: str) -> str:
    """Filter out any pixi/conda directories from a colon-separated path string."""
    if not path_str:
        return ""
    parts = path_str.split(":")
    clean_parts = [p for p in parts if p and not is_pixi_or_conda_path(p)]
    return ":".join(clean_parts)


def get_root_subprocess_env(config: Optional[Any] = None) -> Dict[str, str]:
    """
    Build environment variables for running ROOT macros.
    Ensures ROOT 6.36 and BruFit take precedence and completely filters out
    Pixi's conda library paths (LD_LIBRARY_PATH and PATH) to avoid ABI mismatches.
    """
    env = os.environ.copy()
    rootsys = resolve_rootsys(config)
    brufit = resolve_brufit(config)
    repo_root = get_repo_root()

    env["ROOTSYS"] = rootsys
    env["BRUFIT"] = brufit
    env["MOMS2AMPS"] = str(repo_root)

    # Clean any pixi/conda paths from existing PATH and LD_LIBRARY_PATH
    clean_path = filter_pixi_paths(env.get("PATH", ""))
    clean_ld = filter_pixi_paths(env.get("LD_LIBRARY_PATH", ""))

    # Prepend ROOTSYS/bin to clean PATH
    rootsys_bin = str(Path(rootsys) / "bin")
    env["PATH"] = f"{rootsys_bin}:{clean_path}" if clean_path else rootsys_bin

    # Prepend ROOTSYS/lib and BRUFIT library directories to clean LD_LIBRARY_PATH
    brufit_p = Path(brufit)
    brufit_ld_candidates = [
        brufit_p / "install" / "lib64",
        brufit_p / "install" / "lib",
        brufit_p / "build" / "lib",
        brufit_p / "lib",
    ]
    brufit_ld = [str(p) for p in brufit_ld_candidates if p.is_dir()]
    if not brufit_ld:
        brufit_ld = [str(brufit_p / "lib")]

    ld_paths = [str(Path(rootsys) / "lib")] + brufit_ld
    if clean_ld:
        ld_paths.append(clean_ld)
    env["LD_LIBRARY_PATH"] = ":".join(ld_paths)

    # ROOT_INCLUDE_PATH
    brufit_inc_candidates = [
        brufit_p / "install" / "include",
        brufit_p / "include",
        brufit_p / "core",
        brufit_p / "utility",
    ]
    brufit_inc = [str(p) for p in brufit_inc_candidates if p.is_dir()]
    if not brufit_inc:
        brufit_inc = [str(brufit_p / "include")]

    inc_paths = [
        str(repo_root / "include"),
        str(repo_root / "macros"),
    ] + brufit_inc + [
        str(Path(rootsys) / "include"),
    ]
    existing_inc = filter_pixi_paths(env.get("ROOT_INCLUDE_PATH", ""))
    if existing_inc:
        inc_paths.append(existing_inc)
    env["ROOT_INCLUDE_PATH"] = ":".join(inc_paths)

    env["MPLCONFIGDIR"] = "/tmp"

    return env


def init_environment() -> None:
    """
    Configure the current process environment to prioritize ROOT 6.36 and BruFit.
    Called automatically on package import.
    """
    rootsys = resolve_rootsys()
    brufit = resolve_brufit()
    repo_root = get_repo_root()

    os.environ["ROOTSYS"] = rootsys
    os.environ["BRUFIT"] = brufit
    os.environ["MOMS2AMPS"] = str(repo_root)
    os.environ.setdefault("MPLCONFIGDIR", "/tmp")

    rootsys_bin = str(Path(rootsys) / "bin")
    current_path = os.environ.get("PATH", "")
    if not current_path.startswith(rootsys_bin):
        os.environ["PATH"] = f"{rootsys_bin}:{current_path}"

    ld_paths = [str(Path(rootsys) / "lib"), str(Path(brufit) / "lib")]
    existing_ld = os.environ.get("LD_LIBRARY_PATH", "")
    if existing_ld:
        ld_paths.append(existing_ld)
    os.environ["LD_LIBRARY_PATH"] = ":".join(ld_paths)

    inc_paths = [
        str(repo_root / "include"),
        str(repo_root / "macros"),
        str(Path(brufit) / "include"),
        str(Path(rootsys) / "include"),
    ]
    existing_inc = os.environ.get("ROOT_INCLUDE_PATH", "")
    if existing_inc:
        inc_paths.append(existing_inc)
    os.environ["ROOT_INCLUDE_PATH"] = ":".join(inc_paths)

