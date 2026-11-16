"""
ACLiC compilation management and parallel process synchronization for Moments2Amplitudes.

Guarantees that shared libraries are compiled sequentially in a single process
before any parallel fit workers or t-bin subprocesses start, avoiding ACLiC
concurrency collisions and symbol corruption.
"""

from __future__ import annotations
import os
import sys
import time
import fcntl
import glob
import logging
import threading
import subprocess
from pathlib import Path
from typing import List, Tuple, Optional, Any

from .environment import (
    get_repo_root,
    get_root_binary,
    get_root_subprocess_env,
)

logger = logging.getLogger(__name__)

# Primary C++ translation units that generate ACLiC .so libraries
TARGET_LIBRARIES = [
    "ParameterHelper",
    "Equation",
    "EquationSolver",
    "MomentHelper",
    "MassDependentFitter",
]

_thread_lock = threading.Lock()


class ProcessLock:
    """
    Cross-process file lock using fcntl.flock to synchronize parallel processes.
    Blocks with periodic status updates until the lock is acquired or timeout occurs.
    """

    def __init__(self, lock_file: Path, timeout: float = 300.0, wait_interval: float = 0.5):
        self.lock_file = lock_file
        self.timeout = timeout
        self.wait_interval = wait_interval
        self._fd: Optional[int] = None

    def __enter__(self):
        self.lock_file.parent.mkdir(parents=True, exist_ok=True)
        self._fd = os.open(str(self.lock_file), os.O_RDWR | os.O_CREAT, 0o666)

        start_time = time.time()
        last_log_time = start_time
        logged_waiting = False

        while True:
            try:
                fcntl.flock(self._fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                if logged_waiting:
                    print(f"[ACLiC] Acquired compilation lock ({self.lock_file.name}).")
                return self
            except (BlockingIOError, OSError):
                now = time.time()
                if now - start_time > self.timeout:
                    raise TimeoutError(
                        f"Timed out after {self.timeout:.1f}s waiting for ACLiC compilation lock on {self.lock_file}. "
                        f"Another process may have crashed. If no other fit is running, remove {self.lock_file}."
                    )
                if now - last_log_time >= 5.0:
                    print(f"[ACLiC] Another process is compiling libraries. Waiting... ({now - start_time:.0f}s elapsed)")
                    last_log_time = now
                    logged_waiting = True
                time.sleep(self.wait_interval)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fd is not None:
            try:
                fcntl.flock(self._fd, fcntl.LOCK_UN)
            except Exception:
                pass
            os.close(self._fd)
            self._fd = None


def get_aclic_target_so_paths(repo_root: Optional[Path] = None) -> List[Path]:
    """Return the list of expected .so paths for all target libraries."""
    root = repo_root or get_repo_root()
    src_dir = root / "src"
    return [src_dir / f"{lib}_cpp.so" for lib in TARGET_LIBRARIES]


def get_aclic_source_paths(repo_root: Optional[Path] = None) -> List[Path]:
    """
    Return all source, header, and macro files that affect ACLiC compiled libraries.
    """
    root = repo_root or get_repo_root()
    sources: List[Path] = []

    # C++ sources and headers in src/ and include/
    src_dir = root / "src"
    include_dir = root / "include"
    macros_dir = root / "macros"

    for pattern in ("*.cpp", "*.cxx", "*.cc", "*.h", "*.hpp"):
        sources.extend(src_dir.glob(pattern))
        sources.extend(include_dir.glob(pattern))

    for m in ("Load.C", "Compile.C", "ConfigureAmpsNoValues.C"):
        p = macros_dir / m
        if p.is_file():
            sources.append(p)

    return [
        s for s in sources
        if s.is_file() and "ACLiC" not in s.name and "_tmp_" not in s.name
    ]


def is_aclic_updated(repo_root: Optional[Path] = None) -> Tuple[bool, str]:
    """
    Check if all target ACLiC shared libraries exist, are non-empty, and are newer
    than all source, header, and configuration files.

    Returns:
        (is_updated, reason)
    """
    root = repo_root or get_repo_root()
    so_paths = get_aclic_target_so_paths(root)

    # Check existence and non-zero size
    for p in so_paths:
        if not p.is_file():
            return False, f"Missing target library: {p.name}"
        if p.stat().st_size == 0:
            return False, f"Zero-byte library artifact: {p.name}"

    # Get modification time of the oldest target library
    oldest_so_mtime = min(p.stat().st_mtime for p in so_paths)

    # Compare against modification times of all source/header files
    sources = get_aclic_source_paths(root)
    for src in sources:
        if src.stat().st_mtime > oldest_so_mtime:
            return False, f"Source file '{src.name}' modified after '{p.name}' was built"

    return True, "All ACLiC libraries are up-to-date"


def clean_aclic_artifacts(repo_root: Optional[Path] = None) -> List[str]:
    """
    Remove all stale ACLiC build artifacts (*.so, *.pcm, *.rootmap, *.d, *_tmp_*).
    """
    root = repo_root or get_repo_root()
    src_dir = root / "src"
    removed: List[str] = []

    patterns = [
        "*ACLiC*",
        "*_tmp_*",
        "*.so",
        "*.pcm",
        "*.rootmap",
        "*.d",
    ]

    for pat in patterns:
        for f in src_dir.glob(pat):
            try:
                if f.is_file() or f.is_symlink():
                    f.unlink()
                    removed.append(f.name)
            except OSError as e:
                logger.warning(f"Failed to remove {f}: {e}")

    return removed


def compile_aclic_libraries(
    config: Optional[Any] = None,
    repo_root: Optional[Path] = None,
    timeout: float = 300.0,
) -> bool:
    """
    Execute a single sequential ACLiC compilation of all target libraries using ROOT.
    """
    root = repo_root or get_repo_root()
    compile_macro = root / "macros" / "Compile.C"
    if not compile_macro.is_file():
        raise FileNotFoundError(f"Compile macro not found: {compile_macro}")

    root_bin = get_root_binary(config)
    env = get_root_subprocess_env(config)

    cmd = [
        root_bin,
        "-l",
        "-b",
        "-q",
        str(compile_macro),
    ]

    print(f"[ACLiC] Compiling C++ libraries using {Path(root_bin).name}...")
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        cwd=str(root),
        text=True,
        env=env,
        timeout=timeout,
    )

    if proc.returncode != 0:
        tail = "\n".join(proc.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"ACLiC compilation failed with exit code {proc.returncode}:\n{tail}"
        )

    # Verify that all target .so files now exist and have size > 0
    so_paths = get_aclic_target_so_paths(root)
    for p in so_paths:
        if not p.is_file() or p.stat().st_size == 0:
            raise RuntimeError(
                f"ACLiC compilation finished with exit 0 but expected output '{p.name}' is missing or empty."
            )

    return True


def ensure_aclic_updated(
    config: Optional[Any] = None,
    repo_root: Optional[Path] = None,
    force: bool = False,
    timeout: float = 300.0,
) -> bool:
    """
    Synchronized check and recompile mechanism for parallel execution.

    Uses both thread-level and inter-process file locks so that when multiple
    parallel workers or t-bin jobs start simultaneously:
      1. Only the first process/thread cleans and compiles the libraries.
      2. All other processes/threads wait on the lock until compilation finishes.
      3. Waiting processes then find the libraries up-to-date and proceed immediately.
    """
    root = repo_root or get_repo_root()
    lock_file = root / "src" / ".aclic.lock"

    # Fast-path check without acquiring the lock
    if not force:
        up_to_date, _ = is_aclic_updated(root)
        if up_to_date:
            return True

    with _thread_lock:
        with ProcessLock(lock_file, timeout=timeout):
            # Double-checked locking pattern: re-check after acquiring the lock
            if not force:
                up_to_date, reason = is_aclic_updated(root)
                if up_to_date:
                    return True
            else:
                reason = "Forced recompilation requested"

            print(f"[ACLiC] Libraries require recompilation: {reason}")
            removed = clean_aclic_artifacts(root)
            if removed:
                print(f"[ACLiC] Cleaned {len(removed)} stale build artifacts.")

            start_t = time.time()
            compile_aclic_libraries(config=config, repo_root=root, timeout=timeout)
            duration = time.time() - start_t
            print(f"[ACLiC] Successfully compiled all shared libraries in {duration:.1f}s.")
            return True
