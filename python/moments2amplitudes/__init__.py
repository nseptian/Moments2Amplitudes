"""
Moments2Amplitudes - Mass-dependent partial wave analysis framework using angular moments.
"""

__version__ = "0.2.0"

from .schema import FitConfig, DatasetConfig, Chi2Config, WaveConfig, MinimizerConfig, MultiStartConfig, MCMCPropagationConfig
from .config import load_config, save_config, validate_config
from .rundir import RunDirectory, get_run_directory, resolve_init_file_for_tbin
from .environment import init_environment
from .compiler import is_aclic_updated, clean_aclic_artifacts, compile_aclic_libraries, ensure_aclic_updated

# Ensure process environment is automatically configured for ROOT 6.36 and BruFit in Python 3.12+
import sys
if sys.version_info >= (3, 10):
    init_environment()

__all__ = [
    "FitConfig",
    "DatasetConfig",
    "Chi2Config",
    "WaveConfig",
    "MinimizerConfig",
    "MultiStartConfig",
    "MCMCPropagationConfig",
    "load_config",
    "save_config",
    "validate_config",
    "RunDirectory",
    "get_run_directory",
    "resolve_init_file_for_tbin",
    "is_aclic_updated",
    "clean_aclic_artifacts",
    "compile_aclic_libraries",
    "ensure_aclic_updated",
]

