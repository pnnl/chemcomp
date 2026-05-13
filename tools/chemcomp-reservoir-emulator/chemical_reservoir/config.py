"""
Configuration module for output settings and global parameters.

All tuneable parameters are stored in a YAML config file (default_config.yaml).
Users can override defaults by calling ``load_config("my_config.yaml")`` or by
passing a config dict to individual functions.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Any
import warnings

import yaml

# Suppress warnings globally
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Default config path (bundled with the package)
# ---------------------------------------------------------------------------
_PACKAGE_DIR = Path(__file__).parent
_DEFAULT_CONFIG_PATH = _PACKAGE_DIR.parent / "default_config.yaml"


# ---------------------------------------------------------------------------
# Dataclasses for typed access
# ---------------------------------------------------------------------------

@dataclass
class TrainingConfig:
    n_train: int = 500
    n_test: int = 200
    test_horizon_factor: float = 0.5
    n_train_time: int = 100
    n_test_time: int = 50
    pde_test_horizon_factor: float = 0.3


@dataclass
class ReservoirConfig:
    ridge_alpha: float = 1e-4
    feature_expansion: str = "quadratic"
    delays: List[int] = field(default_factory=lambda: [1, 2, 4, 8])
    ridge_cv_alphas_log_min: float = -6
    ridge_cv_alphas_log_max: float = 2
    ridge_cv_alphas_n: int = 50


@dataclass
class EncodingConfig:
    method: str = "mixed"
    sinusoidal_amplitude: float = 0.3
    mixed_linear_coeff: float = 0.2
    mixed_sine_coeff: float = 0.2
    mixed_gaussian_coeff: float = 0.1
    mixed_gaussian_width: float = 0.1
    multiharmonic_amplitude: float = 0.15
    multiharmonic_n_harmonics: int = 8


@dataclass
class ClosedLoopConfig:
    feedback_scale: float = 0.5
    normalization_eps: float = 1e-8


@dataclass
class OptimizationConfig:
    n_iterations: int = 50
    n_initial: int = 10
    random_state: int = 42
    ei_xi: float = 0.01
    pi_xi: float = 0.01
    ucb_kappa: float = 2.0
    gp_length_scale_bounds: List[float] = field(default_factory=lambda: [0.01, 100.0])
    gp_noise_level: float = 1e-3
    gp_noise_bounds: List[float] = field(default_factory=lambda: [1e-4, 10.0])
    gp_n_restarts: int = 5
    acq_n_restarts: int = 10
    acq_n_samples: int = 1000
    sigma_floor: float = 1e-9
    failure_penalty: float = 1e10


@dataclass
class MetricsConfig:
    eps: float = 1e-10


@dataclass
class SensitivityConfig:
    perturbation_fractions: List[float] = field(
        default_factory=lambda: [0.01, 0.05, 0.1, 0.2, 0.5]
    )
    n_samples: int = 5
    local_epsilon: float = 0.01
    base_val_floor: float = 0.1
    max_measurement_combos: int = 30
    stability_eps: float = 1e-10


@dataclass
class OutputConfig:
    """Configuration for output files and directories."""
    output_dir: str = "test_results"
    figure_format: str = "pdf"
    figure_dpi: int = 300
    save_figures: bool = True
    show_figures: bool = False

    def __post_init__(self):
        """Create output directory if it doesn't exist."""
        if self.save_figures and not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            print(f"Created output directory: {self.output_dir}")

    def get_path(self, filename: str) -> str:
        """Get full path for output file."""
        if not filename.endswith(f".{self.figure_format}"):
            filename = f"{filename}.{self.figure_format}"
        return os.path.join(self.output_dir, filename)

    def get_subdir(self, subdir: str) -> str:
        """Get or create subdirectory."""
        path = os.path.join(self.output_dir, subdir)
        if not os.path.exists(path):
            os.makedirs(path)
        return path


# ---------------------------------------------------------------------------
# Master configuration container
# ---------------------------------------------------------------------------

@dataclass
class ChemCompConfig:
    """Master configuration holding all sub-configs."""
    training: TrainingConfig = field(default_factory=TrainingConfig)
    reservoir: ReservoirConfig = field(default_factory=ReservoirConfig)
    encoding: EncodingConfig = field(default_factory=EncodingConfig)
    closed_loop: ClosedLoopConfig = field(default_factory=ClosedLoopConfig)
    optimization: OptimizationConfig = field(default_factory=OptimizationConfig)
    metrics: MetricsConfig = field(default_factory=MetricsConfig)
    sensitivity: SensitivityConfig = field(default_factory=SensitivityConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


def _build_config_from_dict(raw: Dict[str, Any]) -> ChemCompConfig:
    """Build a ChemCompConfig from a raw dictionary (e.g. parsed YAML)."""
    kwargs = {}
    section_map = {
        "training": TrainingConfig,
        "reservoir": ReservoirConfig,
        "encoding": EncodingConfig,
        "closed_loop": ClosedLoopConfig,
        "optimization": OptimizationConfig,
        "metrics": MetricsConfig,
        "sensitivity": SensitivityConfig,
        "output": OutputConfig,
    }
    for key, cls in section_map.items():
        section_data = raw.get(key, {})
        if section_data:
            kwargs[key] = cls(**section_data)
    return ChemCompConfig(**kwargs)


# ---------------------------------------------------------------------------
# Global config singleton
# ---------------------------------------------------------------------------

_GLOBAL_CONFIG: Optional[ChemCompConfig] = None


def get_config() -> ChemCompConfig:
    """Get the active global configuration (loads defaults if not set)."""
    global _GLOBAL_CONFIG
    if _GLOBAL_CONFIG is None:
        _GLOBAL_CONFIG = load_config()
    return _GLOBAL_CONFIG


def load_config(filepath: Optional[str] = None) -> ChemCompConfig:
    """Load configuration from a YAML file and set it as the global config.

    Args:
        filepath: Path to a YAML config file. If None, loads the bundled
                  default_config.yaml.

    Returns:
        The loaded ChemCompConfig instance (also set as the global config).
    """
    global _GLOBAL_CONFIG

    if filepath is None:
        filepath = str(_DEFAULT_CONFIG_PATH)

    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {filepath}")

    with open(path) as f:
        raw = yaml.safe_load(f) or {}

    _GLOBAL_CONFIG = _build_config_from_dict(raw)
    return _GLOBAL_CONFIG


def update_config(**overrides) -> ChemCompConfig:
    """Update specific config values without reloading from disk.

    Args:
        **overrides: Dot-notation keys mapped to values, e.g.
            update_config(training__n_train=1000, optimization__n_iterations=100)
            (use double-underscore to separate section from field).

    Returns:
        The updated global config.
    """
    cfg = get_config()
    for key, value in overrides.items():
        parts = key.split("__")
        if len(parts) == 2:
            section, field_name = parts
            sub = getattr(cfg, section, None)
            if sub is not None and hasattr(sub, field_name):
                setattr(sub, field_name, value)
            else:
                raise ValueError(f"Unknown config key: {key}")
        else:
            raise ValueError(
                f"Key must be 'section__field' (double underscore), got: {key}"
            )
    return cfg


# ---------------------------------------------------------------------------
# Legacy API (backwards-compatible)
# ---------------------------------------------------------------------------

def get_output_config() -> OutputConfig:
    """Get the global output configuration."""
    return get_config().output


def set_output_directory(path: str, figure_format: str = "pdf",
                         save_figures: bool = True,
                         show_figures: bool = False) -> OutputConfig:
    """Set the global output configuration."""
    cfg = get_config()
    cfg.output = OutputConfig(
        output_dir=path,
        figure_format=figure_format,
        save_figures=save_figures,
        show_figures=show_figures
    )
    return cfg.output