"""
Chemical Reservoir Computing (ChemComp Emulation Module)
========================================================

A framework for solving ODEs and PDEs using chemical reaction
network-based reservoir computing with Bayesian optimization.

Core Components
---------------
- :class:`ConstrainedChemicalReservoir` -- Main reservoir computing engine
- :class:`ODEProblem` / :class:`PDEProblem1D` -- Problem specifications
- :class:`BayesianOptimizer` -- Hyperparameter optimization
- :class:`SensitivityAnalyzer` -- Robustness analysis

Configuration
-------------
All tuneable parameters are stored in ``default_config.yaml`` and managed via:

- :func:`load_config` -- Load a YAML config file
- :func:`get_config` -- Access the active :class:`ChemCompConfig`
- :func:`update_config` -- Modify individual parameters at runtime

Quick Start
-----------
>>> import chemical_reservoir as cr
>>> cr.load_config("my_config.yaml")  # optional: customise parameters
>>> problem = cr.problems.create_vanderpol_problem()
>>> reservoir, config = cr.io.load_reservoir_from_json("network_configs/cascade.json")
>>> results = reservoir.solve_ode(problem)  # uses config defaults for n_train, n_test, etc.
"""

from .config import (
    OutputConfig, set_output_directory, get_output_config,
    ChemCompConfig, get_config, load_config, update_config,
    TrainingConfig, ReservoirConfig, EncodingConfig, ClosedLoopConfig,
    OptimizationConfig, MetricsConfig, SensitivityConfig,
)
from .constraints import (
    SpeciesBounds,
    PhysicalConstraints,
    MeasurementConstraints,
)
from .problems import ODEProblem, PDEProblem1D
from .optimization import (
    AcquisitionFunction,
    ExpectedImprovement,
    ProbabilityOfImprovement,
    UpperConfidenceBound,
    BayesianOptimizer,
)
from .metrics import LearningMetrics, MetricsCalculator
from .sensitivity import SensitivityAnalyzer
from .network_viz import ReactionNetworkVisualizer
from .visualization import ResultsVisualizer, save_figure
from .reservoir import ConstrainedChemicalReservoir

from . import problems  # noqa: F811 — re-export submodule
from . import io
from . import networks
from . import custom_problems
from . import validate

__version__ = "1.0.0"
__author__ = "Chemical Computing Team"

__all__ = [
    # Config
    "OutputConfig",
    "set_output_directory",
    "get_output_config",
    # Constraints
    "SpeciesBounds",
    "PhysicalConstraints",
    "MeasurementConstraints",
    # Problems
    "ODEProblem",
    "PDEProblem1D",
    # Optimization
    "AcquisitionFunction",
    "ExpectedImprovement",
    "ProbabilityOfImprovement",
    "UpperConfidenceBound",
    "BayesianOptimizer",
    # Metrics
    "LearningMetrics",
    "MetricsCalculator",
    # Sensitivity
    "SensitivityAnalyzer",
    # Visualization
    "ReactionNetworkVisualizer",
    "ResultsVisualizer",
    "save_figure",
    # Reservoir
    "ConstrainedChemicalReservoir",
]