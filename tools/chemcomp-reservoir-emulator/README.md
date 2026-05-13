# ChemComp — Emulation Module

The emulation module of the **ChemComp** software package. Solves ordinary and
partial differential equations using chemical reaction networks as a reservoir
computing substrate.

A chemical reaction network is driven by an external input signal derived from
the target ODE/PDE solution.  The resulting species concentrations evolve
according to mass-action kinetics, producing rich, nonlinear dynamics.  A
ridge-regression readout layer is then trained on the recorded concentration
trajectories to map reservoir state back to the target solution — including
forecasting into unseen time intervals.

---

## How it works

```
                ┌──────────────────┐
  ODE / PDE     │  Input encoding  │   constant, sinusoidal,
  solution  ──▶ │  (signal → rates)│   mixed, or multiharmonic
                └────────┬─────────┘
                         │
                         ▼
                ┌──────────────────┐
                │  Chemical        │   Mass-action kinetics
                │  reaction network│   (the "reservoir")
                └────────┬─────────┘
                         │ species concentrations
                         ▼
                ┌──────────────────┐
                │  Feature         │   none, quadratic
                │  expansion       │   (cross-terms), or delayed
                └────────┬─────────┘   (time-delay embedding)
                         │
                         ▼
                ┌──────────────────┐
                │  Ridge regression│   RidgeCV with α ∈ [10⁻⁶, 10²]
                │  readout         │   (50 log-spaced values)
                └────────┬─────────┘
                         │
                         ▼
               Predicted ODE/PDE solution
```

| Concept | Description |
|---------|-------------|
| **Reservoir** | A chemical reaction network whose mass-action kinetics generate nonlinear, high-dimensional dynamics from a low-dimensional input signal. |
| **Input encoding** | The target solution is encoded into species production rates.  Four modes are available: `constant`, `sinusoidal`, `mixed` (default), and `multiharmonic` (8 harmonics). |
| **Feature expansion** | Raw reservoir concentrations can be augmented before readout: `none` (concentrations only), `quadratic` (all pairwise cross-terms), or `delayed` (time-delay embedding with delays [1, 2, 4, 8] — recommended for generalisation). |
| **Readout** | A cross-validated ridge regression (`sklearn.linear_model.RidgeCV`) maps features to the target solution. |
| **Prediction mode** | **Open-loop**: encoding drives the reservoir over the full time domain; predictions in one pass. **Closed-loop** (default): training uses teacher forcing; testing feeds predicted outputs back to modulate the encoding, enabling autonomous attractor tracking. |
| **Bayesian optimisation** | Gaussian-process–based search over input rates and initial concentrations to tune reservoir performance. |
| **Sensitivity analysis** | Perturbation studies over concentrations, input rates, measurement sets, and local kinetic parameters. |

---

## Installation

Requires **Python ≥ 3.9**.

### With conda (recommended)

```bash
conda env create -f environment.yml   # creates the "chemcomp" env
conda activate chemcomp
# The package is installed in editable mode automatically
```

### With pip

```bash
pip install -e ".[notebooks]"
```

### Verify

```bash
python -c "import chemical_reservoir; print(chemical_reservoir.__version__)"
# 1.0.0
```

### Dependencies

numpy ≥ 1.21, scipy ≥ 1.7, scikit-learn ≥ 1.0, matplotlib ≥ 3.4,
networkx ≥ 2.6, pandas ≥ 1.3.  The `[notebooks]` extra adds JupyterLab and
ipywidgets.

---

## Quick start

### Python API

```python
import chemical_reservoir as cr

# 1. Load a reaction network from JSON
reservoir, config = cr.io.load_reservoir_from_json("network_configs/cascade.json")

# 2. Define a target ODE problem
problem = cr.problems.create_vanderpol_problem(mu=1.5)

# 3. Solve — trains the reservoir on the first n_train points,
#    then predicts on n_test unseen points (closed-loop by default)
results = reservoir.solve_ode_closed_loop(problem, n_train=400, n_test=200)

# 4. Evaluate
metrics = cr.MetricsCalculator.calculate_all(
    y_true=results["y_true"],
    y_pred=results["y_pred"],
    t=results["t_all"],
    n_train=len(results["t_train"]),
    variable_names=problem.variable_names,
)
print(f"Test RMSE: {metrics.test_rmse:.4e}")
print(f"Test R²:   {metrics.test_r2:.4f}")

# 5. Visualise
viz = cr.ResultsVisualizer()
viz.plot_learning_metrics_dashboard(results, metrics, problem)
```

### Command-line driver

The `run_reservoir.py` script provides a one-command interface that runs a
problem, generates plots, and saves all outputs to a timestamped directory:

```bash
# Solve the Van der Pol oscillator with the cascade network
python run_reservoir.py --problem vanderpol --network network_configs/cascade.json

# Use open-loop mode instead of the default closed-loop
python run_reservoir.py --problem vanderpol --network network_configs/cascade.json --mode open-loop

# Solve a PDE with a custom output directory
python run_reservoir.py --problem heat_equation --network network_configs/dense.json -o my_results

# List all available problems and networks
python run_reservoir.py --list-problems
python run_reservoir.py --list-networks
```

**Driver options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--problem` | *(required)* | Problem name (see `--list-problems`) |
| `--network` | *(required)* | Path to a reaction-network JSON file |
| `-o` / `--output` | `results/<problem>_<network>_<timestamp>/` | Output directory |
| `--n-train` | 800 | Number of training time-points |
| `--n-test` | 400 | Number of test time-points |
| `--encoding` | `mixed` | Input encoding: `constant`, `sinusoidal`, `mixed`, `multiharmonic` |
| `--features` | `delayed` | Feature expansion: `none`, `quadratic`, `delayed` |
| `--mode` | `closed-loop` | Prediction mode: `open-loop` or `closed-loop` |
| `--feedback-scale` | `0.5` | Output-to-input feedback strength (closed-loop only) |

**Outputs saved:** `prediction_overlay.pdf/.png`, `error_over_time.pdf`,
`phase_portrait.pdf` (2-D+ ODEs), `pde_comparison.pdf/.png` (PDEs),
`reservoir_trajectories.pdf`, `reaction_network*.pdf`, `metrics.json`,
`reservoir_params.json`, `network_config.json`.

---

## Project structure

```
chemical_reservoir/
├── pyproject.toml            # Package metadata & build config (BSD-3-Clause)
├── environment.yml           # Conda environment specification (chemcomp)
├── README.md
├── run_reservoir.py          # CLI driver script (see above)
├── chemical_reservoir/       # Installable Python package
│   ├── __init__.py           # Public API & version (1.0.0)
│   ├── config.py             # OutputConfig dataclass
│   ├── constraints.py        # Physical & measurement constraints
│   ├── io.py                 # JSON load / save / validation for networks
│   ├── metrics.py            # RMSE, R², MAE, correlation metrics
│   ├── networks.py           # Built-in reaction network definitions
│   ├── network_viz.py        # Bipartite reaction-network graph visualisation
│   ├── optimization.py       # Bayesian optimiser (GP + acquisition functions)
│   ├── problems.py           # ODE & PDE problem dataclasses and factories
│   ├── reservoir.py          # ConstrainedChemicalReservoir (core engine)
│   ├── run_analysis.py       # High-level analysis pipelines
│   ├── sensitivity.py        # Sensitivity & perturbation analysis
│   └── visualization.py      # Multi-panel result dashboards
├── network_configs/          # JSON reaction-network definitions
│   ├── cascade.json
│   ├── high_dimensional_random.json
│   ├── minimal.json
│   └── ...  (9 networks total)
├── notebooks/                # Tutorial & demonstration notebooks
│   ├── 01_mathematical_foundations.ipynb
│   ├── 02_quickstart_tutorial.ipynb
│   ├── 03_advanced_examples.ipynb
│   ├── 04_api_reference.ipynb
│   └── 05_open_vs_closed_loop.ipynb
├── scripts/                  # Utility scripts
│   └── gen_highdim_network.py
└── tests/                    # Pytest test suite (9 tests)
    └── test_basic.py
```

---

## Available problems

### ODE problems

| Problem | Factory function | Variables | Description |
|---------|-----------------|-----------|-------------|
| Van der Pol | `create_vanderpol_problem()` | x, v | Relaxation oscillator with limit cycle |
| Schnakenberg | `create_schnakenberg_problem()` | x, y | Morphogenesis / Turing pattern oscillator |
| Brusselator | `create_brusselator_problem()` | x, y | Chemical oscillator model |
| Lotka–Volterra | `create_lotka_volterra_problem()` | prey, predator | Predator–prey population cycles |
| Lorenz | `create_lorenz_problem()` | x, y, z | Chaotic attractor (3-D) |

### PDE problems (1-D, method of lines)

PDEs are semi-discretised in space onto a 1-D grid and solved as a coupled
ODE system.

| Problem | Factory function | Type |
|---------|-----------------|------|
| Heat equation | `create_heat_equation_problem()` | Parabolic diffusion |
| Advection–diffusion | `create_advection_diffusion_problem()` | Parabolic transport |
| Fisher equation | `create_fisher_equation_problem()` | Reaction–diffusion |
| Burgers equation | `create_burgers_equation_problem()` | Nonlinear hyperbolic |

---

## Reaction networks

### Built-in networks (defined in code)

These are available via `chemical_reservoir.networks.list_builtin_networks()` and
can be instantiated without a JSON file.

| Network | Species | Reactions | Topology |
|---------|---------|-----------|----------|
| `minimal` | 2 | 3 | Baseline two-species system |
| `simple_autocatalytic` | 3 | 6 | Autocatalytic feedback cycle |
| `linear_chain` | 5 | 8 | Sequential A → B → C → D → E |
| `cascade` | 5 | 9 | Core reactions + cascade pathway |
| `feedback` | 5 | 8 | Negative feedback loop |
| `oscillatory` | 6 | 6 | Brusselator-like oscillator |
| `dense` | 5 | 12 | Highly connected topology |
| `complex_pathways` | 5 | 10 | Parallel reaction pathways |

### High-dimensional network (JSON only)

| Network | Species | Reactions | Description |
|---------|---------|-----------|-------------|
| `high_dimensional_random` | 20 | 50 | 5 coupled Brusselator oscillator modules with weak inter-module coupling.  Produces limit-cycle dynamics well-suited for challenging problems like the Van der Pol oscillator. |

All networks are stored as JSON in `network_configs/` and loaded via
`cr.io.load_reservoir_from_json()`.  You can also define your own network by
writing a JSON file with the same schema (see the notebooks for examples).

---

## Notebooks

The `notebooks/` directory contains four Jupyter notebooks that progress from
theory to advanced usage:

1. **Mathematical Foundations** — Reservoir computing theory, mass-action
   kinetics, ridge regression, and how they fit together.
2. **Quick-Start Tutorial** — Load a network, solve an ODE, evaluate metrics,
   and visualise predictions step by step.
3. **Advanced Examples** — PDE solving, Bayesian optimisation, sensitivity
   analysis, custom network design, and horizon comparison.
4. **API Reference** — Class-by-class walkthrough with code examples for every
   public component.

---

## Tests

```bash
python -m pytest tests/ -v
```

All 9 tests cover imports, built-in network listing, problem factories, JSON
I/O, and an end-to-end ODE solve.

---

## License

BSD-3-Clause
