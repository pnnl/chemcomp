# ChemComp — Emulation Module

The emulation module of the **ChemComp** software package. Solves ordinary and
partial differential equations using chemical reaction networks as a reservoir
computing substrate.

---

## Algorithm overview

The emulator exploits the rich nonlinear dynamics of chemical reaction networks
(CRNs) to approximate the solutions of differential equations without
traditional numerical integration of the target system during prediction.

### Pipeline

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

### Step-by-step

1. **Target generation** — The true ODE/PDE solution is computed on a training
   time interval using `scipy.integrate.solve_ivp` (ODEs) or the method of
   lines with finite differences (PDEs).

2. **Input encoding** — The known solution is encoded into time-varying species
   production rates that externally drive the CRN.  Four encoding strategies
   are available:
   - `constant` — uniform injection rate
   - `sinusoidal` — harmonic modulation of rates
   - `mixed` (default) — combination of constant and sinusoidal components
   - `multiharmonic` — 8 Fourier harmonics for richer spectral content

3. **Reservoir dynamics** — The CRN is numerically integrated under mass-action
   kinetics.  Each species concentration c_i(t) evolves according to:

       dc_i/dt = Σ_j ν_ij · k_j · ∏_l c_l^α_lj  +  S_i(t)

   where ν_ij are stoichiometric coefficients, k_j are rate constants, and
   S_i(t) is the external input encoding.  The resulting concentration
   trajectories form a high-dimensional, nonlinear expansion of the input
   signal — this is the reservoir state matrix.

4. **Feature expansion** — The raw concentration time series can optionally be
   augmented:
   - `none` — use concentrations directly
   - `quadratic` — append all pairwise cross-products c_i · c_j
   - `delayed` (recommended) — append time-shifted copies at delays [1, 2, 4, 8]

5. **Readout training** — A ridge regression (`sklearn.linear_model.RidgeCV`)
   maps the feature matrix to the target solution values over the training
   window.  The regularisation parameter α is selected by leave-one-out
   cross-validation over 50 log-spaced candidates.

6. **Prediction** — The trained readout is applied to reservoir features on the
   test time interval to produce predictions of the ODE/PDE solution.

### Open-loop vs. closed-loop prediction

The emulator supports two fundamentally different prediction strategies:

#### Open-loop (`solve_ode` / `solve_pde`)

In open-loop mode the full input encoding is pre-computed and used to drive the
reservoir over **both** training and test intervals in a single forward pass.
The readout then maps reservoir states to the target solution.

```
  Known solution ──▶ encode full [0, T] ──▶ integrate CRN ──▶ readout
```

**Characteristics:**
- No feedback — the reservoir never sees its own predictions.
- The input signal for the test interval must still be derived from the true
  solution (or an approximation of it), so open-loop mode is best suited to
  **interpolation** or problems where the driving signal is known ahead of time
  (e.g., PDEs solved via method of lines where spatial modes are known).
- Numerically stable regardless of prediction quality.
- Typically achieves high accuracy for smooth, periodic systems and PDEs.

#### Closed-loop (`solve_ode_closed_loop`)

In closed-loop mode the reservoir is trained using **teacher forcing** (the
true solution drives the encoding during training).  At test time the true
solution is no longer available — instead, the reservoir's own predicted
outputs are fed back to modulate the input encoding at each time step:

```
  Training:  true y(t) ──▶ encode ──▶ CRN ──▶ readout ──▶ train weights

  Testing:   ŷ(t) ──▶ encode ──▶ CRN ──▶ readout ──▶ ŷ(t+Δt) ──┐
              ▲                                                    │
              └────────────────── feedback ─────────────────────────┘
```

**Characteristics:**
- Enables **autonomous prediction** beyond the training window — the
  reservoir generates its own driving signal, making it suitable for
  forecasting chaotic attractors (e.g., Lorenz system).
- A `feedback_scale` parameter (default 0.5) controls how strongly the
  predicted output modulates the input rates.
- More susceptible to error accumulation — small prediction errors are fed
  back and can amplify over time.
- Best suited for **oscillatory** and **chaotic** ODE systems where the
  attractor structure constrains the dynamics.

#### When to use which

| Scenario | Recommended mode |
|----------|-----------------|
| PDE solutions (heat, advection, Burgers, Fisher) | Open-loop |
| Periodic ODEs with known forcing | Open-loop |
| Autonomous ODE forecasting (no future forcing data) | Closed-loop |
| Chaotic attractors (Lorenz) | Closed-loop |
| Limit cycles (Van der Pol, Brusselator) | Closed-loop |

### Multi-subnetwork architecture

For challenging problems (high-dimensional chaotic systems, steep PDE
gradients), the emulator supports a **multi-subnetwork** architecture:

- Multiple weakly-coupled Brusselator modules share a common fuel species.
- Inter-module **leakage** with a tunable rate allows information flow between
  subnetworks without tight coupling.
- A constant **fuel feed** drives all modules, mimicking an open reactor.
- Scaling studies show that 4–6 subnetworks at standard complexity provides
  the best trade-off between expressiveness and trainability.

Network generators in `scripts/` produce these architectures parametrically.

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
networkx ≥ 2.6, pandas ≥ 1.3, sympy ≥ 1.10, pyyaml ≥ 6.0.  The `[notebooks]`
extra adds JupyterLab and ipywidgets.

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
# Solve the Van der Pol oscillator with the cascade network (closed-loop)
python run_reservoir.py --problem vanderpol --network network_configs/cascade.json

# Use open-loop mode instead
python run_reservoir.py --problem vanderpol --network network_configs/cascade.json --mode open-loop

# Solve a PDE (open-loop is automatic for PDEs)
python run_reservoir.py --problem heat_equation --network network_configs/dense.json -o my_results

# Use a custom problem defined in a YAML file
python run_reservoir.py --problem-file problem_configs/lorenz.yaml --network network_configs/cascade.json

# List all available problems and networks
python run_reservoir.py --list-problems
python run_reservoir.py --list-networks
```

**Driver options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--problem` | — | Built-in problem name (see `--list-problems`) |
| `--problem-file` | — | Path to a YAML file defining a custom ODE/PDE |
| `--network` | *(required)* | Path to a reaction-network JSON file |
| `-o` / `--output` | `results/<problem>_<network>_<timestamp>/` | Output directory |
| `--n-train` | 800 | Number of training time-points |
| `--n-test` | 400 | Number of test time-points |
| `--encoding` | `mixed` | Input encoding: `constant`, `sinusoidal`, `mixed`, `multiharmonic` |
| `--features` | `delayed` | Feature expansion: `none`, `quadratic`, `delayed` |
| `--mode` | `closed-loop` | Prediction mode: `open-loop` or `closed-loop` |
| `--feedback-scale` | `0.5` | Output-to-input feedback strength (closed-loop only) |

### Custom problem files (YAML + SymPy)

Instead of using built-in problems, you can define your own ODE or PDE in a
YAML file.  Equations are written as SymPy-compatible expressions and parsed
at load time.

**ODE example** (`problem_configs/lotka_volterra.yaml`):

```yaml
type: ode
name: LotkaVolterra
variables: [prey, predator]
parameters:
  alpha: 1.1
  beta: 0.4
  gamma: 0.4
  delta: 0.1
equations:
  - "alpha*prey - beta*prey*predator"
  - "delta*prey*predator - gamma*predator"
initial_conditions: [10.0, 5.0]
t_span: [0, 15]
```

**PDE example** (`problem_configs/burgers.yaml`):

```yaml
type: pde
name: BurgersEquation
parameters:
  nu: 0.01
equation: "-u*u_x + nu*u_xx"
x_domain: [0, 1]
t_span: [0, 0.5]
n_spatial: 50
initial_condition: "sin(2*pi*x)"
boundary_conditions:
  left: [dirichlet, 0]
  right: [dirichlet, 0]
```

**Supported expression features:** standard arithmetic (`+`, `-`, `*`, `/`,
`**`), functions (`sin`, `cos`, `exp`, `log`, `sqrt`, `tanh`, `abs`),
constants (`pi`, `e`), and implicit multiplication (`2x` → `2*x`).

**Python API:**

```python
from chemical_reservoir.custom_problems import load_problem_from_yaml

problem = load_problem_from_yaml("problem_configs/lorenz.yaml")
results = reservoir.solve_ode_closed_loop(problem, n_train=400, n_test=200)
```

---

## Project structure

```
chemcomp-emulation/
├── pyproject.toml                # Package metadata & build config (BSD-3-Clause)
├── environment.yml               # Conda environment (chemcomp)
├── README.md
├── LICENSE
├── run_reservoir.py              # CLI driver script
├── chemical_reservoir/           # Installable Python package
│   ├── __init__.py               # Public API & version (1.0.0)
│   ├── config.py                 # OutputConfig dataclass
│   ├── constraints.py            # Physical & measurement constraints
│   ├── custom_problems.py        # YAML + SymPy problem loader
│   ├── io.py                     # JSON load / save / validation for networks
│   ├── metrics.py                # RMSE, R², MAE, correlation metrics
│   ├── networks.py               # Built-in reaction network definitions
│   ├── network_viz.py            # Bipartite reaction-network graph visualisation
│   ├── optimization.py           # Bayesian optimiser (GP + acquisition functions)
│   ├── problems.py               # ODE & PDE problem dataclasses and factories
│   ├── reservoir.py              # ConstrainedChemicalReservoir (core engine)
│   ├── run_analysis.py           # High-level analysis pipelines
│   ├── sensitivity.py            # Sensitivity & perturbation analysis
│   └── visualization.py          # Multi-panel result dashboards
├── network_configs/              # JSON reaction-network definitions (18 networks)
│   ├── cascade.json
│   ├── complex_pathways.json
│   ├── dense.json
│   ├── feedback.json
│   ├── high_dimensional_random.json
│   ├── linear_chain.json
│   ├── minimal.json
│   ├── multi_subnet_leaky.json
│   ├── multisubnet_n{4,6,8,10,12}_*.json
│   ├── oscillatory.json
│   └── simple_autocatalytic.json
├── problem_configs/              # YAML custom problem definitions
│   ├── vanderpol.yaml
│   ├── lorenz.yaml
│   ├── lotka_volterra.yaml
│   ├── heat_equation.yaml
│   ├── burgers.yaml
│   └── fisher.yaml
├── notebooks/                    # Tutorial & demonstration notebooks
│   ├── 01_mathematical_foundations.ipynb
│   ├── 02_quickstart_tutorial.ipynb
│   ├── 03_advanced_examples.ipynb
│   ├── 04_api_reference.ipynb
│   └── 05_open_vs_closed_loop.ipynb
├── scripts/                      # Network generation utilities
│   ├── gen_highdim_network.py          # 20-species random Brusselator network
│   ├── gen_multisubnet_network.py      # 6-subnet leaky architecture
│   ├── gen_parameterized_network.py    # Parameterised width/depth generator
│   └── gen_pde_figures.py              # PDE comparison figure generation
└── tests/                        # Pytest test suite (9 tests)
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
ODE system.  The reservoir observes the full spatiotemporal field and the
readout reconstructs the solution at each grid point.

| Problem | Factory function | Type | Key physics |
|---------|-----------------|------|-------------|
| Heat equation | `create_heat_equation_problem()` | Parabolic | Diffusive decay |
| Advection–diffusion | `create_advection_diffusion_problem()` | Parabolic | Traveling wave + dispersion |
| Fisher equation | `create_fisher_equation_problem()` | Reaction–diffusion | Propagating front |
| Burgers equation | `create_burgers_equation_problem()` | Nonlinear hyperbolic | Shock/steep gradient formation |

---

## Reaction networks

### Built-in networks (defined in code)

Available via `chemical_reservoir.networks.list_builtin_networks()`.

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

### High-dimensional and multi-subnetwork networks (JSON)

| Network | Species | Reactions | Description |
|---------|---------|-----------|-------------|
| `high_dimensional_random` | 20 | 50 | 5 coupled Brusselator modules |
| `multi_subnet_leaky` | 25 | 78 | 6 Brusselator modules + leakage + fuel |
| `multisubnet_n4_standard` | 17 | 52 | 4 modules, standard complexity |
| `multisubnet_n6_high` | 31 | 114 | 6 modules, high complexity (extra intermediates) |
| `multisubnet_n10_high` | 51 | 190 | 10 modules, high complexity |
| `multisubnet_n12_standard` | 49 | 156 | 12 modules, standard complexity |

All networks are stored as JSON in `network_configs/` and loaded via
`cr.io.load_reservoir_from_json()`.  Custom networks can be defined by writing
a JSON file with the same schema (see the notebooks for examples).

---

## Notebooks

Five Jupyter notebooks progress from theory to advanced usage:

1. **Mathematical Foundations** — Reservoir computing theory, mass-action
   kinetics, ridge regression, and how they fit together.
2. **Quick-Start Tutorial** — Load a network, solve an ODE, evaluate metrics,
   and visualise predictions step by step.
3. **Advanced Examples** — PDE solving, Bayesian optimisation, sensitivity
   analysis, custom network design, and horizon comparison.
4. **API Reference** — Class-by-class walkthrough with code examples for every
   public component.
5. **Open vs. Closed Loop** — Systematic comparison of prediction modes across
   all ODE/PDE problems, multi-subnetwork scaling study, and training
   efficiency analysis.

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
