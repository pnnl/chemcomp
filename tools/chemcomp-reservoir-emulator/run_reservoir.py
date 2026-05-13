#!/usr/bin/env python
"""
Driver script for chemical reservoir computing.

Provide an ODE/PDE problem name and a reaction network JSON file.
Plots, metrics, and reservoir parameters are saved to an output directory.

Prediction modes
~~~~~~~~~~~~~~~~
- **open-loop** : the reservoir is driven by the input encoding over the
  full time domain; predictions are made in a single forward pass (no
  feedback).
- **closed-loop** (default): training uses teacher forcing; during testing,
  predicted outputs are fed back to modulate the input signal at each step,
  enabling autonomous attractor tracking.  The ``--feedback-scale`` flag
  controls the strength of the output-to-input feedback.

Usage
-----
    python run_reservoir.py --problem vanderpol --network network_configs/cascade.json
    python run_reservoir.py --problem lorenz --network network_configs/high_dimensional_random.json --mode closed-loop
    python run_reservoir.py --problem heat_equation --network network_configs/dense.json -o my_results
    python run_reservoir.py --list-problems
    python run_reservoir.py --list-networks

Options
-------
    --problem         Name of the built-in problem (see --list-problems)
    --network         Path to a reaction network JSON file
    -o / --output     Output directory (default: results/<problem>_<network>_<timestamp>)
    --n-train         Number of training points (default: 800)
    --n-test          Number of test points (default: 400)
    --encoding        Input encoding: constant, sinusoidal, mixed, multiharmonic (default: mixed)
    --features        Feature expansion: none, quadratic, delayed (default: delayed)
    --mode            Prediction mode: open-loop or closed-loop (default: closed-loop)
    --feedback-scale  Output-to-input feedback strength for closed-loop mode (default: 0.5)
    --list-problems   Print available problems and exit
    --list-networks   Print available network configs and exit
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")          # Non-interactive backend so plots save to disk without a display
import matplotlib.pyplot as plt

import chemical_reservoir as cr
from chemical_reservoir.io import load_reservoir_from_json
from chemical_reservoir.network_viz import ReactionNetworkVisualizer
from chemical_reservoir.visualization import ResultsVisualizer
from chemical_reservoir.problems import ODEProblem, PDEProblem1D

# ---- Problem registry --------------------------------------------------------
# Maps human-readable names to (type, factory_function) tuples.
# The type tag ("ODE" or "PDE") controls which solver and plotting path is used.

PROBLEMS = {
    # ODE problems — solved as initial-value problems
    "vanderpol":    ("ODE", cr.problems.create_vanderpol_problem),
    "schnakenberg": ("ODE", cr.problems.create_schnakenberg_problem),
    "lorenz":       ("ODE", cr.problems.create_lorenz_problem),
    "lotka_volterra": ("ODE", cr.problems.create_lotka_volterra_problem),
    "brusselator":  ("ODE", cr.problems.create_brusselator_problem),
    # PDE problems — semi-discretised in space via method of lines, then solved as ODE systems
    "heat_equation":        ("PDE", cr.problems.create_heat_equation_problem),
    "advection_diffusion":  ("PDE", cr.problems.create_advection_diffusion_problem),
    "fisher_equation":      ("PDE", cr.problems.create_fisher_equation_problem),
    "burgers_equation":     ("PDE", cr.problems.create_burgers_equation_problem),
}


def list_problems():
    """Print a formatted table of all registered problem names and types."""
    print("Available problems:")
    print(f"  {'Name':<25s} {'Type':<5s}")
    print("  " + "-" * 32)
    for name, (kind, _) in sorted(PROBLEMS.items()):
        print(f"  {name:<25s} {kind:<5s}")


def list_networks(config_dir: Path):
    """Scan *config_dir* for JSON files and print species/reaction counts."""
    print(f"Available networks in {config_dir}/:")
    for f in sorted(config_dir.glob("*.json")):
        with open(f) as fh:
            data = json.load(fh)
        # Collect unique species names across all reactions
        n_species = len({s for r in data["reactions"] for s in r["species"]})
        n_rxn = len(data["reactions"])
        print(f"  {f.name:<35s} {n_species:>2d} species, {n_rxn:>2d} reactions")


# ---- Plotting helpers --------------------------------------------------------

def save_ode_plots(results, metrics, problem, reservoir, outdir):
    """Generate and save all ODE-related plots.

    Produces:
      - prediction_overlay.pdf/.png — true vs predicted time series per variable
      - error_over_time.pdf — pointwise RMSE across the full time horizon
      - phase_portrait.pdf — 2-D phase plot (only for ≥2 variables)
      - reservoir_trajectories.pdf — individual species concentration traces
    """
    t_all = results["t_all"]
    y_true = results["y_true"]
    y_pred = results["y_pred"]
    n_train = len(results["t_train"])

    # 1. True vs Predicted overlay — one subplot per state variable
    fig, axes = plt.subplots(
        len(problem.variable_names), 1,
        figsize=(12, 3 * len(problem.variable_names)),
        sharex=True,
    )
    if len(problem.variable_names) == 1:
        axes = [axes]
    for i, name in enumerate(problem.variable_names):
        ax = axes[i]
        ax.plot(t_all, y_true[:, i], "b-", lw=1.5, label="True", alpha=0.8)
        ax.plot(t_all, y_pred[:, i], "r--", lw=1.5, label="Predicted", alpha=0.8)
        # Vertical line marks the boundary between training and test regions
        ax.axvline(t_all[n_train - 1], color="k", ls=":", lw=1, label="Train / Test")
        ax.set_ylabel(name)
        ax.legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("t")
    fig.suptitle(f"{problem.name}: True vs Predicted  (Test R²={metrics.test_r2:.4f})")
    plt.tight_layout()
    fig.savefig(outdir / "prediction_overlay.pdf", bbox_inches="tight")
    fig.savefig(outdir / "prediction_overlay.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 2. Error over time — log-scale RMSE at each time step (averaged over variables)
    fig, ax = plt.subplots(figsize=(12, 4))
    errors = np.sqrt(np.mean((y_true - y_pred) ** 2, axis=1))
    ax.semilogy(t_all, errors, "k-", lw=1)
    ax.axvline(t_all[n_train - 1], color="r", ls=":", label="Train / Test")
    ax.set(xlabel="t", ylabel="RMSE(t)", title=f"{problem.name}: Pointwise Error")
    ax.legend()
    plt.tight_layout()
    fig.savefig(outdir / "error_over_time.pdf", bbox_inches="tight")
    plt.close(fig)

    # 3. Phase portrait (if 2-D or 3-D)
    if len(problem.variable_names) >= 2:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.plot(y_true[:, 0], y_true[:, 1], "b-", lw=1.2, label="True", alpha=0.7)
        ax.plot(y_pred[:, 0], y_pred[:, 1], "r--", lw=1.2, label="Predicted", alpha=0.7)
        ax.set(xlabel=problem.variable_names[0], ylabel=problem.variable_names[1],
               title=f"{problem.name}: Phase Portrait")
        ax.legend()
        plt.tight_layout()
        fig.savefig(outdir / "phase_portrait.pdf", bbox_inches="tight")
        plt.close(fig)

    # 4. Reservoir species trajectories — shows the raw chemical concentrations
    #    that the readout layer learns from.  Useful for diagnosing whether the
    #    network dynamics are rich enough (limit cycles, transients, etc.).
    traj = results.get("trajectories")
    if traj is not None:
        n_sp = traj.shape[1]
        n_cols = 4
        n_rows = int(np.ceil(n_sp / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 2.5 * n_rows), sharex=True)
        axes = np.atleast_2d(axes)
        for i in range(n_sp):
            r, c = divmod(i, n_cols)
            axes[r, c].plot(t_all, traj[:, i], lw=0.8)
            axes[r, c].set_title(reservoir.species_names[i], fontsize=9)
        for i in range(n_sp, n_rows * n_cols):
            r, c = divmod(i, n_cols)
            axes[r, c].set_visible(False)
        fig.suptitle("Reservoir Species Trajectories")
        plt.tight_layout()
        fig.savefig(outdir / "reservoir_trajectories.pdf", bbox_inches="tight")
        plt.close(fig)


def save_pde_plots(results, problem, outdir):
    """Generate and save PDE-related plots.

    Produces pde_comparison.pdf/.png — a three-panel heatmap showing the
    true solution, predicted solution, and their pointwise error over the
    (x, t) domain.  A red dashed line marks the train/test time boundary.
    """
    x_grid = results["x_grid"]
    t_all = results["t_all"]
    sol_true = results["sol_true"]
    sol_pred = results["sol_pred"]
    n_train = len(results["t_train"])

    # Three-panel (x, t) heatmap: true | predicted | error
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, data, title in [
        (axes[0], sol_true, "True"),
        (axes[1], sol_pred, "Predicted"),
        (axes[2], sol_true - sol_pred, "Error"),
    ]:
        im = ax.pcolormesh(x_grid, t_all, data, shading="auto")
        ax.axhline(t_all[n_train - 1], color="r", ls=":", lw=1)  # train/test boundary
        ax.set(xlabel="x", ylabel="t", title=f"{problem.name}: {title}")
        plt.colorbar(im, ax=ax)
    plt.tight_layout()
    fig.savefig(outdir / "pde_comparison.pdf", bbox_inches="tight")
    fig.savefig(outdir / "pde_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_network_plots(reservoir, outdir):
    """Save a bipartite graph visualisation of the reaction network to PDF."""
    from chemical_reservoir.config import OutputConfig
    config = OutputConfig(
        output_dir=str(outdir),
        save_figures=True,
        show_figures=False,
        figure_format="pdf",
    )
    reservoir.plot_reaction_network(filename="reaction_network", config=config)


# ---- Main -------------------------------------------------------------------

def main():
    # ---- CLI argument parsing ----
    parser = argparse.ArgumentParser(
        description="Chemical Reservoir Computing Driver",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--problem", type=str, help="Problem name (see --list-problems)")
    parser.add_argument("--problem-file", type=str, default=None,
                        help="Path to a YAML file defining a custom ODE/PDE problem")
    parser.add_argument("--config", type=str, default=None,
                        help="Path to a YAML config file (overrides default_config.yaml)")
    parser.add_argument("--network", type=str, help="Path to reaction network JSON")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output directory (default: results/<problem>_<network>_<timestamp>)")
    parser.add_argument("--n-train", type=int, default=800,
                        help="Number of training time-points (default: 800)")
    parser.add_argument("--n-test", type=int, default=400,
                        help="Number of test time-points (default: 400)")
    parser.add_argument("--encoding", type=str, default="mixed",
                        choices=["constant", "sinusoidal", "mixed", "multiharmonic"],
                        help="How the ODE solution is encoded as reservoir input")
    parser.add_argument("--features", type=str, default="delayed",
                        choices=["none", "quadratic", "delayed"],
                        help="Feature expansion applied to reservoir concentrations")
    parser.add_argument("--mode", type=str, default="closed-loop",
                        choices=["open-loop", "closed-loop"],
                        help="Prediction mode: open-loop or closed-loop (default: closed-loop)")
    parser.add_argument("--feedback-scale", type=float, default=0.5,
                        help="Output-to-input feedback strength for closed-loop mode (default: 0.5)")
    parser.add_argument("--list-problems", action="store_true", help="List available problems")
    parser.add_argument("--list-networks", action="store_true", help="List network configs")
    args = parser.parse_args()

    # Resolve paths relative to the script location so the driver works from
    # any working directory.  network_configs/ lives at the project root.
    script_dir = Path(__file__).resolve().parent
    config_dir = script_dir / "network_configs"

    # Load user config (if supplied) before doing anything else
    if args.config:
        cfg_path = Path(args.config)
        if not cfg_path.is_absolute():
            cfg_path = script_dir / cfg_path
        cr.load_config(str(cfg_path))
        print(f"Config:   {cfg_path}")

    if args.list_problems:
        list_problems()
        return
    if args.list_networks:
        list_networks(config_dir)
        return

    if not args.problem and not args.problem_file:
        parser.error("--problem or --problem-file is required (use --list-problems)")
    if not args.network:
        parser.error("--network is required (use --list-networks)")

    # Resolve the network JSON path; relative paths are anchored to the project root
    network_path = Path(args.network)
    if not network_path.is_absolute():
        network_path = script_dir / network_path
    if not network_path.exists():
        print(f"Network file not found: {network_path}")
        sys.exit(1)

    # ---- Build output directory ----
    # Default: results/<problem>_<network-stem>_<timestamp>/
    if args.output:
        outdir = Path(args.output)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        net_stem = network_path.stem
        prob_name = args.problem or Path(args.problem_file).stem
        outdir = script_dir / "results" / f"{prob_name}_{net_stem}_{ts}"
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Create problem ----
    if args.problem_file:
        from chemical_reservoir.custom_problems import load_problem_from_yaml
        prob_path = Path(args.problem_file)
        if not prob_path.is_absolute():
            prob_path = script_dir / prob_path
        if not prob_path.exists():
            print(f"Problem file not found: {prob_path}")
            sys.exit(1)
        problem = load_problem_from_yaml(prob_path)
        kind = "PDE" if hasattr(problem, 'pde_func') else "ODE"
    else:
        if args.problem not in PROBLEMS:
            print(f"Unknown problem '{args.problem}'.")
            list_problems()
            sys.exit(1)
        kind, factory = PROBLEMS[args.problem]
        problem = factory()
    print(f"Problem:  {problem.name}  ({kind})")
    print(f"Network:  {network_path.name}")
    print(f"Output:   {outdir}")

    # ---- Load reservoir ----
    # Instantiates a ConstrainedChemicalReservoir from the reaction network JSON.
    # The feature expansion (none / quadratic / delayed) is set after construction.
    reservoir, net_config = load_reservoir_from_json(str(network_path))
    reservoir.feature_expansion = args.features
    print(f"Reservoir: {reservoir.n_species} species, {reservoir.n_reactions} reactions, "
          f"{len(reservoir.output_channels)} outputs")
    print(f"Encoding: {args.encoding},  Feature expansion: {args.features},  Mode: {args.mode}")

    # ---- Run reservoir ----
    # solve_ode / solve_pde drives the chemical network with the target solution,
    # collects concentration features, and trains a ridge-regression readout.
    # In closed-loop mode, the test phase feeds predictions back as input.
    if kind == "ODE":
        if args.mode == "closed-loop":
            results = reservoir.solve_ode_closed_loop(
                problem,
                n_train=args.n_train,
                n_test=args.n_test,
                encoding=args.encoding,
                feedback_scale=args.feedback_scale,
            )
        else:
            results = reservoir.solve_ode(
                problem,
                n_train=args.n_train,
                n_test=args.n_test,
                encoding=args.encoding,
            )
        metrics = cr.MetricsCalculator.calculate_all(
            y_true=results["y_true"],
            y_pred=results["y_pred"],
            t=results["t_all"],
            n_train=len(results["t_train"]),
            variable_names=problem.variable_names,
        )
        print(f"\n--- Results ---")
        print(f"Train RMSE: {metrics.train_rmse:.6e}")
        print(f"Test  RMSE: {metrics.test_rmse:.6e}")
        print(f"Test  R²:   {metrics.test_r2:.4f}")

        save_ode_plots(results, metrics, problem, reservoir, outdir)

        # Collect all scalar metrics into a dict for JSON serialisation
        metrics_dict = {
            "problem": problem.name,
            "network": net_config["name"],
            "n_train": args.n_train,
            "n_test": args.n_test,
            "encoding": args.encoding,
            "feature_expansion": args.features,
            "mode": args.mode,
            "feedback_scale": args.feedback_scale if args.mode == "closed-loop" else None,
            "n_species": reservoir.n_species,
            "n_reactions": reservoir.n_reactions,
            "n_features": reservoir.training_info.get("n_features"),
            "ridge_alpha": reservoir.training_info.get("ridge_alpha"),
            "train_rmse": float(metrics.train_rmse),
            "test_rmse": float(metrics.test_rmse),
            "test_r2": float(metrics.test_r2),
            "train_mae": float(metrics.train_mae),
            "test_mae": float(metrics.test_mae),
        }

    else:  # PDE
        results = reservoir.solve_pde(
            problem,
            n_train_time=args.n_train,
            n_test_time=args.n_test,
        )
        print(f"\n--- Results ---")
        print(f"Train RMSE: {results['train_rmse']:.6e}")
        print(f"Test  RMSE: {results['test_rmse']:.6e}")

        save_pde_plots(results, problem, outdir)

        metrics_dict = {
            "problem": problem.name,
            "network": net_config["name"],
            "n_train": args.n_train,
            "n_test": args.n_test,
            "encoding": args.encoding,
            "feature_expansion": args.features,
            "n_species": reservoir.n_species,
            "n_reactions": reservoir.n_reactions,
            "n_features": reservoir.training_info.get("n_features"),
            "ridge_alpha": reservoir.training_info.get("ridge_alpha"),
            "train_rmse": float(results["train_rmse"]),
            "test_rmse": float(results["test_rmse"]),
        }

    # ---- Save outputs ----
    # Network topology graph
    save_network_plots(reservoir, outdir)

    # Scalar metrics (RMSE, R², etc.) for downstream analysis
    with open(outdir / "metrics.json", "w") as f:
        json.dump(metrics_dict, f, indent=2)

    # Copy the exact network config used, so the run is fully reproducible
    with open(network_path) as f:
        net_data = json.load(f)
    with open(outdir / "network_config.json", "w") as f:
        json.dump(net_data, f, indent=2)

    # Reservoir hyperparameters and topology summary
    reservoir_info = {
        "species": reservoir.species_names,
        "n_species": reservoir.n_species,
        "n_reactions": reservoir.n_reactions,
        "input_channels": reservoir.input_channels,
        "output_channels": reservoir.output_channels,
        "feature_expansion": args.features,
        "encoding": args.encoding,
        "ridge_alpha": float(reservoir.training_info.get("ridge_alpha", 0)),
        "n_features": reservoir.training_info.get("n_features", 0),
    }
    with open(outdir / "reservoir_params.json", "w") as f:
        json.dump(reservoir_info, f, indent=2)

    print(f"\nAll outputs saved to: {outdir}")
    print(f"  prediction_overlay.pdf / .png")
    print(f"  error_over_time.pdf")
    if kind == "ODE" and len(problem.variable_names) >= 2:
        print(f"  phase_portrait.pdf")
    if kind == "PDE":
        print(f"  pde_comparison.pdf / .png")
    print(f"  reservoir_trajectories.pdf")
    print(f"  reaction_network*.pdf")
    print(f"  metrics.json")
    print(f"  reservoir_params.json")
    print(f"  network_config.json")


if __name__ == "__main__":
    main()
