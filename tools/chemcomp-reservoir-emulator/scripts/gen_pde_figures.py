"""Generate PDE comparison figures for all four PDE benchmarks.

Produces figures/pde_NAME/pde_comparison.pdf for:
  - heat_equation
  - advection_diffusion
  - fisher_equation
  - burgers_equation

Uses the high-complexity multi-subnet network (best for heat) and the
original network (best for Fisher/Advection), saving figures for both.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import chemical_reservoir as cr
from chemical_reservoir.io import load_reservoir_from_json

# Output directory
FIG_DIR = Path(__file__).parent.parent / "Manuscript" / "figures"

# Networks
NET_ORIGINAL = Path(__file__).parent.parent / "network_configs" / "high_dimensional_random.json"
NET_MULTI_HIGH = Path(__file__).parent.parent / "network_configs" / "multisubnet_n6_high.json"

# PDE problems
PDE_PROBLEMS = [
    ("heat_equation", cr.problems.create_heat_equation_problem(), "Heat Equation"),
    ("advection_diffusion", cr.problems.create_advection_diffusion_problem(), "Advection-Diffusion"),
    ("fisher_equation", cr.problems.create_fisher_equation_problem(), "Fisher Equation"),
    ("burgers_equation", cr.problems.create_burgers_equation_problem(), "Burgers' Equation"),
]

N_TRAIN_TIME = 150
N_TEST_TIME = 50


def plot_pde_comparison(results, problem_title, outpath):
    """Three-panel (x,t) heatmap: true | predicted | error."""
    x_grid = results["x_grid"]
    t_all = results["t_all"]
    sol_true = results["sol_true"]
    sol_pred = results["sol_pred"]
    n_train = len(results["t_train"])

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    
    # Determine common color scale for true/predicted
    vmin = min(sol_true.min(), sol_pred.min())
    vmax = max(sol_true.max(), sol_pred.max())
    
    for ax, data, title in [
        (axes[0], sol_true, "True Solution"),
        (axes[1], sol_pred, "Reservoir Prediction"),
    ]:
        im = ax.pcolormesh(x_grid, t_all, data, shading="auto",
                           vmin=vmin, vmax=vmax, cmap='viridis')
        ax.axhline(t_all[n_train - 1], color="r", ls="--", lw=1.5,
                   label='Train/Test boundary')
        ax.set(xlabel="x", ylabel="t", title=title)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        if ax == axes[0]:
            ax.legend(loc='upper right', fontsize=8)

    # Error panel with symmetric colormap
    error = sol_true - sol_pred
    err_max = max(abs(error.min()), abs(error.max()))
    im = axes[2].pcolormesh(x_grid, t_all, error, shading="auto",
                            vmin=-err_max, vmax=err_max, cmap='RdBu_r')
    axes[2].axhline(t_all[n_train - 1], color="r", ls="--", lw=1.5)
    axes[2].set(xlabel="x", ylabel="t", title="Pointwise Error")
    plt.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04)

    train_rmse = results['train_rmse']
    test_rmse = results['test_rmse']
    fig.suptitle(f"{problem_title}    "
                 f"(Train RMSE = {train_rmse:.2e},  Test RMSE = {test_rmse:.2e})",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath, bbox_inches="tight", dpi=150)
    fig.savefig(outpath.with_suffix('.pdf'), bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def run_pde(network_path, problem, n_train=N_TRAIN_TIME, n_test=N_TEST_TIME):
    """Run PDE prediction with given network."""
    reservoir, _ = load_reservoir_from_json(str(network_path))
    reservoir.feature_expansion = 'delayed'
    results = reservoir.solve_pde(problem, n_train_time=n_train, n_test_time=n_test)
    return results


if __name__ == "__main__":
    # Best network per problem (based on notebook results):
    # Heat: multi_high, Fisher: original, Advection: original, Burgers: multi_high
    best_networks = {
        "heat_equation": NET_MULTI_HIGH,
        "advection_diffusion": NET_ORIGINAL,
        "fisher_equation": NET_ORIGINAL,
        "burgers_equation": NET_MULTI_HIGH,
    }

    for dirname, problem, title in PDE_PROBLEMS:
        print(f"\n{title}:")
        net_path = best_networks[dirname]
        net_name = net_path.stem
        print(f"  Network: {net_name}")
        
        results = run_pde(net_path, problem)
        
        # Compute test R²
        n_train = len(results['t_train'])
        sol_true_test = results['sol_true'][n_train:]
        sol_pred_test = results['sol_pred'][n_train:]
        ss_res = np.sum((sol_true_test - sol_pred_test)**2)
        ss_tot = np.sum((sol_true_test - np.mean(sol_true_test))**2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        print(f"  Train RMSE: {results['train_rmse']:.4e}")
        print(f"  Test RMSE:  {results['test_rmse']:.4e}")
        print(f"  Test R²:    {r2:.4f}")
        
        outpath = FIG_DIR / dirname / "pde_comparison.png"
        plot_pde_comparison(results, title, outpath)
    
    print("\n\nAll PDE figures generated.")
