"""
Main analysis runner for ODE and PDE examples with varying horizons.
"""

import numpy as np
from typing import List, Dict
import os

from .config import set_output_directory, get_output_config
from .constraints import PhysicalConstraints, MeasurementConstraints
from .problems import (
    ODEProblem, PDEProblem1D,
    create_schnakenberg_problem, create_vanderpol_problem, create_lorenz_problem,
    create_brusselator_problem, create_lotka_volterra_problem,
    create_heat_equation_problem, create_advection_diffusion_problem,
    create_fisher_equation_problem, create_burgers_equation_problem
)
from .metrics import MetricsCalculator, LearningMetrics
from .sensitivity import SensitivityAnalyzer
from .visualization import ResultsVisualizer
from .reservoir import ConstrainedChemicalReservoir


def create_default_reactions():
    """Create default reaction network."""
    return [
        {"species": ["A", "B"], "reactants": [1, 1], "products": [2, 0], "rate": 1.2},
        {"species": ["B", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 0.9},
        {"species": ["C", "A"], "reactants": [1, 1], "products": [0, 2], "rate": 0.7},
        {"species": ["A", "B"], "reactants": [2, 0], "products": [1, 1], "rate": 0.15},
        {"species": ["B", "C"], "reactants": [0, 2], "products": [1, 1], "rate": 0.12},
        {"species": ["C", "A"], "reactants": [0, 2], "products": [1, 1], "rate": 0.09},
        {"species": ["A", "D"], "reactants": [1, 0], "products": [0, 1], "rate": 0.4},
        {"species": ["D", "E"], "reactants": [1, 0], "products": [0, 1], "rate": 0.35},
        {"species": ["E", "B"], "reactants": [1, 0], "products": [0, 1], "rate": 0.3},
    ]


def create_default_constraints():
    """Create default physical and measurement constraints."""
    physical = PhysicalConstraints(
        default_min_concentration=0.0,
        default_max_concentration=5.0,
        default_min_input_rate=0.0,
        default_max_input_rate=2.0
    )
    
    physical.add_species_bounds('A', min_conc=0, max_conc=3, min_rate=0.1, max_rate=1.5)
    physical.add_species_bounds('B', min_conc=0, max_conc=4, min_rate=0, max_rate=2.0)
    physical.add_species_bounds('C', min_conc=0, max_conc=2.5, min_rate=0, max_rate=1.0)
    physical.add_species_bounds('D', min_conc=0, max_conc=1.5)
    physical.add_species_bounds('E', min_conc=0, max_conc=0.8)
    
    measurement = MeasurementConstraints(
        max_simultaneous_measurements=3,
        measurable_species={'A', 'B', 'C', 'D'},
        unmeasurable_species={'E'},
        measurement_noise_std=0.01
    )
    
    return physical, measurement


def run_ode_analysis(problem: ODEProblem,
                    output_subdir: str,
                    horizons: List[float] = None,
                    n_train: int = 400,
                    n_test: int = 200,
                    run_optimization: bool = True,
                    run_sensitivity: bool = True,
                    bo_iterations: int = 30) -> Dict:
    """
    Run complete ODE analysis with varying horizons.
    
    Args:
        problem: ODE problem to solve
        output_subdir: Subdirectory for output
        horizons: List of test horizon factors to evaluate
        n_train: Training points
        n_test: Test points
        run_optimization: Whether to run Bayesian optimization
        run_sensitivity: Whether to run sensitivity analysis
        bo_iterations: Bayesian optimization iterations
    """
    if horizons is None:
        horizons = [0.25, 0.5, 0.75, 1.0]
    
    # Setup output
    config = get_output_config()
    subdir = config.get_subdir(output_subdir)
    from .config import OutputConfig
    local_config = OutputConfig(output_dir=subdir, figure_format="pdf",
                                save_figures=True, show_figures=False)
    
    print(f"\n{'='*70}")
    print(f"ODE ANALYSIS: {problem.name}")
    print(f"Output: {subdir}")
    print(f"{'='*70}")
    
    # Create reservoir
    reactions = create_default_reactions()
    physical, measurement = create_default_constraints()
    
    reservoir = ConstrainedChemicalReservoir(
        reactions=reactions,
        input_channels=['A', 'B', 'C'],
        output_channels=['A', 'B', 'D'],
        physical_constraints=physical,
        measurement_constraints=measurement,
        feature_expansion='quadratic'
    )
    
    # Plot reaction network
    print("\n[1] Plotting reaction network...")
    reservoir.plot_reaction_network(filename="reaction_network", config=local_config)
    
    # Run optimization
    opt_results = None
    if run_optimization:
        print(f"\n[2] Running Bayesian optimization ({bo_iterations} iterations)...")
        opt_results = reservoir.optimize_with_bayesian(
            problem, n_iterations=bo_iterations, n_initial=10, verbose=True
        )
        
        visualizer = ResultsVisualizer(config=local_config)
        visualizer.plot_optimization_analysis(opt_results, filename="optimization_analysis")
    else:
        # Set default values
        reservoir.set_initial_concentrations({sp: 1.0 for sp in reservoir.species_names})
        reservoir.set_input_rates({ch: 0.5 for ch in reservoir.input_channels})
    
    # Run for each horizon
    print(f"\n[3] Running analysis for {len(horizons)} horizons...")
    
    horizon_results = {
        'horizons': horizons,
        'train_rmse': [],
        'test_rmse': [],
        'train_r2': [],
        'test_r2': [],
        'all_results': []
    }
    
    visualizer = ResultsVisualizer(config=local_config)
    
    for i, horizon in enumerate(horizons):
        print(f"\n  Horizon {i+1}/{len(horizons)}: {horizon}")
        
        results = reservoir.solve_ode(
            problem, n_train=n_train, n_test=n_test,
            test_horizon_factor=horizon
        )
        
        metrics = MetricsCalculator.calculate_all(
            y_true=results['y_true'],
            y_pred=results['y_pred'],
            t=results['t_all'],
            n_train=len(results['t_train']),
            variable_names=problem.variable_names
        )
        
        horizon_results['train_rmse'].append(metrics.train_rmse)
        horizon_results['test_rmse'].append(metrics.test_rmse)
        horizon_results['train_r2'].append(metrics.train_r2)
        horizon_results['test_r2'].append(metrics.test_r2)
        horizon_results['all_results'].append({
            'results': results,
            'metrics': metrics,
            'horizon': horizon
        })
        
        # Plot for this horizon
        visualizer.plot_learning_metrics_dashboard(
            results, metrics, problem,
            filename=f"learning_metrics_horizon_{horizon:.2f}"
        )
        
        print(f"    Train RMSE: {metrics.train_rmse:.4e}, Test RMSE: {metrics.test_rmse:.4e}")
    
    # Plot horizon comparison
    print("\n[4] Plotting horizon comparison...")
    visualizer.plot_horizon_comparison(horizon_results, filename="horizon_comparison")
    
    # Sensitivity analysis (on best horizon)
    sensitivity_results = None
    if run_sensitivity:
        print("\n[5] Running sensitivity analysis...")
        
        best_horizon_idx = np.argmin(horizon_results['test_rmse'])
        best_horizon = horizons[best_horizon_idx]
        
        analyzer = SensitivityAnalyzer(reservoir, problem)
        
        print("  - Concentration sensitivity...")
        analyzer.analyze_concentration_sensitivity(n_train=200)
        
        print("  - Input rate sensitivity...")
        analyzer.analyze_input_rate_sensitivity(n_train=200)
        
        print("  - Measurement sensitivity...")
        analyzer.analyze_measurement_sensitivity(n_train=200)
        
        print("  - Local sensitivity...")
        analyzer.compute_local_sensitivity(n_train=200)
        
        sensitivity_results = analyzer.sensitivity_results
        visualizer.plot_sensitivity_analysis(sensitivity_results, filename="sensitivity_analysis")
    
    # Summary
    reservoir.report()
    
    print(f"\n{'='*70}")
    print(f"ODE ANALYSIS COMPLETE: {problem.name}")
    print(f"Best horizon: {horizons[np.argmin(horizon_results['test_rmse'])]:.2f}")
    print(f"Best test RMSE: {min(horizon_results['test_rmse']):.4e}")
    print(f"{'='*70}")
    
    return {
        'reservoir': reservoir,
        'horizon_results': horizon_results,
        'opt_results': opt_results,
        'sensitivity_results': sensitivity_results,
        'problem': problem
    }


def run_pde_analysis(problem: PDEProblem1D,
                    output_subdir: str,
                    horizons: List[float] = None,
                    n_train_time: int = 80,
                    n_test_time: int = 40) -> Dict:
    """
    Run complete PDE analysis with varying horizons.
    """
    if horizons is None:
        horizons = [0.15, 0.3, 0.5]
    
    config = get_output_config()
    subdir = config.get_subdir(output_subdir)
    from .config import OutputConfig
    local_config = OutputConfig(output_dir=subdir, figure_format="pdf",
                                save_figures=True, show_figures=False)
    
    print(f"\n{'='*70}")
    print(f"PDE ANALYSIS: {problem.name}")
    print(f"Output: {subdir}")
    print(f"{'='*70}")
    
    reactions = create_default_reactions()
    physical, _ = create_default_constraints()
    
    measurement = MeasurementConstraints(max_simultaneous_measurements=4)
    
    reservoir = ConstrainedChemicalReservoir(
        reactions=reactions,
        input_channels=['A', 'B'],
        output_channels=['A', 'B', 'C', 'D'],
        physical_constraints=physical,
        measurement_constraints=measurement
    )
    
    reservoir.set_initial_concentrations({sp: 1.0 for sp in reservoir.species_names})
    reservoir.set_input_rates({'A': 0.5, 'B': 0.8})
    
    # Plot network
    print("\n[1] Plotting reaction network...")
    reservoir.plot_reaction_network(filename="reaction_network", config=local_config)
    
    # Run for each horizon
    print(f"\n[2] Running analysis for {len(horizons)} horizons...")
    
    horizon_results = {
        'horizons': horizons,
        'train_rmse': [],
        'test_rmse': [],
        'train_r2': [],
        'test_r2': [],
        'all_results': []
    }
    
    visualizer = ResultsVisualizer(config=local_config)
    
    for i, horizon in enumerate(horizons):
        print(f"\n  Horizon {i+1}/{len(horizons)}: {horizon}")
        
        results = reservoir.solve_pde(
            problem, n_train_time=n_train_time, n_test_time=n_test_time,
            test_horizon_factor=horizon
        )
        
        pde_metrics = MetricsCalculator.calculate_pde_metrics(
            sol_true=results['sol_true'],
            sol_pred=results['sol_pred'],
            x=results['x_grid'],
            t=results['t_all'],
            n_train=len(results['t_train'])
        )
        
        horizon_results['train_rmse'].append(pde_metrics['train_rmse'])
        horizon_results['test_rmse'].append(pde_metrics['test_rmse'])
        horizon_results['train_r2'].append(pde_metrics['train_r2'])
        horizon_results['test_r2'].append(pde_metrics.get('test_r2', 0))
        horizon_results['all_results'].append({
            'results': results,
            'metrics': pde_metrics,
            'horizon': horizon
        })
        
        visualizer.plot_pde_metrics(results, pde_metrics,
                                   filename=f"pde_metrics_horizon_{horizon:.2f}")
        
        print(f"    Train RMSE: {pde_metrics['train_rmse']:.4e}, Test RMSE: {pde_metrics['test_rmse']:.4e}")
    
    # Plot comparison
    print("\n[3] Plotting horizon comparison...")
    visualizer.plot_horizon_comparison(horizon_results, filename="horizon_comparison")
    
    print(f"\n{'='*70}")
    print(f"PDE ANALYSIS COMPLETE: {problem.name}")
    print(f"Best horizon: {horizons[np.argmin(horizon_results['test_rmse'])]:.2f}")
    print(f"Best test RMSE: {min(horizon_results['test_rmse']):.4e}")
    print(f"{'='*70}")
    
    return {
        'reservoir': reservoir,
        'horizon_results': horizon_results,
        'problem': problem
    }


def run_complete_analysis(output_dir: str = "test_results",
                         run_ode: bool = True,
                         run_pde: bool = True,
                         ode_problems: List[str] = None,
                         pde_problems: List[str] = None,
                         horizons: List[float] = None):
    """
    Run complete analysis on multiple ODE and PDE problems.
    
    Args:
        output_dir: Base output directory
        run_ode: Whether to run ODE examples
        run_pde: Whether to run PDE examples
        ode_problems: List of ODE problem names to run
        pde_problems: List of PDE problem names to run
        horizons: List of test horizon factors
    """
    # Setup output
    set_output_directory(output_dir, figure_format="pdf", 
                        save_figures=True, show_figures=False)
    
    print("\n" + "#" * 70)
    print("# CHEMICAL RESERVOIR COMPUTING - COMPLETE ANALYSIS")
    print("#" * 70)
    print(f"\nOutput directory: {output_dir}")
    
    all_results = {'ode': {}, 'pde': {}}
    
    # ODE problems
    if run_ode:
        if ode_problems is None:
            ode_problems = ['schnakenberg', 'vanderpol', 'brusselator']
        
        ode_creators = {
            'schnakenberg': create_schnakenberg_problem,
            'vanderpol': create_vanderpol_problem,
            'lorenz': create_lorenz_problem,
            'brusselator': create_brusselator_problem,
            'lotka_volterra': create_lotka_volterra_problem
        }
        
        print(f"\n{'='*70}")
        print("ODE PROBLEMS")
        print(f"{'='*70}")
        
        for prob_name in ode_problems:
            if prob_name in ode_creators:
                problem = ode_creators[prob_name]()
                results = run_ode_analysis(
                    problem,
                    output_subdir=f"ode_{prob_name}",
                    horizons=horizons or [0.25, 0.5, 0.75, 1.0],
                    run_optimization=True,
                    run_sensitivity=True,
                    bo_iterations=25
                )
                all_results['ode'][prob_name] = results
    
    # PDE problems
    if run_pde:
        if pde_problems is None:
            pde_problems = ['heat', 'advection_diffusion']
        
        pde_creators = {
            'heat': create_heat_equation_problem,
            'advection_diffusion': create_advection_diffusion_problem,
            'fisher': create_fisher_equation_problem,
            'burgers': create_burgers_equation_problem
        }
        
        print(f"\n{'='*70}")
        print("PDE PROBLEMS")
        print(f"{'='*70}")
        
        for prob_name in pde_problems:
            if prob_name in pde_creators:
                problem = pde_creators[prob_name]()
                results = run_pde_analysis(
                    problem,
                    output_subdir=f"pde_{prob_name}",
                    horizons=horizons or [0.15, 0.3, 0.5]
                )
                all_results['pde'][prob_name] = results
    
    # Print summary
    print("\n" + "#" * 70)
    print("# ANALYSIS SUMMARY")
    print("#" * 70)
    
    if all_results['ode']:
        print("\nODE Results:")
        print(f"  {'Problem':<20} {'Best Horizon':>15} {'Best Test RMSE':>15}")
        print(f"  {'-'*50}")
        for name, res in all_results['ode'].items():
            hr = res['horizon_results']
            best_idx = np.argmin(hr['test_rmse'])
            print(f"  {name:<20} {hr['horizons'][best_idx]:>15.2f} {hr['test_rmse'][best_idx]:>15.4e}")
    
    if all_results['pde']:
        print("\nPDE Results:")
        print(f"  {'Problem':<20} {'Best Horizon':>15} {'Best Test RMSE':>15}")
        print(f"  {'-'*50}")
        for name, res in all_results['pde'].items():
            hr = res['horizon_results']
            best_idx = np.argmin(hr['test_rmse'])
            print(f"  {name:<20} {hr['horizons'][best_idx]:>15.2f} {hr['test_rmse'][best_idx]:>15.4e}")
    
    # List output files
    print(f"\nOutput files saved to: {output_dir}/")
    for root, dirs, files in os.walk(output_dir):
        level = root.replace(output_dir, '').count(os.sep)
        indent = '  ' * level
        print(f"{indent}{os.path.basename(root)}/")
        subindent = '  ' * (level + 1)
        for file in files:
            print(f"{subindent}{file}")
    
    print("\n" + "#" * 70)
    print("# COMPLETE")
    print("#" * 70)
    
    return all_results


# Main entry point
if __name__ == "__main__":
    results = run_complete_analysis(
        output_dir="test_results",
        run_ode=True,
        run_pde=True,
        ode_problems=['schnakenberg', 'vanderpol'],
        pde_problems=['heat', 'advection_diffusion'],
        horizons=[0.25, 0.5, 1.0]
    )