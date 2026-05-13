"""
Visualization module for reservoir computing results.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from typing import Dict, Optional, Tuple

from .config import get_output_config, OutputConfig
from .metrics import LearningMetrics


def save_figure(fig, filename: str, config: OutputConfig = None):
    """Save figure to configured directory."""
    if config is None:
        config = get_output_config()
    
    if config.save_figures:
        filepath = config.get_path(filename)
        fig.savefig(filepath, format=config.figure_format, 
                   dpi=config.figure_dpi, bbox_inches='tight')
        print(f"  Saved: {filepath}")
    
    if config.show_figures:
        plt.show()
    else:
        plt.close(fig)


class ResultsVisualizer:
    """Comprehensive visualization for reservoir computing results."""
    
    def __init__(self, figsize_scale: float = 1.0, config: OutputConfig = None):
        self.figsize_scale = figsize_scale
        self.colors = plt.cm.tab10.colors
        self.config = config or get_output_config()
    
    def plot_learning_metrics_dashboard(self, 
                                       results: Dict,
                                       metrics: LearningMetrics,
                                       problem,
                                       filename: str = "learning_metrics_dashboard"):
        """Create comprehensive dashboard of learning metrics."""
        fig = plt.figure(figsize=(16 * self.figsize_scale, 12 * self.figsize_scale))
        gs = GridSpec(3, 4, figure=fig, hspace=0.3, wspace=0.3)
        
        t_all = results['t_all']
        y_true = results['y_true']
        y_pred = results['y_pred']
        n_train = len(results['t_train'])
        var_names = problem.variable_names
        
        # 1. Solution comparison
        ax1 = fig.add_subplot(gs[0, :2])
        for i, name in enumerate(var_names):
            ax1.plot(t_all, y_true[:, i], '--', color=self.colors[i], 
                    lw=2, label=f'{name} (true)')
            ax1.plot(t_all, y_pred[:, i], '-', color=self.colors[i], 
                    alpha=0.7, label=f'{name} (pred)')
        ax1.axvline(results['t_train'][-1], color='red', ls=':', lw=2)
        ax1.axvspan(t_all[0], results['t_train'][-1], alpha=0.1, color='blue')
        ax1.set_xlabel('Time', fontsize=12)
        ax1.set_ylabel('Solution', fontsize=12)
        ax1.set_title('Solution Comparison', fontsize=14, fontweight='bold')
        ax1.legend(loc='best', fontsize=9)
        ax1.grid(alpha=0.3)
        
        # 2. Error over time
        ax2 = fig.add_subplot(gs[0, 2:])
        errors = np.sqrt(np.sum((y_true - y_pred)**2, axis=1))
        ax2.semilogy(t_all, errors + 1e-10, 'b-', lw=2)
        ax2.axvline(results['t_train'][-1], color='red', ls=':', lw=2)
        ax2.fill_between(t_all[:n_train], 1e-10, errors[:n_train], alpha=0.3, color='blue', label='Training')
        ax2.fill_between(t_all[n_train:], 1e-10, errors[n_train:], alpha=0.3, color='orange', label='Testing')
        ax2.set_xlabel('Time', fontsize=12)
        ax2.set_ylabel('RMSE', fontsize=12)
        ax2.set_title('Error Evolution', fontsize=14, fontweight='bold')
        ax2.legend()
        ax2.grid(alpha=0.3)
        
        # 3. Metrics bar chart
        ax3 = fig.add_subplot(gs[1, 0])
        metric_names = ['RMSE', 'MAE', 'NRMSE']
        train_vals = [metrics.train_rmse, metrics.train_mae, metrics.train_nrmse]
        test_vals = [metrics.test_rmse, metrics.test_mae, metrics.test_nrmse]
        
        x = np.arange(len(metric_names))
        width = 0.35
        ax3.bar(x - width/2, train_vals, width, label='Train', color='steelblue')
        ax3.bar(x + width/2, test_vals, width, label='Test', color='darkorange')
        ax3.set_xticks(x)
        ax3.set_xticklabels(metric_names)
        ax3.set_ylabel('Value', fontsize=12)
        ax3.set_title('Error Metrics', fontsize=14, fontweight='bold')
        ax3.legend()
        ax3.grid(alpha=0.3, axis='y')
        
        # 4. R² and correlation
        ax4 = fig.add_subplot(gs[1, 1])
        corr_names = ['R²', 'Pearson']
        train_corr = [metrics.train_r2, metrics.train_pearson]
        test_corr = [metrics.test_r2, metrics.test_pearson]
        
        ax4.bar(np.arange(2) - width/2, train_corr, width, label='Train', color='steelblue')
        ax4.bar(np.arange(2) + width/2, test_corr, width, label='Test', color='darkorange')
        ax4.set_xticks(np.arange(2))
        ax4.set_xticklabels(corr_names)
        ax4.set_ylabel('Value', fontsize=12)
        ax4.set_title('Correlation Metrics', fontsize=14, fontweight='bold')
        ax4.set_ylim([0, 1.1])
        ax4.legend()
        ax4.grid(alpha=0.3, axis='y')
        
        # 5. Per-variable R²
        ax5 = fig.add_subplot(gs[1, 2])
        if metrics.per_variable_r2:
            vars_list = list(metrics.per_variable_r2.keys())
            r2_vals = [metrics.per_variable_r2[v] for v in vars_list]
            ax5.bar(vars_list, r2_vals, color=[self.colors[i] for i in range(len(vars_list))])
            ax5.set_ylabel('R²', fontsize=12)
            ax5.set_title('Per-Variable R²', fontsize=14, fontweight='bold')
            ax5.set_ylim([0, 1.1])
            ax5.grid(alpha=0.3, axis='y')
        
        # 6. Scatter
        ax6 = fig.add_subplot(gs[1, 3])
        ax6.scatter(y_true[:n_train].flatten(), y_pred[:n_train].flatten(), 
                   alpha=0.3, s=10, c='blue', label='Train')
        ax6.scatter(y_true[n_train:].flatten(), y_pred[n_train:].flatten(), 
                   alpha=0.3, s=10, c='orange', label='Test')
        lims = [min(y_true.min(), y_pred.min()), max(y_true.max(), y_pred.max())]
        ax6.plot(lims, lims, 'k--', lw=2)
        ax6.set_xlabel('True', fontsize=12)
        ax6.set_ylabel('Predicted', fontsize=12)
        ax6.set_title('Prediction Accuracy', fontsize=14, fontweight='bold')
        ax6.legend()
        ax6.grid(alpha=0.3)
        
        # 7. Phase space
        ax7 = fig.add_subplot(gs[2, 0])
        if len(var_names) >= 2:
            ax7.plot(y_true[:, 0], y_true[:, 1], 'b-', lw=2, alpha=0.7, label='True')
            ax7.plot(y_pred[:, 0], y_pred[:, 1], 'r--', lw=2, alpha=0.7, label='Predicted')
            ax7.set_xlabel(var_names[0], fontsize=12)
            ax7.set_ylabel(var_names[1], fontsize=12)
            ax7.set_title('Phase Space', fontsize=14, fontweight='bold')
            ax7.legend(fontsize=9)
            ax7.grid(alpha=0.3)
        
        # 8. Error distribution
        ax8 = fig.add_subplot(gs[2, 1])
        all_errors = (y_true - y_pred).flatten()
        ax8.hist(all_errors, bins=50, density=True, alpha=0.7, color='steelblue', edgecolor='black')
        ax8.axvline(0, color='red', ls='--', lw=2)
        ax8.set_xlabel('Prediction Error', fontsize=12)
        ax8.set_ylabel('Density', fontsize=12)
        ax8.set_title('Error Distribution', fontsize=14, fontweight='bold')
        ax8.grid(alpha=0.3)
        
        # 9. Per-variable RMSE
        ax9 = fig.add_subplot(gs[2, 2])
        if metrics.per_variable_train_rmse:
            vars_list = list(metrics.per_variable_train_rmse.keys())
            train_rmse = [metrics.per_variable_train_rmse[v] for v in vars_list]
            test_rmse = [metrics.per_variable_test_rmse.get(v, 0) for v in vars_list]
            
            x = np.arange(len(vars_list))
            ax9.bar(x - width/2, train_rmse, width, label='Train', color='steelblue')
            ax9.bar(x + width/2, test_rmse, width, label='Test', color='darkorange')
            ax9.set_xticks(x)
            ax9.set_xticklabels(vars_list)
            ax9.set_ylabel('RMSE', fontsize=12)
            ax9.set_title('Per-Variable RMSE', fontsize=14, fontweight='bold')
            ax9.legend()
            ax9.grid(alpha=0.3, axis='y')
        
        # 10. Summary
        ax10 = fig.add_subplot(gs[2, 3])
        ax10.axis('off')
        
        summary_text = (
            f"SUMMARY\n"
            f"{'='*30}\n\n"
            f"Train Points: {n_train}\n"
            f"Test Points: {len(t_all) - n_train}\n\n"
            f"Train RMSE: {metrics.train_rmse:.4e}\n"
            f"Test RMSE: {metrics.test_rmse:.4e}\n\n"
            f"Train R²: {metrics.train_r2:.4f}\n"
            f"Test R²: {metrics.test_r2:.4f}\n\n"
            f"Stability: {metrics.prediction_stability:.4f}"
        )
        
        ax10.text(0.1, 0.9, summary_text, transform=ax10.transAxes,
                 fontsize=11, verticalalignment='top', fontfamily='monospace',
                 bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.suptitle(f'Learning Metrics: {problem.name}', fontsize=16, fontweight='bold', y=0.98)
        
        save_figure(fig, filename, self.config)
        
        return fig
    
    def plot_sensitivity_analysis(self, sensitivity_results: Dict,
                                 filename: str = "sensitivity_analysis"):
        """Plot sensitivity analysis results."""
        fig = plt.figure(figsize=(16 * self.figsize_scale, 12 * self.figsize_scale))
        gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)
        
        # 1. Concentration sensitivity
        if 'concentration' in sensitivity_results:
            conc_res = sensitivity_results['concentration']
            
            ax1 = fig.add_subplot(gs[0, 0])
            fracs = conc_res['perturbation_fractions']
            rmse = conc_res['rmse_vs_perturbation']
            rmse_std = conc_res['rmse_std_vs_perturbation']
            
            ax1.errorbar(np.array(fracs) * 100, rmse, yerr=rmse_std, 
                        fmt='o-', capsize=5, color='steelblue', lw=2)
            ax1.set_xlabel('Perturbation (%)', fontsize=12)
            ax1.set_ylabel('Test RMSE', fontsize=12)
            ax1.set_title('Global Concentration Sensitivity', fontsize=14, fontweight='bold')
            ax1.grid(alpha=0.3)
            
            ax2 = fig.add_subplot(gs[0, 1])
            species = conc_res['species_names']
            for i, sp in enumerate(species):
                ax2.plot(np.array(fracs) * 100, conc_res['per_species_sensitivity'][sp],
                        'o-', label=sp, color=self.colors[i % len(self.colors)])
            ax2.set_xlabel('Perturbation (%)', fontsize=12)
            ax2.set_ylabel('Test RMSE', fontsize=12)
            ax2.set_title('Per-Species Sensitivity', fontsize=14, fontweight='bold')
            ax2.legend(fontsize=9)
            ax2.grid(alpha=0.3)
        
        # 2. Input rate sensitivity
        if 'input_rate' in sensitivity_results:
            rate_res = sensitivity_results['input_rate']
            
            ax3 = fig.add_subplot(gs[0, 2])
            fracs = rate_res['perturbation_fractions']
            rmse = rate_res['rmse_vs_perturbation']
            rmse_std = rate_res['rmse_std_vs_perturbation']
            
            ax3.errorbar(np.array(fracs) * 100, rmse, yerr=rmse_std,
                        fmt='s-', capsize=5, color='darkorange', lw=2)
            ax3.set_xlabel('Perturbation (%)', fontsize=12)
            ax3.set_ylabel('Test RMSE', fontsize=12)
            ax3.set_title('Input Rate Sensitivity', fontsize=14, fontweight='bold')
            ax3.grid(alpha=0.3)
        
        # 3. Measurement sensitivity
        if 'measurement' in sensitivity_results:
            meas_res = sensitivity_results['measurement']
            
            ax4 = fig.add_subplot(gs[1, 0])
            n_meas = meas_res['n_measurements']
            rmse_vals = meas_res['rmse_values']
            
            ax4.scatter(n_meas, rmse_vals, c='purple', s=80, alpha=0.7)
            ax4.set_xlabel('Number of Measurements', fontsize=12)
            ax4.set_ylabel('Test RMSE', fontsize=12)
            ax4.set_title('Measurement Count vs Performance', fontsize=14, fontweight='bold')
            ax4.grid(alpha=0.3)
            
            if 'best_rmse' in meas_res:
                ax4.axhline(meas_res['best_rmse'], color='green', ls='--')
                ax4.axhline(meas_res['worst_rmse'], color='red', ls='--')
            
            ax5 = fig.add_subplot(gs[1, 1])
            sets_str = [','.join(sorted(s)) for s in meas_res['measurement_sets'][:15]]
            ax5.barh(range(len(sets_str)), meas_res['rmse_values'][:15], color='teal')
            ax5.set_yticks(range(len(sets_str)))
            ax5.set_yticklabels(sets_str, fontsize=8)
            ax5.set_xlabel('Test RMSE', fontsize=12)
            ax5.set_title('Measurement Set Comparison', fontsize=14, fontweight='bold')
            ax5.grid(alpha=0.3, axis='x')
        
        # 4. Local sensitivity
        if 'local' in sensitivity_results:
            local_res = sensitivity_results['local']
            
            ax6 = fig.add_subplot(gs[1, 2])
            params = list(local_res['sensitivities'].keys())
            sens_vals = [local_res['sensitivities'][p] for p in params]
            
            colors_local = ['red' if v > 0 else 'blue' for v in sens_vals]
            ax6.barh(range(len(params)), sens_vals, color=colors_local, alpha=0.7)
            ax6.set_yticks(range(len(params)))
            ax6.set_yticklabels(params, fontsize=9)
            ax6.axvline(0, color='black', lw=1)
            ax6.set_xlabel('Sensitivity', fontsize=12)
            ax6.set_title('Local Sensitivity', fontsize=14, fontweight='bold')
            ax6.grid(alpha=0.3, axis='x')
            
            ax7 = fig.add_subplot(gs[2, 0])
            norm_sens = [local_res['normalized_sensitivities'][p] for p in params]
            abs_norm_sens = np.abs(norm_sens)
            sorted_idx = np.argsort(abs_norm_sens)[::-1]
            
            ax7.bar(range(len(params)), [abs_norm_sens[i] for i in sorted_idx], color='teal', alpha=0.7)
            ax7.set_xticks(range(len(params)))
            ax7.set_xticklabels([params[i] for i in sorted_idx], rotation=45, ha='right', fontsize=8)
            ax7.set_ylabel('|Normalized Sensitivity|', fontsize=12)
            ax7.set_title('Parameter Importance', fontsize=14, fontweight='bold')
            ax7.grid(alpha=0.3, axis='y')
        
        # Summary heatmap
        ax8 = fig.add_subplot(gs[2, 1:])
        
        if 'local' in sensitivity_results:
            local_res = sensitivity_results['local']
            params = list(local_res['sensitivities'].keys())
            
            summary_data = []
            for p in params:
                row = [
                    abs(local_res['sensitivities'].get(p, 0)),
                    abs(local_res['normalized_sensitivities'].get(p, 0))
                ]
                summary_data.append(row)
            
            summary_data = np.array(summary_data)
            for j in range(summary_data.shape[1]):
                col_max = summary_data[:, j].max()
                if col_max > 0:
                    summary_data[:, j] /= col_max
            
            im = ax8.imshow(summary_data, aspect='auto', cmap='YlOrRd')
            ax8.set_xticks(range(2))
            ax8.set_xticklabels(['Local', 'Normalized'], fontsize=10)
            ax8.set_yticks(range(len(params)))
            ax8.set_yticklabels(params, fontsize=9)
            ax8.set_title('Sensitivity Summary', fontsize=14, fontweight='bold')
            
            divider = make_axes_locatable(ax8)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            plt.colorbar(im, cax=cax)
        
        plt.suptitle('Sensitivity Analysis Dashboard', fontsize=16, fontweight='bold', y=0.98)
        
        save_figure(fig, filename, self.config)
        
        return fig
    
    def plot_pde_metrics(self, results: Dict, pde_metrics: Dict,
                        filename: str = "pde_metrics"):
        """Plot PDE-specific metrics."""
        fig = plt.figure(figsize=(16 * self.figsize_scale, 10 * self.figsize_scale))
        gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
        
        x = results['x_grid']
        t = results['t_all']
        n_train = len(results['t_train'])
        sol_true = results['sol_true']
        sol_pred = results['sol_pred']
        
        # Get problem name safely
        problem = results.get('problem')
        if problem is not None:
            problem_name = getattr(problem, 'name', 'PDE')
        else:
            problem_name = results.get('problem_name', 'PDE')
        
        # 1. True solution
        ax1 = fig.add_subplot(gs[0, 0])
        im1 = ax1.pcolormesh(x, t, sol_true, shading='auto', cmap='viridis')
        ax1.axhline(results['t_train'][-1], color='red', ls='--', lw=2)
        ax1.set_xlabel('x', fontsize=12)
        ax1.set_ylabel('t', fontsize=12)
        ax1.set_title('True Solution', fontsize=14, fontweight='bold')
        plt.colorbar(im1, ax=ax1)
        
        # 2. Predicted solution
        ax2 = fig.add_subplot(gs[0, 1])
        im2 = ax2.pcolormesh(x, t, sol_pred, shading='auto', cmap='viridis')
        ax2.axhline(results['t_train'][-1], color='red', ls='--', lw=2)
        ax2.set_xlabel('x', fontsize=12)
        ax2.set_ylabel('t', fontsize=12)
        ax2.set_title('Predicted Solution', fontsize=14, fontweight='bold')
        plt.colorbar(im2, ax=ax2)
        
        # 3. Error
        ax3 = fig.add_subplot(gs[0, 2])
        error = np.abs(sol_true - sol_pred)
        im3 = ax3.pcolormesh(x, t, error, shading='auto', cmap='hot')
        ax3.axhline(results['t_train'][-1], color='white', ls='--', lw=2)
        ax3.set_xlabel('x', fontsize=12)
        ax3.set_ylabel('t', fontsize=12)
        ax3.set_title('Absolute Error', fontsize=14, fontweight='bold')
        plt.colorbar(im3, ax=ax3)
        
        # 4. Spatial error
        ax4 = fig.add_subplot(gs[1, 0])
        spatial_error = pde_metrics['spatial_error_mean']
        ax4.plot(x, spatial_error, 'b-', lw=2)
        ax4.fill_between(x, 0, spatial_error, alpha=0.3)
        ax4.set_xlabel('x', fontsize=12)
        ax4.set_ylabel('Mean |Error|', fontsize=12)
        ax4.set_title('Spatial Error Distribution', fontsize=14, fontweight='bold')
        ax4.grid(alpha=0.3)
        
        # 5. Temporal error
        ax5 = fig.add_subplot(gs[1, 1])
        temporal_error = pde_metrics['temporal_error']
        ax5.semilogy(t, temporal_error + 1e-10, 'b-', lw=2)
        ax5.axvline(results['t_train'][-1], color='red', ls='--', lw=2)
        ax5.set_xlabel('t', fontsize=12)
        ax5.set_ylabel('RMSE', fontsize=12)
        ax5.set_title('Temporal Error Evolution', fontsize=14, fontweight='bold')
        ax5.grid(alpha=0.3)
        
        # 6. Conservation error
        ax6 = fig.add_subplot(gs[1, 2])
        conservation_error = pde_metrics['conservation_error']
        ax6.plot(t, conservation_error, 'g-', lw=2)
        ax6.axvline(results['t_train'][-1], color='red', ls='--', lw=2)
        ax6.set_xlabel('t', fontsize=12)
        ax6.set_ylabel('Conservation Error', fontsize=12)
        ax6.set_title('Conservation Error', fontsize=14, fontweight='bold')
        ax6.grid(alpha=0.3)
        
        plt.suptitle(f'PDE Learning Metrics: {problem_name}', 
                    fontsize=16, fontweight='bold', y=0.98)
        
        save_figure(fig, filename, self.config)
        
        return fig
    
    def plot_optimization_analysis(self, opt_results: Dict,
                                  filename: str = "optimization_analysis"):
        """Plot Bayesian optimization analysis."""
        optimizer = opt_results['optimizer']
        
        fig = plt.figure(figsize=(14 * self.figsize_scale, 8 * self.figsize_scale))
        gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)
        
        iterations = range(1, len(optimizer.y_observed) + 1)
        
        # 1. Convergence
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(iterations, optimizer.y_observed, 'o-', alpha=0.6, label='Observed')
        ax1.axhline(optimizer.best_y, color='red', ls='--', lw=2, label=f'Best: {optimizer.best_y:.4e}')
        ax1.axvline(optimizer.n_initial, color='gray', ls=':', lw=2)
        ax1.set_xlabel('Iteration', fontsize=12)
        ax1.set_ylabel('Objective', fontsize=12)
        ax1.set_title('Convergence', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(alpha=0.3)
        
        # 2. Best so far
        ax2 = fig.add_subplot(gs[0, 1])
        best_so_far = np.minimum.accumulate(optimizer.y_observed)
        ax2.plot(iterations, best_so_far, 'g-', lw=2)
        ax2.set_xlabel('Iteration', fontsize=12)
        ax2.set_ylabel('Best Objective', fontsize=12)
        ax2.set_title('Best Value Progress', fontsize=14, fontweight='bold')
        ax2.grid(alpha=0.3)
        
        # 3. Parameter exploration
        ax3 = fig.add_subplot(gs[0, 2])
        X_obs = np.array(optimizer.X_observed)
        
        for i in range(min(X_obs.shape[1], 5)):
            normalized = (X_obs[:, i] - X_obs[:, i].min()) / (X_obs[:, i].max() - X_obs[:, i].min() + 1e-10)
            ax3.plot(iterations, normalized, '-', alpha=0.7, label=optimizer.variable_names[i][:10])
        ax3.set_xlabel('Iteration', fontsize=12)
        ax3.set_ylabel('Normalized Value', fontsize=12)
        ax3.set_title('Parameter Exploration', fontsize=14, fontweight='bold')
        ax3.legend(fontsize=8)
        ax3.grid(alpha=0.3)
        
        # 4. Optimal parameters
        ax4 = fig.add_subplot(gs[1, 0])
        best_params = optimizer.best_x
        param_names = optimizer.variable_names
        
        normalized_best = (best_params - optimizer.bounds_lower) / (optimizer.bounds_upper - optimizer.bounds_lower + 1e-10)
        
        colors_opt = plt.cm.coolwarm(normalized_best)
        ax4.bar(range(len(param_names)), normalized_best, color=colors_opt)
        ax4.set_xticks(range(len(param_names)))
        ax4.set_xticklabels(param_names, rotation=45, ha='right', fontsize=8)
        ax4.set_ylabel('Normalized Value', fontsize=10)
        ax4.set_title('Optimal Parameters', fontsize=14, fontweight='bold')
        ax4.set_ylim([0, 1])
        ax4.grid(alpha=0.3, axis='y')
        
        # 5. Distribution
        ax5 = fig.add_subplot(gs[1, 1])
        ax5.hist(optimizer.y_observed, bins=20, density=True, alpha=0.7, color='steelblue', edgecolor='black')
        ax5.axvline(optimizer.best_y, color='red', ls='--', lw=2)
        ax5.set_xlabel('Objective Value', fontsize=12)
        ax5.set_ylabel('Density', fontsize=12)
        ax5.set_title('Objective Distribution', fontsize=14, fontweight='bold')
        ax5.grid(alpha=0.3)
        
        # 6. Summary
        ax6 = fig.add_subplot(gs[1, 2])
        ax6.axis('off')
        
        summary_text = (
            f"OPTIMIZATION SUMMARY\n"
            f"{'='*35}\n\n"
            f"Iterations: {len(optimizer.y_observed)}\n"
            f"Initial samples: {optimizer.n_initial}\n"
            f"Parameters: {len(optimizer.variable_names)}\n\n"
            f"Best objective: {optimizer.best_y:.6e}\n"
            f"Median: {np.median(optimizer.y_observed):.6e}\n\n"
            f"TOP PARAMETERS:\n"
        )
        
        for name, val in zip(optimizer.variable_names[:6], optimizer.best_x[:6]):
            summary_text += f"  {name}: {val:.4f}\n"
        
        ax6.text(0.05, 0.95, summary_text, transform=ax6.transAxes,
                fontsize=10, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.suptitle('Bayesian Optimization Analysis', fontsize=16, fontweight='bold', y=0.98)
        
        save_figure(fig, filename, self.config)
        
        return fig
    
    def plot_horizon_comparison(self, horizon_results: Dict,
                               filename: str = "horizon_comparison"):
        """Plot comparison across different train-test horizons."""
        fig = plt.figure(figsize=(14 * self.figsize_scale, 10 * self.figsize_scale))
        gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
        
        horizons = horizon_results['horizons']
        train_rmse = horizon_results['train_rmse']
        test_rmse = horizon_results['test_rmse']
        train_r2 = horizon_results['train_r2']
        test_r2 = horizon_results['test_r2']
        
        # 1. RMSE vs horizon
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.plot(horizons, train_rmse, 'o-', color='steelblue', lw=2, label='Train')
        ax1.plot(horizons, test_rmse, 's-', color='darkorange', lw=2, label='Test')
        ax1.set_xlabel('Test Horizon Factor', fontsize=12)
        ax1.set_ylabel('RMSE', fontsize=12)
        ax1.set_title('RMSE vs Test Horizon', fontsize=14, fontweight='bold')
        ax1.legend()
        ax1.grid(alpha=0.3)
        
        # 2. R² vs horizon
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.plot(horizons, train_r2, 'o-', color='steelblue', lw=2, label='Train')
        ax2.plot(horizons, test_r2, 's-', color='darkorange', lw=2, label='Test')
        ax2.set_xlabel('Test Horizon Factor', fontsize=12)
        ax2.set_ylabel('R²', fontsize=12)
        ax2.set_title('R² vs Test Horizon', fontsize=14, fontweight='bold')
        ax2.set_ylim([0, 1.1])
        ax2.legend()
        ax2.grid(alpha=0.3)
        
        # 3. Error growth
        ax3 = fig.add_subplot(gs[1, 0])
        ax3.semilogy(horizons, test_rmse, 's-', color='darkorange', lw=2)
        ax3.set_xlabel('Test Horizon Factor', fontsize=12)
        ax3.set_ylabel('Test RMSE (log)', fontsize=12)
        ax3.set_title('Error Growth with Horizon', fontsize=14, fontweight='bold')
        ax3.grid(alpha=0.3)
        
        # 4. Summary
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.axis('off')
        
        best_idx = np.argmin(test_rmse)
        
        summary_text = (
            f"HORIZON ANALYSIS SUMMARY\n"
            f"{'='*40}\n\n"
            f"Horizons tested: {len(horizons)}\n"
            f"Range: [{min(horizons):.2f}, {max(horizons):.2f}]\n\n"
            f"Best Test RMSE:\n"
            f"  Horizon: {horizons[best_idx]:.2f}\n"
            f"  RMSE: {test_rmse[best_idx]:.4e}\n"
            f"  R²: {test_r2[best_idx]:.4f}\n\n"
            f"Worst Test RMSE:\n"
            f"  Horizon: {horizons[np.argmax(test_rmse)]:.2f}\n"
            f"  RMSE: {max(test_rmse):.4e}\n"
        )
        
        ax4.text(0.1, 0.9, summary_text, transform=ax4.transAxes,
                fontsize=11, verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        plt.suptitle('Train-Test Horizon Comparison', fontsize=16, fontweight='bold', y=0.98)
        
        save_figure(fig, filename, self.config)
        
        return fig