"""
Learning metrics for evaluating reservoir computing performance.
"""

import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import r2_score
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import warnings

from .config import get_config


@dataclass
class LearningMetrics:
    """Comprehensive metrics for evaluating reservoir learning performance."""
    
    # Basic error metrics
    train_rmse: float = 0.0
    test_rmse: float = 0.0
    train_mae: float = 0.0
    test_mae: float = 0.0
    train_mape: float = 0.0
    test_mape: float = 0.0
    
    # Correlation metrics
    train_r2: float = 0.0
    test_r2: float = 0.0
    train_pearson: float = 0.0
    test_pearson: float = 0.0
    
    # Normalized errors
    train_nrmse: float = 0.0
    test_nrmse: float = 0.0
    
    # Time-dependent metrics
    error_growth_rate: float = 0.0
    lyapunov_time: float = 0.0
    
    # Per-variable metrics
    per_variable_train_rmse: Dict[str, float] = field(default_factory=dict)
    per_variable_test_rmse: Dict[str, float] = field(default_factory=dict)
    per_variable_r2: Dict[str, float] = field(default_factory=dict)
    
    # Stability metrics
    max_error: float = 0.0
    error_std: float = 0.0
    prediction_stability: float = 0.0
    
    def to_dict(self) -> Dict:
        """Convert to dictionary."""
        return {
            'train_rmse': self.train_rmse,
            'test_rmse': self.test_rmse,
            'train_mae': self.train_mae,
            'test_mae': self.test_mae,
            'train_r2': self.train_r2,
            'test_r2': self.test_r2,
            'train_nrmse': self.train_nrmse,
            'test_nrmse': self.test_nrmse,
            'error_growth_rate': self.error_growth_rate,
            'max_error': self.max_error,
            'prediction_stability': self.prediction_stability,
        }
    
    def print_summary(self):
        """Print formatted summary."""
        print("\n" + "=" * 60)
        print("LEARNING METRICS SUMMARY")
        print("=" * 60)
        
        print("\n[Error Metrics]")
        print(f"  {'Metric':<20} {'Training':>12} {'Testing':>12}")
        print(f"  {'-'*44}")
        print(f"  {'RMSE':<20} {self.train_rmse:>12.6f} {self.test_rmse:>12.6f}")
        print(f"  {'MAE':<20} {self.train_mae:>12.6f} {self.test_mae:>12.6f}")
        print(f"  {'NRMSE':<20} {self.train_nrmse:>12.6f} {self.test_nrmse:>12.6f}")
        
        print("\n[Correlation Metrics]")
        print(f"  {'R²':<20} {self.train_r2:>12.6f} {self.test_r2:>12.6f}")
        print(f"  {'Pearson r':<20} {self.train_pearson:>12.6f} {self.test_pearson:>12.6f}")
        
        print("\n[Stability Metrics]")
        print(f"  Max Error: {self.max_error:.6f}")
        print(f"  Prediction Stability: {self.prediction_stability:.4f}")
        
        if self.per_variable_r2:
            print("\n[Per-Variable R²]")
            for var, r2 in self.per_variable_r2.items():
                print(f"  {var}: {r2:.6f}")
        
        print("=" * 60)


class MetricsCalculator:
    """Calculate comprehensive learning metrics."""
    
    @staticmethod
    def calculate_all(y_true: np.ndarray, 
                     y_pred: np.ndarray,
                     t: np.ndarray,
                     n_train: int,
                     variable_names: List[str] = None) -> LearningMetrics:
        """
        Calculate all metrics from predictions.
        
        Args:
            y_true: Ground truth, shape (n_times, n_vars)
            y_pred: Predictions, shape (n_times, n_vars)
            t: Time points
            n_train: Number of training points
            variable_names: Names of variables
        """
        metrics = LearningMetrics()
        
        y_true = np.atleast_2d(y_true)
        y_pred = np.atleast_2d(y_pred)
        
        if y_true.shape[0] == 1:
            y_true = y_true.T
            y_pred = y_pred.T
        
        n_vars = y_true.shape[1]
        if variable_names is None:
            variable_names = [f"y{i}" for i in range(n_vars)]
        
        # Split train/test
        y_true_train = y_true[:n_train]
        y_true_test = y_true[n_train:]
        y_pred_train = y_pred[:n_train]
        y_pred_test = y_pred[n_train:]
        t_test = t[n_train:]
        
        eps = get_config().metrics.eps
        
        # Basic error metrics
        metrics.train_rmse = np.sqrt(np.mean((y_true_train - y_pred_train)**2))
        metrics.test_rmse = np.sqrt(np.mean((y_true_test - y_pred_test)**2)) if len(y_true_test) > 0 else 0
        
        metrics.train_mae = np.mean(np.abs(y_true_train - y_pred_train))
        metrics.test_mae = np.mean(np.abs(y_true_test - y_pred_test)) if len(y_true_test) > 0 else 0
        
        # MAPE
        metrics.train_mape = np.mean(np.abs((y_true_train - y_pred_train) / (np.abs(y_true_train) + eps)))
        metrics.test_mape = np.mean(np.abs((y_true_test - y_pred_test) / (np.abs(y_true_test) + eps))) if len(y_true_test) > 0 else 0
        
        # Normalized RMSE
        y_range = np.max(y_true) - np.min(y_true)
        if y_range > eps:
            metrics.train_nrmse = metrics.train_rmse / y_range
            metrics.test_nrmse = metrics.test_rmse / y_range
        
        # R² score
        metrics.train_r2 = r2_score(y_true_train.flatten(), y_pred_train.flatten())
        if len(y_true_test) > 0:
            metrics.test_r2 = r2_score(y_true_test.flatten(), y_pred_test.flatten())
        
        # Pearson correlation
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            metrics.train_pearson, _ = pearsonr(y_true_train.flatten(), y_pred_train.flatten())
            if len(y_true_test) > 0:
                metrics.test_pearson, _ = pearsonr(y_true_test.flatten(), y_pred_test.flatten())
        
        # Per-variable metrics
        for i, name in enumerate(variable_names):
            metrics.per_variable_train_rmse[name] = np.sqrt(np.mean((y_true_train[:, i] - y_pred_train[:, i])**2))
            if len(y_true_test) > 0:
                metrics.per_variable_test_rmse[name] = np.sqrt(np.mean((y_true_test[:, i] - y_pred_test[:, i])**2))
                try:
                    metrics.per_variable_r2[name] = r2_score(y_true_test[:, i], y_pred_test[:, i])
                except:
                    metrics.per_variable_r2[name] = 0.0
        
        # Stability metrics
        all_errors = np.abs(y_true - y_pred)
        metrics.max_error = np.max(all_errors)
        metrics.error_std = np.std(all_errors)
        mean_error = np.mean(all_errors)
        if mean_error > eps:
            metrics.prediction_stability = 1.0 - metrics.error_std / mean_error
        
        # Error growth rate
        if len(y_true_test) > 1:
            test_errors = np.sqrt(np.mean((y_true_test - y_pred_test)**2, axis=1))
            if len(test_errors) > 1:
                log_errors = np.log(test_errors + eps)
                t_test_shifted = t_test - t_test[0]
                if len(t_test_shifted) > 1:
                    try:
                        coeffs = np.polyfit(t_test_shifted, log_errors, 1)
                        metrics.error_growth_rate = coeffs[0]
                        if metrics.error_growth_rate > eps:
                            metrics.lyapunov_time = 1.0 / metrics.error_growth_rate
                        else:
                            metrics.lyapunov_time = np.inf
                    except:
                        pass
        
        return metrics
    
    @staticmethod
    def calculate_pde_metrics(sol_true: np.ndarray,
                             sol_pred: np.ndarray,
                             x: np.ndarray,
                             t: np.ndarray,
                             n_train: int) -> Dict:
        """Calculate metrics specific to PDE solutions."""
        metrics = {}
        
        sol_true_train = sol_true[:n_train]
        sol_true_test = sol_true[n_train:]
        sol_pred_train = sol_pred[:n_train]
        sol_pred_test = sol_pred[n_train:]
        
        # Global metrics
        metrics['train_rmse'] = np.sqrt(np.mean((sol_true_train - sol_pred_train)**2))
        metrics['test_rmse'] = np.sqrt(np.mean((sol_true_test - sol_pred_test)**2)) if len(sol_true_test) > 0 else 0
        metrics['train_r2'] = r2_score(sol_true_train.flatten(), sol_pred_train.flatten())
        
        if len(sol_true_test) > 0:
            try:
                metrics['test_r2'] = r2_score(sol_true_test.flatten(), sol_pred_test.flatten())
            except:
                metrics['test_r2'] = 0.0
        
        # Spatial error distribution
        spatial_error = np.mean(np.abs(sol_true - sol_pred), axis=0)
        metrics['spatial_error_mean'] = spatial_error
        metrics['spatial_error_max_location'] = x[np.argmax(spatial_error)]
        
        # Temporal error evolution
        temporal_error = np.sqrt(np.mean((sol_true - sol_pred)**2, axis=1))
        metrics['temporal_error'] = temporal_error
        
        # Conservation error
        _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
        true_integral = _trapz(sol_true, x, axis=1)
        pred_integral = _trapz(sol_pred, x, axis=1)
        metrics['conservation_error'] = np.abs(true_integral - pred_integral)
        
        return metrics