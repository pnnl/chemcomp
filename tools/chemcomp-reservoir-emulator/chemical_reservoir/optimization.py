"""
Bayesian optimization for chemical reservoir configuration.
"""

import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize as scipy_minimize
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, RBF, ConstantKernel, WhiteKernel
from abc import ABC, abstractmethod
from typing import Callable, Dict, List, Tuple
import warnings
import time

from .config import get_config


class AcquisitionFunction(ABC):
    """Base class for acquisition functions."""
    
    @abstractmethod
    def __call__(self, X: np.ndarray, gp: GaussianProcessRegressor,
                 y_best: float) -> np.ndarray:
        pass


class ExpectedImprovement(AcquisitionFunction):
    """Expected Improvement acquisition function."""
    
    def __init__(self, xi: float = 0.01):
        self.xi = xi
    
    def __call__(self, X: np.ndarray, gp: GaussianProcessRegressor,
                 y_best: float) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mu, sigma = gp.predict(X, return_std=True)
        sigma = np.maximum(sigma, get_config().optimization.sigma_floor)
        
        imp = y_best - mu - self.xi
        Z = imp / sigma
        ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
        ei[sigma < 1e-9] = 0.0
        
        return ei


class ProbabilityOfImprovement(AcquisitionFunction):
    """Probability of Improvement acquisition function."""
    
    def __init__(self, xi: float = 0.01):
        self.xi = xi
    
    def __call__(self, X: np.ndarray, gp: GaussianProcessRegressor,
                 y_best: float) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mu, sigma = gp.predict(X, return_std=True)
        sigma = np.maximum(sigma, get_config().optimization.sigma_floor)
        
        Z = (y_best - mu - self.xi) / sigma
        return norm.cdf(Z)


class UpperConfidenceBound(AcquisitionFunction):
    """Upper Confidence Bound (Lower for minimization)."""
    
    def __init__(self, kappa: float = 2.0):
        self.kappa = kappa
    
    def __call__(self, X: np.ndarray, gp: GaussianProcessRegressor,
                 y_best: float) -> np.ndarray:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mu, sigma = gp.predict(X, return_std=True)
        return -(mu - self.kappa * sigma)


class BayesianOptimizer:
    """
    Bayesian Optimization for chemical reservoir configuration.
    
    Uses Gaussian Process surrogate model and acquisition function
    to efficiently search the parameter space.
    """
    
    def __init__(self,
                 bounds_lower: np.ndarray,
                 bounds_upper: np.ndarray,
                 variable_names: List[str],
                 acquisition: str = 'EI',
                 kernel: str = 'matern',
                 n_initial: int = None,
                 random_state: int = None):
        """
        Initialize Bayesian optimizer.
        
        Args:
            bounds_lower: Lower bounds for each variable
            bounds_upper: Upper bounds for each variable
            variable_names: Names of variables
            acquisition: 'EI', 'PI', or 'UCB'
            kernel: 'matern' or 'rbf'
            n_initial: Number of initial random samples
            random_state: Random seed
        """
        ocfg = get_config().optimization
        if n_initial is None:
            n_initial = ocfg.n_initial
        if random_state is None:
            random_state = ocfg.random_state

        self.bounds_lower = np.asarray(bounds_lower)
        self.bounds_upper = np.asarray(bounds_upper)
        self.n_dims = len(bounds_lower)
        self.variable_names = variable_names
        self.n_initial = n_initial
        self.rng = np.random.default_rng(random_state)
        
        # Build kernel
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            
            ls_bounds = tuple(ocfg.gp_length_scale_bounds)
            if kernel == 'matern':
                base_kernel = Matern(nu=2.5, length_scale=np.ones(self.n_dims),
                                    length_scale_bounds=ls_bounds)
            else:
                base_kernel = RBF(length_scale=np.ones(self.n_dims),
                                length_scale_bounds=ls_bounds)
            
            self.kernel = ConstantKernel(1.0) * base_kernel + WhiteKernel(
                noise_level=ocfg.gp_noise_level,
                noise_level_bounds=tuple(ocfg.gp_noise_bounds)
            )
        
        # Acquisition function
        if acquisition == 'EI':
            self.acquisition = ExpectedImprovement(xi=ocfg.ei_xi)
        elif acquisition == 'PI':
            self.acquisition = ProbabilityOfImprovement(xi=ocfg.pi_xi)
        elif acquisition == 'UCB':
            self.acquisition = UpperConfidenceBound(kappa=ocfg.ucb_kappa)
        else:
            raise ValueError(f"Unknown acquisition: {acquisition}")
        
        # Storage
        self.X_observed = []
        self.y_observed = []
        self.gp = None
        self.best_x = None
        self.best_y = np.inf
        self.history = []
    
    def _sample_random(self, n: int = 1) -> np.ndarray:
        """Sample random points within bounds."""
        samples = self.rng.uniform(0, 1, size=(n, self.n_dims))
        return samples * (self.bounds_upper - self.bounds_lower) + self.bounds_lower
    
    def _optimize_acquisition(self, n_restarts: int = None,
                             n_samples: int = None) -> np.ndarray:
        ocfg = get_config().optimization
        if n_restarts is None:
            n_restarts = ocfg.acq_n_restarts
        if n_samples is None:
            n_samples = ocfg.acq_n_samples
        """Find point that maximizes acquisition function."""
        X_candidates = self._sample_random(n_samples)
        acq_values = self.acquisition(X_candidates, self.gp, self.best_y)
        best_idx = np.argmax(acq_values)
        x_best = X_candidates[best_idx]
        
        def neg_acquisition(x):
            x = np.clip(x, self.bounds_lower, self.bounds_upper)
            return -self.acquisition(x.reshape(1, -1), self.gp, self.best_y)[0]
        
        best_result = None
        best_acq = -np.inf
        
        start_points = [x_best] + [self._sample_random(1)[0] for _ in range(n_restarts - 1)]
        
        for x0 in start_points:
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = scipy_minimize(
                        neg_acquisition, x0,
                        method='L-BFGS-B',
                        bounds=list(zip(self.bounds_lower, self.bounds_upper))
                    )
                if -result.fun > best_acq:
                    best_acq = -result.fun
                    best_result = result.x
            except Exception:
                continue
        
        if best_result is None:
            return x_best
        
        return np.clip(best_result, self.bounds_lower, self.bounds_upper)
    
    def suggest(self) -> np.ndarray:
        """Suggest next point to evaluate."""
        if len(self.X_observed) < self.n_initial:
            return self._sample_random(1)[0]
        
        X = np.array(self.X_observed)
        y = np.array(self.y_observed)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ocfg = get_config().optimization
            self.gp = GaussianProcessRegressor(
                kernel=self.kernel,
                n_restarts_optimizer=ocfg.gp_n_restarts,
                normalize_y=True,
                random_state=ocfg.random_state
            )
            self.gp.fit(X, y)
        
        return self._optimize_acquisition()
    
    def observe(self, x: np.ndarray, y: float):
        """Record observation."""
        self.X_observed.append(x.copy())
        self.y_observed.append(y)
        
        if y < self.best_y:
            self.best_y = y
            self.best_x = x.copy()
        
        self.history.append({
            'iteration': len(self.X_observed),
            'x': x.copy(),
            'y': y,
            'best_y': self.best_y
        })
    
    def optimize(self, objective_func: Callable,
                n_iterations: int = 50,
                verbose: bool = True) -> Dict:
        """
        Run Bayesian optimization.
        
        Args:
            objective_func: Function to minimize
            n_iterations: Total evaluations
            verbose: Print progress
            
        Returns:
            Dictionary with results
        """
        start_time = time.time()
        
        if verbose:
            print(f"\nBayesian Optimization: {n_iterations} iterations")
            print(f"  Variables: {self.n_dims}")
            print(f"  Initial samples: {self.n_initial}")
            print("-" * 60)
        
        for i in range(n_iterations):
            x = self.suggest()
            y = objective_func(x)
            self.observe(x, y)
            
            if verbose and (i + 1) % max(1, n_iterations // 20) == 0:
                print(f"  Iter {i+1:4d}: y = {y:.6f}, best = {self.best_y:.6f}")
        
        elapsed = time.time() - start_time
        
        if verbose:
            print("-" * 60)
            print(f"Optimization complete in {elapsed:.1f}s")
            print(f"Best objective: {self.best_y:.6f}")
        
        return {
            'best_x': self.best_x,
            'best_y': self.best_y,
            'X_observed': np.array(self.X_observed),
            'y_observed': np.array(self.y_observed),
            'history': self.history,
            'elapsed_time': elapsed,
            'gp': self.gp
        }
    
    def get_best_parameters(self, input_channels: List[str],
                           species_names: List[str]) -> Dict:
        """Extract best parameters as named dictionary."""
        n_inputs = len(input_channels)
        
        input_rates = {ch: self.best_x[i] for i, ch in enumerate(input_channels)}
        concentrations = {sp: self.best_x[n_inputs + i] for i, sp in enumerate(species_names)}
        
        return {
            'input_rates': input_rates,
            'concentrations': concentrations,
            'objective': self.best_y
        }