"""
Sensitivity analysis for chemical reservoir configurations.
"""

import numpy as np
from typing import Dict, List, Optional
from itertools import combinations
import warnings

from .config import get_config


class SensitivityAnalyzer:
    """
    Sensitivity analysis for chemical reservoir configurations.
    
    Analyzes how perturbations affect solution quality.
    """
    
    def __init__(self, reservoir, problem):
        """
        Initialize analyzer.
        
        Args:
            reservoir: ConstrainedChemicalReservoir instance
            problem: ODEProblem or PDEProblem1D instance
        """
        self.reservoir = reservoir
        self.problem = problem
        self.sensitivity_results = {}
    
    def analyze_concentration_sensitivity(self,
                                         perturbation_fractions: List[float] = None,
                                         n_samples: int = None,
                                         n_train: int = None) -> Dict:
        """Analyze sensitivity to initial concentration perturbations."""
        scfg = get_config().sensitivity
        if perturbation_fractions is None:
            perturbation_fractions = scfg.perturbation_fractions
        if n_samples is None:
            n_samples = scfg.n_samples
        if n_train is None:
            n_train = get_config().training.n_train
        
        original_conc = self.reservoir.initial_concentrations.copy()
        
        results = {
            'perturbation_fractions': perturbation_fractions,
            'species_names': self.reservoir.species_names.copy(),
            'rmse_vs_perturbation': [],
            'rmse_std_vs_perturbation': [],
            'per_species_sensitivity': {sp: [] for sp in self.reservoir.species_names}
        }
        
        # Global perturbation
        for frac in perturbation_fractions:
            rmse_samples = []
            
            for _ in range(n_samples):
                perturbed = original_conc * (1 + frac * (2 * np.random.rand(len(original_conc)) - 1))
                perturbed = self.reservoir.physical_constraints.clip_concentrations(
                    perturbed, self.reservoir.species_names
                )
                
                self.reservoir.set_initial_concentrations(
                    dict(zip(self.reservoir.species_names, perturbed))
                )
                
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        res = self.reservoir.solve_ode(self.problem, n_train=n_train)
                    rmse_samples.append(res['test_rmse'])
                except:
                    rmse_samples.append(np.nan)
            
            results['rmse_vs_perturbation'].append(np.nanmean(rmse_samples))
            results['rmse_std_vs_perturbation'].append(np.nanstd(rmse_samples))
        
        # Per-species sensitivity
        for sp_idx, species in enumerate(self.reservoir.species_names):
            species_sensitivity = []
            
            for frac in perturbation_fractions:
                rmse_samples = []
                
                for direction in [-1, 1]:
                    perturbed = original_conc.copy()
                    perturbed[sp_idx] *= (1 + direction * frac)
                    perturbed = self.reservoir.physical_constraints.clip_concentrations(
                        perturbed, self.reservoir.species_names
                    )
                    
                    self.reservoir.set_initial_concentrations(
                        dict(zip(self.reservoir.species_names, perturbed))
                    )
                    
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            res = self.reservoir.solve_ode(self.problem, n_train=n_train)
                        rmse_samples.append(res['test_rmse'])
                    except:
                        rmse_samples.append(np.nan)
                
                species_sensitivity.append(np.nanmean(rmse_samples))
            
            results['per_species_sensitivity'][species] = species_sensitivity
        
        # Restore original
        self.reservoir.set_initial_concentrations(
            dict(zip(self.reservoir.species_names, original_conc))
        )
        
        self.sensitivity_results['concentration'] = results
        return results
    
    def analyze_input_rate_sensitivity(self,
                                      perturbation_fractions: List[float] = None,
                                      n_samples: int = None,
                                      n_train: int = None) -> Dict:
        """Analyze sensitivity to input rate perturbations."""
        scfg = get_config().sensitivity
        if perturbation_fractions is None:
            perturbation_fractions = scfg.perturbation_fractions
        if n_samples is None:
            n_samples = scfg.n_samples
        if n_train is None:
            n_train = get_config().training.n_train
        
        original_rates = self.reservoir.input_rates.copy()
        
        results = {
            'perturbation_fractions': perturbation_fractions,
            'input_channels': self.reservoir.input_channels.copy(),
            'rmse_vs_perturbation': [],
            'rmse_std_vs_perturbation': [],
            'per_channel_sensitivity': {ch: [] for ch in self.reservoir.input_channels}
        }
        
        # Global perturbation
        for frac in perturbation_fractions:
            rmse_samples = []
            
            for _ in range(n_samples):
                perturbed_rates = {}
                for ch in self.reservoir.input_channels:
                    base_rate = original_rates.get(ch, scfg.base_val_floor)
                    perturbed = base_rate * (1 + frac * (2 * np.random.rand() - 1))
                    bounds = self.reservoir.physical_constraints.get_bounds(ch)
                    perturbed_rates[ch] = bounds.clip_input_rate(perturbed)
                
                self.reservoir.set_input_rates(perturbed_rates)
                
                try:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        res = self.reservoir.solve_ode(self.problem, n_train=n_train)
                    rmse_samples.append(res['test_rmse'])
                except:
                    rmse_samples.append(np.nan)
            
            results['rmse_vs_perturbation'].append(np.nanmean(rmse_samples))
            results['rmse_std_vs_perturbation'].append(np.nanstd(rmse_samples))
        
        # Per-channel sensitivity
        for channel in self.reservoir.input_channels:
            channel_sensitivity = []
            
            for frac in perturbation_fractions:
                rmse_samples = []
                
                for direction in [-1, 1]:
                    perturbed_rates = original_rates.copy()
                    base_rate = original_rates.get(channel, scfg.base_val_floor)
                    perturbed = base_rate * (1 + direction * frac)
                    bounds = self.reservoir.physical_constraints.get_bounds(channel)
                    perturbed_rates[channel] = bounds.clip_input_rate(perturbed)
                    
                    self.reservoir.set_input_rates(perturbed_rates)
                    
                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            res = self.reservoir.solve_ode(self.problem, n_train=n_train)
                        rmse_samples.append(res['test_rmse'])
                    except:
                        rmse_samples.append(np.nan)
                
                channel_sensitivity.append(np.nanmean(rmse_samples))
            
            results['per_channel_sensitivity'][channel] = channel_sensitivity
        
        # Restore
        self.reservoir.set_input_rates(original_rates)
        
        self.sensitivity_results['input_rate'] = results
        return results
    
    def analyze_measurement_sensitivity(self, n_train: int = None) -> Dict:
        """Analyze how different measurement selections affect performance."""
        if n_train is None:
            n_train = get_config().training.n_train
        scfg = get_config().sensitivity
        original_outputs = self.reservoir.output_channels.copy()
        measurable = self.reservoir.measurement_constraints.get_measurable_set(
            self.reservoir.species_names
        )
        
        n_max = self.reservoir.measurement_constraints.max_simultaneous_measurements
        
        results = {
            'measurement_sets': [],
            'rmse_values': [],
            'r2_values': [],
            'n_measurements': [],
            'best_set': None,
            'worst_set': None
        }
        
        # Test different combinations
        valid_combos = []
        for r in range(1, n_max + 1):
            for combo in combinations(measurable, r):
                combo_set = set(combo)
                valid, _ = self.reservoir.measurement_constraints.is_valid_measurement_set(
                    combo_set, self.reservoir.species_names
                )
                if valid:
                    valid_combos.append(list(combo_set))
                    if len(valid_combos) >= scfg.max_measurement_combos:
                        break
            if len(valid_combos) >= scfg.max_measurement_combos:
                break
        
        for meas_set in valid_combos:
            self.reservoir.output_channels = meas_set
            
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    res = self.reservoir.solve_ode(self.problem, n_train=n_train)
                
                from sklearn.metrics import r2_score
                results['measurement_sets'].append(meas_set)
                results['rmse_values'].append(res['test_rmse'])
                results['r2_values'].append(
                    r2_score(res['y_true'].flatten(), res['y_pred'].flatten())
                )
                results['n_measurements'].append(len(meas_set))
            except:
                pass
        
        if results['rmse_values']:
            best_idx = np.argmin(results['rmse_values'])
            worst_idx = np.argmax(results['rmse_values'])
            results['best_set'] = results['measurement_sets'][best_idx]
            results['worst_set'] = results['measurement_sets'][worst_idx]
            results['best_rmse'] = results['rmse_values'][best_idx]
            results['worst_rmse'] = results['rmse_values'][worst_idx]
        
        # Restore
        self.reservoir.output_channels = original_outputs
        
        self.sensitivity_results['measurement'] = results
        return results
    
    def compute_local_sensitivity(self, epsilon: float = None, n_train: int = None) -> Dict:
        """Compute local sensitivity at optimal configuration."""
        scfg = get_config().sensitivity
        if epsilon is None:
            epsilon = scfg.local_epsilon
        if n_train is None:
            n_train = get_config().training.n_train
        n_inputs = len(self.reservoir.input_channels)
        
        original_rates = self.reservoir.input_rates.copy()
        original_conc = self.reservoir.initial_concentrations.copy()
        
        param_names = (
            [f"rate_{ch}" for ch in self.reservoir.input_channels] +
            [f"conc_{sp}" for sp in self.reservoir.species_names]
        )
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            baseline_res = self.reservoir.solve_ode(self.problem, n_train=n_train)
        baseline_rmse = baseline_res['test_rmse']
        
        sensitivities = {}
        normalized_sensitivities = {}
        
        # Rate sensitivities
        for ch in self.reservoir.input_channels:
            param_name = f"rate_{ch}"
            base_val = original_rates.get(ch, scfg.base_val_floor)
            
            perturbed_rates = original_rates.copy()
            perturbed_rates[ch] = base_val * (1 + epsilon)
            self.reservoir.set_input_rates(perturbed_rates)
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res_plus = self.reservoir.solve_ode(self.problem, n_train=n_train)
            
            perturbed_rates[ch] = base_val * (1 - epsilon)
            self.reservoir.set_input_rates(perturbed_rates)
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res_minus = self.reservoir.solve_ode(self.problem, n_train=n_train)
            
            sensitivity = (res_plus['test_rmse'] - res_minus['test_rmse']) / (2 * epsilon * base_val + scfg.stability_eps)
            sensitivities[param_name] = sensitivity
            normalized_sensitivities[param_name] = (base_val / (baseline_rmse + scfg.stability_eps)) * sensitivity
        
        self.reservoir.set_input_rates(original_rates)
        
        # Concentration sensitivities
        for i, sp in enumerate(self.reservoir.species_names):
            param_name = f"conc_{sp}"
            base_val = max(original_conc[i], scfg.base_val_floor)
            
            perturbed_conc = original_conc.copy()
            perturbed_conc[i] = base_val * (1 + epsilon)
            perturbed_conc = self.reservoir.physical_constraints.clip_concentrations(
                perturbed_conc, self.reservoir.species_names
            )
            self.reservoir.set_initial_concentrations(
                dict(zip(self.reservoir.species_names, perturbed_conc))
            )
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res_plus = self.reservoir.solve_ode(self.problem, n_train=n_train)
            
            perturbed_conc = original_conc.copy()
            perturbed_conc[i] = base_val * (1 - epsilon)
            perturbed_conc = self.reservoir.physical_constraints.clip_concentrations(
                perturbed_conc, self.reservoir.species_names
            )
            self.reservoir.set_initial_concentrations(
                dict(zip(self.reservoir.species_names, perturbed_conc))
            )
            
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res_minus = self.reservoir.solve_ode(self.problem, n_train=n_train)
            
            sensitivity = (res_plus['test_rmse'] - res_minus['test_rmse']) / (2 * epsilon * base_val + scfg.stability_eps)
            sensitivities[param_name] = sensitivity
            normalized_sensitivities[param_name] = (base_val / (baseline_rmse + scfg.stability_eps)) * sensitivity
        
        # Restore
        self.reservoir.set_initial_concentrations(
            dict(zip(self.reservoir.species_names, original_conc))
        )
        
        results = {
            'param_names': param_names,
            'sensitivities': sensitivities,
            'normalized_sensitivities': normalized_sensitivities,
            'baseline_rmse': baseline_rmse
        }
        
        self.sensitivity_results['local'] = results
        return results