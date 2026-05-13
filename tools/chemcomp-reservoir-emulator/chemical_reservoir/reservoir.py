"""
Main chemical reservoir computing class.

Provides ConstrainedChemicalReservoir with two prediction modes:
  - Open-loop (solve_ode): the reservoir is driven by a predetermined input
    encoding over the full time domain.  Predictions are made in a single
    forward pass with no feedback.
  - Closed-loop (solve_ode_closed_loop): training uses teacher forcing with
    the same encoding; testing feeds the predicted output back to modulate
    the input signal at each step, enabling autonomous attractor tracking.
"""""

import numpy as np
from scipy.integrate import solve_ivp
from sklearn.linear_model import RidgeCV, Ridge
from sklearn.preprocessing import StandardScaler
from typing import List, Dict, Callable, Optional, Tuple, Union
import warnings

from .config import get_config
from .constraints import PhysicalConstraints, MeasurementConstraints
from .optimization import BayesianOptimizer
from .network_viz import ReactionNetworkVisualizer
from .config import get_output_config


class MeasurementSelector:
    """Selects optimal measurement channels."""
    
    def __init__(self, constraints: MeasurementConstraints):
        self.constraints = constraints
    
    def select_optimal(self, trajectories: np.ndarray,
                      target: np.ndarray,
                      species_names: List[str],
                      method: str = 'variance') -> List[str]:
        from itertools import combinations
        
        measurable = self.constraints.get_measurable_set(species_names)
        measurable_list = [s for s in species_names if s in measurable]
        
        n_select = min(self.constraints.max_simultaneous_measurements, len(measurable_list))
        selected = set(self.constraints.required_species & measurable)
        
        if len(selected) >= n_select:
            return list(selected)[:n_select]
        
        candidates = [s for s in measurable_list if s not in selected]
        cand_idx = [species_names.index(s) for s in candidates]
        
        if len(cand_idx) == 0:
            return list(selected)
        
        if method == 'variance':
            scores = np.var(trajectories[:, cand_idx], axis=0)
        else:
            scores = np.ones(len(candidates))
        
        sorted_cands = [c for _, c in sorted(zip(scores, candidates), reverse=True)]
        
        for cand in sorted_cands:
            if len(selected) >= n_select:
                break
            test_set = selected | {cand}
            valid, _ = self.constraints.is_valid_measurement_set(test_set, species_names)
            if valid:
                selected.add(cand)
        
        return list(selected)


class ConstrainedChemicalReservoir:
    """
    Chemical Reservoir Computer with constraints and Bayesian optimization.

    Supports two prediction modes:
      - **Open-loop** (`solve_ode`): input encoding drives the reservoir over
        the entire time domain; the readout predicts in one pass.
      - **Closed-loop** (`solve_ode_closed_loop`): training uses teacher
        forcing; during testing, predicted outputs are fed back as additive
        modulations on the encoding signal, enabling the reservoir to
        self-correct and track attractors autonomously.

    PDE problems are solved via `solve_pde` (open-loop only).
    """
    
    def __init__(self,
                 reactions: List,
                 input_channels: List[str],
                 output_channels: List[str],
                 physical_constraints: Optional[PhysicalConstraints] = None,
                 measurement_constraints: Optional[MeasurementConstraints] = None,
                 ridge_alpha: Optional[float] = None,
                 feature_expansion: Optional[str] = None):
        
        cfg = get_config()
        if ridge_alpha is None:
            ridge_alpha = cfg.reservoir.ridge_alpha
        if feature_expansion is None:
            feature_expansion = cfg.reservoir.feature_expansion

        self._parse_reactions(reactions)
        
        self.input_channels = input_channels
        self.output_channels = output_channels
        
        self.physical_constraints = physical_constraints or PhysicalConstraints()
        self.measurement_constraints = measurement_constraints or MeasurementConstraints()
        
        for sp in self.species_names:
            _ = self.physical_constraints.get_bounds(sp)
        
        self.ridge_alpha = ridge_alpha
        self.feature_expansion = feature_expansion
        
        self.scaler = StandardScaler()
        self.model = None
        
        self.initial_concentrations = None
        self.input_rates = {}
        
        self.measurement_selector = MeasurementSelector(self.measurement_constraints)
        
        self.training_info = {}
        self.optimization_results = None
        self.network_visualizer = None
    
    def _parse_reactions(self, reactions: List):
        """Parse reactions into matrices."""
        species_set = set()
        for rxn in reactions:
            if isinstance(rxn, dict):
                species_set.update(rxn.get('species', []))
            else:
                species_set.update(rxn[3].keys())
                species_set.update(rxn[4].keys())
        
        self.species_names = sorted(list(species_set))
        self.n_species = len(self.species_names)
        self.n_reactions = len(reactions)
        
        self.rate_constants = np.zeros(self.n_reactions)
        self.substrate_matrix = np.zeros((self.n_reactions, self.n_species))
        self.stoich_matrix = np.zeros((self.n_reactions, self.n_species))
        
        for idx, rxn in enumerate(reactions):
            if isinstance(rxn, dict):
                self.rate_constants[idx] = rxn['rate']
                for i, sp in enumerate(rxn['species']):
                    sp_idx = self.species_names.index(sp)
                    self.substrate_matrix[idx, sp_idx] = rxn['reactants'][i]
                    self.stoich_matrix[idx, sp_idx] = rxn['products'][i] - rxn['reactants'][i]
            else:
                self.rate_constants[idx] = rxn[1]
                for sp, stoich in rxn[3].items():
                    sp_idx = self.species_names.index(sp)
                    self.substrate_matrix[idx, sp_idx] = stoich
                for sp, stoich in rxn[4].items():
                    sp_idx = self.species_names.index(sp)
                    self.stoich_matrix[idx, sp_idx] = stoich
    
    def _reaction_rates(self, concentrations: np.ndarray) -> np.ndarray:
        """Compute mass-action reaction rates."""
        rates = self.rate_constants.copy()
        for rxn_idx in range(self.n_reactions):
            for sp_idx in range(self.n_species):
                stoich = self.substrate_matrix[rxn_idx, sp_idx]
                if stoich > 0:
                    rates[rxn_idx] *= np.power(max(concentrations[sp_idx], 0), stoich)
        return rates
    
    def _dynamics(self, t: float, y: np.ndarray,
                  input_func: Optional[Callable] = None) -> np.ndarray:
        """Chemical system dynamics."""
        y = self.physical_constraints.clip_concentrations(y, self.species_names)
        rates = self._reaction_rates(y)
        dydt = self.stoich_matrix.T @ rates
        
        if input_func is not None:
            inputs = input_func(t)
            for sp, rate in inputs.items():
                if sp in self.species_names:
                    dydt[self.species_names.index(sp)] += rate
        
        return dydt
    
    def _simulate(self, t_eval: np.ndarray,
                  input_func: Optional[Callable] = None,
                  y0: Optional[np.ndarray] = None) -> np.ndarray:
        """Simulate reservoir dynamics."""
        if y0 is None:
            if self.initial_concentrations is None:
                y0 = np.array([
                    self.physical_constraints.get_bounds(sp).initial_guess_concentration
                    for sp in self.species_names
                ])
            else:
                y0 = self.initial_concentrations.copy()
        
        y0 = self.physical_constraints.clip_concentrations(y0, self.species_names)
        
        t_eval = np.asarray(t_eval)
        t_eval_unique = np.unique(t_eval)
        t_span = (t_eval_unique[0], t_eval_unique[-1])
        
        def rhs(t, y):
            return self._dynamics(t, y, input_func)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solve_ivp(rhs, t_span, y0, t_eval=t_eval_unique, method='BDF')
        
        result = sol.y.T
        for i in range(result.shape[0]):
            result[i] = self.physical_constraints.clip_concentrations(
                result[i], self.species_names
            )
        
        return result
    
    def _create_input_function(self, t_range: np.ndarray,
                               encoding: str = 'mixed') -> Callable:
        """Create time-dependent input function."""
        cfg = get_config().encoding
        t_min, t_max = t_range.min(), t_range.max()
        t_span = t_max - t_min if t_max > t_min else 1.0
        
        def input_func(t):
            inputs = {}
            tau = (t - t_min) / t_span if t_span > 0 else 0.0
            
            for i, ch in enumerate(self.input_channels):
                bounds = self.physical_constraints.get_bounds(ch)
                base_rate = self.input_rates.get(ch, bounds.initial_guess_input_rate)
                
                if encoding == 'constant':
                    rate = base_rate
                elif encoding == 'sinusoidal':
                    mod = cfg.sinusoidal_amplitude * np.sin(2 * np.pi * tau + i * np.pi / max(len(self.input_channels), 1))
                    rate = base_rate * (1 + mod)
                elif encoding == 'mixed':
                    mod = (cfg.mixed_linear_coeff * tau
                           + cfg.mixed_sine_coeff * np.sin(2*np.pi*tau + i)
                           + cfg.mixed_gaussian_coeff * np.exp(-((tau-0.5)**2)/cfg.mixed_gaussian_width))
                    rate = base_rate * (1 + mod)
                elif encoding == 'multiharmonic':
                    n_ch = max(len(self.input_channels), 1)
                    phase = i * np.pi / n_ch
                    mod = sum(
                        cfg.multiharmonic_amplitude / (n + 1) * np.sin(2 * np.pi * (n + 1) * tau + phase * (n + 1))
                        for n in range(cfg.multiharmonic_n_harmonics)
                    )
                    rate = base_rate * (1 + mod)
                else:
                    rate = base_rate
                
                rate = bounds.clip_input_rate(rate)
                inputs[ch] = rate
            
            return inputs
        
        return input_func
    
    def _extract_features(self, trajectories: np.ndarray) -> np.ndarray:
        """Extract features from reservoir state."""
        output_idx = [self.species_names.index(ch) for ch in self.output_channels]
        measured = trajectories[:, output_idx]
        measured = self.measurement_constraints.add_measurement_noise(measured)

        if self.feature_expansion == 'none':
            return measured
        elif self.feature_expansion == 'quadratic':
            parts = [measured, measured**2]
            n = measured.shape[1]
            for i in range(n):
                for j in range(i+1, n):
                    parts.append((measured[:, i] * measured[:, j]).reshape(-1, 1))
            return np.hstack(parts)
        elif self.feature_expansion == 'delayed':
            # Time-delay embedding: current + several delayed copies
            delays = get_config().reservoir.delays
            parts = [measured]
            for d in delays:
                delayed = np.zeros_like(measured)
                delayed[d:] = measured[:-d]
                delayed[:d] = measured[0]  # pad with initial value
                parts.append(delayed)
            return np.hstack(parts)
        return measured
    
    def _evaluate_configuration(self, problem, t_train: np.ndarray,
                               y_true: np.ndarray) -> float:
        """Evaluate configuration for optimization."""
        try:
            input_func = self._create_input_function(t_train, encoding='mixed')
            trajectories = self._simulate(t_train, input_func)
            
            features = self._extract_features(trajectories)
            
            scaler = StandardScaler()
            features_scaled = scaler.fit_transform(features)
            
            model = Ridge(alpha=self.ridge_alpha)
            model.fit(features_scaled, y_true)
            
            y_pred = model.predict(features_scaled)
            mse = np.mean((y_true - y_pred)**2)
            
            return mse
        except:
            return get_config().optimization.failure_penalty
    
    def set_initial_concentrations(self, concentrations: Union[np.ndarray, Dict]):
        """Set initial concentrations with validation."""
        if isinstance(concentrations, dict):
            conc_array = np.array([
                concentrations.get(sp, self.physical_constraints.get_bounds(sp).initial_guess_concentration)
                for sp in self.species_names
            ])
        else:
            conc_array = np.asarray(concentrations)
        
        conc_array = self.physical_constraints.clip_concentrations(conc_array, self.species_names)
        self.initial_concentrations = conc_array
    
    def set_input_rates(self, rates: Dict[str, float]):
        """Set input rates with validation."""
        validated = {}
        for ch, rate in rates.items():
            if ch in self.input_channels:
                bounds = self.physical_constraints.get_bounds(ch)
                validated[ch] = bounds.clip_input_rate(rate)
        self.input_rates.update(validated)
    
    def optimize_with_bayesian(self, problem,
                              n_iterations: Optional[int] = None,
                              n_initial: Optional[int] = None,
                              acquisition: str = 'EI',
                              verbose: bool = True) -> Dict:
        """Optimize using Bayesian optimization."""
        cfg = get_config().optimization
        if n_iterations is None:
            n_iterations = cfg.n_iterations
        if n_initial is None:
            n_initial = cfg.n_initial
        bounds_lower, bounds_upper = self.physical_constraints.get_bounds_arrays(
            self.species_names, self.input_channels
        )
        
        var_names = (
            [f"rate_{ch}" for ch in self.input_channels] +
            [f"conc_{sp}" for sp in self.species_names]
        )
        
        t_train = np.linspace(problem.t_span[0], problem.t_span[1], 300)
        y_true = problem.solve_numerical(t_train)
        
        n_inputs = len(self.input_channels)
        
        def objective(x):
            input_rates = {ch: x[i] for i, ch in enumerate(self.input_channels)}
            init_conc = {sp: x[n_inputs + i] for i, sp in enumerate(self.species_names)}
            
            self.set_input_rates(input_rates)
            self.set_initial_concentrations(init_conc)
            
            return self._evaluate_configuration(problem, t_train, y_true)
        
        optimizer = BayesianOptimizer(
            bounds_lower=bounds_lower,
            bounds_upper=bounds_upper,
            variable_names=var_names,
            acquisition=acquisition,
            n_initial=n_initial
        )
        
        results = optimizer.optimize(objective, n_iterations=n_iterations, verbose=verbose)
        
        best_params = optimizer.get_best_parameters(self.input_channels, self.species_names)
        
        self.set_input_rates(best_params['input_rates'])
        self.set_initial_concentrations(best_params['concentrations'])
        
        self.optimization_results = {
            'optimizer': optimizer,
            'best_input_rates': best_params['input_rates'],
            'best_concentrations': best_params['concentrations'],
            'best_objective': best_params['objective'],
            'n_iterations': n_iterations,
            'acquisition': acquisition
        }
        
        return self.optimization_results
    
    def solve_ode(self, problem,
                  n_train: Optional[int] = None,
                  n_test: Optional[int] = None,
                  test_horizon_factor: Optional[float] = None,
                  encoding: Optional[str] = None,
                  auto_optimize: bool = False,
                  bo_iterations: Optional[int] = None) -> Dict:
        """Solve ODE in open-loop mode.

        The reservoir is driven by the input encoding over the full time
        domain (train + test).  The readout is fitted on training data
        and predictions are made in a single forward pass with no feedback.

        For closed-loop prediction with output feedback, use
        ``solve_ode_closed_loop`` instead.
        """
        cfg = get_config()
        if n_train is None:
            n_train = cfg.training.n_train
        if n_test is None:
            n_test = cfg.training.n_test
        if test_horizon_factor is None:
            test_horizon_factor = cfg.training.test_horizon_factor
        if encoding is None:
            encoding = cfg.encoding.method
        if bo_iterations is None:
            bo_iterations = cfg.optimization.n_iterations

        if auto_optimize:
            print("Running Bayesian optimization...")
            self.optimize_with_bayesian(problem, n_iterations=bo_iterations)
        
        t_train, t_test, t_all = problem.get_train_test_times(
            n_train, n_test, test_horizon_factor
        )
        
        y_true_train = problem.solve_numerical(t_train)
        y_true_all = problem.solve_numerical(t_all)
        
        input_func = self._create_input_function(t_all, encoding=encoding)
        trajectories = self._simulate(t_all, input_func)
        
        features = self._extract_features(trajectories)
        features_train = self.scaler.fit_transform(features[:n_train])
        features_test = self.scaler.transform(features[n_train:])
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rcfg = cfg.reservoir
            self.model = RidgeCV(alphas=np.logspace(
                rcfg.ridge_cv_alphas_log_min, rcfg.ridge_cv_alphas_log_max, rcfg.ridge_cv_alphas_n
            ))
            self.model.fit(features_train, y_true_train)
        
        y_pred_train = self.model.predict(features_train)
        y_pred_test = self.model.predict(features_test)
        y_pred_all = np.vstack([y_pred_train, y_pred_test])
        
        train_rmse = np.sqrt(np.mean((y_true_train - y_pred_train)**2))
        test_rmse = np.sqrt(np.mean((y_true_all[n_train:] - y_pred_test)**2))
        
        self.training_info = {
            'train_rmse': train_rmse,
            'test_rmse': test_rmse,
            'ridge_alpha': self.model.alpha_,
            'n_features': features.shape[1],
            'output_channels': self.output_channels.copy(),
            'test_horizon_factor': test_horizon_factor
        }
        
        return {
            't_train': t_train,
            't_test': t_test,
            't_all': t_all,
            'y_true': y_true_all,
            'y_pred': y_pred_all,
            'train_rmse': train_rmse,
            'test_rmse': test_rmse,
            'trajectories': trajectories,
            'problem': problem,
            'test_horizon_factor': test_horizon_factor
        }
    
    def solve_ode_closed_loop(self, problem,
                              n_train: Optional[int] = None,
                              n_test: Optional[int] = None,
                              test_horizon_factor: Optional[float] = None,
                              encoding: Optional[str] = None,
                              feedback_scale: Optional[float] = None,
                              auto_optimize: bool = False,
                              bo_iterations: Optional[int] = None) -> Dict:
        """Solve ODE with closed-loop (autonomous) prediction.
        
        Training phase uses teacher forcing: the reservoir is driven by the
        standard input encoding and the readout is fitted to the true solution.
        
        Test phase uses output feedback: at each step the reservoir's own
        predicted output is fed back to modulate input rates, enabling
        autonomous tracking of the attractor without access to the true
        solution.
        
        Args:
            problem: ODE problem instance
            n_train: Number of training time points
            n_test: Number of test time points
            test_horizon_factor: Fraction of training span for test horizon
            encoding: Input encoding for training phase
            feedback_scale: Strength of output-to-input feedback (0 to 1)
            auto_optimize: Whether to run Bayesian optimization first
            bo_iterations: Number of BO iterations
        """
        cfg = get_config()
        if n_train is None:
            n_train = cfg.training.n_train
        if n_test is None:
            n_test = cfg.training.n_test
        if test_horizon_factor is None:
            test_horizon_factor = cfg.training.test_horizon_factor
        if encoding is None:
            encoding = cfg.encoding.method
        if feedback_scale is None:
            feedback_scale = cfg.closed_loop.feedback_scale
        if bo_iterations is None:
            bo_iterations = cfg.optimization.n_iterations

        if auto_optimize:
            print("Running Bayesian optimization...")
            self.optimize_with_bayesian(problem, n_iterations=bo_iterations)
        
        t_train, t_test, t_all = problem.get_train_test_times(
            n_train, n_test, test_horizon_factor
        )
        
        y_true_train = problem.solve_numerical(t_train)
        y_true_all = problem.solve_numerical(t_all)
        n_out = y_true_train.shape[1]
        
        # ---- Training phase (teacher forcing) ----
        # Drive reservoir with standard input encoding, fit readout to true solution.
        # The encoding is created over t_all so it extends into the test region.
        input_func = self._create_input_function(t_all, encoding=encoding)
        traj_train = self._simulate(t_train, input_func)
        
        features_train = self._extract_features(traj_train)
        features_train_scaled = self.scaler.fit_transform(features_train)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rcfg = cfg.reservoir
            self.model = RidgeCV(alphas=np.logspace(
                rcfg.ridge_cv_alphas_log_min, rcfg.ridge_cv_alphas_log_max, rcfg.ridge_cv_alphas_n
            ))
            self.model.fit(features_train_scaled, y_true_train)
        
        y_pred_train = self.model.predict(features_train_scaled)
        
        # Normalization statistics for the feedback signal so that the
        # feedback modulation is O(1) regardless of problem scale
        y_train_mean = y_true_train.mean(axis=0)
        y_train_std = y_true_train.std(axis=0) + cfg.closed_loop.normalization_eps
        
        # ---- Test phase (closed-loop with output feedback) ----
        # Step through test time points one at a time.
        # The input is the SAME encoding used in training, plus a feedback
        # correction from the predicted output.  This keeps the reservoir
        # in a consistent input regime while the feedback steers it.
        max_delay = max(cfg.reservoir.delays) if self.feature_expansion == 'delayed' else 0
        buffer_size = max_delay + 1
        
        # Seed the state buffer with the tail of the training trajectory
        # so that time-delay features are available from the first test step
        n_buf_init = min(buffer_size, len(traj_train))
        state_buffer = list(traj_train[-n_buf_init:])
        
        y_pred_test = np.zeros((n_test, n_out))
        current_state = traj_train[-1].copy()
        last_pred = y_pred_train[-1].copy()
        
        for k in range(n_test):
            t_start = t_train[-1] if k == 0 else t_test[k - 1]
            t_end = t_test[k]
            
            # Normalize the last prediction for feedback modulation
            pred_norm = (last_pred - y_train_mean) / y_train_std
            
            # Combined input: original time-dependent encoding + feedback
            # correction.  The encoding provides the base driving signal
            # (consistent with training), the feedback nudges the reservoir.
            captured_pn = pred_norm.copy()
            def combined_func(t, _pn=captured_pn, _enc=input_func):
                base_inputs = _enc(t)
                for i, ch in enumerate(self.input_channels):
                    bounds = self.physical_constraints.get_bounds(ch)
                    j = i % len(_pn)
                    fb = feedback_scale * float(_pn[j])
                    base_inputs[ch] *= (1.0 + fb)
                    base_inputs[ch] = bounds.clip_input_rate(base_inputs[ch])
                return base_inputs
            
            # Integrate reservoir one time step with combined input
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sol = solve_ivp(
                    lambda t, y: self._dynamics(t, y, combined_func),
                    (t_start, t_end),
                    current_state,
                    method='BDF',
                    t_eval=[t_end],
                    max_step=(t_end - t_start),
                )
            
            if sol.success and sol.y.shape[1] > 0:
                current_state = sol.y[:, -1].copy()
            current_state = self.physical_constraints.clip_concentrations(
                current_state, self.species_names
            )
            
            # Append to buffer and extract features from the sliding window
            state_buffer.append(current_state.copy())
            buf_arr = np.array(state_buffer[-buffer_size:])
            features_step = self._extract_features(buf_arr)
            feature_vec = self.scaler.transform(features_step[-1:])
            
            # Predict and feed back
            pred = self.model.predict(feature_vec)[0]
            y_pred_test[k] = pred
            last_pred = pred.copy()
        
        y_pred_all = np.vstack([y_pred_train, y_pred_test])
        
        # Build full trajectory array for plotting
        traj_test = np.array(state_buffer[-n_test:])
        trajectories_all = np.vstack([traj_train, traj_test])
        
        train_rmse = np.sqrt(np.mean((y_true_train - y_pred_train) ** 2))
        test_rmse = np.sqrt(np.mean((y_true_all[n_train:] - y_pred_test) ** 2))
        
        self.training_info = {
            'train_rmse': train_rmse,
            'test_rmse': test_rmse,
            'ridge_alpha': self.model.alpha_,
            'n_features': features_train.shape[1],
            'output_channels': self.output_channels.copy(),
            'test_horizon_factor': test_horizon_factor,
            'mode': 'closed-loop',
            'feedback_scale': feedback_scale,
        }
        
        return {
            't_train': t_train,
            't_test': t_test,
            't_all': t_all,
            'y_true': y_true_all,
            'y_pred': y_pred_all,
            'train_rmse': train_rmse,
            'test_rmse': test_rmse,
            'trajectories': trajectories_all,
            'problem': problem,
            'test_horizon_factor': test_horizon_factor,
            'mode': 'closed-loop',
        }
    
    def solve_pde(self, problem,
                  n_train_time: Optional[int] = None,
                  n_test_time: Optional[int] = None,
                  test_horizon_factor: Optional[float] = None) -> Dict:
        """Solve PDE with configurable horizon."""
        cfg = get_config().training
        if n_train_time is None:
            n_train_time = cfg.n_train_time
        if n_test_time is None:
            n_test_time = cfg.n_test_time
        if test_horizon_factor is None:
            test_horizon_factor = cfg.pde_test_horizon_factor

        t_train, t_test, t_all = problem.get_train_test_times(
            n_train_time, n_test_time, test_horizon_factor
        )
        
        x_grid, sol_train = problem.solve_numerical(t_train)
        _, sol_all = problem.solve_numerical(t_all)
        
        input_func = self._create_input_function(t_all)
        trajectories = self._simulate(t_all, input_func)
        
        features = self._extract_features(trajectories)
        features_train = self.scaler.fit_transform(features[:n_train_time])
        features_all = self.scaler.transform(features)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.model = RidgeCV(alphas=np.logspace(-6, 2, 50))
            self.model.fit(features_train, sol_train)
        
        sol_pred = self.model.predict(features_all)
        
        train_rmse = np.sqrt(np.mean((sol_train - sol_pred[:n_train_time])**2))
        test_rmse = np.sqrt(np.mean((sol_all[n_train_time:] - sol_pred[n_train_time:])**2))
        
        return {
            'x_grid': x_grid,
            't_train': t_train,
            't_test': t_test,
            't_all': t_all,
            'sol_true': sol_all,
            'sol_pred': sol_pred,
            'train_rmse': train_rmse,
            'test_rmse': test_rmse,
            'problem': problem,
            'test_horizon_factor': test_horizon_factor
        }
    
    def plot_reaction_network(self, filename: str = "reaction_network",
                             config=None):
        """Plot reaction network as bipartite graph."""
        if self.network_visualizer is None:
            self.network_visualizer = ReactionNetworkVisualizer(self)
        
        if config is None:
            config = get_output_config()
        
        self.network_visualizer.plot_bipartite_network(filename=filename, config=config)
        self.network_visualizer.plot_detailed_bipartite(filename=f"{filename}_detailed", config=config)
        self.network_visualizer.plot_reaction_equations(filename=f"{filename}_equations", config=config)
    
    def report(self):
        """Print configuration summary."""
        print("\n" + "=" * 70)
        print("CONSTRAINED CHEMICAL RESERVOIR")
        print("=" * 70)
        
        print(f"\nSpecies: {self.species_names}")
        print(f"Reactions: {self.n_reactions}")
        print(f"Input channels: {self.input_channels}")
        print(f"Output channels: {self.output_channels}")
        
        self.physical_constraints.print_bounds_table(
            self.species_names, self.input_channels
        )
        
        if self.training_info:
            print(f"\nTraining Results:")
            print(f"  Train RMSE: {self.training_info.get('train_rmse', 'N/A'):.6f}")
            print(f"  Test RMSE: {self.training_info.get('test_rmse', 'N/A'):.6f}")
        
        print("=" * 70)