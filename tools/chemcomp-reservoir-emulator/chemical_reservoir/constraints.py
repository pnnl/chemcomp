"""
Physical and measurement constraints for chemical reservoir computing.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Set
import warnings


@dataclass
class SpeciesBounds:
    """Bounds for a single chemical species."""
    name: str
    min_concentration: float = 0.0
    max_concentration: float = 10.0
    min_input_rate: float = 0.0
    max_input_rate: float = 5.0
    initial_guess_concentration: Optional[float] = None
    initial_guess_input_rate: Optional[float] = None
    
    def __post_init__(self):
        # Enforce concentration >= 0
        if self.min_concentration < 0:
            self.min_concentration = 0.0
        
        if self.max_concentration < self.min_concentration:
            raise ValueError(f"{self.name}: max_concentration < min_concentration")
        if self.max_input_rate < self.min_input_rate:
            raise ValueError(f"{self.name}: max_input_rate < min_input_rate")
        
        if self.initial_guess_concentration is None:
            self.initial_guess_concentration = 0.5 * (self.min_concentration + self.max_concentration)
        if self.initial_guess_input_rate is None:
            self.initial_guess_input_rate = 0.5 * (self.min_input_rate + self.max_input_rate)
    
    def clip_concentration(self, value: float) -> float:
        """Clip concentration to valid bounds."""
        return np.clip(value, self.min_concentration, self.max_concentration)
    
    def clip_input_rate(self, value: float) -> float:
        """Clip input rate to valid bounds."""
        return np.clip(value, self.min_input_rate, self.max_input_rate)
    
    def sample_concentration(self, rng: np.random.Generator = None) -> float:
        """Sample random concentration within bounds."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.uniform(self.min_concentration, self.max_concentration)
    
    def sample_input_rate(self, rng: np.random.Generator = None) -> float:
        """Sample random input rate within bounds."""
        if rng is None:
            rng = np.random.default_rng()
        return rng.uniform(self.min_input_rate, self.max_input_rate)


@dataclass
class PhysicalConstraints:
    """Physical constraints with per-species bounds."""
    species_bounds: Dict[str, SpeciesBounds] = field(default_factory=dict)
    default_min_concentration: float = 0.0
    default_max_concentration: float = 10.0
    default_min_input_rate: float = 0.0
    default_max_input_rate: float = 5.0
    max_total_concentration: Optional[float] = None
    
    def __post_init__(self):
        if self.default_min_concentration < 0:
            self.default_min_concentration = 0.0
    
    def add_species_bounds(self, species_name: str,
                          min_conc: float = None,
                          max_conc: float = None,
                          min_rate: float = None,
                          max_rate: float = None,
                          init_conc: float = None,
                          init_rate: float = None):
        """Add or update bounds for a specific species."""
        self.species_bounds[species_name] = SpeciesBounds(
            name=species_name,
            min_concentration=min_conc if min_conc is not None else self.default_min_concentration,
            max_concentration=max_conc if max_conc is not None else self.default_max_concentration,
            min_input_rate=min_rate if min_rate is not None else self.default_min_input_rate,
            max_input_rate=max_rate if max_rate is not None else self.default_max_input_rate,
            initial_guess_concentration=init_conc,
            initial_guess_input_rate=init_rate
        )
    
    def get_bounds(self, species_name: str) -> SpeciesBounds:
        """Get bounds for a species (creates default if not exists)."""
        if species_name not in self.species_bounds:
            self.species_bounds[species_name] = SpeciesBounds(
                name=species_name,
                min_concentration=self.default_min_concentration,
                max_concentration=self.default_max_concentration,
                min_input_rate=self.default_min_input_rate,
                max_input_rate=self.default_max_input_rate
            )
        return self.species_bounds[species_name]
    
    def get_concentration_bounds(self, species: str) -> Tuple[float, float]:
        """Get (min, max) concentration for a species."""
        bounds = self.get_bounds(species)
        return (bounds.min_concentration, bounds.max_concentration)
    
    def get_input_rate_bounds(self, species: str) -> Tuple[float, float]:
        """Get (min, max) input rate for a species."""
        bounds = self.get_bounds(species)
        return (bounds.min_input_rate, bounds.max_input_rate)
    
    def validate_concentrations(self, concentrations: np.ndarray,
                               species_names: List[str]) -> Tuple[bool, str]:
        """Check if concentrations satisfy all constraints."""
        for conc, species in zip(concentrations, species_names):
            lb, ub = self.get_concentration_bounds(species)
            if conc < lb - 1e-10 or conc > ub + 1e-10:
                return False, f"{species}: {conc:.6f} outside [{lb}, {ub}]"
        
        if self.max_total_concentration is not None:
            total = np.sum(concentrations)
            if total > self.max_total_concentration:
                return False, f"Total {total:.4f} > max {self.max_total_concentration}"
        
        return True, "OK"
    
    def clip_concentrations(self, concentrations: np.ndarray,
                           species_names: List[str]) -> np.ndarray:
        """Clip concentrations to valid bounds."""
        clipped = concentrations.copy()
        for i, species in enumerate(species_names):
            bounds = self.get_bounds(species)
            clipped[i] = bounds.clip_concentration(clipped[i])
        return clipped
    
    def get_bounds_arrays(self, species_names: List[str],
                         input_channels: List[str]) -> Tuple[np.ndarray, np.ndarray]:
        """Get bounds as arrays for optimization."""
        n_inputs = len(input_channels)
        n_species = len(species_names)
        
        lower = np.zeros(n_inputs + n_species)
        upper = np.zeros(n_inputs + n_species)
        
        for i, ch in enumerate(input_channels):
            lb, ub = self.get_input_rate_bounds(ch)
            lower[i] = lb
            upper[i] = ub
        
        for i, sp in enumerate(species_names):
            lb, ub = self.get_concentration_bounds(sp)
            lower[n_inputs + i] = lb
            upper[n_inputs + i] = ub
        
        return lower, upper
    
    def print_bounds_table(self, species_names: List[str], input_channels: List[str]):
        """Print formatted table of all bounds."""
        print("\n" + "=" * 80)
        print("SPECIES BOUNDS")
        print("=" * 80)
        print(f"{'Species':<12} {'Conc Min':>10} {'Conc Max':>10} {'Rate Min':>10} {'Rate Max':>10} {'Input?':>8}")
        print("-" * 80)
        
        for sp in species_names:
            bounds = self.get_bounds(sp)
            is_input = "Yes" if sp in input_channels else "No"
            print(f"{sp:<12} {bounds.min_concentration:>10.4f} {bounds.max_concentration:>10.4f} "
                  f"{bounds.min_input_rate:>10.4f} {bounds.max_input_rate:>10.4f} {is_input:>8}")
        
        print("=" * 80)


@dataclass
class MeasurementConstraints:
    """Constraints on which species can be measured."""
    max_simultaneous_measurements: int = 4
    measurable_species: Optional[Set[str]] = None
    unmeasurable_species: Set[str] = field(default_factory=set)
    required_species: Set[str] = field(default_factory=set)
    measurement_noise_std: float = 0.01
    interfering_groups: List[Set[str]] = field(default_factory=list)
    
    def get_measurable_set(self, all_species: List[str]) -> Set[str]:
        """Get the set of species that can be measured."""
        if self.measurable_species is not None:
            candidates = self.measurable_species & set(all_species)
        else:
            candidates = set(all_species)
        return candidates - self.unmeasurable_species
    
    def is_valid_measurement_set(self, selected: Set[str],
                                all_species: List[str]) -> Tuple[bool, str]:
        """Check if a measurement selection is valid."""
        measurable = self.get_measurable_set(all_species)
        
        if not selected.issubset(measurable):
            return False, f"Cannot measure: {selected - measurable}"
        
        if len(selected) > self.max_simultaneous_measurements:
            return False, f"Too many: {len(selected)} > {self.max_simultaneous_measurements}"
        
        missing = self.required_species - selected
        if missing & measurable:
            return False, f"Missing required: {missing & measurable}"
        
        for group in self.interfering_groups:
            if len(selected & group) > 1:
                return False, f"Interfering: {selected & group}"
        
        return True, "OK"
    
    def add_measurement_noise(self, measurements: np.ndarray) -> np.ndarray:
        """Add realistic measurement noise."""
        if self.measurement_noise_std > 0:
            noise = np.random.normal(0, self.measurement_noise_std, measurements.shape)
            return np.maximum(measurements + noise, 0)
        return measurements