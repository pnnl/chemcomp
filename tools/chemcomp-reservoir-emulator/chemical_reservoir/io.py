"""
I/O utilities for loading, saving, and validating chemical reaction networks.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from .constraints import PhysicalConstraints, MeasurementConstraints, SpeciesBounds
from .reservoir import ConstrainedChemicalReservoir


def validate_reaction(rxn: Dict, index: int) -> Dict:
    """
    Validate and normalize a reaction dictionary.

    Parameters
    ----------
    rxn : dict
        Reaction dictionary with keys ``species``, ``reactants``, ``products``, ``rate``.
    index : int
        Reaction index (for error messages).

    Returns
    -------
    dict
        Validated reaction dictionary.
    """
    required_fields = ["species", "reactants", "products", "rate"]
    for field in required_fields:
        if field not in rxn:
            raise ValueError(f"Reaction {index} missing required field: '{field}'")

    if not isinstance(rxn["species"], list) or len(rxn["species"]) < 1:
        raise ValueError(f"Reaction {index}: 'species' must be a non-empty list")

    n_species = len(rxn["species"])

    if not isinstance(rxn["reactants"], list) or len(rxn["reactants"]) != n_species:
        raise ValueError(f"Reaction {index}: 'reactants' must be a list of length {n_species}")

    if not isinstance(rxn["products"], list) or len(rxn["products"]) != n_species:
        raise ValueError(f"Reaction {index}: 'products' must be a list of length {n_species}")

    if not isinstance(rxn["rate"], (int, float)) or rxn["rate"] < 0:
        raise ValueError(f"Reaction {index}: 'rate' must be a non-negative number")

    return {
        "species": [str(s) for s in rxn["species"]],
        "reactants": [int(r) for r in rxn["reactants"]],
        "products": [int(p) for p in rxn["products"]],
        "rate": float(rxn["rate"]),
        "name": rxn.get("name", f"R{index}"),
    }


def load_reaction_network(filepath: Union[str, Path]) -> Dict:
    """
    Load a reaction network configuration from a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Path to the JSON configuration file.

    Returns
    -------
    dict
        Network configuration containing ``reactions``, ``input_channels``,
        ``output_channels``, ``name``, ``all_species``, and optional fields.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Reaction network file not found: {filepath}")

    with open(filepath, "r") as f:
        config = json.load(f)

    if "reactions" not in config:
        raise ValueError("JSON must contain 'reactions' field")

    validated_reactions = []
    all_species = set()
    for i, rxn in enumerate(config["reactions"]):
        validated_rxn = validate_reaction(rxn, i)
        validated_reactions.append(validated_rxn)
        all_species.update(validated_rxn["species"])

    config["reactions"] = validated_reactions
    config["all_species"] = sorted(list(all_species))

    if "name" not in config:
        config["name"] = filepath.stem
    if "input_channels" not in config:
        config["input_channels"] = config["all_species"][:2]
    if "output_channels" not in config:
        config["output_channels"] = config["all_species"][:min(3, len(config["all_species"]))]

    for ch in config["input_channels"]:
        if ch not in all_species:
            raise ValueError(f"Input channel '{ch}' not found in species")
    for ch in config["output_channels"]:
        if ch not in all_species:
            raise ValueError(f"Output channel '{ch}' not found in species")

    return config


def save_reaction_network(
    filepath: Union[str, Path],
    reactions: List[Dict],
    input_channels: List[str],
    output_channels: List[str],
    name: str = "custom_network",
    constraints: Optional[Dict] = None,
    initial_conditions: Optional[Dict] = None,
    input_rates: Optional[Dict] = None,
    description: str = "",
    metadata: Optional[Dict] = None,
):
    """
    Save a reaction network configuration to a JSON file.

    Parameters
    ----------
    filepath : str or Path
        Output file path.
    reactions : list of dict
        Reaction definitions.
    input_channels : list of str
        Input species names.
    output_channels : list of str
        Output species names.
    name : str
        Network name.
    constraints, initial_conditions, input_rates : dict, optional
        Optional configuration sections.
    description : str
        Human-readable description.
    metadata : dict, optional
        Additional metadata.
    """
    config = {
        "name": name,
        "description": description,
        "reactions": reactions,
        "input_channels": input_channels,
        "output_channels": output_channels,
    }
    if constraints is not None:
        config["constraints"] = constraints
    if initial_conditions is not None:
        config["initial_conditions"] = initial_conditions
    if input_rates is not None:
        config["input_rates"] = input_rates
    if metadata is not None:
        config["metadata"] = metadata

    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(config, f, indent=2)


def _create_physical_constraints(constraints_dict: Dict, all_species: List[str]) -> PhysicalConstraints:
    """Create PhysicalConstraints from a config dictionary."""
    physical = PhysicalConstraints(
        default_min_concentration=constraints_dict.get("default_min_concentration", 0.0),
        default_max_concentration=constraints_dict.get("default_max_concentration", 5.0),
        default_min_input_rate=constraints_dict.get("default_min_input_rate", 0.0),
        default_max_input_rate=constraints_dict.get("default_max_input_rate", 2.0),
    )
    if "species_bounds" in constraints_dict:
        for species, bounds in constraints_dict["species_bounds"].items():
            if species in all_species:
                physical.add_species_bounds(
                    species,
                    min_conc=bounds.get("min_conc", 0.0),
                    max_conc=bounds.get("max_conc", 5.0),
                    min_rate=bounds.get("min_rate", 0.0),
                    max_rate=bounds.get("max_rate", 2.0),
                )
    return physical


def _create_measurement_constraints(constraints_dict: Dict, all_species: List[str]) -> MeasurementConstraints:
    """Create MeasurementConstraints from a config dictionary."""
    measurable = set(constraints_dict.get("measurable_species", all_species))
    unmeasurable = set(constraints_dict.get("unmeasurable_species", []))
    return MeasurementConstraints(
        max_simultaneous_measurements=constraints_dict.get("max_simultaneous_measurements", 3),
        measurable_species=measurable,
        unmeasurable_species=unmeasurable,
        measurement_noise_std=constraints_dict.get("measurement_noise_std", 0.01),
    )


def _default_physical_constraints(species: List[str]) -> PhysicalConstraints:
    """Create default physical constraints for a list of species."""
    physical = PhysicalConstraints(
        default_min_concentration=0.0,
        default_max_concentration=5.0,
        default_min_input_rate=0.0,
        default_max_input_rate=2.0,
    )
    for i, sp in enumerate(species):
        physical.add_species_bounds(
            sp,
            min_conc=0.0,
            max_conc=3.0 + (i % 3),
            min_rate=0.0,
            max_rate=1.5 + 0.5 * (i % 2),
        )
    return physical


def _default_measurement_constraints(species: List[str]) -> MeasurementConstraints:
    """Create default measurement constraints for a list of species."""
    n_measurable = min(len(species) - 1, 4) if len(species) > 1 else len(species)
    measurable = set(species[:n_measurable])
    unmeasurable = set(species[n_measurable:]) if len(species) > n_measurable else set()
    return MeasurementConstraints(
        max_simultaneous_measurements=min(3, len(measurable)),
        measurable_species=measurable,
        unmeasurable_species=unmeasurable,
        measurement_noise_std=0.01,
    )


def load_reservoir_from_json(
    filepath: Union[str, Path],
    feature_expansion: str = "quadratic",
    physical_constraints: Optional[PhysicalConstraints] = None,
    measurement_constraints: Optional[MeasurementConstraints] = None,
) -> Tuple[ConstrainedChemicalReservoir, Dict]:
    """
    Create a fully-configured reservoir from a JSON network file.

    Parameters
    ----------
    filepath : str or Path
        Path to JSON network configuration.
    feature_expansion : str
        Feature expansion type (``'none'`` or ``'quadratic'``).
    physical_constraints : PhysicalConstraints, optional
        Override constraints from JSON.
    measurement_constraints : MeasurementConstraints, optional
        Override constraints from JSON.

    Returns
    -------
    reservoir : ConstrainedChemicalReservoir
        Configured reservoir instance.
    config : dict
        Parsed network configuration.
    """
    config = load_reaction_network(filepath)

    if physical_constraints is not None:
        physical = physical_constraints
    elif "constraints" in config and "physical" in config["constraints"]:
        physical = _create_physical_constraints(config["constraints"]["physical"], config["all_species"])
    else:
        physical = _default_physical_constraints(config["all_species"])

    if measurement_constraints is not None:
        measurement = measurement_constraints
    elif "constraints" in config and "measurement" in config["constraints"]:
        measurement = _create_measurement_constraints(config["constraints"]["measurement"], config["all_species"])
    else:
        measurement = _default_measurement_constraints(config["all_species"])

    reservoir = ConstrainedChemicalReservoir(
        reactions=config["reactions"],
        input_channels=config["input_channels"],
        output_channels=config["output_channels"],
        physical_constraints=physical,
        measurement_constraints=measurement,
        feature_expansion=feature_expansion,
    )

    if "initial_conditions" in config:
        reservoir.set_initial_concentrations(config["initial_conditions"])
    else:
        reservoir.set_initial_concentrations({sp: 1.0 for sp in reservoir.species_names})

    if "input_rates" in config:
        reservoir.set_input_rates(config["input_rates"])
    else:
        reservoir.set_input_rates({ch: 0.5 for ch in config["input_channels"]})

    return reservoir, config


def load_all_networks(config_dir: Union[str, Path]) -> Dict[str, Dict]:
    """
    Load all ``*.json`` network files from a directory.

    Parameters
    ----------
    config_dir : str or Path
        Directory containing JSON network files.

    Returns
    -------
    dict
        Mapping of network name → parsed configuration dict.
    """
    config_dir = Path(config_dir)
    if not config_dir.exists():
        raise FileNotFoundError(f"Config directory not found: {config_dir}")

    configs = {}
    for json_file in sorted(config_dir.glob("*.json")):
        try:
            config = load_reaction_network(json_file)
            configs[config["name"]] = config
        except Exception as e:
            import warnings
            warnings.warn(f"Failed to load {json_file.name}: {e}")
    return configs
