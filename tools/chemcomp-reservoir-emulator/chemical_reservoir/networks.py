"""
Built-in reaction network definitions.

Each network is a list of reaction dictionaries with keys:
``species``, ``reactants``, ``products``, ``rate``, and optional ``name``.
"""

from typing import Dict, List, Tuple


def get_builtin_networks() -> Dict[str, List[Dict]]:
    """
    Return all built-in reaction network topologies.

    Returns
    -------
    dict
        Mapping of network name → list of reaction dicts.

    Available Networks
    ------------------
    - ``simple_autocatalytic`` : 3-species autocatalytic cycle (baseline)
    - ``cascade`` : 5-species with A→B→C core + D→E→B cascade
    - ``dense`` : 5-species highly connected topology
    - ``oscillatory`` : Brusselator-like oscillatory network
    - ``linear_chain`` : A→B→C→D→E linear chain with reverse
    - ``feedback`` : Negative feedback loop C→D→A
    - ``minimal`` : 2-species minimal network
    - ``complex_pathways`` : Parallel pathways with cross-coupling
    """
    networks = {}

    networks["simple_autocatalytic"] = [
        {"species": ["A", "B"], "reactants": [1, 1], "products": [2, 0], "rate": 1.0, "name": "R0"},
        {"species": ["B", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 0.8, "name": "R1"},
        {"species": ["C", "A"], "reactants": [1, 1], "products": [0, 2], "rate": 0.6, "name": "R2"},
        {"species": ["A", "B"], "reactants": [2, 0], "products": [1, 1], "rate": 0.1, "name": "R3"},
        {"species": ["B", "C"], "reactants": [0, 2], "products": [1, 1], "rate": 0.08, "name": "R4"},
        {"species": ["C", "A"], "reactants": [0, 2], "products": [1, 1], "rate": 0.06, "name": "R5"},
    ]

    networks["cascade"] = [
        {"species": ["A", "B"], "reactants": [1, 1], "products": [2, 0], "rate": 1.2, "name": "R0"},
        {"species": ["B", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 0.9, "name": "R1"},
        {"species": ["C", "A"], "reactants": [1, 1], "products": [0, 2], "rate": 0.7, "name": "R2"},
        {"species": ["A", "B"], "reactants": [2, 0], "products": [1, 1], "rate": 0.15, "name": "R3"},
        {"species": ["B", "C"], "reactants": [0, 2], "products": [1, 1], "rate": 0.12, "name": "R4"},
        {"species": ["C", "A"], "reactants": [0, 2], "products": [1, 1], "rate": 0.09, "name": "R5"},
        {"species": ["A", "D"], "reactants": [1, 0], "products": [0, 1], "rate": 0.4, "name": "R6"},
        {"species": ["D", "E"], "reactants": [1, 0], "products": [0, 1], "rate": 0.35, "name": "R7"},
        {"species": ["E", "B"], "reactants": [1, 0], "products": [0, 1], "rate": 0.3, "name": "R8"},
    ]

    networks["dense"] = [
        {"species": ["A", "B"], "reactants": [1, 1], "products": [2, 0], "rate": 1.0, "name": "R0"},
        {"species": ["B", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 0.9, "name": "R1"},
        {"species": ["C", "D"], "reactants": [1, 1], "products": [0, 2], "rate": 0.8, "name": "R2"},
        {"species": ["D", "E"], "reactants": [1, 1], "products": [0, 2], "rate": 0.7, "name": "R3"},
        {"species": ["E", "A"], "reactants": [1, 1], "products": [0, 2], "rate": 0.6, "name": "R4"},
        {"species": ["A", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 0.5, "name": "R5"},
        {"species": ["B", "D"], "reactants": [1, 1], "products": [0, 2], "rate": 0.4, "name": "R6"},
        {"species": ["C", "E"], "reactants": [1, 1], "products": [0, 2], "rate": 0.3, "name": "R7"},
        {"species": ["A", "B"], "reactants": [2, 0], "products": [1, 1], "rate": 0.1, "name": "R8"},
        {"species": ["B", "C"], "reactants": [0, 2], "products": [1, 1], "rate": 0.09, "name": "R9"},
        {"species": ["C", "D"], "reactants": [0, 2], "products": [1, 1], "rate": 0.08, "name": "R10"},
        {"species": ["D", "E"], "reactants": [0, 2], "products": [1, 1], "rate": 0.07, "name": "R11"},
    ]

    networks["oscillatory"] = [
        {"species": ["A", "X"], "reactants": [1, 0], "products": [0, 1], "rate": 1.0, "name": "R0"},
        {"species": ["X", "Y"], "reactants": [2, 1], "products": [3, 0], "rate": 1.0, "name": "R1"},
        {"species": ["B", "X"], "reactants": [1, 1], "products": [0, 2], "rate": 1.0, "name": "R2"},
        {"species": ["X", "C"], "reactants": [1, 0], "products": [0, 1], "rate": 1.0, "name": "R3"},
        {"species": ["Y", "D"], "reactants": [1, 0], "products": [0, 1], "rate": 0.5, "name": "R4"},
        {"species": ["D", "X"], "reactants": [1, 0], "products": [0, 1], "rate": 0.3, "name": "R5"},
    ]

    networks["linear_chain"] = [
        {"species": ["A", "B"], "reactants": [1, 0], "products": [0, 1], "rate": 1.0, "name": "R0"},
        {"species": ["B", "C"], "reactants": [1, 0], "products": [0, 1], "rate": 0.9, "name": "R1"},
        {"species": ["C", "D"], "reactants": [1, 0], "products": [0, 1], "rate": 0.8, "name": "R2"},
        {"species": ["D", "E"], "reactants": [1, 0], "products": [0, 1], "rate": 0.7, "name": "R3"},
        {"species": ["B", "A"], "reactants": [1, 0], "products": [0, 1], "rate": 0.1, "name": "R4"},
        {"species": ["C", "B"], "reactants": [1, 0], "products": [0, 1], "rate": 0.09, "name": "R5"},
        {"species": ["D", "C"], "reactants": [1, 0], "products": [0, 1], "rate": 0.08, "name": "R6"},
        {"species": ["E", "D"], "reactants": [1, 0], "products": [0, 1], "rate": 0.07, "name": "R7"},
    ]

    networks["feedback"] = [
        {"species": ["A", "B"], "reactants": [1, 1], "products": [2, 0], "rate": 1.2, "name": "R0"},
        {"species": ["B", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 1.0, "name": "R1"},
        {"species": ["C", "A"], "reactants": [1, 0], "products": [0, 1], "rate": 0.8, "name": "R2"},
        {"species": ["A", "B"], "reactants": [2, 0], "products": [1, 1], "rate": 0.2, "name": "R3"},
        {"species": ["B", "C"], "reactants": [0, 2], "products": [1, 1], "rate": 0.15, "name": "R4"},
        {"species": ["C", "D"], "reactants": [1, 0], "products": [0, 1], "rate": 0.5, "name": "R5"},
        {"species": ["D", "A"], "reactants": [1, 1], "products": [0, 0], "rate": 0.3, "name": "R6"},
        {"species": ["D", "E"], "reactants": [1, 0], "products": [0, 1], "rate": 0.4, "name": "R7"},
    ]

    networks["minimal"] = [
        {"species": ["A", "B"], "reactants": [1, 1], "products": [2, 0], "rate": 1.0, "name": "R0"},
        {"species": ["B", "A"], "reactants": [1, 0], "products": [0, 1], "rate": 0.5, "name": "R1"},
        {"species": ["A", "B"], "reactants": [2, 0], "products": [1, 1], "rate": 0.1, "name": "R2"},
    ]

    networks["complex_pathways"] = [
        {"species": ["A", "B"], "reactants": [1, 1], "products": [2, 0], "rate": 1.0, "name": "R0"},
        {"species": ["B", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 0.8, "name": "R1"},
        {"species": ["A", "D"], "reactants": [1, 0], "products": [0, 1], "rate": 0.6, "name": "R2"},
        {"species": ["D", "C"], "reactants": [1, 1], "products": [0, 2], "rate": 0.5, "name": "R3"},
        {"species": ["B", "D"], "reactants": [1, 1], "products": [0, 2], "rate": 0.3, "name": "R4"},
        {"species": ["C", "E"], "reactants": [1, 0], "products": [0, 1], "rate": 0.4, "name": "R5"},
        {"species": ["A", "B"], "reactants": [2, 0], "products": [1, 1], "rate": 0.1, "name": "R6"},
        {"species": ["B", "C"], "reactants": [0, 2], "products": [1, 1], "rate": 0.08, "name": "R7"},
        {"species": ["D", "A"], "reactants": [1, 0], "products": [0, 1], "rate": 0.06, "name": "R8"},
        {"species": ["E", "C"], "reactants": [1, 0], "products": [0, 1], "rate": 0.04, "name": "R9"},
    ]

    return networks


def get_network_io_config(network_name: str) -> Tuple[List[str], List[str]]:
    """
    Return the default (input_channels, output_channels) for a built-in network.

    Parameters
    ----------
    network_name : str
        Name of a built-in network.

    Returns
    -------
    input_channels : list of str
    output_channels : list of str
    """
    configs = {
        "simple_autocatalytic": (["A", "B"], ["A", "B", "C"]),
        "cascade": (["A", "B", "C"], ["A", "B", "D"]),
        "dense": (["A", "B"], ["A", "B", "C", "D"]),
        "oscillatory": (["A", "B"], ["X", "Y", "D"]),
        "linear_chain": (["A"], ["A", "B", "C", "D"]),
        "feedback": (["A", "B"], ["A", "C", "D"]),
        "minimal": (["A"], ["A", "B"]),
        "complex_pathways": (["A", "B"], ["B", "C", "D", "E"]),
    }
    return configs.get(network_name, (["A", "B"], ["A", "B", "C"]))


def list_builtin_networks() -> List[str]:
    """Return names of all built-in reaction networks."""
    return list(get_builtin_networks().keys())
