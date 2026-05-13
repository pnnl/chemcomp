"""
Parameterized multi-subnetwork CRN generator.

Generates networks with configurable:
  - Number of subnetworks (scaling width)
  - Complexity per subnetwork (scaling depth: extra intermediate species + reactions)
  - Leakage rate between compartments
  - Fuel feed-in rate

Usage:
    from gen_parameterized_network import generate_network
    network = generate_network(n_sub=8, complexity='high', leak_rate=0.005)
"""
import json
import numpy as np
import os


def generate_network(n_sub=6, complexity='standard', leak_rate=0.005,
                     fuel_feed=2.0, fuel_consumption=0.05, seed=123):
    """Generate a multi-subnetwork CRN with configurable architecture.

    Parameters
    ----------
    n_sub : int
        Number of Brusselator-like subnetworks.
    complexity : str
        'minimal' — standard Brusselator (4 species/module: X, Y, B, D)
        'standard' — Brusselator + fuel coupling (4 species/module + 1 fuel)
        'high' — adds intermediate species Z per module (5 species/module + fuel)
                  with extra reactions: X -> Z -> Y pathway and Z autocatalysis
    leak_rate : float
        Base leakage rate between adjacent compartments.
    fuel_feed : float
        Constant production rate of shared fuel F.
    fuel_consumption : float
        Rate at which fuel activates each subnetwork.
    seed : int
        Random seed for rate perturbations.

    Returns
    -------
    dict
        Network configuration dict ready for JSON serialization.
    """
    rng = np.random.RandomState(seed)

    # Generate diverse oscillation parameters
    base_params = [
        (1.0, 3.0, 1.0), (0.8, 2.8, 1.2), (1.2, 3.5, 0.9),
        (0.7, 2.2, 1.1), (1.1, 3.2, 1.0), (0.9, 2.6, 1.3),
        (0.6, 2.0, 1.0), (1.3, 3.8, 0.8), (0.75, 2.4, 1.15),
        (1.05, 3.1, 1.05), (0.85, 2.7, 1.25), (1.15, 3.4, 0.95),
    ]
    # Cycle through if more subnetworks than base params
    params = [base_params[i % len(base_params)] for i in range(n_sub)]

    species_all = []
    modules = []

    for m in range(n_sub):
        mod = {"X": f"X{m+1}", "Y": f"Y{m+1}", "B": f"B{m+1}", "D": f"D{m+1}"}
        species_all.extend([mod["X"], mod["Y"], mod["B"], mod["D"]])
        if complexity == 'high':
            mod["Z"] = f"Z{m+1}"
            species_all.append(mod["Z"])
        modules.append(mod)

    # Shared fuel
    FUEL = "F"
    species_all.append(FUEL)

    reactions = []

    def add_rxn(sp_list, reactants, products, rate):
        reactions.append({
            "species": sp_list,
            "reactants": [int(x) for x in reactants],
            "products": [int(x) for x in products],
            "rate": round(float(rate), 6),
        })

    # --- Intra-subnetwork reactions ---
    for m, (a, b_target, c_rate) in enumerate(params):
        mod = modules[m]

        # Core Brusselator
        add_rxn([mod["X"]], [0], [1], a)
        add_rxn([mod["X"], mod["Y"]], [2, 1], [3, 0], c_rate)
        add_rxn([mod["B"], mod["X"], mod["Y"], mod["D"]], [1, 1, 0, 0], [0, 0, 1, 1], 1.0)
        add_rxn([mod["X"], mod["D"]], [1, 0], [0, 1], 1.0)

        # Buffered B_m
        add_rxn([mod["B"]], [0], [1], b_target * 2.0)
        add_rxn([mod["B"]], [1], [0], 2.0)

        # Waste management
        add_rxn([mod["D"], mod["B"]], [1, 0], [0, 1], 0.08)
        add_rxn([mod["D"]], [1], [0], 0.04)

        # Fuel coupling
        add_rxn([FUEL, mod["X"]], [1, 1], [0, 2], fuel_consumption)

        # High-complexity: intermediate Z species
        if complexity == 'high':
            Z = mod["Z"]
            # X -> Z (slow conversion to intermediate)
            add_rxn([mod["X"], Z], [1, 0], [0, 1], 0.3 + 0.1 * rng.rand())
            # Z + Y -> 2Y (Z catalyzes Y production — amplification loop)
            add_rxn([Z, mod["Y"]], [1, 1], [0, 2], 0.2 + 0.1 * rng.rand())
            # Z -> D (intermediate decay)
            add_rxn([Z, mod["D"]], [1, 0], [0, 1], 0.15 + 0.05 * rng.rand())
            # Z autocatalysis: 2Z -> 3Z (weak)
            add_rxn([Z], [2], [3], 0.05 + 0.02 * rng.rand())
            # Z decay to prevent blowup
            add_rxn([Z], [1], [0], 0.5)

    # --- Fuel dynamics ---
    add_rxn([FUEL], [0], [1], fuel_feed)
    add_rxn([FUEL], [1], [0], fuel_feed / 3.0)

    # --- Inter-subnetwork leakage ---
    for m1 in range(n_sub):
        m2 = (m1 + 1) % n_sub
        rate_x = leak_rate * (1 + 0.3 * rng.rand())
        add_rxn([modules[m1]["X"], modules[m2]["X"]], [1, 0], [0, 1], rate_x)
        rate_y = leak_rate * (1 + 0.3 * rng.rand())
        add_rxn([modules[m1]["Y"], modules[m2]["Y"]], [1, 0], [0, 1], rate_y)

        if complexity == 'high':
            # Z leakage too
            rate_z = leak_rate * 0.5 * (1 + 0.3 * rng.rand())
            add_rxn([modules[m1]["Z"], modules[m2]["Z"]], [1, 0], [0, 1], rate_z)

    # Skip connections
    for m1 in range(n_sub):
        m2 = (m1 + 2) % n_sub
        if m2 <= m1:
            continue
        rate_skip = leak_rate * 0.5 * (1 + 0.2 * rng.rand())
        add_rxn([modules[m1]["X"], modules[m2]["X"]], [1, 0], [0, 1], rate_skip)

    # Cross-inhibition
    for m1 in range(n_sub):
        m2 = (m1 + 3) % n_sub
        if n_sub <= 3:
            continue  # skip for very small networks
        rate_inhib = leak_rate * 0.3 * (1 + 0.5 * rng.rand())
        add_rxn(
            [modules[m1]["Y"], modules[m2]["X"], modules[m1]["X"], modules[m2]["Y"]],
            [1, 1, 0, 0], [0, 0, 1, 1], rate_inhib,
        )

    # --- Build config ---
    input_channels = [modules[m]["B"] for m in range(n_sub)] + [FUEL]
    output_channels = species_all[:]

    name = f"multisubnet_n{n_sub}_{complexity}"
    description = (
        f"Multi-subnetwork CRN ({len(species_all)} species, {len(reactions)} reactions): "
        f"{n_sub} Brusselator modules, complexity={complexity}, "
        f"leak_rate={leak_rate}, fuel_feed={fuel_feed}"
    )

    return {
        "name": name,
        "description": description,
        "reactions": reactions,
        "input_channels": input_channels,
        "output_channels": output_channels,
        "_meta": {
            "n_sub": n_sub,
            "complexity": complexity,
            "n_species": len(species_all),
            "n_reactions": len(reactions),
        }
    }


def save_network(network, outdir=None):
    """Save network to JSON file."""
    if outdir is None:
        outdir = os.path.join(os.path.dirname(__file__), "..", "network_configs")
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"{network['name']}.json")
    # Remove _meta before saving
    to_save = {k: v for k, v in network.items() if k != '_meta'}
    with open(outpath, "w") as f:
        json.dump(to_save, f, indent=2)
    return outpath


if __name__ == "__main__":
    configs = [
        # Width scaling: vary n_sub with standard complexity
        (4, 'standard'),
        (6, 'standard'),
        (8, 'standard'),
        (10, 'standard'),
        (12, 'standard'),
        # Depth scaling: vary complexity with fixed n_sub=6
        (6, 'minimal'),
        (6, 'standard'),
        (6, 'high'),
        # Large complex network
        (10, 'high'),
    ]

    for n_sub, complexity in configs:
        net = generate_network(n_sub=n_sub, complexity=complexity)
        path = save_network(net)
        meta = net['_meta']
        print(f"  {meta['n_species']:3d} species, {meta['n_reactions']:3d} reactions  ->  {path}")
