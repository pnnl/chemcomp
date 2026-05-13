"""
Generate a high-dimensional oscillatory reaction network for reservoir computing.

Uses 5 coupled Brusselator oscillators.  The standard Brusselator produces
limit-cycle oscillations via mass-action kinetics:

    (1) 0 -> X           rate = a
    (2) 2X + Y -> 3X     rate = c           (trimolecular autocatalysis)
    (3) B + X  -> Y + D  rate = b_eff       (linear conversion)
    (4) X -> D            rate = 1           (decay)

giving dx/dt = a - (b_eff + 1)x + c x² y,  dy/dt = b_eff x - c x² y.
Oscillation condition: b_eff > 1 + a².

Each module has different (a, b_eff, c) to produce different periods.
Modules are weakly coupled by Y-X exchange reactions.
"""
import json
import numpy as np
import os

rng = np.random.RandomState(42)

# ---- 5 Brusselator modules ----
# Species per module: X_m, Y_m, B_m (held approximately constant), D_m (waste)
# We keep B_m constant via a fast source + decay pair.

n_modules = 5
modules = []
species_all = []
for m in range(n_modules):
    X = f"X{m+1}"
    Y = f"Y{m+1}"
    B = f"B{m+1}"
    D = f"D{m+1}"
    modules.append({"X": X, "Y": Y, "B": B, "D": D})
    species_all.extend([X, Y, B, D])

reactions = []


def add_rxn(sp_list, reactants, products, rate):
    reactions.append({
        "species": sp_list,
        "reactants": [int(x) for x in reactants],
        "products": [int(x) for x in products],
        "rate": round(float(rate), 6),
    })


# Brusselator parameters per module
# (a, b_target, c_rate)
# b_eff = b_target (we'll set B_ss ≈ b_target and k_conv = 1)
# Need b_eff > 1 + a² for oscillation
# Period ≈ 2π / sqrt(b - 1 - a²) * factor
module_params = [
    (1.0, 3.0, 1.0),    # Module 1: b=3 > 1+1=2 ✓, period ~6s
    (1.0, 2.5, 1.0),    # Module 2: b=2.5 > 2 ✓, period ~9s
    (0.8, 2.5, 1.0),    # Module 3: b=2.5 > 1.64 ✓, period ~7s
    (1.2, 3.5, 1.0),    # Module 4: b=3.5 > 2.44 ✓, period ~6s
    (0.6, 2.0, 1.0),    # Module 5: b=2.0 > 1.36 ✓, period ~10s
]

for m, (a, b_target, c_rate) in enumerate(module_params):
    mod = modules[m]

    # (1) Source: 0 -> X_m   (rate = a)
    add_rxn([mod["X"]], [0], [1], a)

    # (2) Autocatalysis: 2 X_m + Y_m -> 3 X_m   (rate = c)
    add_rxn([mod["X"], mod["Y"]], [2, 1], [3, 0], c_rate)

    # (3) Conversion: B_m + X_m -> Y_m + D_m   (rate = 1.0)
    #     With B_ss ≈ b_target, effective rate = b_target * x
    add_rxn([mod["B"], mod["X"], mod["Y"], mod["D"]], [1, 1, 0, 0], [0, 0, 1, 1], 1.0)

    # (4) Decay: X_m -> D_m   (rate = 1.0)
    add_rxn([mod["X"], mod["D"]], [1, 0], [0, 1], 1.0)

    # (5) Buffered B_m: 0 -> B_m (fast source) and B_m -> 0 (fast decay)
    #     Equilibrium B_ss = source_rate / decay_rate = b_target
    #     Use moderate rates so B responds slowly to perturbation
    b_source = b_target * 2.0
    b_decay = 2.0
    add_rxn([mod["B"]], [0], [1], b_source)
    add_rxn([mod["B"]], [1], [0], b_decay)

    # (6) Slow waste recycling: D_m -> B_m  (prevents D from accumulating forever)
    add_rxn([mod["D"], mod["B"]], [1, 0], [0, 1], 0.1)

    # (7) Slow D decay: D_m -> 0
    add_rxn([mod["D"]], [1], [0], 0.05)


# ---- Inter-module coupling ----
# Weak Y-X exchange between neighboring modules
for m1 in range(n_modules):
    m2 = (m1 + 1) % n_modules
    # Y_m1 + X_m2 -> X_m1 + Y_m2  (inhibitor from m1 converts activator of m2)
    rate = 0.02 + 0.01 * rng.rand()
    add_rxn(
        [modules[m1]["Y"], modules[m2]["X"], modules[m1]["X"], modules[m2]["Y"]],
        [1, 1, 0, 0], [0, 0, 1, 1], rate,
    )

# Skip coupling (modules i <-> i+2)
for m1 in range(n_modules):
    m2 = (m1 + 2) % n_modules
    rate = 0.01 + 0.005 * rng.rand()
    add_rxn(
        [modules[m1]["X"], modules[m2]["Y"], modules[m1]["Y"], modules[m2]["X"]],
        [1, 1, 0, 0], [0, 0, 1, 1], rate,
    )


print(f"Total reactions: {len(reactions)}")
print(f"Total species:   {len(species_all)}")

# Input channels: B_m species (control the bifurcation parameter)
input_channels = [modules[m]["B"] for m in range(n_modules)]

# Output channels: all species for maximum feature diversity
output_channels = species_all[:]

n_out = len(output_channels)
n_feat = n_out + n_out + n_out * (n_out - 1) // 2
print(f"Output channels: {n_out} -> {n_feat} quadratic features")
print(f"Input channels:  {input_channels}")
print(f"Output channels: {output_channels}")

network = {
    "name": "high_dimensional_random",
    "description": (
        f"High-dimensional Brusselator network ({len(species_all)} species, "
        f"{len(reactions)} reactions) with {n_modules} coupled oscillator "
        "modules for reservoir computing"
    ),
    "reactions": reactions,
    "input_channels": input_channels,
    "output_channels": output_channels,
}

outdir = os.path.join(os.path.dirname(__file__), "..", "network_configs")
outpath = os.path.join(outdir, "high_dimensional_random.json")
with open(outpath, "w") as f:
    json.dump(network, f, indent=2)

print(f"Saved to {outpath}")
