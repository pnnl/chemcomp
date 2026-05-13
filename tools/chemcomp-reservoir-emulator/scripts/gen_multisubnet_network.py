"""
Generate a multi-subnetwork CRN reservoir with inter-compartment leakage
and constant reagent feed-in.

Architecture:
  - N_SUB independent Brusselator-like oscillator subnetworks, each with
    different oscillation frequencies.
  - Small leakage reactions between adjacent compartments model diffusive
    coupling (X_i -> X_j at low rate).
  - A shared "fuel" species F is continuously produced (constant feed-in)
    and slowly consumed by every subnetwork.  This provides a constant
    energy source that keeps the reservoir in an excitable regime and
    couples all subnetworks globally.

The leakage and fuel feed give the reservoir richer transient dynamics
and longer memory horizons compared to strongly-coupled networks.
"""
import json
import numpy as np
import os

rng = np.random.RandomState(123)

# ==============================================================================
# Configuration
# ==============================================================================

N_SUB = 6           # Number of subnetworks
LEAK_RATE = 0.005   # Leakage rate between adjacent compartments
FUEL_FEED = 2.0     # Constant production rate of shared fuel F
FUEL_CONSUMPTION = 0.05  # Rate at which each subnetwork consumes fuel

# Brusselator parameters per subnetwork: (a, b_eff, c)
# Oscillation requires b_eff > 1 + a²
SUBNETWORK_PARAMS = [
    (1.0, 3.0, 1.0),   # Sub 1: standard
    (0.8, 2.8, 1.2),   # Sub 2: faster nonlinearity
    (1.2, 3.5, 0.9),   # Sub 3: high b
    (0.7, 2.2, 1.1),   # Sub 4: near bifurcation
    (1.1, 3.2, 1.0),   # Sub 5: intermediate
    (0.9, 2.6, 1.3),   # Sub 6: high nonlinearity
]

# ==============================================================================
# Build species list
# ==============================================================================

species_all = []
modules = []

for m in range(N_SUB):
    X = f"X{m+1}"
    Y = f"Y{m+1}"
    B = f"B{m+1}"
    D = f"D{m+1}"
    modules.append({"X": X, "Y": Y, "B": B, "D": D})
    species_all.extend([X, Y, B, D])

# Shared fuel species
FUEL = "F"
species_all.append(FUEL)

print(f"Species ({len(species_all)}): {species_all}")

# ==============================================================================
# Build reactions
# ==============================================================================

reactions = []


def add_rxn(sp_list, reactants, products, rate):
    reactions.append({
        "species": sp_list,
        "reactants": [int(x) for x in reactants],
        "products": [int(x) for x in products],
        "rate": round(float(rate), 6),
    })


# --- Intra-subnetwork reactions (Brusselator) ---
for m, (a, b_target, c_rate) in enumerate(SUBNETWORK_PARAMS):
    mod = modules[m]

    # (1) Source: 0 -> X_m   (rate = a)
    add_rxn([mod["X"]], [0], [1], a)

    # (2) Autocatalysis: 2X_m + Y_m -> 3X_m   (rate = c)
    add_rxn([mod["X"], mod["Y"]], [2, 1], [3, 0], c_rate)

    # (3) Conversion: B_m + X_m -> Y_m + D_m   (rate = 1.0)
    add_rxn([mod["B"], mod["X"], mod["Y"], mod["D"]], [1, 1, 0, 0], [0, 0, 1, 1], 1.0)

    # (4) Decay: X_m -> D_m   (rate = 1.0)
    add_rxn([mod["X"], mod["D"]], [1, 0], [0, 1], 1.0)

    # (5) Buffered B_m: source + decay to maintain B_ss ≈ b_target
    b_source = b_target * 2.0
    b_decay = 2.0
    add_rxn([mod["B"]], [0], [1], b_source)
    add_rxn([mod["B"]], [1], [0], b_decay)

    # (6) Waste recycling and decay
    add_rxn([mod["D"], mod["B"]], [1, 0], [0, 1], 0.08)
    add_rxn([mod["D"]], [1], [0], 0.04)

    # (7) Fuel-dependent activation: F + X_m -> 2X_m
    #     Fuel catalyzes activator production — couples all subnetworks
    add_rxn([FUEL, mod["X"]], [1, 1], [0, 2], FUEL_CONSUMPTION)


# --- Constant fuel feed-in: 0 -> F ---
add_rxn([FUEL], [0], [1], FUEL_FEED)

# --- Fuel natural decay: F -> 0 ---
# Steady state F_ss = FUEL_FEED / fuel_decay_rate
fuel_decay_rate = FUEL_FEED / 3.0  # F_ss ≈ 3.0
add_rxn([FUEL], [1], [0], fuel_decay_rate)


# --- Inter-subnetwork leakage (adjacent compartments) ---
# Small diffusive exchange: X_i -> X_j and Y_i -> Y_j
for m1 in range(N_SUB):
    m2 = (m1 + 1) % N_SUB  # circular topology

    # X leakage: X_m1 -> X_m2  (and reverse)
    rate_x = LEAK_RATE * (1 + 0.3 * rng.rand())
    add_rxn([modules[m1]["X"], modules[m2]["X"]], [1, 0], [0, 1], rate_x)

    # Y leakage: Y_m1 -> Y_m2
    rate_y = LEAK_RATE * (1 + 0.3 * rng.rand())
    add_rxn([modules[m1]["Y"], modules[m2]["Y"]], [1, 0], [0, 1], rate_y)

# --- Skip connections (every other subnetwork) ---
for m1 in range(N_SUB):
    m2 = (m1 + 2) % N_SUB
    if m2 <= m1:
        continue  # avoid duplicates
    rate_skip = LEAK_RATE * 0.5 * (1 + 0.2 * rng.rand())
    add_rxn([modules[m1]["X"], modules[m2]["X"]], [1, 0], [0, 1], rate_skip)


# --- Cross-inhibition via Y (weak) ---
# Y from one subnetwork converts X of a non-adjacent one
for m1 in range(N_SUB):
    m2 = (m1 + 3) % N_SUB
    rate_inhib = LEAK_RATE * 0.3 * (1 + 0.5 * rng.rand())
    add_rxn(
        [modules[m1]["Y"], modules[m2]["X"], modules[m1]["X"], modules[m2]["Y"]],
        [1, 1, 0, 0], [0, 0, 1, 1], rate_inhib,
    )


# ==============================================================================
# Network configuration
# ==============================================================================

print(f"\nTotal reactions: {len(reactions)}")
print(f"Total species:   {len(species_all)}")

# Input channels: B_m species (control bifurcation parameter) + Fuel
input_channels = [modules[m]["B"] for m in range(N_SUB)] + [FUEL]

# Output channels: all species
output_channels = species_all[:]

print(f"Input channels ({len(input_channels)}):  {input_channels}")
print(f"Output channels ({len(output_channels)}): {output_channels}")

network = {
    "name": "multi_subnet_leaky",
    "description": (
        f"Multi-subnetwork CRN reservoir ({len(species_all)} species, "
        f"{len(reactions)} reactions) with {N_SUB} weakly-coupled Brusselator "
        f"oscillators, inter-compartment leakage (rate={LEAK_RATE}), "
        f"and constant fuel feed-in (F_ss≈3.0)."
    ),
    "reactions": reactions,
    "input_channels": input_channels,
    "output_channels": output_channels,
}

outdir = os.path.join(os.path.dirname(__file__), "..", "network_configs")
outpath = os.path.join(outdir, "multi_subnet_leaky.json")
with open(outpath, "w") as f:
    json.dump(network, f, indent=2)

print(f"\nSaved to {outpath}")
print(f"\nArchitecture summary:")
print(f"  {N_SUB} Brusselator subnetworks (circular coupling)")
print(f"  Leakage rate: {LEAK_RATE} (X and Y diffusion between neighbours)")
print(f"  Skip connections: every-other at {LEAK_RATE*0.5:.4f}")
print(f"  Cross-inhibition: Y→X at rate ~{LEAK_RATE*0.3:.4f}")
print(f"  Shared fuel F: feed={FUEL_FEED}, steady-state≈3.0")
print(f"  Fuel activation: F + X_m -> 2X_m at rate {FUEL_CONSUMPTION}")
