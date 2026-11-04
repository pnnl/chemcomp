#!/usr/bin/env python3
"""Check that a compiled reaction network reproduces its ODE system.

The compiler's ODE-to-network flow claims that the mass-action kinetics of the
network it emits are exactly the input ODEs. This script tests that claim
numerically:

  1. solve the reference ODE (``chemcomp-translate -ode-to-json``);
  2. simulate the network (``chemcomp-translate -mlir-to-json``) with the
     emulator's mass-action dynamics (``ConstrainedChemicalReservoir``);
  3. map species back to ODE states with the JSON ``species_map``: a dual-rail
     pair contributes ``x = x_p - x_m``;
  4. report the mean squared error (MSE) and the maximum absolute error per
     state, and fail when any MSE exceeds ``--tol``.

Both systems are integrated with tight tolerances, so the error measures the
realization, not the integrator.

Example::

    chemcomp-opt model.mlir --ode-to-crn | chemcomp-translate -mlir-to-json > crn.json
    chemcomp-translate -ode-to-json model.mlir > ode.json
    python scripts/check_realization.py --ode ode.json --crn crn.json --t-end 10
"""

import argparse
import json
import os
import sys
import tempfile

import numpy as np
from scipy.integrate import solve_ivp

from chemical_reservoir.constraints import PhysicalConstraints
from chemical_reservoir.io import load_reaction_network
from chemical_reservoir.reservoir import ConstrainedChemicalReservoir

RTOL = 1e-10
ATOL = 1e-12


def reference_rhs(ode):
    """Right-hand side of the reference ODE as a function of the state vector."""
    index = {state["name"]: i for i, state in enumerate(ode["states"])}
    terms = [
        (index[t["state"]], t["coeff"], [(index[s], e) for s, e in t["monomial"].items()])
        for t in ode["terms"]
    ]

    def rhs(_t, y):
        dydt = np.zeros(len(index))
        for target, coeff, monomial in terms:
            value = coeff
            for i, exponent in monomial:
                value *= y[i] ** exponent
            dydt[target] += value
        return dydt

    return rhs


def simulate_network(crn_path, t_eval):
    """Simulate the network; return species names and trajectories (time x species)."""
    config = load_reaction_network(crn_path)
    # The defaults clip concentrations at 3-5; the check must not clip.
    physical = PhysicalConstraints(default_max_concentration=1e12)
    reservoir = ConstrainedChemicalReservoir(
        reactions=config["reactions"],
        input_channels=[],
        output_channels=[],
        physical_constraints=physical,
    )
    populations = config.get("populations", {})
    y0 = np.array([float(populations.get(sp, 0.0)) for sp in reservoir.species_names])
    sol = solve_ivp(
        lambda t, y: reservoir._dynamics(t, y),
        (t_eval[0], t_eval[-1]),
        y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=RTOL,
        atol=ATOL,
    )
    if not sol.success:
        raise RuntimeError(f"network simulation failed: {sol.message}")
    return reservoir.species_names, sol.y.T, config.get("species_map", {})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--ode", required=True, help="reference ODE (ode-to-json)")
    parser.add_argument("--crn", required=True, help="network (mlir-to-json)")
    parser.add_argument("--t-end", type=float, default=10.0, help="simulated time")
    parser.add_argument("--points", type=int, default=201, help="samples compared")
    parser.add_argument("--tol", type=float, default=1e-6, help="maximum MSE per state")
    args = parser.parse_args(argv)

    with open(args.ode) as f:
        ode = json.load(f)
    t_eval = np.linspace(0.0, args.t_end, args.points)

    y0 = np.array([s["initial"] for s in ode["states"]], dtype=float)
    ref = solve_ivp(
        reference_rhs(ode),
        (0.0, args.t_end),
        y0,
        t_eval=t_eval,
        method="LSODA",
        rtol=RTOL,
        atol=ATOL,
    )
    if not ref.success:
        print(f"reference ODE failed: {ref.message}", file=sys.stderr)
        return 2

    species, traj, species_map = simulate_network(args.crn, t_eval)

    print(f"{ode['name']}: {len(species)} species, t in [0, {args.t_end:g}], {args.points} samples")
    print(f"  {'state':<12} {'species':<24} {'MSE':>10} {'max |err|':>10}")
    worst = 0.0
    for i, state in enumerate(ode["states"]):
        name = state["name"]
        members = [
            (j, species_map.get(sp, {}).get("sign", 1))
            for j, sp in enumerate(species)
            if species_map.get(sp, {}).get("state", sp) == name
        ]
        if not members:
            print(f"  {name:<12} no species encodes this state", file=sys.stderr)
            return 2
        value = sum(sign * traj[:, j] for j, sign in members)
        err = value - ref.y[i]
        mse = float(np.mean(err**2))
        worst = max(worst, mse)
        encoding = " - ".join(species[j] for j, _ in sorted(members, key=lambda m: -m[1]))
        print(f"  {name:<12} {encoding:<24} {mse:10.3e} {np.max(np.abs(err)):10.3e}")

    status = "PASS" if worst <= args.tol else "FAIL"
    print(f"  {status}: worst MSE {worst:.3e} (tolerance {args.tol:g})")
    return 0 if status == "PASS" else 1


def run_in_scratch_dir(argv=None):
    """Run main() from a scratch directory: the emulator creates its figure
    directory in the working directory, and the check should leave nothing
    behind. Paths given on the command line are made absolute first."""
    argv = list(sys.argv[1:] if argv is None else argv)
    for i, arg in enumerate(argv[:-1]):
        if arg in ("--ode", "--crn"):
            argv[i + 1] = os.path.abspath(argv[i + 1])
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as scratch:
        os.chdir(scratch)
        try:
            return main(argv)
        finally:
            os.chdir(cwd)


if __name__ == "__main__":
    sys.exit(run_in_scratch_dir())
