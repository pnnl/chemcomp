"""
Validation utilities for problem YAML files and reaction network JSON files.

Users can call these functions to check that their custom input files are
well-formed and will work with the reservoir before running a full solve.

CLI usage:
    python -m chemical_reservoir.validate problem_configs/my_ode.yaml
    python -m chemical_reservoir.validate network_configs/my_network.json

Python API:
    from chemical_reservoir.validate import validate_problem_file, validate_network_file
    errors = validate_problem_file("my_problem.yaml")
    errors = validate_network_file("my_network.json")
"""

import json
import sys
from pathlib import Path
from typing import List

import numpy as np


def validate_network_file(filepath: str) -> List[str]:
    """
    Validate a reaction network JSON file.

    Checks:
    - File exists and is valid JSON
    - Required top-level keys present (reactions, input_channels, output_channels)
    - Each reaction has species, reactants, products, rate fields
    - Stoichiometry vectors have consistent lengths
    - Rate constants are positive
    - Input/output channels reference species that appear in reactions
    - Network can be loaded and integrated (basic smoke test)

    Args:
        filepath: Path to a reaction network JSON file.

    Returns:
        List of error/warning strings. Empty list means valid.
    """
    errors = []
    filepath = Path(filepath)

    # --- File existence ---
    if not filepath.exists():
        return [f"File not found: {filepath}"]

    # --- JSON parsing ---
    try:
        with open(filepath) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return [f"Invalid JSON: {e}"]

    # --- Top-level structure ---
    required_keys = ["reactions"]
    for key in required_keys:
        if key not in data:
            errors.append(f"Missing required top-level key: '{key}'")

    if errors:
        return errors

    # --- Reactions validation ---
    reactions = data["reactions"]
    if not isinstance(reactions, list) or len(reactions) == 0:
        errors.append("'reactions' must be a non-empty list")
        return errors

    all_species = set()
    for i, rxn in enumerate(reactions):
        # Required fields
        for field in ["species", "reactants", "products", "rate"]:
            if field not in rxn:
                errors.append(f"Reaction {i}: missing field '{field}'")
                continue

        if "species" not in rxn:
            continue

        species = rxn["species"]
        if not isinstance(species, list) or len(species) == 0:
            errors.append(f"Reaction {i}: 'species' must be a non-empty list")
            continue

        all_species.update(species)
        n = len(species)

        if "reactants" in rxn:
            if not isinstance(rxn["reactants"], list):
                errors.append(f"Reaction {i}: 'reactants' must be a list")
            elif len(rxn["reactants"]) != n:
                errors.append(
                    f"Reaction {i}: 'reactants' length ({len(rxn['reactants'])}) "
                    f"!= 'species' length ({n})"
                )
            else:
                for j, val in enumerate(rxn["reactants"]):
                    if not isinstance(val, (int, float)) or val < 0:
                        errors.append(
                            f"Reaction {i}: reactants[{j}] must be a non-negative number"
                        )

        if "products" in rxn:
            if not isinstance(rxn["products"], list):
                errors.append(f"Reaction {i}: 'products' must be a list")
            elif len(rxn["products"]) != n:
                errors.append(
                    f"Reaction {i}: 'products' length ({len(rxn['products'])}) "
                    f"!= 'species' length ({n})"
                )
            else:
                for j, val in enumerate(rxn["products"]):
                    if not isinstance(val, (int, float)) or val < 0:
                        errors.append(
                            f"Reaction {i}: products[{j}] must be a non-negative number"
                        )

        if "rate" in rxn:
            if not isinstance(rxn["rate"], (int, float)):
                errors.append(f"Reaction {i}: 'rate' must be a number")
            elif rxn["rate"] <= 0:
                errors.append(f"Reaction {i}: 'rate' must be positive (got {rxn['rate']})")

    # --- Input/output channels ---
    for channel_key in ["input_channels", "output_channels"]:
        if channel_key in data:
            channels = data[channel_key]
            if not isinstance(channels, list):
                errors.append(f"'{channel_key}' must be a list")
            else:
                for ch in channels:
                    if ch not in all_species:
                        errors.append(
                            f"'{channel_key}' references species '{ch}' "
                            f"not found in any reaction"
                        )

    # --- Try loading ---
    if not errors:
        try:
            from .io import load_reservoir_from_json
            reservoir, _ = load_reservoir_from_json(str(filepath))

            # Quick integration test
            t_test = np.linspace(0, 1, 20)
            from scipy.integrate import solve_ivp

            y0 = reservoir.initial_concentrations
            if y0 is None:
                y0 = np.ones(reservoir.n_species) * 0.1

            def rhs(t, y):
                return reservoir._dynamics(t, y)

            sol = solve_ivp(rhs, [0, 1], y0, t_eval=t_test)
            if not sol.success:
                errors.append(f"Integration test failed: {sol.message}")
            elif np.any(np.isnan(sol.y)) or np.any(np.isinf(sol.y)):
                errors.append(
                    "Integration test produced NaN/Inf — check rate constants "
                    "and initial concentrations"
                )
        except Exception as e:
            errors.append(f"Failed to load/integrate network: {e}")

    return errors


def validate_problem_file(filepath: str) -> List[str]:
    """
    Validate a custom problem YAML file.

    Checks:
    - File exists and is valid YAML
    - Required fields present for the declared type (ode/pde)
    - Equations parse correctly via SymPy
    - Number of equations matches number of variables (ODE)
    - Initial conditions have correct length (ODE)
    - t_span is a valid 2-element list with t_start < t_end
    - Boundary conditions are well-formed (PDE)
    - A short numerical integration succeeds without errors

    Args:
        filepath: Path to a YAML problem definition file.

    Returns:
        List of error/warning strings. Empty list means valid.
    """
    errors = []
    filepath = Path(filepath)

    if not filepath.exists():
        return [f"File not found: {filepath}"]

    # --- YAML parsing ---
    try:
        import yaml
        with open(filepath) as f:
            spec = yaml.safe_load(f)
    except ImportError:
        return ["PyYAML is not installed: pip install pyyaml"]
    except Exception as e:
        return [f"Invalid YAML: {e}"]

    if not isinstance(spec, dict):
        return ["YAML file must contain a mapping (dictionary) at the top level"]

    # --- Type ---
    problem_type = spec.get("type", "").lower()
    if problem_type not in ("ode", "pde"):
        errors.append(
            f"'type' must be 'ode' or 'pde' (got '{spec.get('type', '<missing>')}')"
        )
        return errors

    # --- Common fields ---
    if "name" not in spec:
        errors.append("Missing 'name' field (recommended for identification)")

    if "variables" not in spec:
        errors.append("Missing 'variables' field — list the dependent variable names")

    t_span = spec.get("t_span")
    if t_span is None:
        errors.append("Missing 't_span' field (e.g., [0, 10])")
    elif not isinstance(t_span, list) or len(t_span) != 2:
        errors.append("'t_span' must be a 2-element list [t_start, t_end]")
    else:
        try:
            t0, t1 = float(t_span[0]), float(t_span[1])
            if t0 >= t1:
                errors.append(f"'t_span': t_start ({t0}) must be < t_end ({t1})")
        except (TypeError, ValueError):
            errors.append("'t_span' elements must be numbers")

    # --- Parameters ---
    params = spec.get("parameters", {})
    if not isinstance(params, dict):
        errors.append("'parameters' must be a mapping (key: value pairs)")
        params = {}

    for k, v in params.items():
        try:
            float(v)
        except (TypeError, ValueError):
            errors.append(f"Parameter '{k}' value must be a number (got {v!r})")

    # --- Type-specific validation ---
    if problem_type == "ode":
        errors.extend(_validate_ode_spec(spec, params))
    else:
        errors.extend(_validate_pde_spec(spec, params))

    # --- Try full load if no blocking errors ---
    if not errors:
        try:
            from .custom_problems import load_problem_from_yaml
            problem = load_problem_from_yaml(str(filepath))

            # Short integration test
            t_test = np.linspace(
                float(spec["t_span"][0]),
                float(spec["t_span"][1]),
                min(30, 10),
            )
            if problem_type == "ode":
                sol = problem.solve_numerical(t_test)
                if np.any(np.isnan(sol)) or np.any(np.isinf(sol)):
                    errors.append(
                        "Integration produced NaN/Inf — check equations and parameters"
                    )
            else:
                x, sol = problem.solve_numerical(t_test)
                if np.any(np.isnan(sol)) or np.any(np.isinf(sol)):
                    errors.append(
                        "Integration produced NaN/Inf — check equation and parameters"
                    )
        except Exception as e:
            errors.append(f"Load/integration failed: {e}")

    return errors


def _validate_ode_spec(spec: dict, params: dict) -> List[str]:
    """Validate ODE-specific fields."""
    errors = []

    variables = spec.get("variables", [])
    if not isinstance(variables, list) or len(variables) == 0:
        errors.append("'variables' must be a non-empty list of variable names")
        return errors

    equations = spec.get("equations")
    if equations is None:
        errors.append("Missing 'equations' field — list of RHS expressions")
        return errors
    if not isinstance(equations, list):
        errors.append("'equations' must be a list of strings")
        return errors

    if len(equations) != len(variables):
        errors.append(
            f"Number of equations ({len(equations)}) must match "
            f"number of variables ({len(variables)})"
        )

    # Check each equation parses
    try:
        import sympy
        from sympy.parsing.sympy_parser import (
            parse_expr, standard_transformations,
            implicit_multiplication_application,
        )
        transformations = standard_transformations + (implicit_multiplication_application,)

        local_dict = {
            'pi': sympy.pi, 'e': sympy.E, 't': sympy.Symbol('t'),
            'sin': sympy.sin, 'cos': sympy.cos, 'tan': sympy.tan,
            'exp': sympy.exp, 'log': sympy.log, 'sqrt': sympy.sqrt,
            'tanh': sympy.tanh, 'cosh': sympy.cosh, 'sinh': sympy.sinh,
            'abs': sympy.Abs,
        }
        for v in variables:
            local_dict[v] = sympy.Symbol(v)
        for p in params:
            local_dict[p] = sympy.Symbol(p)

        for i, eq in enumerate(equations):
            if not isinstance(eq, str):
                errors.append(f"equations[{i}]: must be a string (got {type(eq).__name__})")
                continue
            try:
                parse_expr(eq, local_dict=local_dict, transformations=transformations)
            except Exception as e:
                errors.append(f"equations[{i}]: failed to parse '{eq}' — {e}")

    except ImportError:
        errors.append("SymPy not installed — cannot validate equation syntax")

    # Initial conditions
    ic = spec.get("initial_conditions")
    if ic is None:
        errors.append("Missing 'initial_conditions' field")
    elif not isinstance(ic, list):
        errors.append("'initial_conditions' must be a list of numbers")
    elif len(ic) != len(variables):
        errors.append(
            f"'initial_conditions' length ({len(ic)}) must match "
            f"number of variables ({len(variables)})"
        )

    return errors


def _validate_pde_spec(spec: dict, params: dict) -> List[str]:
    """Validate PDE-specific fields."""
    errors = []

    # Equation
    equation = spec.get("equation")
    if equation is None:
        errors.append("Missing 'equation' field — RHS expression for u_t = ...")
    elif not isinstance(equation, str):
        errors.append("'equation' must be a string")
    else:
        try:
            import sympy
            from sympy.parsing.sympy_parser import (
                parse_expr, standard_transformations,
                implicit_multiplication_application,
            )
            transformations = standard_transformations + (implicit_multiplication_application,)
            local_dict = {
                'u': sympy.Symbol('u'), 'x': sympy.Symbol('x'),
                't': sympy.Symbol('t'),
                'u_x': sympy.Symbol('u_x'), 'u_xx': sympy.Symbol('u_xx'),
                'pi': sympy.pi, 'e': sympy.E,
                'sin': sympy.sin, 'cos': sympy.cos, 'exp': sympy.exp,
                'sqrt': sympy.sqrt, 'tanh': sympy.tanh, 'abs': sympy.Abs,
            }
            for p in params:
                local_dict[p] = sympy.Symbol(p)
            parse_expr(equation, local_dict=local_dict, transformations=transformations)
        except ImportError:
            errors.append("SymPy not installed — cannot validate equation syntax")
        except Exception as e:
            errors.append(f"'equation': failed to parse '{equation}' — {e}")

    # x_domain
    x_domain = spec.get("x_domain")
    if x_domain is None:
        errors.append("Missing 'x_domain' field (e.g., [0, 1])")
    elif not isinstance(x_domain, list) or len(x_domain) != 2:
        errors.append("'x_domain' must be a 2-element list [x_min, x_max]")

    # initial_condition
    ic = spec.get("initial_condition")
    if ic is None:
        errors.append("Missing 'initial_condition' field — expression in x")
    elif not isinstance(ic, str):
        errors.append("'initial_condition' must be a string expression in x")
    else:
        try:
            import sympy
            from sympy.parsing.sympy_parser import (
                parse_expr, standard_transformations,
                implicit_multiplication_application,
            )
            transformations = standard_transformations + (implicit_multiplication_application,)
            local_dict = {
                'x': sympy.Symbol('x'), 'pi': sympy.pi, 'e': sympy.E,
                'sin': sympy.sin, 'cos': sympy.cos, 'exp': sympy.exp,
                'sqrt': sympy.sqrt, 'tanh': sympy.tanh, 'abs': sympy.Abs,
            }
            for p in params:
                local_dict[p] = sympy.Symbol(p)
            parse_expr(ic, local_dict=local_dict, transformations=transformations)
        except ImportError:
            pass
        except Exception as e:
            errors.append(f"'initial_condition': failed to parse '{ic}' — {e}")

    # boundary_conditions
    bc = spec.get("boundary_conditions")
    if bc is None:
        errors.append("Missing 'boundary_conditions' field")
    elif not isinstance(bc, dict):
        errors.append("'boundary_conditions' must be a mapping")
    else:
        valid_types = ("dirichlet", "neumann")
        for side in ("left", "right"):
            if side not in bc:
                errors.append(f"'boundary_conditions' missing '{side}' entry")
            else:
                entry = bc[side]
                if not isinstance(entry, list) or len(entry) != 2:
                    errors.append(
                        f"boundary_conditions.{side} must be [type, value] "
                        f"(e.g., [dirichlet, 0])"
                    )
                elif entry[0] not in valid_types:
                    errors.append(
                        f"boundary_conditions.{side}[0] must be one of "
                        f"{valid_types} (got '{entry[0]}')"
                    )

    return errors


# ---- CLI entry point ----

def main():
    """Validate a problem YAML or network JSON file from the command line."""
    if len(sys.argv) < 2:
        print("Usage: python -m chemical_reservoir.validate <file> [<file2> ...]")
        print()
        print("Validates problem YAML files and reaction network JSON files.")
        print("Exit code 0 = all valid, 1 = errors found.")
        sys.exit(1)

    all_ok = True
    for fpath in sys.argv[1:]:
        path = Path(fpath)
        print(f"\n{'='*60}")
        print(f"Validating: {path}")
        print(f"{'='*60}")

        if path.suffix in (".yaml", ".yml"):
            errors = validate_problem_file(str(path))
            file_type = "problem"
        elif path.suffix == ".json":
            errors = validate_network_file(str(path))
            file_type = "network"
        else:
            print(f"  Unknown file type: {path.suffix}")
            print(f"  Expected .yaml/.yml (problem) or .json (network)")
            all_ok = False
            continue

        if errors:
            all_ok = False
            print(f"  INVALID — {len(errors)} error(s) found:")
            for err in errors:
                print(f"    ✗ {err}")
        else:
            print(f"  ✓ Valid {file_type} file")

    print()
    if all_ok:
        print("All files valid.")
        sys.exit(0)
    else:
        print("Validation failed — see errors above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
