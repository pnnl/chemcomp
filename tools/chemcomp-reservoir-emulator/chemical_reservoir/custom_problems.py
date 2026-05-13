"""
Load custom ODE/PDE problems from YAML files using SymPy for symbolic parsing.

YAML format for an ODE:

    type: ode
    name: MyOscillator
    variables: [x, v]
    parameters:
      mu: 1.5
      omega: 1.0
    equations:
      - "v"
      - "mu*(1 - x**2)*v - omega**2*x"
    initial_conditions: [2.0, 0.0]
    t_span: [0, 20]

YAML format for a PDE:

    type: pde
    name: CustomDiffusion
    parameters:
      D: 0.1
      k: 1.0
    equation: "D*u_xx - k*u"
    x_domain: [0, 1]
    t_span: [0, 2]
    n_spatial: 50
    initial_condition: "sin(pi*x)"
    boundary_conditions:
      left: [dirichlet, 0]
      right: [dirichlet, 0]
"""

import numpy as np
import yaml
from pathlib import Path
from typing import Union

from .problems import ODEProblem, PDEProblem1D

try:
    import sympy
    from sympy import lambdify, symbols, sympify, pi, sin, cos, exp, sqrt, tanh
    from sympy.parsing.sympy_parser import (
        parse_expr,
        standard_transformations,
        implicit_multiplication_application,
    )
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False


def _require_sympy():
    if not HAS_SYMPY:
        raise ImportError(
            "SymPy is required to load custom problem files. "
            "Install it with: pip install sympy"
        )


def _parse_ode_expr(expr_str: str, var_names: list, param_names: list):
    """Parse a single ODE RHS expression string into a callable."""
    _require_sympy()

    transformations = standard_transformations + (implicit_multiplication_application,)

    # Build the local namespace with known functions and constants
    local_dict = {
        'pi': sympy.pi, 'e': sympy.E,
        'sin': sympy.sin, 'cos': sympy.cos, 'tan': sympy.tan,
        'exp': sympy.exp, 'log': sympy.log, 'sqrt': sympy.sqrt,
        'tanh': sympy.tanh, 'cosh': sympy.cosh, 'sinh': sympy.sinh,
        'abs': sympy.Abs, 'Abs': sympy.Abs,
    }

    # Create symbols for variables and parameters
    var_syms = tuple(symbols(var_names)) if var_names else ()
    param_syms = tuple(symbols(param_names)) if param_names else ()

    for v in var_syms:
        local_dict[str(v)] = v
    for p in param_syms:
        local_dict[str(p)] = p

    # Add a time symbol
    t_sym = symbols('t')
    local_dict['t'] = t_sym

    expr = parse_expr(expr_str, local_dict=local_dict, transformations=transformations)
    return expr, var_syms, param_syms, t_sym


def _parse_pde_expr(expr_str: str, param_names: list):
    """Parse a PDE RHS expression (function of u, x, t, u_x, u_xx, params)."""
    _require_sympy()

    transformations = standard_transformations + (implicit_multiplication_application,)

    local_dict = {
        'pi': sympy.pi, 'e': sympy.E,
        'sin': sympy.sin, 'cos': sympy.cos, 'tan': sympy.tan,
        'exp': sympy.exp, 'log': sympy.log, 'sqrt': sympy.sqrt,
        'tanh': sympy.tanh, 'cosh': sympy.cosh, 'sinh': sympy.sinh,
        'abs': sympy.Abs, 'Abs': sympy.Abs,
    }

    # PDE built-in symbols
    u_sym, x_sym, t_sym = symbols('u x t')
    u_x_sym, u_xx_sym = symbols('u_x u_xx')
    local_dict.update({
        'u': u_sym, 'x': x_sym, 't': t_sym,
        'u_x': u_x_sym, 'u_xx': u_xx_sym,
    })

    param_syms = tuple(symbols(param_names)) if param_names else ()
    for p in param_syms:
        local_dict[str(p)] = p

    expr = parse_expr(expr_str, local_dict=local_dict, transformations=transformations)
    return expr, u_sym, x_sym, t_sym, u_x_sym, u_xx_sym, param_syms


def load_problem_from_yaml(filepath: Union[str, Path]) -> Union[ODEProblem, PDEProblem1D]:
    """
    Load a custom ODE or PDE problem from a YAML file.

    The YAML file uses SymPy-compatible expressions for the equations.
    See module docstring for format specification.

    Args:
        filepath: Path to a YAML problem definition file.

    Returns:
        An ODEProblem or PDEProblem1D instance ready for use with the reservoir.

    Raises:
        ImportError: If sympy is not installed.
        ValueError: If the YAML format is invalid.
    """
    _require_sympy()

    filepath = Path(filepath)
    with open(filepath, 'r') as f:
        spec = yaml.safe_load(f)

    problem_type = spec.get('type', 'ode').lower()

    if problem_type == 'ode':
        return _build_ode_problem(spec)
    elif problem_type == 'pde':
        return _build_pde_problem(spec)
    else:
        raise ValueError(f"Unknown problem type: '{problem_type}'. Must be 'ode' or 'pde'.")


def _build_ode_problem(spec: dict) -> ODEProblem:
    """Build an ODEProblem from a parsed YAML spec."""
    var_names = spec['variables']
    params_dict = spec.get('parameters', {})
    param_names = list(params_dict.keys())
    param_values = tuple(float(v) for v in params_dict.values())
    equations = spec['equations']

    if len(equations) != len(var_names):
        raise ValueError(
            f"Number of equations ({len(equations)}) must match "
            f"number of variables ({len(var_names)})"
        )

    # Parse all expressions and build a single vectorised function
    all_exprs = []
    for eq_str in equations:
        expr, var_syms, param_syms, t_sym = _parse_ode_expr(eq_str, var_names, param_names)
        all_exprs.append(expr)

    # Build a single lambdified function: f(y, t, *params) -> list
    # odeint expects func(y, t, *args)
    all_symbols = list(var_syms) + [t_sym] + list(param_syms)
    funcs = [lambdify(all_symbols, expr, modules='numpy') for expr in all_exprs]

    def ode_func(y, t, *params):
        # Unpack variables
        args = list(y) + [t] + list(params)
        return [f(*args) for f in funcs]

    y0 = np.array(spec['initial_conditions'], dtype=float)
    t_span = tuple(spec['t_span'])

    return ODEProblem(
        func=ode_func,
        params=param_values,
        y0=y0,
        t_span=t_span,
        variable_names=var_names,
        name=spec.get('name', 'CustomODE'),
    )


def _build_pde_problem(spec: dict) -> PDEProblem1D:
    """Build a PDEProblem1D from a parsed YAML spec."""
    params_dict = spec.get('parameters', {})
    param_names = list(params_dict.keys())
    param_values = tuple(float(v) for v in params_dict.values())

    # Parse main PDE equation
    expr, u_sym, x_sym, t_sym, u_x_sym, u_xx_sym, param_syms = _parse_pde_expr(
        spec['equation'], param_names
    )

    pde_symbols = [u_sym, x_sym, t_sym, u_x_sym, u_xx_sym] + list(param_syms)
    pde_lambda = lambdify(pde_symbols, expr, modules='numpy')

    def pde_func(u, x, t, u_x, u_xx, *params):
        return pde_lambda(u, x, t, u_x, u_xx, *params)

    # Parse initial condition expression
    ic_expr_str = spec['initial_condition']
    transformations = standard_transformations + (implicit_multiplication_application,)
    local_dict = {
        'pi': sympy.pi, 'e': sympy.E,
        'sin': sympy.sin, 'cos': sympy.cos, 'exp': sympy.exp,
        'sqrt': sympy.sqrt, 'tanh': sympy.tanh,
        'abs': sympy.Abs, 'x': x_sym,
    }
    # Also allow parameters in IC
    for p in param_syms:
        local_dict[str(p)] = p

    ic_expr = parse_expr(ic_expr_str, local_dict=local_dict, transformations=transformations)
    ic_func = lambdify([x_sym] + list(param_syms), ic_expr, modules='numpy')

    def initial_condition(x):
        return ic_func(x, *param_values)

    # Boundary conditions
    bc_spec = spec.get('boundary_conditions', {})
    boundary_conditions = {}
    for side in ('left', 'right'):
        if side in bc_spec:
            bc = bc_spec[side]
            boundary_conditions[side] = (bc[0], float(bc[1]))

    return PDEProblem1D(
        pde_func=pde_func,
        params=param_values,
        x_domain=tuple(spec['x_domain']),
        t_span=tuple(spec['t_span']),
        initial_condition=initial_condition,
        boundary_conditions=boundary_conditions,
        n_spatial=spec.get('n_spatial', 50),
        name=spec.get('name', 'CustomPDE'),
    )
