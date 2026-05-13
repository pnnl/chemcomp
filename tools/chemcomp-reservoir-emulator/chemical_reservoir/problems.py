"""
ODE and PDE problem definitions.
"""

import numpy as np
from scipy.integrate import odeint
from dataclasses import dataclass, field
from typing import Callable, List, Tuple, Dict
import warnings


@dataclass
class ODEProblem:
    """Specification for an ODE problem."""
    func: Callable  # dy/dt = func(y, t, *params)
    params: tuple
    y0: np.ndarray
    t_span: Tuple[float, float]
    variable_names: List[str] = field(default_factory=list)
    name: str = "ODE"
    
    def __post_init__(self):
        self.y0 = np.atleast_1d(np.array(self.y0, dtype=float))
        if not self.variable_names:
            self.variable_names = [f"y{i}" for i in range(len(self.y0))]
    
    def solve_numerical(self, t_eval: np.ndarray) -> np.ndarray:
        """Get ground truth solution."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return odeint(self.func, self.y0, t_eval, args=self.params)
    
    def get_train_test_times(self, n_train: int, n_test: int, 
                            test_horizon_factor: float = 0.5) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate train and test time arrays.
        
        Args:
            n_train: Number of training points
            n_test: Number of test points
            test_horizon_factor: Fraction of training span to use for testing
            
        Returns:
            t_train, t_test, t_all
        """
        t_start = self.t_span[0]
        t_end = self.t_span[1]
        train_span = t_end - t_start
        
        t_train = np.linspace(t_start, t_end, n_train)
        
        dt = train_span / n_train
        t_test_start = t_end + dt
        t_test_end = t_end + test_horizon_factor * train_span
        t_test = np.linspace(t_test_start, t_test_end, n_test)
        
        t_all = np.concatenate([t_train, t_test])
        
        return t_train, t_test, t_all


@dataclass
class PDEProblem1D:
    """Specification for a 1D PDE problem (Method of Lines)."""
    pde_func: Callable  # u_t = f(u, x, t, u_x, u_xx, *params)
    params: tuple
    x_domain: Tuple[float, float]
    t_span: Tuple[float, float]
    initial_condition: Callable  # u(x, 0) = ic(x)
    boundary_conditions: Dict  # {'left': ('dirichlet', val), 'right': ...}
    n_spatial: int = 50
    name: str = "PDE"
    
    def solve_numerical(self, t_eval: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Get ground truth solution using method of lines."""
        x = np.linspace(self.x_domain[0], self.x_domain[1], self.n_spatial)
        dx = x[1] - x[0]
        
        def ode_system(u, t):
            dudt = np.zeros_like(u)
            u_x = np.zeros_like(u)
            u_xx = np.zeros_like(u)
            
            # Central differences for interior
            u_x[1:-1] = (u[2:] - u[:-2]) / (2 * dx)
            u_xx[1:-1] = (u[2:] - 2*u[1:-1] + u[:-2]) / (dx**2)
            
            # Boundary handling
            bc_l = self.boundary_conditions.get('left', ('dirichlet', 0))
            bc_r = self.boundary_conditions.get('right', ('dirichlet', 0))
            
            if bc_l[0] == 'dirichlet':
                u_xx[0] = 2 * (u[1] - u[0]) / (dx**2)
            if bc_r[0] == 'dirichlet':
                u_xx[-1] = 2 * (u[-2] - u[-1]) / (dx**2)
            
            for i in range(len(u)):
                dudt[i] = self.pde_func(u[i], x[i], t, u_x[i], u_xx[i], *self.params)
            
            if bc_l[0] == 'dirichlet':
                dudt[0] = 0
            if bc_r[0] == 'dirichlet':
                dudt[-1] = 0
            
            return dudt
        
        u0 = self.initial_condition(x)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = odeint(ode_system, u0, t_eval)
        return x, sol
    
    def get_train_test_times(self, n_train: int, n_test: int,
                            test_horizon_factor: float = 0.3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Generate train and test time arrays."""
        t_start = self.t_span[0]
        t_end = self.t_span[1]
        train_span = t_end - t_start
        
        t_train = np.linspace(t_start, t_end, n_train)
        
        dt = train_span / n_train
        t_test_start = t_end + dt
        t_test_end = t_end + test_horizon_factor * train_span
        t_test = np.linspace(t_test_start, t_test_end, n_test)
        
        t_all = np.concatenate([t_train, t_test])
        
        return t_train, t_test, t_all


# ============================================================================
# PREDEFINED TEST PROBLEMS
# ============================================================================

def create_schnakenberg_problem(t_end: float = 15.0) -> ODEProblem:
    """Schnakenberg oscillator ODE."""
    def schnakenberg(y, t, a, b, c, d):
        x, z = y
        return [d*z*x**2 - b*x + a, -d*z*x**2 + c]
    
    return ODEProblem(
        func=schnakenberg,
        params=(0.5, 2.2, 1.5, 1.0),
        y0=[2.5, 1.1],
        t_span=(0, t_end),
        variable_names=['x', 'y'],
        name="Schnakenberg"
    )


def create_lorenz_problem(t_end: float = 10.0) -> ODEProblem:
    """Lorenz chaotic system."""
    def lorenz(y, t, sigma, rho, beta):
        x, z, w = y
        return [
            sigma * (z - x),
            x * (rho - w) - z,
            x * z - beta * w
        ]
    
    return ODEProblem(
        func=lorenz,
        params=(10.0, 28.0, 8/3),
        y0=[1.0, 1.0, 1.0],
        t_span=(0, t_end),
        variable_names=['x', 'y', 'z'],
        name="Lorenz"
    )


def create_vanderpol_problem(mu: float = 1.5, t_end: float = 20.0) -> ODEProblem:
    """Van der Pol oscillator."""
    def vanderpol(y, t, mu):
        x, v = y
        return [v, mu*(1 - x**2)*v - x]
    
    return ODEProblem(
        func=vanderpol,
        params=(mu,),
        y0=[2.0, 0.0],
        t_span=(0, t_end),
        variable_names=['x', 'v'],
        name="VanDerPol"
    )


def create_lotka_volterra_problem(t_end: float = 15.0) -> ODEProblem:
    """Lotka-Volterra predator-prey model."""
    def lotka_volterra(y, t, alpha, beta, gamma, delta):
        prey, predator = y
        return [
            alpha * prey - beta * prey * predator,
            delta * prey * predator - gamma * predator
        ]
    
    return ODEProblem(
        func=lotka_volterra,
        params=(1.1, 0.4, 0.4, 0.1),
        y0=[10.0, 5.0],
        t_span=(0, t_end),
        variable_names=['prey', 'predator'],
        name="LotkaVolterra"
    )


def create_brusselator_problem(t_end: float = 20.0) -> ODEProblem:
    """Brusselator chemical oscillator."""
    def brusselator(y, t, a, b):
        x, z = y
        return [
            a + x**2 * z - (b + 1) * x,
            b * x - x**2 * z
        ]
    
    return ODEProblem(
        func=brusselator,
        params=(1.0, 3.0),
        y0=[1.0, 1.0],
        t_span=(0, t_end),
        variable_names=['x', 'y'],
        name="Brusselator"
    )


def create_heat_equation_problem(t_end: float = 2.0, alpha: float = 0.1) -> PDEProblem1D:
    """1D Heat equation."""
    def heat_pde(u, x, t, u_x, u_xx, alpha):
        return alpha * u_xx
    
    return PDEProblem1D(
        pde_func=heat_pde,
        params=(alpha,),
        x_domain=(0, 1),
        t_span=(0, t_end),
        initial_condition=lambda x: np.sin(np.pi * x),
        boundary_conditions={'left': ('dirichlet', 0), 'right': ('dirichlet', 0)},
        n_spatial=30,
        name="HeatEquation"
    )


def create_advection_diffusion_problem(t_end: float = 1.0) -> PDEProblem1D:
    """1D Advection-diffusion equation."""
    def advection_diffusion(u, x, t, u_x, u_xx, v, D):
        return -v * u_x + D * u_xx
    
    return PDEProblem1D(
        pde_func=advection_diffusion,
        params=(0.5, 0.05),  # velocity, diffusion
        x_domain=(0, 1),
        t_span=(0, t_end),
        initial_condition=lambda x: np.exp(-50*(x - 0.3)**2),
        boundary_conditions={'left': ('dirichlet', 0), 'right': ('dirichlet', 0)},
        n_spatial=50,
        name="AdvectionDiffusion"
    )


def create_fisher_equation_problem(t_end: float = 5.0) -> PDEProblem1D:
    """Fisher's equation (reaction-diffusion)."""
    def fisher(u, x, t, u_x, u_xx, D, r):
        return D * u_xx + r * u * (1 - u)
    
    return PDEProblem1D(
        pde_func=fisher,
        params=(0.1, 1.0),  # diffusion, growth rate
        x_domain=(0, 1),
        t_span=(0, t_end),
        initial_condition=lambda x: np.where(x < 0.3, 1.0, 0.0),
        boundary_conditions={'left': ('neumann', 0), 'right': ('neumann', 0)},
        n_spatial=40,
        name="FisherEquation"
    )


def create_burgers_equation_problem(t_end: float = 0.5) -> PDEProblem1D:
    """Burgers' equation."""
    def burgers(u, x, t, u_x, u_xx, nu):
        return -u * u_x + nu * u_xx
    
    return PDEProblem1D(
        pde_func=burgers,
        params=(0.01,),  # viscosity
        x_domain=(0, 1),
        t_span=(0, t_end),
        initial_condition=lambda x: np.sin(2 * np.pi * x),
        boundary_conditions={'left': ('dirichlet', 0), 'right': ('dirichlet', 0)},
        n_spatial=50,
        name="BurgersEquation"
    )


# Dictionary of all available problems
ODE_PROBLEMS = {
    'schnakenberg': create_schnakenberg_problem,
    'lorenz': create_lorenz_problem,
    'vanderpol': create_vanderpol_problem,
    'lotka_volterra': create_lotka_volterra_problem,
    'brusselator': create_brusselator_problem,
}

PDE_PROBLEMS = {
    'heat': create_heat_equation_problem,
    'advection_diffusion': create_advection_diffusion_problem,
    'fisher': create_fisher_equation_problem,
    'burgers': create_burgers_equation_problem,
}