"""Tests for the validation module and custom problem loading."""

import json
import numpy as np
import pytest
from pathlib import Path

import chemical_reservoir as cr
from chemical_reservoir.validate import validate_problem_file, validate_network_file
from chemical_reservoir.custom_problems import load_problem_from_yaml


PROJECT_ROOT = Path(__file__).parent.parent


# ============================================================================
# Network validation tests
# ============================================================================

class TestNetworkValidation:
    """Test validate_network_file against valid and invalid inputs."""

    def test_valid_minimal_network(self):
        errors = validate_network_file(str(PROJECT_ROOT / "network_configs" / "minimal.json"))
        assert errors == [], f"Unexpected errors: {errors}"

    def test_valid_cascade_network(self):
        errors = validate_network_file(str(PROJECT_ROOT / "network_configs" / "cascade.json"))
        assert errors == []

    def test_valid_highdim_network(self):
        errors = validate_network_file(
            str(PROJECT_ROOT / "network_configs" / "high_dimensional_random.json")
        )
        assert errors == []

    def test_valid_multisubnet_network(self):
        path = PROJECT_ROOT / "network_configs" / "multi_subnet_leaky.json"
        if path.exists():
            errors = validate_network_file(str(path))
            assert errors == []

    def test_all_builtin_networks(self):
        """Every JSON in network_configs/ should validate cleanly."""
        config_dir = PROJECT_ROOT / "network_configs"
        for f in sorted(config_dir.glob("*.json")):
            errors = validate_network_file(str(f))
            assert errors == [], f"{f.name} has errors: {errors}"

    def test_missing_file(self):
        errors = validate_network_file("/nonexistent/path.json")
        assert len(errors) == 1
        assert "not found" in errors[0].lower()

    def test_invalid_json(self, tmp_path):
        bad = tmp_path / "bad.json"
        bad.write_text("{ not valid json")
        errors = validate_network_file(str(bad))
        assert any("invalid json" in e.lower() for e in errors)

    def test_missing_reactions_key(self, tmp_path):
        f = tmp_path / "no_rxn.json"
        f.write_text(json.dumps({"name": "test"}))
        errors = validate_network_file(str(f))
        assert any("reactions" in e.lower() for e in errors)

    def test_empty_reactions(self, tmp_path):
        f = tmp_path / "empty_rxn.json"
        f.write_text(json.dumps({"reactions": []}))
        errors = validate_network_file(str(f))
        assert any("non-empty" in e.lower() for e in errors)

    def test_reaction_missing_rate(self, tmp_path):
        f = tmp_path / "no_rate.json"
        f.write_text(json.dumps({
            "reactions": [
                {"species": ["A", "B"], "reactants": [1, 0], "products": [0, 1]}
            ],
            "input_channels": ["A"],
            "output_channels": ["B"],
        }))
        errors = validate_network_file(str(f))
        assert any("rate" in e.lower() for e in errors)

    def test_reaction_negative_rate(self, tmp_path):
        f = tmp_path / "neg_rate.json"
        f.write_text(json.dumps({
            "reactions": [
                {"species": ["A", "B"], "reactants": [1, 0], "products": [0, 1], "rate": -1.0}
            ],
            "input_channels": ["A"],
            "output_channels": ["B"],
        }))
        errors = validate_network_file(str(f))
        assert any("positive" in e.lower() for e in errors)

    def test_stoichiometry_length_mismatch(self, tmp_path):
        f = tmp_path / "mismatch.json"
        f.write_text(json.dumps({
            "reactions": [
                {"species": ["A", "B", "C"], "reactants": [1, 0], "products": [0, 1, 0], "rate": 1.0}
            ],
            "input_channels": ["A"],
            "output_channels": ["B"],
        }))
        errors = validate_network_file(str(f))
        assert any("length" in e.lower() for e in errors)

    def test_invalid_channel_species(self, tmp_path):
        f = tmp_path / "bad_channel.json"
        f.write_text(json.dumps({
            "reactions": [
                {"species": ["A", "B"], "reactants": [1, 0], "products": [0, 1], "rate": 1.0}
            ],
            "input_channels": ["X"],
            "output_channels": ["A"],
        }))
        errors = validate_network_file(str(f))
        assert any("not found" in e.lower() for e in errors)


# ============================================================================
# Problem file validation tests
# ============================================================================

class TestProblemValidation:
    """Test validate_problem_file against valid and invalid inputs."""

    def test_valid_vanderpol(self):
        errors = validate_problem_file(str(PROJECT_ROOT / "problem_configs" / "vanderpol.yaml"))
        assert errors == [], f"Unexpected errors: {errors}"

    def test_valid_lorenz(self):
        errors = validate_problem_file(str(PROJECT_ROOT / "problem_configs" / "lorenz.yaml"))
        assert errors == []

    def test_valid_lotka_volterra(self):
        errors = validate_problem_file(
            str(PROJECT_ROOT / "problem_configs" / "lotka_volterra.yaml")
        )
        assert errors == []

    def test_valid_heat_equation(self):
        errors = validate_problem_file(
            str(PROJECT_ROOT / "problem_configs" / "heat_equation.yaml")
        )
        assert errors == []

    def test_valid_burgers(self):
        errors = validate_problem_file(str(PROJECT_ROOT / "problem_configs" / "burgers.yaml"))
        assert errors == []

    def test_valid_fisher(self):
        errors = validate_problem_file(str(PROJECT_ROOT / "problem_configs" / "fisher.yaml"))
        assert errors == []

    def test_all_problem_configs(self):
        """Every YAML in problem_configs/ should validate cleanly."""
        config_dir = PROJECT_ROOT / "problem_configs"
        for f in sorted(config_dir.glob("*.yaml")):
            errors = validate_problem_file(str(f))
            assert errors == [], f"{f.name} has errors: {errors}"

    def test_missing_file(self):
        errors = validate_problem_file("/nonexistent/path.yaml")
        assert len(errors) == 1
        assert "not found" in errors[0].lower()

    def test_invalid_yaml(self, tmp_path):
        bad = tmp_path / "bad.yaml"
        bad.write_text("{{{{invalid yaml")
        errors = validate_problem_file(str(bad))
        assert len(errors) >= 1

    def test_missing_type(self, tmp_path):
        f = tmp_path / "no_type.yaml"
        f.write_text("name: Test\nvariables: [x]\n")
        errors = validate_problem_file(str(f))
        assert any("type" in e.lower() for e in errors)

    def test_ode_missing_equations(self, tmp_path):
        f = tmp_path / "no_eq.yaml"
        f.write_text(
            "type: ode\nname: Test\nvariables: [x]\n"
            "initial_conditions: [1.0]\nt_span: [0, 10]\n"
        )
        errors = validate_problem_file(str(f))
        assert any("equations" in e.lower() for e in errors)

    def test_ode_equation_count_mismatch(self, tmp_path):
        f = tmp_path / "mismatch.yaml"
        f.write_text(
            "type: ode\nname: Test\nvariables: [x, y]\n"
            "equations:\n  - 'y'\n"
            "initial_conditions: [1.0, 0.0]\nt_span: [0, 10]\n"
        )
        errors = validate_problem_file(str(f))
        assert any("number of equations" in e.lower() for e in errors)

    def test_ode_ic_length_mismatch(self, tmp_path):
        f = tmp_path / "bad_ic.yaml"
        f.write_text(
            "type: ode\nname: Test\nvariables: [x, y]\n"
            "equations:\n  - 'y'\n  - '-x'\n"
            "initial_conditions: [1.0]\nt_span: [0, 10]\n"
        )
        errors = validate_problem_file(str(f))
        assert any("initial_conditions" in e.lower() for e in errors)

    def test_ode_bad_expression(self, tmp_path):
        f = tmp_path / "bad_expr.yaml"
        f.write_text(
            "type: ode\nname: Test\nvariables: [x]\n"
            "equations:\n  - '***invalid***'\n"
            "initial_conditions: [1.0]\nt_span: [0, 10]\n"
        )
        errors = validate_problem_file(str(f))
        assert any("parse" in e.lower() for e in errors)

    def test_ode_bad_tspan(self, tmp_path):
        f = tmp_path / "bad_tspan.yaml"
        f.write_text(
            "type: ode\nname: Test\nvariables: [x]\n"
            "equations:\n  - '-x'\n"
            "initial_conditions: [1.0]\nt_span: [10, 0]\n"
        )
        errors = validate_problem_file(str(f))
        assert any("t_start" in e.lower() or "t_span" in e.lower() for e in errors)

    def test_pde_missing_equation(self, tmp_path):
        f = tmp_path / "no_eq.yaml"
        f.write_text(
            "type: pde\nname: Test\nvariables: [u]\n"
            "x_domain: [0, 1]\nt_span: [0, 1]\n"
            "initial_condition: 'sin(pi*x)'\n"
            "boundary_conditions:\n  left: [dirichlet, 0]\n  right: [dirichlet, 0]\n"
        )
        errors = validate_problem_file(str(f))
        assert any("equation" in e.lower() for e in errors)

    def test_pde_missing_bc(self, tmp_path):
        f = tmp_path / "no_bc.yaml"
        f.write_text(
            "type: pde\nname: Test\nvariables: [u]\n"
            "parameters:\n  D: 0.1\n"
            "equation: 'D*u_xx'\n"
            "x_domain: [0, 1]\nt_span: [0, 1]\n"
            "initial_condition: 'sin(pi*x)'\n"
        )
        errors = validate_problem_file(str(f))
        assert any("boundary" in e.lower() for e in errors)

    def test_pde_invalid_bc_type(self, tmp_path):
        f = tmp_path / "bad_bc.yaml"
        f.write_text(
            "type: pde\nname: Test\nvariables: [u]\n"
            "parameters:\n  D: 0.1\n"
            "equation: 'D*u_xx'\n"
            "x_domain: [0, 1]\nt_span: [0, 1]\n"
            "initial_condition: 'sin(pi*x)'\n"
            "boundary_conditions:\n  left: [periodic, 0]\n  right: [dirichlet, 0]\n"
        )
        errors = validate_problem_file(str(f))
        assert any("dirichlet" in e.lower() or "neumann" in e.lower() for e in errors)


# ============================================================================
# Custom problem loading tests (functional)
# ============================================================================

class TestCustomProblemLoading:
    """Test that loaded custom problems produce correct numerical results."""

    def test_vanderpol_matches_builtin(self):
        custom = load_problem_from_yaml(str(PROJECT_ROOT / "problem_configs" / "vanderpol.yaml"))
        builtin = cr.problems.create_vanderpol_problem()

        t = np.linspace(0, 5, 50)
        sol_custom = custom.solve_numerical(t)
        sol_builtin = builtin.solve_numerical(t)
        np.testing.assert_allclose(sol_custom, sol_builtin, rtol=1e-4)

    def test_lorenz_integration(self):
        problem = load_problem_from_yaml(str(PROJECT_ROOT / "problem_configs" / "lorenz.yaml"))
        t = np.linspace(0, 1, 50)
        sol = problem.solve_numerical(t)
        assert sol.shape == (50, 3)
        assert not np.any(np.isnan(sol))

    def test_lotka_volterra_matches_builtin(self):
        custom = load_problem_from_yaml(
            str(PROJECT_ROOT / "problem_configs" / "lotka_volterra.yaml")
        )
        builtin = cr.problems.create_lotka_volterra_problem()

        t = np.linspace(0, 5, 50)
        sol_custom = custom.solve_numerical(t)
        sol_builtin = builtin.solve_numerical(t)
        np.testing.assert_allclose(sol_custom, sol_builtin, rtol=1e-4)

    def test_heat_equation_integration(self):
        problem = load_problem_from_yaml(
            str(PROJECT_ROOT / "problem_configs" / "heat_equation.yaml")
        )
        t = np.linspace(0, 0.5, 20)
        x, sol = problem.solve_numerical(t)
        assert sol.shape == (20, 30)
        assert not np.any(np.isnan(sol))
        # Heat equation should decay
        assert sol[-1].max() < sol[0].max()

    def test_burgers_integration(self):
        problem = load_problem_from_yaml(str(PROJECT_ROOT / "problem_configs" / "burgers.yaml"))
        t = np.linspace(0, 0.1, 10)
        x, sol = problem.solve_numerical(t)
        assert sol.shape == (10, 50)
        assert not np.any(np.isnan(sol))

    def test_fisher_integration(self):
        problem = load_problem_from_yaml(str(PROJECT_ROOT / "problem_configs" / "fisher.yaml"))
        t = np.linspace(0, 1, 10)
        x, sol = problem.solve_numerical(t)
        assert sol.shape == (10, 40)
        assert not np.any(np.isnan(sol))


# ============================================================================
# End-to-end reservoir tests with custom problems
# ============================================================================

class TestEndToEnd:
    """Test full reservoir pipeline with custom-loaded problems."""

    @pytest.fixture
    def reservoir(self):
        json_path = PROJECT_ROOT / "network_configs" / "minimal.json"
        if not json_path.exists():
            pytest.skip("minimal.json not found")
        r, _ = cr.io.load_reservoir_from_json(str(json_path))
        return r

    def test_custom_ode_open_loop(self, reservoir):
        problem = load_problem_from_yaml(str(PROJECT_ROOT / "problem_configs" / "vanderpol.yaml"))
        results = reservoir.solve_ode(problem, n_train=80, n_test=40)
        assert "y_pred" in results
        assert results["y_pred"].shape[0] == 120

    def test_custom_ode_closed_loop(self, reservoir):
        problem = load_problem_from_yaml(str(PROJECT_ROOT / "problem_configs" / "vanderpol.yaml"))
        results = reservoir.solve_ode_closed_loop(problem, n_train=80, n_test=40)
        assert "y_pred" in results
        assert results["y_pred"].shape[0] == 120

    def test_custom_pde(self, reservoir):
        problem = load_problem_from_yaml(
            str(PROJECT_ROOT / "problem_configs" / "heat_equation.yaml")
        )
        results = reservoir.solve_pde(problem, n_train_time=30, n_test_time=15)
        assert "train_rmse" in results
        assert results["train_rmse"] >= 0
