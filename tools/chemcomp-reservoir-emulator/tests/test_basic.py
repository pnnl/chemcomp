"""Basic smoke tests for chemical_reservoir package."""

import numpy as np
import pytest

import chemical_reservoir as cr
from chemical_reservoir import io, networks, problems


class TestImports:
    def test_version(self):
        assert cr.__version__ == "1.0.0"

    def test_core_classes(self):
        assert cr.ConstrainedChemicalReservoir is not None
        assert cr.ODEProblem is not None
        assert cr.PDEProblem1D is not None
        assert cr.BayesianOptimizer is not None


class TestNetworks:
    def test_list_builtin(self):
        names = networks.list_builtin_networks()
        assert len(names) == 8
        assert "cascade" in names

    def test_get_io_config(self):
        inp, out = networks.get_network_io_config("cascade")
        assert isinstance(inp, list)
        assert isinstance(out, list)


class TestProblems:
    def test_vanderpol(self):
        p = problems.create_vanderpol_problem()
        assert p.name == "VanDerPol"
        t = np.linspace(*p.t_span, 100)
        y = p.solve_numerical(t)
        assert y.shape == (100, 2)

    def test_heat_equation(self):
        p = problems.create_heat_equation_problem()
        assert p.name == "HeatEquation"


class TestIO:
    def test_load_network(self, tmp_path):
        cr.io.save_reaction_network(
            filepath=tmp_path / "test.json",
            reactions=[
                {"species": ["A", "B"], "reactants": [1, 0], "products": [0, 1], "rate": 1.0},
            ],
            input_channels=["A"],
            output_channels=["A", "B"],
            name="test_net",
        )
        config = cr.io.load_reaction_network(tmp_path / "test.json")
        assert config["name"] == "test_net"
        assert len(config["reactions"]) == 1

    def test_load_reservoir_from_json(self):
        from pathlib import Path
        json_path = Path(__file__).parent.parent / "network_configs" / "minimal.json"
        if json_path.exists():
            reservoir, config = cr.io.load_reservoir_from_json(json_path)
            assert config["name"] == "minimal"


class TestReservoir:
    def test_solve_ode(self):
        from pathlib import Path
        json_path = Path(__file__).parent.parent / "network_configs" / "minimal.json"
        if not json_path.exists():
            pytest.skip("minimal.json not found")

        reservoir, _ = cr.io.load_reservoir_from_json(json_path)
        problem = problems.create_vanderpol_problem(t_end=5.0)
        results = reservoir.solve_ode(problem, n_train=100, n_test=50)

        assert "y_true" in results
        assert "y_pred" in results
        assert results["y_true"].shape[0] > 0
