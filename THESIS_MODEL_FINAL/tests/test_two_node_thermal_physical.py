from __future__ import annotations

import inspect

import numpy as np
import pytest
from scipy.linalg import expm

from models.two_node_thermal_physical import (
    TwoNodeThermalParameters,
    simulate_two_node_physical,
)


def test_exact_nonuniform_update_matches_matrix_exponential() -> None:
    parameters = TwoNodeThermalParameters(80.0, 20.0, 0.4, 0.1)
    time_s = np.asarray([0.0, 0.13, 2.7, 12.0])
    heat_w = np.asarray([1.2, -0.2, 0.7, 0.7])
    ambient_k = np.asarray([298.0, 298.1, 297.9, 297.9])
    initial = np.asarray([299.0, 298.0])
    matrix = np.asarray([[-0.4 / 80.0, 0.4 / 80.0], [0.4 / 20.0, -0.5 / 20.0]])
    reference = [initial.copy()]
    state = initial.copy()
    for index, dt_s in enumerate(np.diff(time_s)):
        transition = expm(matrix * dt_s)
        forcing = np.asarray([heat_w[index] / 80.0, 0.1 * ambient_k[index] / 20.0])
        state = transition @ state + np.linalg.solve(matrix, (transition - np.eye(2))) @ forcing
        reference.append(state.copy())
    result = simulate_two_node_physical(
        time_s, heat_w, ambient_k, parameters, initial[0], initial[1]
    )
    np.testing.assert_allclose(
        np.column_stack((result.t_core_k, result.t_surface_k)),
        np.asarray(reference),
        rtol=0.0,
        atol=1.0e-11,
    )


def test_equilibrium_is_invariant() -> None:
    parameters = TwoNodeThermalParameters(90.0, 25.0, 0.5, 0.2)
    ambient_k = 298.0
    heat_w = 1.0
    surface_equilibrium = ambient_k + heat_w / parameters.hA_w_per_k
    core_equilibrium = surface_equilibrium + heat_w / parameters.g_cs_w_per_k
    result = simulate_two_node_physical(
        np.asarray([0.0, 0.2, 10.0, 1000.0]),
        np.full(4, heat_w),
        ambient_k,
        parameters,
        core_equilibrium,
        surface_equilibrium,
    )
    np.testing.assert_allclose(result.t_core_k, core_equilibrium, atol=1.0e-12)
    np.testing.assert_allclose(result.t_surface_k, surface_equilibrium, atol=1.0e-12)


@pytest.mark.parametrize("field", range(4))
def test_nonpositive_parameter_is_rejected(field: int) -> None:
    values = [80.0, 20.0, 0.4, 0.1]
    values[field] = 0.0
    with pytest.raises(ValueError, match="positive"):
        simulate_two_node_physical(
            np.asarray([0.0, 1.0]),
            np.zeros(2),
            298.0,
            TwoNodeThermalParameters(*values),
            298.0,
            298.0,
        )


def test_clean_api_has_no_empirical_correction_arguments() -> None:
    names = set(inspect.signature(simulate_two_node_physical).parameters)
    forbidden = {"soc", "current_a", "correction", "multiplier", "bias", "offset"}
    assert names.isdisjoint(forbidden)
