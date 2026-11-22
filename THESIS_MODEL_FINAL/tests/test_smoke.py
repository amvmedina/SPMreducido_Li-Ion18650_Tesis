"""Minimal execution tests; these do not perform identification or optimization."""

import numpy as np

from models import GroupedSPMParameters, PARAMETER_ORDER, simulate_grouped_spm


def smoke_parameters() -> GroupedSPMParameters:
    # Numerical fixture only; this is not a calibrated or synthetic-truth dataset.
    return GroupedSPMParameters(
        alpha_n_s=3000.0,
        alpha_p_s=1800.0,
        q_n_c=11000.0,
        q_p_c=11000.0,
        d_n=5.0e-5,
        d_p=7.0e-4,
        soc0_n=0.60,
        soc0_p=0.40,
        r0_ohm=0.02,
    )


def test_model_imports_and_executes_without_synthetic_machinery() -> None:
    time_s = np.linspace(0.0, 20.0, 21)
    current_a = np.full(time_s.size, 0.2)
    result = simulate_grouped_spm(time_s, current_a, smoke_parameters())

    assert result.domain_valid
    assert np.all(np.isfinite(result.voltage_v))
    assert result.x_n_avg[-1] < result.x_n_avg[0]
    assert result.x_p_avg[-1] > result.x_p_avg[0]


def test_parameter_vector_order_is_explicit_and_round_trips() -> None:
    parameters = smoke_parameters()
    restored = GroupedSPMParameters.from_vector(parameters.as_vector())

    assert len(PARAMETER_ORDER) == 9
    np.testing.assert_array_equal(restored.as_vector(), parameters.as_vector())

