from __future__ import annotations

import numpy as np
import pytest

from models.panasonic_entropic import (
    DU_EQ_DT_MV_PER_K,
    DU_EQ_DT_V_PER_K,
    SOC_AVERAGE,
    SUPPORTED_SOC_RANGE,
    dUeq_dT_panasonic,
)


def test_exact_tabulated_values_are_returned_in_volts_per_kelvin() -> None:
    expected_mv_per_k = np.asarray(
        [0.069, -0.082, 0.058, 0.141, 0.169, 0.143, 0.065, 0.083, 0.083]
    )
    np.testing.assert_array_equal(DU_EQ_DT_MV_PER_K, expected_mv_per_k)
    np.testing.assert_allclose(
        DU_EQ_DT_V_PER_K,
        expected_mv_per_k * 1.0e-3,
        rtol=0.0,
        atol=1.0e-18,
    )
    np.testing.assert_allclose(
        dUeq_dT_panasonic(SOC_AVERAGE),
        expected_mv_per_k * 1.0e-3,
        rtol=0.0,
        atol=1.0e-18,
    )


def test_linear_interpolation() -> None:
    assert dUeq_dT_panasonic(0.55) == pytest.approx(0.156e-3, abs=1.0e-18)
    # Midway between +0.069 and -0.082 mV/K.
    assert dUeq_dT_panasonic(0.15) == pytest.approx(-0.0065e-3, abs=1.0e-18)


def test_units_are_volts_per_kelvin_not_millivolts_per_kelvin() -> None:
    assert dUeq_dT_panasonic(0.90) == pytest.approx(0.000083, abs=1.0e-18)
    assert dUeq_dT_panasonic(0.50) == pytest.approx(0.000169, abs=1.0e-18)


def test_negative_coefficient_at_soc_point_two() -> None:
    assert dUeq_dT_panasonic(0.20) == pytest.approx(-0.000082, abs=1.0e-18)
    assert dUeq_dT_panasonic(0.20) < 0.0


@pytest.mark.parametrize(
    "invalid_soc",
    [0.0, 0.099999, 0.900001, 1.0, np.nan, np.inf, -np.inf],
)
def test_domain_errors_are_explicit_and_never_clip(invalid_soc: float) -> None:
    with pytest.raises(ValueError, match="finite|no clipping or extrapolation"):
        dUeq_dT_panasonic(invalid_soc)


def test_array_domain_error_rejects_the_entire_call() -> None:
    with pytest.raises(ValueError, match="no clipping or extrapolation"):
        dUeq_dT_panasonic(np.asarray([SUPPORTED_SOC_RANGE[0], 0.5, 0.91]))


def test_supported_domain_is_finite_and_inclusive() -> None:
    grid = np.linspace(*SUPPORTED_SOC_RANGE, 1001)
    assert np.all(np.isfinite(dUeq_dT_panasonic(grid)))
    assert np.isfinite(dUeq_dT_panasonic(SUPPORTED_SOC_RANGE[0]))
    assert np.isfinite(dUeq_dT_panasonic(SUPPORTED_SOC_RANGE[1]))
