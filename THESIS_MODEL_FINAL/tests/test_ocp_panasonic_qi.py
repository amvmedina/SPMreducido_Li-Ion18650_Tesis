from __future__ import annotations

import numpy as np
import pytest

from models.ocp_panasonic_qi import (
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
    N_UP_KNOTS,
    UP_KNOT_STOICHIOMETRIES,
    U_n,
    U_p,
    validate_positive_knot_voltages,
)


def test_guo_mcmb_reference_values() -> None:
    x = np.asarray([0.80, 0.40, 0.20, 0.10, 0.05, 0.01])
    expected = np.asarray([0.11526, 0.15258, 0.18262, 0.23117, 0.27967, 0.63407])
    np.testing.assert_allclose(U_n(x), expected, rtol=0.0, atol=6.0e-5)


def test_guo_mcmb_is_finite_through_supported_domain_and_continuous_at_knots() -> None:
    values = U_n(np.linspace(*MCMB_SUPPORTED_STOICHIOMETRY_RANGE, 20001))
    assert np.all(np.isfinite(values))
    for knot in (0.47423, 0.16833, 0.078218, 0.015397):
        epsilon = 1.0e-10
        assert U_n(knot - epsilon) == pytest.approx(U_n(knot), abs=1.0e-8)
        assert U_n(knot + epsilon) == pytest.approx(U_n(knot), abs=1.0e-8)


@pytest.mark.parametrize("invalid", [0.0, 0.009999, 0.800001, 1.0, np.nan, np.inf])
def test_guo_mcmb_rejects_unsupported_values_without_clipping(invalid: float) -> None:
    with pytest.raises(ValueError, match="No clipping|finite"):
        U_n(invalid)


def test_positive_ocp_interpolates_exactly_at_fixed_knots() -> None:
    voltages = np.linspace(4.35, 2.75, N_UP_KNOTS)
    np.testing.assert_allclose(U_p(UP_KNOT_STOICHIOMETRIES, voltages), voltages)
    midpoint = 0.5 * (UP_KNOT_STOICHIOMETRIES[4] + UP_KNOT_STOICHIOMETRIES[5])
    assert U_p(midpoint, voltages) == pytest.approx(0.5 * (voltages[4] + voltages[5]))


def test_positive_ocp_is_finite_and_nonincreasing_through_working_domain() -> None:
    voltages = np.linspace(4.3, 2.8, N_UP_KNOTS)
    values = U_p(np.linspace(*NCA_WORKING_STOICHIOMETRY_RANGE, 20001), voltages)
    assert np.all(np.isfinite(values))
    assert np.all(np.diff(values) <= 0.0)


@pytest.mark.parametrize("invalid", [0.0, 0.179999, 0.960001, 1.0, np.nan, np.inf])
def test_positive_ocp_rejects_unsupported_values_without_clipping(invalid: float) -> None:
    with pytest.raises(ValueError, match="No clipping|finite"):
        U_p(invalid, np.linspace(4.3, 2.8, N_UP_KNOTS))


def test_positive_knot_constraints_are_explicit() -> None:
    valid = np.linspace(4.3, 2.8, N_UP_KNOTS)
    np.testing.assert_array_equal(validate_positive_knot_voltages(valid), valid)

    increasing = valid.copy()
    increasing[6] = increasing[5] + 0.01
    with pytest.raises(ValueError, match="non-increasing"):
        validate_positive_knot_voltages(increasing)

    below_bound = valid.copy()
    below_bound[-1] = 2.49
    with pytest.raises(ValueError, match="no clipping"):
        validate_positive_knot_voltages(below_bound)

    with pytest.raises(ValueError, match="Exactly 13"):
        validate_positive_knot_voltages(valid[:-1])
