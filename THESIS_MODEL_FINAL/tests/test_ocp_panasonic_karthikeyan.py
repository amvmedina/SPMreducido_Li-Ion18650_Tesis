from __future__ import annotations

import numpy as np
import pytest

from models.ocp_panasonic_karthikeyan import (
    DEFAULT_TEMPERATURE_K,
    FARADAY_C_PER_MOL,
    MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL,
    MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE,
    MCMB_GRAPHITE_U0_V,
    NCA_COEFFICIENTS_J_PER_MOL,
    NCA_LITERATURE_STOICHIOMETRY_RANGE,
    NCA_U0_V,
    U_n,
    U_p,
    redlich_kister_equilibrium_potential,
)


def _direct_finite_expression(
    x: np.ndarray,
    coefficients: tuple[float, ...],
) -> np.ndarray:
    t = 2.0 * x - 1.0
    result = coefficients[0] * t
    for k, coefficient in enumerate(coefficients[1:], start=1):
        result += coefficient * (
            t ** (k + 1) - 2.0 * x * k * (1.0 - x) * t ** (k - 1)
        )
    return result


@pytest.mark.parametrize(
    ("u0_v", "coefficients"),
    [
        (MCMB_GRAPHITE_U0_V, MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL),
        (NCA_U0_V, NCA_COEFFICIENTS_J_PER_MOL),
    ],
)
def test_low_level_evaluator_is_finite_through_mathematical_domain(
    u0_v: float,
    coefficients: tuple[float, ...],
) -> None:
    x = np.linspace(1.0e-8, 1.0 - 1.0e-8, 20001)
    values = redlich_kister_equilibrium_potential(
        x,
        DEFAULT_TEMPERATURE_K,
        u0_v=u0_v,
        coefficients_j_per_mol=coefficients,
    )
    assert values.shape == x.shape
    assert np.all(np.isfinite(values))


@pytest.mark.parametrize(
    ("function", "supported_range"),
    [
        (U_n, MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE),
        (U_p, NCA_LITERATURE_STOICHIOMETRY_RANGE),
    ],
)
def test_panasonic_functions_are_finite_through_supported_ranges(
    function,
    supported_range: tuple[float, float],
) -> None:
    x = np.linspace(*supported_range, 20001)
    values = function(x)
    assert values.shape == x.shape
    assert np.all(np.isfinite(values))


@pytest.mark.parametrize(
    ("function", "u0_v", "coefficients"),
    [
        (U_n, MCMB_GRAPHITE_U0_V, MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL),
        (U_p, NCA_U0_V, NCA_COEFFICIENTS_J_PER_MOL),
    ],
)
def test_continuous_and_exact_at_half_stoichiometry(
    function,
    u0_v: float,
    coefficients: tuple[float, ...],
) -> None:
    epsilon = 1.0e-9
    left = function(0.5 - epsilon)
    center = function(0.5)
    right = function(0.5 + epsilon)
    expected_center = u0_v - coefficients[1] / (2.0 * FARADAY_C_PER_MOL)
    assert center == pytest.approx(expected_center, rel=0.0, abs=1.0e-13)
    assert left == pytest.approx(center, rel=0.0, abs=1.0e-7)
    assert right == pytest.approx(center, rel=0.0, abs=1.0e-7)


@pytest.mark.parametrize(
    ("u0_v", "coefficients"),
    [
        (MCMB_GRAPHITE_U0_V, MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL),
        (NCA_U0_V, NCA_COEFFICIENTS_J_PER_MOL),
    ],
)
def test_stable_basis_matches_direct_k_minus_one_form_away_from_half(
    u0_v: float,
    coefficients: tuple[float, ...],
) -> None:
    x = np.asarray([0.03, 0.15, 0.35, 0.65, 0.85, 0.97])
    expected = (
        u0_v
        + 8.314
        * DEFAULT_TEMPERATURE_K
        / FARADAY_C_PER_MOL
        * np.log((1.0 - x) / x)
        + _direct_finite_expression(x, coefficients) / FARADAY_C_PER_MOL
    )
    actual = redlich_kister_equilibrium_potential(
        x,
        DEFAULT_TEMPERATURE_K,
        u0_v=u0_v,
        coefficients_j_per_mol=coefficients,
    )
    np.testing.assert_allclose(actual, expected, rtol=1.0e-13, atol=1.0e-13)


def test_mcmb_graphite_curve_has_plausible_scale_and_direction() -> None:
    values = U_n(np.asarray([0.01, 0.05, 0.10, 0.20, 0.50, 0.80]))
    assert np.all((values > 0.0) & (values < 1.0))
    assert np.all(np.diff(values) < 0.0)


def test_nca_curve_has_plausible_scale_and_direction() -> None:
    low, middle, high = U_p(np.asarray([0.18, 0.5, 0.96]))
    assert 3.0 < middle < 4.5
    assert low > middle > high
    assert 2.0 < high < 4.0
    assert 3.5 < low < 5.0


@pytest.mark.parametrize(
    ("u0_v", "coefficients"),
    [
        (MCMB_GRAPHITE_U0_V, MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL),
        (NCA_U0_V, NCA_COEFFICIENTS_J_PER_MOL),
    ],
)
@pytest.mark.parametrize("invalid_x", [-0.1, 0.0, 1.0, 1.1, np.nan, np.inf])
def test_low_level_evaluator_rejects_outside_mathematical_domain_without_clipping(
    u0_v: float,
    coefficients: tuple[float, ...],
    invalid_x: float,
) -> None:
    with pytest.raises(ValueError, match="stoichiometry|Stoichiometry"):
        redlich_kister_equilibrium_potential(
            invalid_x,
            DEFAULT_TEMPERATURE_K,
            u0_v=u0_v,
            coefficients_j_per_mol=coefficients,
        )


@pytest.mark.parametrize(
    ("function", "invalid_x"),
    [
        (U_n, 0.009999),
        (U_n, 0.800001),
        (U_p, 0.179999),
        (U_p, 0.960001),
    ],
)
def test_panasonic_functions_reject_literature_domain_extrapolation(
    function,
    invalid_x: float,
) -> None:
    with pytest.raises(ValueError, match="literature-supported"):
        function(invalid_x)


def test_low_level_and_panasonic_facing_domains_are_separate() -> None:
    extrapolated = redlich_kister_equilibrium_potential(
        0.90,
        DEFAULT_TEMPERATURE_K,
        u0_v=MCMB_GRAPHITE_U0_V,
        coefficients_j_per_mol=MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL,
    )
    assert np.isfinite(extrapolated)
    with pytest.raises(ValueError, match="No extrapolation"):
        U_n(0.90)


@pytest.mark.parametrize("function", [U_n, U_p])
def test_invalid_array_is_not_partially_clipped(function) -> None:
    with pytest.raises(ValueError, match="No extrapolation"):
        function(np.asarray([0.2, 0.0, 0.8]))


@pytest.mark.parametrize(
    ("function", "stoichiometry"),
    [
        (U_n, np.asarray([0.01, 0.50, 0.80])),
        (U_p, np.asarray([0.18, 0.50, 0.96])),
    ],
)
def test_temperature_broadcast_and_validation(function, stoichiometry: np.ndarray) -> None:
    result = function(stoichiometry, np.asarray([273.15, 298.15, 323.15]))
    assert result.shape == (3,)
    assert np.all(np.isfinite(result))
    with pytest.raises(ValueError, match="temperature"):
        function(0.5, 0.0)
