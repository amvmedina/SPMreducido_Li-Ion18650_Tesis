"""Regression and physical smoke tests for the coupled model."""

import numpy as np

from models import (
    CellSOCParameters,
    GroupedSPMParameters,
    LumpedThermalParameters,
    simulate_electrothermal_spm,
    simulate_grouped_spm,
    simulate_lumped_thermal,
)


def electrical_parameters() -> GroupedSPMParameters:
    # Numerical fixture only; it is not a calibrated battery parameter set.
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


def thermal_parameters(
    *, t_initial_k: float = 298.15, t_ambient_k: float = 298.15
) -> LumpedThermalParameters:
    return LumpedThermalParameters(
        c_th_j_per_k=500.0,
        g_th_w_per_k=0.5,
        t_initial_k=t_initial_k,
        t_ambient_k=t_ambient_k,
    )


def test_electrochemical_prescribed_temperature_regression() -> None:
    time_s = np.linspace(0.0, 20.0, 21)
    current_a = np.full(time_s.size, 0.2)
    result = simulate_grouped_spm(
        time_s, current_a, electrical_parameters(), temperature=301.15
    )

    np.testing.assert_allclose(
        result.voltage_v[[0, -1]], [3.33360684, 3.33337151], rtol=2e-9
    )
    np.testing.assert_allclose(
        result.x_n_avg[[0, -1]], [0.6, 0.5996363636363636], rtol=1e-14
    )
    np.testing.assert_allclose(
        result.x_p_avg[[0, -1]], [0.4, 0.4003636363636364], rtol=1e-14
    )


def test_constant_discharge_current_decreases_cell_soc_by_coulomb_counting() -> None:
    time_s = np.array([0.0, 180.0, 360.0])
    current_a = np.full(time_s.size, 1.0)
    result = simulate_electrothermal_spm(
        time_s,
        current_a,
        electrical_parameters(),
        CellSOCParameters(q_cell_ah=2.0, soc_cell_0=0.8),
        thermal_parameters(),
    )

    np.testing.assert_allclose(result.cell_soc, [0.8, 0.775, 0.75])
    assert result.cell_soc_domain_valid
    assert result.first_soc_violation is None


def test_soc_is_not_clipped_and_reports_violation() -> None:
    result = simulate_electrothermal_spm(
        np.array([0.0, 3600.0]),
        np.array([1.0, 1.0]),
        electrical_parameters(),
        CellSOCParameters(q_cell_ah=1.0, soc_cell_0=0.5),
        thermal_parameters(),
    )

    assert result.cell_soc[-1] == -0.5
    assert not result.cell_soc_domain_valid
    assert result.first_soc_violation is not None
    assert result.first_soc_violation.sample_index == 1


def test_thermal_equilibrium_at_zero_heat() -> None:
    time_s = np.array([0.0, 10.0, 100.0])
    result = simulate_lumped_thermal(
        time_s, np.zeros(time_s.size), thermal_parameters()
    )
    np.testing.assert_array_equal(result.temperature_k, np.full(3, 298.15))


def test_coupled_zero_current_at_ambient_remains_at_equilibrium() -> None:
    time_s = np.array([0.0, 10.0, 100.0])
    result = simulate_electrothermal_spm(
        time_s,
        np.zeros(time_s.size),
        electrical_parameters(),
        CellSOCParameters(q_cell_ah=3.0, soc_cell_0=0.7),
        thermal_parameters(),
    )

    np.testing.assert_array_equal(result.temperature_k, np.full(3, 298.15))
    np.testing.assert_array_equal(
        result.heat_generation_total_w, np.zeros(time_s.size)
    )


def test_passive_cooling_decays_toward_ambient() -> None:
    time_s = np.array([0.0, 100.0, 200.0])
    result = simulate_lumped_thermal(
        time_s,
        np.zeros(time_s.size),
        thermal_parameters(t_initial_k=318.15, t_ambient_k=298.15),
    )

    assert 298.15 < result.temperature_k[-1] < result.temperature_k[1] < 318.15


def test_positive_heat_generation_increases_temperature() -> None:
    time_s = np.array([0.0, 10.0, 20.0])
    result = simulate_lumped_thermal(
        time_s, np.full(time_s.size, 10.0), thermal_parameters()
    )

    assert result.temperature_k[-1] > result.temperature_k[1] > 298.15


def test_short_coupled_profile_has_finite_consistent_outputs() -> None:
    time_s = np.arange(0.0, 31.0)
    current_a = np.r_[np.full(10, 0.4), np.full(10, -0.2), np.zeros(11)]
    result = simulate_electrothermal_spm(
        time_s,
        current_a,
        electrical_parameters(),
        CellSOCParameters(q_cell_ah=3.0, soc_cell_0=0.7),
        thermal_parameters(),
    )

    requested = (
        result.time_s,
        result.current_a,
        result.terminal_voltage_v,
        result.cell_soc,
        result.temperature_k,
        result.heat_generation_total_w,
        result.heat_irreversible_w,
        result.heat_ohmic_w,
        result.heat_kinetic_w,
        result.heat_reversible_w,
        result.u_eq_v,
        result.x_n_avg,
        result.x_n_surface,
        result.x_p_avg,
        result.x_p_surface,
    )
    assert all(values.shape == time_s.shape for values in requested)
    assert all(np.all(np.isfinite(values)) for values in requested)
    assert result.electrochemical_domain_valid
    assert result.cell_soc_domain_valid
    assert np.all(result.heat_irreversible_nonnegative)
    np.testing.assert_array_equal(result.heat_reversible_w, np.zeros(time_s.size))
    assert not result.reversible_heat_enabled


def test_predicted_temperature_is_fed_back_to_voltage() -> None:
    time_s = np.array([0.0, 100.0, 200.0])
    current_a = np.full(time_s.size, 1.0)
    coupled = simulate_electrothermal_spm(
        time_s,
        current_a,
        electrical_parameters(),
        CellSOCParameters(q_cell_ah=3.0, soc_cell_0=0.8),
        LumpedThermalParameters(1.0, 0.0, 298.15, 298.15),
    )
    prescribed = simulate_grouped_spm(
        time_s, current_a, electrical_parameters(), temperature=298.15
    )

    assert coupled.temperature_k[1] > coupled.temperature_k[0]
    assert coupled.terminal_voltage_v[1] != prescribed.voltage_v[1]
