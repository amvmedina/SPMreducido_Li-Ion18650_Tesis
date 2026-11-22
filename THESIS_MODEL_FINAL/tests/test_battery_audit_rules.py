from __future__ import annotations

import numpy as np

from analysis.battery_audit_rules import (
    TIME_CONFIRMED_INTERVAL,
    TIME_ESTIMATED,
    CommonAuditRecord,
    evaluate_common_audit,
)
from analysis.battery_audit_app_service import ROLE_DYNAMIC, ROLE_LOW_RATE, ROLE_VALIDATION, SETUP_UNKNOWN


def _record(
    *,
    record_id: str = "record",
    name: str = "Ensayo",
    current_a: np.ndarray | None = None,
    voltage_v: np.ndarray | None = None,
    capacity_ah: np.ndarray | None = None,
    roles: tuple[str, ...] = (),
    time_s: np.ndarray | None = None,
    measured_time: bool = False,
    discharge_cutoff_v: float | None = 3.0,
    pulse_candidate: bool = False,
    time_basis: str = "unavailable",
    cell_temperature_c: np.ndarray | None = None,
    ambient_temperature_c: np.ndarray | None = None,
) -> CommonAuditRecord:
    current = np.asarray([0.0, -0.15, -0.15, -0.15, 0.0]) if current_a is None else current_a
    voltage = np.asarray([4.2, 4.1, 3.6, 3.0, 3.1]) if voltage_v is None else voltage_v
    capacity = np.asarray([0.0, 0.0, 1.0, 3.0, 3.0]) if capacity_ah is None else capacity_ah
    return CommonAuditRecord(
        record_id=record_id,
        name=name,
        voltage_v=voltage,
        current_a=current,
        capacity_ah=capacity,
        roles=roles,
        setup_role=SETUP_UNKNOWN,
        time_s=time_s,
        measured_time=measured_time,
        time_basis=time_basis,
        cell_temperature_c=cell_temperature_c,
        ambient_temperature_c=ambient_temperature_c,
        discharge_cutoff_v=discharge_cutoff_v,
        pulse_candidate=pulse_candidate,
    )


def _states(record: CommonAuditRecord, nominal_capacity_ah: float | None = 3.0) -> dict[str, str]:
    return {item.code: item.status for item in evaluate_common_audit((record,), expected_record_count=1, nominal_capacity_ah=nominal_capacity_ah)}


def test_same_signals_give_same_states_regardless_of_names() -> None:
    first = _record(record_id="first", name="Marca A")
    second = _record(record_id="second", name="Modalidad B")
    assert _states(first) == _states(second)


def test_incomplete_capacity_is_not_invalid_data() -> None:
    record = _record(voltage_v=np.asarray([4.2, 4.1, 3.9, 3.7, 3.8]), discharge_cutoff_v=3.0)
    states = _states(record)
    assert states["E02"] == "FALTANTE"
    assert states["E03"] == "CONFORME"


def test_e01_requires_measurements_not_simulations() -> None:
    assert _states(_record())["E01"] == "CONFORME"


def test_e04_is_not_evaluated_without_comparable_mounts() -> None:
    assert _states(_record())["E04"] == "NO EVALUADO"


def test_declared_low_rate_and_hppc_roles_do_not_guarantee_conformity() -> None:
    high_current = _record(current_a=np.asarray([0.0, -3.0, -3.0, -3.0, 0.0]), roles=(ROLE_LOW_RATE,))
    assert _states(high_current)["E05"] == "NO APTO"

    pulse_without_time = _record(
        current_a=np.asarray([0.0, -3.0, 0.0, -3.0, 0.0]),
        roles=(ROLE_DYNAMIC,),
        pulse_candidate=True,
    )
    assert _states(pulse_without_time)["E06"] == "NO APTO"


def test_removing_required_capacity_changes_e01_only_as_expected() -> None:
    complete = _record(roles=(ROLE_VALIDATION,))
    without_capacity = _record(capacity_ah=np.asarray([]), roles=(ROLE_VALIDATION,))
    assert _states(complete)["E01"] == "CONFORME"
    assert _states(without_capacity)["E01"] == "FALTANTE"


def test_confirmed_sampling_interval_satisfies_thermal_time_requirement() -> None:
    temperature = np.asarray([25.0, 25.1, 25.2, 25.3, 25.4])
    confirmed = _record(
        time_s=np.arange(5, dtype=float) * 3.0,
        time_basis=TIME_CONFIRMED_INTERVAL,
        cell_temperature_c=temperature,
        ambient_temperature_c=np.full(5, 25.0),
    )
    assert _states(confirmed)["E07"] == "CONFORME"


def test_unverified_estimated_interval_requires_review_and_confirmation_changes_result() -> None:
    temperature = np.asarray([25.0, 25.1, 25.2, 25.3, 25.4])
    estimated = _record(
        time_s=np.arange(5, dtype=float) * 3.0,
        time_basis=TIME_ESTIMATED,
        cell_temperature_c=temperature,
        ambient_temperature_c=np.full(5, 25.0),
    )
    confirmed = _record(
        time_s=np.arange(5, dtype=float) * 3.0,
        time_basis=TIME_CONFIRMED_INTERVAL,
        cell_temperature_c=temperature,
        ambient_temperature_c=np.full(5, 25.0),
    )
    estimated_outcomes = {item.code: item for item in evaluate_common_audit((estimated,), expected_record_count=1, nominal_capacity_ah=3.0)}
    assert estimated_outcomes["E07"].status == "REVISAR"
    assert "Confirmar el intervalo de muestreo" in estimated_outcomes["E07"].action
    assert _states(confirmed)["E07"] == "CONFORME"
