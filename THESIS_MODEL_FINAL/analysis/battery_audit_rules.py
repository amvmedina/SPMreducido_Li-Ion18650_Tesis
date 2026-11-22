"""Brand-independent E01–E08 rules for experimental battery records."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from analysis.battery_dataset_assessment import NonidealityFlag, evaluate_relative_nonideality_gate


STATUS_CONFORME = "CONFORME"
STATUS_REVISAR = "REVISAR"
STATUS_FALTANTE = "FALTANTE"
STATUS_NO_APTO = "NO APTO"
STATUS_NO_EVALUADO = "NO EVALUADO"

ROLE_LOW_RATE = "referencia de baja tasa"
ROLE_DYNAMIC = "pulsos para identificación dinámica"
ROLE_THERMAL = "identificación térmica"
ROLE_VALIDATION = "validación independiente"
SETUP_REFERENCE = "montaje de referencia"

REST_CURRENT_A = 0.05
LOW_RATE_MAX_C = 0.10
MINIMUM_PULSE_TRANSITIONS = 2
TIME_RECORDED = "recorded"
TIME_CONFIRMED_INTERVAL = "confirmed_interval"
TIME_ESTIMATED = "estimated_current_capacity"
TIME_UNAVAILABLE = "unavailable"


@dataclass(frozen=True)
class CommonAuditRecord:
    record_id: str
    name: str
    voltage_v: np.ndarray
    current_a: np.ndarray
    capacity_ah: np.ndarray
    roles: tuple[str, ...]
    setup_role: str
    include_in_mount_comparison: bool = True
    time_s: np.ndarray | None = None
    measured_time: bool = False
    cell_temperature_c: np.ndarray | None = None
    ambient_temperature_c: np.ndarray | None = None
    discharge_cutoff_v: float | None = None
    charge_cutoff_v: float | None = None
    capacity_test_eligible: bool = True
    pulse_candidate: bool = False
    time_basis: str = TIME_UNAVAILABLE


@dataclass(frozen=True)
class RuleOutcome:
    code: str
    status: str
    evidence: str
    action: str = ""


def _arrays(record: CommonAuditRecord) -> tuple[tuple[str, np.ndarray], ...]:
    output = [
        ("voltaje", np.asarray(record.voltage_v, dtype=float)),
        ("corriente", np.asarray(record.current_a, dtype=float)),
        ("capacidad", np.asarray(record.capacity_ah, dtype=float)),
    ]
    if record.time_s is not None:
        output.append(("tiempo", np.asarray(record.time_s, dtype=float)))
    if record.cell_temperature_c is not None:
        output.append(("temperatura de la celda", np.asarray(record.cell_temperature_c, dtype=float)))
    if record.ambient_temperature_c is not None:
        output.append(("temperatura ambiente", np.asarray(record.ambient_temperature_c, dtype=float)))
    return tuple(output)


def _required_available(record: CommonAuditRecord) -> bool:
    sizes = [np.asarray(values).size for _, values in _arrays(record)[:3]]
    return min(sizes, default=0) >= 2 and len(set(sizes)) == 1


def _invalid_details(record: CommonAuditRecord) -> list[str]:
    details: list[str] = []
    required_size = np.asarray(record.voltage_v).size
    for label, values in _arrays(record):
        if values.size != required_size:
            details.append(f"{label}: longitud diferente")
            continue
        count = int(np.count_nonzero(~np.isfinite(values)))
        if count:
            details.append(f"{label}: {count} muestras")
    return details


def _direction_runs(current_a: np.ndarray, direction: str) -> list[tuple[int, int]]:
    current = np.asarray(current_a, dtype=float)
    selected = current < -REST_CURRENT_A if direction == "discharge" else current > REST_CURRENT_A
    indices = np.flatnonzero(selected)
    if not indices.size:
        return []
    breaks = np.flatnonzero(np.diff(indices) > 1)
    starts = np.concatenate(([0], breaks + 1))
    stops = np.concatenate((breaks, [indices.size - 1]))
    return [(int(indices[start]), int(indices[stop])) for start, stop in zip(starts, stops, strict=True)]


def _capacity_complete(record: CommonAuditRecord) -> tuple[bool, bool]:
    """Return (complete, assessable) using record-specific declared cutoffs."""

    if not record.capacity_test_eligible:
        return False, False
    voltage = np.asarray(record.voltage_v, dtype=float)
    if not _required_available(record) or not np.all(np.isfinite(voltage)):
        return False, False
    assessable = False
    if record.discharge_cutoff_v is not None:
        assessable = True
        for _, stop in _direction_runs(record.current_a, "discharge"):
            if voltage[stop] <= record.discharge_cutoff_v:
                return True, True
    if record.charge_cutoff_v is not None:
        assessable = True
        for _, stop in _direction_runs(record.current_a, "charge"):
            if voltage[stop] >= record.charge_cutoff_v:
                return True, True
    return False, assessable


def _rest_to_discharge(record: CommonAuditRecord) -> np.ndarray:
    voltage = np.asarray(record.voltage_v, dtype=float)
    current = np.asarray(record.current_a, dtype=float)
    if voltage.size != current.size or voltage.size < 2:
        return np.asarray([], dtype=float)
    transitions = np.flatnonzero((np.abs(current[:-1]) <= REST_CURRENT_A) & (current[1:] < -REST_CURRENT_A)) + 1
    values: list[float] = []
    for index in transitions:
        delta_i = current[index] - current[index - 1]
        if abs(delta_i) <= REST_CURRENT_A:
            continue
        value = abs((voltage[index - 1] - voltage[index]) / delta_i)
        if np.isfinite(value):
            values.append(float(value))
    return np.asarray(values, dtype=float)


def _iqr(values: np.ndarray) -> float:
    return float(np.percentile(values, 75.0) - np.percentile(values, 25.0)) if values.size else float("nan")


def _time_valid(record: CommonAuditRecord) -> bool:
    verified_basis = record.measured_time or record.time_basis in {TIME_RECORDED, TIME_CONFIRMED_INTERVAL}
    if not verified_basis or record.time_s is None:
        return False
    time = np.asarray(record.time_s, dtype=float)
    differences = np.diff(time)
    return time.size == record.voltage_v.size and np.all(np.isfinite(time)) and np.all(differences >= 0.0) and np.any(differences > 0.0)


def _time_pending_confirmation(record: CommonAuditRecord) -> bool:
    if record.time_basis != TIME_ESTIMATED or record.time_s is None:
        return False
    time = np.asarray(record.time_s, dtype=float)
    return time.size == record.voltage_v.size and np.all(np.isfinite(time)) and np.all(np.diff(time) >= 0.0)


def evaluate_common_audit(
    records: tuple[CommonAuditRecord, ...],
    *,
    expected_record_count: int,
    load_error_names: tuple[str, ...] = (),
    nominal_capacity_ah: float | None = None,
) -> tuple[RuleOutcome, ...]:
    """Apply one set of rules irrespective of source format or cell name."""

    outcomes: list[RuleOutcome] = []
    electrically_available = [record for record in records if _required_available(record)]
    if len(electrically_available) == expected_record_count:
        e01_status = STATUS_CONFORME
        e01_evidence = f"Los {expected_record_count} ensayos contienen voltaje, corriente y capacidad."
        e01_action = ""
    elif electrically_available:
        e01_status = STATUS_REVISAR
        e01_evidence = f"{len(electrically_available)} de {expected_record_count} ensayos contienen las tres variables eléctricas requeridas."
        e01_action = "Revisar los ensayos que no contienen voltaje, corriente y capacidad."
    else:
        e01_status = STATUS_FALTANTE
        e01_evidence = "Ningún ensayo contiene simultáneamente voltaje, corriente y capacidad."
        e01_action = "Incluir voltaje, corriente y capacidad en al menos un ensayo."
    if load_error_names:
        e01_evidence += f" No pudieron cargarse: {', '.join(load_error_names)}."
    outcomes.append(RuleOutcome("E01", e01_status, e01_evidence, e01_action))

    completeness = {record.record_id: _capacity_complete(record) for record in electrically_available}
    complete = [record for record in electrically_available if completeness[record.record_id][0]]
    assessable = [record for record in electrically_available if completeness[record.record_id][1]]
    if complete:
        names = ", ".join(record.name for record in complete)
        outcomes.append(RuleOutcome("E02", STATUS_CONFORME, f"Existe al menos un ensayo completo según sus límites declarados: {names}."))
    elif not assessable:
        outcomes.append(RuleOutcome(
            "E02", STATUS_FALTANTE,
            "No se declararon límites de terminación suficientes para confirmar una carga o descarga completa.",
            "Indicar el límite de voltaje de terminación correspondiente al ensayo.",
        ))
    else:
        outcomes.append(RuleOutcome(
            "E02", STATUS_FALTANTE,
            "Ningún ensayo alcanza el límite de terminación declarado.",
            "Incluir una carga o descarga que alcance el límite declarado de la celda.",
        ))

    invalid_rows: list[str] = []
    for record in records:
        details = _invalid_details(record)
        if details:
            invalid_rows.append(f"{record.name} ({'; '.join(details)})")
    if not records:
        outcomes.append(RuleOutcome("E03", STATUS_FALTANTE, "No hay ensayos disponibles para comprobar su integridad.", "Cargar al menos un ensayo."))
    elif invalid_rows:
        outcomes.append(RuleOutcome(
            "E03", STATUS_REVISAR,
            f"Se encontraron datos faltantes o inválidos en: {'; '.join(invalid_rows)}.",
            "Corregir únicamente las variables y muestras identificadas.",
        ))
    else:
        outcomes.append(RuleOutcome("E03", STATUS_CONFORME, f"Los {len(records)} ensayos no contienen datos faltantes o inválidos."))

    comparison_records = [record for record in electrically_available if record.include_in_mount_comparison]
    references = [record for record in comparison_records if record.setup_role == SETUP_REFERENCE]
    transition_values = {record.record_id: _rest_to_discharge(record) for record in comparison_records}
    comparable = [record for record in comparison_records if transition_values[record.record_id].size]
    elevated_ids: set[str] = set()
    if len(references) != 1 or len(comparable) < 2:
        outcomes.append(RuleOutcome("E04", STATUS_NO_EVALUADO, "No hay dos montajes comparables con una referencia y transiciones reposo–descarga."))
    else:
        rows = []
        for record in comparable:
            values = transition_values[record.record_id]
            rows.append({
                "dataset_id": record.record_id,
                "discharge_r_app_count": int(values.size),
                "discharge_r_app_median_ohm": float(np.median(values)),
                "discharge_r_app_iqr_ohm": _iqr(values),
            })
        try:
            comparison = evaluate_relative_nonideality_gate(
                rows,
                references[0].record_id,
                transition_count_field="discharge_r_app_count",
                median_field="discharge_r_app_median_ohm",
                iqr_field="discharge_r_app_iqr_ohm",
                minimum_transition_count=1,
                metric_label="R_app,desc",
            )
            by_id = {record.record_id: record for record in comparable}
            parts = []
            for item in comparison:
                suffix = "referencia" if item.dataset == references[0].record_id else f"{item.ratio_to_reference:.3g} veces la referencia"
                parts.append(f"{by_id[item.dataset].name}: {item.median_r_app_ohm:.6g} Ω ({suffix})")
                if item.flag == NonidealityFlag.ELEVATED_SERIES_PATH:
                    elevated_ids.add(item.dataset)
            status = STATUS_REVISAR if elevated_ids else STATUS_CONFORME
            action = "Revisar el montaje y las conexiones; la comparación no identifica una causa física única." if elevated_ids else ""
            outcomes.append(RuleOutcome("E04", status, "; ".join(parts) + ".", action))
        except ValueError:
            outcomes.append(RuleOutcome("E04", STATUS_REVISAR, "La referencia seleccionada no permite una comparación relativa consistente.", "Verificar que exista una única referencia comparable."))

    low_rate_records = [record for record in records if ROLE_LOW_RATE in record.roles]
    if not low_rate_records:
        outcomes.append(RuleOutcome("E05", STATUS_FALTANTE, "No se incluyó un ensayo de baja corriente.", "Realizar una descarga C/20."))
    else:
        suitable_low_rate: list[CommonAuditRecord] = []
        metadata_missing = False
        for record in low_rate_records:
            active = np.abs(record.current_a) > REST_CURRENT_A
            if nominal_capacity_ah is None or nominal_capacity_ah <= 0.0 or not completeness.get(record.record_id, (False, False))[1]:
                metadata_missing = True
                continue
            c_rate = float(np.median(np.abs(record.current_a[active])) / nominal_capacity_ah) if np.any(active) else float("inf")
            if c_rate <= LOW_RATE_MAX_C and completeness.get(record.record_id, (False, False))[0] and bool(_direction_runs(record.current_a, "discharge")):
                suitable_low_rate.append(record)
        if suitable_low_rate:
            outcomes.append(RuleOutcome("E05", STATUS_CONFORME, f"El ensayo de baja corriente es apto: {', '.join(record.name for record in suitable_low_rate)}."))
        elif metadata_missing:
            outcomes.append(RuleOutcome("E05", STATUS_REVISAR, "El ensayo de baja corriente está declarado, pero faltan la capacidad nominal o el límite de terminación para comprobarlo.", "Indicar la capacidad nominal y el límite de descarga."))
        else:
            outcomes.append(RuleOutcome("E05", STATUS_NO_APTO, "El ensayo declarado no combina baja corriente y descarga completa.", "Realizar una descarga C/20 completa."))

    pulse_records = [record for record in records if ROLE_DYNAMIC in record.roles or record.pulse_candidate]
    if not pulse_records:
        outcomes.append(RuleOutcome("E06", STATUS_FALTANTE, "No se incluyó un ensayo de pulsos con periodos de reposo.", "Realizar un ensayo HPPC con reposos, pulsos de corriente y tiempo medido."))
    else:
        structurally_valid = [record for record in pulse_records if _rest_to_discharge(record).size >= MINIMUM_PULSE_TRANSITIONS]
        suitable_pulses = [record for record in structurally_valid if _time_valid(record)]
        if len(suitable_pulses) == len(pulse_records):
            affected = [record.name for record in suitable_pulses if record.record_id in elevated_ids]
            evidence = f"El ensayo de pulsos contiene tiempo medido, reposos y al menos {MINIMUM_PULSE_TRANSITIONS} pulsos de descarga."
            action = ""
            status = STATUS_CONFORME
            if affected:
                status = STATUS_REVISAR
                evidence += f" E04 detectó una respuesta elevada en: {', '.join(affected)}."
                action = "Revisar el montaje antes de utilizar los pulsos."
            outcomes.append(RuleOutcome("E06", status, evidence, action))
        elif len(structurally_valid) == len(pulse_records) and all(
            _time_pending_confirmation(record) or record.time_basis == TIME_UNAVAILABLE
            for record in pulse_records
        ):
            outcomes.append(RuleOutcome(
                "E06", STATUS_REVISAR,
                f"El ensayo contiene reposos y al menos {MINIMUM_PULSE_TRANSITIONS} pulsos de descarga, pero la escala temporal no está confirmada.",
                "Confirmar el intervalo de muestreo y la continuidad del registro.",
            ))
        else:
            outcomes.append(RuleOutcome("E06", STATUS_NO_APTO, "Los pulsos existen, pero la secuencia temporal, los reposos o el número de pulsos no permiten evaluarlos.", "Corregir el protocolo HPPC o la información temporal indicada."))

    thermal_records = [record for record in records if record.cell_temperature_c is not None or record.ambient_temperature_c is not None]
    complete_thermal = [
        record for record in thermal_records
        if record.cell_temperature_c is not None
        and record.ambient_temperature_c is not None
    ]
    valid_thermal = [record for record in complete_thermal if not _invalid_details(record)]
    suitable_thermal = [record for record in valid_thermal if _time_valid(record)]
    pending_thermal = [
        record for record in valid_thermal
        if not _time_valid(record)
        and (record.time_basis in {TIME_ESTIMATED, TIME_UNAVAILABLE})
    ]
    if not thermal_records:
        outcomes.append(RuleOutcome("E07", STATUS_FALTANTE, "No se encontraron mediciones de temperatura.", "Medir la temperatura de la celda y la temperatura ambiente."))
    elif suitable_thermal:
        outcomes.append(RuleOutcome("E07", STATUS_CONFORME, f"Hay mediciones térmicas aptas en: {', '.join(record.name for record in suitable_thermal)}."))
    elif pending_thermal:
        outcomes.append(RuleOutcome(
            "E07", STATUS_REVISAR,
            "Se dispone de mediciones de temperatura de la celda y del ambiente.",
            "Confirmar el intervalo de muestreo del equipo y la continuidad del registro.",
        ))
    elif not complete_thermal:
        missing_cell = not any(record.cell_temperature_c is not None for record in thermal_records)
        missing_ambient = not any(record.ambient_temperature_c is not None for record in thermal_records)
        missing = []
        if missing_cell:
            missing.append("temperatura de la celda")
        if missing_ambient:
            missing.append("temperatura ambiente")
        if not missing:
            missing.append("temperatura de la celda y ambiente asociadas al mismo ensayo")
        outcomes.append(RuleOutcome("E07", STATUS_FALTANTE, f"Falta: {', '.join(missing)}.", "Incluir las mediciones de temperatura indicadas."))
    else:
        outcomes.append(RuleOutcome("E07", STATUS_NO_APTO, "Las mediciones térmicas contienen datos inválidos o una escala temporal inconsistente.", "Corregir las mediciones identificadas antes de utilizarlas."))

    validation_records = [record for record in records if ROLE_VALIDATION in record.roles]
    calibration_roles = {ROLE_LOW_RATE, ROLE_DYNAMIC, ROLE_THERMAL}
    independent = [
        record for record in validation_records
        if not calibration_roles.intersection(record.roles)
        and _required_available(record)
        and not _invalid_details(record)
    ]
    if independent:
        outcomes.append(RuleOutcome("E08", STATUS_CONFORME, f"Existe validación independiente: {', '.join(record.name for record in independent)}."))
    elif validation_records:
        outcomes.append(RuleOutcome("E08", STATUS_NO_APTO, "Los ensayos declarados para validación también se usan en calibración o no contienen las variables requeridas.", "Reservar un ensayo válido que no se utilice para calibrar."))
    else:
        outcomes.append(RuleOutcome("E08", STATUS_FALTANTE, "No se reservó un ensayo independiente para validación.", "Reservar un ensayo que no se utilice para calibrar."))

    return tuple(outcomes)


__all__ = ["CommonAuditRecord", "RuleOutcome", "evaluate_common_audit"]
