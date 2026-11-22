"""Minimal evidence-gated automatic audit for the TRC18650 case study.

The module orchestrates existing data-quality and transfer gates.  It does not
run an optimizer, alter raw measurements, correct voltage, or claim that the
apparent rest-to-load response is intrinsic cell resistance.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from analysis.battery_dataset_assessment import (
    DataQualityGateResult,
    NonidealityFlag,
    RelativeNonidealityResult,
    SuitabilityStatus,
    TransferDecision,
    TransferEvidence,
    evaluate_data_quality_gate,
    evaluate_relative_nonideality_gate,
    evaluate_transfer_gate,
)
from analysis.experimental_data_audit import StandardExperiment


REFERENCE_DATASET_ID = "cycle_cell1_ambient"
TRC_REFERENCE_CAPACITY_AH = 2.8999
TRC_REFERENCE_R_APP_DESC_OHM = 0.036383074111115


@dataclass(frozen=True)
class PanasonicModelSnapshot:
    q_cell_ah: float
    q_n_c: float
    q_p_c: float
    x_n_start: float
    x_n_end: float
    x_p_start: float
    x_p_end: float
    r0_ohm: float
    alpha_n_s: float
    alpha_p_s: float
    d_n_per_s: float
    d_p_per_s: float
    c_th_j_per_k: float
    g_th_w_per_k: float
    source_files: tuple[Path, ...]


@dataclass(frozen=True)
class ReferenceDataset:
    role: str
    path: Path


@dataclass(frozen=True)
class AuditBundle:
    classifications: tuple[dict[str, object], ...]
    recommendations: tuple[dict[str, object], ...]
    transfer_decision: TransferDecision
    model_snapshot: PanasonicModelSnapshot
    reference_datasets: tuple[ReferenceDataset, ...]


def _finite_missing_count(experiment: StandardExperiment) -> int:
    return int(
        sum(np.count_nonzero(~np.isfinite(values)) for values in experiment.numeric.values())
    )


def _temperature_channels(experiment: StandardExperiment) -> list[str]:
    return sorted(name for name in experiment.numeric if name.startswith("Temp"))


def _temperature_roles(experiment: StandardExperiment) -> str:
    return "; ".join(
        f"{name}={experiment.metadata.temperature_roles.get(name, 'ubicación no documentada')}"
        for name in _temperature_channels(experiment)
    )


def _format_capacity(value: object) -> float | str:
    number = float(value)
    return number if np.isfinite(number) else "not_available"


def _anomaly_text(
    experiment: StandardExperiment,
    audit: Mapping[str, object],
    nonideality: RelativeNonidealityResult,
) -> str:
    findings: list[str] = []
    if not experiment.header_present:
        findings.append("encabezado ausente; esquema heredado explícitamente")
    findings.append("tiempo medido ausente; base temporal provisional")
    missing = _finite_missing_count(experiment)
    if missing:
        findings.append(f"{missing} valores numéricos no finitos")
    if nonideality.flag == NonidealityFlag.ELEVATED_SERIES_PATH:
        findings.append(
            "respuesta serie aparente elevada respecto a la referencia; atribución física no única"
        )
    elif nonideality.flag == NonidealityFlag.INSUFFICIENT_EVIDENCE:
        findings.append("transiciones insuficientes para clasificar la respuesta serie")
    else:
        findings.append("grupo relativo de referencia; no certifica R0 intrínseco")
    if int(audit["discharge_r_app_count"]) == 1:
        findings.append(
            "una sola transición reposo–descarga; no demuestra repetibilidad dentro del registro"
        )
    unknown_roles = [
        name
        for name in _temperature_channels(experiment)
        if "unknown" in experiment.metadata.temperature_roles.get(name, "unknown")
        or "desconoc" in experiment.metadata.temperature_roles.get(name, "").lower()
    ]
    if unknown_roles:
        findings.append("ubicación no documentada para " + ", ".join(unknown_roles))
    if experiment.test_type == "pulse discharge/rest":
        findings.append("secuencia parcial; no representa capacidad completa")
    if experiment.metadata.dataset_id == "cycle_cell2_undeclared":
        findings.append(
            "Temp7 es compartido con cell1_ambient, sin sincronización muestra a muestra"
        )
    return "; ".join(findings)


def decide_dataset_use(
    *,
    dataset_id: str,
    test_type: str,
    capacity_complete: bool,
    time_measured: bool,
    nonideality_flag: NonidealityFlag,
) -> dict[str, str]:
    """Return conservative presentation-level decisions from audited evidence."""

    is_reference = dataset_id == REFERENCE_DATASET_ID
    elevated = nonideality_flag == NonidealityFlag.ELEVATED_SERIES_PATH
    if is_reference and not elevated and capacity_complete:
        characterization = "SUITABLE WITH CAVEAT"
        recommended = (
            "caracterización de capacidad y diagnóstico de transferencia voltaje–capacidad"
        )
        electrical_validation = "SUITABLE WITH CAVEAT"
    elif test_type == "pulse discharge/rest" or elevated:
        characterization = "QUALITATIVE / DIAGNOSTIC ONLY"
        recommended = "diagnóstico de no idealidad experimental y del montaje"
        electrical_validation = "QUALITATIVE / DIAGNOSTIC ONLY"
    else:
        characterization = "QUALITATIVE / DIAGNOSTIC ONLY"
        recommended = "caracterización descriptiva"
        electrical_validation = "QUALITATIVE / DIAGNOSTIC ONLY"

    electrical_calibration = "NOT SUITABLE"
    thermal_calibration = "NOT SUITABLE"
    thermal_validation = (
        "QUALITATIVE / DIAGNOSTIC ONLY"
        if not time_measured or elevated
        else "SUITABLE WITH CAVEAT"
    )
    return {
        "characterization_status": characterization,
        "electrical_calibration_status": electrical_calibration,
        "electrical_validation_status": electrical_validation,
        "thermal_calibration_status": thermal_calibration,
        "thermal_validation_status": thermal_validation,
        "recommended_use": recommended,
    }


def evaluate_trc_discharge_nonideality(
    audit_rows: list[Mapping[str, object]],
    *,
    reference_dataset: str = REFERENCE_DATASET_ID,
) -> dict[str, RelativeNonidealityResult]:
    """Evaluate only rest-to-discharge R_app for the TRC comparison."""

    results = {
        item.dataset: item
        for item in evaluate_relative_nonideality_gate(
            audit_rows,
            reference_dataset=reference_dataset,
            transition_count_field="discharge_r_app_count",
            median_field="discharge_r_app_median_ohm",
            iqr_field="discharge_r_app_iqr_ohm",
            minimum_transition_count=1,
            metric_label="R_app,desc",
        )
    }
    reference_value = results[reference_dataset].median_r_app_ohm
    if not np.isclose(
        reference_value,
        TRC_REFERENCE_R_APP_DESC_OHM,
        rtol=0.0,
        atol=5.0e-15,
    ):
        raise ValueError(
            "The directional TRC reference no longer matches the verified "
            f"R_app,desc value: {reference_value:.15g} ohm."
        )
    return results


def classify_experiments(
    experiments: Sequence[StandardExperiment],
    *,
    reference_dataset: str = REFERENCE_DATASET_ID,
) -> tuple[dict[str, object], ...]:
    """Audit and classify experiments without modifying their arrays or files."""

    if len({item.metadata.dataset_id for item in experiments}) != len(experiments):
        raise ValueError("Dataset identifiers must be unique.")
    gates: dict[str, DataQualityGateResult] = {}
    audit_rows: list[Mapping[str, object]] = []
    by_id = {item.metadata.dataset_id: item for item in experiments}
    if reference_dataset not in by_id:
        raise ValueError(f"Required reference dataset is missing: {reference_dataset}")
    for experiment in experiments:
        gate = evaluate_data_quality_gate(experiment)
        gates[experiment.metadata.dataset_id] = gate
        audit_rows.append(gate.audit_metrics)
    relative = evaluate_trc_discharge_nonideality(
        audit_rows, reference_dataset=reference_dataset
    )

    rows: list[dict[str, object]] = []
    for experiment in experiments:
        dataset_id = experiment.metadata.dataset_id
        gate = gates[dataset_id]
        audit = gate.audit_metrics
        nonideality = relative[dataset_id]
        decisions = decide_dataset_use(
            dataset_id=dataset_id,
            test_type=experiment.test_type,
            capacity_complete=bool(audit["capacity_sequence_appears_complete"]),
            time_measured=gate.timebase_status.value == "measured",
            nonideality_flag=nonideality.flag,
        )
        rows.append(
            {
                "dataset_id": dataset_id,
                "source_file": str(experiment.metadata.path),
                "cell_id": experiment.metadata.cell_id,
                "thermal_condition": experiment.metadata.nominal_temperature,
                "test_type": experiment.test_type,
                "sample_count": int(experiment.voltage_v.size),
                "header_present": experiment.header_present,
                "available_columns": ";".join(experiment.columns),
                "missing_numeric_values": _finite_missing_count(experiment),
                "measured_time_available": False,
                "timebase_status": gate.timebase_status.value,
                "provisional_dt_s": experiment.estimated_dt_s,
                "provisional_duration_s": float(experiment.estimated_time_s[-1]),
                "raw_current_sign": experiment.raw_current_sign,
                "discharge_capacity_ah": _format_capacity(
                    audit["main_discharge_counter_capacity_ah"]
                ),
                "charge_capacity_ah": _format_capacity(
                    audit["main_charge_counter_capacity_ah"]
                ),
                "capacity_sequence_appears_complete": bool(
                    audit["capacity_sequence_appears_complete"]
                ),
                "initial_condition_status": (
                    "main discharge preceded by recorded rest; independent initial SOC unavailable"
                    if bool(audit["main_discharge_preceded_by_rest"])
                    else "independent initial SOC and preceding rest unavailable"
                ),
                "temperature_channels": ";".join(_temperature_channels(experiment)),
                "temperature_roles": _temperature_roles(experiment),
                "r_app_desc_transition_count": nonideality.transition_count,
                "r_app_desc_median_ohm": nonideality.median_r_app_ohm,
                "r_app_desc_iqr_ohm": nonideality.iqr_r_app_ohm,
                "r_app_desc_ratio_to_reference": nonideality.ratio_to_reference,
                "r_app_desc_interpretation": (
                    "R_app,desc: apparent observed-path rest-to-discharge response; "
                    "not intrinsic cell resistance, isolated setup resistance, or model R0"
                ),
                "nonideality_flag": nonideality.flag.value,
                "anomalies_and_caveats": _anomaly_text(
                    experiment, audit, nonideality
                ),
                **decisions,
                "traceability": (
                    "analysis.experimental_data_audit.audit_experiment; "
                    "analysis.experimental_data_audit.extract_rest_to_load_transitions; "
                    "analysis.battery_dataset_assessment.evaluate_data_quality_gate; "
                    "analysis.battery_dataset_assessment.evaluate_relative_nonideality_gate"
                ),
            }
        )
    return tuple(rows)


def build_transfer_decision() -> TransferDecision:
    """Evaluate the already established Panasonic-to-TRC evidence hierarchy."""

    return evaluate_transfer_gate(
        TransferEvidence(
            model_supplied=True,
            direct_transfer_adequate=False,
            capacity_mapping_compatible=False,
            voltage_reference_compatible=False,
            dynamic_parameters_compatible=False,
            capacity_information_available=True,
            voltage_reference_information_available=True,
            measured_dynamic_time_available=False,
        )
    )


def build_recommendations(
    classifications: Sequence[Mapping[str, object]],
) -> tuple[dict[str, object], ...]:
    """Build a small deterministic action plan from the observed evidence gaps."""

    missing_time = any(not bool(row["measured_time_available"]) for row in classifications)
    elevated = [
        str(row["dataset_id"])
        for row in classifications
        if row["nonideality_flag"] == NonidealityFlag.ELEVATED_SERIES_PATH.value
    ]
    any_calibration_ready = any(
        row["electrical_calibration_status"] != "NOT SUITABLE"
        for row in classifications
    )
    recommendations: list[dict[str, object]] = []

    if missing_time or elevated:
        recommendations.append(
            {
                "priority": 0,
                "action_id": "recover_documentation_and_verify_setup",
                "detected_problem": (
                    "Measured time is absent in all Chroma exports; six datasets show an elevated apparent series-path group."
                ),
                "engineer_action": (
                    "Recover the Chroma recipe/sample clock and document force/sense wiring, voltage-tap locations, contacts, cable path, cell identity, and sensor locations; verify the fixture with a known load."
                ),
                "why_needed": (
                    "Dynamic and thermal identification require measured time, and the observed voltage step cannot be assigned uniquely to cell, fixture, wiring, sensing, or fast polarization."
                ),
                "preliminary_protocol": (
                    "Inventory existing configuration files and logs first. Before another cell test, verify current and voltage channels with a traceable load and confirm whether voltage sensing excludes force-lead drops."
                ),
                "safety_constraint": (
                    "Do not execute a new cell recipe until the TRC18650 datasheet and tester limits are verified."
                ),
                "active_execution_time_h": "not_applicable",
                "rest_time_h": "not_applicable",
                "preparation_time_h": "not_estimable_from_available_records",
                "total_time_status": "not estimable; documentation/access dependent",
                "enables": (
                    "Measured-time audit, defensible R_app,desc interpretation, cross-file synchronization decision, and go/no-go decision for reuse of existing runs."
                ),
                "trigger_evidence": (
                    "timebase_status=inferred/provisional for all datasets; elevated group="
                    + ",".join(elevated)
                ),
            }
        )

    recommendations.extend(
        [
            {
                "priority": 1,
                "action_id": "trc_low_rate_reference",
                "detected_problem": "No TRC18650 low-rate equilibrium-voltage reference is available.",
                "engineer_action": (
                    "Acquire one complete C/20 discharge with measured time and, preferably, a separate C/20 charge using the verified reference fixture."
                ),
                "why_needed": (
                    "The loaded 1C curve cannot identify a cell-compatible OCP/voltage reference or effective electrode windows."
                ),
                "preliminary_protocol": (
                    "Condition the cell under documented manufacturer limits; record current, voltage, Ah, measured time, surface temperature and ambient; retain initial/intermediate/final rests and cutoff events as separate segments. For Q=2.8999 Ah, C/20 is approximately 0.145 A."
                ),
                "safety_constraint": (
                    "Use manufacturer voltage, current and temperature limits; the observed historical cutoffs are not a substitute for a datasheet."
                ),
                "active_execution_time_h": "20 discharge only; 40 for discharge plus charge",
                "rest_time_h": "not estimable until the recipe is fixed",
                "preparation_time_h": "not estimable until conditioning is specified",
                "total_time_status": "active time calculable; total laboratory time not yet established",
                "enables": (
                    "TRC capacity confirmation at low rate, voltage-reference/OCP compatibility gate, and effective state-window adaptation."
                ),
                "trigger_evidence": "no low-rate TRC dataset in the seven-record inventory",
            },
            {
                "priority": 2,
                "action_id": "clean_timed_hppc",
                "detected_problem": (
                    "The existing pulse test is partial, lacks measured time, and belongs to the elevated apparent-response group."
                ),
                "engineer_action": (
                    "After fixture verification and the low-rate gate, acquire a clean timed HPPC profile across documented SOC levels."
                ),
                "why_needed": (
                    "R0, diffusion/relaxation and kinetic groups cannot be separated from setup effects or identified without actual time."
                ),
                "preliminary_protocol": (
                    "Use measured timestamps, SOC set by coulomb counting from a documented initial state, and at least 0.5C, 1C and 2C pulses with approximately 10 s pulse duration and approximately 20 min rest, following the supported Panasonic protocol structure. Reserve pulses/SOC levels for validation."
                ),
                "safety_constraint": (
                    "Current levels require manufacturer approval and tester thermal/voltage limits; omit higher rates unless explicitly supported."
                ),
                "active_execution_time_h": "not fully estimable; pulse active time is minor relative to conditioning",
                "rest_time_h": "8.67 h minimum for 26 verified 20 min rest intervals",
                "preparation_time_h": "not estimable; SOC stepping and conditioning recipe required",
                "total_time_status": (
                    "Panasonic 25 °C reference span is 27.11 h, but it is not a directly transferable TRC duration"
                ),
                "enables": "TRC-specific conditional dynamic identification and held-out pulse validation.",
                "trigger_evidence": (
                    "no electrically calibration-ready dataset; existing pulse median R_app,desc is elevated"
                    if not any_calibration_ready
                    else "existing dynamic evidence remains conditional"
                ),
            },
            {
                "priority": 3,
                "action_id": "clean_timed_thermal_discharge",
                "detected_problem": (
                    "Surface and ambient temperatures exist, but there is no measured time and affected fixtures confound heat generation."
                ),
                "engineer_action": (
                    "First try to recover the clock for the clean ambient cycle; otherwise repeat one 1C discharge and cooling rest with the verified fixture and documented sensors."
                ),
                "why_needed": "Thermal storage and heat-transfer rate cannot be identified in capacity coordinates alone.",
                "preliminary_protocol": (
                    "Record measured time, current, voltage, Ah, central surface temperature, end surface temperature and chamber air continuously through discharge and cooling."
                ),
                "safety_constraint": (
                    "Use manufacturer limits and stop criteria; monitor both cell surface and chamber temperature."
                ),
                "active_execution_time_h": "approximately 1 h at nominal 1C",
                "rest_time_h": "not estimable until a cooling criterion is specified",
                "preparation_time_h": "not estimable from available records",
                "total_time_status": "active time calculable; rest/preparation unresolved",
                "enables": "TRC thermal calibration or validation after electrical heat input is defensible.",
                "trigger_evidence": "all thermal calibration statuses are NOT SUITABLE",
            },
            {
                "priority": 4,
                "action_id": "independent_validation_profile",
                "detected_problem": "No independent TRC profile supports quantitative validation of an adapted model.",
                "engineer_action": (
                    "Reserve a complete cycle or acquire a dynamic profile not used in OCP, HPPC or thermal calibration."
                ),
                "why_needed": "Calibration data cannot also establish independent predictive validity.",
                "preliminary_protocol": (
                    "Start from documented SOC, preserve bidirectional current if present, and record measured time, voltage, current, Ah, surface and ambient temperatures without fitting to the profile."
                ),
                "safety_constraint": "Scale the profile to verified TRC current, voltage and temperature limits.",
                "active_execution_time_h": "6.2355 h for the existing Panasonic UDDS reference only",
                "rest_time_h": "not estimable for a future TRC recipe",
                "preparation_time_h": "not estimable for a future TRC recipe",
                "total_time_status": "reference duration only; TRC total not established",
                "enables": "Independent electrical or electro-thermal validation of the adapted TRC model.",
                "trigger_evidence": "current files are characterization/diagnostic records, not an untouched adapted-model validation set",
            },
        ]
    )
    return tuple(recommendations)


def validate_inputs(
    model_snapshot: PanasonicModelSnapshot,
    reference_datasets: Sequence[ReferenceDataset],
    experiments: Sequence[StandardExperiment],
) -> None:
    missing = [path for path in model_snapshot.source_files if not path.is_file()]
    missing.extend(item.path for item in reference_datasets if not item.path.is_file())
    missing.extend(item.metadata.path for item in experiments if not item.metadata.path.is_file())
    if missing:
        raise FileNotFoundError("Missing audit inputs:\n" + "\n".join(str(path) for path in missing))


def run_audit(
    experiments: Sequence[StandardExperiment],
    model_snapshot: PanasonicModelSnapshot,
    reference_datasets: Sequence[ReferenceDataset],
) -> AuditBundle:
    validate_inputs(model_snapshot, reference_datasets, experiments)
    classifications = classify_experiments(experiments)
    recommendations = build_recommendations(classifications)
    transfer = build_transfer_decision()
    if transfer.automatic_fit_triggered:
        raise RuntimeError("The automatic audit must never trigger parameter fitting.")
    return AuditBundle(
        classifications=classifications,
        recommendations=recommendations,
        transfer_decision=transfer,
        model_snapshot=model_snapshot,
        reference_datasets=tuple(reference_datasets),
    )


def write_csv(path: Path, rows: Iterable[Mapping[str, object]]) -> None:
    materialized = list(rows)
    if not materialized:
        raise ValueError("Cannot write an empty audit CSV.")
    fields: list[str] = []
    for row in materialized:
        for field in row:
            if field not in fields:
                fields.append(field)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(materialized)


def write_engineer_report(path: Path, bundle: AuditBundle) -> None:
    rows = bundle.classifications
    reference = next(row for row in rows if row["dataset_id"] == REFERENCE_DATASET_ID)
    elevated = [
        row
        for row in rows
        if row["nonideality_flag"] == NonidealityFlag.ELEVATED_SERIES_PATH.value
    ]
    model = bundle.model_snapshot

    def capacity_text(value: object) -> str:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return "no disponible"
        return f"{number:.4f}" if np.isfinite(number) else "no disponible"

    test_labels = {
        "combined cycle": "ciclo combinado",
        "pulse discharge/rest": "pulsos de descarga/reposo",
        "discharge": "descarga",
        "charge": "carga",
        "rest": "reposo",
    }
    recommendation_es = {
        "recover_documentation_and_verify_setup": {
            "title": "Recuperar documentación y verificar el montaje",
            "detected_problem": "Los siete CSV carecen de tiempo medido y seis registros pertenecen al grupo de respuesta serie aparente elevada.",
            "engineer_action": "Recuperar la receta y el reloj Chroma; documentar fuerza/sensado, puntos de voltaje, contactos, cables, celda y sensores; verificar el montaje con una carga conocida.",
            "why_needed": "La identificación dinámica y térmica requiere tiempo medido y el salto de voltaje no puede atribuirse de forma única a celda, montaje, cableado, sensado o polarización rápida.",
            "preliminary_protocol": "Inventariar configuraciones y registros existentes. Antes de otra celda, verificar corriente y voltaje con una carga trazable y confirmar si el sensado excluye las caídas en conductores de fuerza.",
            "safety_constraint": "No ejecutar una receta nueva hasta verificar la ficha técnica TRC18650 y los límites del equipo.",
            "enables": "Auditoría temporal, interpretación defendible de R_app,desc, decisión de sincronización y decisión de reutilización de ensayos existentes.",
        },
        "trc_low_rate_reference": {
            "title": "Obtener una referencia TRC de baja tasa",
            "detected_problem": "No existe una referencia TRC18650 de voltaje de equilibrio a baja tasa.",
            "engineer_action": "Adquirir una descarga C/20 completa y, preferiblemente, una carga C/20 separada con tiempo medido y montaje verificado.",
            "why_needed": "La curva cargada a 1C no identifica una referencia OCP/voltaje compatible ni ventanas efectivas de electrodo.",
            "preliminary_protocol": "Acondicionar dentro de límites documentados; registrar corriente, voltaje, Ah, tiempo y temperaturas; conservar reposos y cortes como segmentos. Para 2.8999 Ah, C/20 es aproximadamente 0.145 A.",
            "safety_constraint": "Usar límites de voltaje, corriente y temperatura del fabricante; los cortes históricos no sustituyen una ficha técnica.",
            "enables": "Capacidad TRC a baja tasa, puerta de compatibilidad OCP/voltaje y adaptación de ventanas efectivas.",
        },
        "clean_timed_hppc": {
            "title": "Adquirir un HPPC limpio y temporizado",
            "detected_problem": "El pulso existente es parcial, no tiene tiempo medido y pertenece al grupo de respuesta aparente elevada.",
            "engineer_action": "Después de verificar montaje y baja tasa, adquirir HPPC limpio en niveles de SOC documentados.",
            "why_needed": "R0 y los grupos de difusión, relajación y cinética no pueden separarse del montaje ni identificarse sin tiempo real.",
            "preliminary_protocol": "Usar tiempo medido, SOC por conteo desde un estado documentado y pulsos de 0.5C, 1C y 2C de unos 10 s, con reposos cercanos a 20 min. Reservar pulsos o SOC para validación.",
            "safety_constraint": "Los niveles de corriente requieren soporte del fabricante y límites térmicos/eléctricos verificados; omitir tasas mayores sin soporte explícito.",
            "enables": "Identificación dinámica condicional TRC y validación con pulsos reservados.",
        },
        "clean_timed_thermal_discharge": {
            "title": "Obtener una descarga térmica limpia y temporizada",
            "detected_problem": "Hay temperaturas superficial y ambiente, pero falta tiempo medido y los montajes afectados confunden el calor.",
            "engineer_action": "Intentar recuperar el reloj del ciclo ambiente; si no es posible, repetir una descarga 1C y enfriamiento con montaje y sensores verificados.",
            "why_needed": "La capacidad térmica y la transferencia de calor no pueden identificarse solo en coordenadas de capacidad.",
            "preliminary_protocol": "Registrar tiempo, corriente, voltaje, Ah, temperatura superficial central, temperatura cerca del extremo y aire durante descarga y enfriamiento.",
            "safety_constraint": "Usar límites y criterios de parada del fabricante; vigilar superficie y aire de cámara.",
            "enables": "Calibración o validación térmica TRC cuando el calor eléctrico sea defendible.",
        },
        "independent_validation_profile": {
            "title": "Reservar un perfil para validación independiente",
            "detected_problem": "No existe un perfil TRC independiente para validar cuantitativamente un modelo adaptado.",
            "engineer_action": "Reservar un ciclo o adquirir un perfil dinámico no usado en OCP, HPPC ni calibración térmica.",
            "why_needed": "Los datos de calibración no pueden establecer también validez predictiva independiente.",
            "preliminary_protocol": "Partir de SOC documentado, conservar corriente bidireccional si existe y registrar tiempo, voltaje, corriente, Ah y temperaturas sin ajustar al perfil.",
            "safety_constraint": "Escalar el perfil a límites TRC verificados de corriente, voltaje y temperatura.",
            "enables": "Validación eléctrica o electro-térmica independiente del modelo TRC adaptado.",
        },
    }

    def time_text(value: object) -> str:
        translations = {
            "not_applicable": "no aplica",
            "not_estimable_from_available_records": "no estimable con los registros disponibles",
            "20 discharge only; 40 for discharge plus charge": "20 h para descarga; 40 h para descarga y carga",
            "not estimable until the recipe is fixed": "no estimable hasta fijar la receta",
            "not estimable until conditioning is specified": "no estimable hasta definir el acondicionamiento",
            "not fully estimable; pulse active time is minor relative to conditioning": "no estimable completamente; los pulsos son breves frente al acondicionamiento",
            "8.67 h minimum for 26 verified 20 min rest intervals": "mínimo 8.67 h para 26 reposos verificados de 20 min",
            "not estimable; SOC stepping and conditioning recipe required": "no estimable; requiere receta de acondicionamiento y cambio de SOC",
            "approximately 1 h at nominal 1C": "aproximadamente 1 h a 1C nominal",
            "not estimable until a cooling criterion is specified": "no estimable hasta definir el criterio de enfriamiento",
            "not estimable from available records": "no estimable con los registros disponibles",
            "6.2355 h for the existing Panasonic UDDS reference only": "6.2355 h solo para la referencia UDDS Panasonic existente",
            "not estimable for a future TRC recipe": "no estimable para una futura receta TRC",
        }
        return translations.get(str(value), str(value))
    lines = [
        "# Informe automático de auditoría TRC18650",
        "",
        "Este informe aplica reglas deterministas de disponibilidad, calidad, no idealidad experimental y transferencia. No ejecuta calibración, no corrige voltajes y no interpreta `R_app,desc` como resistencia interna de la celda o como `R0`.",
        "",
        "## A. Resumen del estado actual",
        "",
        f"- Registros TRC auditados: **{len(rows)}**.",
        f"- Referencia experimental relativa reconocida: **`{REFERENCE_DATASET_ID}`**, capacidad de descarga `{float(reference['discharge_capacity_ah']):.4f} Ah`.",
        f"- Registros en el grupo de respuesta serie aparente elevada: **{len(elevated)} de {len(rows)}**.",
        "- Tiempo medido: **ausente en los siete CSV**. Las duraciones son provisionales y no autorizan calibración dinámica o térmica.",
        "- Referencia TRC de baja tasa: **no disponible**.",
        "- Comparación Panasonic→TRC: permitida únicamente como diagnóstico de transferencia, no como validación de un modelo TRC adaptado.",
        "",
        "Entradas Panasonic verificadas:",
        "",
        f"- `Q_cell={model.q_cell_ah:.12f} Ah`, `Q_n={model.q_n_c:.12f} C`, `Q_p={model.q_p_c:.12f} C`.",
        f"- Dinámica condicional: `R0={model.r0_ohm:.13g} Ω`, `alpha_p={model.alpha_p_s:.12g} s`, `d_n={model.d_n_per_s:.12g} s⁻¹`; `alpha_n={model.alpha_n_s:.12g} s`, `d_p={model.d_p_per_s:.12g} s⁻¹`.",
        f"- Térmica efectiva: `Cth={model.c_th_j_per_k:.12g} J/K`, `Gth={model.g_th_w_per_k:.12g} W/K`.",
        "",
        "## B. Diagnóstico experimental",
        "",
        "| Ensayo | Condición | Tipo | Q descarga/carga [Ah] | Tiempo | R_app,desc [Ω] | N descarga | Relación | Aptitud principal |",
        "|---|---|---|---:|---|---:|---:|---:|---|",
    ]
    for row in rows:
        capacities = (
            f"{capacity_text(row['discharge_capacity_ah'])}/"
            f"{capacity_text(row['charge_capacity_ah'])}"
        )
        lines.append(
            f"| `{row['dataset_id']}` | {row['thermal_condition']} | {test_labels.get(str(row['test_type']), row['test_type'])} | {capacities} | provisional inferida | {float(row['r_app_desc_median_ohm']):.6f} | {int(row['r_app_desc_transition_count'])} | {float(row['r_app_desc_ratio_to_reference']):.2f}× | {row['recommended_use']} |"
        )
    lines.extend(
        [
            "",
            "`R_app,desc` utiliza exclusivamente transiciones reposo → descarga. Los ciclos completos contienen una sola transición de descarga y, por tanto, no demuestran repetibilidad dentro de cada registro; los archivos de pulsos contienen siete y tres transiciones. La separación es evidencia de una respuesta elevada en el camino observado, pero no localiza la causa en contactos, cables, sensado, celda, instrumento o polarización rápida. Los registros eléctricos y térmicos de archivos distintos no se combinan sin una base temporal común.",
            "",
            "## C. Estado de transferencia Panasonic → TRC18650",
            "",
            f"La puerta jerárquica devuelve: **{bundle.transfer_decision.summary}**.",
            "",
            "| Componente | Estado para TRC |",
            "|---|---|",
            "| Capacidad de celda | Disponible: 2.8999 Ah en la referencia ambiente; debe sustituir la capacidad Panasonic para TRC. |",
            "| OCP y ventanas efectivas | No identificadas para TRC; la curva 1C cargada no es OCV. |",
            "| Parámetros dinámicos | No identificables con pulsos afectados y tiempo provisional. |",
            "| Parámetros térmicos | No identificables sin tiempo medido y sin partición defendible del calor. |",
            "| Comparación directa Panasonic | Diagnóstico de transferencia solamente. |",
            "| Validación cuantitativa TRC | No autorizada con la evidencia actual. |",
            "",
            "## D. Plan experimental recomendado",
            "",
        ]
    )
    for item in bundle.recommendations:
        display = recommendation_es[str(item["action_id"])]
        lines.extend(
            [
                f"### Prioridad {item['priority']} — {display['title']}",
                "",
                f"- **Carencia:** {display['detected_problem']}",
                f"- **Acción:** {display['engineer_action']}",
                f"- **Justificación:** {display['why_needed']}",
                f"- **Receta preliminar:** {display['preliminary_protocol']}",
                f"- **Seguridad:** {display['safety_constraint']}",
                f"- **Tiempo activo:** {time_text(item['active_execution_time_h'])}.",
                f"- **Reposos:** {time_text(item['rest_time_h'])}.",
                f"- **Preparación:** {time_text(item['preparation_time_h'])}.",
                f"- **Resultado habilitado:** {display['enables']}",
                "",
            ]
        )
    lines.extend(
        [
            "## E. Estado final",
            "",
            "- **¿Es posible caracterizar la celda con los datos actuales?** Sí, parcialmente: capacidad de la referencia ambiente, curvas voltaje–capacidad, señales térmicas descriptivas y diagnóstico relativo del camino serie.",
            "- **¿Es posible realizar una comparación diagnóstica con el modelo Panasonic?** Sí. Debe rotularse como diagnóstico de transferencia y no como precisión de un modelo TRC.",
            "- **¿Es posible validar cuantitativamente un modelo adaptado a la TRC18650?** No.",
            "- **¿Qué falta?** Tiempo medido, trazabilidad del montaje y sensado, referencia TRC de baja tasa, pulsos limpios, evidencia térmica temporizada y un perfil reservado para validación.",
            "- **¿Cuál es la siguiente acción recomendada?** Recuperar primero la receta/reloj Chroma y documentar/verificar el camino de fuerza y sensado antes de repetir ensayos.",
            "",
            "## Trazabilidad",
            "",
            "- Calidad y segmentación: `analysis/experimental_data_audit.py`.",
            "- Puertas de calidad/no idealidad/transferencia: `analysis/battery_dataset_assessment.py`.",
            "- Configuración de los siete registros: `scripts/audit_own_trc18650_data.py`.",
            "- Diagnóstico directo previo: `scripts/audit_panasonic_to_trc18650_transfer.py` y `thesis_data/own_model_transfer_metrics.csv`.",
            "- La clasificación completa y las recomendaciones se exportan en los CSV asociados a este informe.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


__all__ = [
    "AuditBundle",
    "PanasonicModelSnapshot",
    "REFERENCE_DATASET_ID",
    "ReferenceDataset",
    "TRC_REFERENCE_CAPACITY_AH",
    "build_recommendations",
    "build_transfer_decision",
    "classify_experiments",
    "decide_dataset_use",
    "run_audit",
    "validate_inputs",
    "write_csv",
    "write_engineer_report",
]
