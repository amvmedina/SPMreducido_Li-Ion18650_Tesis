"""Application adapters for the brand-independent battery-data audit."""

from __future__ import annotations

import csv
import io
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from analysis.battery_audit_rules import (
    TIME_CONFIRMED_INTERVAL,
    TIME_ESTIMATED,
    TIME_RECORDED,
    CommonAuditRecord,
    RuleOutcome,
    evaluate_common_audit,
)
from analysis.battery_dataset_assessment import NonidealityFlag, evaluate_data_quality_gate, evaluate_relative_nonideality_gate
from analysis.experimental_data_audit import ExperimentMetadata, StandardExperiment, load_chroma_experiment
from scripts.audit_own_trc18650_data import PART2_COLUMNS, experiment_metadata
from scripts.run_trc18650_automatic_audit import model_snapshot, reference_datasets
from runtime_paths import require_data_file


STATUS_CONFORME = "CONFORME"
STATUS_REVISAR = "REVISAR"
STATUS_FALTANTE = "FALTANTE"
STATUS_NO_APTO = "NO APTO"
STATUS_NO_EVALUADO = "NO EVALUADO"

SETUP_REFERENCE = "montaje de referencia"
SETUP_OTHER = "otro montaje"
SETUP_UNKNOWN = "desconocido"

ROLE_LOW_RATE = "referencia de baja tasa"
ROLE_DYNAMIC = "pulsos para identificación dinámica"
ROLE_THERMAL = "identificación térmica"
ROLE_VALIDATION = "validación independiente"
ROLE_CHARACTERIZATION = "caracterización"
ROLE_UNKNOWN = "desconocido"

STANDARD_CHROMA_COLUMNS = PART2_COLUMNS

AUDIT_CRITERIA = (
    {"código": "E01", "nombre": "Variables eléctricas", "qué verifica": "Disponibilidad de voltaje, corriente y capacidad."},
    {"código": "E02", "nombre": "Capacidad experimental", "qué verifica": "Disponibilidad de una carga o descarga completa."},
    {"código": "E03", "nombre": "Integridad de los datos", "qué verifica": "Ausencia de datos faltantes o inválidos."},
    {"código": "E04", "nombre": "Montaje eléctrico", "qué verifica": "Comparación de la resistencia aparente entre montajes."},
    {"código": "E05", "nombre": "Ensayo de baja corriente", "qué verifica": "Disponibilidad de una descarga de caracterización a baja corriente."},
    {"código": "E06", "nombre": "Ensayo de pulsos", "qué verifica": "Disponibilidad de pulsos de corriente y reposos."},
    {"código": "E07", "nombre": "Información térmica", "qué verifica": "Disponibilidad y aptitud de las mediciones de temperatura."},
    {"código": "E08", "nombre": "Validación independiente", "qué verifica": "Disponibilidad de un ensayo no utilizado para calibración."},
)


@dataclass(frozen=True)
class CellContext:
    chemistry: str = "desconocida"
    nominal_capacity_ah: float | None = None
    operating_limits: str = "desconocidos"
    discharge_cutoff_v: float | None = None
    charge_cutoff_v: float | None = None

    @property
    def safety_documented(self) -> bool:
        return self.operating_limits.strip().lower() not in {"", "desconocido", "desconocidos", "unknown"}


@dataclass(frozen=True)
class UploadedFileSpec:
    filename: str
    content: bytes
    dataset_id: str
    cell_id: str = "desconocida"
    condition: str = "desconocida"
    setup_role: str = SETUP_UNKNOWN
    experimental_roles: tuple[str, ...] = (ROLE_UNKNOWN,)
    temperature_roles: Mapping[str, str] | None = None
    fallback_columns: tuple[str, ...] | None = None
    include_in_rapp_comparison: bool = True
    sampling_interval_s: float | None = None
    sampling_interval_constant_confirmed: bool = False
    record_continuous_confirmed: bool = False


@dataclass(frozen=True)
class ApplicationAuditResult:
    mode: str
    file_summaries: tuple[dict[str, object], ...]
    checks: tuple[dict[str, object], ...]
    comparison_rows: tuple[dict[str, object], ...]
    report_markdown: str

    def checks_csv(self) -> bytes:
        return rows_to_csv_bytes(self.checks)

    def summaries_csv(self) -> bytes:
        return rows_to_csv_bytes(self.file_summaries)


def rows_to_csv_bytes(rows: Iterable[Mapping[str, object]]) -> bytes:
    materialized = list(rows)
    if not materialized:
        return b""
    fields: list[str] = []
    for row in materialized:
        for key in row:
            if key not in fields:
                fields.append(key)
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields)
    writer.writeheader()
    writer.writerows(materialized)
    return buffer.getvalue().encode("utf-8-sig")


def audit_criteria_payload() -> tuple[dict[str, str], ...]:
    return tuple(dict(item) for item in AUDIT_CRITERIA)


def _decode_csv(content: bytes) -> str:
    for encoding in ("utf-8-sig", "utf-8", "cp1252"):
        try:
            return content.decode(encoding)
        except UnicodeDecodeError:
            continue
    raise ValueError("El archivo no puede decodificarse como UTF-8 o Windows-1252.")


def inspect_csv_bytes(filename: str, content: bytes) -> dict[str, object]:
    text = _decode_csv(content)
    rows = list(csv.reader(io.StringIO(text)))
    required = {"Voltage(V)", "Current(A)", "Capacity(Ah)"}
    header_index = next((index for index, row in enumerate(rows) if required.issubset(set(row))), None)
    columns = tuple(rows[header_index]) if header_index is not None else tuple()
    return {
        "filename": Path(filename).name,
        "header_detected": header_index is not None,
        "columns": columns,
        "column_count": len(columns),
        "sample_count_estimate": max(0, len(rows) - header_index - 1) if header_index is not None else len(rows),
        "temperature_columns": tuple(name for name in columns if re.fullmatch(r"Temp\d+", name)),
        "required_electrical_columns_detected": required.issubset(set(columns)),
    }


def reference_model_payload() -> dict[str, tuple[dict[str, object], ...]]:
    snapshot = model_snapshot()
    parameters = (
        {"Parámetro": "Q_celda", "Descripción": "Capacidad de la celda usada para el balance de carga", "Valor": snapshot.q_cell_ah, "Unidad": "Ah"},
        {"Parámetro": "Qₙ", "Descripción": "Escala efectiva de capacidad del electrodo negativo", "Valor": snapshot.q_n_c, "Unidad": "C"},
        {"Parámetro": "Qₚ", "Descripción": "Escala efectiva de capacidad del electrodo positivo", "Valor": snapshot.q_p_c, "Unidad": "C"},
        {"Parámetro": "xₙ,ini", "Descripción": "Estado estequiométrico efectivo inicial del electrodo negativo", "Valor": snapshot.x_n_start, "Unidad": "–"},
        {"Parámetro": "xₙ,fin", "Descripción": "Estado estequiométrico efectivo final del electrodo negativo", "Valor": snapshot.x_n_end, "Unidad": "–"},
        {"Parámetro": "xₚ,ini", "Descripción": "Estado estequiométrico efectivo inicial del electrodo positivo", "Valor": snapshot.x_p_start, "Unidad": "–"},
        {"Parámetro": "xₚ,fin", "Descripción": "Estado estequiométrico efectivo final del electrodo positivo", "Valor": snapshot.x_p_end, "Unidad": "–"},
        {"Parámetro": "R₀", "Descripción": "Resistencia óhmica efectiva en la ecuación terminal", "Valor": snapshot.r0_ohm, "Unidad": "Ω"},
        {"Parámetro": "αₙ", "Descripción": "Escala temporal efectiva de difusión del electrodo negativo", "Valor": snapshot.alpha_n_s, "Unidad": "s"},
        {"Parámetro": "αₚ", "Descripción": "Escala temporal efectiva de difusión del electrodo positivo", "Valor": snapshot.alpha_p_s, "Unidad": "s"},
        {"Parámetro": "dₙ", "Descripción": "Grupo cinético efectivo del electrodo negativo", "Valor": snapshot.d_n_per_s, "Unidad": "s⁻¹"},
        {"Parámetro": "dₚ", "Descripción": "Grupo cinético efectivo del electrodo positivo", "Valor": snapshot.d_p_per_s, "Unidad": "s⁻¹"},
        {"Parámetro": "Cₜₕ", "Descripción": "Capacitancia térmica efectiva de la celda", "Valor": snapshot.c_th_j_per_k, "Unidad": "J/K"},
        {"Parámetro": "Gₜₕ", "Descripción": "Conductancia térmica efectiva hacia el ambiente", "Valor": snapshot.g_th_w_per_k, "Unidad": "W/K"},
    )
    tests = tuple({"función": item.role, "archivo": item.path.name, "ruta": str(item.path), "disponible": item.path.is_file()} for item in reference_datasets())
    sources = tuple({"archivo_fuente": str(path), "disponible": path.is_file()} for path in snapshot.source_files)
    return {"parameters": parameters, "reference_tests": tests, "sources": sources}


def _safe_dataset_id(value: str, index: int) -> str:
    normalized = re.sub(r"[^A-Za-z0-9_]+", "_", value.strip()).strip("_")
    return normalized or f"dataset_{index + 1}"


def _capacity_text(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "no disponible"
    return f"{number:.6g} Ah" if np.isfinite(number) else "no disponible"


def _check(outcome: RuleOutcome) -> dict[str, object]:
    criterion = next(item for item in AUDIT_CRITERIA if item["código"] == outcome.code)
    return {
        "conjunto": "conjunto completo",
        "código": outcome.code,
        "requisito evaluado": criterion["nombre"],
        "qué verifica": criterion["qué verifica"],
        "estado": outcome.status,
        "evidencia": outcome.evidence,
        "acción requerida": outcome.action,
    }


def _generic_report(summaries: Sequence[Mapping[str, object]], checks: Sequence[Mapping[str, object]], context: CellContext) -> str:
    counts = {status: sum(row["estado"] == status for row in checks) for status in (STATUS_CONFORME, STATUS_REVISAR, STATUS_FALTANTE, STATUS_NO_APTO, STATUS_NO_EVALUADO)}
    next_action = next((str(row["acción requerida"]) for row in checks if row["estado"] not in {STATUS_CONFORME, STATUS_NO_EVALUADO} and row["acción requerida"]), "No hay acciones pendientes.")
    lines = [
        "# Informe técnico de auditoría de datos de batería", "",
        "El informe se genera mediante reglas determinísticas aplicadas a una estructura experimental común. Los datos no se modifican.", "",
        "## Contexto declarado", "",
        f"- Química: **{context.chemistry or 'desconocida'}**.",
        f"- Capacidad nominal: **{context.nominal_capacity_ah if context.nominal_capacity_ah is not None else 'desconocida'}**.",
        f"- Límites de operación: **{context.operating_limits or 'desconocidos'}**.", "",
        "## Resumen", "",
        f"- CONFORME: **{counts[STATUS_CONFORME]}**.", f"- REVISAR: **{counts[STATUS_REVISAR]}**.",
        f"- FALTANTE: **{counts[STATUS_FALTANTE]}**.", f"- NO APTO: **{counts[STATUS_NO_APTO]}**.",
        f"- NO EVALUADO: **{counts[STATUS_NO_EVALUADO]}**.", f"- Siguiente acción: **{next_action}**.", "",
        "## Verificaciones", "", "| Código | Verificación | Estado | Resultado | Acción |", "|---|---|---|---|---|",
    ]
    for row in checks:
        values = [str(row[key]).replace("|", "/") for key in ("código", "requisito evaluado", "estado", "evidencia", "acción requerida")]
        lines.append("| " + " | ".join(values) + " |")
    lines.extend(["", "La resistencia aparente reposo–descarga es un diagnóstico relativo y no identifica por sí sola una causa física."])
    return "\n".join(lines) + "\n"


def _temperature_arrays(experiment: StandardExperiment, spec: UploadedFileSpec) -> tuple[np.ndarray | None, np.ndarray | None]:
    surface = [name for name, role in (spec.temperature_roles or {}).items() if "surface" in role and name in experiment.numeric]
    ambient = [name for name, role in (spec.temperature_roles or {}).items() if ("ambient" in role or "chamber" in role) and name in experiment.numeric]
    return (
        np.asarray(experiment.numeric[surface[0]], dtype=float) if surface else None,
        np.asarray(experiment.numeric[ambient[0]], dtype=float) if ambient else None,
    )


def _summary_nonideality(loaded: Mapping[str, StandardExperiment], specs: Mapping[str, UploadedFileSpec], gates: Mapping[str, object]) -> dict[str, object]:
    references = [dataset_id for dataset_id, spec in specs.items() if spec.setup_role == SETUP_REFERENCE and spec.include_in_rapp_comparison and dataset_id in loaded]
    finite = [dataset_id for dataset_id, gate in gates.items() if specs[dataset_id].include_in_rapp_comparison and np.isfinite(float(gate.audit_metrics["discharge_r_app_median_ohm"]))]
    if len(references) != 1 or len(finite) < 2:
        return {}
    try:
        return {item.dataset: item for item in evaluate_relative_nonideality_gate(
            [gates[item].audit_metrics for item in finite], references[0],
            transition_count_field="discharge_r_app_count", median_field="discharge_r_app_median_ohm",
            iqr_field="discharge_r_app_iqr_ohm", minimum_transition_count=1, metric_label="R_app,desc",
        )}
    except ValueError:
        return {}


def audit_uploaded_files(specs: Sequence[UploadedFileSpec], context: CellContext | None = None) -> ApplicationAuditResult:
    context = context or CellContext()
    if not specs:
        raise ValueError("Debe cargarse al menos un archivo CSV.")
    dataset_ids = [_safe_dataset_id(spec.dataset_id, index) for index, spec in enumerate(specs)]
    if len(set(dataset_ids)) != len(dataset_ids):
        raise ValueError("Los identificadores de conjunto deben ser únicos.")
    loaded: dict[str, StandardExperiment] = {}
    gates: dict[str, object] = {}
    errors: dict[str, str] = {}
    spec_by_id = dict(zip(dataset_ids, specs, strict=True))
    with tempfile.TemporaryDirectory(prefix="battery_audit_upload_") as directory:
        root = Path(directory)
        for index, (dataset_id, spec) in enumerate(zip(dataset_ids, specs, strict=True)):
            path = root / f"{index:03d}_{Path(spec.filename).name}"
            path.write_bytes(spec.content)
            metadata = ExperimentMetadata(dataset_id, path, spec.cell_id, spec.condition, "no suministrada", dict(spec.temperature_roles or {}))
            try:
                experiment = load_chroma_experiment(metadata, spec.fallback_columns)
                loaded[dataset_id] = experiment
                gates[dataset_id] = evaluate_data_quality_gate(experiment)
            except Exception as exc:
                errors[dataset_id] = f"{type(exc).__name__}: {exc}"

        relative = _summary_nonideality(loaded, spec_by_id, gates)
        summaries: list[dict[str, object]] = []
        common: list[CommonAuditRecord] = []
        for dataset_id, spec in spec_by_id.items():
            name = f"{spec.condition} · {spec.setup_role}"
            if dataset_id not in loaded:
                summaries.append({"dataset_id": dataset_id, "nombre visible": name, "estado": "no evaluable", "uso recomendado": "corregir el archivo"})
                continue
            experiment = loaded[dataset_id]
            gate = gates[dataset_id]
            audit = gate.audit_metrics
            response = relative.get(dataset_id)
            ratio: float | str = response.ratio_to_reference if response is not None else "no disponible"
            rapp = float(audit["discharge_r_app_median_ohm"])
            summaries.append({
                "dataset_id": dataset_id, "nombre visible": name, "archivo": Path(spec.filename).name,
                "condición": spec.condition, "montaje declarado": spec.setup_role, "muestras": int(experiment.voltage_v.size),
                "tipo de ensayo": experiment.test_type, "capacidad descarga": _capacity_text(audit["main_discharge_counter_capacity_ah"]),
                "capacidad carga": _capacity_text(audit["main_charge_counter_capacity_ah"]),
                "R_app,desc [ohm]": rapp if np.isfinite(rapp) else "no disponible", "razón relativa": ratio,
                "uso recomendado": "diagnóstico experimental según las verificaciones E01–E08",
            })
            cell_temp, ambient_temp = _temperature_arrays(experiment, spec)
            if "Time" in experiment.numeric:
                audit_time = np.asarray(experiment.numeric["Time"], dtype=float)
                measured_time = True
                time_basis = TIME_RECORDED
            elif (
                spec.sampling_interval_s is not None
                and spec.sampling_interval_s > 0.0
                and spec.sampling_interval_constant_confirmed
                and spec.record_continuous_confirmed
            ):
                audit_time = np.arange(experiment.voltage_v.size, dtype=float) * spec.sampling_interval_s
                measured_time = False
                time_basis = TIME_CONFIRMED_INTERVAL
            else:
                interval = spec.sampling_interval_s if spec.sampling_interval_s is not None and spec.sampling_interval_s > 0.0 else experiment.estimated_dt_s
                audit_time = np.arange(experiment.voltage_v.size, dtype=float) * interval
                measured_time = False
                time_basis = TIME_ESTIMATED
            common.append(CommonAuditRecord(
                record_id=dataset_id, name=name, voltage_v=np.asarray(experiment.voltage_v), current_a=np.asarray(experiment.raw_current_a),
                capacity_ah=np.asarray(experiment.numeric["Capacity(Ah)"]), roles=spec.experimental_roles, setup_role=spec.setup_role,
                include_in_mount_comparison=spec.include_in_rapp_comparison, time_s=audit_time, measured_time=measured_time,
                cell_temperature_c=cell_temp, ambient_temperature_c=ambient_temp,
                discharge_cutoff_v=context.discharge_cutoff_v, charge_cutoff_v=context.charge_cutoff_v,
                capacity_test_eligible=experiment.test_type != "pulse discharge/rest",
                pulse_candidate=experiment.test_type == "pulse discharge/rest",
                time_basis=time_basis,
            ))
    outcomes = evaluate_common_audit(tuple(common), expected_record_count=len(specs),
                                     load_error_names=tuple(f"{spec_by_id[item].condition} · {spec_by_id[item].setup_role}" for item in errors),
                                     nominal_capacity_ah=context.nominal_capacity_ah)
    checks = tuple(_check(outcome) for outcome in outcomes)
    return ApplicationAuditResult("archivos cargados", tuple(summaries), checks, tuple(), _generic_report(summaries, checks, context))


def build_trc_case_study(root: Path | None = None) -> ApplicationAuditResult:
    del root
    metadata = {item.dataset_id: item for item in experiment_metadata()}
    selections = (
        ("cycle_cell1_ambient", "rionegro_referencia", "Rionegro · Aire 22–24 °C", SETUP_REFERENCE, (ROLE_CHARACTERIZATION, ROLE_THERMAL), True),
        ("cycle_cell2_undeclared", "rionegro_resistencia_externa", "Rionegro · Aire 22–24 °C", "montaje con resistencia externa", (ROLE_CHARACTERIZATION,), True),
        ("cycle_cell2_30c", "cali_electrico", "Cali · Aire 29–31 °C", "montaje con resistencia externa", (ROLE_CHARACTERIZATION,), False),
        ("cycle_cell1_30c", "cali_termico_valledupar_descarga", "Cali / Valledupar · procedencia separada", "montaje con resistencia externa", (ROLE_CHARACTERIZATION, ROLE_THERMAL), False),
        ("cycle_cell1_35c", "valledupar_carga_termica", "Valledupar · Aire 33–37 °C", "montaje con resistencia externa", (ROLE_CHARACTERIZATION, ROLE_THERMAL), False),
    )
    specs = [UploadedFileSpec(metadata[source].path.name, require_data_file(metadata[source].path).read_bytes(), public, "identificador conservado", condition, setup, roles, metadata[source].temperature_roles, None, include)
             for source, public, condition, setup, roles, include in selections]
    calculated = audit_uploaded_files(specs, CellContext("no documentada", 2.8999, "no documentados"))
    visible = (
        {"condición": "Rionegro · Aire 22–24 °C", "ensayo": "Descarga", "montaje": "Referencia", "variables disponibles": "Voltaje, corriente y capacidad"},
        {"condición": "Rionegro · Aire 22–24 °C", "ensayo": "Descarga", "montaje": "Resistencia externa", "variables disponibles": "Voltaje, corriente y capacidad"},
        {"condición": "Rionegro · Aire 22–24 °C", "ensayo": "Carga CC–CV", "montaje": "Referencia", "variables disponibles": "Voltaje, corriente, capacidad y temperatura"},
        {"condición": "Rionegro · Aire 22–24 °C", "ensayo": "Carga CC–CV", "montaje": "Resistencia externa", "variables disponibles": "Voltaje, corriente y capacidad"},
        {"condición": "Cali · Aire 29–31 °C", "ensayo": "Descarga", "montaje": "Resistencia externa", "variables disponibles": "Voltaje, corriente y capacidad"},
        {"condición": "Cali · Aire 29–31 °C", "ensayo": "Carga CC–CV", "montaje": "Resistencia externa", "variables disponibles": "Voltaje, corriente y capacidad"},
        {"condición": "Valledupar · Aire 33–37 °C", "ensayo": "Descarga", "montaje": "Resistencia externa", "variables disponibles": "Voltaje y capacidad"},
        {"condición": "Valledupar · Aire 33–37 °C", "ensayo": "Carga CC–CV", "montaje": "Resistencia externa", "variables disponibles": "Voltaje, corriente, capacidad y temperatura"},
    )
    return ApplicationAuditResult("caso TRC18650", visible, calculated.checks, tuple(), calculated.report_markdown)


def _panasonic_common_records() -> tuple[CommonAuditRecord, ...]:
    from scripts.audit_panasonic_hppc import load_hppc
    from scripts.validate_panasonic_1c_forward import DISCHARGE_CURRENT_THRESHOLD_A, load_raw_discharge, longest_contiguous_true
    from scripts.validate_panasonic_cccv_charge import RAW_SOURCE_PATH as CCCV_PATH, load_mat_channels, remove_duplicate_times
    from scripts.validate_panasonic_udds_forward import load_udds, remove_exact_duplicate_times
    from scripts.verify_panasonic_c20_qi_frozen import DATA_PATH, DISCHARGE_END_INDEX_INCLUSIVE, DISCHARGE_START_INDEX

    table = np.genfromtxt(require_data_file(DATA_PATH), delimiter=",", names=True, dtype=float, encoding="utf-8")
    selection = slice(DISCHARGE_START_INDEX, DISCHARGE_END_INDEX_INCLUSIVE + 1)
    c20_time = np.asarray(table["Time"][selection]); c20_current = np.asarray(table["Current"][selection]); c20_voltage = np.asarray(table["Voltage"][selection])
    c20_capacity = np.concatenate(([0.0], np.cumsum(0.5 * (np.abs(c20_current[:-1]) + np.abs(c20_current[1:])) * np.diff(c20_time) / 3600.0)))

    hppc = load_hppc()
    _, raw_1c = load_raw_discharge(); idx = longest_contiguous_true(raw_1c["Current"] <= DISCHARGE_CURRENT_THRESHOLD_A)
    _, raw_cccv = load_mat_channels(CCCV_PATH); cccv, _ = remove_duplicate_times(raw_cccv)
    _, raw_udds = load_udds(); udds, _ = remove_exact_duplicate_times(raw_udds)
    return (
        CommonAuditRecord("panasonic_c20", "C/20", c20_voltage, c20_current, c20_capacity, (ROLE_LOW_RATE,), SETUP_UNKNOWN, False, c20_time, True, discharge_cutoff_v=2.5),
        CommonAuditRecord("panasonic_hppc", "HPPC", hppc.voltage_v, hppc.current_a, hppc.exported_ah, (ROLE_DYNAMIC,), SETUP_UNKNOWN, False, hppc.time_s, True, hppc.battery_temperature_c, hppc.chamber_temperature_c, 2.5, 4.2, False, True),
        CommonAuditRecord("panasonic_1c", "Descarga a 1C", raw_1c["Voltage"][idx], raw_1c["Current"][idx], raw_1c["Ah"][idx], (ROLE_THERMAL, ROLE_VALIDATION), SETUP_UNKNOWN, False, raw_1c["Time"][idx], True, raw_1c["Battery_Temp_degC"][idx], raw_1c["Chamber_Temp_degC"][idx], 2.5, None),
        CommonAuditRecord("panasonic_cccv", "Carga CC–CV", cccv["Voltage"], cccv["Current"], cccv["Ah"], (ROLE_VALIDATION,), SETUP_UNKNOWN, False, cccv["Time"], True, cccv["Battery_Temp_degC"], cccv["Chamber_Temp_degC"], None, 4.2),
        CommonAuditRecord("panasonic_udds", "UDDS", udds["Voltage"], udds["Current"], udds["Ah"], (ROLE_VALIDATION,), SETUP_UNKNOWN, False, udds["Time"], True, udds["Battery_Temp_degC"], udds["Chamber_Temp_degC"], 2.5, 4.2, False, False),
    )


def build_panasonic_reference_audit() -> ApplicationAuditResult:
    records = _panasonic_common_records()
    outcomes = evaluate_common_audit(records, expected_record_count=len(records), nominal_capacity_ah=2.994979138519)
    checks = tuple(_check(outcome) for outcome in outcomes)
    summaries = tuple({
        "nombre visible": record.name,
        "función experimental": ", ".join(record.roles),
        "variables disponibles": "Voltaje, corriente, capacidad" + (", temperatura de celda y ambiente" if record.cell_temperature_c is not None and record.ambient_temperature_c is not None else ""),
        "uso recomendado": "evaluar según E01–E08",
    } for record in records)
    context = CellContext("NCR18650PF", 2.994979138519, "límites declarados por ensayo")
    return ApplicationAuditResult("datos de referencia Panasonic", summaries, checks, tuple(), _generic_report(summaries, checks, context))


__all__ = [
    "ApplicationAuditResult", "CellContext", "ROLE_CHARACTERIZATION", "ROLE_DYNAMIC", "ROLE_LOW_RATE", "ROLE_THERMAL",
    "ROLE_UNKNOWN", "ROLE_VALIDATION", "SETUP_OTHER", "SETUP_REFERENCE", "SETUP_UNKNOWN", "STANDARD_CHROMA_COLUMNS",
    "UploadedFileSpec", "audit_criteria_payload", "audit_uploaded_files", "build_panasonic_reference_audit", "build_trc_case_study",
    "inspect_csv_bytes", "reference_model_payload", "rows_to_csv_bytes",
]
