"""High-level, evidence-gated battery experimental-dataset assessment.

This module consolidates data quality, experimental nonideality, model
transferability, and model-accuracy eligibility.  It deliberately contains no
optimizer, voltage-correction routine, or automatic parameter-fitting path.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Mapping

import numpy as np

from analysis.experimental_data_audit import (
    StandardExperiment,
    audit_experiment,
    extract_rest_to_load_transitions,
    largest_log_gap_reference_group,
)


class TimebaseStatus(str, Enum):
    MEASURED = "measured"
    INFERRED_PROVISIONAL = "inferred/provisional"
    UNAVAILABLE = "unavailable"


class SuitabilityStatus(str, Enum):
    SUITABLE = "SUITABLE"
    SUITABLE_WITH_CAVEAT = "SUITABLE WITH CAVEAT"
    DIAGNOSTIC_ONLY = "QUALITATIVE / DIAGNOSTIC ONLY"
    NOT_SUITABLE = "NOT SUITABLE"


class NonidealityFlag(str, Enum):
    NO_RELATIVE_ELEVATION = "NO RELATIVE ELEVATION DETECTED"
    ELEVATED_SERIES_PATH = "EVIDENCE OF ELEVATED SERIES-PATH RESPONSE"
    INSUFFICIENT_EVIDENCE = "INSUFFICIENT EVIDENCE"


class MetricRole(str, Enum):
    MODEL_ACCURACY = "MODEL ACCURACY"
    TRANSFER_DIAGNOSTIC = "TRANSFER DIAGNOSTIC"
    DATA_QUALITY_DIAGNOSTIC = "DATA-QUALITY DIAGNOSTIC"
    DESCRIPTIVE_ONLY = "DESCRIPTIVE ONLY"


class TransferLevel(str, Enum):
    LEVEL_0 = "LEVEL 0: DIRECT PARAMETER TRANSFER SUFFICIENT"
    LEVEL_1 = "LEVEL 1: CAPACITY / STATE-MAPPING ADAPTATION REQUIRED"
    LEVEL_2 = "LEVEL 2: CELL-SPECIFIC VOLTAGE-REFERENCE / OCP ADAPTATION REQUIRED"
    LEVEL_3 = "LEVEL 3: DYNAMIC PARAMETER ADAPTATION ALSO REQUIRED"
    LEVEL_4 = "LEVEL 4: AVAILABLE DATA INSUFFICIENT FOR DEFENSIBLE QUANTITATIVE MODEL IDENTIFICATION"


@dataclass(frozen=True)
class EvidenceDecision:
    status: SuitabilityStatus
    evidence: str

    def __post_init__(self) -> None:
        if not self.evidence.strip():
            raise ValueError("Every suitability decision requires explicit evidence.")


@dataclass(frozen=True)
class MetricRecord:
    name: str
    value: float | int | str
    unit: str
    metric_role: MetricRole
    evidence_basis: str
    domain: str = "electrical"
    timebase_status: TimebaseStatus = TimebaseStatus.UNAVAILABLE
    independent_soc_reference: bool = False

    def __post_init__(self) -> None:
        if not self.evidence_basis.strip():
            raise ValueError("Every metric requires an evidence basis.")
        normalized = self.name.lower().replace("_", " ")
        if "r app" in normalized and ("intrinsic" in normalized or "r0" in normalized):
            raise ValueError("R_app cannot be labelled intrinsic battery R0.")
        if self.domain == "thermal_dynamic" and self.timebase_status != TimebaseStatus.MEASURED:
            raise ValueError("THERMAL DYNAMIC ACCURACY NOT ASSESSABLE without measured time.")
        if self.domain == "soc_accuracy" and not self.independent_soc_reference:
            raise ValueError("SOC accuracy requires an independent SOC reference.")


@dataclass(frozen=True)
class OcpDomainScreen:
    values: np.ndarray
    valid_mask: np.ndarray
    lower: float
    upper: float

    @property
    def invalid_count(self) -> int:
        return int(np.count_nonzero(~self.valid_mask))


@dataclass(frozen=True)
class TransferEvidence:
    model_supplied: bool
    direct_transfer_adequate: bool | None
    capacity_mapping_compatible: bool | None
    voltage_reference_compatible: bool | None
    dynamic_parameters_compatible: bool | None
    capacity_information_available: bool
    voltage_reference_information_available: bool
    measured_dynamic_time_available: bool


@dataclass(frozen=True)
class TransferDecision:
    levels: tuple[TransferLevel, ...]
    stopped_because: str
    automatic_fit_triggered: bool
    next_action: str

    @property
    def summary(self) -> str:
        return "; ".join(level.value for level in self.levels)


@dataclass(frozen=True)
class DatasetAssessment:
    dataset: str
    cell: str
    nominal_temperature: str
    test_type: str
    timebase_status: TimebaseStatus
    capacity_supported: bool
    apparent_r_median_ohm: float
    apparent_r_ratio: float
    nonideality_flag: NonidealityFlag
    electrical_calibration: EvidenceDecision
    electrical_validation: EvidenceDecision
    thermal_calibration: EvidenceDecision
    thermal_validation: EvidenceDecision
    model_applicability: str
    parameter_transferability: str
    model_accuracy_status: str
    quantitative_voltage_metrics_supported: str
    quantitative_thermal_metrics_supported: str
    recommended_use: str

    def to_row(self) -> dict[str, object]:
        return {
            "dataset": self.dataset,
            "cell": self.cell,
            "nominal_temperature": self.nominal_temperature,
            "test_type": self.test_type,
            "measured_time_available": self.timebase_status == TimebaseStatus.MEASURED,
            "timebase_status": self.timebase_status.value,
            "capacity_supported": self.capacity_supported,
            "apparent_R_median": self.apparent_r_median_ohm,
            "apparent_R_ratio": self.apparent_r_ratio,
            "apparent_R_metric_role": MetricRole.DATA_QUALITY_DIAGNOSTIC.value,
            "nonideality_flag": self.nonideality_flag.value,
            "electrical_calibration_status": self.electrical_calibration.status.value,
            "electrical_calibration_evidence": self.electrical_calibration.evidence,
            "electrical_validation_status": self.electrical_validation.status.value,
            "electrical_validation_evidence": self.electrical_validation.evidence,
            "thermal_calibration_status": self.thermal_calibration.status.value,
            "thermal_calibration_evidence": self.thermal_calibration.evidence,
            "thermal_validation_status": self.thermal_validation.status.value,
            "thermal_validation_evidence": self.thermal_validation.evidence,
            "model_applicability": self.model_applicability,
            "transfer_level_required": self.parameter_transferability,
            "model_accuracy_status": self.model_accuracy_status,
            "quantitative_voltage_metrics_supported": self.quantitative_voltage_metrics_supported,
            "voltage_metric_role": MetricRole.TRANSFER_DIAGNOSTIC.value if self.dataset == "cycle_cell1_ambient" else MetricRole.DATA_QUALITY_DIAGNOSTIC.value,
            "quantitative_thermal_metrics_supported": self.quantitative_thermal_metrics_supported,
            "thermal_metric_role": MetricRole.DESCRIPTIVE_ONLY.value,
            "recommended_use": self.recommended_use,
            "automatic_voltage_correction": False,
            "automatic_model_fitting": False,
        }


@dataclass(frozen=True)
class DataQualityGateResult:
    audit_metrics: Mapping[str, object]
    timebase_status: TimebaseStatus
    electrical_calibration: EvidenceDecision
    electrical_validation: EvidenceDecision
    thermal_calibration: EvidenceDecision
    thermal_validation: EvidenceDecision


@dataclass(frozen=True)
class RelativeNonidealityResult:
    dataset: str
    flag: NonidealityFlag
    transition_count: int
    median_r_app_ohm: float
    iqr_r_app_ohm: float
    reference_dataset: str
    reference_median_r_app_ohm: float
    reference_iqr_r_app_ohm: float
    ratio_to_reference: float
    reference_group_provenance: str
    evidence: str


def parse_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    if str(value).strip().lower() in {"true", "1", "yes"}:
        return True
    if str(value).strip().lower() in {"false", "0", "no", ""}:
        return False
    raise ValueError(f"Cannot parse boolean value {value!r}.")


def evaluate_data_quality_gate(experiment: StandardExperiment) -> DataQualityGateResult:
    """Run the reusable low-level audit and issue conservative preliminary decisions.

    These Stage-A decisions precede the cross-experiment nonideality gate and
    may be downgraded later. They never authorize parameter fitting.
    """

    transitions = extract_rest_to_load_transitions(experiment)
    audit = audit_experiment(experiment, transitions)
    timebase = resolve_timebase(
        measured_available=parse_bool(audit["raw_time_available"]),
        provisional_available="provisional" in str(audit["time_basis"]),
    )
    capacity_complete = parse_bool(audit["capacity_sequence_appears_complete"])
    temperature_available = bool(str(audit["temperature_channels"]).strip())
    ambient_available = any(
        (key.endswith("_role") and ("ambient" in str(value) or "chamber" in str(value)))
        for key, value in audit.items()
    )
    transition_count = int(audit["r_app_transition_count"])

    if timebase == TimebaseStatus.MEASURED and transition_count >= 2:
        electrical_calibration = EvidenceDecision(
            SuitabilityStatus.SUITABLE_WITH_CAVEAT,
            "Measured time and repeatable transitions exist; cross-experiment nonideality and excitation sufficiency remain to be checked.",
        )
    else:
        electrical_calibration = EvidenceDecision(
            SuitabilityStatus.NOT_SUITABLE,
            "Measured time or sufficient transition excitation is absent; calibration is not authorized by Stage A.",
        )
    electrical_validation = EvidenceDecision(
        SuitabilityStatus.SUITABLE_WITH_CAVEAT if capacity_complete else SuitabilityStatus.DIAGNOSTIC_ONLY,
        "Complete capacity coverage supports capacity-coordinate comparison, pending nonideality screening."
        if capacity_complete
        else "Incomplete capacity coverage limits use to qualitative/diagnostic comparison.",
    )
    if timebase == TimebaseStatus.MEASURED and temperature_available and ambient_available:
        thermal_calibration = EvidenceDecision(
            SuitabilityStatus.SUITABLE_WITH_CAVEAT,
            "Measured time plus surface and ambient channels exist; heat-input identifiability still requires review.",
        )
        thermal_validation = EvidenceDecision(
            SuitabilityStatus.SUITABLE_WITH_CAVEAT,
            "Measured time and thermal boundary channels exist, pending sensor-integrity review.",
        )
    else:
        thermal_calibration = EvidenceDecision(
            SuitabilityStatus.NOT_SUITABLE,
            "Thermal calibration requires measured time and defensible surface/ambient boundary data.",
        )
        thermal_validation = EvidenceDecision(
            SuitabilityStatus.DIAGNOSTIC_ONLY,
            "Temperature channels may be described, but dynamic thermal validation is unsupported without measured time and boundary evidence.",
        )
    return DataQualityGateResult(
        audit_metrics=audit,
        timebase_status=timebase,
        electrical_calibration=electrical_calibration,
        electrical_validation=electrical_validation,
        thermal_calibration=thermal_calibration,
        thermal_validation=thermal_validation,
    )


def evaluate_relative_nonideality_gate(
    audit_rows: list[Mapping[str, object]],
    reference_dataset: str,
    *,
    transition_count_field: str = "r_app_transition_count",
    median_field: str = "r_app_median_ohm",
    iqr_field: str = "r_app_iqr_ohm",
    minimum_transition_count: int = 2,
    metric_label: str = "R_app",
) -> list[RelativeNonidealityResult]:
    """Use the audit's data-driven largest-log-gap grouping; no fixed threshold.

    Field selection is explicit so a directional audit cannot silently mix
    charge and discharge transitions.  Existing callers retain the original
    all-transition fields and two-transition evidence requirement.
    """

    if minimum_transition_count < 1:
        raise ValueError("minimum_transition_count must be at least one.")

    by_id = {str(row["dataset_id"]): row for row in audit_rows}
    if reference_dataset not in by_id:
        raise ValueError(f"Reference dataset {reference_dataset!r} is absent.")
    medians = {
        dataset: float(row[median_field])
        for dataset, row in by_id.items()
        if np.isfinite(float(row[median_field]))
    }
    reference_group, _ = largest_log_gap_reference_group(medians)
    if reference_dataset not in reference_group:
        raise ValueError("Requested reference is not in the data-driven lowest group.")
    reference = by_id[reference_dataset]
    reference_median = float(reference[median_field])
    reference_iqr = float(reference[iqr_field])
    provenance = (
        f"largest adjacent gap in log experiment-median {metric_label}; "
        "reference selected from lowest own-data group"
    )
    results: list[RelativeNonidealityResult] = []
    for dataset, row in by_id.items():
        count = int(row[transition_count_field])
        median = float(row[median_field])
        iqr = float(row[iqr_field])
        ratio = median / reference_median if np.isfinite(median) else float("nan")
        if count < minimum_transition_count or not np.isfinite(median):
            flag = NonidealityFlag.INSUFFICIENT_EVIDENCE
            evidence = (
                f"Fewer than {minimum_transition_count} comparable transitions "
                "or no finite experiment median."
            )
        elif dataset in reference_group:
            flag = NonidealityFlag.NO_RELATIVE_ELEVATION
            evidence = "Experiment belongs to the data-driven lowest relative-response group; this is not intrinsic-R0 certification."
        else:
            flag = NonidealityFlag.ELEVATED_SERIES_PATH
            evidence = "Experiment belongs to the data-driven elevated group; evidence is consistent with an elevated series-path response, with no unique physical attribution."
        results.append(
            RelativeNonidealityResult(
                dataset=dataset,
                flag=flag,
                transition_count=count,
                median_r_app_ohm=median,
                iqr_r_app_ohm=iqr,
                reference_dataset=reference_dataset,
                reference_median_r_app_ohm=reference_median,
                reference_iqr_r_app_ohm=reference_iqr,
                ratio_to_reference=ratio,
                reference_group_provenance=provenance,
                evidence=evidence,
            )
        )
    return results


def resolve_timebase(*, measured_available: bool, provisional_available: bool) -> TimebaseStatus:
    if measured_available:
        return TimebaseStatus.MEASURED
    if provisional_available:
        return TimebaseStatus.INFERRED_PROVISIONAL
    return TimebaseStatus.UNAVAILABLE


def convert_current_for_model(raw_current_a: np.ndarray, raw_current_sign: str) -> np.ndarray:
    """Return a discharge-positive copy without mutating or aliasing raw data."""

    raw = np.asarray(raw_current_a, dtype=float)
    if raw_current_sign == "discharge_negative":
        converted = -raw
    elif raw_current_sign == "discharge_positive":
        converted = raw.copy()
    else:
        raise ValueError("Ambiguous current sign cannot be converted automatically.")
    return np.array(converted, copy=True)


def validate_voltage_curve_label(source_kind: str, requested_label: str) -> str:
    if source_kind.lower() in {"loaded_1c", "loaded terminal voltage", "loaded"} and "ocv" in requested_label.lower():
        raise ValueError("Loaded 1C voltage cannot be labelled OCV.")
    return requested_label


def screen_ocp_domain(values: np.ndarray, lower: float, upper: float) -> OcpDomainScreen:
    """Return explicit support flags and preserve unsupported values unchanged."""

    array = np.array(values, dtype=float, copy=True)
    valid = np.isfinite(array) & (array >= lower) & (array <= upper)
    return OcpDomainScreen(values=array, valid_mask=valid, lower=lower, upper=upper)


def request_automatic_voltage_correction() -> None:
    raise RuntimeError("Automatic voltage correction is prohibited by the assessment methodology.")


def require_thermal_dynamic_support(timebase_status: TimebaseStatus, heat_boundary_supported: bool) -> str:
    if timebase_status != TimebaseStatus.MEASURED or not heat_boundary_supported:
        return "THERMAL DYNAMIC ACCURACY NOT ASSESSABLE"
    return "THERMAL DYNAMIC ACCURACY ASSESSABLE"


def translate_nonideality_flag(value: str) -> NonidealityFlag:
    normalized = value.strip().upper()
    if normalized == NonidealityFlag.NO_RELATIVE_ELEVATION.value:
        return NonidealityFlag.NO_RELATIVE_ELEVATION
    if normalized == NonidealityFlag.ELEVATED_SERIES_PATH.value:
        return NonidealityFlag.ELEVATED_SERIES_PATH
    if normalized in {"INSUFFICIENT TRANSITIONS", "INSUFFICIENT EVIDENCE"}:
        return NonidealityFlag.INSUFFICIENT_EVIDENCE
    raise ValueError(f"Unknown nonideality flag: {value}")


def evaluate_transfer_gate(evidence: TransferEvidence) -> TransferDecision:
    """Evaluate transfer hierarchically and never trigger optimization."""

    if not evidence.model_supplied:
        return TransferDecision(
            levels=(TransferLevel.LEVEL_4,),
            stopped_because="No pre-existing calibrated model was supplied.",
            automatic_fit_triggered=False,
            next_action="Report model-transfer gate as not entered; do not fit automatically.",
        )
    if evidence.direct_transfer_adequate is True:
        return TransferDecision(
            levels=(TransferLevel.LEVEL_0,),
            stopped_because="Direct transfer met the evidence requirements.",
            automatic_fit_triggered=False,
            next_action="Proceed only to evidence-supported validation metrics.",
        )

    levels: list[TransferLevel] = []
    if not evidence.capacity_information_available or evidence.capacity_mapping_compatible is None:
        return TransferDecision(
            levels=(TransferLevel.LEVEL_4,),
            stopped_because="Capacity/state-mapping evidence is unavailable.",
            automatic_fit_triggered=False,
            next_action="Report limitation; do not infer or fit the missing mapping.",
        )
    if evidence.capacity_mapping_compatible is False:
        levels.append(TransferLevel.LEVEL_1)

    if not evidence.voltage_reference_information_available or evidence.voltage_reference_compatible is None:
        levels.append(TransferLevel.LEVEL_4)
        return TransferDecision(
            levels=tuple(levels),
            stopped_because="Cell-compatible voltage-reference evidence is unavailable.",
            automatic_fit_triggered=False,
            next_action="Report limitation; do not construct an OCP from loaded voltage.",
        )
    if evidence.voltage_reference_compatible is False:
        levels.append(TransferLevel.LEVEL_2)

    if evidence.dynamic_parameters_compatible is False:
        levels.append(TransferLevel.LEVEL_3)
    if not evidence.measured_dynamic_time_available:
        levels.append(TransferLevel.LEVEL_4)
        return TransferDecision(
            levels=tuple(dict.fromkeys(levels)),
            stopped_because="Measured dynamic time is unavailable; dynamic identification is not defensible.",
            automatic_fit_triggered=False,
            next_action="Report provisional transfer diagnostics only; do not fit dynamics.",
        )
    return TransferDecision(
        levels=tuple(dict.fromkeys(levels)) or (TransferLevel.LEVEL_0,),
        stopped_because="All supplied compatibility evidence was evaluated.",
        automatic_fit_triggered=False,
        next_action="Proceed only with separately authorized, evidence-supported work.",
    )


def decision_from_row(row: Mapping[str, object], prefix: str) -> EvidenceDecision:
    status = SuitabilityStatus(str(row[f"{prefix}_suitability"]))
    evidence = str(row[f"{prefix}_evidence"])
    return EvidenceDecision(status=status, evidence=evidence)


def build_dataset_assessment(
    audit_row: Mapping[str, object],
    case_row: Mapping[str, object],
    *,
    clean_reference_id: str,
    transfer_decision: TransferDecision,
) -> DatasetAssessment:
    dataset = str(audit_row["dataset_id"])
    timebase = resolve_timebase(
        measured_available=parse_bool(audit_row["raw_time_available"]),
        provisional_available="provisional" in str(audit_row["time_basis"]),
    )
    nonideality = translate_nonideality_flag(str(case_row["detector_flag"]))
    is_reference = dataset == clean_reference_id
    if is_reference:
        model_applicability = "CONDITIONAL ELECTRICAL TRANSFER DIAGNOSTICS IN CAPACITY/SOC COORDINATES"
        transferability = transfer_decision.summary
        accuracy_status = "NO TRC MODEL ACCURACY CLAIM: INCOMPATIBLE DIRECT TRANSFER"
        quantitative_voltage = "CONDITIONAL: transfer-diagnostic voltage metrics only"
    else:
        model_applicability = "DIAGNOSTIC CASE STUDY ONLY AFTER NONIDEALITY GATE"
        transferability = "NOT ENTERED: experimental nonideality/data-quality gate"
        accuracy_status = "NOT ASSESSED"
        quantitative_voltage = "NO: relative voltage/R_app diagnostics only"
    return DatasetAssessment(
        dataset=dataset,
        cell=str(audit_row["cell_id"]),
        nominal_temperature=str(audit_row["nominal_test_temperature"]),
        test_type=str(audit_row["test_type_signal_verified"]),
        timebase_status=timebase,
        capacity_supported=parse_bool(audit_row["capacity_sequence_appears_complete"]),
        apparent_r_median_ohm=float(audit_row["r_app_median_ohm"]),
        apparent_r_ratio=float(audit_row["r_app_ratio_vs_reference"]),
        nonideality_flag=nonideality,
        electrical_calibration=decision_from_row(audit_row, "electrical_calibration"),
        electrical_validation=decision_from_row(audit_row, "electrical_validation"),
        thermal_calibration=decision_from_row(audit_row, "thermal_calibration"),
        thermal_validation=decision_from_row(audit_row, "thermal_validation"),
        model_applicability=model_applicability,
        parameter_transferability=transferability,
        model_accuracy_status=accuracy_status,
        quantitative_voltage_metrics_supported=quantitative_voltage,
        quantitative_thermal_metrics_supported=require_thermal_dynamic_support(timebase, False),
        recommended_use=str(audit_row["recommended_thesis_use"]),
    )


def methodology_flow() -> dict[str, object]:
    nodes = [
        {"id": "new_data", "label": "NEW EXPERIMENTAL DATA"},
        {"id": "quality", "label": "DATA QUALITY AUDIT"},
        {"id": "suitability", "label": "SUITABILITY CLASSIFICATION"},
        {"id": "diagnostic", "label": "DIAGNOSTIC CASE STUDY"},
        {"id": "transfer", "label": "MODEL TRANSFER TEST"},
        {"id": "compatibility", "label": "PARAMETER-COMPATIBILITY AUDIT"},
        {"id": "limitation", "label": "REPORT LIMITATION"},
        {"id": "validation", "label": "QUANTITATIVE VALIDATION"},
    ]
    edges = [
        {"from": "new_data", "to": "quality", "condition": "always"},
        {"from": "quality", "to": "suitability", "condition": "always"},
        {"from": "suitability", "to": "diagnostic", "condition": "nonideal or insufficient"},
        {"from": "suitability", "to": "transfer", "condition": "suitable evidence and pre-existing model"},
        {"from": "transfer", "to": "compatibility", "condition": "always after transfer test"},
        {"from": "compatibility", "to": "limitation", "condition": "adaptation required but unsupported"},
        {"from": "compatibility", "to": "validation", "condition": "all required experimental evidence supports metrics"},
    ]
    return {"schema_version": 1, "nodes": nodes, "edges": edges, "automatic_fit": False, "automatic_voltage_correction": False}


__all__ = [
    "DatasetAssessment",
    "DataQualityGateResult",
    "EvidenceDecision",
    "MetricRecord",
    "MetricRole",
    "NonidealityFlag",
    "OcpDomainScreen",
    "RelativeNonidealityResult",
    "SuitabilityStatus",
    "TimebaseStatus",
    "TransferDecision",
    "TransferEvidence",
    "TransferLevel",
    "build_dataset_assessment",
    "convert_current_for_model",
    "evaluate_data_quality_gate",
    "evaluate_relative_nonideality_gate",
    "evaluate_transfer_gate",
    "methodology_flow",
    "request_automatic_voltage_correction",
    "require_thermal_dynamic_support",
    "resolve_timebase",
    "screen_ocp_domain",
    "validate_voltage_curve_label",
]
