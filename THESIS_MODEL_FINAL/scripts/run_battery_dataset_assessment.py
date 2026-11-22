"""Consolidate completed own-data artifacts into the final assessment pipeline.

The script reads existing audit, case-study, and transfer outputs. It does not
rerun a battery model, optimize parameters, alter voltage, or modify raw data.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from analysis.battery_dataset_assessment import (  # noqa: E402
    MetricRecord,
    MetricRole,
    TimebaseStatus,
    TransferEvidence,
    build_dataset_assessment,
    evaluate_transfer_gate,
    methodology_flow,
    validate_voltage_curve_label,
)


DATA_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
QUALITY_PATH = DATA_DIRECTORY / "own_data_quality_audit.csv"
CASE_METRICS_PATH = DATA_DIRECTORY / "own_case_study_cross_experiment_metrics.csv"
TRANSITIONS_PATH = DATA_DIRECTORY / "own_data_apparent_resistance_transitions.csv"
TRANSFER_METRICS_PATH = DATA_DIRECTORY / "own_model_transfer_metrics.csv"
TRANSFER_REPORT_PATH = DATA_DIRECTORY / "own_model_transfer_audit.md"
CASE_REPORT_PATH = DATA_DIRECTORY / "own_case_study_nonideality_report.md"

SUMMARY_PATH = DATA_DIRECTORY / "battery_dataset_assessment_summary.csv"
FINAL_REPORT_PATH = DATA_DIRECTORY / "uniandes_trc18650_case_study_final.md"
FLOW_MD_PATH = DATA_DIRECTORY / "battery_dataset_assessment_flow.md"
FLOW_PNG_PATH = DATA_DIRECTORY / "battery_dataset_assessment_flow.png"
REFERENCE_ID = "cycle_cell1_ambient"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required completed artifact is missing: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def transfer_value(rows: list[dict[str, str]], test: str, subset: str, metric: str) -> float:
    matches = [row for row in rows if row["test"] == test and row["subset"] == subset and row["metric"] == metric]
    if len(matches) != 1:
        raise ValueError(f"Expected one transfer metric for {(test, subset, metric)}, found {len(matches)}.")
    return float(matches[0]["value"])


def format_mv(value_v: float) -> str:
    return f"{1000.0 * value_v:.3f}"


def metric_record(
    name: str,
    value: float,
    unit: str,
    role: MetricRole,
    basis: str,
) -> MetricRecord:
    return MetricRecord(
        name=name,
        value=value,
        unit=unit,
        metric_role=role,
        evidence_basis=basis,
        timebase_status=TimebaseStatus.INFERRED_PROVISIONAL,
    )


def write_flow_markdown() -> None:
    flow = methodology_flow()
    lines = [
        "# Battery experimental-dataset assessment flow",
        "",
        "The fenced JSON block is the machine-readable source for the rendered flow. Branches explicitly prevent unsuitable datasets from being promoted automatically to calibration or quantitative validation.",
        "",
        "```json",
        json.dumps(flow, indent=2),
        "```",
        "",
        "```mermaid",
        "flowchart TD",
        "    A[NEW EXPERIMENTAL DATA] --> B[DATA QUALITY AUDIT]",
        "    B --> C[SUITABILITY CLASSIFICATION]",
        "    C -->|nonideal or insufficient| D[DIAGNOSTIC CASE STUDY]",
        "    C -->|evidence suitable and model supplied| E[MODEL TRANSFER TEST]",
        "    E --> F[PARAMETER-COMPATIBILITY AUDIT]",
        "    F -->|adaptation required but unsupported| G[REPORT LIMITATION]",
        "    F -->|all evidence supports metrics| H[QUANTITATIVE VALIDATION]",
        "    D --> G",
        "```",
        "",
        "Quantitative validation is conditional, not the default endpoint. The pipeline never triggers model fitting or voltage correction automatically.",
    ]
    FLOW_MD_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def draw_flow_figure() -> None:
    positions = {
        "new": (0.50, 0.93),
        "quality": (0.50, 0.80),
        "suitability": (0.50, 0.67),
        "diagnostic": (0.17, 0.50),
        "transfer": (0.50, 0.50),
        "compatibility": (0.50, 0.34),
        "limitation": (0.18, 0.16),
        "validation": (0.73, 0.16),
    }
    labels = {
        "new": "NEW EXPERIMENTAL DATA",
        "quality": "DATA QUALITY AUDIT",
        "suitability": "SUITABILITY CLASSIFICATION",
        "diagnostic": "DIAGNOSTIC CASE STUDY",
        "transfer": "MODEL TRANSFER TEST",
        "compatibility": "PARAMETER-COMPATIBILITY AUDIT",
        "limitation": "REPORT LIMITATION",
        "validation": "QUANTITATIVE VALIDATION\nonly when evidence supports it",
    }
    colors = {
        "new": "#D9EAF7",
        "quality": "#D9EAF7",
        "suitability": "#FFF2CC",
        "diagnostic": "#FCE4D6",
        "transfer": "#D9EAD3",
        "compatibility": "#FFF2CC",
        "limitation": "#F4CCCC",
        "validation": "#D9EAD3",
    }
    figure, axis = plt.subplots(figsize=(10.5, 8.3), constrained_layout=True)
    axis.set_xlim(0.0, 1.0)
    axis.set_ylim(0.05, 1.0)
    axis.axis("off")
    width, height = 0.31, 0.075
    for node, (x, y) in positions.items():
        box = FancyBboxPatch(
            (x - width / 2.0, y - height / 2.0),
            width,
            height,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            facecolor=colors[node],
            edgecolor="0.25",
            linewidth=1.2,
        )
        axis.add_patch(box)
        axis.text(x, y, labels[node], ha="center", va="center", fontsize=9.5, weight="bold")

    def arrow(start: str, end: str, label: str = "", curve: float = 0.0) -> None:
        x1, y1 = positions[start]
        x2, y2 = positions[end]
        patch = FancyArrowPatch(
            (x1, y1 - height / 2.0),
            (x2, y2 + height / 2.0),
            arrowstyle="-|>",
            mutation_scale=13,
            linewidth=1.2,
            color="0.25",
            connectionstyle=f"arc3,rad={curve}",
        )
        axis.add_patch(patch)
        if label:
            axis.text((x1 + x2) / 2.0, (y1 + y2) / 2.0 + 0.018, label, ha="center", va="center", fontsize=8, color="0.25")

    arrow("new", "quality")
    arrow("quality", "suitability")
    arrow("suitability", "transfer", "suitable evidence")
    arrow("suitability", "diagnostic", "nonideal / insufficient", 0.12)
    arrow("transfer", "compatibility")
    arrow("compatibility", "limitation", "adaptation unsupported", 0.12)
    arrow("compatibility", "validation", "evidence supports metrics", -0.10)
    arrow("diagnostic", "limitation", "retain and report")
    axis.text(0.50, 0.07, "No automatic fitting • No automatic voltage correction", ha="center", fontsize=10, weight="bold")
    axis.set_title("Evidence-gated battery experimental-dataset assessment", fontsize=15, pad=10)
    figure.savefig(FLOW_PNG_PATH, dpi=190)
    plt.close(figure)


def write_final_report(
    assessments,
    quality_rows,
    case_rows,
    transition_rows,
    transfer_rows,
    transfer_decision,
) -> None:
    quality_by_id = {row["dataset_id"]: row for row in quality_rows}
    case_by_id = {row["dataset"]: row for row in case_rows}
    assessment_by_id = {item.dataset: item for item in assessments}
    clean_quality = quality_by_id[REFERENCE_ID]
    clean_case = case_by_id[REFERENCE_ID]

    backbone_rmse = transfer_value(transfer_rows, "A_voltage_backbone", "overall", "rmse_v")
    backbone_bias = transfer_value(transfer_rows, "A_voltage_backbone", "overall", "bias_v")
    full_rmse = transfer_value(transfer_rows, "B_literal_full_transfer", "valid_ocp_support_only", "rmse_v")
    full_bias = transfer_value(transfer_rows, "B_literal_full_transfer", "valid_ocp_support_only", "bias_v")
    invalid_count = int(transfer_value(transfer_rows, "audit", "transfer_integrity", "discharge_surface_invalid_count"))
    full_valid = int(transfer_value(transfer_rows, "B_literal_full_transfer", "valid_ocp_support_only", "n_samples"))
    charge_rmse = transfer_value(transfer_rows, "directional_recharge", "all_recharge", "rmse_v")
    charge_cc_rmse = transfer_value(transfer_rows, "directional_recharge", "cc", "rmse_v")
    charge_cv_rmse = transfer_value(transfer_rows, "directional_recharge", "cv", "rmse_v")
    regional = {
        region: transfer_value(transfer_rows, "B_literal_full_transfer", f"soc_{region}_valid_support_only", "rmse_v")
        for region in ("high", "middle", "low")
    }

    reported_metrics = [
        metric_record("clean_reference_Q_discharge", float(clean_case["measured_discharge_capacity_ah"]), "Ah", MetricRole.DATA_QUALITY_DIAGNOSTIC, "Complete Chroma discharge counter."),
        metric_record("clean_reference_apparent_R", float(clean_case["median_r_app_ohm"]), "ohm", MetricRole.DATA_QUALITY_DIAGNOSTIC, "Adjacent rest-to-load DeltaV/DeltaI; not intrinsic R0."),
        metric_record("Panasonic_backbone_RMSE", backbone_rmse, "V", MetricRole.TRANSFER_DIAGNOSTIC, "Panasonic Ueq versus loaded TRC 1C curve; not OCP validation."),
        metric_record("literal_full_transfer_RMSE", full_rmse, "V", MetricRole.TRANSFER_DIAGNOSTIC, "Valid frozen-OCP-support samples using provisional time."),
        metric_record("charge_transfer_RMSE", charge_rmse, "V", MetricRole.TRANSFER_DIAGNOSTIC, "Held-out directional transfer using provisional time."),
    ]

    lines = [
        "# Universidad de los Andes TRC18650 final case-study assessment",
        "",
        "## Methodological outcome",
        "",
        "The case study is closed as an evidence-gated dataset assessment, not as a TRC parameter calibration. Data quality, model applicability, parameter transferability, and model accuracy eligibility are kept separate.",
        "",
        "**Direct Panasonic-to-TRC transfer is not adequate. Q_cell is directly supported by the clean own experiment. Capacity/state mapping requires cell-specific adaptation. A cell-compatible equilibrium-voltage/OCP representation is required before rigorous TRC reduced-model validation. Available own data do not identify that OCP or the dynamic/thermal parameters uniquely. Therefore, NO FULL TRC PARAMETER REFIT IS PERFORMED.**",
        "",
        "This is a methodological decision, not a failed optimization.",
        "",
        "## Metric-role register",
        "",
        "| Metric | Value | Unit | metric_role | Evidence basis |",
        "|---|---:|---|---|---|",
    ]
    for metric in reported_metrics:
        value = f"{metric.value:.9f}" if isinstance(metric.value, float) else str(metric.value)
        lines.append(f"| {metric.name} | {value} | {metric.unit} | {metric.metric_role.value} | {metric.evidence_basis} |")
    lines.extend([
        "",
        "No metric in this final TRC case study is labelled MODEL ACCURACY because the direct model-transfer gate failed. Voltage comparisons remain transfer or data-quality diagnostics.",
        "",
        "## A. Clean reference",
        "",
        f"`{REFERENCE_ID}` is the clean relative own-data reference. Its complete discharge capacity is `{float(clean_case['measured_discharge_capacity_ah']):.4f} Ah`; median apparent rest-to-load series resistance is `{float(clean_case['median_r_app_ohm']):.9f} ohm`. The latter is a **DATA-QUALITY DIAGNOSTIC**, not intrinsic battery R0.",
        "",
        f"Electrical calibration: **{clean_quality['electrical_calibration_suitability']}** — {clean_quality['electrical_calibration_evidence']}. Electrical validation: **{clean_quality['electrical_validation_suitability']}** — {clean_quality['electrical_validation_evidence']}.",
        "",
        "## B. Affected complete-cycle datasets",
        "",
        "Voltage deviations are affected-minus-clean at matched normalized reference SOC. They combine temperature, cell, state, electrochemical, and experimental/setup effects; no voltage correction is applied.",
        "",
        "| Dataset | Cell / nominal T | Capacity [Ah] | Median R_app [ohm] | Ratio | Mean / RMS voltage deviation [mV] | Nonideality flag | Interpretation | metric_role |",
        "|---|---|---:|---:|---:|---:|---|---|---|",
    ])
    for dataset in ("cycle_cell1_30c", "cycle_cell1_35c", "cycle_cell2_undeclared", "cycle_cell2_30c"):
        row = case_by_id[dataset]
        lines.append(
            f"| {dataset} | {row['cell']} / {row['nominal_temperature']} | {float(row['measured_discharge_capacity_ah']):.4f} | {float(row['median_r_app_ohm']):.9f} | {float(row['r_app_ratio_to_clean_reference']):.3f}x | {1000*float(row['mean_voltage_deviation_v']):.2f} / {1000*float(row['rms_voltage_deviation_v']):.2f} | {row['detector_flag']} | {row['interpretation']} | DATA-QUALITY DIAGNOSTIC |"
        )
    lines.extend([
        "",
        "These experiments provide evidence consistent with an elevated series-path response and experimental/setup nonideality. Fixture resistance was not uniquely identified. The datasets remain valuable because the methodology detects that they should not be treated as calibration ground truth.",
        "",
        "## C. Repeated pulse evidence",
        "",
        "| Dataset | Transitions | Median [ohm] | IQR [ohm] | Range [ohm] | Ratio | Finding | metric_role |",
        "|---|---:|---:|---:|---|---:|---|---|",
    ])
    for dataset in ("pulse_cell1_part1", "pulse_cell1_part2"):
        q = quality_by_id[dataset]
        values = [float(row["r_app_ohm"]) for row in transition_rows if row["dataset_id"] == dataset]
        lines.append(
            f"| {dataset} | {int(q['r_app_transition_count'])} | {float(q['r_app_median_ohm']):.9f} | {float(q['r_app_iqr_ohm']):.9f} | {min(values):.9f}-{max(values):.9f} | {float(q['r_app_ratio_vs_reference']):.3f}x | Persistent elevated group | DATA-QUALITY DIAGNOSTIC |"
        )
    lines.extend([
        "",
        "All ten retained pulse transitions remain in the elevated group. The elevation is persistent across both parts, not an isolated transition; its physical source remains non-unique.",
        "",
        "## D. Direct Panasonic model transfer",
        "",
        "| Diagnostic | Result | metric_role | Interpretation |",
        "|---|---:|---|---|",
        f"| Voltage-backbone RMSE | {format_mv(backbone_rmse)} mV | TRANSFER DIAGNOSTIC | Panasonic Ueq versus loaded TRC 1C voltage; not OCP validation |",
        f"| Voltage-backbone bias | {format_mv(backbone_bias)} mV | TRANSFER DIAGNOSTIC | Frozen backbone is systematically above loaded TRC voltage |",
        f"| Literal full-transfer RMSE | {format_mv(full_rmse)} mV | TRANSFER DIAGNOSTIC | Valid support only; provisional-time integration |",
        f"| Literal full-transfer bias | {format_mv(full_bias)} mV | TRANSFER DIAGNOSTIC | Predicted minus measured |",
        f"| High/middle/low-SOC RMSE | {format_mv(regional['high'])} / {format_mv(regional['middle'])} / {format_mv(regional['low'])} mV | TRANSFER DIAGNOSTIC | Deterioration is concentrated at low SOC |",
        f"| Surface-state support | {full_valid} valid; {invalid_count} invalid | TRANSFER DIAGNOSTIC | 22 low-SOC positive-surface violations; no clipping |",
        f"| Recharge RMSE | {format_mv(charge_rmse)} mV | TRANSFER DIAGNOSTIC | Held-out directional check |",
        f"| Recharge CC / CV RMSE | {format_mv(charge_cc_rmse)} / {format_mv(charge_cv_rmse)} mV | TRANSFER DIAGNOSTIC | No direction-specific correction |",
        "",
        "These results are not a calibrated TRC model and cannot be reported as model accuracy. The dynamic terms partially compensate the voltage-backbone mismatch, but capacity mapping remains incompatible and the discharge enters unsupported OCP surface states.",
        "",
        "## Four separate assessment dimensions",
        "",
        "| Dataset | Data quality | Model applicability | Parameter transferability | Model accuracy eligibility | Recommended use |",
        "|---|---|---|---|---|---|",
    ])
    for assessment in assessments:
        data_quality = f"E-cal {assessment.electrical_calibration.status.value}; E-val {assessment.electrical_validation.status.value}; T-cal {assessment.thermal_calibration.status.value}; T-val {assessment.thermal_validation.status.value}"
        lines.append(f"| {assessment.dataset} | {data_quality} | {assessment.model_applicability} | {assessment.parameter_transferability} | {assessment.model_accuracy_status} | {assessment.recommended_use} |")
    lines.extend([
        "",
        "## Final transfer decision",
        "",
        f"The hierarchical gate returns: `{transfer_decision.summary}`.",
        "",
        "- LEVEL 1 is required because the TRC 2.8999 Ah discharge does not traverse the frozen Panasonic capacity/state mapping completely.",
        "- LEVEL 2 is required because the voltage-backbone residual is large and SOC-dependent.",
        "- LEVEL 3 is indicated by residual dynamic mismatch but cannot be calibrated with provisional time.",
        "- LEVEL 4 records that full OCP, dynamic, and thermal identification is unsupported by the available own data.",
        "",
        "## Reporting policy",
        "",
        "- Voltage-vs-reference-SOC errors are quantitative only as **TRANSFER DIAGNOSTICS** until a compatible TRC parameterization exists.",
        "- Capacity and R_app are quantitative **DATA-QUALITY DIAGNOSTICS**; R_app is never intrinsic R0.",
        "- SOC accuracy is not reported because no independent SOC ground truth exists.",
        "- Dynamic time-domain results are provisional and not primary.",
        "- **THERMAL DYNAMIC ACCURACY NOT ASSESSABLE** for every own dataset because measured time is unavailable.",
        "- Temperature versus capacity/SOC remains **DESCRIPTIVE ONLY**.",
        "",
        "## Provenance and stop point",
        "",
        f"This consolidation reads `{QUALITY_PATH.name}`, `{CASE_METRICS_PATH.name}`, `{TRANSITIONS_PATH.name}`, and `{TRANSFER_METRICS_PATH.name}`. Detailed source analyses remain in `{CASE_REPORT_PATH.name}` and `{TRANSFER_REPORT_PATH.name}`.",
        "",
        "No model was rerun, no parameter was fitted, no OCP/window was rescaled, no voltage was corrected, and affected data were not used as calibration truth.",
    ])
    FINAL_REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    for path in (TRANSFER_REPORT_PATH, CASE_REPORT_PATH):
        if not path.is_file():
            raise FileNotFoundError(f"Required completed report is missing: {path}")
    quality_rows = read_csv(QUALITY_PATH)
    case_rows = read_csv(CASE_METRICS_PATH)
    transition_rows = read_csv(TRANSITIONS_PATH)
    transfer_rows = read_csv(TRANSFER_METRICS_PATH)
    quality_by_id = {row["dataset_id"]: row for row in quality_rows}
    case_by_id = {row["dataset"]: row for row in case_rows}
    if set(quality_by_id) != set(case_by_id):
        raise RuntimeError("Audit and case-study dataset identities differ.")

    validate_voltage_curve_label("loaded_1c", "clean ambient 1C discharge reference curve")
    transfer_decision = evaluate_transfer_gate(
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
    if transfer_decision.automatic_fit_triggered:
        raise RuntimeError("Assessment methodology must never trigger fitting automatically.")

    assessments = [
        build_dataset_assessment(
            quality_by_id[dataset],
            case_by_id[dataset],
            clean_reference_id=REFERENCE_ID,
            transfer_decision=transfer_decision,
        )
        for dataset in quality_by_id
    ]
    write_csv(SUMMARY_PATH, [assessment.to_row() for assessment in assessments])
    write_final_report(
        assessments,
        quality_rows,
        case_rows,
        transition_rows,
        transfer_rows,
        transfer_decision,
    )
    write_flow_markdown()
    draw_flow_figure()
    print(f"Consolidated {len(assessments)} own datasets without rerunning or fitting a model.")
    print(f"Transfer gate: {transfer_decision.summary}")
    print(f"Outputs written under {DATA_DIRECTORY}")


if __name__ == "__main__":
    main()
