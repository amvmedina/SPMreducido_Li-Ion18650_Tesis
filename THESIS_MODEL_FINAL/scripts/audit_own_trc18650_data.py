"""Reproducible audit of all own TRC18650/Chroma battery experiments.

No electrochemical or thermal model is run. Raw files are read-only and any
timebase reconstructed from capacity/current is explicitly provisional.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from analysis.experimental_data_audit import (  # noqa: E402
    ExperimentMetadata,
    StandardExperiment,
    audit_experiment,
    extract_rest_to_load_transitions,
    largest_log_gap_reference_group,
    load_chroma_experiment,
)


SOURCE_DIRECTORY = Path(r"C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data")
CONTROL_DIRECTORY = Path(r"C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos")
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
FIGURE_DIRECTORY = OUTPUT_DIRECTORY / "own_data_audit_figures"
MANIFEST_PATH = OUTPUT_DIRECTORY / "own_data_manifest.csv"
QUALITY_PATH = OUTPUT_DIRECTORY / "own_data_quality_audit.csv"
TRANSITION_PATH = OUTPUT_DIRECTORY / "own_data_apparent_resistance_transitions.csv"
REPORT_PATH = OUTPUT_DIRECTORY / "own_data_audit.md"
CROSS_RESISTANCE_PATH = FIGURE_DIRECTORY / "apparent_resistance_cross_experiment.png"

PART2_COLUMNS = (
    "Description",
    "Step No.",
    "Voltage(V)",
    "Current(A)",
    "Capacity(Ah)",
    "Mode",
    "Total Capacity(Ah)",
    "Temp1",
    "Temp2",
    "Temp7",
)


def experiment_metadata() -> list[ExperimentMetadata]:
    """File identity metadata only; test type and conclusions come from signals."""

    cell1_roles = {
        "Temp1": "cell surface middle",
        "Temp2": "cell surface end",
        "Temp7": "chamber/ambient",
    }
    cell2_roles = {"Temp3": "temperature location unknown", "Temp4": "temperature location unknown"}
    definitions = [
        ("pulse_cell1_part1", "amvTRC18650LiDisPulse1C cell1 part1.csv", "cell 1", "ambient/undeclared", "not encoded", cell1_roles),
        ("pulse_cell1_part2", "amvTRC18650LiDisPulse1C cell1 part2.csv", "cell 1", "ambient/undeclared", "not encoded", cell1_roles),
        ("cycle_cell1_ambient", "amvTRC18650LiIonCycle1C Tambient cell1.csv", "cell 1", "Aire 22–24 °C", "not encoded", cell1_roles),
        ("cycle_cell1_30c", "amvTRC18650LiIonCycle1C T30C cell1.csv", "cell 1", "Aire 29–31 °C", "not encoded", cell1_roles),
        ("cycle_cell1_35c", "amvTRC18650LiIonCycle1C T35C cell1.csv", "cell 1", "Aire 33–37 °C", "not encoded", cell1_roles),
        ("cycle_cell2_undeclared", "amvTRC18650LiIonCycle1C cell2 07232026.csv", "cell 2", "Aire 22–24 °C", "2026-07-23", cell2_roles),
        ("cycle_cell2_30c", "amvTRC18650LiIonCycle1C T30C cell2 07242026.csv", "cell 2", "Aire 29–31 °C", "2026-07-24", cell2_roles),
    ]
    return [
        ExperimentMetadata(
            dataset_id=dataset_id,
            path=SOURCE_DIRECTORY / filename,
            cell_id=cell_id,
            nominal_temperature=temperature,
            encoded_date=date,
            temperature_roles=dict(roles),
        )
        for dataset_id, filename, cell_id, temperature, date, roles in definitions
    ]


def load_all() -> list[StandardExperiment]:
    expected = {metadata.path.name for metadata in experiment_metadata()}
    discovered = {path.name for path in SOURCE_DIRECTORY.glob("*.csv")}
    if expected != discovered:
        raise RuntimeError(
            f"Own-data discovery changed. Missing={sorted(expected-discovered)}, unexpected={sorted(discovered-expected)}"
        )
    experiments: list[StandardExperiment] = []
    for metadata in experiment_metadata():
        fallback = PART2_COLUMNS if metadata.dataset_id == "pulse_cell1_part2" else None
        experiments.append(load_chroma_experiment(metadata, fallback))
    return experiments


def classify_and_augment(
    audits: list[dict[str, object]],
    transitions: list[dict[str, object]],
) -> tuple[set[str], float, float, float]:
    medians = {
        str(row["dataset_id"]): float(row["r_app_median_ohm"])
        for row in audits
        if np.isfinite(float(row["r_app_median_ohm"]))
    }
    reference_group, split_threshold = largest_log_gap_reference_group(medians)
    reference_values = np.asarray(
        [
            float(row["r_app_ohm"])
            for row in transitions
            if str(row["dataset_id"]) in reference_group
        ]
    )
    reference_median = float(np.median(reference_values))
    reference_iqr = float(np.quantile(reference_values, 0.75) - np.quantile(reference_values, 0.25))
    for row in audits:
        dataset_id = str(row["dataset_id"])
        count = int(row["r_app_transition_count"])
        median = float(row["r_app_median_ohm"])
        has_ambient = any(
            "ambient" in str(value) or "chamber" in str(value)
            for key, value in row.items()
            if key.endswith("_role")
        )
        temp_roles_known = has_ambient or any(
            "surface" in str(value) for key, value in row.items() if key.endswith("_role")
        )
        if count == 0 or not np.isfinite(median):
            resistance_class = "uncertain: no reliable comparable transition"
            excess = float("nan")
            ratio = float("nan")
        elif dataset_id in reference_group:
            resistance_class = "lowest internally consistent reference group"
            excess = median - reference_median
            ratio = median / reference_median
        else:
            resistance_class = "substantially elevated relative to reference group"
            excess = median - reference_median
            ratio = median / reference_median
        values = [
            float(item["r_app_ohm"])
            for item in transitions
            if str(item["dataset_id"]) == dataset_id
        ]
        if count <= 2:
            variability = "insufficient transitions for distribution stability"
        elif min(values) > split_threshold:
            variability = "internally consistent elevated band"
        elif max(values) < split_threshold:
            variability = "internally consistent low band"
        else:
            variability = "highly variable / straddles data-driven group split"
        row.update(
            {
                "resistance_group": resistance_class,
                "r_app_variability": variability,
                "r_app_excess_vs_reference_ohm": excess,
                "r_app_ratio_vs_reference": ratio,
                "reference_group_median_ohm": reference_median,
                "reference_group_iqr_ohm": reference_iqr,
                "data_driven_group_split_ohm": split_threshold,
                "evidence_of_parasitic_setup_effects": (
                    "no relative elevation detected; still not an intrinsic-R0 certification"
                    if dataset_id in reference_group
                    else "evidence consistent with parasitic series resistance; source not uniquely identifiable"
                ),
            }
        )

        if dataset_id in reference_group:
            row["electrical_calibration_suitability"] = "NOT SUITABLE"
            row["electrical_calibration_evidence"] = (
                "lowest relative apparent-resistance group and complete cycle, but no measured time channel and insufficient dynamic excitation"
            )
            row["electrical_validation_suitability"] = "SUITABLE WITH CAVEAT"
            row["electrical_validation_evidence"] = (
                "cleanest current/voltage transitions and complete capacity coverage; timebase remains provisional"
            )
            row["thermal_calibration_suitability"] = "NOT SUITABLE"
            row["thermal_calibration_evidence"] = (
                "surface/ambient temperatures exist, but thermal dynamics cannot be calibrated defensibly without measured timestamps"
            )
            row["thermal_validation_suitability"] = "QUALITATIVE / DIAGNOSTIC ONLY"
            row["thermal_validation_evidence"] = (
                "stable surface/ambient channels, but duration and gaps are reconstructed rather than measured"
            )
            row["recommended_thesis_use"] = (
                "strongest own quantitative voltage/capacity comparison with timing caveat; qualitative thermal comparison"
            )
        else:
            row["electrical_calibration_suitability"] = "NOT SUITABLE"
            row["electrical_calibration_evidence"] = (
                f"apparent rest-to-load resistance is {ratio:.2f}x the clean reference median; no measured time channel"
            )
            row["electrical_validation_suitability"] = "QUALITATIVE / DIAGNOSTIC ONLY"
            row["electrical_validation_evidence"] = (
                "voltage response contains a systematic elevated transition drop consistent with setup/path nonideality"
            )
            row["thermal_calibration_suitability"] = "NOT SUITABLE"
            row["thermal_calibration_evidence"] = (
                "unknown partition of extra electrical loss inside/outside cell plus absent measured timestamps"
            )
            row["thermal_validation_suitability"] = "QUALITATIVE / DIAGNOSTIC ONLY"
            if dataset_id == "cycle_cell2_undeclared":
                thermal_context = (
                    " and Temp3/Temp4 placement is undocumented; the shared measured "
                    "Temp7 air channel is confirmed but cross-file sample alignment is not"
                )
            elif not temp_roles_known:
                thermal_context = " and sensor placement/ambient alignment are undocumented"
            else:
                thermal_context = ""
            row["thermal_validation_evidence"] = (
                "temperature traces remain descriptive, but heat input is setup-confounded"
                + thermal_context
            )
            row["recommended_thesis_use"] = (
                "experimental-nonideality case study; qualitative temperature-response diagnostic"
            )
    return reference_group, split_threshold, reference_median, reference_iqr


def write_csv_union(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def manifest_rows(experiments: list[StandardExperiment], audits: list[dict[str, object]]) -> list[dict[str, object]]:
    audit_by_id = {str(row["dataset_id"]): row for row in audits}
    rows: list[dict[str, object]] = []
    for experiment in experiments:
        audit = audit_by_id[experiment.metadata.dataset_id]
        rows.append(
            {
                "dataset_id": experiment.metadata.dataset_id,
                "full_path": str(experiment.metadata.path.resolve()),
                "filename": experiment.metadata.path.name,
                "cell_id": experiment.metadata.cell_id,
                "nominal_test_temperature": experiment.metadata.nominal_temperature,
                "test_type_signal_verified": experiment.test_type,
                "date_encoded_or_metadata": experiment.metadata.encoded_date,
                "sample_count": audit["sample_count"],
                "duration_s": audit["estimated_duration_s"],
                "duration_status": "provisional: no raw time channel",
                "available_channels": audit["available_channels"],
                "header_present": experiment.header_present,
                "raw_current_sign": experiment.raw_current_sign,
                "model_loader_mapping": "I_model=-I_raw; raw file unchanged",
            }
        )
    return rows


def save_experiment_figure(experiment: StandardExperiment) -> Path:
    time_minutes = experiment.estimated_time_s / 60.0
    stride = max(1, time_minutes.size // 12000)
    shown = np.arange(0, time_minutes.size, stride)
    figure, axes = plt.subplots(4, 1, figsize=(12.0, 9.5), sharex=True, constrained_layout=True)
    axes[0].plot(time_minutes[shown], experiment.numeric["Current(A)"][shown], color="#0072B2", lw=0.8)
    axes[0].axhline(0.0, color="black", lw=0.6)
    axes[0].set_ylabel("Raw current [A]")
    axes[1].plot(time_minutes[shown], experiment.numeric["Voltage(V)"][shown], color="#D55E00", lw=0.8)
    axes[1].set_ylabel("Voltage [V]")
    temp_columns = [name for name in experiment.numeric if name.startswith("Temp")]
    ambient = [name for name in temp_columns if "ambient" in experiment.metadata.temperature_roles.get(name, "") or "chamber" in experiment.metadata.temperature_roles.get(name, "")]
    battery = [name for name in temp_columns if name not in ambient]
    for name in battery:
        axes[2].plot(time_minutes[shown], experiment.numeric[name][shown], lw=0.8, label=f"{name}: {experiment.metadata.temperature_roles.get(name, 'unknown')}")
    if battery:
        axes[2].legend(fontsize=8)
    else:
        axes[2].text(0.5, 0.5, "Battery temperature unavailable", transform=axes[2].transAxes, ha="center")
    axes[2].set_ylabel("Battery T [degC]")
    for name in ambient:
        axes[3].plot(time_minutes[shown], experiment.numeric[name][shown], color="#009E73", lw=0.8, label=name)
    if ambient:
        axes[3].legend(fontsize=8)
    else:
        axes[3].text(0.5, 0.5, "Chamber/ambient temperature unavailable", transform=axes[3].transAxes, ha="center")
    axes[3].set(xlabel="Provisional elapsed time from capacity/current [min]", ylabel="Ambient T [degC]")
    for axis in axes:
        axis.grid(alpha=0.22)
    figure.suptitle(f"{experiment.metadata.dataset_id}: raw signals (timebase provisional)")
    path = FIGURE_DIRECTORY / f"{experiment.metadata.dataset_id}.png"
    figure.savefig(path, dpi=170)
    plt.close(figure)
    return path


def save_resistance_figure(
    experiments: list[StandardExperiment],
    transitions: list[dict[str, object]],
    reference_group: set[str],
    split_threshold: float,
) -> None:
    ids = [experiment.metadata.dataset_id for experiment in experiments]
    values = [
        [float(row["r_app_ohm"]) for row in transitions if row["dataset_id"] == dataset_id]
        for dataset_id in ids
    ]
    figure, axis = plt.subplots(figsize=(12.5, 6.5), constrained_layout=True)
    plot = axis.boxplot(values, tick_labels=ids, patch_artist=True, showmeans=True)
    for patch, dataset_id in zip(plot["boxes"], ids, strict=True):
        patch.set_facecolor("#56B4E9" if dataset_id in reference_group else "#E69F00")
        patch.set_alpha(0.65)
    for index, dataset_id in enumerate(ids, start=1):
        rows = [row for row in transitions if row["dataset_id"] == dataset_id]
        for row in rows:
            marker = "v" if row["direction"] == "discharge" else "^"
            axis.scatter(index, float(row["r_app_ohm"]), marker=marker, color="black", s=28, zorder=4)
    axis.axhline(split_threshold, color="#CC79A7", ls="--", label="Data-driven largest-log-gap split")
    axis.set(ylabel="Apparent rest-to-load series resistance [ohm]", title="Own TRC18650 rest-to-load diagnostic (not intrinsic cell R0)")
    axis.tick_params(axis="x", rotation=28)
    axis.grid(alpha=0.22, axis="y")
    axis.legend()
    figure.savefig(CROSS_RESISTANCE_PATH, dpi=180)
    plt.close(figure)


def format_value(value: object, digits: int = 6) -> str:
    if isinstance(value, (float, np.floating)):
        return "not assessable" if not np.isfinite(float(value)) else f"{float(value):.{digits}f}"
    return str(value)


def write_report(
    experiments: list[StandardExperiment],
    audits: list[dict[str, object]],
    transitions: list[dict[str, object]],
    reference_group: set[str],
    split_threshold: float,
    reference_median: float,
    reference_iqr: float,
    supporting_files: list[Path],
) -> None:
    audit_by_id = {str(row["dataset_id"]): row for row in audits}
    strongest = next(iter(reference_group))
    affected = [row for row in audits if row["dataset_id"] not in reference_group]
    lines = [
        "# Own TRC18650 experimental-data quality audit",
        "",
        "This is a read-only data-quality audit. No SPM/OCP/electrical/thermal parameter was fitted, no Panasonic parameter was applied, and no raw voltage was corrected.",
        "",
        "## Discovery and standardized representation",
        "",
        f"Seven primary Chroma/TRC18650 CSVs were found under `{SOURCE_DIRECTORY}` and audited individually. {len(supporting_files)} chamber/controller CSVs were also found under `{CONTROL_DIRECTORY}`; they are supporting temperature logs, not distinct battery experiments, and the available columns do not provide a defensible absolute alignment to the Chroma files.",
        "",
        "Every primary CSV contains voltage, current, step capacity, signed total capacity, step number/mode, and at least two temperature channels. The reusable loader exposes canonical voltage/current views while preserving the raw fields and their sign. None contains elapsed time or a timestamp. Therefore all duration, dt, current integration, and time-domain plots are explicitly **provisional**, using a uniform sampling interval inferred from active-step capacity/current. Raw timestamp duplicates, reversals, and long gaps are **not assessable**. The raw data remain unchanged.",
        "",
        "Because the provisional timebase is inferred from the same capacity/current records, agreement between reconstructed current integration and the Ah counters is a consistency check, **not an independent capacity measurement**.",
        "",
        "Signal evidence establishes `I_raw<0` for discharge and `I_raw>0` for charge; a future model loader must record `I_model=-I_raw` without overwriting source current.",
        "",
        "The `cycle_cell1_ambient` and `cycle_cell2_undeclared` experiments were simultaneous and shared the measured air thermocouple `Temp7` (22.4896–23.4919 degC). This channel is associated with both tests as experimental context but is not copied into the cell-2 arrays. Different record counts and absent timestamps prevent sample-by-sample alignment. `Temp3` and `Temp4` remain separate measurements of undocumented placement.",
        "",
        "## Suitability summary",
        "",
        "| Dataset | Electrical calibration | Electrical validation | Thermal calibration | Thermal validation | Recommended thesis use |",
        "|---|---|---|---|---|---|",
    ]
    for experiment in experiments:
        row = audit_by_id[experiment.metadata.dataset_id]
        lines.append(
            f"| {row['dataset_id']} | {row['electrical_calibration_suitability']} | {row['electrical_validation_suitability']} | {row['thermal_calibration_suitability']} | {row['thermal_validation_suitability']} | {row['recommended_thesis_use']} |"
        )
    lines.extend(
        [
            "",
            "## Apparent rest-to-load series resistance",
            "",
            "Common rule: use adjacent records with `abs(I_pre)<=0.02 A`, `abs(I_loaded)>=2.5 A`, and `abs(delta I)>=2.5 A`; then calculate `R_app=abs(delta V)/abs(delta I)`. Both discharge and charge onsets are retained and labelled. This is the apparent response of the observed cell/fixture/wiring/measurement/electrochemical path, **not intrinsic cell resistance or model R0**.",
            "",
            f"The lowest group was selected from experiment medians by the largest adjacent gap in log(R_app), not by a universal resistance threshold. The data-driven split is `{split_threshold:.9f} ohm`. Reference-group median/IQR are `{reference_median:.9f}/{reference_iqr:.9f} ohm`.",
            "",
            "| Dataset | N | Median [ohm] | Mean [ohm] | Std [ohm] | Range [ohm] | IQR [ohm] | Excess [ohm] | Ratio | Interpretation |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|",
        ]
    )
    for experiment in experiments:
        row = audit_by_id[experiment.metadata.dataset_id]
        lines.append(
            f"| {row['dataset_id']} | {row['r_app_transition_count']} | {format_value(row['r_app_median_ohm'], 9)} | {format_value(row['r_app_mean_ohm'], 9)} | {format_value(row['r_app_std_ohm'], 9)} | {format_value(row['r_app_min_ohm'], 9)}-{format_value(row['r_app_max_ohm'], 9)} | {format_value(row['r_app_iqr_ohm'], 9)} | {format_value(row['r_app_excess_vs_reference_ohm'], 9)} | {format_value(row['r_app_ratio_vs_reference'], 3)}x | {row['resistance_group']}; {row['r_app_variability']} |"
        )
    lines.extend(
        [
            "",
            "The elevated group is separated from the reference by a factor of roughly five to six. Within the repeated pulse records, every transition remains on the elevated side of the data-driven split; this is evidence of a persistent path condition rather than one isolated jump. Attribution is not unique, so the report uses **evidence consistent with parasitic series resistance** and does not subtract any inferred excess.",
            "",
            "## Per-dataset report cards",
        ]
    )
    for experiment in experiments:
        row = audit_by_id[experiment.metadata.dataset_id]
        temp_columns = [name for name in experiment.numeric if name.startswith("Temp")]
        temp_descriptions = ", ".join(
            f"{name}: initial/final {format_value(row.get(name+'_initial_c'))}/{format_value(row.get(name+'_final_c'))}, range {format_value(row.get(name+'_min_c'))}-{format_value(row.get(name+'_max_c'))} degC, apparent step {format_value(row.get(name+'_quantization_c'), 4)} degC, max adjacent jump {format_value(row.get(name+'_max_adjacent_jump_c'))} degC, robust jump flags {row.get(name+'_jump_flag_count')} above {format_value(row.get(name+'_jump_screen_threshold_c'))} degC, missing {row.get(name+'_nonfinite_count')} ({experiment.metadata.temperature_roles.get(name, 'unknown')})"
            for name in temp_columns
        )
        lines.extend(
            [
                "",
                f"### {row['dataset_id']}",
                "",
                f"- **Dataset:** `{row['full_path']}`",
                f"- **Cell / nominal test temperature / signal-verified type:** {row['cell_id']} / {row['nominal_test_temperature']} / {row['test_type_signal_verified']}.",
                f"- **Time integrity:** no raw time; provisional dt `{row['estimated_dt_s']:.9f} s`, duration `{row['estimated_duration_s']:.3f} s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.",
                f"- **Current sign:** {row['raw_current_sign']} ({row['current_sign_confidence']} confidence); range `{row['current_min_a']:.5f}` to `{row['current_max_a']:.5f} A`.",
                f"- **Capacity consistency:** main discharge/charge counters `{format_value(row['main_discharge_counter_capacity_ah'])}/{format_value(row['main_charge_counter_capacity_ah'])} Ah`; provisional integrations `{format_value(row['main_discharge_integrated_capacity_ah'])}/{format_value(row['main_charge_integrated_capacity_ah'])} Ah`. Whole-file provisional signed/absolute integrals are `{format_value(row['integrated_signed_charge_ah_provisional'])}/{format_value(row['integrated_absolute_charge_ah_provisional'])} Ah`, versus signed/absolute counter changes `{format_value(row['counter_signed_change_ah'])}/{format_value(row['counter_absolute_throughput_ah'])} Ah`; differences `{format_value(row['signed_integration_minus_counter_ah'])}/{format_value(row['absolute_integration_minus_counter_ah'])} Ah`. Complete sequence `{row['capacity_sequence_appears_complete']}`; SOC status: `{row['soc_label']}`.",
                f"- **Voltage quality:** range `{row['voltage_min_v']:.4f}-{row['voltage_max_v']:.4f} V`; main discharge start/end `{format_value(row['main_discharge_start_voltage_v'], 4)}/{format_value(row['main_discharge_end_voltage_v'], 4)} V` with preceding/following rest `{row['main_discharge_preceded_by_rest']}/{row['main_discharge_followed_by_rest']}`; largest adjacent jump `{row['max_voltage_jump_v']:.4f} V`; largest jump with `abs(delta I)<0.05 A` `{row['max_voltage_jump_at_small_current_change_v']:.4f} V`; near-extreme repeat count `{row['voltage_extreme_near_repeat_count']}`; inconsistent transition signs `{row['transition_sign_inconsistent_count']}`.",
                f"- **Thermal quality:** {temp_descriptions}; active runs with net heating `{row['temperature_active_runs_showing_net_heating']}`, rest runs with net cooling `{row['temperature_rest_runs_showing_net_cooling']}`. Sensor-pair offset/correlation/lag are `{format_value(row['temperature_pair_mean_offset_c'])} degC / {format_value(row['temperature_pair_correlation'])} / {row['temperature_pair_lag_samples']} samples`.",
                f"- **Apparent series resistance:** N `{row['r_app_transition_count']}`, median `{format_value(row['r_app_median_ohm'], 9)} ohm`; `{row['resistance_group']}`. {row['evidence_of_parasitic_setup_effects']}.",
                f"- **Electrical calibration:** {row['electrical_calibration_suitability']} — {row['electrical_calibration_evidence']}.",
                f"- **Electrical validation:** {row['electrical_validation_suitability']} — {row['electrical_validation_evidence']}.",
                f"- **Thermal calibration:** {row['thermal_calibration_suitability']} — {row['thermal_calibration_evidence']}.",
                f"- **Thermal validation:** {row['thermal_validation_suitability']} — {row['thermal_validation_evidence']}.",
                f"- **Recommended thesis use:** {row['recommended_thesis_use']}.",
            ]
        )
    affected_names = ", ".join(str(row["dataset_id"]) for row in affected)
    cell1_thermal = ", ".join(
        str(row["dataset_id"])
        for row in affected
        if row["cell_id"] == "cell 1" and "Temp7" in str(row["temperature_channels"])
    )
    lines.extend(
        [
            "",
            "## Methodological interpretation",
            "",
            f"1. **Strongest own experiment for quantitative model comparison:** `{strongest}`. It is the sole lowest-resistance-group experiment and has complete discharge/charge capacity plus documented cell-surface and ambient channels. Quantitative use should be restricted to voltage/capacity until an actual acquisition timebase is recovered; thermal use remains qualitative.",
            f"2. **Experiments showing setup-related voltage evidence:** {affected_names}. Their relative separation is consistent with added series-path resistance, but the data cannot uniquely locate it in contacts, fixture, wiring, sensing, or cell state.",
            "3. **Consistency of affected behavior:** repeated pulse transitions are internally consistent; heated cell-1 cycles occupy the same elevated band. Cell-2 values are also elevated but show more between-experiment spread, so attribution across physical cells remains cautious.",
            f"4. **Affected datasets retaining useful thermal information:** {cell1_thermal} retain two documented cell-surface channels plus ambient and show heating/cooling structure. They are useful descriptively, not as quantitative heat-ground-truth. Cell-2 temperature traces remain descriptive; Temp3/Temp4 placement is undocumented. The simultaneous Aire 22–24 °C cell-2 run shares the measured Temp7 air record, but cross-file sample alignment is not established.",
            "5. **Best experimental-nonideality examples:** the two pulse parts provide repeated evidence across ten transitions; the 30/35 degC cell-1 cycles show that the elevated step also appears in nominally steady cycle tests.",
            "6. **Framework reuse:** yes. Channel aliases/configuration, explicit current-sign metadata, time-source status, integrity/capacity/temperature checks, adjacent rest-to-load extraction, and data-driven cross-experiment grouping can be reused. New datasets with measured time improve the audit automatically without changing the methodology.",
            "",
            "## Clean, affected, and ambiguous evidence",
            "",
            f"- **Relative clean reference:** `{strongest}` is the strongest clean candidate in this collection, not a certified intrinsic-resistance experiment. Only two comparable onsets are present, and the absence of raw time prevents a full dynamic-quality certification.",
            "- **Clearly affected for electrical use:** the two cell-1 pulse parts and the 30/35 degC cell-1 cycles show persistent, approximately five-to-six-fold elevated apparent transition resistance relative to the reference. They are rejected for electrical identification.",
            "- **Attribution remains ambiguous:** both cell-2 cycles are also elevated and therefore unsuitable for electrical identification, but comparison with cell 1 cannot separate fixture resistance from genuine cell-to-cell or temperature differences. Temp3/Temp4 placement is undocumented. The Aire 22–24 °C test has a confirmed shared measured Temp7 air record, without sample-level alignment.",
            "",
            "The reported temperature `apparent step` is the modal nonzero adjacent digital increment, not a certified sensor resolution. The robust jump screen flags adjacent changes above `max(0.5 degC, median + 10 x 1.4826 x MAD)` for review; a flag alone is not an automatic rejection.",
            "",
            "## Supporting controller logs",
            "",
            "These files were discovered but are not independently classified as battery experiments because voltage/current/cell identity and a defensible alignment key are absent:",
            "",
        ]
    )
    lines.extend(f"- `{path}`" for path in supporting_files)
    lines.extend(
        [
            "",
            "## Stop point",
            "",
            "No SPM, OCP, R0, diffusion, or thermal parameter was fitted; no voltage correction or model-accuracy calculation was performed.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    FIGURE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    experiments = load_all()
    transitions_by_id = {
        experiment.metadata.dataset_id: extract_rest_to_load_transitions(experiment)
        for experiment in experiments
    }
    transitions = [
        row
        for experiment in experiments
        for row in transitions_by_id[experiment.metadata.dataset_id]
    ]
    audits = [
        audit_experiment(experiment, transitions_by_id[experiment.metadata.dataset_id])
        for experiment in experiments
    ]
    reference_group, split_threshold, reference_median, reference_iqr = classify_and_augment(
        audits, transitions
    )
    write_csv_union(MANIFEST_PATH, manifest_rows(experiments, audits))
    write_csv_union(QUALITY_PATH, audits)
    write_csv_union(TRANSITION_PATH, transitions)
    for experiment in experiments:
        save_experiment_figure(experiment)
    save_resistance_figure(experiments, transitions, reference_group, split_threshold)
    supporting_files = sorted(CONTROL_DIRECTORY.glob("*.csv")) if CONTROL_DIRECTORY.is_dir() else []
    write_report(
        experiments,
        audits,
        transitions,
        reference_group,
        split_threshold,
        reference_median,
        reference_iqr,
        supporting_files,
    )
    print(f"Audited {len(experiments)} own TRC18650 files; {len(transitions)} comparable transitions.")
    print(f"Reference group={sorted(reference_group)}; median={reference_median:.9f} ohm; split={split_threshold:.9f} ohm")
    print(f"Outputs written under {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
