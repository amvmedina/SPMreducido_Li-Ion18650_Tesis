from __future__ import annotations

from analysis.battery_dataset_assessment import NonidealityFlag, TransferLevel
from analysis.trc18650_automatic_audit import (
    REFERENCE_DATASET_ID,
    TRC_REFERENCE_R_APP_DESC_OHM,
    build_recommendations,
    build_transfer_decision,
    decide_dataset_use,
    evaluate_trc_discharge_nonideality,
)


def test_verified_directional_reference_is_frozen_for_the_audit() -> None:
    assert TRC_REFERENCE_R_APP_DESC_OHM == 0.036383074111115


def test_trc_route_uses_discharge_values_not_combined_values() -> None:
    rows = [
        {
            "dataset_id": REFERENCE_DATASET_ID,
            "r_app_transition_count": 2,
            "r_app_median_ohm": 0.020,
            "r_app_iqr_ohm": 0.010,
            "discharge_r_app_count": 1,
            "discharge_r_app_median_ohm": TRC_REFERENCE_R_APP_DESC_OHM,
            "discharge_r_app_iqr_ohm": 0.0,
        },
        {
            "dataset_id": "cycle_cell1_30c",
            "r_app_transition_count": 2,
            "r_app_median_ohm": 0.050,
            "r_app_iqr_ohm": 0.010,
            "discharge_r_app_count": 1,
            "discharge_r_app_median_ohm": 0.24060721062618612,
            "discharge_r_app_iqr_ohm": 0.0,
        },
    ]
    result = evaluate_trc_discharge_nonideality(rows)
    elevated = result["cycle_cell1_30c"]
    assert elevated.median_r_app_ohm == 0.24060721062618612
    assert elevated.ratio_to_reference == (
        0.24060721062618612 / TRC_REFERENCE_R_APP_DESC_OHM
    )
    assert elevated.transition_count == 1
    assert elevated.flag == NonidealityFlag.ELEVATED_SERIES_PATH


def _classification(
    dataset_id: str,
    *,
    elevated: bool,
    measured_time: bool = False,
    calibration: str = "NOT SUITABLE",
) -> dict[str, object]:
    return {
        "dataset_id": dataset_id,
        "measured_time_available": measured_time,
        "nonideality_flag": (
            NonidealityFlag.ELEVATED_SERIES_PATH.value
            if elevated
            else NonidealityFlag.NO_RELATIVE_ELEVATION.value
        ),
        "electrical_calibration_status": calibration,
    }


def test_clean_reference_is_characterization_not_calibration() -> None:
    decision = decide_dataset_use(
        dataset_id=REFERENCE_DATASET_ID,
        test_type="combined cycle",
        capacity_complete=True,
        time_measured=False,
        nonideality_flag=NonidealityFlag.NO_RELATIVE_ELEVATION,
    )
    assert decision["characterization_status"] == "SUITABLE WITH CAVEAT"
    assert decision["electrical_calibration_status"] == "NOT SUITABLE"
    assert "diagnóstico de transferencia" in decision["recommended_use"]


def test_elevated_pulse_is_diagnostic_only() -> None:
    decision = decide_dataset_use(
        dataset_id="pulse_cell1_part1",
        test_type="pulse discharge/rest",
        capacity_complete=False,
        time_measured=False,
        nonideality_flag=NonidealityFlag.ELEVATED_SERIES_PATH,
    )
    assert decision["characterization_status"] == "QUALITATIVE / DIAGNOSTIC ONLY"
    assert decision["electrical_validation_status"] == "QUALITATIVE / DIAGNOSTIC ONLY"
    assert decision["thermal_calibration_status"] == "NOT SUITABLE"


def test_transfer_gate_never_fits_and_reports_all_required_levels() -> None:
    decision = build_transfer_decision()
    assert decision.automatic_fit_triggered is False
    assert decision.levels == (
        TransferLevel.LEVEL_1,
        TransferLevel.LEVEL_2,
        TransferLevel.LEVEL_3,
        TransferLevel.LEVEL_4,
    )


def test_recommendations_prioritize_recovery_before_new_experiments() -> None:
    rows = (
        _classification(REFERENCE_DATASET_ID, elevated=False),
        _classification("pulse_cell1_part1", elevated=True),
    )
    recommendations = build_recommendations(rows)
    assert recommendations[0]["priority"] == 0
    assert recommendations[0]["action_id"] == "recover_documentation_and_verify_setup"
    assert recommendations[1]["action_id"] == "trc_low_rate_reference"
    assert "40" in str(recommendations[1]["active_execution_time_h"])


def test_recommendations_do_not_call_rapp_intrinsic_resistance() -> None:
    rows = (
        _classification(REFERENCE_DATASET_ID, elevated=False),
        _classification("cycle_cell1_30c", elevated=True),
    )
    text = " ".join(str(value) for row in build_recommendations(rows) for value in row.values()).lower()
    assert "r_app" in text or "apparent" in text
    assert "intrinsic r0" not in text
    assert "internal resistance estimate" not in text
