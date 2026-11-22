from __future__ import annotations

import numpy as np
import pytest

from analysis.battery_dataset_assessment import (
    MetricRecord,
    MetricRole,
    NonidealityFlag,
    TimebaseStatus,
    TransferEvidence,
    TransferLevel,
    convert_current_for_model,
    evaluate_relative_nonideality_gate,
    evaluate_transfer_gate,
    request_automatic_voltage_correction,
    require_thermal_dynamic_support,
    resolve_timebase,
    screen_ocp_domain,
    validate_voltage_curve_label,
)


def test_raw_current_is_never_overwritten_or_aliased() -> None:
    raw = np.asarray([0.0, -2.9, 1.4])
    original = raw.copy()
    model = convert_current_for_model(raw, "discharge_negative")
    np.testing.assert_array_equal(raw, original)
    np.testing.assert_array_equal(model, -original)
    assert not np.shares_memory(raw, model)


def test_provisional_time_cannot_be_marked_measured() -> None:
    status = resolve_timebase(measured_available=False, provisional_available=True)
    assert status == TimebaseStatus.INFERRED_PROVISIONAL
    assert status != TimebaseStatus.MEASURED


def test_r_app_cannot_be_labelled_intrinsic_r0() -> None:
    with pytest.raises(ValueError, match="R_app cannot be labelled"):
        MetricRecord(
            name="R_app intrinsic R0",
            value=0.03,
            unit="ohm",
            metric_role=MetricRole.DATA_QUALITY_DIAGNOSTIC,
            evidence_basis="Adjacent transition.",
        )


def test_voltage_correction_is_never_automatic() -> None:
    with pytest.raises(RuntimeError, match="prohibited"):
        request_automatic_voltage_correction()


def test_thermal_dynamic_metrics_blocked_without_measured_time() -> None:
    assert (
        require_thermal_dynamic_support(TimebaseStatus.INFERRED_PROVISIONAL, True)
        == "THERMAL DYNAMIC ACCURACY NOT ASSESSABLE"
    )
    with pytest.raises(ValueError, match="THERMAL DYNAMIC ACCURACY NOT ASSESSABLE"):
        MetricRecord(
            name="temperature_time_rmse",
            value=0.5,
            unit="degC",
            metric_role=MetricRole.MODEL_ACCURACY,
            evidence_basis="Provisional time only.",
            domain="thermal_dynamic",
            timebase_status=TimebaseStatus.INFERRED_PROVISIONAL,
        )


def test_loaded_1c_voltage_cannot_be_labelled_ocv() -> None:
    with pytest.raises(ValueError, match="cannot be labelled OCV"):
        validate_voltage_curve_label("loaded_1c", "TRC OCV curve")


def test_unsupported_ocp_states_are_flagged_not_clipped() -> None:
    original = np.asarray([0.17, 0.18, 0.50, 0.96, 0.97])
    screen = screen_ocp_domain(original, 0.18, 0.96)
    np.testing.assert_array_equal(screen.values, original)
    np.testing.assert_array_equal(screen.valid_mask, [False, True, True, True, False])
    assert screen.invalid_count == 2
    assert screen.values[0] == 0.17 and screen.values[-1] == 0.97


def test_failed_transfer_never_triggers_model_fitting() -> None:
    decision = evaluate_transfer_gate(
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
    assert decision.automatic_fit_triggered is False
    assert decision.levels == (
        TransferLevel.LEVEL_1,
        TransferLevel.LEVEL_2,
        TransferLevel.LEVEL_3,
        TransferLevel.LEVEL_4,
    )
    assert "do not fit" in decision.next_action.lower()


def test_nonideality_gate_uses_relative_data_grouping() -> None:
    rows = [
        {"dataset_id": "clean", "r_app_transition_count": 3, "r_app_median_ohm": 0.030, "r_app_iqr_ohm": 0.002},
        {"dataset_id": "peer", "r_app_transition_count": 3, "r_app_median_ohm": 0.032, "r_app_iqr_ohm": 0.002},
        {"dataset_id": "elevated", "r_app_transition_count": 4, "r_app_median_ohm": 0.200, "r_app_iqr_ohm": 0.010},
    ]
    results = {item.dataset: item for item in evaluate_relative_nonideality_gate(rows, "clean")}
    assert results["clean"].flag.value == "NO RELATIVE ELEVATION DETECTED"
    assert results["peer"].flag.value == "NO RELATIVE ELEVATION DETECTED"
    assert results["elevated"].flag.value == "EVIDENCE OF ELEVATED SERIES-PATH RESPONSE"
    assert "largest adjacent gap" in results["elevated"].reference_group_provenance


def test_nonideality_gate_can_select_discharge_only_fields() -> None:
    rows = [
        {
            "dataset_id": "clean",
            "r_app_transition_count": 2,
            "r_app_median_ohm": 0.025,
            "r_app_iqr_ohm": 0.010,
            "discharge_r_app_count": 1,
            "discharge_r_app_median_ohm": 0.036,
            "discharge_r_app_iqr_ohm": 0.0,
        },
        {
            "dataset_id": "elevated",
            "r_app_transition_count": 2,
            "r_app_median_ohm": 0.040,
            "r_app_iqr_ohm": 0.010,
            "discharge_r_app_count": 1,
            "discharge_r_app_median_ohm": 0.240,
            "discharge_r_app_iqr_ohm": 0.0,
        },
    ]
    results = {
        item.dataset: item
        for item in evaluate_relative_nonideality_gate(
            rows,
            "clean",
            transition_count_field="discharge_r_app_count",
            median_field="discharge_r_app_median_ohm",
            iqr_field="discharge_r_app_iqr_ohm",
            minimum_transition_count=1,
            metric_label="R_app,desc",
        )
    }
    assert results["clean"].median_r_app_ohm == 0.036
    assert results["elevated"].median_r_app_ohm == 0.240
    assert results["elevated"].ratio_to_reference == 0.240 / 0.036
    assert results["elevated"].transition_count == 1
    assert results["elevated"].flag == NonidealityFlag.ELEVATED_SERIES_PATH
    assert "R_app,desc" in results["elevated"].reference_group_provenance
