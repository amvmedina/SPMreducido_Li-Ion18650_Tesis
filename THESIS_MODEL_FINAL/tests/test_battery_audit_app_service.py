from __future__ import annotations

from pathlib import Path

import pytest

from analysis.battery_audit_app_service import (
    CellContext,
    ROLE_CHARACTERIZATION,
    ROLE_DYNAMIC,
    ROLE_LOW_RATE,
    SETUP_OTHER,
    SETUP_REFERENCE,
    UploadedFileSpec,
    audit_uploaded_files,
    build_panasonic_reference_audit,
    build_trc_case_study,
    inspect_csv_bytes,
    reference_model_payload,
)


def _chroma_csv(*, loaded_voltage: float) -> bytes:
    rows = [
        "Description,Step No.,Voltage(V),Current(A),Capacity(Ah),Mode,Total Capacity(Ah),Temp1,Temp7",
        ",1,4.2,0.0,0.0,Rest,0.0,25.0,25.0",
        f",2,{loaded_voltage},-3.0,0.0,CC,0.0,25.1,25.0",
        ",2,3.9,-3.0,0.02,CC,-0.02,25.2,25.0",
        ",2,3.7,-3.0,0.04,CC,-0.04,25.4,25.0",
        ",2,3.5,-3.0,0.06,CC,-0.06,25.6,25.1",
        ",2,3.3,-3.0,0.08,CC,-0.08,25.8,25.1",
        ",2,3.0,-3.0,0.10,CC,-0.10,26.0,25.1",
        ",3,3.1,0.0,0.10,Rest,-0.10,26.0,25.1",
    ]
    return ("\n".join(rows) + "\n").encode("utf-8")


def _spec(dataset_id: str, voltage: float, setup: str) -> UploadedFileSpec:
    return UploadedFileSpec(
        filename=f"{dataset_id}.csv",
        content=_chroma_csv(loaded_voltage=voltage),
        dataset_id=dataset_id,
        cell_id=dataset_id,
        condition="condición declarada",
        setup_role=setup,
        experimental_roles=(ROLE_CHARACTERIZATION,),
        temperature_roles={
            "Temp1": "cell surface middle",
            "Temp7": "chamber/ambient",
        },
    )


def test_uploaded_csv_is_structurally_inspected_without_project_paths() -> None:
    preview = inspect_csv_bytes("new_cell.csv", _chroma_csv(loaded_voltage=4.08))
    assert preview["header_detected"] is True
    assert preview["required_electrical_columns_detected"] is True
    assert preview["sample_count_estimate"] == 8
    assert preview["temperature_columns"] == ("Temp1", "Temp7")


def test_reference_parameter_table_contains_descriptions_and_physical_units() -> None:
    parameters = reference_model_payload()["parameters"]
    assert len(parameters) == 14
    assert all(tuple(row) == ("Parámetro", "Descripción", "Valor", "Unidad") for row in parameters)
    assert all(str(row["Descripción"]).strip() for row in parameters)
    resistance = next(row for row in parameters if row["Parámetro"] == "R₀")
    assert resistance["Unidad"] == "Ω"


def test_uploaded_audit_uses_calculated_directional_response() -> None:
    result = audit_uploaded_files(
        (
            _spec("reference", 4.08, SETUP_REFERENCE),
            _spec("comparison", 3.48, SETUP_OTHER),
        ),
        CellContext(),
    )
    summaries = {str(row["dataset_id"]): row for row in result.file_summaries}
    assert summaries["reference"]["R_app,desc [ohm]"] == pytest.approx(0.04)
    assert summaries["comparison"]["R_app,desc [ohm]"] == pytest.approx(0.24)
    assert summaries["comparison"]["razón relativa"] == pytest.approx(6.0)
    check = next(row for row in result.checks if row["código"] == "E04")
    assert check["estado"] == "REVISAR"
    assert "otro montaje: 0.24 Ω (6 veces la referencia)" in str(check["evidencia"])
    assert "no identifica por sí sola" in result.report_markdown.lower()


def test_uploaded_audit_does_not_invent_a_reference() -> None:
    result = audit_uploaded_files((_spec("only", 4.08, SETUP_OTHER),), CellContext())
    check = next(row for row in result.checks if row["código"] == "E04")
    assert check["estado"] == "NO EVALUADO"
    assert "no hay dos montajes comparables" in str(check["evidencia"]).lower()
    assert check["acción requerida"] == ""


def test_case_study_report_is_generated_from_engine_results() -> None:
    source_directory = Path(r"C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data")
    if not source_directory.is_dir():
        pytest.skip("Configured TRC case-study files are not installed.")
    result = build_trc_case_study()
    assert len(result.file_summaries) == 8
    assert {row["condición"] for row in result.file_summaries} == {
        "Rionegro · Aire 22–24 °C", "Cali · Aire 29–31 °C", "Valledupar · Aire 33–37 °C"
    }
    visible = str(result.file_summaries)
    assert "cycle_cell" not in visible
    assert "pulse_cell" not in visible
    statuses = {row["código"]: row["estado"] for row in result.checks}
    assert statuses["E06"] == "FALTANTE"
    for row in result.checks:
        message = f"{row['evidencia']} {row['acción requerida']}"
        assert "[" not in message
        assert "valores no finitos" not in message
        assert "conflictos de reutilización" not in message
        assert "partición defendible" not in message


def test_low_rate_and_pulse_requirements_follow_declared_inputs() -> None:
    ordinary = _spec("ordinary", 4.08, SETUP_OTHER)
    baseline = audit_uploaded_files((ordinary,), CellContext())
    baseline_states = {row["código"]: row["estado"] for row in baseline.checks}
    assert baseline_states["E05"] == "FALTANTE"
    assert baseline_states["E06"] == "FALTANTE"
    low_rate_check = next(row for row in baseline.checks if row["código"] == "E05")
    pulse_check = next(row for row in baseline.checks if row["código"] == "E06")
    assert low_rate_check["acción requerida"] == "Realizar una descarga C/20."
    assert "No se incluyó un ensayo de pulsos" in str(pulse_check["evidencia"])
    assert "HPPC" in str(pulse_check["acción requerida"])

    low_rate = UploadedFileSpec(
        **{**ordinary.__dict__, "dataset_id": "low_rate", "experimental_roles": (ROLE_LOW_RATE,)}
    )
    with_low_rate = audit_uploaded_files((ordinary, low_rate), CellContext())
    states = {row["código"]: row["estado"] for row in with_low_rate.checks}
    assert states["E05"] == "REVISAR"

    pulse = UploadedFileSpec(
        **{**ordinary.__dict__, "dataset_id": "pulse", "experimental_roles": (ROLE_DYNAMIC,)}
    )
    with_pulse = audit_uploaded_files((pulse,), CellContext())
    pulse_check = next(row for row in with_pulse.checks if row["código"] == "E06")
    assert pulse_check["estado"] == "NO APTO"


def test_commercial_cell_name_does_not_change_decisions() -> None:
    first = _spec("same_data", 4.08, SETUP_OTHER)
    second = UploadedFileSpec(**{**first.__dict__, "cell_id": "another commercial name"})
    states_first = [(row["código"], row["estado"]) for row in audit_uploaded_files((first,), CellContext()).checks]
    states_second = [(row["código"], row["estado"]) for row in audit_uploaded_files((second,), CellContext()).checks]
    assert states_first == states_second


def test_panasonic_and_trc_use_the_same_check_identifiers(monkeypatch: pytest.MonkeyPatch) -> None:
    source_directory = Path(r"C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data")
    if not source_directory.is_dir():
        pytest.skip("Configured reference files are not installed.")
    import analysis.battery_audit_visualization as visualization

    monkeypatch.setattr(visualization, "load_panasonic_plot", lambda experiment, variable: object())
    monkeypatch.setattr(visualization, "plot_is_finite", lambda plot: True)
    trc = build_trc_case_study()
    panasonic = build_panasonic_reference_audit()
    assert [row["código"] for row in trc.checks] == [f"E{index:02d}" for index in range(1, 9)]
    assert [row["código"] for row in panasonic.checks] == [f"E{index:02d}" for index in range(1, 9)]
    panasonic_e04 = next(row for row in panasonic.checks if row["código"] == "E04")
    assert panasonic_e04["estado"] == "NO EVALUADO"
    assert panasonic_e04["acción requerida"] == ""
