"""Interfaz Streamlit para evaluar conjuntos de datos de baterías."""

from __future__ import annotations

import html
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

from analysis.battery_audit_app_service import (
    CellContext, ROLE_CHARACTERIZATION, ROLE_DYNAMIC, ROLE_LOW_RATE, ROLE_THERMAL,
    ROLE_UNKNOWN, ROLE_VALIDATION, SETUP_OTHER, SETUP_REFERENCE, SETUP_UNKNOWN,
    STANDARD_CHROMA_COLUMNS, UploadedFileSpec, audit_uploaded_files,
    audit_criteria_payload, build_panasonic_reference_audit, build_trc_case_study, inspect_csv_bytes,
    reference_model_payload,
)
from analysis.battery_audit_visualization import (
    PANASONIC_CATALOG, TRC_OPTIONS, ExperimentPlot, load_panasonic_plot,
    load_trc_plot, load_uploaded_plot, panasonic_catalog_rows, trc_catalog_rows,
    uploaded_plot_options,
)


ROOT = Path(__file__).resolve().parent
STATUS_COLORS = {
    "CONFORME": ("#e8f5e9", "#1b5e20"),
    "REVISAR": ("#fff8e1", "#8a5a00"),
    "FALTANTE": ("#ffebee", "#b71c1c"),
    "NO APTO": ("#ffebee", "#b71c1c"),
    "NO EVALUADO": ("#eceff1", "#455a64"),
}
st.set_page_config(page_title="Auditoría de datos de batería", page_icon="🔋", layout="wide")
st.markdown(
    """
    <style>
    /* Theme-independent contrast for the soft audit cards. */
    .audit-card,
    .audit-card .audit-title,
    .audit-card .audit-description,
    .audit-card .audit-evidence {
        color: #1f2937 !important;
        -webkit-text-fill-color: #1f2937 !important;
    }
    .audit-card .audit-status {
        color: var(--audit-status-color) !important;
        -webkit-text-fill-color: var(--audit-status-color) !important;
    }
    .audit-action,
    .audit-action strong {
        color: #374151 !important;
        -webkit-text-fill-color: #374151 !important;
    }
    .audit-action {
        background-color: #f3f4f6 !important;
        border-left: 3px solid #9ca3af;
        border-radius: 0.25rem;
        padding: 0.42rem 0.65rem;
        margin: -0.15rem 0 0.65rem 0;
        font-size: 0.875rem;
    }
    .audit-description {
        display: block;
        margin: 0.3rem 0 0.22rem 0;
    }
    .audit-evidence {
        display: block;
    }
    .audit-legend {
        display: flex;
        flex-wrap: wrap;
        gap: 0.4rem 0.8rem;
        margin: 0.35rem 0 1rem 0;
    }
    .audit-legend-item {
        color: #1f2937 !important;
        -webkit-text-fill-color: #1f2937 !important;
        padding: 0.3rem 0.55rem;
        border-radius: 0.3rem;
        font-size: 0.86rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


def _result() -> object | None:
    return st.session_state.get("battery_audit_result")


def _frame(rows: object, columns: tuple[str, ...] | None = None) -> pd.DataFrame:
    frame = pd.DataFrame(rows)
    if columns:
        frame = frame[[name for name in columns if name in frame.columns]]
    for column in frame.select_dtypes(include="object").columns:
        frame[column] = frame[column].map(lambda value: "" if value is None else str(value))
    return frame


@st.cache_data(show_spinner=False)
def _panasonic_plot(experiment: str, variable: str) -> ExperimentPlot:
    return load_panasonic_plot(experiment, variable)


@st.cache_data(show_spinner=False)
def _trc_plot(condition: str, test: str, setup: str, variable: str) -> ExperimentPlot:
    return load_trc_plot(condition, test, setup, variable)


def _render_plot(plot: ExperimentPlot) -> None:
    figure, axes = plt.subplots(len(plot.panels), 1, figsize=(10.2, 4.3 * len(plot.panels)), squeeze=False, constrained_layout=True)
    color_map = {
        "Datos experimentales": "#202124", "Modelo": "#2166ac",
        "Superficie central": "#b2182b", "Superficie cerca del extremo": "#ef8a62",
        "Aire de la cámara": "#4d4d4d",
    }
    for axis, panel in zip(axes[:, 0], plot.panels, strict=True):
        for index, series in enumerate(panel.series):
            color = color_map.get(series.label, f"C{index}")
            if series.style == "points":
                step = max(1, len(panel.x) // 1200)
                axis.scatter(panel.x[::step], series.values[::step], s=8, color=color, alpha=0.7, label=series.label)
            else:
                axis.plot(panel.x, series.values, color=color, linewidth=1.7,
                          linestyle="--" if series.style == "dashed" else "-", label=series.label)
        axis.set_title(panel.title, loc="left", fontsize=11, fontweight="semibold")
        axis.set_xlabel(panel.x_label)
        axis.set_ylabel(panel.y_label)
        axis.grid(alpha=0.18, linewidth=0.7)
        axis.legend(frameon=False, ncol=min(3, len(panel.series)))
    st.pyplot(figure, width="stretch")
    plt.close(figure)
    st.caption(f"Función experimental: {plot.role}")
    if plot.note:
        st.info(plot.note)


def _show_panasonic_explorer(key_prefix: str) -> None:
    experiment = st.selectbox("Experimento", tuple(PANASONIC_CATALOG), key=f"{key_prefix}_experiment")
    variable = st.selectbox("Variable", tuple(PANASONIC_CATALOG[experiment]["variables"]), key=f"{key_prefix}_variable")
    try:
        with st.spinner("Cargando comparación..."):
            _render_plot(_panasonic_plot(experiment, variable))
    except Exception as exc:
        st.error(f"No fue posible cargar esta comparación: {type(exc).__name__}: {exc}")


def _show_model_reference() -> None:
    st.header("Modelo de referencia")
    st.write("Modelo Panasonic NCR18650PF y ensayos experimentales asociados.")
    try:
        payload = reference_model_payload()
    except Exception as exc:
        st.error(f"No fue posible leer la configuración del modelo: {exc}")
        return
    left, right = st.columns((1.05, 1.15), gap="large")
    with left:
        st.subheader("Parámetros del modelo")
        st.dataframe(
            _frame(payload["parameters"]),
            hide_index=True,
            width="stretch",
            height=520,
            column_config={
                "Parámetro": st.column_config.TextColumn(width="small"),
                "Descripción": st.column_config.TextColumn(width="large"),
                "Valor": st.column_config.NumberColumn(format="%.9g", width="medium"),
                "Unidad": st.column_config.TextColumn(width="small"),
            },
        )
    with right:
        st.subheader("Explorador del modelo")
        _show_panasonic_explorer("model")


def _temperature_roles(index: int, columns: tuple[str, ...]) -> dict[str, str]:
    choices = {
        "Superficie central": "cell surface middle",
        "Superficie cerca del extremo": "cell surface end",
        "Aire de cámara o ambiente": "chamber/ambient",
        "Ubicación no documentada": "temperature location unknown",
    }
    roles: dict[str, str] = {}
    if columns:
        st.markdown("**Ubicación de las mediciones de temperatura**")
    for number, name in enumerate(columns, start=1):
        label = st.selectbox(f"Sensor de temperatura {number}", tuple(choices), index=3, key=f"temp_role_{index}_{name}")
        roles[name] = choices[label]
    return roles


def _uploaded_data_controls() -> None:
    uploads = st.file_uploader("Seleccione uno o varios archivos CSV", type=("csv",), accept_multiple_files=True)
    if not uploads:
        st.info("Seleccione los registros que desea evaluar.")
        return
    chemistry = st.text_input("Química de la celda", value="desconocida")
    capacity_known = st.checkbox("Se conoce la capacidad nominal")
    capacity = st.number_input("Capacidad nominal (Ah)", min_value=0.001, value=3.0, step=0.1, disabled=not capacity_known)
    limits = st.text_input("Límites de operación documentados", value="desconocidos")
    limit_columns = st.columns(2)
    with limit_columns[0]:
        discharge_limit_known = st.checkbox("Se conoce el límite de descarga")
        discharge_cutoff = st.number_input("Voltaje final de descarga (V)", min_value=0.1, value=2.5, step=0.05, disabled=not discharge_limit_known)
    with limit_columns[1]:
        charge_limit_known = st.checkbox("Se conoce el límite de carga")
        charge_cutoff = st.number_input("Voltaje final de carga (V)", min_value=0.1, value=4.2, step=0.05, disabled=not charge_limit_known)
    role_options = (ROLE_CHARACTERIZATION, ROLE_LOW_RATE, ROLE_DYNAMIC, ROLE_THERMAL, ROLE_VALIDATION, ROLE_UNKNOWN)
    specs: list[UploadedFileSpec] = []
    for index, upload in enumerate(uploads):
        content = upload.getvalue()
        try:
            preview = inspect_csv_bytes(upload.name, content)
        except Exception as exc:
            preview = {"header_detected": False, "column_count": 0, "sample_count_estimate": 0, "temperature_columns": (), "preview_error": str(exc)}
        sampling_interval: float | None = None
        interval_constant_confirmed = False
        record_continuous_confirmed = False
        with st.expander(f"Registro {index + 1}", expanded=index == 0):
            first, second, third = st.columns(3)
            first.metric("Muestras", int(preview["sample_count_estimate"]))
            second.metric("Variables", int(preview.get("column_count", 0)))
            third.metric("Sensores térmicos", len(preview["temperature_columns"]))
            if "preview_error" in preview:
                st.warning(str(preview["preview_error"]))
            dataset_id = st.text_input("Identificador", value=f"registro_{index + 1}", key=f"dataset_{index}")
            condition = st.text_input("Condición experimental", value="sin especificar", key=f"condition_{index}")
            cell_id = st.text_input("Identificación de la celda", value="sin especificar", key=f"cell_{index}")
            setup = st.selectbox("Montaje", (SETUP_UNKNOWN, SETUP_REFERENCE, SETUP_OTHER), key=f"setup_{index}")
            roles = tuple(st.multiselect("Función experimental", role_options, default=[ROLE_UNKNOWN], key=f"roles_{index}"))
            temperature_roles = _temperature_roles(index, tuple(preview["temperature_columns"]))
            if "Time" not in tuple(preview.get("columns", ())):
                st.markdown("**Información temporal**")
                sampling_known = st.checkbox("Conozco el intervalo de muestreo", key=f"sampling_known_{index}")
                interval_value = st.number_input(
                    "Intervalo de muestreo (s)", min_value=0.001, value=1.0, step=0.1,
                    disabled=not sampling_known, key=f"sampling_interval_{index}",
                )
                interval_constant_confirmed = st.checkbox(
                    "Confirmo que el intervalo fue constante", disabled=not sampling_known,
                    key=f"sampling_constant_{index}",
                )
                record_continuous_confirmed = st.checkbox(
                    "Confirmo que el registro fue continuo", disabled=not sampling_known,
                    key=f"sampling_continuous_{index}",
                )
                sampling_interval = float(interval_value) if sampling_known else None
            fallback = None
            if not bool(preview["header_detected"]) and st.checkbox("El registro usa el formato estándar del equipo", key=f"fallback_{index}"):
                fallback = STANDARD_CHROMA_COLUMNS
            specs.append(
                UploadedFileSpec(
                    filename=upload.name,
                    content=content,
                    dataset_id=dataset_id,
                    cell_id=cell_id,
                    condition=condition,
                    setup_role=setup,
                    experimental_roles=roles or (ROLE_UNKNOWN,),
                    temperature_roles=temperature_roles,
                    fallback_columns=fallback,
                    sampling_interval_s=sampling_interval,
                    sampling_interval_constant_confirmed=interval_constant_confirmed,
                    record_continuous_confirmed=record_continuous_confirmed,
                )
            )
    if st.button("Cargar y auditar registros", type="primary"):
        context = CellContext(
            chemistry,
            float(capacity) if capacity_known else None,
            limits,
            float(discharge_cutoff) if discharge_limit_known else None,
            float(charge_cutoff) if charge_limit_known else None,
        )
        try:
            result = audit_uploaded_files(specs, context)
            st.session_state.update(battery_audit_result=result, battery_data_mode="propios", uploaded_specs=tuple(specs))
            st.success("Registros cargados y auditados.")
        except Exception as exc:
            st.error(f"No fue posible completar la auditoría: {type(exc).__name__}: {exc}")


def _show_data() -> None:
    st.header("Datos experimentales")
    mode = st.radio("Conjunto de datos", ("Datos experimentales propios", "Caso TRC18650", "Datos de referencia Panasonic"), horizontal=True)
    if mode == "Datos experimentales propios":
        _uploaded_data_controls()
    elif mode == "Caso TRC18650":
        st.dataframe(_frame(trc_catalog_rows()), hide_index=True, width="stretch")
        if st.button("Cargar caso TRC18650", type="primary"):
            try:
                result = build_trc_case_study(ROOT)
                st.session_state.update(battery_audit_result=result, battery_data_mode="trc")
                st.success("Caso cargado y auditado.")
            except Exception as exc:
                st.error(f"No fue posible cargar el caso: {type(exc).__name__}: {exc}")
    else:
        st.dataframe(_frame(panasonic_catalog_rows()), hide_index=True, width="stretch")
        if st.button("Cargar datos Panasonic", type="primary"):
            try:
                with st.spinner("Verificando los ensayos de referencia..."):
                    result = build_panasonic_reference_audit()
                st.session_state.update(battery_audit_result=result, battery_data_mode="panasonic")
                st.success("Datos Panasonic cargados y auditados.")
            except Exception as exc:
                st.error(f"No fue posible cargar los datos: {type(exc).__name__}: {exc}")


def _show_trc_explorer() -> None:
    conditions = tuple(dict.fromkeys(item[0] for item in TRC_OPTIONS))
    condition = st.selectbox("Condición experimental", conditions)
    tests = tuple(dict.fromkeys(item[1] for item in TRC_OPTIONS if item[0] == condition))
    test = st.selectbox("Tipo de ensayo", tests)
    setups = tuple(dict.fromkeys(item[2] for item in TRC_OPTIONS if item[:2] == (condition, test)))
    setup = st.selectbox("Montaje", setups)
    variables = next(item[3] for item in TRC_OPTIONS if item[:3] == (condition, test, setup))
    variable = st.selectbox("Variable", variables)
    try:
        _render_plot(_trc_plot(condition, test, setup, variable))
    except Exception as exc:
        st.error(f"No fue posible cargar la gráfica: {type(exc).__name__}: {exc}")


def _show_uploaded_explorer() -> None:
    specs = tuple(st.session_state.get("uploaded_specs", ()))
    if not specs:
        st.info("Cargue y audite sus registros antes de visualizarlos.")
        return
    labels = [f"{spec.condition} · Registro {index + 1}" for index, spec in enumerate(specs)]
    selected = st.selectbox("Condición experimental", range(len(specs)), format_func=lambda index: labels[index])
    spec = specs[selected]
    try:
        tests, variables = uploaded_plot_options(spec)
        test = st.selectbox("Tipo de ensayo", tests)
        variable = st.selectbox("Variable", variables)
        _render_plot(load_uploaded_plot(spec, test, variable))
    except Exception as exc:
        st.error(f"No fue posible visualizar este registro: {type(exc).__name__}: {exc}")


def _show_experiments() -> None:
    st.header("Visualización de ensayos")
    mode = st.session_state.get("battery_data_mode")
    if mode == "trc":
        _show_trc_explorer()
    elif mode == "panasonic":
        _show_panasonic_explorer("data")
    elif mode == "propios":
        _show_uploaded_explorer()
    else:
        st.info("Seleccione y cargue un conjunto en Datos experimentales.")


def _audit_card(row: dict[str, object]) -> None:
    status = str(row["estado"])
    background, foreground = STATUS_COLORS.get(status, ("#f5f5f5", "#333333"))
    code = html.escape(str(row["código"]))
    requirement = html.escape(str(row["requisito evaluado"]))
    description = html.escape(str(row["qué verifica"]))
    evidence = html.escape(str(row["evidencia"]))
    safe_status = html.escape(status)
    st.markdown(
        f"<div class='audit-card' style='--audit-status-color:{foreground};border-left:5px solid {foreground};"
        f"background-color:{background};padding:0.8rem 1rem;border-radius:0.35rem;margin-bottom:0.45rem'>"
        f"<b class='audit-title'>{code} · {requirement}</b>"
        f"<span class='audit-status' style='float:right;font-weight:700'>{safe_status}</span><br>"
        f"<span class='audit-description'>{description}</span>"
        f"<span class='audit-evidence'><strong>Resultado:</strong> {evidence}</span></div>", unsafe_allow_html=True)
    action = str(row.get("acción requerida", "")).strip()
    if status != "CONFORME" and action:
        st.markdown(
            f"<div class='audit-action'><strong>Acción requerida:</strong> {html.escape(action)}</div>",
            unsafe_allow_html=True,
        )


def _show_audit_guide() -> None:
    with st.expander("Criterios de evaluación"):
        for criterion in audit_criteria_payload():
            st.markdown(f"**{criterion['código']} · {criterion['nombre']}**")
            st.write(criterion["qué verifica"])
    legend = (
        ("CONFORME", "el requisito se cumple"),
        ("REVISAR", "se detectó una condición que requiere comprobación"),
        ("FALTANTE", "no se dispone de la información o del ensayo requerido"),
        ("NO APTO", "los datos existen, pero no son adecuados para el uso evaluado"),
        ("NO EVALUADO", "no se cumplen las condiciones para realizar una comprobación adicional"),
    )
    items = []
    for status, explanation in legend:
        background, _ = STATUS_COLORS[status]
        items.append(
            f"<span class='audit-legend-item' style='background-color:{background}'><strong>{status}:</strong> {explanation}</span>"
        )
    st.markdown(f"<div class='audit-legend'>{''.join(items)}</div>", unsafe_allow_html=True)


def _show_audit() -> None:
    st.header("Auditoría")
    result = _result()
    if result is None:
        st.info("Cargue un conjunto de datos para ejecutar las verificaciones.")
        return
    _show_audit_guide()
    counts = pd.Series([row["estado"] for row in result.checks]).value_counts()
    metric_columns = st.columns(5)
    for column, status in zip(metric_columns, ("CONFORME", "REVISAR", "FALTANTE", "NO APTO", "NO EVALUADO"), strict=True):
        column.metric(status, int(counts.get(status, 0)))
    for row in result.checks:
        _audit_card(dict(row))
    pending = [row for row in result.checks if row["estado"] not in {"CONFORME", "NO EVALUADO"}]
    next_action = next((str(row["acción requerida"]) for row in pending if str(row.get("acción requerida", "")).strip()), "No hay acciones pendientes.")
    st.subheader("Resumen")
    first, second = st.columns(2)
    with first:
        st.markdown("**Datos aprovechables**")
        for row in result.file_summaries:
            name = row.get("nombre visible") or " · ".join(str(row.get(key, "")) for key in ("condición", "ensayo", "montaje") if row.get(key))
            use = row.get("uso recomendado") or row.get("función experimental") or row.get("variables disponibles", "")
            st.write(f"• {name}: {use}")
    with second:
        st.metric("Requisitos pendientes", len(pending))
        st.markdown("**Siguiente acción recomendada**")
        st.write(next_action)
    buttons = st.columns(3)
    buttons[0].download_button("Descargar informe técnico", result.report_markdown.encode("utf-8"), "auditoria_bateria_informe.md", "text/markdown", width="stretch")
    buttons[1].download_button("Descargar verificaciones", result.checks_csv(), "auditoria_bateria_verificaciones.csv", "text/csv", width="stretch")
    buttons[2].download_button("Descargar resumen", result.summaries_csv(), "auditoria_bateria_resumen.csv", "text/csv", width="stretch")


st.title("Evaluación de datos para modelos de batería")
tab_model, tab_data, tab_plots, tab_audit = st.tabs(("Modelo de referencia", "Datos experimentales", "Visualización de ensayos", "Auditoría"))
with tab_model:
    _show_model_reference()
with tab_data:
    _show_data()
with tab_plots:
    _show_experiments()
with tab_audit:
    _show_audit()
