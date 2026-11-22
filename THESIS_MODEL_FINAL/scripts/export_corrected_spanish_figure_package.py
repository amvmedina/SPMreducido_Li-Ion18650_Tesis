"""Regenerate and export the corrected final Spanish thesis figure package.

This script performs frozen forward evaluations only. It does not fit, optimize,
or modify any accepted electrical or thermal parameter. Metrics displayed by the
four corrected figures are read from the authoritative numerical CSV outputs;
the thermal builders additionally verify their plotted arrays against them.
"""

from __future__ import annotations

import csv
import hashlib
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import build_final_thesis_package as package  # noqa: E402


FINAL_ES = ROOT / "thesis_figures" / "final_es"
LOCALIZED = ROOT / "figuras_finales_es"
LOCALIZED_THESIS = LOCALIZED / "tesis"
LOCALIZED_PRESENTATION = LOCALIZED / "presentacion"
AUDIT_CSV = ROOT / "thesis_data" / "final_external_model_results_audit_2026-11-12.csv"

FIGURES = (
    ("01_ajuste_ocp_panasonic", "panasonic_ocp_fit_final"),
    ("02_1_calibracion_hppc_validacion_bajo_soc", "hppc_calibracion_validacion_final"),
    ("02_2_definicion_fases_pulso_hppc", "hppc_fases_pulso_final"),
    ("03_validacion_electrica_1c_panasonic", "panasonic_1c_voltage_validation_final"),
    ("04_validacion_electrica_udds_panasonic", "panasonic_udds_voltage_validation_final"),
    ("05_transferencia_electrica_carga_cccv", "panasonic_cccv_voltage_transfer_final"),
    ("06_separacion_direccional_carga_descarga", "panasonic_directional_voltage_separation_final"),
    ("07_calibracion_termica_1c", "panasonic_1c_thermal_calibration_final"),
    ("08_transferencia_termica_udds", "panasonic_udds_thermal_transfer_final"),
    ("09_transferencia_termica_carga_cccv", "panasonic_cccv_thermal_transfer_final"),
    ("10_seleccion_modelo_termico", "thermal_model_selection_summary"),
    ("11_respuesta_serie_aparente_uniandes", "uniandes_apparent_series_path_final"),
    ("12_1_voltaje_contactos_defectuosos_uniandes", "uniandes_defective_contacts_discharge_final"),
    ("12_2_comparacion_montajes_contacto_uniandes", "uniandes_contact_arrangements_comparison_final"),
    ("13_transferencia_panasonic_a_trc", "panasonic_to_trc_transfer_final"),
    ("14_flujo_evaluacion_dataset_bateria", "battery_dataset_assessment_flow_final"),
    ("15_contexto_geografico_intervalos_uniandes", "uniandes_geographical_context_final"),
    ("16_carga_cccv_resistencia_externa_uniandes", "uniandes_cccv_external_conditions_final"),
    ("17_comparacion_carga_cccv_montajes_uniandes", "uniandes_cccv_contact_comparison_final"),
    ("18_temperaturas_descarga_uniandes", "uniandes_temperature_discharge_final"),
    ("19_temperaturas_carga_uniandes", "uniandes_temperature_charge_final"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_expected_metrics() -> None:
    udds_path = package.DATA / "panasonic_udds_thermal_validation_metrics.csv"
    cccv_path = package.DATA / "panasonic_cccv_thermal_validation_metrics.csv"
    expected = {
        "udds_b_rmse": 0.6746964046939614,
        "udds_b_bias": -0.6511176255833355,
        "cccv_b_rmse": 0.8643869967745215,
    }
    actual = {
        "udds_b_rmse": package.metric_value(
            udds_path, validation="validation_b", subset="overall", metric="rmse_k"
        ),
        "udds_b_bias": package.metric_value(
            udds_path,
            validation="validation_b",
            subset="overall",
            metric="mean_signed_error_k",
        ),
        "cccv_b_rmse": package.metric_value(
            cccv_path,
            scope="temperature_metrics",
            subset="validation_b",
            metric="rmse_k",
        ),
    }
    for name, value in actual.items():
        if value != expected[name]:
            raise ValueError(f"Authoritative metric {name} changed: {value!r}.")


def regenerate_corrected_sources() -> None:
    package.configure_style()
    package.THESIS_FIGURES.mkdir(parents=True, exist_ok=True)
    package.PRESENTATION.mkdir(parents=True, exist_ok=True)

    electrical = package.electrical_forward_arrays()
    package.figure_1c(electrical["one_c"])
    package.figure_udds(electrical["udds"])
    package.figure_cccv(electrical["cccv"])

    thermal = package.thermal_forward_arrays()
    package.figure_thermal_udds(thermal["udds"])
    package.figure_thermal_cccv(thermal["cccv"])

    _, _, experiments = package.own_case_arrays()
    package.figure_uniandes_discharge_comparisons(experiments)
    package.figure_uniandes_geographical_context()
    package.figure_uniandes_charge_temperature()


def export_package() -> list[dict[str, str]]:
    for directory in (FINAL_ES, LOCALIZED_THESIS, LOCALIZED_PRESENTATION):
        directory.mkdir(parents=True, exist_ok=True)

    manifest: list[dict[str, str]] = []
    for spanish_stem, source_stem in FIGURES:
        thesis_sources = {
            extension: package.THESIS_FIGURES / f"{source_stem}.{extension}"
            for extension in ("pdf", "png")
        }
        presentation_sources = {
            extension: package.PRESENTATION / f"{source_stem}.{extension}"
            for extension in ("pdf", "png")
        }
        for source in (*thesis_sources.values(), *presentation_sources.values()):
            if not source.is_file() or source.stat().st_size == 0:
                raise FileNotFoundError(source)

        for extension, source in thesis_sources.items():
            shutil.copy2(source, FINAL_ES / f"{spanish_stem}.{extension}")
            shutil.copy2(source, LOCALIZED_THESIS / f"{spanish_stem}.{extension}")
        for extension, source in presentation_sources.items():
            shutil.copy2(source, LOCALIZED_PRESENTATION / f"{spanish_stem}.{extension}")

        final_png = FINAL_ES / f"{spanish_stem}.png"
        manifest.append(
            {
                "spanish_stem": spanish_stem,
                "source_stem": source_stem,
                "thesis_pdf": f"thesis_figures/final_es/{spanish_stem}.pdf",
                "thesis_png": f"thesis_figures/final_es/{spanish_stem}.png",
                "sha256_thesis_png": sha256(final_png),
                "data_policy": "accepted frozen forward-only configuration; no fit",
            }
        )

    for destination in (
        FINAL_ES / "resultados_verificados_modelo_final.csv",
        LOCALIZED / "resultados_verificados_modelo_final.csv",
    ):
        shutil.copy2(AUDIT_CSV, destination)

    for destination in (
        FINAL_ES / "manifest_paquete_es_verificado.csv",
        LOCALIZED / "manifest_paquete_es_verificado.csv",
    ):
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
            writer.writeheader()
            writer.writerows(manifest)
    return manifest


def contact_sheet() -> None:
    columns = 3
    rows = (len(FIGURES) + columns - 1) // columns
    figure, axes = plt.subplots(rows, columns, figsize=(16, 4 * rows), constrained_layout=True)
    for axis, (spanish_stem, _) in zip(axes.flat, FIGURES, strict=False):
        axis.imshow(plt.imread(LOCALIZED_PRESENTATION / f"{spanish_stem}.png"))
        axis.set_title(spanish_stem, fontsize=10)
        axis.axis("off")
    for axis in axes.flat[len(FIGURES) :]:
        axis.axis("off")
    figure.suptitle("Paquete final de figuras en espa\u00f1ol", fontsize=22, weight="bold")
    figure.savefig(LOCALIZED / "contact_sheet_figuras.png", dpi=180, bbox_inches="tight")
    plt.close(figure)


def verify_exports(manifest: list[dict[str, str]]) -> None:
    if len(manifest) != len(FIGURES):
        raise RuntimeError("The exported Spanish package is incomplete.")
    for row in manifest:
        spanish_stem = row["spanish_stem"]
        for directory in (FINAL_ES, LOCALIZED_THESIS, LOCALIZED_PRESENTATION):
            for extension in ("pdf", "png"):
                path = directory / f"{spanish_stem}.{extension}"
                if not path.is_file() or path.stat().st_size == 0:
                    raise RuntimeError(f"Missing or empty export: {path}")


def main() -> None:
    if not AUDIT_CSV.is_file():
        raise FileNotFoundError(AUDIT_CSV)
    verify_expected_metrics()
    regenerate_corrected_sources()
    manifest = export_package()
    contact_sheet()
    verify_exports(manifest)
    print("Corrected Spanish figure package exported from accepted frozen results.")


if __name__ == "__main__":
    main()
