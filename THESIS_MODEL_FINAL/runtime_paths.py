"""Portable paths for datasets required by the Streamlit application."""

from __future__ import annotations

from pathlib import Path


APP_ROOT = Path(__file__).resolve().parent
DATA_DIRECTORY = APP_ROOT / "data"
PANASONIC_DATA_DIRECTORY = DATA_DIRECTORY / "panasonic"
TRC18650_DATA_DIRECTORY = DATA_DIRECTORY / "trc18650"

PANASONIC_C20_PATH = PANASONIC_DATA_DIRECTORY / "c20_ocv_25C_converted.csv"
PANASONIC_HPPC_PATH = PANASONIC_DATA_DIRECTORY / "hppc_25C_converted.csv"
PANASONIC_1C_PATH = PANASONIC_DATA_DIRECTORY / "discharge_1C_25C_converted.csv"
PANASONIC_PRECEDING_1C_PATH = PANASONIC_1C_PATH
PANASONIC_CCCV_PATH = PANASONIC_DATA_DIRECTORY / "charge_3349_2_converted.csv"
PANASONIC_UDDS_PATH = PANASONIC_DATA_DIRECTORY / "udds_25C_converted.csv"


def require_data_file(path: Path) -> Path:
    """Return *path* or raise a deployment-oriented missing-file error."""

    if path.is_file():
        return path
    try:
        expected = path.relative_to(APP_ROOT).as_posix()
    except ValueError:
        expected = path.name
    raise FileNotFoundError(
        f"Required application dataset is missing: {path.name}. "
        f"Expected repository location: THESIS_MODEL_FINAL/{expected}"
    )

