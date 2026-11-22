from __future__ import annotations

from pathlib import Path

import pytest

from runtime_paths import (
    PANASONIC_1C_PATH,
    PANASONIC_C20_PATH,
    PANASONIC_CCCV_PATH,
    PANASONIC_HPPC_PATH,
    PANASONIC_UDDS_PATH,
    TRC18650_DATA_DIRECTORY,
    require_data_file,
)


def test_all_bundled_runtime_datasets_are_available() -> None:
    panasonic = (
        PANASONIC_C20_PATH,
        PANASONIC_HPPC_PATH,
        PANASONIC_1C_PATH,
        PANASONIC_CCCV_PATH,
        PANASONIC_UDDS_PATH,
    )
    assert all(path.is_file() for path in panasonic)
    assert len(tuple(TRC18650_DATA_DIRECTORY.glob("*.csv"))) == 7


def test_missing_dataset_error_identifies_repository_filename(tmp_path: Path) -> None:
    missing = tmp_path / "required_measurement.csv"
    with pytest.raises(FileNotFoundError) as error:
        require_data_file(missing)
    message = str(error.value)
    assert "required_measurement.csv" in message
    assert "C:\\Users" not in message

