from __future__ import annotations

import numpy as np
import pytest

from models.ocp_panasonic_qi_frozen import (
    Q_CELL_AH,
    SOURCE_MULTISTART_START_ID,
    UP_KNOT_STOICHIOMETRIES,
    UP_KNOT_VOLTAGES_V,
    U_n,
    U_p,
    V_eq,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
    electrode_stoichiometries,
)


def test_frozen_provenance_and_nominal_values() -> None:
    assert SOURCE_MULTISTART_START_ID == 3
    assert Q_CELL_AH == 2.994979138519
    assert X_N_START == 0.720451136448
    assert X_N_END == 0.251849704066
    assert X_P_START == 0.224027933455
    assert X_P_END == 0.953477823621


def test_frozen_arrays_are_read_only_and_monotonic() -> None:
    assert not UP_KNOT_STOICHIOMETRIES.flags.writeable
    assert not UP_KNOT_VOLTAGES_V.flags.writeable
    assert UP_KNOT_STOICHIOMETRIES.shape == (13,)
    assert UP_KNOT_VOLTAGES_V.shape == (13,)
    assert np.all(np.diff(UP_KNOT_VOLTAGES_V) <= 0.0)


def test_nominal_effective_trajectories_and_voltage_are_finite() -> None:
    z = np.linspace(0.0, 1.0, 1001)
    x_n, x_p = electrode_stoichiometries(z)
    assert x_n[0] == pytest.approx(X_N_START)
    assert x_n[-1] == pytest.approx(X_N_END)
    assert x_p[0] == pytest.approx(X_P_START)
    assert x_p[-1] == pytest.approx(X_P_END)
    np.testing.assert_allclose(V_eq(z), U_p(x_p) - U_n(x_n))
    assert np.all(np.isfinite(V_eq(z)))


@pytest.mark.parametrize("invalid_z", [-0.001, 1.001, np.nan, np.inf])
def test_nominal_trajectory_rejects_invalid_z_without_clipping(invalid_z: float) -> None:
    with pytest.raises(ValueError, match="finite|no clipping"):
        electrode_stoichiometries(invalid_z)
