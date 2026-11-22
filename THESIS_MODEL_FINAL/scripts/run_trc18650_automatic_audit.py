"""Run the first functional automatic audit for the seven TRC18650 records."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.trc18650_automatic_audit import (  # noqa: E402
    PanasonicModelSnapshot,
    ReferenceDataset,
    run_audit,
    write_csv,
    write_engineer_report,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from models.panasonic_thermal_parameters import (  # noqa: E402
    C_TH_J_PER_K,
    G_TH_W_PER_K,
)
from scripts.audit_own_trc18650_data import load_all  # noqa: E402
from scripts.validate_panasonic_1c_forward import (  # noqa: E402
    NOMINAL_PARAMETERS,
    Q_N_C,
    Q_P_C,
)
from runtime_paths import (  # noqa: E402
    PANASONIC_1C_PATH,
    PANASONIC_C20_PATH,
    PANASONIC_CCCV_PATH,
    PANASONIC_HPPC_PATH,
    PANASONIC_UDDS_PATH,
)


DEFAULT_OUTPUT_DIRECTORY = ROOT / "thesis_data"
REPORT_NAME = "trc18650_automatic_audit_report.md"
EXPERIMENTS_NAME = "trc18650_automatic_audit_experiments.csv"
RECOMMENDATIONS_NAME = "trc18650_automatic_audit_recommendations.csv"


def reference_datasets() -> tuple[ReferenceDataset, ...]:
    return (
        ReferenceDataset(
            "C/20 capacity and OCP reference",
            PANASONIC_C20_PATH,
        ),
        ReferenceDataset(
            "HPPC dynamic calibration",
            PANASONIC_HPPC_PATH,
        ),
        ReferenceDataset(
            "1C electrical validation and thermal calibration",
            PANASONIC_1C_PATH,
        ),
        ReferenceDataset(
            "UDDS independent validation",
            PANASONIC_UDDS_PATH,
        ),
        ReferenceDataset(
            "CC-CV charge validation",
            PANASONIC_CCCV_PATH,
        ),
    )


def model_snapshot() -> PanasonicModelSnapshot:
    return PanasonicModelSnapshot(
        q_cell_ah=Q_CELL_AH,
        q_n_c=Q_N_C,
        q_p_c=Q_P_C,
        x_n_start=X_N_START,
        x_n_end=X_N_END,
        x_p_start=X_P_START,
        x_p_end=X_P_END,
        r0_ohm=NOMINAL_PARAMETERS["R0"],
        alpha_n_s=NOMINAL_PARAMETERS["alpha_n"],
        alpha_p_s=NOMINAL_PARAMETERS["alpha_p"],
        d_n_per_s=NOMINAL_PARAMETERS["d_n"],
        d_p_per_s=NOMINAL_PARAMETERS["d_p"],
        c_th_j_per_k=C_TH_J_PER_K,
        g_th_w_per_k=G_TH_W_PER_K,
        source_files=(
            ROOT / "models" / "ocp_panasonic_qi_frozen.py",
            ROOT / "scripts" / "validate_panasonic_1c_forward.py",
            ROOT / "models" / "panasonic_thermal_parameters.py",
        ),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the evidence-gated TRC18650 automatic audit."
    )
    parser.add_argument(
        "--output-directory",
        type=Path,
        default=DEFAULT_OUTPUT_DIRECTORY,
        help="Directory for the new report and CSV files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    experiments = load_all()
    bundle = run_audit(experiments, model_snapshot(), reference_datasets())
    output = args.output_directory.resolve()
    report_path = output / REPORT_NAME
    experiments_path = output / EXPERIMENTS_NAME
    recommendations_path = output / RECOMMENDATIONS_NAME
    write_engineer_report(report_path, bundle)
    write_csv(experiments_path, bundle.classifications)
    write_csv(recommendations_path, bundle.recommendations)
    print(f"Audited {len(bundle.classifications)} TRC18650 datasets.")
    print(f"Elevated apparent-response group: {sum(row['nonideality_flag'].startswith('EVIDENCE') for row in bundle.classifications)} datasets.")
    print(f"Transfer gate: {bundle.transfer_decision.summary}")
    print(f"Report: {report_path}")
    print(f"Experiment table: {experiments_path}")
    print(f"Recommendations: {recommendations_path}")


if __name__ == "__main__":
    main()
