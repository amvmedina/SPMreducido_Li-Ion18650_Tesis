from __future__ import annotations

from pathlib import Path

import pytest


def _relative_luminance(hex_color: str) -> float:
    values = [int(hex_color[index:index + 2], 16) / 255.0 for index in (1, 3, 5)]
    linear = [value / 12.92 if value <= 0.04045 else ((value + 0.055) / 1.055) ** 2.4 for value in values]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def _contrast(first: str, second: str) -> float:
    bright, dark = sorted((_relative_luminance(first), _relative_luminance(second)), reverse=True)
    return (bright + 0.05) / (dark + 0.05)


@pytest.mark.parametrize("background", ("#e8f5e9", "#fff8e1", "#ffebee"))
def test_audit_card_palette_has_theme_independent_text_contrast(background: str) -> None:
    assert _contrast("#1f2937", background) >= 7.0


def test_streamlit_case_button_renders_calculated_results() -> None:
    pytest.importorskip("streamlit")
    source_directory = Path(r"C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data")
    if not source_directory.is_dir():
        pytest.skip("Configured TRC case-study files are not installed.")

    from streamlit.testing.v1 import AppTest

    app_path = Path(__file__).resolve().parents[1] / "app.py"
    app = AppTest.from_file(str(app_path), default_timeout=30).run()
    assert not app.exception
    assert [tab.label for tab in app.tabs] == [
        "Modelo de referencia",
        "Datos experimentales",
        "Visualización de ensayos",
        "Auditoría",
    ]

    app.radio[0].set_value("Caso TRC18650").run()
    case_button = next(button for button in app.button if button.label == "Cargar caso TRC18650")
    case_button.click().run()
    assert not app.exception
    rendered_markdown = "\n".join(str(element.value) for element in app.markdown)
    assert "class='audit-card'" in rendered_markdown
    assert "class='audit-title'" in rendered_markdown
    assert "class='audit-evidence'" in rendered_markdown
    assert "class='audit-action'" in rendered_markdown
    assert "Criterios de evaluación" in [expander.label for expander in app.expander]
    assert "NO EVALUADO:" in rendered_markdown
    assert "<strong>Resultado:</strong>" in rendered_markdown
    assert "color: #1f2937 !important" in rendered_markdown
    assert "-webkit-text-fill-color: #1f2937 !important" in rendered_markdown
    assert "cycle_cell1" not in rendered_markdown
    assert "pulse_cell1" not in rendered_markdown
    assert "Prototipo local" not in rendered_markdown
