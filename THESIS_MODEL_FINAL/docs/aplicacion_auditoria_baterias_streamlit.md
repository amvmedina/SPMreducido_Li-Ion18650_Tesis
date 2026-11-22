# Aplicación de evaluación de datos de batería

## Alcance

La aplicación organiza el trabajo en cuatro secciones: **Modelo de referencia**, **Datos experimentales**, **Visualización de ensayos** y **Auditoría**. Reutiliza el modelo Panasonic y el motor determinístico de evaluación existentes; no ajusta parámetros ni altera mediciones.

Los cargadores convierten CSV y archivos Panasonic a una misma estructura experimental. E01–E08 se ejecutan después de esa conversión mediante `evaluate_common_audit`; no existen reglas de aprobación por marca, nombre de archivo o modalidad.

Admite tres modalidades:

- archivos CSV propios;
- caso TRC18650, presentado por condición térmica y tipo de ensayo;
- datos Panasonic de referencia: C/20, HPPC, descarga a 1C, carga CC–CV y UDDS.

El explorador muestra únicamente variables respaldadas por la correspondencia experimental. Los gráficos TRC usan capacidad cuando esa es la coordenada experimental aceptada y no presentan un tiempo reconstruido como tiempo medido.

Cuando un CSV no contiene tiempo, el formulario permite indicar el intervalo de muestreo y confirmar que fue constante y que el registro fue continuo. Solo con ambas confirmaciones la escala reconstruida se acepta como verificada. Una escala estimada a partir de corriente y capacidad permanece provisional y conduce a REVISAR cuando el tiempo es el único requisito pendiente.

## Instalación y ejecución

Desde la raíz del proyecto:

```powershell
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

## Uso

1. Consulte el modelo y seleccione un ensayo Panasonic en **Modelo de referencia**.
2. Seleccione una de las tres modalidades en **Datos experimentales** y cargue el conjunto.
3. Elija condición, ensayo y variable en **Visualización de ensayos**.
4. Revise los estados E01–E08 en **Auditoría** y descargue el informe o los CSV.

Los estados son `CONFORME`, `REVISAR`, `FALTANTE`, `NO APTO` y `NO EVALUADO`. Este último se reserva para comprobaciones adicionales que no pueden realizarse, como una comparación entre montajes cuando el conjunto no fue diseñado para ello; no sustituye un requisito obligatorio faltante. Las acciones se muestran solo para estados que requieren atención.

La sección desplegable **Criterios de evaluación** muestra únicamente el nombre y una descripción breve de cada verificación. Los criterios exactos se documentan en `docs/tabla_tecnica_verificaciones_auditoria.md`. La comparación eléctrica entre montajes usa una respuesta aparente de la transición reposo → descarga; no separa celda, montaje, contactos, cables, instrumentación o polarización rápida.

## Caso TRC18650

La vista presenta Rionegro (22–24 °C), Cali (29–31 °C) y Valledupar (33–37 °C), con descarga y carga CC–CV cuando existe correspondencia respaldada. Los archivos fuente conservan su procedencia interna, pero no se presentan como siete experimentos finales. Los pulsos excluidos no participan: si no se carga otro HPPC apto, la verificación correspondiente queda como FALTANTE.

## Pruebas

```powershell
python -m pytest tests/test_battery_audit_app_service.py tests/test_battery_audit_visualization.py tests/test_battery_audit_streamlit_app.py -q
```

La tabla `docs/tabla_tecnica_verificaciones_auditoria.md` documenta el criterio exacto de E01–E08. La inspección de ausencia de IA se registra en `docs/verificacion_sin_ia_auditoria.md`.
