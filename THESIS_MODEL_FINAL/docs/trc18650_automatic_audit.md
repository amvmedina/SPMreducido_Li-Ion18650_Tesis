# Ejecución de la auditoría automática TRC18650

## Objetivo

La primera versión funcional audita los siete CSV Chroma, verifica la disponibilidad de los archivos Panasonic aceptados y genera decisiones conservadoras de caracterización, diagnóstico, calibración, validación y transferencia.

El programa no modifica datos originales, no ejecuta optimización, no recalibra el modelo Panasonic y no corrige voltajes.

## Ejecución

Desde la raíz del repositorio:

```powershell
python scripts/run_trc18650_automatic_audit.py
```

Para escribir las tres salidas en otro directorio:

```powershell
python scripts/run_trc18650_automatic_audit.py --output-directory ruta/de/salida
```

## Salidas predeterminadas

- `thesis_data/trc18650_automatic_audit_report.md`
- `thesis_data/trc18650_automatic_audit_experiments.csv`
- `thesis_data/trc18650_automatic_audit_recommendations.csv`

## Entradas y trazabilidad

- Los siete archivos TRC y sus roles se obtienen mediante `scripts/audit_own_trc18650_data.py`.
- La carga, segmentación, capacidad y transiciones se calculan con `analysis/experimental_data_audit.py`.
- Las puertas de calidad, respuesta serie relativa y transferencia proceden de `analysis/battery_dataset_assessment.py`.
- Los parámetros Panasonic se leen de los módulos y scripts aceptados; no se copian a los archivos TRC.
- La presencia de los cinco MAT Panasonic de referencia se verifica antes de generar resultados.

La configuración actual conserva las rutas absolutas de los CSV originales. Si el conjunto se mueve, debe actualizarse únicamente la configuración de metadatos; nunca los archivos medidos.

## Pruebas

```powershell
python -m pytest tests/test_trc18650_automatic_audit.py -q
```

Las pruebas verifican que la referencia ambiente no se promueva a calibración, que los pulsos afectados permanezcan como diagnóstico, que la puerta de transferencia no active ajustes y que la recuperación de documentación preceda a nuevos experimentos.
