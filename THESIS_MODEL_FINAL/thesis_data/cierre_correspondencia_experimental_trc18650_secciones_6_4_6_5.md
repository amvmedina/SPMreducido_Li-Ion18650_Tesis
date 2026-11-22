# Cierre de correspondencia experimental TRC18650 para las secciones 6.4 y 6.5

## 1. Alcance y fuentes autoritativas

Los siete CSV Chroma constituyen el **inventario original auditado**, pero no representan siete experimentos finales independientes en la presentación de resultados. Las figuras de la sección 6.4 emplean una selección y unas correspondencias aceptadas, documentadas en:

- `thesis_data/uniandes_corrected_discharge_mapping_audit.csv` para las curvas de descarga;
- `thesis_data/uniandes_charge_temperature_correspondence.csv` y `thesis_data/uniandes_charge_audit.csv` para carga CC–CV;
- `thesis_data/uniandes_temperature_ranges.csv` para las señales térmicas;
- `thesis_data/trc18650_automatic_audit_experiments.csv` para `R_app,desc` y la clasificación de calidad.

La corrección aceptada de descarga reasigna exclusivamente `Voltage(V)` y la capacidad descargada calculada desde `Capacity(Ah)`. No transfiere automáticamente corriente, carga, capacidad total ni temperatura entre archivos.

## 2. Correspondencia definitiva de los resultados presentados

### 2.1 Curvas eléctricas de descarga

| Curva presentada | Configuración | Archivo original | Variables utilizadas | Observación de correspondencia |
|---|---|---|---|---|
| Rionegro · Aire 22–24 °C | Montaje de referencia, celda 1 | `amvTRC18650LiIonCycle1C Tambient cell1.csv` (`cycle_cell1_ambient`) | `Voltage(V)` y `Capacity(Ah)` del tramo de descarga | Correspondencia en el mismo archivo. |
| Rionegro · Aire 22–24 °C | Montaje con resistencia externa, celda 2 | `amvTRC18650LiIonCycle1C cell2 07232026.csv` (`cycle_cell2_undeclared`) | `Voltage(V)` y `Capacity(Ah)` del tramo de descarga | Temp7 es un canal de aire medido y compartido con el ensayo simultáneo de celda 1; no existe alineación muestra a muestra entre archivos. |
| Cali · Aire 29–31 °C | Montaje con resistencia externa, celda 2 | `amvTRC18650LiIonCycle1C T30C cell2 07242026.csv` (`cycle_cell2_30c`) | `Voltage(V)` y `Capacity(Ah)` del tramo de descarga | La temperatura de aire no se combina muestra a muestra con este archivo. |
| Valledupar · Aire 33–37 °C | Montaje con resistencia externa, curva eléctrica aceptada | `amvTRC18650LiIonCycle1C T30C cell1.csv` (`cycle_cell1_30c`) | Solo `Voltage(V)` y `Capacity(Ah)` del tramo de descarga | La temperatura 33–37 °C permanece en `cycle_cell1_35c`. Los dos registros no se combinan muestra a muestra y no se reasigna la corriente. |

La curva eléctrica de descarga de `cycle_cell1_35c` está excluida de estas comparaciones, aunque ese archivo conserva su papel como registro térmico y de carga para la condición 33–37 °C.

El eje de las figuras es `Capacidad descargada (Ah)`, calculada como `abs(Capacity(Ah)_i - Capacity(Ah)_inicio)` dentro del tramo registrado. No es SOC, no se normaliza entre 0 y 1 y no se extrapola.

### 2.2 Curvas eléctricas de carga CC–CV

La corrección especial de descarga **no se aplica a la carga**. Cada curva de carga usa voltaje, corriente, capacidad, paso y modo del mismo CSV:

| Curva presentada | Configuración | Archivo original | Variables utilizadas | Estado |
|---|---|---|---|---|
| Montaje de referencia, Rionegro · Aire 22–24 °C | Celda 1 | `amvTRC18650LiIonCycle1C Tambient cell1.csv` (`cycle_cell1_ambient`) | `Voltage(V)`, `Current(A)`, `Capacity(Ah)`, `Step No.`, `Mode` | Correspondencia confirmada en el mismo archivo. |
| Montaje con resistencia externa, Rionegro · Aire 22–24 °C | Celda 2 | `amvTRC18650LiIonCycle1C cell2 07232026.csv` (`cycle_cell2_undeclared`) | Las mismas cinco columnas | Correspondencia confirmada en el mismo archivo. |
| Montaje con resistencia externa, Cali · Aire 29–31 °C | Celda 2 | `amvTRC18650LiIonCycle1C T30C cell2 07242026.csv` (`cycle_cell2_30c`) | Las mismas cinco columnas | Correspondencia confirmada en el mismo archivo. |
| Montaje con resistencia externa, Valledupar · Aire 33–37 °C | Celda 1 | `amvTRC18650LiIonCycle1C T35C cell1.csv` (`cycle_cell1_35c`) | Las mismas cinco columnas | Correspondencia confirmada en el mismo archivo. |

La capacidad cargada se calcula desde el inicio del tramo de carga con `abs(Capacity(Ah)_i - Capacity(Ah)_inicio)`. Las etapas CC y CV se identifican con los pasos y las señales del mismo registro; no se reconstruye una carga con columnas de archivos distintos.

### 2.3 Resultados térmicos presentados

Las figuras térmicas de descarga y carga usan exclusivamente los tres registros de celda 1 con ubicaciones documentadas:

| Condición presentada | Archivo | Señales |
|---|---|---|
| Rionegro · Aire 22–24 °C | `cycle_cell1_ambient` | `Temp1` superficie central, `Temp2` cerca del extremo y `Temp7` aire de la cámara |
| Cali · Aire 29–31 °C | `cycle_cell1_30c` | `Temp1`, `Temp2`, `Temp7` |
| Valledupar · Aire 33–37 °C | `cycle_cell1_35c` | `Temp1`, `Temp2`, `Temp7` |

`Temp3` y `Temp4` de la celda 2 no se presentan como ubicaciones superficiales específicas porque esas ubicaciones no están documentadas. El Temp7 compartido de 22–24 °C no se copia ni se alinea muestra a muestra con el archivo de celda 2.

Las ciudades son referencias climáticas asociadas con los intervalos experimentales, no lugares declarados de realización de los ensayos.

## 3. Comparaciones entre montajes ya implementadas

### 3.1 Descarga

La comparación solicitada ya está implementada por `scripts/build_uniandes_context_figures.py` y usa:

- montaje de referencia: `cycle_cell1_ambient`, celda 1;
- montaje con resistencia externa: `cycle_cell2_undeclared`, celda 2;
- condición compartida: Rionegro · Aire 22–24 °C.

Archivos finales para tesis:

- `thesis_figures/final_es/12_2_comparacion_montajes_contacto_uniandes.pdf`;
- `thesis_figures/final_es/12_2_comparacion_montajes_contacto_uniandes.png`.

Archivos fuente equivalentes:

- `thesis_figures/final/uniandes_contact_arrangements_comparison_final.pdf`;
- `thesis_figures/final/uniandes_contact_arrangements_comparison_final.png`.

### 3.2 Carga CC–CV

La comparación solicitada ya está implementada por `scripts/build_uniandes_charge_temperature_figures.py`. Usa los mismos dos archivos de Rionegro, pero toma en cada caso `Voltage(V)`, `Current(A)` y `Capacity(Ah)` del tramo de carga del propio archivo.

Archivos finales para tesis:

- `thesis_figures/final_es/17_comparacion_carga_cccv_montajes_uniandes.pdf`;
- `thesis_figures/final_es/17_comparacion_carga_cccv_montajes_uniandes.png`.

Archivos fuente equivalentes:

- `thesis_figures/final/uniandes_cccv_contact_comparison_final.pdf`;
- `thesis_figures/final/uniandes_cccv_contact_comparison_final.png`.

La inspección de las salidas confirma títulos, ejes y leyendas correctos. No es necesario regenerar estas dos figuras.

Ambas comparaciones usan celdas diferentes. No son ensayos antes/después sobre una misma celda y no permiten aislar cuantitativamente una resistencia de contacto ni atribuir toda la diferencia de voltaje a una única causa.

## 4. Uso defendible de `R_app,desc` en la sección 6.5

Los valores direccionales verificados permanecen ligados a las transiciones de voltaje y corriente de sus archivos originales:

| Registro original | `R_app,desc` (Ω) | N descarga | Razón | Uso en la discusión |
|---|---:|---:|---:|---|
| `cycle_cell1_ambient` | 0.036383074 | 1 | 1.0000 | Referencia relativa del montaje de referencia a 22–24 °C. |
| `cycle_cell2_undeclared` | 0.185261851 | 1 | 5.0920 | Respuesta aparente del archivo usado como montaje con resistencia externa a 22–24 °C. |
| `cycle_cell2_30c` | 0.226546596 | 1 | 6.2267 | Respuesta aparente del archivo de celda 2 usado en las curvas eléctricas de Cali. |
| `cycle_cell1_30c` | 0.240607211 | 1 | 6.6132 | Diagnóstico del archivo original. No debe relabelarse automáticamente como `R_app,desc` de Valledupar, porque la corrección aceptada no reasignó corriente. |
| `cycle_cell1_35c` | 0.237533623 | 1 | 6.5287 | Diagnóstico del archivo usado para carga y temperatura de Valledupar; su descarga eléctrica no aparece en la comparación corregida. |
| `pulse_cell1_part1` | 0.221348508 | 7 | 6.0838 | Evidencia repetida del grupo elevado; no es un experimento final de las figuras de carga/descarga. |
| `pulse_cell1_part2` | 0.228938929 | 3 | 6.2925 | Continuación parcial del ensayo de pulsos; no es un experimento final independiente. |

Para la comparación directa entre montajes a 22–24 °C pueden citarse `0.036383 Ω` y `0.185262 Ω`, aclarando que corresponden a celdas distintas y a una sola transición de descarga por ciclo. Los pulsos aportan repetibilidad al diagnóstico del grupo elevado, pero sus partes 1 y 2 no deben contarse como dos experimentos finales completos.

`R_app,desc` es una respuesta aparente del camino observado. No es la resistencia interna efectiva de la celda, una medición aislada de la resistencia externa del montaje ni el `R0` del modelo Panasonic. La respuesta puede incluir contribuciones de celda, polarización rápida, contactos, cableado, sensado, instrumentación y resolución de muestreo.

## 5. Recomendaciones que deben conservarse en la sección 6.5

La discusión puede mantenerse breve y centrarse en estas conclusiones:

1. Verificar y documentar el montaje de fuerza y sensado, las conexiones de voltaje, los contactos, el cableado, la identidad de las celdas y la ubicación de los sensores antes de usar los registros para identificación.
2. Usar una referencia TRC de baja tasa y tiempo medido para adaptar capacidad y representación de voltaje/OCP; una descarga cargada a 1C no sustituye esa referencia.
3. Reservar para identificación dinámica solo pulsos limpios, temporizados y obtenidos con el montaje verificado; los pulsos actuales son diagnósticos.
4. Conservar separadas calibración y validación, reservando un perfil independiente.
5. Tratar los resultados térmicos disponibles como descriptivos mientras no exista una correspondencia temporal y una partición del calor defendibles.

La ausencia de tiempo medido es una limitación para identificación dinámica y térmica, pero no debe desplazar el resultado central: los registros seleccionados permiten documentar la influencia del montaje sobre el voltaje observado, sin identificar de manera única la contribución física responsable.

## 6. Advertencias y salidas que no deben reutilizarse sin corrección

No requieren regeneración las dos figuras comparativas identificadas en la sección 3 ni las figuras eléctricas y térmicas generadas por los scripts de correspondencia aceptados.

Sin embargo, las siguientes salidas anteriores conservan la mediana conjunta carga/descarga y no deben usarse para citar los valores definitivos de `R_app,desc`:

- `thesis_figures/final/uniandes_apparent_series_path_final.*`;
- `thesis_figures/final_es/11_respuesta_serie_aparente_uniandes.*`;
- `thesis_data/battery_dataset_assessment_summary.csv`;
- el párrafo Uniandes de `thesis_data/final_results_narrative.md`;
- la sección 3.1 de `thesis_data/diagnostico_auditoria_transferencia_trc18650.md`.

Si la figura 11 de resistencia aparente se mantiene en la tesis, debe regenerarse desde `trc18650_automatic_audit_experiments.csv`, usando `R_app,desc`, y debe titular los siete elementos como **registros originales auditados**, no como siete experimentos finales. Si no es indispensable, la opción más clara es omitirla y citar la tabla compacta de esta sección.

Para redactar la sección 6.5 deben considerarse autoritativos:

- `thesis_data/verificacion_resistencia_aparente_trc18650.md`;
- `thesis_data/trc18650_automatic_audit_report.md`;
- `thesis_data/trc18650_automatic_audit_experiments.csv`;
- el presente informe de correspondencia.

No se modificaron datos originales, cálculos de `R_app,desc`, parámetros Panasonic, curvas aceptadas ni correspondencias experimentales.
