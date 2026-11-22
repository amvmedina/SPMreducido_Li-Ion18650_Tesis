# Verificación de la resistencia aparente en los experimentos TRC18650

## 1. Resultado de la verificación

La diferencia entre `0.240607 Ω` y `0.232494 Ω` para `cycle_cell1_30c` no procede de una fórmula distinta, un cambio de los datos ni una versión diferente del extractor. Ambos valores se obtienen de las mismas dos transiciones guardadas en `own_data_apparent_resistance_transitions.csv`:

- reposo → descarga: `R_app = 0.240607210626186 Ω`;
- reposo → carga: `R_app = 0.224380834280597 Ω`.

El diagnóstico preliminar mostró la transición de **descarga**. La auditoría automática mostró la **mediana de descarga y carga**; al existir dos valores, esa mediana es su promedio:

`median(0.240607210626186, 0.224380834280597) = 0.232494022453391 Ω`.

Por tanto, ninguno de los dos números contiene un error aritmético. La inconsistencia consiste en haber usado agregaciones direccionales diferentes bajo una misma etiqueta general de `R_app`.

## 2. Trazabilidad del cálculo

La función `extract_rest_to_load_transitions()` de `analysis/experimental_data_audit.py` identifica el índice previo `i` y el primer índice cargado `i+1` mediante estas condiciones:

- `abs(I[i]) <= 0.02 A`;
- `abs(I[i+1]) >= 2.5 A`;
- `abs(I[i+1] - I[i]) >= 2.5 A`.

No se utiliza una ventana temporal ni un promedio de muestras. Para cada par adyacente calcula:

`Delta I = I[i+1] - I[i]`

`Delta V = V[i+1] - V[i]`

`R_app = abs(Delta V) / abs(Delta I)`.

En todas las transiciones conservadas la muestra previa tiene `I=0 A`. Para descarga, esta expresión coincide con `abs(V_reposo - V_carga)/abs(I_carga)`. Las 20 transiciones almacenadas satisfacen exactamente la fórmula, con diferencia numérica máxima igual a cero a la precisión guardada, y todas tienen el signo de respuesta de voltaje esperado.

`audit_experiment()` calcula y conserva dos resúmenes diferentes:

- `r_app_median_ohm`: mediana de **todas** las transiciones reposo → carga activa, sin distinguir descarga y carga;
- `discharge_r_app_median_ohm`: mediana de las transiciones cuya dirección es **descarga**.

La ruta de la auditoría automática es:

1. `evaluate_data_quality_gate()` llama al mismo extractor y a `audit_experiment()`;
2. `evaluate_relative_nonideality_gate()` selecciona actualmente `r_app_median_ohm` y `r_app_transition_count`;
3. `analysis/trc18650_automatic_audit.py` exporta esa mediana general como `r_app_median_ohm`.

El valor preliminar de cada ciclo coincide exactamente con `discharge_r_app_median_ohm` del archivo `own_data_quality_audit.csv`. No se encontró un generador independiente del informe `diagnostico_auditoria_transferencia_trc18650.md`; sus valores se pueden rastrear directamente a esa columna y a las transiciones individuales. El propio informe preliminar mezcla ambos criterios: la tabla por experimento usa descarga, mientras que su sección 3.1 usa la mediana general `0.035875537 Ω` como referencia y las razones generales `4.95–6.48`.

La regeneración en memoria desde los siete CSV originales produjo las mismas 20 claves de transición y diferencias máximas de `0` en `Delta V`, `Delta I` y `R_app` respecto al CSV guardado. Además, `own_data_quality_audit.csv`, creado antes de la auditoría automática, ya contiene simultáneamente ambas columnas. La evidencia descarta una diferencia causada por versiones de los datos o del extractor.

## 3. Comparación de los siete registros

| Registro | Valor del diagnóstico inicial [Ω] | Valor de la auditoría final [Ω] | Causa de la diferencia | Valor recomendado [Ω] |
|---|---:|---:|---|---:|
| `pulse_cell1_part1` | 0.221348508 | 0.221348508 | Las siete transiciones son de descarga; ambos criterios coinciden. | 0.221348508 |
| `pulse_cell1_part2` | 0.228938929 | 0.228938929 | Las tres transiciones son de descarga; ambos criterios coinciden. | 0.228938929 |
| `cycle_cell1_ambient` | 0.036383074 | 0.035875537 | El primero usa la única transición de descarga; el segundo es la mediana de una descarga y una carga. | 0.036383074 |
| `cycle_cell1_30c` | 0.240607211 | 0.232494022 | El primero usa la única transición de descarga; el segundo combina 0.240607211 Ω en descarga y 0.224380834 Ω en carga. | 0.240607211 |
| `cycle_cell1_35c` | 0.237533623 | 0.222850133 | El primero usa la única transición de descarga; el segundo combina 0.237533623 Ω en descarga y 0.208166644 Ω en carga. | 0.237533623 |
| `cycle_cell2_undeclared` | 0.185261851 | 0.177540167 | El primero usa la única transición de descarga; el segundo combina 0.185261851 Ω en descarga y 0.169818483 Ω en carga. | 0.185261851 |
| `cycle_cell2_30c` | 0.226546596 | 0.215820169 | El primero usa la única transición de descarga; el segundo combina 0.226546596 Ω en descarga y 0.205093741 Ω en carga. | 0.226546596 |

Las diferencias aparecen en los cinco ciclos combinados. No aparecen en los archivos de pulsos porque estos solo contienen transiciones de descarga comparables.

## 4. Criterio definitivo para la tesis

Para la sección 6.5 se recomienda usar exclusivamente las transiciones **reposo → descarga cercana a 1C** y reportar la mediana de dichas transiciones dentro de cada registro:

`R_app,desc = median(abs(Delta V_j)/abs(Delta I_j)), j perteneciente a reposo → descarga`.

Este criterio es preferible porque:

1. coincide con el objetivo original del diagnóstico, que solicitó comparar la entrada desde reposo a descarga cercana a 1C;
2. mantiene fija la dirección de corriente entre todos los experimentos;
3. no mezcla transiciones de carga tomadas en estados de carga y voltajes distintos con transiciones de descarga;
4. evita que los archivos de pulsos se resuman solo con descargas mientras los ciclos se resumen con una mezcla de descarga y carga.

No se aplicaron filtros adicionales, exclusiones posteriores ni correcciones de voltaje. La única selección adicional respecto al extractor general es `direction == "discharge"`.

Los ciclos combinados contienen una sola transición de descarga. Por ello, su valor recomendado es un **estimado de una transición**, aunque numéricamente la mediana de una observación sea ese mismo valor. No debe presentarse como evidencia de repetibilidad dentro del registro. La repetibilidad interna sí puede evaluarse en los pulsos, con siete y tres transiciones, respectivamente.

### Valores y razones definitivas

Usando `cycle_cell1_ambient = 0.036383074111115 Ω` como referencia direccional de descarga:

| Registro | `R_app,desc` [Ω] | N de transiciones de descarga | Razón respecto a `cycle_cell1_ambient` |
|---|---:|---:|---:|
| `cycle_cell1_ambient` | 0.036383074 | 1 | 1.0000 |
| `pulse_cell1_part1` | 0.221348508 | 7 | 6.0838 |
| `pulse_cell1_part2` | 0.228938929 | 3 | 6.2925 |
| `cycle_cell1_30c` | 0.240607211 | 1 | 6.6132 |
| `cycle_cell1_35c` | 0.237533623 | 1 | 6.5287 |
| `cycle_cell2_undeclared` | 0.185261851 | 1 | 5.0920 |
| `cycle_cell2_30c` | 0.226546596 | 1 | 6.2267 |

Al aplicar el mismo criterio de separación por la mayor brecha adyacente en escala logarítmica, `cycle_cell1_ambient` continúa siendo el único registro del grupo inferior. Por tanto, la conclusión cualitativa de separación entre la referencia y los otros seis registros no cambia, aunque sí cambian los valores y las razones que deben citarse.

## 5. Interpretación física y limitaciones

`R_app` es una respuesta aparente del camino eléctrico observado entre dos muestras adyacentes. Puede contener, de forma no separable con estos datos, contribuciones de la celda, polarización electroquímica rápida, contactos, cables, geometría y puntos de sensado, resolución temporal e instrumentación.

En consecuencia:

- `R_app` no es automáticamente la resistencia interna efectiva de la celda;
- no identifica por separado una resistencia externa del montaje;
- no permite atribuir la anomalía exclusivamente a contactos, cables o a la celda;
- no es el parámetro `R0` del modelo Panasonic y no debe transferirse a ese modelo;
- la ausencia de tiempo experimental medido limita la interpretación de cuán "instantánea" es la muestra cargada.

Los resultados solo justifican describir una **respuesta aparente elevada del camino observado** respecto al registro de referencia.

## 6. Corrección pendiente en la auditoría automática

Sí es necesario corregir, en una revisión posterior, las columnas numéricas `r_app_median_ohm`, `r_app_iqr_ohm` y `r_app_ratio_to_reference` de la auditoría automática para que la tabla destinada a la tesis use el criterio direccional de descarga. También debe reconciliarse la sección 3.1 del diagnóstico preliminar, que actualmente cita la referencia y las razones del criterio general.

La corrección de implementación sería hacer explícita la selección de `discharge_r_app_median_ohm` y `discharge_r_app_count` en la ruta TRC, en lugar de reutilizar silenciosamente `r_app_median_ohm`. No debe sustituirse el conteo sin revisar la regla de evidencia: los ciclos quedarían con una sola transición direccional y esa limitación debe conservarse como advertencia separada de la magnitud observada.

No es necesario cambiar los datos originales, las transiciones guardadas, el modelo Panasonic ni ningún parámetro. La clasificación cualitativa de separación relativa permanece, pero los valores definitivos que deben usarse en la tesis son los de descarga indicados en este informe. De acuerdo con la restricción de esta verificación, no se modificó todavía ningún código ni se sobrescribió ningún resultado previamente aceptado.
