# Contexto experimental y presentación de datos Uniandes

## Canal de aire compartido

`cycle_cell1_ambient` y `cycle_cell2_undeclared` se realizaron simultáneamente y compartieron el canal medido de aire `Temp7` (22.4896–23.4919 °C). El canal no se copia dentro de las señales de la celda 2. `Temp3` y `Temp4` permanecen como mediciones independientes de ubicación no documentada.

Los archivos tienen diferentes cantidades de registros y no incluyen tiempo medido. Por ello, la simultaneidad no demuestra alineación muestra a muestra y no se calculan diferencias de temperatura entre archivos.

## Etiquetas de condición

Las etiquetas `Aire 22–24 °C`, `Aire 29–31 °C` y `Aire 33–37 °C` resumen rangos medidos de `Temp7`; no son referencias programadas ni afirmaciones climáticas permanentes.

## Mapeo eléctrico corregido para las figuras de descarga

La corrección se aplica exclusivamente a `Voltage(V)` y a la capacidad descargada acumulada del intervalo seleccionado. Rionegro/montaje con resistencia externa usa `cycle_cell2_undeclared`; Cali/montaje con resistencia externa usa `cycle_cell2_30c`; Valledupar/montaje con resistencia externa usa el voltaje y la capacidad de `cycle_cell1_30c`; Rionegro/referencia usa `cycle_cell1_ambient`. La curva eléctrica de `cycle_cell1_35c` queda excluida.

Para Valledupar, las temperaturas permanecen en `cycle_cell1_35c` (`Temp1`, `Temp2`, `Temp7`). No se transfieren columnas térmicas hacia `cycle_cell1_30c` ni se combinan registros muestra a muestra. Carga, corriente y capacidades totales permanecen asociadas con sus archivos originales hasta que su procedencia sea verificada.

La auditoría confirmó que el eje anterior forzaba el primer registro a `SOC_ref=1` y el último a `SOC_ref=0` usando la propia capacidad final de cada archivo como denominador. Los CSV no demuestran que esos extremos sean 100 % y 0 % de SOC. Las figuras corregidas muestran `Q_desc = abs(Q-Q_inicio)` en Ah. Cada curva comienza en 0 Ah y termina en su capacidad descargada realmente registrada, sin normalización, sin extrapolación y sin afirmar una descarga completa.

## Auditoría de apariencia y método visual

Las descargas seleccionadas no contienen coordenadas de capacidad descargada duplicadas. El aspecto dentado procede principalmente de la cuantización/variación de voltaje y de saltos reales, no del trazado. Para mejorar la lectura se usa una mediana móvil centrada de 7 muestras, solamente dentro de la descarga continua. No se filtran ventanas que crucen saltos de voltaje ≥10 mV, cambios de corriente ≥0.05 A, extremos ni el punto de corte. Después se usa interpolación lineal; no se emplean splines y no puede aparecer sobreimpulso.

Los puntos originales se dibujan debajo de las líneas. Todos los indicadores de resistencia aparente, capacidad y desviación continúan calculándose con las mediciones originales; el panel de desviación se retira únicamente de la figura de presentación.

## Fuentes geográficas

- Límite de Colombia: Natural Earth, `ne_110m_admin_0_countries.geojson`, versión 5.1.2, SHA-256 `6866c877d39cba9c357620878839b336d569f8c662d3cfab4cb1dbe2d39c977f`.
- Rionegro: https://www.minvivienda.gov.co/node/38738
- Cali: https://ventanilla.minvivienda.gov.co/node/39659
- Valledupar: https://www.minvivienda.gov.co/node/39057

Las ciudades proporcionan contexto geográfico. No se incorporan referencias meteorológicas ni se afirma que un intervalo experimental sea habitual a partir de un solo día.
