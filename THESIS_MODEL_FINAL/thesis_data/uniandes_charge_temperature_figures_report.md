# Figuras de carga y temperatura de los ensayos propios

## Correspondencia de registros

La reasignación eléctrica anterior se mantiene limitada a voltaje y capacidad de **descarga**. No se trasladó a la carga, la corriente ni las temperaturas. Cada curva de carga usa voltaje, corriente y `Capacity(Ah)` del mismo CSV.

Las figuras térmicas usan solamente archivos de celda 1 que contienen Temp1, Temp2 y Temp7 en el mismo registro. Temp3 y Temp4 de celda 2 no se muestran porque sus ubicaciones no están documentadas. El Temp7 compartido por los dos ensayos simultáneos de Aire 22–24 °C no se alinea muestra a muestra con el archivo de celda 2. En los paneles térmicos, `Capacity(Ah)` se usa solo como coordenada del mismo archivo para las temperaturas; esto no modifica la reasignación de la curva eléctrica de descarga.

## Carga CC–CV

La capacidad cargada se calcula como `Q_chg,i = abs(Capacity(Ah)_i-Capacity(Ah)_inicio_carga)`. En todos los registros seleccionados, CC y CV pertenecen al mismo paso 5, las filas son contiguas y el contador es monótono; no se detecta reinicio entre etapas.

| Curva | Capacidad cargada (Ah) | Voltaje final (V) | Corriente final (A) | Terminación |
|---|---:|---:|---:|---|
| Rionegro · Aire 22–24 °C | 2.6896 | 4.2027 | 0.1378 | cambio de carga CV a reposo a 4.2027 V y 0.1378 A; la corriente venía disminuyendo, pero no alcanza un umbral documentado; el CSV no contiene receta, código de terminación, evento manual ni bandera de error; la causa definitiva no puede establecerse |
| Cali · Aire 29–31 °C | 2.5053 | 4.2042 | 0.2917 | cambio de carga CV a reposo a 4.2042 V y 0.2917 A; la corriente venía disminuyendo, pero no alcanza un umbral documentado; el CSV no contiene receta, código de terminación, evento manual ni bandera de error; la causa definitiva no puede establecerse |
| Valledupar · Aire 33–37 °C | 2.5185 | 4.2034 | 0.3553 | cambio de carga CV a reposo a 4.2034 V y 0.3553 A; la corriente venía disminuyendo, pero no alcanza un umbral documentado; el CSV no contiene receta, código de terminación, evento manual ni bandera de error; la causa definitiva no puede establecerse |
| Montaje de referencia | 2.8999 | 4.2033 | 0.1333 | cambio de carga CV a reposo cerca de 2.900 Ah; compatible con un límite de capacidad; el CSV no contiene receta, código de terminación, evento manual ni bandera de error; la causa definitiva no puede establecerse |

## Temperaturas por fase

| Fase | Condición | Sensor | Mínimo (°C) | Máximo (°C) |
|---|---|---|---:|---:|
| discharge | Rionegro · Aire 22–24 °C | Temp1 | 23.8919 | 33.0483 |
| discharge | Rionegro · Aire 22–24 °C | Temp2 | 23.6380 | 32.1979 |
| discharge | Rionegro · Aire 22–24 °C | Temp7 | 22.4896 | 23.4467 |
| discharge | Cali · Aire 29–31 °C | Temp1 | 30.0534 | 37.8710 |
| discharge | Cali · Aire 29–31 °C | Temp2 | 29.8657 | 36.3229 |
| discharge | Cali · Aire 29–31 °C | Temp7 | 29.4847 | 30.9245 |
| discharge | Valledupar · Aire 33–37 °C | Temp1 | 34.6002 | 44.4022 |
| discharge | Valledupar · Aire 33–37 °C | Temp2 | 34.4398 | 43.0099 |
| discharge | Valledupar · Aire 33–37 °C | Temp7 | 33.3432 | 36.9376 |
| charge | Rionegro · Aire 22–24 °C | Temp1 | 26.3592 | 31.8844 |
| charge | Rionegro · Aire 22–24 °C | Temp2 | 25.9556 | 31.1914 |
| charge | Rionegro · Aire 22–24 °C | Temp7 | 22.5108 | 23.4919 |
| charge | Cali · Aire 29–31 °C | Temp1 | 30.2842 | 36.0649 |
| charge | Cali · Aire 29–31 °C | Temp2 | 30.0027 | 34.9878 |
| charge | Cali · Aire 29–31 °C | Temp7 | 29.4393 | 30.7055 |
| charge | Valledupar · Aire 33–37 °C | Temp1 | 35.5307 | 42.7163 |
| charge | Valledupar · Aire 33–37 °C | Temp2 | 35.1656 | 41.7465 |
| charge | Valledupar · Aire 33–37 °C | Temp7 | 33.6238 | 36.2619 |

## Procesamiento gráfico

Las cuatro figuras se construyen a partir de las muestras originales. No se usan tiempos reconstruidos, no se incluyen reposos sobre ejes de capacidad y no se extrapolan registros. Como única excepción visual, Temp7 de Valledupar usa una mediana móvil centrada de 11 muestras para reducir fluctuaciones rápidas sin sobreimpulso. En carga, el filtro se aplica por separado dentro de CC y CV; las primeras y últimas cinco muestras de cada tramo permanecen sin alterar. Temp1 y Temp2 no se suavizan. Todos los rangos y resultados se calculan con las muestras originales. Las ciudades se emplean solo como referencias de las condiciones térmicas.
