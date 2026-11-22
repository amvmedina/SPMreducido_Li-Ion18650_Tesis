# Diagnóstico técnico para la auditoría de datos experimentales y transferencia del modelo

## Alcance y decisión de esta etapa

Este documento corresponde a la inspección previa solicitada para la sección **6.5. Auditoría de los datos experimentales y transferencia del modelo**. Se revisaron los modelos, scripts, datos, tablas y reportes numéricos disponibles. No se ejecutó ninguna calibración, no se modificaron parámetros Panasonic y no se implementó una auditoría nueva.

La conclusión principal es:

> Con los datos actuales puede implementarse directamente una auditoría automática de disponibilidad, integridad básica, capacidad, canales térmicos, respuesta serie aparente, aptitud de uso y prerrequisitos de transferencia. No puede realizarse una identificación electroquímica-térmica completa de la TRC18650 porque faltan una base de tiempo medida, una referencia de voltaje de equilibrio/OCP compatible con la celda y una caracterización inequívoca del montaje eléctrico.

El archivo `cycle_cell1_ambient` es la referencia propia más sólida para comparaciones en coordenadas de capacidad. Los demás ensayos muestran una respuesta serie aparente aproximadamente cinco a seis veces mayor. Esa observación justifica revisar el camino eléctrico, pero no identifica por sí sola una resistencia de contacto ni permite sustituir el parámetro interno `R0`.

## 1. Inventario verificado del proyecto

### 1.1 Modelo Panasonic aceptado

| Componente | Archivo autoritativo | Estado y contenido verificado |
|---|---|---|
| OCP Panasonic | `models/ocp_panasonic_qi_frozen.py` | OCP negativa MCMB Guo–White, OCP positiva de 13 nodos, capacidad y ventanas estequiométricas nominales efectivas; evaluación directa, sin optimizador. |
| Dinámica agrupada | `models/grouped_spm.py` y las rutinas Panasonic de validación | Modelo agrupado de cuatro estados, discretización ZOH exacta, corriente positiva en descarga, sin recorte de estados u OCP. |
| Acoplamiento electro-térmico | `models/electrothermal_spm.py` | Integra SOC de celda, estados eléctricos, calor y temperatura; no redefine los parámetros eléctricos. |
| Coeficiente entrópico | `models/panasonic_entropic.py` | Coeficiente de celda completa válido solamente en `0.10 <= SOC <= 0.90`; no extrapola. |
| Parámetros térmicos | `models/panasonic_thermal_parameters.py` | Modelo térmico efectivo de un nodo, calibrado con descarga Panasonic 1C a 25 °C en `SOC=0.90–0.10`. |

Parámetros finales verificados:

| Grupo | Valor | Procedencia |
|---|---:|---|
| `Q_cell` | 2.994979138519 Ah | Integración de la descarga Panasonic C/20. |
| `x_n,start`, `x_n,end` | 0.720451136448, 0.251849704066 | Solución nominal 3 del ajuste C/20; ventanas efectivas, no estequiometrías físicas únicas. |
| `x_p,start`, `x_p,end` | 0.224027933455, 0.953477823621 | Misma identificación C/20. |
| `Q_n`, `Q_p` | 23008.732269258333 C, 14780.898652565113 C | Derivados analíticamente de `Q_cell` y las ventanas efectivas. |
| `R0` | 0.0259594900586 Ω | Calibración HPPC condicionada a pulsos con SOC inicial `>=1/3`. |
| `alpha_p`, `d_n` | 3204.22765511 s, 0.000183029692927 s⁻¹ | Misma calibración HPPC. |
| `alpha_n`, `d_p` | 3000 s, 0.0007 s⁻¹ | Fijos durante la calibración HPPC. |
| `Cth`, `Gth` | 117.52317800099644 J/K, 0.10911208974639978 W/K | Ajuste térmico 1C a 25 °C, usando calor irreversible basado en voltaje medido. |

Las constantes eléctricas aceptadas son `R=8.314472 J mol⁻¹ K⁻¹` y `F=96487 C mol⁻¹`.

### 1.2 Datos Panasonic usados por la rama final

El repositorio contiene una copia local bajo `Panasonic-18650PF-Data-master/Panasonic-18650PF-Data-master/Panasonic 18650PF Data`. El manifiesto conserva además las rutas originales externas. Los cinco archivos principales son:

| Uso | Archivo MAT seleccionado | Papel real en el flujo |
|---|---|---|
| Referencia de baja tasa | `25degC/C20 OCV and 1C discharge tests_start_of_tests/05-08-17_13.26 C20 OCV Test_C20_25dC.mat` | Capacidad C/20 y ajuste de la representación OCP de celda completa. |
| Calibración dinámica | `25degC/5 pulse disch/03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat` | HPPC a 0.5, 1, 2, 4 y 6C; el conjunto nominal final usa 40 pulsos con SOC inicial `>=1/3`. |
| Calibración térmica y validación eléctrica 1C | `25degC/C20 OCV and 1C discharge tests_start_of_tests/03-09-17_17.59 3349_Dis1C_1.mat` | Ajuste de `Cth,Gth` y validación eléctrica directa; son usos distintos sobre el mismo registro y deben declararse como tales. |
| Validación dinámica independiente | `25degC/Drive cycles/03-21-17_00.29 25degC_UDDS_Pan18650PF.mat` | Validación directa eléctrica y térmica, sin reajuste. |
| Validación de carga | `25degC/C20 OCV and 1C discharge tests_start_of_tests/03-09-17_19.03 3349_Charge2.mat` | Validación CC–CV directa, sin controlador simulado ni reajuste. |

Los MAT Panasonic contienen `TimeStamp`, `Voltage`, `Current`, `Ah`, `Wh`, `Power`, `Battery_Temp_degC`, `Time` y `Chamber_Temp_degC`. La corriente original es negativa en descarga y se invierte una sola vez en la frontera de carga del modelo.

La colección local también contiene HPPC y otros perfiles a 10, 0, −10 y −20 °C. Existe una auditoría multitemperatura ejecutable, pero concluyó que `R0(T)`, `alpha_p(T)` y `d_n(T)` no eran suficientemente recuperables; no se implementó una ley térmica o de Arrhenius. Por tanto, el modelo final debe considerarse condicionado principalmente a 25 °C.

### 1.3 Scripts Panasonic y resultados numéricos existentes

| Etapa | Script verificable | Salida numérica principal |
|---|---|---|
| C/20 y cierre OCP | `scripts/verify_panasonic_c20_qi_frozen.py` | `panasonic_ocp_identification_closure.md`; RMSE 5–95 % = 1.085749161 mV. |
| Auditoría HPPC | `scripts/audit_panasonic_hppc.py` | `panasonic_hppc_pulse_summary.csv` y `panasonic_hppc_audit.md`. |
| Calibración HPPC | `scripts/calibrate_panasonic_hppc_primary.py` y `scripts/diagnose_panasonic_hppc_primary_no_low_soc.py` | Parámetros nominales condicionados a SOC `>=1/3`; RMSE balanceado 36.835198 mV. |
| Validación HPPC de bajo SOC | `scripts/validate_panasonic_hppc_low_soc.py` | Evaluación directa de pulsos reservados; sin reajuste. |
| Validación 1C | `scripts/validate_panasonic_1c_forward.py` | RMSE global 114.362188 mV; deterioro marcado a bajo SOC. |
| Validación UDDS | `scripts/validate_panasonic_udds_forward.py` | RMSE global 30.988192 mV; error mayor a bajo SOC y en escalones grandes. |
| Validación CC–CV | `scripts/validate_panasonic_cccv_charge.py` | RMSE global 61.033423 mV; sesgos distintos en reposo, CC y CV. |
| Calibración térmica 1C | `scripts/calibrate_panasonic_1c_thermal.py` | RMSE 0.143622065 K, `Cth` y `Gth` finales. |
| Validaciones térmicas | `scripts/validate_panasonic_udds_thermal.py`, `scripts/validate_panasonic_cccv_thermal.py` | Evaluaciones directas A/B; no reajuste. |

Estas cifras demuestran que se realizaron cálculos numéricos. No se infirió la existencia de un análisis a partir de una figura aislada.

### 1.4 Datos TRC18650 y dependencia de rutas externas

Los siete CSV originales existen actualmente, pero están fuera del repositorio, bajo `C:/Users/amvel/Documents/Andes/Tesis/TRC18650experimental data`. El código depende de esas rutas absolutas. Para una auditoría portátil futura conviene recibir las rutas mediante configuración, sin copiar ni alterar los originales.

Todos los archivos contienen voltaje, corriente, `Capacity(Ah)`, `Total Capacity(Ah)`, paso y modo. Ninguno contiene tiempo transcurrido ni marca temporal. Los cinco archivos de celda 1 contienen `Temp1`, `Temp2` y `Temp7`; los dos de celda 2 contienen `Temp3` y `Temp4`.

Roles documentados:

- `Temp1`: superficie central de la celda.
- `Temp2`: superficie cerca de un extremo.
- `Temp7`: aire de cámara/ambiente.
- `Temp3` y `Temp4`: temperaturas medidas, pero su ubicación no está documentada.

`cycle_cell1_ambient` y `cycle_cell2_undeclared` se realizaron simultáneamente y compartieron el `Temp7` medido, 22.4896–23.4919 °C. Los registros tienen longitudes distintas y no hay tiempo medido; no está permitido alinearlos muestra a muestra.

## 2. Inventario y aptitud de los registros TRC18650

Las duraciones de la tabla son **provisionales**. Se obtuvieron a partir de capacidad/corriente, no de una base de tiempo registrada. Por ello no permiten auditar duplicados temporales, inversión del tiempo, pausas perdidas ni sincronización entre archivos.

| Registro | Ensayo y condición | Variables / muestras | Duración provisional | Condición inicial disponible | Hallazgos e interrupciones | Uso defendible actual |
|---|---|---|---:|---|---|---|
| `pulse_cell1_part1` | Pulsos/restos de descarga ≈1C, celda 1, ambiente no declarado | V, I, capacidades, paso/modo, Temp1/2/7; 59 990 | 5986.15 s (1.663 h) | No hay SOC inicial certificado; inicia en una secuencia ya activa | Siete transiciones; archivo incompleto. `R_app` mediana 0.221349 Ω. | Diagnóstico de no idealidad y respuesta térmica cualitativa; no calibración. |
| `pulse_cell1_part2` | Continuación del anterior | Mismos canales inferidos por posición; 14 275 | 1424.09 s (0.396 h) | Continuación presumida, sin tiempo absoluto | No tiene encabezado. La frontera conserva V/I/capacidad, pero no excluye un intervalo perdido. Tres transiciones; `R_app` 0.228939 Ω. | Igual que parte 1. |
| `cycle_cell1_ambient` | Descarga ≈1C y carga CC–CV, celda 1, Aire 22–24 °C, montaje de referencia | V, I, capacidades, paso/modo, Temp1/2/7; 4241 | 12697.54 s (3.527 h) | La descarga está precedida por reposo y una carga parcial visible; no hay SOC inicial independiente | Descarga 2.8999 Ah y recarga 2.8999 Ah. Diez valores faltantes en Temp7. `R_app` de descarga 0.036383 Ω. | Mejor referencia propia: capacidad y comparación V–capacidad; térmica solo cualitativa mientras falte tiempo medido. |
| `cycle_cell1_30c` | Descarga ≈1C y carga CC–CV, celda 1, Aire 29–31 °C | V, I, capacidades, paso/modo, Temp1/2/7; 4372 | 13085.03 s (3.635 h) | Descarga precedida por reposo; no hay SOC inicial independiente | Descarga/carga 2.4894/2.4653 Ah. `R_app` de descarga 0.240607 Ω. | Caso de montaje con respuesta serie elevada; diagnóstico cualitativo. |
| `cycle_cell1_35c` | Descarga ≈1C y carga CC–CV, celda 1, Aire 33–37 °C | V, I, capacidades, paso/modo, Temp1/2/7; 4937 | 14776.82 s (4.105 h) | Descarga precedida por reposo; no hay SOC inicial independiente | Descarga/carga 2.6171/2.5185 Ah. Seis cambios de Temp7 superan el cribado robusto de 0.5 °C. `R_app` de descarga 0.237534 Ω. | Temperaturas descriptivas y caso de montaje; no identificación electro-térmica. |
| `cycle_cell2_undeclared` | Descarga ≈1C y carga CC–CV, celda 2, Aire 22–24 °C compartido | V, I, capacidades, paso/modo, Temp3/4; 6532 | 19546.34 s (5.430 h) | Descarga precedida por reposo; no hay SOC inicial independiente | Descarga/carga 2.7133/2.6896 Ah. Temp3/4 sin ubicación; Temp7 está en otro archivo sin alineación. `R_app` de descarga 0.185262 Ω. | Comparación entre montajes/celdas y diagnóstico cualitativo. |
| `cycle_cell2_30c` | Descarga ≈1C y carga CC–CV, celda 2, Aire 29–31 °C | V, I, capacidades, paso/modo, Temp3/4; 4501 | 13471.71 s (3.742 h) | Descarga precedida por reposo; no hay SOC inicial independiente | Descarga/carga 2.5994/2.5053 Ah. No contiene ambiente. `R_app` de descarga 0.226547 Ω. | Diagnóstico cualitativo; no calibración. |

La corriente Chroma es negativa en descarga y positiva en carga. La conversión requerida es `I_model=-I_raw`, únicamente en el cargador de datos.

### 2.1 Correspondencia experimental que no debe generalizarse

La corrección confirmada para las figuras de descarga reasigna únicamente voltaje y capacidad descargada:

- Rionegro, montaje de referencia: `cycle_cell1_ambient`.
- Rionegro, montaje con resistencia externa: `cycle_cell2_undeclared`.
- Cali, montaje con resistencia externa: `cycle_cell2_30c`.
- Valledupar, montaje con resistencia externa: voltaje/capacidad de `cycle_cell1_30c`; temperaturas permanecen en `cycle_cell1_35c`.
- La curva eléctrica de descarga de `cycle_cell1_35c` queda excluida de esas comparaciones.

Esta reasignación no autoriza a transferir corriente, carga, capacidad total o temperatura entre archivos. En particular, el registro eléctrico y el térmico de Valledupar no pueden combinarse muestra a muestra.

## 3. Estado real de los análisis del montaje

### 3.1 Evidencia experimental comprobada

La rutina `extract_rest_to_load_transitions()` aplica de forma reproducible el mismo criterio local a todos los CSV y calcula:

`R_app = abs(ΔV)/abs(ΔI)`.

El archivo `own_data_apparent_resistance_transitions.csv` conserva cada transición. El grupo de referencia tiene mediana/IQR 0.035875537/0.000507537 Ω. Los pulsos y ciclos restantes caen en un grupo elevado, con razones aproximadas de 4.95 a 6.48 respecto a la referencia. Las diez transiciones de los dos archivos de pulsos se mantienen entre 0.218014 y 0.237146 Ω.

Esto demuestra una separación reproducible de la respuesta terminal observada y descarta que el pulso sea un único valor atípico.

### 3.2 Anomalías que justifican revisar el montaje

- A corriente de descarga cercana a −2.90 A y voltajes de reposo próximos a 4.18–4.19 V, los ciclos elevados presentan caídas iniciales de aproximadamente 0.54–0.70 V, frente a 0.1055 V en la referencia.
- Los archivos de 30 y 35 °C de celda 1 y ambos archivos de celda 2 pertenecen al mismo grupo elevado que los pulsos.
- Las diferencias de capacidad y forma de voltaje no pueden atribuirse de manera aislada a temperatura, porque coexisten diferencias de montaje, celda y terminación.

### 3.3 Hipótesis que los registros no confirman

Los CSV no documentan dónde se conectaron los terminales de sensado de voltaje, si se usó sensado Kelvin, qué resistencia tenían cables/contactos ni dónde se disipó la potencia adicional. La denominación “montaje con resistencia externa” proviene del contexto experimental confirmado por el usuario, no de una identificación eléctrica completa en los archivos.

Por tanto:

- `R_app` no es la resistencia interna de la celda ni el `R0` del modelo.
- No puede afirmarse que toda la diferencia esté en contactos, cables o celda.
- No puede calcularse el calor interno como si toda la caída adicional ocurriera dentro de la celda.
- La comparación entre celda 1 y celda 2 no es un ensayo antes/después sobre la misma celda.

## 4. Requisitos de transferencia Panasonic → TRC18650

| Información/parámetro | Evidencia Panasonic | Evidencia TRC disponible | Decisión para TRC |
|---|---|---|---|
| Capacidad de celda | C/20 medida: 2.994979138519 Ah | Descarga/carga completa de referencia: 2.8999 Ah | Adaptar a 2.8999 Ah; no transferir la capacidad Panasonic. |
| Estado inicial | Ventanas efectivas y Ah acumulados con tiempo medido | Reposos y contadores Chroma, pero sin SOC inicial independiente | Puede inicializarse solo de forma condicional en coordenadas de capacidad; hace falta protocolo de acondicionamiento documentado. |
| OCP negativa y positiva | MCMB fija + OCP positiva de 13 nodos aceptada contra C/20 Panasonic | No hay química, curvas de media celda ni descarga de baja tasa TRC | No transferible como representación validada; se requiere referencia TRC de baja tasa y una hipótesis química explícita. |
| Ventanas y `Q_n,Q_p` | Ajustadas/derivadas para Panasonic | No identificadas | Deben adaptarse conjuntamente con la referencia de voltaje; no basta cambiar `Q_cell`. |
| `R0, alpha_p, d_n` | HPPC Panasonic 25 °C, condicionado a SOC `>=1/3` | Pulsos propios contaminados por respuesta serie y sin tiempo medido | No identificables con los datos actuales. |
| `alpha_n,d_p` | Fijos en el flujo Panasonic | Sin evidencia TRC | No deben asumirse parámetros TRC identificados. |
| `Cth,Gth` | Ajuste Panasonic 1C con tiempo y ambiente medidos | Temperaturas TRC, pero sin tiempo medido; calentamiento confuso en montajes elevados | No transferibles como parámetros físicos TRC ni identificables actualmente. |
| Coeficiente entrópico | Coeficiente de celda Panasonic, SOC 0.10–0.90 | Ninguna medición/referencia TRC | No transferible sin justificar química/celda. |
| Dependencia con temperatura | No se aceptó ninguna ley multitemperatura Panasonic | Tres intervalos térmicos propios, pero montajes y tiempo confunden el efecto | No puede identificarse una ley térmica de parámetros. |

El diagnóstico de transferencia ya ejecutado (`scripts/audit_panasonic_to_trc18650_transfer.py`) produjo comparaciones directas, pero sus RMSE están rotulados correctamente como **diagnósticos de transferencia**, no como validación de un modelo TRC. El propio análisis concluye que se requieren al menos adaptación de capacidad/estado y de referencia de voltaje/OCP; la adaptación dinámica está indicada, pero no es identificable sin tiempo medido.

## 5. Experimentos e información faltantes, priorizados

Antes de repetir ensayos debe intentarse recuperar la receta Chroma, la base de tiempo original y la documentación del cableado. Si esos elementos existen, pueden evitar repeticiones parciales.

### Prioridad 0 — Trazabilidad del montaje y del tiempo

**Objetivo:** determinar qué incluye el voltaje registrado y convertir las duraciones provisionales en tiempo medido.

**Procedimiento preliminar:** recuperar configuración Chroma, periodo de adquisición, esquema de fuerza/sensado, punto de conexión de los cables de voltaje, calibre/longitud de cables, identificación de celda y ubicación de sensores. Verificar el montaje con una carga o patrón conocido antes de usar otra celda.

**Tiempo:** no es estimable con una receta disponible; depende de la documentación y del acceso al equipo. No debe inventarse una duración.

### Prioridad 1 — Referencia TRC de baja tasa

**Objetivo:** medir capacidad y forma de voltaje cuasiestática específica de la TRC, con el montaje de referencia.

**Procedimiento preliminar:** acondicionar la celda con límites documentados; registrar tiempo, corriente, voltaje, Ah, temperaturas superficial y ambiente; ejecutar al menos una descarga C/20 completa y, preferiblemente, una carga C/20 separada para evaluar dependencia con la dirección. Conservar reposos y cortes como segmentos explícitos.

**Tiempo activo calculable:** con `Q=2.8999 Ah`, C/20 equivale aproximadamente a 0.145 A y una rama ideal completa requiere unas **20 h**. Descarga y carga requieren unas **40 h activas**. Preparación, CV y reposos no pueden estimarse hasta fijar la receta; deben reportarse por separado.

Este experimento no repite la capacidad 1C ya disponible: añade la referencia de baja tasa que falta.

### Prioridad 2 — HPPC limpio y temporizado

**Objetivo:** identificar respuesta óhmica efectiva, difusión/relajación y cinética sin incorporar la resistencia externa del montaje.

**Procedimiento preliminar:** usar el montaje verificado, tiempo medido y estados de carga establecidos por coulomb counting desde una condición inicial documentada. Aplicar pulsos y reposos en varios SOC, manteniendo al menos los niveles 0.5, 1 y 2C que tienen soporte común en la referencia Panasonic. Reservar una fracción de pulsos o niveles de SOC para validación.

**Tiempo verificable de referencia:** el archivo HPPC Panasonic 25 °C ocupa 97 599.4 s, **27.11 h**, e incluye grandes separaciones entre grupos. Los pulsos son aproximadamente de 10 s y muchos reposos son de aproximadamente 1200 s. Como las descargas de ajuste de SOC están en archivos/intervalos separados, ese valor no es una duración completa transferible a una receta TRC. Para 26 combinaciones comunes, solo 26 reposos de 20 min ya representan 8.67 h; deben añadirse acondicionamiento y desplazamientos de SOC.

El pulso Chroma existente no sustituye este ensayo: está afectado y carece de tiempo medido.

### Prioridad 3 — Ensayo térmico limpio

**Objetivo:** identificar o validar el modelo térmico con una fuente de calor eléctricamente interpretable.

**Procedimiento preliminar:** descarga 1C con tiempo medido, montaje de referencia, temperatura central, temperatura cerca del extremo y aire de cámara, seguida de un reposo de enfriamiento registrado. Antes de repetirla, intentar recuperar el reloj del ensayo ambiente existente.

**Tiempo activo calculable:** una descarga nominal 1C de 2.8999 Ah requiere aproximadamente **1 h**. Como referencia, la descarga Panasonic aceptada duró 3474.37 s (57.91 min), y su intervalo térmico SOC 0.90–0.10 duró 2974.93 s (49.58 min). El reposo de enfriamiento y la preparación no tienen duración TRC verificada.

### Prioridad 4 — Validación independiente

**Objetivo:** evaluar transferencia sin reutilizar el perfil de calibración.

**Procedimiento preliminar:** reservar una descarga/carga completa o ejecutar un perfil dinámico diferente al HPPC. Debe comenzar desde un SOC documentado, conservar corriente regenerativa si existe y registrar tiempo y temperaturas.

**Tiempo de referencia verificable:** el UDDS Panasonic usado dura 22 447.866 s, **6.2355 h**. Es solo una referencia de tamaño del registro; no incluye la preparación de una celda TRC ni garantiza que el mismo perfil sea seguro o apropiado.

### Lo que no debe repetirse automáticamente

- No repetir una medición de capacidad 1C únicamente para obtener `Q_cell`: la referencia ambiente ya soporta 2.8999 Ah en coordenadas de capacidad.
- No repetir todos los ensayos a tres temperaturas antes de corregir el montaje y el tiempo; se reproduciría la misma ambigüedad.
- No intentar identificar OCP, dinámica y parámetros térmicos simultáneamente con un único ciclo CC–CV.

## 6. Propuesta de auditoría automática mínima

La primera implementación puede limitarse a siete reglas deterministas:

1. **Descubrimiento y esquema.** Verificar existencia, encabezado, columnas obligatorias, unidades declaradas, identificador de celda y procedencia. No aceptar `part2` sin el esquema de `part1` explícitamente configurado.
2. **Base de tiempo y signo.** Clasificar el tiempo como medido, provisional o ausente; detectar el signo de corriente y convertirlo solo en una copia destinada al modelo. Prohibir métricas dinámicas si el tiempo no es medido.
3. **Segmentación y capacidad.** Segmentar por paso/modo/corriente, reportar carga, descarga, reposos, capacidad recorrida, corte observado y completitud. No llamar SOC a una fracción normalizada del registro sin referencia independiente.
4. **Canales térmicos y alineación.** Registrar rol, faltantes, saltos y rango de cada termocupla. Prohibir operaciones muestra a muestra entre archivos sin tiempo común. No asignar ubicación a Temp3/Temp4.
5. **Respuesta serie relativa.** Extraer transiciones adyacentes con una regla común, calcular `R_app`, comparar medianas entre ensayos y marcar grupos relativos. Rotular siempre el resultado como diagnóstico del camino observado, nunca como `R0`.
6. **Matriz de aptitud.** Emitir por registro: calibración eléctrica, validación eléctrica, calibración térmica, validación térmica y uso diagnóstico, cada decisión con evidencia textual.
7. **Puerta de transferencia y siguiente acción.** Comparar disponibilidad de capacidad, OCP/referencia de voltaje, tiempo dinámico, temperatura/ambiente y experimento independiente. Si falta un prerrequisito, detener la transferencia cuantitativa y recomendar el experimento mínimo correspondiente; nunca iniciar un ajuste automáticamente.

La salida debe incluir una tabla de inventario, una tabla de aptitud, banderas de revisión, prerrequisitos faltantes y una recomendación priorizada. No necesita aprendizaje automático, corrección automática de voltaje ni optimización.

## 7. Archivos y funciones reutilizables

### Carga y auditoría básica

- `analysis/experimental_data_audit.py`
  - `load_chroma_experiment()`
  - `infer_raw_current_sign()`
  - `estimate_sample_interval_s()` — debe conservar el rótulo provisional.
  - `segment_capacities()`
  - `extract_rest_to_load_transitions()`
  - `audit_experiment()`
  - `largest_log_gap_reference_group()`
- `scripts/audit_own_trc18650_data.py`
  - `experiment_metadata()` y `load_all()` como configuración actual de los siete archivos.

### Segmentación y comparaciones propias

- `analysis/own_case_study_reference.py`
  - `segment_complete_cycle()`
  - `discharge_curve()` y `charge_curve()`
  - `relative_series_path_detector()`
- `scripts/build_own_trc18650_case_study.py`
  - generación de tablas y reportes, no como fuente primaria de datos.
- `scripts/build_uniandes_context_figures.py` y `scripts/build_uniandes_charge_temperature_figures.py`
  - contienen el mapeo corregido y la correspondencia térmica; las transformaciones visuales no deben entrar en métricas.

### Reglas de decisión ya implementadas

- `analysis/battery_dataset_assessment.py`
  - `evaluate_data_quality_gate()`
  - `evaluate_relative_nonideality_gate()`
  - `convert_current_for_model()`
  - `screen_ocp_domain()`
  - `evaluate_transfer_gate()`
  - `require_thermal_dynamic_support()`
- `scripts/run_battery_dataset_assessment.py`
  - ejemplo de ensamblaje del informe final.
- `tests/test_battery_dataset_assessment.py`
  - pruebas de las prohibiciones y puertas metodológicas.

### Simulación y métricas

- `models/grouped_spm.py`: propagación ZOH, estados superficiales, voltaje y auditoría de dominio.
- `models/electrothermal_spm.py`: SOC y acoplamiento térmico.
- `scripts/audit_panasonic_to_trc18650_transfer.py`: comparación jerárquica sin reajuste.
- Las funciones `voltage_metrics()` y `region_masks()` de ese script pueden reutilizarse solo después de que la puerta de aptitud autorice métricas cuantitativas.

### Salidas numéricas que deben tratarse como evidencia

- `thesis_data/own_data_manifest.csv`
- `thesis_data/own_data_quality_audit.csv`
- `thesis_data/own_data_apparent_resistance_transitions.csv`
- `thesis_data/own_case_study_cross_experiment_metrics.csv`
- `thesis_data/own_model_transfer_metrics.csv`
- `thesis_data/uniandes_corrected_discharge_mapping_audit.csv`
- `thesis_data/uniandes_charge_temperature_correspondence.csv`

Las figuras son productos de presentación; no sustituyen estas tablas ni los scripts que las generan.

## 8. Implementable ahora y dependencias pendientes

### Implementable directamente con los datos actuales

- Descubrimiento y normalización de los siete CSV.
- Inventario de canales, muestras, rangos, pasos y modos.
- Segmentación de descarga, carga CC–CV y reposos.
- Capacidades recorridas basadas en el contador Chroma.
- Conversión explícita del signo de corriente.
- Auditoría de temperaturas, faltantes y saltos.
- Diagnóstico relativo `R_app` y clasificación del grupo elevado.
- Tabla de aptitud y puerta jerárquica de transferencia.
- Comparaciones de voltaje frente a capacidad del ensayo ambiente, rotuladas como diagnóstico de transferencia y no como validación TRC.
- Recomendación automática del siguiente requisito faltante.

### Requiere información adicional

- Tiempo medido o registro del periodo real de adquisición y continuidad entre archivos.
- Receta Chroma completa, cortes, reposos, eventos manuales y estado de acondicionamiento.
- Diagrama del camino de fuerza y sensado de voltaje; ubicación de terminales Kelvin si existieron.
- Identificación/datasheet y química de la TRC18650.
- Referencia TRC de baja tasa para OCP/voltaje de equilibrio y ventanas efectivas.
- Pulsos limpios para parámetros dinámicos.
- Posición documentada de Temp3/Temp4 y una base temporal para alinear el `Temp7` compartido.
- Ensayo térmico con tiempo medido y calor eléctricamente defendible.
- Perfil reservado para validación independiente.

## Stop point

El diagnóstico termina aquí. Se creó solamente este informe. No se modificó ningún archivo existente, no se ejecutó el modelo, no se recalibró ningún parámetro y no se implementó todavía el sistema automático.
