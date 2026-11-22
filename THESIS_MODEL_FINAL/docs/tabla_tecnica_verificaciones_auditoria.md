# Criterios técnicos de la auditoría E01–E08

Todos los formatos se convierten primero a `CommonAuditRecord`. La función única `evaluate_common_audit` aplica las reglas siguientes sin consultar marca, nombre comercial, nombre de archivo o modalidad de carga. Los adaptadores CSV y Panasonic solo leen señales y metadatos.

| ID | Función responsable | Entradas | Criterio exacto | Estados posibles | Acción generada |
|---|---|---|---|---|---|
| E01 | `evaluate_common_audit` | Voltaje, corriente y capacidad | CONFORME si todos los ensayos contienen las tres series, con al menos dos muestras y la misma longitud; REVISAR si solo una parte; FALTANTE si ninguno | CONFORME, REVISAR, FALTANTE | Incluir o corregir las variables ausentes. No intervienen simulaciones. |
| E02 | `evaluate_common_audit` / `_capacity_complete` | Voltaje, corriente, capacidad y límites declarados por ensayo | CONFORME si al menos una carga alcanza su límite superior declarado o una descarga alcanza su límite inferior declarado; FALTANTE si no hay límites suficientes o ningún ensayo los alcanza. Los ensayos de pulsos y perfiles dinámicos pueden marcarse como no destinados a esta comprobación | CONFORME, FALTANTE | Declarar el límite correspondiente o aportar una carga/descarga completa. No existen límites universales codificados para archivos nuevos. |
| E03 | `evaluate_common_audit` / `_invalid_details` | Todas las series suministradas | CONFORME si no hay muestras inválidas ni longitudes distintas; REVISAR si se identifica cualquier variable afectada; FALTANTE sin ensayos | CONFORME, REVISAR, FALTANTE | Informa ensayo, variable y número de muestras. No depende de E02. |
| E04 | `evaluate_common_audit` / `_rest_to_discharge` | Una referencia y otro montaje con transición reposo–descarga | Calcula la mediana de `abs(ΔV/ΔI)` y aplica la agrupación relativa por mayor separación logarítmica ya verificada. NO EVALUADO sin dos montajes comparables; REVISAR ante respuesta relativa elevada; CONFORME si no se detecta elevación | CONFORME, REVISAR, NO EVALUADO | Revisar montaje y conexiones sin atribuir una causa física única. |
| E05 | `evaluate_common_audit` | Ensayo declarado de baja corriente, capacidad nominal, corriente y límite de descarga | CONFORME si la mediana de corriente activa es ≤0.10C y la descarga alcanza su límite; REVISAR si faltan capacidad nominal o límite; NO APTO si existe pero no combina baja corriente y descarga completa; FALTANTE si no existe | CONFORME, REVISAR, FALTANTE, NO APTO | Completar metadatos o realizar una descarga C/20 completa. |
| E06 | `evaluate_common_audit` | Ensayo HPPC/pulsos, escala temporal, corriente y voltaje | CONFORME con tiempo registrado o reconstruido desde un intervalo constante y continuo confirmado, reposos y al menos dos transiciones reposo–descarga; REVISAR si la secuencia de pulsos es suficiente pero la escala temporal requiere confirmación; NO APTO ante una secuencia temporal o de pulsos inconsistente; FALTANTE sin ensayo. Si E04 eleva el mismo registro, E06 pasa a REVISAR | CONFORME, REVISAR, FALTANTE, NO APTO | Confirmar el intervalo, corregir el protocolo o revisar el montaje señalado por E04. La independencia se evalúa solo en E08. |
| E07 | `evaluate_common_audit` | Temperatura de celda, temperatura ambiente y escala temporal | CONFORME si un ensayo contiene ambas temperaturas y tiempo registrado o reconstruido desde un intervalo constante/continuo confirmado; REVISAR si ambas temperaturas existen pero el intervalo solo fue estimado o no está confirmado; FALTANTE si falta alguna temperatura requerida; NO APTO ante datos térmicos o temporales inconsistentes | CONFORME, REVISAR, FALTANTE, NO APTO | Indica la medición faltante o solicita confirmar el intervalo de muestreo, sin exigir repetir el ensayo si la información puede recuperarse. |
| E08 | `evaluate_common_audit` | Funciones experimentales declaradas y variables E01 | CONFORME si existe un ensayo de validación con variables válidas y sin funciones de baja tasa, pulsos o calibración térmica; NO APTO si los candidatos se usan también en calibración o carecen de variables; FALTANTE si no existen | CONFORME, FALTANTE, NO APTO | Reservar un ensayo válido que no se utilice en calibración. |

## Constantes del criterio

- Umbral de reposo: `|I| ≤ 0.05 A`.
- Máxima tasa admitida como baja corriente: `0.10C`.
- Mínimo para E06: dos transiciones reposo–descarga.
- Los límites de voltaje de E02 y E05 proceden de los metadatos del ensayo; no se asignan por marca.

## Clasificación temporal

- `recorded`: tiempo incluido directamente en el registro.
- `confirmed_interval`: tiempo reconstruido con un intervalo indicado por el usuario, confirmado como constante y con registro continuo.
- `estimated_current_capacity`: escala provisional estimada; puede orientar la revisión, pero no produce CONFORME por sí sola.
- `unavailable`: no existe información temporal utilizable.
