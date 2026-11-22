# Objetivos aprobados y evidencia producida

## OBJETIVO GENERAL

> Estudiar la implementación, calibración y validación de un modelo electroquímico–térmico reducido para una batería de ion-litio tipo 18650, considerando la identificación de parámetros a partir de bases de datos reportadas y mediciones experimentales propias, con el propósito de representar la respuesta eléctrica y térmica de la celda bajo diferentes condiciones de operación.

### Evidencia producida

- Implementación documentada del modelo agrupado SPM, generación de calor y dinámica térmica de un nodo.
- Calibración eléctrica y térmica en la base Panasonic NCR18650PF.
- Validaciones independientes 1C, UDDS y CC-CV, con limitaciones declaradas.
- Auditoría experimental y transferencia al caso TRC18650 de Uniandes.

## OBJETIVOS ESPECÍFICOS

### 1.

> Implementar un modelo electroquímico–térmico reducido, basado en el Single Particle Model (SPM) o en una variante afín, para representar variables relevantes de una celda 18650, tales como voltaje terminal, estado de carga, generación de calor y temperatura.

#### Evidencia producida

- Ecuaciones implementadas y auditadas en `thesis_data/final_model_equations.md`.
- Figuras 01, 07 y 10: representación OCP, respuesta térmica y selección de arquitectura.
- Salidas explícitas de voltaje, SOC, estados de electrodo, calor y temperatura sin recorte silencioso.

### 2.

> Formular una metodología de selección y calibración de parámetros del modelo, considerando parámetros eléctricos, electroquímicos y térmicos obtenidos de la literatura, especificaciones del fabricante, bases de datos reportadas o mediciones propias.

#### Evidencia producida

- Figura 14: flujo reproducible de auditoría, clasificación, transferencia y validación.
- Figuras 01, 02 y 07: identificación OCP, calibración HPPC condicional y calibración térmica.
- Registros trazables de procedencia, dominio de soporte, identificabilidad y rol de cada conjunto de datos.

### 3.

> Contrastar el desempeño del modelo implementado frente a bases de datos experimentales reportadas y mediciones propias bajo condiciones controladas de operación, utilizando perfiles de carga, descarga, pulsos de corriente o perfiles dinámicos disponibles.

#### Evidencia producida

- Figuras 02–09: contraste con pulsos HPPC, descarga 1C, UDDS y carga CC-CV.
- Figuras 11 y 12: mediciones propias, referencia experimental y ensayos afectados.
- Métricas RMSE, MAE, sesgo y R² reportadas según el rol de evidencia de cada ensayo.

### 4.

> Analizar la transferencia de los parámetros calibrados entre diferentes conjuntos de datos o celdas 18650, evaluando cuándo los parámetros identificados son suficientes para representar la respuesta del sistema y cuándo es necesario considerar un reajuste o adaptación parcial de los mismos.

#### Evidencia producida

- Figuras 02–05, 08 y 09: transferencia entre regiones de SOC, perfiles y dirección de corriente en Panasonic.
- Figura 13: transferencia directa Panasonic a TRC18650 y diagnóstico de incompatibilidad.
- Decisión jerárquica de adaptación: mapeo de capacidad/estado, referencia de voltaje compatible y evidencia dinámica adicional cuando exista base de tiempo medida.
