# Pie de figura propuesto

Definición de las fases de evaluación sobre un pulso HPPC real. La corriente se presenta positiva en descarga.

- **Respuesta al inicio del pulso:** error asociado al cambio inicial de voltaje después de aplicar el escalón de corriente.
- **Durante el pulso:** error mientras la corriente permanece aplicada.
- **Relajación:** error durante la recuperación del voltaje después de retirar la corriente.
- **RMSE balanceado:** raíz cuadrada del promedio de los tres errores cuadráticos medios de fase. Para cada fase, primero se obtiene el error cuadrático medio de cada pulso y luego se promedian los pulsos con el mismo peso; así se evita que una fase o un pulso domine solamente por contener más muestras.
