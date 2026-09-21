# CAMBIO EN EL INFORME — Conclusiones por sección, síntesis y conclusión revisada

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

---

## 0. Alcance: SOLO REDACCIÓN

Este pedido **no agrega ningún análisis**. No se corre ningún test nuevo, no se
recalcula nada, no se modifica ningún script de análisis (`02` a `11`). Se toca
únicamente `12_informe` (R y Python) para insertar texto.

**Regla de procedencia para cada número del texto:** los textos de abajo son
borradores escritos a partir del informe. Antes de insertar cada uno, verificá
**cada número** contra el CSV correspondiente en `outputs/tables/`. Si un número
no coincide, **no lo corrijas en silencio**: usá el valor del CSV y listame al
final qué números cambiaste y de qué tabla salió el correcto. Si una afirmación
no se puede verificar contra una tabla existente, redactala como observación
visual de la figura, sin afirmar más de lo que se ve.

Ningún número del texto se escribe a mano en el script: se lee del CSV al
generar el informe, igual que el resto de las secciones. Así, si los datos
cambian, el texto cambia con ellos.

---

## 1. Conclusión breve al final de cada sección

Agregar al final de cada sección un bloque titulado **"Conclusión de la sección"**,
de 2 a 5 oraciones, con el siguiente contenido.

### Acto 1.1 — ELISA de IL-6

> El LPS indujo una respuesta inflamatoria sistémica: la IL-6 fue detectable en
> 9/9 madres LPS frente a 1/5 control (Fisher p = 0.005), lo que valida el
> modelo. En líquido amniótico ningún contraste alcanzó significancia
> (Peto-Peto ♀ p = 0.30; ♂ p = 0.33), aunque la mediana de los valores
> detectados fue mayor bajo LPS en ambos sexos. Con 5 sacos control por sexo,
> la ausencia de significancia no permite concluir que la IL-6 no llegue al
> compartimento fetal.

### Acto 1.2 — Cuantificación relativa

Sin conclusión: es una sección de método.

### Acto 1.3 — Modelos de expresión

> **Placenta.** Ningún gen mostró interacción SEXO×TTO. El LPS modificó la
> expresión de il6, glut3, slc38a2, glut1 y CD36 de forma equivalente en ambos
> sexos (efecto principal de tratamiento).
>
> **Cerebro fetal.** Siete genes mostraron interacción SEXO×TTO significativa,
> todas robustas a la corrección BH. En glut1, slc38a2 y fatp1 el post hoc
> localiza el efecto en hembras: ♀Control difiere de ♀LPS y ♀LPS difiere de
> ♂LPS, sin cambios en machos. En fatcd36, fatp4, gp130 y slc38a1 la
> interacción no se explica por ninguna comparación de medias (todos los p de
> Holm ≥ 0.07). La sección 2.3 muestra que en esos cuatro genes el efecto está
> en la dispersión y no en la media.

### Acto 1.4 — pSTAT3

> El LPS aumentó la abundancia de fosfo-STAT3 en placenta solo en hembras
> (♀Control 1.98 → ♀LPS 4.81; Holm p = 1.3×10⁻⁷). En machos no cambió
> (2.70 → 2.75; p = 0.88). Es el mismo patrón que los genes con efecto en media
> en cerebro. Como los transportadores placentarios responden igual en ambos
> sexos, sus cambios no dependen de la activación de STAT3, o los machos los
> alcanzan por otra vía. Recordar D9: se mide abundancia de fosfo-STAT3, no
> fracción fosforilada.

### Acto 2.1–2.2 — Correlación placenta↔cerebro y co-expresión

> Descriptivamente, en hembras control la correlación placenta–cerebro es alta
> en varios genes y cae cerca de cero con LPS; en machos no hay correlación en
> ningún grupo. Estas figuras describen y no testean (prohibición 4). En los
> diagramas triangulares de cerebro de hembras, la distribución ♀Control es
> ancha con una cola hacia valores bajos, y la ♀LPS es un pico angosto: la
> reducción de dispersión de la sección 2.3 es visible directamente.

### Acto 2.3–2.4 — Dispersión y test de Δρ

> Ningún Δρ Control vs LPS es significativo, ni agrupando sexos ni dentro de
> cada sexo (0 de 28; mínimo fatcd36 en hembras, p = 0.075). El test de
> interacción SEXO×TTO sobre la dispersión es el resultado positivo del Acto 2:
> en cerebro, cinco genes sobreviven a BH (gp130, fatcd36, fatp4, fatp1,
> slc38a2; slc38a1 en tendencia), con un patrón cruzado: el LPS reduce la
> dispersión en hembras y la aumenta en machos. En placenta ninguno sobrevive a
> BH. Como los transportadores varían mayormente juntos (el eigengene explica
> el 82 % de la varianza en cerebro; sección 2.6), estos genes no son efectos
> independientes: reflejan un mismo patrón compartido por el conjunto de
> transportadores.
>
> La figura de SD por ítem agrupa los sexos y por eso no muestra este efecto:
> los cambios opuestos de hembras y machos se cancelan al promediarlos. Es un
> ejemplo directo de lo que oculta agrupar los sexos.

### Acto 2.5 — Simulación de restricción de rango

> La caída de correlación observada en hembras es compatible con la
> compactación de la expresión bajo LPS: cuando el rango de una variable se
> reduce, la correlación cae aunque la relación biológica no haya cambiado. Las
> excepciones (fatcd36 y score compuesto en hembras) quedan en el límite del
> intervalo y se toman como pista, no como hallazgo. La pérdida aparente de
> acoplamiento placenta–cerebro en hembras no es un efecto independiente: es la
> reducción de dispersión vista desde otro ángulo.

### Acto 2.6 — Sensibilidad

> El eigengene (PC1) explica el 76 % de la varianza de los transportadores en
> placenta y el 82 % en cerebro, con cargas similares para los siete genes: los
> transportadores varían mayormente juntos, como un único eje por feto.
> Reemplazar el score compuesto por el eigengene no cambia ninguna conclusión.
> Excluir el feto más influyente de cada ítem tampoco cambia veredictos, pero
> aproximadamente duplica los p de los ítems que estaban cerca del umbral
> (fatcd36 0.107 → 0.211; score 0.079 → 0.175): esas señales son sensibles a la
> exclusión de un solo feto.

---

## 2. Sección nueva: Síntesis del eje madre → placenta → cerebro

Al final del Acto 1, antes del Acto 2. Texto:

> **Madre.** El LPS produjo una respuesta inflamatoria sistémica inequívoca
> (IL-6 sérica detectable en 9/9 madres tratadas frente a 1/5 control).
>
> **Líquido amniótico.** Sin resultado concluyente: con 5 sacos control por
> sexo no se detectó diferencia, lo que no descarta que la IL-6 llegue al
> compartimento fetal.
>
> **Placenta.** Dos respuestas que no coinciden. La señalización IL-6/STAT3 se
> activa solo en placentas de fetos hembra. La expresión de transportadores de
> nutrientes, en cambio, se modifica en ambos sexos por igual.
>
> **Cerebro fetal.** La respuesta depende del sexo del feto. En hembras, el
> LPS modifica la expresión de glut1, slc38a2 y fatp1, con el mismo patrón que
> pSTAT3 en placenta; en machos, esos genes no cambian.
>
> **Lectura del Acto 1.** El LPS materno induce en la descendencia hembra una
> respuesta coherente a lo largo del eje: activación de STAT3 en placenta y
> cambio de la expresión de transportadores en cerebro. En la descendencia
> macho, la respuesta en cerebro no es detectable como cambio de nivel.

---

## 3. Sección nueva: Conclusión revisada (Acto 1 frente a Acto 2)

Al final del Acto 2, antes de Limitaciones. Una tabla y un párrafo.

| Afirmación del Acto 1 | Qué la puso a prueba en el Acto 2 | Resultado |
|---|---|---|
| Siete genes de cerebro responden con dimorfismo sexual | Test de interacción SEXO×TTO sobre la dispersión | Se precisa: 3 genes con desplazamiento de la media en hembras (glut1, slc38a2, fatp1); en los otros 4 (fatcd36, fatp4, gp130, slc38a1) lo que cambia es la variabilidad, como parte de un patrón compartido por los transportadores |
| (Lectura intuitiva de las figuras) El LPS desacopla placenta y cerebro en hembras | Test formal de Δρ + simulación de restricción de rango | Se descarta: ningún Δρ significativo; la caída de correlación se explica por la compactación de la expresión |
| La respuesta resumida por el score compuesto es robusta | Eigengene PC1 y exclusión del feto extremo | Se sostiene: el eigengene no cambia conclusiones; las señales cercanas al umbral son sensibles a la exclusión de un solo feto |

> **Conclusión revisada.** En hembras, el LPS materno produce una respuesta
> direccional y homogénea: activa STAT3 en placenta, desplaza la expresión
> cerebral de glut1, slc38a2 y fatp1, y reduce la variabilidad entre
> individuos. En machos no hay respuesta direccional, pero aumenta la
> heterogeneidad: cada individuo responde distinto. La aparente pérdida de
> acoplamiento placenta–cerebro en hembras no es un segundo hallazgo, sino la
> misma reducción de dispersión vista a través de la correlación.

---

## 4. Actualizar el Resumen

El Resumen del inicio es anterior a la sesión 17 y no menciona el test de
interacción sobre la dispersión, que es el resultado principal del Acto 2.
Reescribirlo de modo que refleje la conclusión revisada de la sección 3: el
párrafo de "Coordinación placenta↔cerebro" tiene que decir que no hay cambio de
coordinación **y** que sí hay un efecto sexo-dependiente sobre la dispersión en
cerebro.

---

## 5. OPCIONAL — Correcciones de consistencia

No son análisis nuevos: son figuras que no reflejan resultados que ya están
calculados. Aplicalas solo si confirmo.

a) **Las figuras de 2.3, 2.4 y 2.5 muestran solo el estrato agrupado.** Las
tablas tienen la columna `ESTRATO` desde la sesión 17, pero
`acto2_dispersion_sd.png`, `acto2_test_delta_rho.png` y
`acto2_simulacion_delta_rho.png` siguen mostrando solo `AMBOS_SEXOS`. Agregar
los estratos `HEMBRA` y `MACHO` a esas tres figuras.

b) **`il6R` sigue apareciendo en el Acto 2** (figuras de SD, Δρ, simulación y
sensibilidad, y filas de las tablas), aunque se pidió excluirlo. La figura de
correlación global dice en el título que está excluido, así que hay una
inconsistencia visible.

c) **El triángulo superior de los SPLOM por sexo muestra un único rho** que
junta Control y LPS del mismo sexo, en vez de un rho por tratamiento como en la
versión anterior.

---

## 6. Cierre

Regenerá el informe HTML y el PDF. Un commit para los puntos 1 a 4 y, si lo
confirmo, otro para el 5. Actualizá `ESTADO.md`.

Antes de escribir código decime en dos líneas qué vas a hacer, y al terminar
listame los números que no coincidían con los CSV.
