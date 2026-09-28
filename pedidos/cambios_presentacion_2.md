# CAMBIO — Presentación, segunda tanda

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.
Todo esto es sobre `R/13_presentacion.R` y sus dos salidas.

---

## 0. El hilo de la presentación

La charla tiene una tesis, y las diapositivas de resultados tienen que servirla:

> **Sola habría hecho boxplots. Con el agente pude hacer análisis que no hacía, y
> esos análisis explicaron cosas que el boxplot dejaba sin resolver.**

Por eso las diapositivas 6 y 7 son "el análisis convencional" y las 8 y 9 son "lo que
agregó explorar con el agente". Ese contraste es el argumento, no un detalle de
orden.

La presentación pasa de 11 a **12 diapositivas**.

---

## 1. Diapositiva 1 — Título

Título del trabajo:

> **Transportadores de nutrientes en el eje placenta–cerebro fetal en un modelo de
> activación inmune materna**

(Se acortó respecto de "Caracterización de los transportadores...": en una
diapositiva de título, el sustantivo solo funciona mejor que el nominalizado.)

Subtítulo, como está: reanálisis bioestadístico dirigido por un agente de IA.

---

## 2. Diapositiva 2 — "Antecedentes"

- **Título: "Antecedentes"**, no "El problema".
- **Los rayos tienen que apuntar al vientre**, no salir hacia afuera. Hoy se leen
  como algo que emana de la figura; tienen que leerse como algo que incide sobre la
  gestación. Ubicalos apuntando hacia el vientre y, si hace falta, con una punta de
  flecha sutil.
- **La frase de cierre va dentro de un recuadro**, con fondo suave y borde fino,
  alineada al ancho del contenido. Hoy queda suelta al pie.
- **Sacar los `--` del texto.** Aparecen en varias diapositivas como dos guiones
  literales. Reemplazalos por puntuación normal (coma, punto o dos puntos) o por un
  guion tipográfico. Revisá **todas** las diapositivas, no solo esta.

---

## 3. Diapositiva 3 — "Modelo experimental"

- **Título: "Modelo experimental"**, no "El experimento".
- La imagen `assets/modelo-experimental.png` va a ser reemplazada por una versión sin
  las anotaciones de qué se mide en cada tejido, para no repetir lo que ya dice el
  texto. Usá el mismo nombre de archivo; no cambies el código por esto.
- El texto de las tres mediciones se mantiene.

---

## 4. Diapositiva 4 — Lo que se le pidió al agente

Hoy muestra solo las decisiones fijas, y entonces la diapositiva 5 aparece
descolgada: no se entiende de dónde salieron 26 figuras si lo único que se mostró
fueron tres decisiones metodológicas.

Agregá al final de la diapositiva, después de "No fue «analiza mis datos»", una
línea breve que cierre el arco: sobre esas decisiones se pidió después **el análisis
convencional** (modelos por gen y tejido, con sus boxplots) y, encima, una
**exploración de los datos** que el análisis convencional no incluye (correlaciones
entre tejidos, co-expresión entre genes, dispersión).

---

## 5. Diapositiva 5 — Lo que salió de ahí

Sin cambios, salvo los `--`.

---

## 6. Diapositivas 6 y 7 — El análisis convencional

Encabezado de sección de estas dos: **"El análisis convencional"** en vez de "Lo que
produjo el agente".

**Diapositiva 6 — Placenta.** Hoy muestra solo pSTAT3. Pasa a mostrar:
- Boxplots de **tres genes** de placenta con efecto de tratamiento, recortados del
  panel completo o generados como subconjunto: `il6`, `glut3`, `slc38a2`.
- `acto1_pstat3.png`, más chico, al costado.
- La idea, en una línea: el LPS modifica la expresión en placenta y activa la vía
  IL-6/STAT3, pero el boxplot no muestra interacción sexo × tratamiento en ningún
  gen.

**Diapositiva 7 — Cerebro fetal.** Hoy muestra el panel completo y no se lee nada.
Pasa a mostrar **tres genes elegidos**: `glut1`, `slc38a2`, `fatp1`, que son los que
tienen post hoc significativo.
- La idea: acá sí hay interacción, y el post hoc la ubica en las hembras.

Para los dos: si recortar el panel existente no da buena calidad, generá una figura
nueva con solo esos genes, reusando las funciones de `07_figuras_acto1.R`, sin
recalcular nada. Anotalo en `procedencia.csv`.

---

## 7. Diapositivas 8 y 9 — Lo que agregó explorar con el agente

Encabezado de sección de estas dos: **"Lo que agregó explorar con el agente"**.
Quitar el rótulo "la diapositiva clave" del encabezado: si una diapositiva necesita
anunciarse como importante, no lo es.

**Diapositiva 8 — Co-expresión entre genes (nueva).** Un diagrama triangular por
tejido o por sexo, el que se lea mejor proyectado. La idea: con el boxplot cada gen
se mira por separado; el diagrama triangular muestra cómo se mueven los genes entre
sí, y eso es un análisis que no estaba en el plan original.

**Diapositiva 9 — La dispersión.** Esta es la que hoy no se entiende, y el problema
es el texto, no la figura. Tiene que decir el razonamiento completo, en tres pasos:

1. En cerebro, siete genes mostraron interacción sexo × tratamiento en el boxplot.
2. En cuatro de ellos, **ninguna comparación entre grupos la explicaba**: el post hoc
   no encontraba entre qué pares estaba la diferencia. Mirando el boxplot quedaba
   "hay algo sexo-dependiente, pero no se ve dónde".
3. El análisis de dispersión lo explicó: la diferencia no está en dónde se ubican los
   datos sino en **cuánto se dispersan**. El LPS **compacta** la expresión en hembras
   y la **dispersa** en machos. Es un dimorfismo que el boxplot no puede mostrar.

Figura: `acto2_dispersion_sd.png` con los paneles por sexo. Si hay una figura que
muestre mejor el contraste (por ejemplo la diagonal del diagrama triangular de
cerebro de hembras, donde Control es ancho y LPS es un pico angosto), usá esa y
decime cuál elegiste.

**Nada de esto puede incluir números reales en la versión pública** (ver punto 10).

---

## 8. Diapositivas 10 y 11 — Errores

Sin cambios de contenido.

**Diapositiva 10**: mover el texto a una columna lateral y sumar una figura que
muestre el desacople, si hay alguna que lo haga sin exponer números inéditos.
Si no la hay, dejala como está y avisame.

---

## 9. Diapositiva 12 — Cierre

Dejala como está por ahora. Se reescribe cuando termine la segunda prueba de
reproducción con el agente externo.

---

## 10. Reglas que siguen valiendo

- Ningún número que dependa de los datos se escribe a mano: todo se lee de
  `outputs/tables/`.
- La versión pública no lleva ninguna afirmación biológica ni ningún número real, y
  las figuras que muestre son las generadas con datos sintéticos, acompañadas del
  aviso. **Confirmame que en la versión pública las figuras efectivamente aparecen**:
  en el HTML actual, las diapositivas 6, 7 y 8 parecen mostrar solo el aviso.
- Toda figura nueva tiene su fila en `procedencia.csv`.
- La verificación de que la versión pública no contiene frases interpretativas tiene
  que seguir pasando, y cubrir las diapositivas nuevas.

---

## 11. Cierre

Generá las dos versiones, decime tamaño y ruta de cada una, y confirmame que
`git check-ignore` sigue excluyendo `outputs/presentacion_real/`. Un commit por
bloque. Actualizá `ESTADO.md` y `procedencia.csv`.

Antes de escribir código decime en dos líneas qué vas a hacer, y avisame si alguna
diapositiva no se puede armar con lo que hay en el repo.
