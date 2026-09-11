# CAMBIO EN T7 — Correlaciones y co-expresión separadas por sexo

Pegar como mensaje al agente. Aplica a `R/08_acto2_correlaciones.R` **y** a su gemelo
`python/08_acto2_correlaciones.py`. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

---

## 0. Qué NO cambia

- Las decisiones D1–D12 siguen fijas.
- Sigue valiendo la prohibición 4: **no se compara rho entre estratos** (ni entre
  tratamientos, ni entre sexos, ni entre paneles). Las figuras describen; el test
  formal es el de pendientes de T8.
- Sigue valiendo la prohibición 5: un cambio de rho no se interpreta como cambio
  de coordinación sin descartar restricción de rango.
- Emparejamiento por `FETO`: un par entra solo si el gen está detectado en
  **ambos** tejidos.
- Todo se calcula sobre `-ΔΔCt`. Los ejes se muestran como `2^(-ΔΔCt)` en escala
  log2 (D2).
- **Solo Spearman.** No agregar Pearson en ninguna figura ni tabla.

---

## 1. Extender los estratos

`ESTRATOS` pasa de 3 a 7:

```
GLOBAL, CONTROL, LPS, HEMBRA_CONTROL, HEMBRA_LPS, MACHO_CONTROL, MACHO_LPS
```

Los tres actuales **se conservan**: la vista sin separar por sexo sigue siendo
parte del informe.

`cargar()` hoy no extrae el sexo. Agregarlo: usar la columna `SEXO` de
`qpcr_cuantificacion_long.tsv` si existe; si no existe, derivarla de `GRUPO`
(`HEMBRA_CONTROL` → `HEMBRA`). Verificá cuál de los dos casos es el real y
dejalo registrado.

`pares()` tiene que filtrar por sexo y tratamiento según el estrato pedido.

**Piso de n.** `PISO_PAR = 5` sigue igual. Con los estratos por sexo, n ≤ 9 por
celda, y en genes con muchos no detectados (`il6R` en cerebro) varias celdas van
a quedar por debajo. En esos casos: se grafican los puntos igual, pero **no** se
dibuja la línea de tendencia, no se reporta rho ni p, y la leyenda dice
`n = X (sin rho)`. No fuerces un rho con n < 5.

---

## 2. Figura nueva: una por gen, dos paneles por sexo

**Una figura por cada item** de la lista nueva de `ITEMS`:

```
fatcd36, fatp1, fatp4, glut1, glut3, slc38a1, slc38a2, gp130, score_compuesto
```

Nueve figuras. Cambios respecto de lo que hay hoy:

- `il6` sigue fuera: no tiene `-ΔΔCt` en cerebro (calibrador ♀Control 0/9, D7).
- **`il6R` sale** de las correlaciones por pedido explícito: en cerebro tiene
  3/7/3/3 detectados y al estratificar por sexo casi todas las celdas caen bajo
  el piso de 5 pares. `il6R` **se conserva** en los boxplots del Acto 1 y en las
  tablas de modelos: este cambio es solo para las figuras de correlación.

Estructura de cada figura:

- Título: el nombre del gen, en itálica.
- Subtítulo: `Placenta–brain correlation at E15`.
- **Dos paneles lado a lado**: `Females` (izquierda) y `Males` (derecha), como
  facetas con recuadro visible.
- Eje X: `Placenta E15 — <gen>/rsp29 relative expression`.
- Eje Y: `Brain E15 — <gen>/rsp29 relative expression`.
- Ambos ejes en escala log2, mostrando `2^(-ΔΔCt)`.
- Dentro de cada panel, los dos tratamientos superpuestos: **Control** y **LPS**,
  con forma y color distintos por grupo.
- Línea de tendencia lineal por tratamiento, con **banda de dispersión
  transparente al 95%**.
- Leyenda al pie, con una entrada por combinación sexo × tratamiento (4 en
  total), cada una mostrando: `Spearman rho = X.XX ; p = X.XXX` y `n = X`.
  **Quitar Pearson.**
- Nota al pie: `Axes: 2^-ΔΔCt (log2 display) | Spearman on -ΔΔCt | Linear fit
  with 95% CI | Same fetus pairing`.

Nombre de archivo: `outputs/figures/acto2_corr_placenta_cerebro_<item>.png`

La figura global actual (`acto2_corr_*` con facetas por gen, sin separar sexo)
**se conserva** tal como está.

---

## 3. SPLOM por sexo

Los dos SPLOM actuales (`acto2_coexpresion_SPLOM_PLACENTA_E15.png` y
`..._BRAIN_E15.png`) **se conservan sin cambios**, incluida la diagonal con
densidad suavizada y sombreada por tratamiento, que es el estilo pedido.

**Agregar cuatro SPLOM nuevos**, uno por tejido × sexo:

```
acto2_coexpresion_SPLOM_PLACENTA_E15_HEMBRA.png
acto2_coexpresion_SPLOM_PLACENTA_E15_MACHO.png
acto2_coexpresion_SPLOM_BRAIN_E15_HEMBRA.png
acto2_coexpresion_SPLOM_BRAIN_E15_MACHO.png
```

**Además, el conjunto de genes del SPLOM cambia y es distinto por tejido.** Hoy
son los 7 transportadores en los dos. Pasa a ser:

- **Placenta E15**: los 7 transportadores + `il6` + `gp130` (9 variables).
- **Cerebro E15**: los 7 transportadores + `gp130` (8 variables).

`il6` entra solo en placenta porque ahí sí es cuantificable; en cerebro no lo es
(D7). `il6R` no entra en ningún SPLOM.

Esto vale para las seis figuras: las dos que ya existen (ambos sexos juntos) y
las cuatro nuevas por sexo. Como los dos tejidos ya no tienen el mismo conjunto
de variables, declaralo en el pie de figura y en el informe: las matrices de
placenta y cerebro **no son comparables celda por celda**.

Mismo formato que los actuales: diagonal con densidad KDE de Control vs LPS
superpuestas, triángulo inferior con scatter de individuos, triángulo superior
con Spearman. Dentro de cada figura por sexo el n baja a ~9 por tratamiento, así
que **poné el n en cada panel del triángulo superior** junto al rho, y respetá el
piso de 5.

En el triángulo superior, `GGally::wrap("cor", method = "spearman")` etiqueta el
valor como `Corr:`, lo cual es ambiguo. Cambiar la etiqueta a `rho:` para que no
se confunda con Pearson.

`tabla_coexpresion()` tiene que ganar una columna `SEXO` con los valores
`AMBOS`, `HEMBRA`, `MACHO`, y recorrer el conjunto de genes que corresponda a
cada tejido, de modo que la tabla acompañe a las seis figuras. La constante
`TRANSP` deja de alcanzar: definí en `00_config` dos listas, `GENES_SPLOM_PLACENTA`
y `GENES_SPLOM_BRAIN`, para que el conjunto no quede escrito a mano dentro del
script de figuras.

---

## 4. Verificaciones a actualizar

La verificación `acto2_estratos` afirma hoy que los estratos son exactamente
`GLOBAL;CONTROL;LPS` y va a fallar. Reescribirla para que compruebe los siete
estratos nuevos.

Agregar estas verificaciones:

1. Cada panel sexo × tratamiento de la figura por gen tiene **n ≤ 9**, y el n
   reportado en la leyenda coincide con la cantidad de puntos dibujados.
2. Los estratos por sexo particionan a los globales: `n(HEMBRA_CONTROL) +
   n(MACHO_CONTROL) = n(CONTROL)` para cada item. Si no cierra, hay un error de
   filtrado.
3. Ninguna figura ni tabla del Acto 2 contiene valores de Pearson.
4. Los items con n < 5 en algún estrato aparecen sin rho y sin línea de
   tendencia, no con un rho calculado sobre n chico.
5. `il6R` no aparece en ninguna figura ni tabla de correlación del Acto 2, y sí
   sigue apareciendo en las figuras y tablas del Acto 1.
6. El SPLOM de placenta tiene 9 variables y el de cerebro 8, con los conjuntos
   exactos especificados arriba.

---

## 5. Procedencia

Una fila nueva en `procedencia.csv` por cada figura nueva (9 por gen/score + 4
SPLOM), y actualizá las filas de las dos figuras SPLOM que ya existían, porque
les cambia el conjunto de variables. En la columna de decisión, para las figuras por sexo:

> Estratificación por sexo pedida explícitamente. Alternativa descartada: dejar
> solo los estratos GLOBAL/CONTROL/LPS (n mayor, IC más angostos), descartada
> porque el dimorfismo sexual es la pregunta del proyecto y los estratos
> agregados lo promedian. Limitación declarada: con n ≤ 9 los IC de rho son
> anchos y la comparación entre paneles no está testeada.

Agregar también en `analisis_descartados.md` la corrección del criterio anterior:
el script había fijado `GLOBAL/CONTROL/LPS` argumentando que las celdas
SEXO×TTO darían IC inútiles; esa decisión se revierte por pedido explícito, con
la limitación declarada en el informe en lugar de omitir el análisis.

Y una entrada nueva: `il6R` se excluye de las correlaciones del Acto 2 por
detección insuficiente en cerebro (3/7/3/3), que al estratificar por sexo deja
casi todas las celdas bajo el piso de 5 pares. Se conserva en el Acto 1.

---

## 6. Orden de trabajo

1. Cambios en R, correr y revisar las figuras.
2. Espejar en Python.
3. Correr la comparación R↔Python y `99_verificar`.
4. Regenerar el informe HTML.
5. Un commit por paso, con mensaje descriptivo.

Antes de escribir código, decime en dos líneas qué vas a hacer y avisame si algo
de esto entra en conflicto con `AGENTS.md`.
