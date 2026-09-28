# AGENTS.md — Reanálisis reproducible MIA-LPS (placenta E15 / cerebro fetal E15)

> **Todo agente que abra este proyecto debe leer este archivo y `ESTADO.md` completos
> antes de ejecutar nada.** Este archivo fija las convenciones y decisiones; `ESTADO.md`
> es la bitácora de avance entre sesiones. No se reabren las decisiones D1–D13.

---

## 1. Objetivo verificable y checklist de finalización

**El proyecto está terminado cuando**, partiendo de un clon limpio del repositorio y sin
acceso a los datos crudos, un tercero puede correr un único comando de arranque
documentado en el README (`.\run_all.ps1`) y obtener, sobre **datos sintéticos**, todas
las salidas listadas abajo, y `R/99_verificar.R` (o `python/99_verificar.py`) imprime
`TODAS LAS VERIFICACIONES PASARON`.

Checklist (el script `99_verificar` chequea existencia **y no-vacuidad** de cada ítem):

- [x] `AGENTS.md` y `CLAUDE.md` en la raíz, con convenciones y decisiones fijas.
- [x] `README.md` con instalación del entorno y ejecución completa desde cero.
- [x] `ESTADO.md` con el registro de avance entre sesiones.
- [x] Scripts numerados en `R/` y en `python/`, funcionalmente equivalentes.
- [x] Generador de datos sintéticos que reproduce la estructura exacta de los tres archivos de entrada.
- [x] Figuras del ACTO 1: boxplots de expresión (por gen × tejido), boxplot de pSTAT3, figuras de ELISA.
- [x] Figuras del ACTO 2: correlaciones placenta–cerebro por gen, pair plots, figuras de dispersión, figura de la simulación.
- [x] `outputs/tables/procedencia.csv` — una fila por figura y por tabla.
- [x] `outputs/tables/verificaciones.csv` — una fila por resultado principal.
- [x] `outputs/tables/analisis_descartados.md` — qué se probó, por qué no funcionó, qué se hizo en su lugar.
- [x] `outputs/tables/comparacion_R_python.csv` — concordancia numérica entre ambas implementaciones.
- [x] `docs/informe.html` — informe completo autocontenido.
- [x] `docs/informe.pdf` — versión imprimible del informe.
- [x] `logs/corrida_<fecha>.txt` — fecha, entorno, versiones de paquetes, cuántas verificaciones pasaron.

Estado de avance del checklist: ver la tabla de la Sección 6 (plan de tareas).

---

## 2. Rutas y parámetros fijos

| Parámetro | Valor |
|---|---|
| `RAIZ_REPO` | `C:\Users\Usuario\Documents\Doctorado\Resultados\LPS 100\qPCR\MIA_LPS_reanalisis` |
| `ARCHIVO_QPCR` | `data/raw/Raw data CTs.xlsx` (copiado del crudo `qPCR\Raw data CTs.xlsx`) |
| `HOJA_QPCR` | `Sheet1` |
| `ARCHIVO_ELISA` | `data/raw/ELISA IL6 2026 Dosis 100.xlsx` |
| `HOJA_ELISA` | `Sueros y LA` |
| `ARCHIVO_PSTAT3` | `data/raw/pstat3 placenta.xlsx` |
| `HOJA_PSTAT3` | `Sheet1` |
| `CARPETA_SALIDA` | `outputs/` (relativa a `RAIZ_REPO`; **nunca** rutas absolutas en el código) |
| `SEMILLA` | `20260101` |

**Regla de acceso a datos:** el único acceso a datos crudos dentro de los scripts es
`data/raw/`. Nunca leer los originales desde su ubicación fuera del repo. `data/raw/`
(contenido) va en `.gitignore`: los datos son inéditos y no se versionan.

### 2.1 Notas de decodificación de los archivos crudos (verificadas al inicio)

- **qPCR (`Raw data CTs.xlsx` / `Sheet1`)**: formato largo, columnas
  `MADRE, NOMINACION, FETO, GRUPO, SEXO, TTO, TEJIDO, rsp29, GEN, CT_CRUDO`.
  Tejidos presentes: `PLACENTA_E15`, `BRAIN_E15`, `BRAIN_P1`.
  **`BRAIN_P1` se EXCLUYE en `02_ingesta_qc`** (fuera del alcance de la Sección 3;
  se analiza en otro informe) y se registra en `analisis_descartados.md`.
  E15 = 36 fetos exactos. `CT_CRUDO`: 129 celdas `== 40` + 30 celdas ya vacías (NA).
  Ambas se tratan igual: **no detectado → NA** (ver D3). `rsp29` no tiene ningún 40 ni NA.
  `GRUPO` viene con espaciado inconsistente (`"♀ LPS"` vs `"♀Control"`); `SEXO`/`TTO`
  con capitalización distinta entre archivos. La ingesta **normaliza** derivando los 4
  grupos de `SEXO` + `TTO`. El generador sintético reproduce la inconsistencia.
- **ELISA (`Sueros y LA`)**: fila de encabezado en índice 1, datos desde índice 2.
  Columnas: `0`=tratamiento(+sexo en LA), `1`=`Abs 450`, `2`=`Conc` (medición real, con
  signo), `3`=`IL-6` (= `Conc` con negativos pisados a 0 — **NO usar esta columna**),
  `4`=ID (madre en MS, feto en LA), `5`=`TEJIDO` (`Suero materno` / `Líquido amniótico`).
  Columnas 6–13 son celdas sueltas de scratch: se ignoran. La ingesta parte la hoja en
  dos tablas largas por la columna 5. MS: n=14 madres (4 con `Conc`<0). LA: n=28 sacos
  (10 con `Conc`<0). **Los 4 negativos de MS se tratan como censura a izquierda igual
  que los 10 de LA** (D10 es general). `LOD` registrado = `Conc = 0` (el blanco);
  alternativa documentada = menor estándar de la hoja `CURVA IL6`.
- **pSTAT3 (`Sheet1`)**: columnas `MEMBRANA, MADRE, NOMINACION, FETO, GRUPO, SEXO, TTO,
  TEJIDO, PSTAT3`. 36 filas, 3 membranas × 4 grupos × 3 réplicas, balanceado.

### 2.2 Discrepancias conocidas de n (resueltas explícitamente)

Diseño: 18 madres, 36 fetos. ELISA: 14 madres en suero, 28 sacos en líquido amniótico.
Es esperable (no todo se pudo colectar). **No se completa, no se imputa, no se asume.**
El reporte de QC (`02_ingesta_qc`) informa: n real por archivo × grupo × sexo; qué
madres/fetos faltan en cada bloque; y deja constancia de que el ELISA se analiza con su
n propio, sin forzar al n=36 del diseño de qPCR.

---

## 3. Estructura de carpetas y convención de nombres

```
MIA_LPS_reanalisis/
├─ AGENTS.md  CLAUDE.md  README.md  ESTADO.md  .gitignore
├─ requirements.txt        <- entorno Python (pip)
├─ renv.lock               <- entorno R (renv::snapshot); se genera al cerrar T1
├─ run_all.ps1             <- corre todo en orden (PowerShell)
├─ data/
│  ├─ raw/         <- 3 Excel crudos. CONTENIDO IGNORADO por git (.gitkeep versionado)
│  ├─ synthetic/   <- datos sintéticos. VERSIONADO
│  └─ processed/   <- intermedios regenerables. CONTENIDO IGNORADO por git
├─ R/    00_config.R 01_generar_sinteticos.R 02_ingesta_qc.R 03_elisa.R
│        04_qpcr_cuantificacion.R 05_qpcr_modelos.R 06_pstat3.R 07_figuras_acto1.R
│        08_acto2_correlaciones.R 09_acto2_dispersion.R 10_acto2_simulacion.R
│        11_sensibilidad.R 98_comparacion.R 12_informe.R 99_verificar.R
│        13_presentacion.R   <- SOLO R, ver excepcion mas abajo
│        14_informe_breve.R  <- SOLO R, ver cuarta excepcion mas abajo
├─ python/  <- mismos números y nombres, extensión .py (SIN 13_presentacion.py ni
│             14_informe_breve.py)
├─ assets/  <- ilustracion-mia.png, modelo-experimental.png, epidemiologia.png
│             (autoria propia, VERSIONADO; epidemiologia.png no se usa, ver 13_presentacion)
├─ outputs/
│  ├─ figures/       <- <fig>_<detalle>.png (300 dpi)
│  ├─ tables/        <- procedencia.csv, verificaciones.csv, comparacion_R_python.csv,
│  │                    comparacion_reporte.md, analisis_descartados.md, y CSV de resultados
│  │  ├─ R/          <- salidas numéricas de la implementación R (nombres idénticos a python/)
│  │  └─ python/     <- salidas numéricas de la implementación Python
│  ├─ intermediate/  <- CONTENIDO IGNORADO por git
│  ├─ presentacion_real/  <- 13_presentacion sobre datos reales. NUNCA se versiona
│  └─ informe_breve_real/ <- 14_informe_breve sobre datos reales. NUNCA se versiona
├─ docs/   informe.html  informe.pdf   (listo para GitHub Pages, CONTENIDO IGNORADO)
│          index.html  presentacion.pdf  referencias.md   <- 13_presentacion sobre datos
│          sinteticos: presentacion PUBLICA, estos 3 SI se versionan (excepcion a docs/*)
│          informe_breve.html  informe_breve.pdf  <- 14_informe_breve sobre datos
│          sinteticos: informe breve PUBLICO, estos 2 SI se versionan (misma excepcion)
└─ logs/   corrida_<AAAA-MM-DD>.txt
```

**Excepcion a "todo se implementa dos veces" -- `13_presentacion` (solo R):** es capa de
presentacion (`pedidos/cambios_presentacion.md`), no analisis -- arma la pagina/PDF de
exposicion a partir de lo que ya escribieron 02..12, sin calcular ningun resultado nuevo
que comparar entre lenguajes. Por eso no tiene gemelo en `python/`. Genera DOS versiones
con el mismo codigo, decidiendo el destino por `fuente_datos()` (igual mecanismo que el
aviso sintetico de `12_informe`): sintetica -> `docs/index.html` (PUBLICA, se versiona);
real -> `outputs/presentacion_real/` (NUNCA se versiona, cubierto por `.gitignore`). El
grafico epidemiologico de contexto (diapositiva 2) es el unico grafico del repo que no
sale de `outputs/figures/`: no es un resultado del proyecto, es contexto bibliografico
(`pedidos/referencias_epidemiologia.md` -> `docs/referencias.md`, citas copiadas tal cual).

**Segunda excepcion (pedido 2, `pedidos/cambios_presentacion_2.md`):** dos diapositivas
(placenta y cerebro fetal) necesitan boxplots de un SUBCONJUNTO de 3 genes por tejido, no
el panel completo. `13_presentacion.R` reusa `panel_gen()` y sus dependencias de
`07_figuras_acto1.R` (pedido explicito: "sin recalcular nada") sourceando ese archivo
COMPLETO pero con `local = <environment nuevo>`, para no chocar con los nombres propios de
`13_presentacion` (07 redefine `ESTE_SCRIPT`, `registrar_procedencia`, etc. con firmas
distintas). Las dos figuras resultantes (`outputs/figures/acto1_expresion_{PLACENTA,
BRAIN}_E15_subset3.png`) llevan `tipo = "figura_presentacion"` en `procedencia.csv`, no
`"figura"`: ese segundo valor las deja fuera del chequeo de `12_informe` que exige que
toda fila `tipo == "figura"` este embebida en `docs/informe.html` -- estas dos son solo
para la presentacion, nunca para el informe.

**Tercera figura reusada (pedido del usuario en chat, sin archivo de pedido):**
`13_presentacion.R` tambien regenera pSTAT3 para la diapositiva 6, reusando
`figura_pstat3()` del mismo `.ENV07()` -- guardada como
`outputs/figures/acto1_pstat3_presentacion.png` (NO `acto1_pstat3.png`, el que usa la
corrida principal), para no pisar el archivo que lee `12_informe`. Motivo: en la version
PUBLICA (sintetica), los brackets de TENDENCIA (D11, 0.05<=p<0.1, estilo "punteada") no
imprimen el numero de p -- se sobreescribe `d11_texto` dentro de `.ENV07()`
(`aplicar_override_tendencia()`, nunca en `07_figuras_acto1.R` mismo) para que, cuando el
estilo sea "punteada", el texto quede vacio; el bracket punteado se sigue dibujando. Los
brackets SOLIDOS (`*`, `**`, `***`) se muestran igual en las dos versiones -- son datos
simulados en la publica, no exponen nada. La version REAL nunca aplica el override.

**Cuarta excepcion (`pedidos/pedido_informe_breve_final.md`, reemplaza a
`pedidos/pedido_informe_breve.md`):** `14_informe_breve.R` (SOLO R, mismo criterio que
`13_presentacion.R`: capa de presentacion, no analisis, sin gemelo en `python/`) arma un
informe breve de maximo 5 paginas para lectura rapida. Sourcea `13_presentacion.R`
COMPLETO con `local = <environment nuevo>` (`.ENV13()`), que a su vez tiene su propio
`.ENV07()` -- por eso este informe hereda, sin recalcular nada, el subconjunto de
boxplots de cerebro de la segunda excepcion y, a traves de el, `panel_deteccion()` de
`07_figuras_acto1.R`. Agrega TRES figuras propias (todas con `tipo =
"figura_presentacion"` en `procedencia.csv`, no `"figura"`, para no disparar el chequeo
de `12_informe` que exige que toda fila `tipo == "figura"` este embebida en
`docs/informe.html`): un subconjunto de boxplots de placenta con lista de genes propia
(il6, fatp1, slc38a2, glut1 -- distinta de la de la segunda excepcion), un panel
standalone de deteccion de il6 en cerebro (reusa `panel_deteccion()`, D7), y un recorte
de la diagonal del SPLOM de co-expresion (densidades de -ddCt por tratamiento, SOLO
hembras, cerebro, en los genes con interaccion SEXOxTTO significativa sobre la
dispersion, BH<0.05, leidos de `acto2_dispersion_interaccion.csv` -- con resguardo a un
conjunto fijo de 5 genes si esa corrida puntual no deja ninguno bajo el umbral, ya que
los efectos sinteticos son arbitrarios). Genera DOS versiones por el mismo mecanismo de
`fuente_datos()`: sintetica -> `docs/informe_breve.{html,pdf}` (PUBLICA, se versiona);
real -> `outputs/informe_breve_real/` (NUNCA se versiona, cubierto por `.gitignore`).
**Mecanismo de aviso mas granular que el de la version anterior de este pedido:** en vez
de reemplazar la pagina 3/4/5 ENTERA por un aviso generico, se conservan siempre la
narrativa fija de cada pagina, la descripcion de que analisis se hicieron, las figuras
(siempre de la corrida vigente) y las limitaciones; solo la frase puntual que afirma que
DIO un analisis sobre datos reales se reemplaza por un parrafo corto de clase `sim`
(`sim_bloque()`). Las cifras y figuras siempre se leen en vivo de `outputs/tables/` y
`outputs/figures/` tal como esten pobladas al momento de correr el script -- por eso,
para armar la version PUBLICA correctamente, es obligatorio regenerar figuras/tablas
desde datos sinteticos ANTES (`run_all.ps1 -Only R -FromSynthetic`) y restaurar el
estado real DESPUES (`run_all.ps1` completo + `14_informe_breve.R` de nuevo), exactamente
la misma secuencia que exige `13_presentacion.R` (ver leccion documentada en `ESTADO.md`,
sesion de la presentacion: no alcanza con bloquear el texto, las figuras tambien tienen
que venir de la corrida sintetica).

**Correcciones de redaccion (`pedidos/pedido_informe_breve_correcciones.md`, reemplaza
al pedido anterior en lo que difiere):** el texto de las 5 paginas lleva acentuacion
correcta (nunca caracteres acentuados literales en el codigo fuente de R -- entidades
HTML en el cuerpo de cada pagina, que no pasa por `.esc()`; escapes `\uXXXX` en los
titulos de pagina, que si pasan por `.esc()` y esa funcion escaparia una entidad
literal; y `\uXXXX` tambien para el texto que va dentro de una figura PNG, que no es
HTML), p-valores formateados con una unica funcion `.p_fmt()` (nunca notacion
cientifica cruda), sin codigos internos del repositorio en el texto (D7/D13 explicados
en palabras) ni remision a `docs/referencias.md` (reemplazada por citas numeradas
normales, copiadas tal cual de ese archivo, al pie de la primera pagina), y figuras
numeradas ("Figura N.") referenciadas por numero en el texto. Se agrega una CUARTA
figura reusada -- `acto2_corr_placenta_cerebro_fatp4.png`, que ya genera
`08_acto2_correlaciones.R` con su propia fila `tipo = "figura"` en `procedencia.csv`
(este script no le agrega una fila propia) -- con un epigrafe que aclara que la
diferencia entre grupos no fue significativa por test formal, para que el lector no
interprete el patron visual (correlacion aparente en control, ninguna en LPS, en
hembras) como una perdida de correlacion real.

**Convención de nombres de salidas:**

- Figuras: `outputs/figures/actoN_<tema>_<gen|tejido|detalle>.png`, 300 dpi, texto legible.
- Tablas de resultados: `outputs/tables/<R|python>/<contenido>.csv`, **nombres idénticos
  entre R y Python** para que `comparacion_R_python` los cruce automáticamente.
- Toda figura y toda tabla tiene **una fila** en `outputs/tables/procedencia.csv`.
- Todo resultado principal tiene **una fila** en `outputs/tables/verificaciones.csv`.

---

## 4. Decisiones metodológicas no negociables (D1–D13)

> Si hay una razón técnica de peso para que alguna sea inaplicable en un caso puntual,
> **parar y preguntar**. No se cambian por cuenta propia. No se reabren D1–D12.

| # | Decisión | Justificación (1 línea) |
|---|---|---|
| **D1** | Cuantificación relativa: `ΔCt = Ct_gen − Ct_rsp29`. Calibrador = **HEMBRA CONTROL**, por gen × tejido, promediando **solo valores detectados**. `ΔΔCt = ΔCt_muestra − ΔCt_calibrador`. | Método estándar de cuantificación relativa; el calibrador fija el cero biológico. |
| **D2** | Analizar sobre **−ΔΔCt** (log2: simétrica y aditiva). Graficar `FC = 2^(−ΔΔCt)` con eje Y log. | El fold-change está acotado en 0 y es asimétrico → falla los supuestos de modelos lineales casi siempre. |
| **D3** | No detectados → **NA**. Se evaluó imputación MNAR (`nondetects`) y **se descartó**. No imputar por ningún método. | La imputación MNAR introdujo estructura artificial en el baseline Control; ver `analisis_descartados.md`. |
| **D4** | Población de análisis: **los 36 fetos**. Variabilidad entre individuos = biológica (ARN verificado, 260/280 ∈ [1.8, 2.0]). **Sin exclusión de outliers** en el análisis principal. | Excluir individuos "extremos" sesga hacia la hipótesis; los controles de sensibilidad (D8.6) cubren la robustez. |
| **D5** | Modelo: `−ΔΔCt ~ SEXO * TTO`, ajustado por gen × tejido. Cascada de supuestos sobre residuos (ver 4.1). Registrar en procedencia qué rama se usó y por qué. | Un solo modelo con selección de método guiada por diagnóstico, no por conveniencia. |
| **D6** | Post hoc **solo si `SEXO×TTO` es significativa**. 4 comparaciones fijas: ♀Control–♀LPS · ♂Control–♂LPS · ♀LPS–♂LPS · ♀Control–♂Control. Corrección **Holm** dentro de esas 4. | Comparaciones preplanificadas; Holm controla el FWER sin ser tan conservador como Bonferroni. |
| **D7** | `il6` en cerebro E15 **no es cuantificable**: calibrador ♀Control tiene 0/9 detectados. Se **excluye del modelo** y se analiza solo como **proporción de detección** (Fisher exacto 2×2 dentro de cada sexo, y tabla 2×4). Verificar programáticamente "0 detectados en calibrador" y aplicar la misma regla a cualquier otro gen × tejido donde ocurra. | El fold-change se calcularía contra un ΔCt calibrador inexistente. |
| **D8** | Score compuesto de transportadores: **z-score de cada gen por separado** (dentro de cada tejido, sobre los 36 fetos) y **promedio de los 7 z por feto** (usando los z disponibles si el feto tiene <7 detectados). Nombre: *composite z-score / module score / gene set score*; variante PCA = *eigengene*. **Prohibido llamarlo "TONE"** en cualquier archivo, gráfico o texto. | Score robusto y estándar; "TONE" es nomenclatura interna no publicable. |
| **D9** | pSTAT3: `PSTAT3 ~ SEXO * TTO + MEMBRANA` (MEMBRANA = bloque fijo; diseño balanceado). Misma cascada D5 y mismo post hoc D6. **Limitación obligatoria en el informe:** normalizado a proteína total sin STAT3 total → refleja **abundancia de fosfo-STAT3**, no fracción fosforilada. | El bloque absorbe la variación entre membranas; sin STAT3 total no hay fracción. |
| **D10** | ELISA con censura a izquierda: columna indicadora `censurado = TRUE`, valor guardado como `NA`, `LOD` registrado aparte. Reportar **% de censura por grupo antes** de cualquier estadístico. Métodos que manejan censura: KM/ROS (`NADA` en R) o no paramétricos con censurados como empates en el rango más bajo (Peto-Peto / Gehan). Si un grupo tiene censura tan alta que ningún estimador es defendible, decirlo y reportar solo proporción de detección. | Un valor "< blanco" no es un número negativo ni un 0; es información parcial. |
| **D11** | Anotación de boxplots de expresión y de pSTAT3, en **cascada de 3 ramas** (ver 4.2 — **ampliada por pedido explícito del usuario, 2026-09-11**, revierte la restricción previa "solo si la interacción es significativa"): (a) `SEXO×TTO` significativa → brackets por par del post hoc D6, como antes. (b) interacción NO significativa y efecto principal de `TTO` significativo/tendencia → **un solo bracket** que abarca los 4 grupos, etiqueta `Control vs LPS` + estrellas/`p`. (c) interacción NO significativa y efecto principal de `SEXO` significativo/tendencia → bracket entre los centros de cada sexo, etiqueta `♀ vs ♂` + estrellas/`p`. (b) y (c) no son excluyentes entre sí. Símbolos (las 3 ramas): `p<0.001` → `***`; `p<0.01` → `**`; `p<0.05` → `*` (bracket línea llena); `0.05≤p<0.1` → tendencia: bracket punteado + `p` (3 decimales); `p≥0.1` → sin anotar. **Una sola función decide qué anotar** (una en R, una en Python), usada en todas las figuras. | El modelo no sostiene que el efecto difiera entre sexos cuando la interacción no es significativa: marcar los 4 pares sugeriría un dimorfismo no sostenido; el efecto principal es la conclusión que el modelo sí sostiene. |
| **D12** | Corrección entre genes: criterio **primario = sin corrección**. Agregar igualmente, como **columna suplementaria rotulada**, el `p` ajustado por **Benjamini-Hochberg dentro de cada tejido** a través de los genes. **No cambiar ninguna conclusión** por esa columna; mencionar en el informe cuántos resultados sobreviven. | Decisión abierta declarada; la columna BH queda disponible sin dirigir la inferencia. |
| **D13** | Sin término de camada. Los modelos no incluyen `MADRE`, ni como efecto fijo ni aleatorio. En análisis previos de este mismo modelo experimental se evaluó incluirla y no modificaba los resultados. Se asume independencia entre los fetos para el análisis; la limitación se declara en el informe. | En análisis previos de este mismo modelo experimental, agregar `MADRE` no modificaba los resultados; con 2 fetos por camada (1 hembra + 1 macho) no hay grados de libertad para estimarla como aleatoria dentro de cada celda SEXO×TTO. |

### 4.1 Cascada de supuestos (D5)

| Situación | Método | R | Python |
|---|---|---|---|
| Shapiro-Wilk **y** Levene se cumplen | ANOVA tipo III (contrastes suma-cero) | `car::Anova(m, type=3)` con `options(contrasts=c("contr.sum","contr.poly"))` | `statsmodels` sum coding + `anova_lm(typ=3)` |
| Falla homocedasticidad (Levene) | OLS con errores robustos HC3 | `car::Anova(m, white.adjust="hc3")` o `sandwich::vcovHC(type="HC3")` + `lmtest::coeftest` | `fit(cov_type="HC3")` |
| Falla normalidad (Shapiro-Wilk) | ART (Aligned Rank Transform) | `ARTool::art` + `ARTool::art.con` | ART **implementado a mano** (Wobbrock et al.: alinear restando efectos no relevantes estimados, rankear, ANOVA sobre rangos alineados), **verificado contra la salida de R** (verificación obligatoria) |

- Post hoc de ART: **ART-C** (`art.con`), **nunca `emmeans` directo** sobre el modelo ART (infla el error tipo I).
- El orden de la cascada es: normalidad y homocedasticidad OK → ANOVA-III; si falla solo Levene → HC3; si falla Shapiro (con o sin Levene) → ART.

### 4.2 Cascada de anotación de brackets (D11 ampliada)

| Situación | Anotación |
|---|---|
| `SEXO×TTO` significativa (post hoc D6 corrido) | Un bracket por cada comparación de D6 cuyo `p_holm` cruce el umbral (como antes de esta sesión) |
| `SEXO×TTO` NO significativa y `TTO` (efecto principal) significativo/tendencia | Un único bracket que abarca los 4 grupos, etiqueta `Control vs LPS` |
| `SEXO×TTO` NO significativa y `SEXO` (efecto principal) significativo/tendencia | Un único bracket entre los centros de ♀ y ♂, etiqueta `♀ vs ♂` |

- Los dos brackets de efecto principal **no son excluyentes entre sí**: si `TTO` y `SEXO`
  cruzan el umbral a la vez (sin interacción significativa), se apilan ambos.
- Símbolos: mismo criterio en las 3 filas — `p<0.001` → `***`, `p<0.01` → `**`, `p<0.05` →
  `*` (bracket sólido); `0.05≤p<0.1` → bracket punteado + `p = 0.NNN`; `p≥0.1` → sin anotar.
- **Cambio del 2026-09-11** (pedido explícito del usuario, tratando
  `pedidos/boxplots_acto1_base_R.R` como especificación de estilo y de lógica de
  anotación, no como código a copiar): antes, sin interacción significativa, el panel
  quedaba sin ninguna marca — aunque hubiera un efecto principal fuerte (p. ej.
  `il6@PLACENTA_E15` con `p_TTO` = 6.0e-06). Esa restricción se revierte con esta cascada.
- Los boxplots de expresión (`07_figuras_acto1`) pasan a **R base** (no `ggplot2`) porque
  el estilo pedido necesita bigote y tope del bigote con trazo distinto del borde de la
  caja, algo que `geom_boxplot` no expone y `boxplot()`/`bxp()` de R base sí. El Acto 2
  (08) sigue en ggplot2/GGally, y pSTAT3 (mismo script que los boxplots de expresión)
  también sigue en ggplot2. matplotlib (Python) ya permitía estilar los tres trazos por
  separado, así que ese lenguaje no cambió de librería.

---

## 5. Prohibiciones explícitas

1. **No imputar** valores no detectados por ningún método (ni MNAR, ni LOD/2, ni mínimo observado, ni k-NN).
2. **No tratar los ELISA negativos como números negativos**, ni ponerlos en 0, ni borrarlos silenciosamente.
3. **No usar datos crudos fuera de `data/raw/`**, y **no versionarlos** (`data/raw/` en `.gitignore`). Datos inéditos.
4. **No reportar "significativo en Control y no en LPS" como prueba de que las correlaciones difieren.** La ausencia de significancia no es evidencia de diferencia. Para comparar correlaciones: test formal (interacción de pendientes, Fisher z, o permutación), y **ese** test es el que se reporta.
5. **No interpretar cambios de correlación como cambios de coordinación biológica** sin descartar antes cambio de dispersión (restricción de rango). La simulación de la Sección 8.4 existe para eso.
6. **No correr el análisis principal sobre el fold-change** (ver D2).
7. **No usar `emmeans` directo sobre modelos ART** (ver D5).
8. **No llamar "TONE" al score compuesto** (ver D8).
9. **No reabrir las decisiones D1–D12.**
10. **No inventar datos**, ni completar celdas faltantes, ni "arreglar" un archivo de entrada. Si algo no cierra: parar y preguntar.
11. **No hardcodear rutas absolutas** en los scripts.
12. **No dejar figuras sin su fila** en `procedencia.csv`.

---

## 6. Plan de tareas (Sección 10 del brief) y estado

| # | Tarea | Entregable | Cierre sesión | Estado |
|---|---|---|---|---|
| T0 | `AGENTS.md`, `CLAUDE.md`, estructura, `.gitignore`, `README.md` esqueleto | repo inicializado | Sí | **HECHO** (2026-09-07) |
| T1 | `00_config` + `01_generar_sinteticos` (R y Python) | datos sintéticos que corren | Sí | **HECHO** (2026-09-07) |
| T2 | `02_ingesta_qc`: lectura, 40→NA, censura ELISA, tabla de n real | reporte de QC | Sí | **HECHO** (2026-09-07) |
| T3 | `03_elisa` + figura de validación | Acto 1.1 | No | **HECHO** (2026-09-08) |
| T4 | `04_qpcr_cuantificacion`: ΔCt, calibrador, ΔΔCt, z-scores | tablas intermedias | Sí | **HECHO** (2026-09-08) |
| T5 | `05_qpcr_modelos`: cascada D5, post hoc D6, Fisher para D7 | tabla de clasificación | Sí | **HECHO** (2026-09-08) |
| T6 | `06_pstat3` + `07_figuras_acto1` (función de brackets D11) | Acto 1 completo | Sí | **HECHO** (2026-09-09) |
| T7 | `08_acto2_correlaciones`: correlaciones por gen + pair plots | Acto 2.1–2.2 | Sí | **HECHO** (2026-09-09) |
| T8 | `09_acto2_dispersion` + `10_acto2_simulacion` + test de pendientes | Acto 2.3–2.5 | Sí | **HECHO** (2026-09-10) |
| T9 | `11_sensibilidad`: eigengene y exclusión del extremo | Acto 2.6 | Sí | **HECHO** (2026-09-10) |
| T10 | Procedencia, verificaciones, análisis descartados, comparación R/Python | tablas de auditoría | Sí | **HECHO** (2026-09-10) |
| T11 | `12_informe` + `99_verificar` + `run_all.ps1` + log de corrida | informe HTML y PDF | — | **HECHO** (2026-09-10) |

> **T11 cierra el proyecto.** `run_all.ps1` corre en pasadas separadas por lenguaje
> (PRIME python 00–11 → R 00–11/98/12 → Python 00–11/98/12 → VERIFY R/99 y python/99):
> es la única forma de que `12_informe` deje en disco los `.md` y el `informe.html`
> de CADA lenguaje para que `99_verificar` los byte-compare. `12_informe` arma el
> HTML a mano (sin rmarkdown/pandoc: no están en el toolchain y no darían salida
> byte-idéntica); las figuras van incrustadas en base64, así que `informe.html` NO
> es byte-idéntico R/Python (los PNG difieren, ver §7) y la paridad se chequea sobre
> `informe.textonly.html` (sin los `data:`). El PDF es best-effort por Edge/Chrome
> headless y **nunca** frena el pipeline.

**Antes de empezar cada tarea:** decir en dos líneas qué se va a hacer y esperar
confirmación si implica una decisión no fijada acá. **Al terminar cada tarea:**
actualizar `ESTADO.md`, hacer commit descriptivo, e indicar "conviene cerrar la sesión
acá" donde corresponde. **Al abrir sesión nueva:** leer `AGENTS.md` y `ESTADO.md` primero.

---

## 7. Convenciones de implementación (Sección 9 del brief)

- Todo se implementa **dos veces**: R y Python, funcionalmente equivalentes.
- Cada script arranca con secciones tituladas por un comentario breve que explica
  **por qué** se corre ese bloque (no qué hace la función).
- Semilla fija (`SEMILLA = 20260101`), declarada en `00_config`.
- **Solo rutas relativas.** La raíz se resuelve desde la ubicación del script.
- En cada bloque, marcar `# LIBRERÍA: <paquete>::<función>` / `# PROPIO: ...` — alimenta
  `procedencia.csv` (columna `origen_codigo`).
- Salidas numéricas con **nombres idénticos** entre R y Python, en `outputs/tables/R/` y
  `outputs/tables/python/`.
- Comparación R↔Python (`99_verificar` o `09_comparar`): tolerancia declarada
  (`1e-6` estadísticos, `1e-4` p de tests iterativos). Discrepancia por encima de la
  tolerancia → se investiga y se documenta en `analisis_descartados.md`, no se ignora.
- Datos sintéticos (`01_generar_sinteticos`): misma estructura, nombres de columna, tipos
  y **patologías** que los reales — 18 camadas / 36 fetos, `CT_CRUDO` con valores 40, un
  gen con 0 detectados en ♀Control en cerebro (dispara D7), bloque de ELISA con valores
  negativos, pSTAT3 con 3 membranas balanceadas. Los efectos simulados son arbitrarios.

---

## 8. Toolchain (verificado al inicio de T0)

- **Python** 3.13 con numpy, pandas, scipy, statsmodels, matplotlib, seaborn, scikit-learn,
  pingouin, patsy, openpyxl. **A instalar:** `lifelines`, `nbconvert`, `nbformat`, `jinja2`.
- **R** 4.6.1 (`C:\Program Files\R\R-4.6.1` en esta maquina; `run_all.ps1` resuelve
  `Rscript` solo, sin rutas fijas -- ver README 2.0) con car, ARTool, emmeans, survival, sandwich,
  lmtest, nlme, rmarkdown, ggplot2, readxl, writexl. **A instalar:** `NADA`, `GGally`,
  y **pandoc** (para render del informe).
- Si no hubiera LaTeX: `docs/informe.pdf` vía impresión headless del HTML
  (`msedge.exe --headless --print-to-pdf`) — comando documentado en el README.
- El pipeline **no debe fallar** por la generación del PDF; si falla, deja el HTML y avisa.
