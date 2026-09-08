# ESTADO.md — Bitácora de handoff entre sesiones

> Después de cada tarea se escribe acá: qué se completó, qué archivos se generaron, qué
> quedó pendiente, y cuál es el siguiente paso concreto. Al abrir una sesión nueva, lo
> primero que se hace es leer `AGENTS.md` y este archivo.

---

## Sesión 1 — 2026-09-07 — T0

### Qué se completó

- **Confirmación de la Sección 2** con el usuario (respuestas registradas):
  - `ARCHIVO_QPCR = data/raw/Raw data CTs.xlsx`, `HOJA_QPCR = Sheet1` (el "valor probable"
    del brief apuntaba a `RESULTADOS_COMPLETOS_qPCR_AYELEN.xlsx/BD_FINAL`, que ya está
    post-imputación MNAR → incompatible con D1/D3; se descarta como fuente).
  - `RAIZ_REPO = C:\Users\Usuario\Documents\Doctorado\Resultados\LPS 100\qPCR\MIA_LPS_reanalisis`.
  - Instalaciones permitidas (pandoc, lifelines, nbconvert, jinja2, NADA, GGally).
  - `BRAIN_P1` se excluye en la ingesta; se registrará en `analisis_descartados.md`.
- **Inspección de los 3 crudos** (ver AGENTS.md §2.1). Verificado: E15 = 36 fetos;
  `CT_CRUDO` 129 celdas `==40` + 30 NA; `il6`/`BRAIN_E15`/`♀Control` = 0/9 detectados
  (D7 dispara); ELISA MS n=14 (4 `Conc`<0), LA n=28 (10 `Conc`<0).
- **Estructura de carpetas** creada (árbol de AGENTS.md §3), con `.gitkeep` en cada dir.
- **`AGENTS.md`** (+ copia idéntica **`CLAUDE.md`**, md5 verificado): checklist S1,
  estructura + convención de nombres, D1–D12 con justificación, cascada D5, prohibiciones
  S5, plan de tareas S10 con estado, toolchain verificado.
- **`.gitignore`**: `data/raw/*`, `data/processed/*`, `outputs/intermediate/*`, `*.xlsx`
  fuera de `data/synthetic/`, artefactos R/Python/OS.
- **`README.md`** esqueleto: instalación R + Python, obtención de datos, `.\run_all.ps1`,
  salidas esperadas, advertencia de validez solo con datos reales.
- **`requirements.txt`** con versiones ancladas de lo instalado + `lifelines`, `nbconvert`,
  `nbformat`, `jinja2`.
- **`run_all.ps1`** esqueleto (orden de scripts, aún sin cuerpo real).
- `git init`, copia de los 3 crudos a `data/raw/`, verificación `git status`.

### Archivos generados

```
MIA_LPS_reanalisis/
├─ AGENTS.md  CLAUDE.md  README.md  ESTADO.md  .gitignore
├─ requirements.txt  run_all.ps1
├─ data/{raw,synthetic,processed}/.gitkeep
├─ data/raw/  <- 3 Excel crudos copiados (IGNORADOS por git)
├─ R/.gitkeep  python/.gitkeep
├─ outputs/{figures,tables,intermediate}/.gitkeep
└─ docs/.gitkeep  logs/.gitkeep
```

### Pendiente / siguiente paso concreto

- **T1**: `R/00_config.R` + `python/00_config.py` (rutas relativas, semilla, constantes:
  lista de genes, 7 transportadores, 3 vía IL-6, housekeeping, mapa de grupos, tolerancias
  de comparación) y `R/01_generar_sinteticos.R` + `python/01_generar_sinteticos.py`
  (generador que reproduce estructura + patologías de los 3 crudos → `data/synthetic/`).
  Cerrar sesión al terminar T1.
- Al cerrar T1: `renv::init()` / `renv::snapshot()` para generar `renv.lock`.

### Notas para la próxima sesión

- Los supuestos menores (4 negativos de MS como censura; 30 NA de `CT_CRUDO` = no
  detectado; `LOD = Conc 0`) están anotados en AGENTS.md §2.1 y pueden revisarse en T2/T3.
- pandoc / paquetes R faltantes todavía **no** instalados: hacerlo al inicio de T11
  (informe) o antes si T3 necesita `NADA`/`lifelines` para el ELISA con censura.

---

## Sesión 2 — 2026-09-07 — T1

### Qué se completó

- **`python/00_config.py` + `R/00_config.R`** (equivalentes):
  - Raíz del repo resuelta desde la ubicación del script (`__file__` / `--file=` /
    `sys.frame()$ofile`), **sin rutas absolutas**. Todas las rutas de datos y
    salidas derivan de ahí y las carpetas regenerables se crean si faltan.
  - `ruta_datos(archivo)`: devuelve `data/raw/<archivo>` si existe, si no
    `data/synthetic/<archivo>`. **Nunca escribe en `data/raw/`.** Un clon limpio
    (sin crudos) cae solo al sintético. `fuente_datos()` informa cuál se usaría.
  - Constantes: `SEMILLA=20260101`; `GENES` (10, orden canónico del crudo);
    `GENES_TRANSPORTADORES` (7), `GENES_VIA_IL6` (3 = gp130/il6/il6R),
    `GEN_HOUSEKEEPING="rsp29"`; `TEJIDOS_E15`, `TEJIDO_EXCLUIDO="BRAIN_P1"`;
    `NIVELES_SEXO/TTO`, `CALIBRADOR` (HEMBRA/CONTROL, D1), `etiqueta_grupo()` que
    reproduce el espaciado inconsistente real (`♀Control` pegado, resto con
    espacio), `normalizar_sexo/tto()`; `TOL_ESTADISTICO=1e-6`,
    `TOL_P_ITERATIVO=1e-4`.
  - **RNG propio determinista** (`nuevo_rng`): LCG 32-bit (Numerical Recipes) +
    normal por método polar de Marsaglia con caché (solo `sqrt`/`log` → bits
    idénticos en R y Python; se evita `sin/cos` a propósito). Verificado:
    `unif1() x3 = [0.0873843804, 0.7218535838, 0.5726731312]` idéntico en ambos.
- **`python/01_generar_sinteticos.py` + `R/01_generar_sinteticos.R`** (equivalentes):
  generan `data/synthetic/` con la estructura y las patologías de los 3 crudos.
  Consumo del RNG en el mismo orden (qpcr → elisa → curva → pstat3).
  - `Raw data CTs.xlsx` / `Sheet1`: 1000 filas, 10 columnas en el orden real
    (`MADRE…TEJIDO, rsp29, GEN, CT_CRUDO`). 18 madres (C1–C9, L1–L9), 36 fetos
    E15 (1 hembra + 1 macho c/u), `PLACENTA_E15`+`BRAIN_E15` para los 36 y
    `BRAIN_P1` para C1–C7/L1–L7 (28 fetos, 6 con `FETO` vacío → 22 únicos + 6
    vacíos = 280 filas). No-detectados: `CT_CRUDO == 40` (86) **y** celda vacía
    (30, exactamente 6 en c/u de fatcd36/glut1/glut3/slc38a1/slc38a2), ambas =
    no detectado (D3). `rsp29` sin 40 ni NA, constante por feto×tejido. **D7:**
    `il6`/`BRAIN_E15`/♀Control = 0/9 detectados forzado (dispara la regla); ~11
    detectados en el resto de `il6`/cerebro (el crudo tiene 7).
  - `ELISA IL6 2026 Dosis 100.xlsx`: hojas `CURVA IL6` + `Sueros y LA` (orden del
    crudo). `Sueros y LA` = grilla 44×14 sin encabezado real (rótulo `IL6` en
    fila 0, encabezados con espacios finales `'Conc '`/`'IL-6 '` en fila 1,
    scratch `Abs/m/oo` + `0.0016`/`0.1057` a la derecha). Bloque MS n=14 (4
    `Conc`<0) + bloque LA n=28 (8 `Conc`<0), partidos por la columna 5. Relación
    real `Conc = (Abs−0.1057)/0.0016`, columna `IL-6` = `max(Conc,0)` (a NO usar).
  - `pstat3 placenta.xlsx` / `Sheet1`: 36 filas, 9 columnas, `MEMBRANA` int
    1/2/3 balanceada (3×12; 3 réplicas por membrana×sexo×tto), `SEXO`/`TTO` en
    Title case (`Hembra`/`Control`, distinto del qPCR en mayúsculas).
  - Mirrors canónicos `.tsv` + `MANIFEST.tsv` (sha256 por hoja) + `PATOLOGIAS.tsv`
    (contadores). Cada script hace `assert`/`stopifnot` de los invariantes duros.

### Paridad R ↔ Python

Los **6 `.tsv` canónicos** (`Raw data CTs`, `pstat3 placenta`, `ELISA … Sueros y
LA`, `ELISA … CURVA IL6`, `MANIFEST`, `PATOLOGIAS`) salen **byte a byte idénticos**
entre `Rscript R/01…` y `python python/01…` (verificado con `cmp`). El único ajuste
necesario fue no redondear la columna scratch `media` de la curva (caso "mitad
exacta" que `round()` resuelve distinto en cada lenguaje).

### Archivos generados / modificados

```
python/00_config.py          R/00_config.R
python/01_generar_sinteticos.py   R/01_generar_sinteticos.R
data/synthetic/Raw data CTs.xlsx  + .tsv
data/synthetic/pstat3 placenta.xlsx  + .tsv
data/synthetic/ELISA IL6 2026 Dosis 100.xlsx  + __Sueros y LA.tsv + __CURVA IL6.tsv
data/synthetic/MANIFEST.tsv   data/synthetic/PATOLOGIAS.tsv
ESTADO.md (esta entrada)
```

### Desviaciones menores respecto del crudo (documentar en `analisis_descartados.md` en T2/T10)

- **`ELISA … CURVA IL6`**: se reproduce solo el bloque inferior de la curva
  (`Abs 450` / `pg/ml`, 8 estándares 0…1000 en filas 15–22); el crudo además
  trae arriba una matriz de réplicas de absorbancia (filas 1–8) que no se replica
  (no la usa el pipeline E15/placenta; solo se referencia como alternativa de LOD
  en D10). El sintético queda 22×4 en vez de 22×10.
- **`il6`/`BRAIN_E15` detectados fuera del calibrador**: sintético ~11 vs 7 en el
  crudo; el invariante que importa (0/9 en el calibrador, D7) se cumple exacto.
- **`CT_CRUDO == 40`**: 86 en el sintético vs 129 en el crudo (mismo fenómeno y
  misma forma —literal 40 + celda vacía—, proporción del mismo orden).
- **Conc<0 en LA**: 8 vs 10 en el crudo.

### Dependencias nuevas usadas (ya instaladas; las capturará `renv::snapshot()`)

- R: `openxlsx` (4.2.8.1, escritura de `.xlsx` con tipo por celda) y `digest`
  (sha256 de los `.tsv`). Python: `openpyxl` (ya en `requirements.txt`).

### Cierre de T1

- `renv` 1.2.4 instalado y **`renv::init()` ejecutado** desde `RAIZ_REPO`:
  crea `renv/` (`activate.R`, `settings.json`, `.gitignore`), `.Rprofile`
  (`source("renv/activate.R")`) y **`renv.lock`**. La lockfile captura por ahora
  solo lo que usan los scripts (`openxlsx` 4.2.8.1, `digest` 0.6.39 + deps
  `cli`, `Rcpp`, `stringi`, `zip`, `renv`) y R 4.6.1; se ampliará con cada tarea
  que agregue `library(...)` y un nuevo `renv::snapshot()`.
- `R/01_generar_sinteticos.R` re-corrido bajo renv → hashes idénticos. OK.
- Tabla de estado en AGENTS.md/CLAUDE.md: **T1 → HECHO (2026-09-07)**.
- **T1 cerrado.** Conviene cerrar la sesión acá.

### Siguiente paso concreto

- **T2** (`02_ingesta_qc`): lectura de los 3 archivos vía `ruta_datos()`, `40→NA`
  y celda vacía→NA unificados, exclusión de `BRAIN_P1` (registrar en
  `analisis_descartados.md`), censura ELISA (`Conc<0` → `censurado=TRUE`, valor
  `NA`, `LOD` aparte) con **detección de la fila de encabezado** y coerción
  numérica (para leer igual el crudo y el sintético), y tabla de n real por
  archivo × grupo × sexo.

### Notas para la próxima sesión

- `run_all.ps1`: previsto que regenere el sintético **con un solo lenguaje**
  (p. ej. `python python/01…`) para no alternar el contenido de `data/synthetic/`;
  la equivalencia con R se chequea aparte en T10/`99_verificar` comparando los
  `.tsv`/`MANIFEST`.
- Emojis/VS16: `NOMINACION` = `<madre>` + `♀`/`♂` + `U+FE0F`; `GRUPO` sin VS16;
  LA col0 = `"Control ♂"` (espacio + emoji, sin VS16). Ya reproducido.

---

## Sesión 3 — 2026-09-07 — T2

### Qué se completó

- **`python/02_ingesta_qc.py` + `R/02_ingesta_qc.R`** (equivalentes): única puerta
  de entrada de datos al pipeline. Leen los 3 Excel vía `ruta_datos()` (crudo real
  si está en `data/raw/`, si no el sintético), formato largo + saneo, sin tocar
  otras decisiones.
  - **qPCR**: detección de encabezado por nombres exactos; `CT_CRUDO == 40` y celda
    vacía unificadas → `no_detectado=TRUE`, `CT=NA` (D3, sin imputar); `CT_CRUDO`
    original se conserva para auditar. `BRAIN_P1` se aparta (`qpcr_brain_p1_excluido.tsv`)
    y se registra en `analisis_descartados.md`. Invariantes duros: 36 fetos E15,
    720 filas, tejidos = {PLACENTA_E15, BRAIN_E15}, `rsp29` sin NA/40, `CT_CRUDO`
    sin no-numéricos fuera de vacío. `MADRE_ID` = `FETO` sin sufijo `.NN`.
  - **ELISA**: encabezado localizado por la celda `TEJIDO`; hoja partida en MS/LA
    por esa columna (ignora scratch cols 6–13); `TTO` de la col 0, `SEXO` sólo en
    LA (del emoji), MS sin sexo fetal. `Conc < 0` → censura a izquierda D10:
    `censurado=TRUE`, `IL6_pgml=NA`, `LOD=0` aparte. Col `IL-6` del crudo NO se usa
    (registrado en `analisis_descartados.md`: es `max(Conc,0)` salvo MS/LPS que es
    `Conc×4`).
  - **pSTAT3**: `SEXO`/`TTO` normalizados (venían en Title case); verifica 3
    membranas balanceadas (12/12/12).
- **Salidas** (todas byte-idénticas R↔Python, verificado con `cmp`):
  - `data/processed/`: `qpcr_e15_long.tsv`, `qpcr_brain_p1_excluido.tsv`,
    `elisa_long.tsv`, `pstat3_long.tsv` (intermedios para T3–T6).
  - `outputs/tables/{R,python}/`: `qc_n_por_grupo.csv` (n real por archivo×grupo×sexo),
    `qc_no_detectados_qpcr.csv` (80 filas gen×tejido×grupo; marca D7),
    `qc_censura_elisa.csv` (% censura por bloque×grupo, D10),
    `qc_faltantes.csv` (IDs de referencia sin dato y viceversa),
    `qc_resumen.csv` (chequeos máquina-legibles).
  - `outputs/tables/`: `qc_reporte.md` (reporte legible), y alta de
    `analisis_descartados.md` (sección `02_ingesta_qc`, delimitada por marcadores
    HTML para reemplazo idempotente), `procedencia.csv` y `verificaciones.csv`
    (merge por columna `script`, id **sin extensión** `02_ingesta_qc` para que R y
    Python reemplacen la misma fila → archivo único).

### Números (datos reales, para referencia; NO se versionan)

qPCR E15: 36 fetos / 720 filas; BRAIN_P1 excluido 280 filas; CT no detectado
129 (`==40`) + 30 (vacío) = 159. D7 dispara en `il6@BRAIN_E15` (0/9 en calibrador).
ELISA MS n=14 (4 censurados), LA n=28 (10 censurados). pSTAT3 36 filas, membranas
1:12|2:12|3:12. Faltan en ELISA-MS 4 madres (`C_240719_1`, `C_240724_1`,
`C_241029_1`, `C_241031_1`); ELISA-LA 8 sacos. pSTAT3 cruza 36/36 con qPCR E15.
(Sobre sintético: CT=40 → 86, LA censurados → 8; el resto igual.)

### Decisiones / desviaciones de esta sesión

- **`leer_grid()` bifurca por fuente:** datos **reales** → `.xlsx` de `data/raw/`
  (Python: `openpyxl`; R: `readxl` con `col_types="text"` → los números vuelven como
  texto de precisión completa y `as.numeric()`/`float()` recuperan el mismo `double`
  → `%.10g` coincide). Datos **sintéticos** → el **mirror `.tsv` canónico** de
  `data/synthetic/` (`<archivo sin .xlsx>.tsv` o `…__<hoja>.tsv`), que es la forma
  versionada, determinista y byte-idéntica R/Python creada en T1. Fallback al `.xlsx`
  sintético sólo si no hubiera `.tsv`.
  - Motivos para NO leer el `.xlsx` sintético: (a) `openxlsx::read.xlsx` (R)
    mal-parsea cadenas con `xml:space="preserve"` → los encabezados `"Conc "` /
    `"IL-6 "` del ELISA salían como `xml:space="preserve">Conc `; (b) `openpyxl`
    (Python) **no podía re-leer** el `.xlsx` que escribe `openxlsx` (`KeyError:
    xl/drawings/drawing1.xml`); (c) el `.xlsx` de `openpyxl` **no es
    byte-determinista** (mete timestamps en el zip). El `.tsv` esquiva las tres.
  - Costo asumido: el `rsp29`/CT del sintético leído del `.tsv` tiene 10 cifras
    (`%.10g`) en vez de la precisión completa del `.xlsx`. Es consistente entre R y
    Python (ambos usan el `.tsv`) y está muy dentro de `TOL_ESTADISTICO = 1e-6`.
- **`data/synthetic/*.xlsx` dejan de versionarse** (`.gitignore`: `*.xlsx` sin
  excepción; `git rm --cached` de los tres). `01_generar_sinteticos` los sigue
  generando para inspección manual; la forma canónica del fixture son los `.tsv` +
  `MANIFEST.tsv` (hash del `.tsv`, sin cambios).
- **Lector R:** se instaló `readxl` 1.5.0 + 17 deps en renv y se corrió
  `renv::snapshot()` (renv.lock actualizado). Sólo se usa en la rama de datos reales.
- **`.gitignore`**: además de los `.xlsx`, se agregan `outputs/tables/*` y
  `outputs/figures/*` (salvo `.gitkeep` y los subdirs `R/`, `python/`). Son
  regenerables por `run_all` y, sobre datos reales, contendrían números de datos
  inéditos (AGENTS §2/§5). `99_verificar` igual chequea su existencia/no-vacuidad.

### Verificado

- R↔Python byte-idénticos en las 4 tablas `data/processed/*.tsv`, las 5
  `qc_*.csv` (× R/ y python/), `qc_reporte.md`, `procedencia.csv`,
  `verificaciones.csv`, `analisis_descartados.md` — sobre datos reales y sobre
  sintético.
- Idempotente: 2ª corrida de cada implementación no cambia ningún archivo.
- Ambas rutas de `leer_grid` verificadas: real (`.xlsx`) y sintético (`.tsv`, con
  los `.xlsx` sintéticos borrados) dan salidas byte-idénticas R↔Python.
- `renv::status()` limpio tras el snapshot (readxl + deps en el lock).

### Siguiente paso concreto

- **T3** (`03_elisa`): análisis del ELISA IL-6 con censura a izquierda (D10) sobre
  `data/processed/elisa_long.tsv` — % de detección/censura por grupo primero,
  luego KM/ROS (`NADA` en R) o no paramétrico con censurados como empates en el
  rango más bajo; figura de validación (Acto 1.1). MS por tratamiento, LA por
  tratamiento×sexo. Instalar `NADA` (R) y `lifelines` (Python) si hacen falta.
  T3 NO cierra sesión (según plan).

### Notas para la próxima sesión

- Los `.tsv` de `data/processed/` los escriben ambas implementaciones (cada corrida
  pisa con contenido idéntico). En `run_all` (R y luego Python) el estado final lo
  deja Python; la paridad se chequea en T10.
- `qc_no_detectados_qpcr.csv` ya trae la columna que necesita T5 para D7: filtrar
  `GRUPO == "HEMBRA_CONTROL" & n_total > 0 & n_detectado == 0`.

---

## Sesión 4 — 2026-09-08 — T3

### Decisiones confirmadas con el usuario (dentro del menú de D10)

- **Contraste primario con censura = Peto-Peto** (Fleming-Harrington G-rho = 1),
  con los censurados empatados en el rango más bajo. NO ROS/KM-NADA.
- **Suero materno (MS): solo detección + Fisher exacto** Control vs LPS. Control
  tiene 4/5 censurados (80%) → ningún estimador de ubicación es defendible
  (cláusula final de D10). Regla fija para MS (no depende de umbral ni de la
  fuente de datos).
- **Líquido amniótico (LA): estratificado por sexo** — Peto-Peto Control vs LPS
  dentro de ♀ y dentro de ♂, + Fisher de detección por sexo. Sin modelo
  factorial SEXO×TTO (censura 0–60% + n 5–9 por celda).

### Qué se completó

- **`python/03_elisa.py` + `R/03_elisa.R`** (equivalentes). Leen
  `data/processed/elisa_long.tsv`. No tocan ninguna decisión fuera de D10.
  - **Detección/censura por grupo** (D10 paso 1, antes de cualquier estadístico):
    `elisa_deteccion.csv` (6 grupos: MS×2 + LA×4).
  - **Descriptivo** `elisa_descriptivo.csv`: n, n_det, %censura, min/mediana/max
    de detectados, y **mediana Kaplan-Meier** del dato reflejado
    `t' = M - t` (`M = ceil(max Conc detectada) + 1`; sobre real M=1407). Vacío
    cuando `S(t)` no baja de 0.5 (censura alta) → MS Control y LA Control ♀.
  - **Peto-Peto** `peto_peto_2grupos()` — **PROPIO**, misma fórmula G-rho que
    `survival::survdiff`. Reflexión censura izquierda → derecha
    (`M - Conc` evento, `M - LOD` censura a derecha) para que survdiff, que
    maneja censura a derecha, deje a los "< LOD" en el rango más bajo.
    **R verifica en corrida** `abs(hand − survdiff(rho=1)) < 1e-8` (`stopifnot`).
    `elisa_petopeto_la.csv` (♀, ♂).
  - **Fisher exacto 2×2** `fisher_2x2()` — **PROPIO**, misma regla de dos colas
    que `stats::fisher.test` (suma de tablas con `prob ≤ prob(obs)·(1+1e-7)`);
    OR con corrección de Haldane (`+0.5`). **R verifica en corrida** contra
    `fisher.test` (`< 1e-9`). `elisa_fisher_deteccion.csv` (MS, LA♀, LA♂).
  - **Figuras Acto 1.1** (300 dpi): `outputs/figures/acto1_elisa_ms.png`
    (Conc por tratamiento; censurados = símbolo abierto en el LOD; box de
    detectados si n≥3; Fisher p en subtítulo) y `acto1_elisa_la.png` (facet por
    sexo; Peto-Peto p por panel). R con ggplot2, Python con matplotlib —
    equivalentes, no byte-idénticas (son PNG).
  - `elisa_reporte.md` legible; sección `03_elisa` en `analisis_descartados.md`
    (marcadores `<!-- 03_elisa:inicio/fin -->`); filas en `procedencia.csv` y
    `verificaciones.csv` (merge por `script`, id sin extensión `03_elisa`).

### Resultados (datos reales; NO se versionan)

- MS: Control 1/5 detectados vs LPS 9/9. Fisher p = 4.995e-03. (Descriptivo:
  mediana detectados MS-LPS = 567.5 pg/mL; MS-Control no estimable.)
- LA ♀: Peto-Peto χ² = 1.0664, p = 0.302. LA ♂: χ² = 0.9375, p = 0.333.
  Fisher detección LA♀ p = 1.0, LA♂ p = 0.258.
- (Sobre sintético M = 1132; MS igual patrón 4/5 censura en Control; LA con otra
  censura → LA♀ Peto-Peto p = 0.031, LA♂ p = 0.795. Los números cambian con la
  fuente pero R==Python para cada fuente.)

### Paridad y robustez verificadas

- **R ↔ Python byte-idénticos** en `elisa_deteccion.csv`,
  `elisa_descriptivo.csv`, `elisa_fisher_deteccion.csv`,
  `elisa_petopeto_la.csv` (× `R/` y `python/`), `elisa_reporte.md`,
  `procedencia.csv`, `verificaciones.csv`, `analisis_descartados.md` — **sobre
  datos reales y sobre sintético**.
- Idempotente: 2ª corrida de cada implementación no cambia ningún archivo.
- Truco de paridad: los p-valores y el `chisq` que dependen de funciones
  trascendentes (`lgamma`, `exp`, `erfc`/`pchisq`) se guardan como texto
  `"%.6e"` ya formateado → mismo string aunque la libm difiera en el último bit
  (documentado en `analisis_descartados.md §03_elisa`). Los conteos, %, medianas
  y OR (aritmética exacta) van con `%.10g`.

### Entorno

- **ggplot2 4.0.3** instalado en renv (+ deps: scales, farver, gtable, isoband,
  labeling, RColorBrewer, viridisLite, withr, S7) y **`survival` 3.8-6** +
  `Matrix`/`lattice` capturados. `renv::snapshot()` corrido → `renv.lock` (38
  paquetes), `renv::status()` limpio.
- **`renv.lock` Repositories y `renv/settings.json`**: CRAN pasó de
  `packagemanager.posit.co` (fallaba acá con error SSL 60) a
  `https://cloud.r-project.org`; `ppm.enabled: null → false`. Mejora de
  reproducibilidad en redes donde Posit PM está bloqueado.
- **`NADA` (R) y `lifelines` (Python) NO se instalaron**: la elección de
  Peto-Peto + KM-reflejado no los necesita.

### Siguiente paso concreto

- **T4** (`04_qpcr_cuantificacion`): sobre `data/processed/qpcr_e15_long.tsv`,
  calcular `ΔCt = CT_gen − CT_rsp29` por muestra; calibrador **HEMBRA_CONTROL**
  por gen×tejido promediando **solo detectados** (D1); `ΔΔCt` y analizar sobre
  **−ΔΔCt** (D2); z-score por gen dentro de tejido sobre los 36 fetos y score
  compuesto = promedio de los 7 z de transportadores por feto (D8, **prohibido
  "TONE"**). `il6@BRAIN_E15` se deja fuera de la cuantificación (D7, 0/9 en
  calibrador) — solo entra como proporción de detección en T5. Tablas
  intermedias en `data/processed/` + `outputs/tables/{R,python}/`. T4 **sí**
  cierra sesión.

### Notas para la próxima sesión

- La regla D7 se detecta programáticamente con
  `qc_no_detectados_qpcr.csv` (filtrar `GRUPO=="HEMBRA_CONTROL" & n_total>0 &
  n_detectado==0`) — hoy solo `il6@BRAIN_E15`.
- Las figuras del ELISA las escriben ambas implementaciones a la misma ruta; en
  `run_all` (R y luego Python) queda la de Python. No son byte-idénticas (PNG);
  la equivalencia es visual/estructural.
- Helpers de escritura/merge (`.fmt`, `escribir_csv`, `merge_por_script`,
  `registrar_procedencia/verificaciones`) están duplicados en cada script
  numerado (patrón del proyecto desde 01/02); si crecen mucho, evaluar
  factorizarlos en T10.

---

## Sesión 5 — 2026-09-08 — T4

### Encuadre confirmado con el usuario

- **T4 = solo cuantificación.** NO se corre `nondetects`, NO se implementa EM, NO
  se imputa por ningún camino (D3 está tomada, no se re-litiga).
- z-score con **desvío estándar muestral (n−1)**.
- **Eigengene / PCA diferido a T9** (Acto 2.6). T4 entrega solo el score
  compuesto (promedio de los 7 z de transportadores).
- "Documentar la MNAR con números" = una **tabla descriptiva de no detectados**
  (n y % de NA por gen×tejido×grupo + marca de calibrador ♀Control 0/detectados).
  Esa tabla es el sustento numérico de D3 **y** de D7. La redacción del descarte
  MNAR se escribe en T10 apoyándose en ella (cualitativo en la lógica,
  cuantitativo en los conteos).

### Qué se completó

- **`python/04_qpcr_cuantificacion.py` + `R/04_qpcr_cuantificacion.R`**
  (equivalentes). Leen `data/processed/qpcr_e15_long.tsv`.
  - **D1** `dCt = CT_gen − CT_rsp29`; calibrador = promedio de `dCt` en
    `HEMBRA_CONTROL` por gen×tejido, **solo detectados**; `ddCt = dCt − dCt_cal`.
  - **D2** medida de análisis = `neg_ddCt = −ddCt`. `FC = 2^(−ddCt)` **no se
    guarda** (se calcula al graficar en T6/T7) — así se evita el redondeo de
    `2^x` entre libm.
  - **D7** `il6@BRAIN_E15` no cuantificable (calibrador ♀Control 0/9):
    `cuantificable = FALSE`, `ddCt/neg_ddCt/z` NA. Único caso; se detecta
    programáticamente (`HEMBRA_CONTROL n_detectado==0 & n_total>0`).
  - **D8** z-score por gen dentro de tejido sobre los 36 fetos, solo detectados,
    **sd muestral (n−1)**. **Score compuesto** = promedio de los z de los 7
    transportadores por feto×tejido (los disponibles si <7). **No se usa "TONE"**
    en ningún archivo.
  - Tabla `qpcr_no_detectados_descriptivo.csv`: n y % de NA por
    gen×tejido×{4 grupos + TODOS} + `calibrador_cero_detectados` (sustento D3/D7).

### Salidas

- `data/processed/qpcr_cuantificacion_long.tsv` (720 filas): …, `dCt`,
  `dCt_calibrador`, `ddCt`, `neg_ddCt`, `z`, flags `no_detectado`,
  `cuantificable`, `es_transportador`.
- `data/processed/qpcr_score_compuesto_long.tsv` (72 filas): `n_z_disponibles`,
  `score_compuesto` por feto×tejido.
- `outputs/tables/{R,python}/`: `qpcr_calibradores.csv` (20),
  `qpcr_no_detectados_descriptivo.csv` (100), `qpcr_cuantificacion_resumen.csv`
  (20: n_det/36, media/sd/mediana de `neg_ddCt`),
  `qpcr_score_compuesto_resumen.csv` (8: media/sd del score por grupo×tejido).
- `outputs/tables/qpcr_cuantificacion_reporte.md`; sección
  `04_qpcr_cuantificacion` en `analisis_descartados.md`; filas en
  `procedencia.csv` y `verificaciones.csv`.

### Números (datos reales; NO se versionan)

- Calibradores: 20 gen×tejido, todos con ≥5 detectados en ♀Control salvo
  `il6R@BRAIN_E15` (3/9, cuantificable) y `il6@BRAIN_E15` (0/9, **no**
  cuantificable, D7).
- `dCt` NA = 77 = exactamente los no detectados E15 → nada imputado.
  (El "159" de `qc_resumen` de T2 es el total E15+BRAIN_P1; E15 solo = 77, todos
  `CT_CRUDO==40`; las 30 celdas vacías están todas en BRAIN_P1.)
- Score compuesto (media por grupo): PLACENTA ♀Control −0.93, ♀LPS +0.18,
  ♂Control +0.11, ♂LPS +0.41 · BRAIN ♀Control −0.50, ♀LPS +0.68,
  ♂Control −0.09, ♂LPS −0.34. (Descriptivo; el modelo es T5.)

### Paridad y robustez

- **R ↔ Python byte-idénticos** en las 4 `qpcr_*.csv` (× `R/` y `python/`), los
  2 `.tsv` de `data/processed/`, `qpcr_cuantificacion_reporte.md`,
  `procedencia.csv`, `verificaciones.csv`, `analisis_descartados.md` — **sobre
  datos reales y sobre sintético**. Idempotente.
- Truco de paridad clave: **promedios y desvíos con acumulador `double`
  explícito** (`a <- a + v` en un loop), NO `sum()`/`mean()`/`sd()` de R (usan
  `long double` en Windows/MinGW → divergían del `double` de Python). Con el loop
  explícito y el mismo orden de suma (fetos ordenados por `MADRE_ID, FETO`
  radix), `dCt_calibrador`, `z` y `score_compuesto` salen bit-idénticos → todo
  con `%.10g`, sin necesidad de formateo `%.6e` como en T3.
- R/04 corre ~40–50 s (lookups por clave-string en loops); aceptable.

### Cierre de T4

- Tabla de estado en AGENTS.md/CLAUDE.md: **T4 → HECHO (2026-09-08)**.
- **T4 cerrado. Conviene cerrar la sesión acá.**

### Siguiente paso concreto

- **T5** (`05_qpcr_modelos`): sobre `qpcr_cuantificacion_long.tsv`, modelo
  `neg_ddCt ~ SEXO * TTO` por gen×tejido con la **cascada de supuestos D5/§4.1**
  (Shapiro + Levene → ANOVA III `car::Anova(type=3)` con `contr.sum`; falla
  Levene → HC3; falla Shapiro → **ART**, en Python implementado a mano y
  verificado contra `ARTool`). Post hoc **solo si `SEXO×TTO` significativa**: 4
  comparaciones fijas, corrección **Holm** (D6). Para `il6@BRAIN_E15` (D7):
  fuera del modelo, solo **Fisher exacto 2×2 de detección dentro de cada sexo** +
  tabla 2×4. Columna suplementaria BH entre genes por tejido (D12), sin cambiar
  conclusiones. Entregable: tabla de clasificación por gen×tejido (qué rama de la
  cascada, p de interacción, post hoc). T5 **sí** cierra sesión.

### Notas para la próxima sesión

- Instalar en renv, al abrir T5: `car`, `ARTool`, `emmeans`, `sandwich`,
  `lmtest` (están en el toolchain de AGENTS §8 pero **no** en la lockfile
  mínima; hoy la lock tiene 38 paquetes: readxl/openxlsx/digest/ggplot2/survival
  + deps). Python: `statsmodels`, `scipy`, `pingouin`, `patsy` ya están.
- El z-score y el score compuesto de T4 **no** entran al modelo de T5 (ese modela
  `neg_ddCt` por gen); el score compuesto se usa en el Acto 2 (T7–T9).
- `qpcr_cuantificacion_long.tsv` ya trae `cuantificable` y `es_transportador`
  para filtrar directo en T5.
