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

---

## Sesión 6 — 2026-09-08 — T5

### Decisiones confirmadas con el usuario (dentro del menú de D5/D6/D12)

- **alfa = 0.05** para la compuerta D6 (interacción `SEXO×TTO` que habilita el
  post hoc) y para la cascada de supuestos. **Shapiro-Wilk sobre los residuos del
  modelo conjunto** del gen×tejido (no por celda); Levene por celda a 0.05.
- **Piso de celda = 5 detectados.** Si alguna de las 4 celdas `SEXO×TTO` tiene
  <5 valores detectados tras descartar no-detectados, no se ajusta el modelo
  factorial: vía `descriptivo_n_bajo` (solo n + descriptivo de `neg_ddCt`).
- **D12: tres columnas BH** — BH entre los genes modelados de cada tejido, por
  separado para `p_SEXO`, `p_TTO`, `p_SEXO×TTO`.

### Qué se completó

- **`python/05_qpcr_modelos.py` + `R/05_qpcr_modelos.R`** (equivalentes). Leen
  `data/processed/qpcr_cuantificacion_long.tsv`. Una fila de clasificación por
  gen×tejido (20). Vías: `modelo` / `D7_deteccion` / `descriptivo_n_bajo`.
  - **Cascada D5 / §4.1**: modelo `neg_ddCt ~ SEXO * TTO` (OLS, contr.sum).
    Shapiro≥.05 & Levene≥.05 → `anova3` (SS tipo III por comparación de modelos);
    Shapiro≥.05 & Levene<.05 → `hc3` (Wald III 1 gl con sándwich HC3);
    Shapiro<.05 → `art` (ART Wobbrock et al. 2011). La rama usada queda en
    `rama_cascada` y en `procedencia.csv`.
  - **Post hoc D6** (solo si interacción p<.05): 4 comparaciones fijas
    (`♀C-♀L`, `♂C-♂L`, `♀L-♂L`, `♀C-♂C`), corrección **Holm** dentro de las 4.
    `anova3` → contraste de medias marginales (vcov OLS); `hc3` → ídem con vcov
    HC3; `art` → **ART-C** (Elkin et al. 2021), nunca `emmeans` directo sobre el
    modelo ART.
  - **D7** (`cuantificable == FALSE`, calibrador ♀Control 0/detectados):
    `il6@BRAIN_E15` fuera del modelo → Fisher exacto 2×2 de detección Control vs
    LPS dentro de cada sexo + tabla 2×4. Se detecta programáticamente.
  - **D12**: columnas `p_*_BH` suplementarias, no dirigen la inferencia.
  - **Núcleo numérico PROPIO idéntico R/Python**: solver 4×4 (Gauss-Jordan),
    SS III por comparación de modelos, sándwich HC3, Levene (Brown-Forsythe,
    centro = mediana), ART y ART-C, Holm, BH, Fisher 2×2 (misma regla que
    `fisher.test`). **R cruza-verifica en corrida** cada rama contra
    `car::Anova` / `emmeans` (vcov OLS o HC3) / `ARTool::art` + `art.con`
    (`stopifnot`, tol 1e-6). **Única dependencia de librería en el resultado:**
    Shapiro-Wilk (`shapiro.test` / `scipy.stats.shapiro`) — coinciden a ~1e-12.

### Resultados (datos reales; NO se versionan)

- Vías: 18 modelados, 1 D7 (`il6@BRAIN_E15`), 1 descriptivo
  (`il6R@BRAIN_E15`, celdas 3/7/3/3 < 5).
- Reparto de ramas: **anova3 = 5, hc3 = 6, art = 7**.
- Interacción `SEXO×TTO` p<.05 en **7 gen×tejido, todos en BRAIN_E15**
  (fatcd36, fatp1, fatp4, glut1, gp130, slc38a1, slc38a2); ninguna en placenta
  (mín. gp130@PLA p=.085). Post hoc corrido en esos 7.
- **D12**: dentro de BRAIN los **7/7** siguen p_BH<.05; en placenta no había
  ninguna cruda significativa. **Ninguna conclusión cambia por la columna BH.**
- `il6@BRAIN_E15` (D7): detectado solo bajo LPS (♀ 4/9, ♂ 3/9), 0/9 en ambos
  Control. Fisher ♀ p=8.24e-02, ♂ p=2.06e-01 → **0/2 alcanzan p<.05**.
- (Sobre sintético: 18 anova3, 0 hc3, 0 art; D7 y descriptivo_n_bajo iguales;
  4 post hoc. Los números cambian con la fuente pero R==Python para cada fuente.)

### Paridad y robustez verificadas

- **R ↔ Python byte-idénticos** en `qpcr_modelos_clasificacion.csv` (× R/ y
  python/), `qpcr_modelos_posthoc.csv`, `qpcr_modelos_nomodelo_descriptivo.csv`,
  `qpcr_il6_brain_fisher.csv`, `qpcr_il6_brain_tabla2x4.csv`,
  `qpcr_modelos_reporte.md`, `procedencia.csv`, `verificaciones.csv`,
  `analisis_descartados.md`, y `data/processed/qpcr_modelos_clasificacion.tsv`
  — **sobre datos reales y sobre sintético**.
- Cruza-verificación R vs `car`/`emmeans`/`ARTool`: peor |dif| = 4.2e-14 (real),
  7.9e-14 (sintético) — muy por debajo de 1e-6 (`stopifnot`).
- Idempotente: 2ª corrida de cada implementación no cambia ningún archivo;
  secuencia completa `R 02→05` seguida de `python 02→05` deja los merge files
  (`verificaciones/procedencia/analisis_descartados`) byte-idénticos.
- Truco de paridad: todo estadístico / p dependiente de trascendentes
  (`pf`, `pt`, `lgamma`, W de Shapiro) se guarda como texto `%.6e` (igual que
  T3); conteos, n y aritmética exacta con `%.10g`. Como el núcleo es el mismo
  código en ambos lenguajes, en la práctica salen byte-idénticos y no solo
  dentro de tolerancia.
- `art.con` (ARTool) para el `SEXO:TTO` de un 2×2 se reduce a
  `rank(round(y − media_global, 8))` → modelo de una vía sobre el factor de 4
  celdas → contrastes t con MSE combinado (gl = n−4); verificado numéricamente
  contra `ARTool::art.con` (rel < 5e-14).

### Entorno

- **Instalados en renv** (estaban en toolchain §8, no en la lock mínima):
  `car` 3.1-5, `ARTool` 0.11.2, `emmeans` 2.0.4, `sandwich` 3.1-3,
  `lmtest` 0.9-40 + deps (lme4, pbkrtest, quantreg, doBy, broom, dplyr, …).
  `renv::snapshot()` corrido → `renv.lock` pasa de 38 a **104 paquetes**;
  `renv::status()` limpio.

### Cierre de T5

- Tabla de estado en AGENTS.md/CLAUDE.md: **T5 → HECHO (2026-09-08)**.
- **T5 cerrado. Conviene cerrar la sesión acá.**

### Siguiente paso concreto

- **T6** (`06_pstat3` + `07_figuras_acto1`): modelo pSTAT3
  `PSTAT3 ~ SEXO * TTO + MEMBRANA` (D9, misma cascada D5 y mismo post hoc D6;
  MEMBRANA = bloque fijo), con la limitación obligatoria del informe (abundancia
  de fosfo-STAT3, no fracción). Después `07_figuras_acto1`: boxplots de
  expresión por gen×tejido (FC = 2^(−ΔΔCt), eje Y log), boxplot de pSTAT3 y
  figuras de ELISA reunidas, con **una sola función de brackets D11** (una en R,
  una en Python) usada en todas. Anotar solo si interacción y post hoc
  significativos. T6 **sí** cierra sesión.

### Notas para la próxima sesión

- El modelo de T5 reusa el núcleo hand-rolled (`anova3_terminos`, `hc3_terminos`,
  `levene_bf`, `art_terminos`, `emmeans_pares`, `holm`, `bh`, `fila_diseno`,
  `resolver`, `invertir`): para pSTAT3 en T6 basta agregar la columna `MEMBRANA`
  al diseño (bloque fijo, 2 columnas contr.sum para 3 niveles) y repetir la
  cascada. Conviene evaluar factorizar ese núcleo a un módulo compartido en T10.
- La función de brackets D11 (T7) necesita, por figura, saber si la interacción
  fue significativa y el `p_holm` de la comparación concreta: sale de
  `qpcr_modelos_posthoc.csv` (`contraste`, `p_holm`) + `interaccion_significativa`
  de `qpcr_modelos_clasificacion.csv`.
- `run_all.ps1` ya lista `05_qpcr_modelos`. Las figuras las escriben ambas
  implementaciones a la misma ruta; en `run_all` queda la de Python.

---

## Sesión 7 — 2026-09-09 — T6

### Decisiones confirmadas con el usuario (fuera del menú fijado en AGENTS)

- **Boxplots de expresión**: un PNG por tejido, faceteado por gen (no un PNG por
  gen×tejido).
- **Figura de pSTAT3**: valores **crudos** por SEXO×TTO, con la MEMBRANA como
  **forma de punto** (no residualizar por membrana).
- **il6 @ BRAIN_E15** en el Acto 1: panel de **proporción de detección** Control
  vs LPS por sexo (desde `qpcr_il6_brain_tabla2x4.csv` + p de Fisher), no un
  boxplot de FC (no es cuantificable, D7).
- **Rama `art` de la cascada D5 para pSTAT3**: `ARTool::art()` **rechaza**
  `y ~ SEXO*TTO + MEMBRANA` (`parse.art.formula` exige diseño completamente
  cruzado). Fallback fijado: **ART aditivo hand-rolled** (alinear cada término
  restando el ajuste OLS de los demás términos —incl. efecto principal de
  MEMBRANA—, rankear, ANOVA III sobre los rangos con el diseño completo; post hoc
  ART-C hand-rolled). Sobre datos **reales y sintéticos** la cascada cae en
  `anova3` → **la rama `art` no se ejecuta**; queda definida y auto-verificada
  (R recomputa el alineado por dos vías con `stopifnot`).

### Qué se completó

- **`python/06_pstat3.py` + `R/06_pstat3.R`** (equivalentes). Leen
  `data/processed/pstat3_long.tsv`. Una fila de clasificación (pSTAT3 @
  PLACENTA_E15). Modelo D9 `PSTAT3 ~ SEXO * TTO + MEMBRANA`, OLS contr.sum,
  **diseño de 6 columnas** `[1, s, t, s*t, m1, m2]` (MEMBRANA = bloque fijo, 3
  niveles → 2 columnas contr.sum), `df_resid = n − 6`.
  - **Núcleo numérico PROPIO** generalizado a `p = 6` (Gauss-Jordan, SS tipo III
    por comparación de modelos, sándwich HC3 con Wald de `k` gl por término,
    Levene BF sobre las 4 celdas SEXO×TTO, ART y ART-C aditivos, Holm).
  - **Cascada D5** (α = 0.05): Shapiro-Wilk sobre los **residuos del modelo
    conjunto** (con MEMBRANA); Levene BF sobre las 4 celdas. `anova3` si ambos
    ≥ .05; `hc3` si falla Levene; `art` (aditivo hand-rolled) si falla Shapiro.
  - **Post hoc D6** (solo si interacción `SEXO×TTO` p < .05): 4 comparaciones
    fijas, Holm. `anova3`/`hc3` → contraste de medias marginales promediando
    sobre MEMBRANA (vcov OLS / HC3); `art` → ART-C hand-rolled.
  - **F_MEMBRANA / p_MEMBRANA** se reportan como diagnóstico del bloque (2 gl),
    no gatean nada.
  - **Limitación obligatoria D9** escrita en `pstat3_reporte.md` (§6) y en
    `analisis_descartados.md`: pSTAT3 normalizado a proteína total sin STAT3
    total → **abundancia de fosfo-STAT3, no fracción fosforilada**.
- **`python/07_figuras_acto1.py` + `R/07_figuras_acto1.R`**. **UNA función de
  anotación D11** (`d11_anotacion(p, interaccion_sig)`, idéntica en R y Python)
  + `anotar_comparaciones` / `brackets_df` (una por lenguaje), usada en las 3
  figuras. Se anota una comparación solo si la interacción del gen×tejido es
  significativa **y** su `p_holm` cruza el umbral D11
  (`***`<.001, `**`<.01, `*`<.05 sólido; `.05≤p<.1` punteado + `p = 0.NNN`;
  `p≥.1` nada). 4 comparaciones D6 → pares de cajas (0=HC,1=HL,2=MC,3=ML).
  - `acto1_expresion_PLACENTA_E15.png`: 10 paneles FC = 2^(−ΔΔCt), eje log,
    4 cajas SEXO×TTO + puntos, línea en FC = 1. Sin brackets (placenta no tiene
    interacciones significativas, T5).
  - `acto1_expresion_BRAIN_E15.png`: 9 paneles FC + panel de **detección** de
    il6 (barras % detectado + p de Fisher). Brackets D11: `fatp1` (♀C–♀L `*`,
    ♀L–♂L `*`), `fatp4` (♀L–♂L tendencia `p = 0.074`), `glut1` (♀C–♀L `***`,
    ♀L–♂L `*`), `slc38a2` (♀C–♀L `**`, ♀L–♂L `**`).
  - `acto1_pstat3.png`: 4 cajas SEXO×TTO, PSTAT3 **crudo** (eje lineal), forma
    de punto = membrana; brackets `pSTAT3` ♀C–♀L `***`, ♀L–♂L `***`.
  - Las figuras de ELISA (`acto1_elisa_ms.png`, `acto1_elisa_la.png`) las hizo
    T3; 07 **no las regenera**, solo completan el set del Acto 1.
  - R: figura de BRAIN compuesta con `cowplot::ggdraw` + `draw_plot` (el panel
    de detección de il6 va en el hueco libre de la 3ª fila del facet).

### Salidas

- `data/processed/pstat3_modelo_clasificacion.tsv` (1 fila).
- `outputs/tables/{R,python}/`: `pstat3_modelo_clasificacion.csv` (F/p de
  SEXO/TTO/SEXOxTTO/MEMBRANA, rama, flags), `pstat3_posthoc.csv` (4 filas si hay
  interacción), `pstat3_descriptivo.csv` (n/media/sd/mediana por grupo).
- `outputs/tables/pstat3_reporte.md`; secciones `06_pstat3` y `07_figuras_acto1`
  en `analisis_descartados.md`; filas en `procedencia.csv` y `verificaciones.csv`.
- `outputs/figures/acto1_expresion_{PLACENTA_E15,BRAIN_E15}.png`,
  `acto1_pstat3.png` (300 dpi).

### Resultados (datos reales; NO se versionan)

- pSTAT3: rama **`anova3`** (Shapiro p residuos = 0.1006, Levene p = 0.2686).
  SEXO p = 1.98e-02 · TTO p = 9.30e-06 · **SEXO×TTO p = 1.68e-05** · MEMBRANA
  (bloque) p = 4.14e-05. Interacción significativa → post hoc D6 (OLS marginal,
  Holm): **♀Control–♀LPS `p_holm` = 1.27e-07 (`***`)**, ♂Control–♂LPS 0.883,
  **♀LPS–♂LPS 2.55e-05 (`***`)**, ♀Control–♂Control 0.141. Medias: ♀Control 1.98,
  ♀LPS 4.81, ♂Control 2.70, ♂LPS 2.75 → el efecto LPS sobre pSTAT3 es
  **restringido a hembras**.
- (Sobre sintético: rama `anova3`; SEXO×TTO p = 0.153 → **sin post hoc**;
  0 brackets en las figuras. R == Python para cada fuente.)

### Paridad y robustez verificadas

- **R ↔ Python byte-idénticos** en `pstat3_modelo_clasificacion.csv` (× R/ y
  python/), `pstat3_posthoc.csv`, `pstat3_descriptivo.csv`, `pstat3_reporte.md`,
  `data/processed/pstat3_modelo_clasificacion.tsv`, y las filas nuevas de
  `procedencia.csv` / `verificaciones.csv` / `analisis_descartados.md` (secciones
  06 y 07) — **sobre datos reales y sobre sintético**. Idempotente.
- **Cruza-verificación (solo `06_pstat3.R`)** contra `car::Anova(type=3)`,
  `car::Anova(white.adjust="hc3")` y `emmeans` (vcov OLS / HC3): peor `|dif|` =
  **9.99e-16** (real), 6.44e-15 (sintético) — muy por debajo de 1e-6
  (`stopifnot`). ART/ART-C: auto-verificación (dos vías de alineado) porque
  ARTool no ajusta el bloque aditivo; no se ejecuta con estos datos.
- Truco de paridad igual que 03/05: estadísticos / p que pasan por trascendentes
  (`pf`, `pt`, W de Shapiro) → texto `%.6e` (`p6e`); conteos y aritmética exacta
  → `%.10g`. Núcleo idéntico en ambos lenguajes → en la práctica byte-idénticos.
- La decisión D11 (`d11_anotacion`) produce **exactamente los mismos brackets**
  en R y Python (7 en expresión + 2 en pSTAT3 sobre datos reales; 0 + 0 sobre
  sintético) — verificado por la fila `figuras_d11_gate` de `verificaciones.csv`.
- Figuras PNG: equivalentes, no byte-idénticas (ggplot2 vs matplotlib).

### Entorno

- **`cowplot` 1.2.0** ya estaba en `renv.lock` (entró como dependencia en el
  snapshot de T5); `renv::status()` limpio, **no hizo falta re-snapshot**.
  Python: `matplotlib` ya en `requirements.txt`.
- `run_all.ps1` ya listaba `06_pstat3` y `07_figuras_acto1` en `$Steps` (T0):
  sin cambios.

### Cierre de T6

- Tabla de estado en AGENTS.md/CLAUDE.md: **T6 → HECHO (2026-09-09)**.
- **T6 cerrado. Conviene cerrar la sesión acá.**

### Siguiente paso concreto

- **T7** (`08_acto2_correlaciones`): sobre el score compuesto de T4
  (`qpcr_score_compuesto_long.tsv`) y/o los `neg_ddCt` por gen, correlaciones
  **placenta ↔ cerebro** por gen (y del score) con pair plots. Acto 2.1–2.2.
  **Prohibición 4**: no reportar "significativo en Control y no en LPS" como
  prueba de diferencia — para comparar correlaciones entre grupos hace falta un
  test formal (Fisher z / interacción de pendientes / permutación) y **ese** es
  el que se reporta (eso es T8). T7 **sí** cierra sesión.

### Notas para la próxima sesión

- El emparejamiento placenta–cerebro es **por feto** (`FETO`): cada uno tiene una
  fila `PLACENTA_E15` y una `BRAIN_E15` en `qpcr_cuantificacion_long.tsv` y en
  `qpcr_score_compuesto_long.tsv`. il6 no tiene `neg_ddCt` en cerebro (D7) → queda
  fuera de las correlaciones por gen de ese tejido.
- `GGally` (toolchain AGENTS §8) **no** está en la lockfile; instalarlo al abrir
  T7 si se usa para los pair plots, y `renv::snapshot()`. Python: `seaborn` ya
  está.
- El núcleo hand-rolled de 05/06 (`ajustar`, `resolver`, `anova3_terminos`,
  `emmeans_pares`, `holm`, …) puede reusarse si T8 modela pendientes; evaluar
  factorizarlo a un módulo compartido en T10 (ya anotado en T4/T5).

---

## Sesión 8 — 2026-09-09 — T7

### Decisiones confirmadas con el usuario (fuera del menú fijado en AGENTS)

- **Qué correlacionar**: `-ΔΔCt` por gen (9 genes: todos menos `il6`, que no tiene
  `-ΔΔCt` en cerebro por D7) **y** el score compuesto de 7 transportadores (D8).
- **Coeficiente**: **Spearman ρ** (no Pearson). Robusto a outliers de qPCR.
- **Estratos**: `GLOBAL` (n≤36) + por `TTO` (`CONTROL` / `LPS`, n≤18). Las 4
  celdas SEXO×TTO (n~9) NO se usan (IC inútiles).
- **Pair plots**: SPLOM de co-expresión de los 7 transportadores por tejido +
  grilla de dispersión placenta↔cerebro por gen (y score).

### Qué se completó

- **`python/08_acto2_correlaciones.py` + `R/08_acto2_correlaciones.R`**
  (equivalentes). Leen `qpcr_cuantificacion_long.tsv` + `qpcr_score_compuesto_long.tsv`.
  **T7 SOLO describe**: no compara correlaciones entre grupos (prohibición 4) ni
  las interpreta como coordinación biológica sin descartar restricción de rango
  (prohibición 5) — eso es **T8**. El reporte y `analisis_descartados.md` lo dicen
  explícitamente.
  - **Emparejamiento por FETO**: un par entra si el feto tiene `-ΔΔCt` (o score)
    detectado en **ambos** tejidos (PLACENTA_E15 y BRAIN_E15).
  - **Núcleo Spearman PROPIO** (idéntico R/Python): `ρ` = Pearson sobre rangos
    promedio (corrige empates), con sumas de acumulador `double` explícito (mismo
    orden). `p` por t-aproximación `t = ρ·√((n-2)/(1-ρ²))`, df = n-2, 2 colas —
    **la misma fórmula que `cor.test(method="spearman", exact=FALSE)` y
    `scipy.stats.spearmanr` por defecto**. IC 95% Bonett-Wright:
    `SE_z = √((1+ρ²/2)/(n-3))`, `z = atanh(ρ)`, `IC = tanh(z ± 1.959963984540054·SE_z)`.
    Piso: `n_par < 5` → solo `n` (sin ρ/IC/p).
  - **`08_acto2_correlaciones.R` cruza-verifica** cada ρ y p contra
    `cor.test(..., exact=FALSE)` (`stopifnot`, tol 1e-9). **Nunca** se llama al
    método exacto de `cor.test` (AS 89): **segfaultea** en este build de R
    (documentado en `analisis_descartados.md`).
  - **Co-expresión**: ρ de Spearman entre los 7 transportadores dentro de cada
    tejido (21 pares × 2 tejidos) → tabla + SPLOM.

### Salidas

- `outputs/tables/{R,python}/acto2_correlaciones.csv` (30 filas: 9 genes + score
  × {GLOBAL, CONTROL, LPS}; `n_par`, `rho_spearman`, `ic95_low/high`, `p_valor`).
- `outputs/tables/{R,python}/acto2_coexpresion_transportadores.csv` (42 filas).
- `outputs/tables/acto2_correlaciones_reporte.md`; sección `08_acto2_correlaciones`
  en `analisis_descartados.md`; filas en `procedencia.csv` y `verificaciones.csv`.
- `outputs/figures/acto2_dispersion_placenta_cerebro.png` (10 paneles, ρ global
  + IC y ρ por TTO anotados, coloreado por TTO, con el caveat de prohibición 4 en
  el subtítulo); `acto2_coexpresion_SPLOM_{PLACENTA_E15,BRAIN_E15}.png` (7×7,
  GGally::ggpairs con ρ overall/Control/LPS en el triángulo superior). 300 dpi.

### Resultados (datos reales; NO se versionan)

- Correlación placenta↔cerebro **GLOBAL** (Spearman ρ): score compuesto ρ = 0.24
  [-0.10, 0.53] p = 0.17; por gen todas débiles-moderadas y **solo `slc38a1`
  alcanza p < .05** (ρ = 0.47 [0.13, 0.71] p = 6.4e-3); `gp130` ρ = 0.32 p = .058.
- Por TTO (descriptivo, **NO se compara**): en varios genes y en el score la ρ
  de Control es mayor que la de LPS (p. ej. score Control ρ = 0.48 vs LPS
  ρ = -0.14; fatcd36 0.47 vs -0.13; fatp4 0.42 vs -0.11). **Interpretar esa
  diferencia es T8** (test formal + restricción de rango).
- Co-expresión intra-tejido de los transportadores: fuerte en ambos tejidos
  (ρ ~ 0.5–0.87), contexto para la restricción de rango de T8.
- (Sobre sintético: R == Python para cada fuente; los números cambian.)

### Paridad y robustez verificadas

- **R ↔ Python byte-idénticos** en `acto2_correlaciones.csv`,
  `acto2_coexpresion_transportadores.csv` (× R/ y python/),
  `acto2_correlaciones_reporte.md`, y las filas nuevas de `procedencia.csv` /
  `verificaciones.csv` / `analisis_descartados.md` — **sobre datos reales y
  sobre sintético**. Idempotente.
- Cruza-verificación R vs `cor.test(exact=FALSE)`: peor `|dif|` = **8.9e-16**
  (real) — muy por debajo de 1e-9 (`stopifnot`). Además el núcleo PROPIO de
  Python reproduce `scipy.stats.spearmanr` a `%.10g` / `%.6e`.
- Truco de paridad: ρ y sumas con acumulador `double` explícito → `%.10g`
  bit-idéntico; `p` (vía `pt`) e IC (vía `tanh`/`atanh`) → texto `%.6e`.
- Figuras PNG: equivalentes, no byte-idénticas (ggplot2/GGally vs matplotlib).

### Entorno

- **`GGally` 2.4.0** instalado en renv (+ `forcats` 1.0.1, `ggstats` 0.14.0,
  `patchwork` 1.3.2). `renv::snapshot()` → `renv.lock` pasa a **111 paquetes**;
  `renv::status()` limpio. Python: `seaborn` / `scipy` ya estaban.
- `run_all.ps1` ya listaba `08_acto2_correlaciones`: sin cambios.

### Cierre de T7

- Tabla de estado en AGENTS.md/CLAUDE.md: **T7 → HECHO (2026-09-09)**.
- **T7 cerrado. Conviene cerrar la sesión acá.**

### Siguiente paso concreto

- **T8** (`09_acto2_dispersion` + `10_acto2_simulacion` + test de pendientes):
  Acto 2.3–2.5. (a) Comparar formalmente las correlaciones placenta↔cerebro
  entre Control y LPS con un test que **es el que se reporta** (Fisher z sobre
  Spearman, interacción de pendientes en un modelo, o permutación — elegir con el
  usuario). (b) **Antes** de interpretar cualquier cambio de correlación como
  cambio de coordinación biológica, descartar cambio de dispersión / restricción
  de rango: `10_acto2_simulacion` (Sección 8.4 del brief) genera datos con la
  misma correlación verdadera pero distinta dispersión por grupo y muestra cuánto
  se mueve la ρ observada. T8 **sí** cierra sesión.

### Notas para la próxima sesión

- La entrada de T8 son los mismos pares por feto de T7 (`_pares(D, item, estrato)`
  en Python / `pares()` en R) — conviene factorizar esa construcción y el núcleo
  Spearman a un helper compartido si T8 los reusa mucho.
- Para el Fisher z sobre Spearman: usar `z = atanh(ρ)`, `SE = 1/√(n-3)` (o la
  variante Bonett-Wright ya implementada en `spearman_ci`), y el estadístico
  `(z1 - z2)/√(SE1² + SE2²)`. Confirmar con el usuario cuál de las tres vías
  (Fisher z / pendientes / permutación) es **la** que se reporta (prohibición 4).
- `il6R@BRAIN` tiene `n_par` chico por estrato (Control 6, LPS 9): el test de
  diferencia para ese gen tendrá potencia casi nula — dejarlo explícito, no
  omitirlo en silencio.

---

## Sesión 9 — 2026-09-10 — T8

### Decisiones confirmadas con el usuario (fuera del menú fijado en AGENTS)

- **Test reportado de diferencia de correlaciones (prohibición 4)**: **Fisher z
  sobre ρ de Spearman**. `z = atanh(ρ)`; estadístico
  `(z_control − z_lps)/√(SE_control² + SE_lps²)`; `p` normal a dos colas.
- **Error estándar del Fisher z**: **Bonett-Wright primario**
  `SE_i = √((1 + ρ_i²/2)/(n_i − 3))` (el mismo del IC de T7) **+** SE clásico
  `1/√(n_i − 3)` como columna al lado (no cambia conclusiones). Además `p_bw_bh`
  (BH entre ítems) como columna suplementaria en el espíritu de D12.
- **Simulación de restricción de rango (prohibición 5)**: **normal bivariada en
  −ΔΔCt**. `r` verdadera común a ambos grupos, SD marginales = SD observadas por
  grupo. Inversión `r = 2·sin(π·ρ_S/6)` (exacta para la normal bivariada).
  **Dos escenarios de `r` verdadera como rango**: `GLOBAL` (ρ de T7, pooled) y
  `CONTROL` (ρ de Control); el Δρ observado se compara contra la distribución
  simulada de cada uno. B = 2000, semilla 20260101, RNG PROPIO.

### Qué se completó

- **`python/09_acto2_dispersion.py` + `R/09_acto2_dispersion.R`** (Acto 2.3–2.4).
  - **2.4 test reportado (prohibición 4)**: Fisher z sobre ρ de Spearman
    Control vs LPS, sobre los mismos pares por feto de T7. Piso: se testea solo
    si Control **y** LPS tienen `n_par ≥ 5`. Salida
    `acto2_test_correlaciones.csv` (10 filas: 9 genes + score): `rho_control`,
    `rho_lps`, `delta_rho`, `z_*`, `se_bw_*`, `stat_z_bw`, `p_bw`, `p_bw_bh`,
    `se_clasico_*`, `stat_z_clasico`, `p_clasico`.
  - **2.3 dispersión (insumo de la prohibición 5)**: `acto2_dispersion.csv`
    (20 filas: ítem × tejido): `sd_control`, `sd_lps`, `ratio_var_lps_control`,
    Levene Brown-Forsythe por lado (`levene_bf_F`, `levene_bf_p`). SD (n−1) sobre
    los mismos pares por feto que la correlación.
  - **09.R cruza-verifica** cada ρ contra `cor.test(exact=FALSE)` (`stopifnot`
    tol 1e-9) y cada `F` de Levene contra `car::leveneTest(center=median)`
    (tol 1e-8). Peores |dif| reales: ρ = 1.1e-16, F Levene = 5.1e-15.
- **`python/10_acto2_simulacion.py` + `R/10_acto2_simulacion.R`** (Acto 2.5,
  prohibición 5). `acto2_simulacion.csv` (20 filas: ítem × {GLOBAL, CONTROL}):
  `rho_true`, `r_pearson_gen`, SD por grupo/lado, `delta_rho_obs`,
  `sim_mean_delta_rho`, `sim_sd_delta_rho`, `sim_q025`, `sim_q975`,
  `sim_frac_abs_ge_obs` (p de simulación), `sim_tasa_fisher_sig` (falsos
  positivos del Fisher z de 09 bajo pura restricción de rango), `veredicto`
  DENTRO/FUERA. Núcleo PROPIO: `cuantil_tipo7` (== `quantile` default de R),
  `rho_s_a_r`, sorteos con `norm1()` del RNG de `00_config` (un único RNG
  recorre ítems × escenarios en orden).

### Salidas

- `outputs/tables/{R,python}/acto2_dispersion.csv`,
  `acto2_test_correlaciones.csv`, `acto2_simulacion.csv`.
- `outputs/tables/acto2_dispersion_reporte.md`,
  `acto2_simulacion_reporte.md`; secciones `09_acto2_dispersion` y
  `10_acto2_simulacion` en `analisis_descartados.md`; filas en
  `procedencia.csv` (8) y `verificaciones.csv` (17).
- `outputs/figures/acto2_dispersion_sd.png`, `acto2_test_delta_rho.png`
  (forest ρ_control vs ρ_lps con Δρ y `p_bw`), `acto2_simulacion_delta_rho.png`
  (Δρ observado vs IC95 simulado por ítem y escenario). 300 dpi.

### Resultados (datos reales; NO se versionan)

- **Test formal (prohibición 4): NINGÚN ítem alcanza `p < .05`.** El más chico
  es `score_compuesto` `p_bw = 0.079` (Δρ = 0.62), luego `fatcd36` `p_bw = 0.107`
  (Δρ = 0.60) y `fatp4` `p_bw = 0.140` (Δρ = 0.53). Con BH entre ítems, el
  mínimo `p_bw_bh = 0.47`. → **No hay evidencia formal de que las correlaciones
  placenta↔cerebro difieran entre Control y LPS**; el patrón descriptivo de T7
  (Control coordinado, LPS desacoplado) no se sostiene como diferencia.
- **Simulación (prohibición 5): 18/20 celdas DENTRO del IC95 simulado.** Solo
  `fatcd36` y `score_compuesto` caen FUERA, y únicamente bajo el escenario
  `CONTROL` (más exigente), con `sim_frac_abs_ge_obs` ≈ 0.03–0.04 (bajo `GLOBAL`
  están DENTRO). → **No se puede descartar restricción de rango**: la sola
  diferencia de dispersión por grupo produce Δρ del tamaño observado.
- La tasa de falsos positivos del Fisher z bajo pura restricción de rango es
  ≈ 0.032–0.047 (bien calibrada al 5% nominal).
- (Sobre sintético: R == Python byte a byte; los números cambian.)

### Paridad y robustez verificadas

- **R ↔ Python byte-idénticos** en `acto2_dispersion.csv`,
  `acto2_test_correlaciones.csv`, `acto2_simulacion.csv` (× R/ y python/),
  los dos `*_reporte.md`, y las filas nuevas de `procedencia.csv` /
  `verificaciones.csv` / `analisis_descartados.md` — **sobre datos reales y
  sobre sintético** (probado moviendo `data/raw/` y reprocesando 02→10).
  Idempotente (R y Python, 2 corridas cada uno).
- Truco de paridad: ρ, SD, cociente de varianzas, medias/SD/cuantiles de la
  simulación y todos los conteos → acumulador `double` explícito, texto
  `%.10g` bit-idéntico. `z` (atanh), `p` (normal), `F` de Levene (pf) y
  `r_pearson_gen` (sin) → texto `%.6e`. Los sorteos de la simulación usan
  `norm1()` (LCG + polar de Marsaglia de T1): byte-idéntico R/Python.
- Figuras PNG: equivalentes, no byte-idénticas (ggplot2 vs matplotlib;
  R usa ASCII en labels, Python usa Unicode).

### Entorno

- Sin dependencias nuevas: `09` usa `car` (ya en renv por 05) + `ggplot2`;
  `10` usa `ggplot2`. `renv.lock` sin cambios (111 paquetes).
- `run_all.ps1` ya listaba `09_acto2_dispersion` y `10_acto2_simulacion`:
  sin cambios.

### Cierre de T8

- Tabla de estado en AGENTS.md/CLAUDE.md: **T8 → HECHO (2026-09-10)**.
- **T8 cerrado. Conviene cerrar la sesión acá.**

### Siguiente paso concreto

- **T9** (`11_sensibilidad`, Acto 2.6): controles de sensibilidad de la
  correlación / del test — (a) **eigengene** (variante PCA del score compuesto,
  D8) en lugar del score promedio de z, rehaciendo la correlación
  placenta↔cerebro y el Fisher z; (b) **exclusión del feto extremo** (el más
  influyente por distancia de Cook o por |Δρ| al quitarlo, uno por ítem) y ver
  si mueve ρ y el veredicto. Doble implementación R/Python byte-idéntica.
  Confirmar con el usuario el criterio exacto de "feto extremo" (Cook vs
  leave-one-out sobre ρ) antes de implementar. T9 cierra sesión.

### Notas para la próxima sesión

- El núcleo Spearman + `pares()` + `cargar()` está ahora replicado en 08, 09 y
  10 (idéntico). Si T9 lo vuelve a usar, sigue siendo preferible re-inlinear
  (la convención del repo es que cada script importe solo `00_config`) antes que
  tocar `00_config` y arriesgar sus salidas.
- El eigengene de T9 necesita PCA PROPIO (no `prcomp`/`sklearn` en el resultado)
  para paridad byte a byte: eigendescomposición de la matriz de correlación
  7×7 por Jacobi, o SVD hand-rolled. Verificar signo del PC1 (convención:
  cargas mayoritariamente positivas) igual que en `informe-e15-reproducible`.
- `10_acto2_simulacion` tarda ~15 s en R (B=2000 × 20 celdas × ~30 `norm1`).
  Si T9 simula encima, considerar bajar B o reusar las corridas de 10.
- La conclusión sustantiva del Acto 2 para el informe: la coordinación
  placenta↔cerebro **no cambia de forma demostrable** entre Control y LPS —
  ni el test formal (prohibición 4) ni la simulación (prohibición 5) permiten
  afirmar una diferencia. Es un resultado negativo limpio.

---

## Sesión 10 — 2026-09-10 — T9

### Decisiones confirmadas con el usuario (fuera del menú fijado en AGENTS)

- **Feto extremo (control B)**: **leave-one-out sobre ρ**, no distancia de Cook
  (Cook mide influencia sobre un ajuste lineal; el estadístico del Acto 2 es de
  rango). **Blanco = ρ GLOBAL por ítem**: el feto de mayor
  `|ρ_full − ρ_sin_i|` sobre la ρ pooled; se lo excluye de los 3 estratos y se
  rehacen ρ + IC + Fisher z.
- **Eigengene con fetos incompletos (<7 transportadores detectados)**:
  **proyectar con los z disponibles** y reescalar por `sqrt(Σ carga_g²)` de las
  cargas usadas (análogo al "promedio de los z disponibles" de D8). No casos
  completos (descartaría 4–6 fetos y cambiaría el módulo para todos).

### Qué se completó

- **`python/11_sensibilidad.py` + `R/11_sensibilidad.R`** (equivalentes). Leen
  `qpcr_cuantificacion_long.tsv` (columna `z`, D8) + `qpcr_score_compuesto_long.tsv`.
  Dos controles de sensibilidad del Acto 2, cada uno rehaciendo la correlación
  placenta↔cerebro (Spearman ρ, GLOBAL/CONTROL/LPS) y el Fisher z Control vs LPS
  de T8:
  - **(A) Eigengene** = proyección sobre **PC1 PROPIO** del módulo de 7
    transportadores, **por tejido**. Matriz de entrada: correlación de Pearson
    pairwise-complete de las 7 columnas `z` (se documenta el descarte de la
    matriz de Spearman y de casos completos). **Jacobi clásico SIN
    trigonometría** (`t = 1/(θ ± √(θ²+1))`, `c = 1/√(t²+1)`, `s = t·c` → sólo
    `+ − × ÷ √`; barrido cíclico p<q, umbral `1e-15`, ≤100 sweeps) → **bit-idéntico
    R/Python**, mismo truco que el RNG polar de `00_config`. Signo de PC1 fijado
    a `Σ cargas > 0`. `R` cruza-verifica contra `eigen()` (tol 1e-8), Python
    contra `numpy.linalg.eigh` (tol 1e-6): peor |dif| = 4.0e-13.
  - **(B) Exclusión del feto extremo** (LOO sobre ρ GLOBAL, uno por ítem) para
    **11 ítems** (9 genes de T7/T8 + score compuesto + eigengene). Se excluye de
    los 3 estratos y se rehace ρ (IC Bonett-Wright) + Fisher z; `veredicto_cambia`
    marca si `p_bw<.05` cambia. **NO se re-estima PC1** al quitar el feto (cargas
    fijas sobre los 36). `R` cruza-verifica cada ρ de estrato vs
    `cor.test(exact=FALSE)` (tol 1e-9): peor |dif| = 1.1e-16.

### Salidas

- `outputs/tables/{R,python}/`: `acto2_sensibilidad_eigengene_correlacion.csv`
  (6 filas: eigengene + score × 3 estratos), `acto2_sensibilidad_eigengene_test.csv`
  (2 filas: Fisher z eigengene vs score), `acto2_sensibilidad_pca_loadings.csv`
  (14 filas: carga PC1 × 7 genes × 2 tejidos), `acto2_sensibilidad_pca_varianza.csv`
  (14 filas: autovalor/prop_var de los 7 PC × 2 tejidos),
  `acto2_sensibilidad_excl_extremo.csv` (11 filas).
- `outputs/tables/acto2_sensibilidad_reporte.md`; sección `11_sensibilidad` en
  `analisis_descartados.md`; **8 filas** en `procedencia.csv` y **11** en
  `verificaciones.csv` (todas TRUE).
- `outputs/figures/acto2_sensibilidad_eigengene.png` (cargas PC1 por tejido +
  dispersión eigengene + ρ score vs eigengene por estrato),
  `acto2_sensibilidad_excl_extremo.png` (ρ GLOBAL con/sin feto extremo por ítem,
  con el cambio de `p_bw`). 300 dpi.

### Resultados (datos reales; NO se versionan)

- **(A) Eigengene**: PC1 explica **76.3 %** (placenta) / **81.6 %** (cerebro) de
  la varianza, **7/7 cargas positivas** en ambos tejidos → eje de "tono de
  expresión". El eigengene **replica el resultado del score compuesto**: ρ GLOBAL
  placenta↔cerebro = 0.244 (n.s., score 0.236); Fisher z Control vs LPS
  `p_bw` = 0.109 (eigengene) vs 0.079 (score) → **ninguno alcanza `p<.05`**,
  consistente con T8.
- **(B) Exclusión del feto extremo**: **0/11 ítems cambian el veredicto**
  `p_bw<.05`. El feto más influyente para el score, el eigengene y varios genes
  es `L_091025_2.08` (LPS); excluirlo no vuelve significativo ningún test
  (p. ej. score `p_bw` 0.079 → 0.175). `il6R` sí llega al piso n≥5/grupo (Control
  6, LPS 9), se testea.
- **Nota de método**: la matriz de correlación pairwise-complete es levemente
  indefinida (7º autovalor en cerebro ≈ −0.11 sobre traza 7); PC1 domina y
  coincide con `eigen()` a 1e-13 → no afecta el eigengene. Registrado en
  `verificaciones.csv` (`sens_eigengene_matriz_psd`) y `analisis_descartados.md`.
- **Conclusión del Acto 2 reforzada**: el resultado negativo (sin diferencia
  demostrable en la coordinación placenta↔cerebro entre Control y LPS) es
  **robusto** al resumen elegido (promedio-de-z vs eigengene PC1) y al feto
  individual más influyente.
- (Sobre sintético: R == Python byte a byte; PC1 var. expl. 88.1 % / 80.9 %;
  0/11 cambian veredicto. Los números cambian con la fuente.)

### Paridad y robustez verificadas

- **R ↔ Python byte-idénticos** en los 5 `acto2_sensibilidad_*.csv` (× `R/` y
  `python/`), `acto2_sensibilidad_reporte.md`, y las filas nuevas de
  `procedencia.csv` / `verificaciones.csv` / `analisis_descartados.md` — **sobre
  datos reales y sobre sintético** (probado moviendo `data/raw/` y reprocesando
  02→11). Idempotente (R y Python, 2 corridas cada uno).
- Cadena completa **R 02→11** vs **Python 02→11**: los 3 merge files quedan
  byte-idénticos.
- Truco de paridad: Jacobi sin trig + Spearman/Pearson/SD/normas con acumulador
  `double` explícito → `%.10g` bit-idéntico para ρ, cargas, autovalores y
  eigengene; `p` (pt), IC (tanh/atanh) y `p` del Fisher z (normal) → texto
  `%.6e`. Las filas `score_compuesto` de las tablas (A) **reproducen T7/T8**
  exactamente (ρ, IC, z, SE, stat, p) — chequeo de consistencia interna.
- Figuras PNG: equivalentes, no byte-idénticas (ggplot2/patchwork vs matplotlib).

### Entorno

- Sin dependencias nuevas: `11` usa `ggplot2` + `patchwork` (ya en `renv.lock`
  desde T7) y, sólo en la cruza-verificación, `eigen()` base. Python:
  `numpy` (ya en `requirements.txt`). `renv.lock` sin cambios (111 paquetes).
- `run_all.ps1` ya listaba `11_sensibilidad` en `$Steps` (T0): sin cambios.
  **Nota para T11**: `run_all.ps1` tiene un bug preexistente — bajo
  `Set-StrictMode -Version Latest`, `(Get-Command Rscript.exe -EA
  SilentlyContinue).Source` explota si el comando no está. Arreglar al finalizar
  `run_all` en T11 (guardar el resultado de `Get-Command` antes de `.Source`).

### Cierre de T9

- Tabla de estado en AGENTS.md/CLAUDE.md: **T9 → HECHO (2026-09-10)** (md5 de
  ambos idéntico tras el cambio).
- **T9 cerrado. Conviene cerrar la sesión acá.**

### Siguiente paso concreto

- **T10** (tablas de auditoría): consolidar `procedencia.csv`,
  `verificaciones.csv`, `analisis_descartados.md` y generar
  `comparacion_R_python.csv` (concordancia numérica automática cruzando
  `outputs/tables/R/` vs `outputs/tables/python/` archivo por archivo, con las
  tolerancias declaradas 1e-6 / 1e-4). Evaluar factorizar a un módulo compartido
  el núcleo hand-rolled (`spearman_rho`, `pares`/`cargar`, merges) hoy replicado
  en 08–11. T10 cierra sesión.

### Notas para la próxima sesión

- El núcleo Spearman + `pares()` + `cargar()` está replicado idéntico en 08, 09,
  10 y 11; `11` agrega Jacobi sin trig + `pearson_pairwise`. Si T10 factoriza,
  mantener la regla del repo (cada script importa sólo `00_config`) o mover el
  núcleo a `00_config` con mucho cuidado de no tocar sus salidas.
- `comparacion_R_python.csv` de T10 debería cruzar TODOS los
  `outputs/tables/{R,python}/*.csv` (ya son byte-idénticos por construcción en
  02–11; el archivo lo deja documentado y lo chequea `99_verificar`).
- La conclusión del Acto 2 para el informe (T11) ya está firme: negativo limpio
  y **robusto** (T9). Mencionar var. expl. de PC1 (76/82 %) y que el eigengene
  no cambia nada.

---

## Sesión 11 — 2026-09-10 — T10

### Decisiones tomadas (ninguna toca D1–D12 ni salidas existentes)

- **Script nuevo `98_comparacion` (R y Python)**, numerado 98 = utilitario
  hermano de `99_verificar`. En `run_all.ps1` va **después de `11_sensibilidad`
  y antes de `12_informe`** (el informe podrá citar sus tablas). Ambos lenguajes
  producen EL MISMO `comparacion_R_python.csv` + `comparacion_reporte.md` (control
  cruzado uno del otro).
- **No se factoriza** el núcleo hand-rolled de 08–11 (`spearman_rho` /
  `rangos_promedio` / `pares()` / `cargar()`). Se evaluó mover a `00_config` /
  módulo compartido y se descartó: rompe la regla "cada script importa sólo
  `00_config`", arriesga perturbar salidas byte-idénticas ya cerradas justo antes
  del informe, y la réplica ya está verificada byte a byte + cruzada contra
  `cor.test`/`car`/`eigen()` en cada script. Registrado en
  `analisis_descartados.md` (sección `98_comparacion`). Revisable post-T11.
- `procedencia.csv` columna `entradas` → `data/real/…` **no es un bug**: es
  `data/<fuente>/<archivo>` con `fuente ∈ {real, sintetico}` (de
  `cfg.fuente_datos`), etiqueta de procedencia, no ruta de FS. Sin cambios.

### Qué se completó

- **`python/98_comparacion.py` + `R/98_comparacion.R`** (equivalentes):
  1. **Concordancia numérica R↔Python** archivo por archivo de los 31 CSV de
     `outputs/tables/{R,python}/`. Cada celda que parsea a número finito en ambos
     lados: `|a−b| ≤ 1e-6` **o** `|a−b|/max(|a|,|b|) ≤ 1e-6` (columna `tol` pasa a
     `1e-4` si alguna celda necesita el margen laxo — reservado a `p` de tests
     iterativos); el resto, texto exacto. Salida
     `outputs/tables/comparacion_R_python.csv` (31 filas + fila `__TOTAL__`),
     columnas: `filas_R/py`, `cols_R/py`, `header_igual`, `celdas`,
     `celdas_numericas/texto`, `n_dif_texto`, `max_dif_abs/rel`, `peor_celda`,
     `tol`, `n_fuera_tol`, `byte_identico`, `ok`.
  2. **Auditoría de las 3 tablas** → 7 filas en `verificaciones.csv`
     (script `98_comparacion`): cobertura de gemelos R/python (0 huérfanos),
     headers iguales (31/31), dentro de tolerancia (0 fuera), toda figura con
     fila en `procedencia.csv` (13/13), toda tabla `{R,python}/*.csv` con fila
     (31/31), sección en `analisis_descartados.md` por script 02–11 (10/10),
     ninguna `verificacion` distinta de TRUE (79/79). El conteo de
     `verificaciones` **excluye las propias filas de `98_comparacion`** para ser
     estable ante el orden de ejecución (R antes/después de Python) y re-corridas.
  3. **Sección `98_comparacion` en `analisis_descartados.md`** con las
     justificaciones que 02/03/04 diferían a T10:
     - **Descarte de imputación MNAR** (D3), con números de la corrida: 77/720
       obs. no detectadas (10.7 %); celdas D7 (calibrador ♀Control 0/n en
       il6@BRAIN_E15). Tres razones: ancla inexistente en celdas D7 → FC contra
       valor inventado (prohibición 10); estructura artificial en el baseline
       Control (exploración previa); y no hace falta (D1/D8 sólo sobre
       detectados, D7 analiza detección).
     - **LOD alternativo del ELISA** (menor estándar de `CURVA IL6`): **no
       aplicado**. Mover el LOD de 0 al menor estándar no reordena los empates
       (censura = `Conc < 0`, no el valor del LOD) → no puede cambiar la
       conclusión del Peto-Peto/Fisher por construcción. Se documenta, no se
       corre.
     - **No factorizar** el núcleo 08–11 (arriba).
  4. **`outputs/tables/comparacion_reporte.md`**: reporte legible (tabla de las
     31 comparaciones + resultados de las 4 auditorías).
  5. **5 filas en `procedencia.csv`** (script `98_comparacion`):
     `comparacion_R_python.csv`, `comparacion_reporte.md` y las 3 tablas de
     auditoría registradas como `tipo = auditoria`.
- **`run_all.ps1`**: `'98_comparacion'` agregado a `$Steps` entre
  `11_sensibilidad` y `12_informe`.
- **AGENTS.md / CLAUDE.md** (md5 idéntico, `6bd87dee…`): §3 árbol de carpetas
  (`98_comparacion.R` en `R/`, `comparacion_reporte.md` en `tables/`); §6 tabla:
  **T10 → HECHO (2026-09-10)**.

### Resultados (datos reales; NO se versionan — `outputs/tables/*` está en .gitignore)

- **31/31 CSV byte-idénticos** R↔Python; peor |dif| absoluta = 0; 0 celdas fuera
  de tolerancia; 0 diferencias de texto; `tol = 1e-6` en todos (nunca hizo falta
  el fallback 1e-4).
- **Auditoría**: 0 figuras sin fila (13), 0 tablas sin fila (31), 0 scripts sin
  sección en `analisis_descartados.md` (10/10), 79/79 `verificaciones` en TRUE.
- `comparacion_R_python.csv`, `comparacion_reporte.md`, `procedencia.csv`,
  `verificaciones.csv` y `analisis_descartados.md` quedan **byte-idénticos**
  tras correr Python y después R (probado), y **estables ante re-corridas**
  (2×Python, 2×R → mismos md5).

### Entorno

- Sin dependencias nuevas. `98_comparacion` sólo usa base R / stdlib de Python
  (parser CSV propio, ya en 02–11). `renv.lock` sin cambios. `requirements.txt`
  sin cambios.

### Cierre de T10

- **T10 cerrado. Conviene cerrar la sesión acá.**

### Siguiente paso concreto

- **T11** (último): `R/12_informe.R` + `python/12_informe.py` (informe HTML
  autocontenido + PDF por impresión headless si no hay LaTeX; el pipeline no debe
  fallar por el PDF), `R/99_verificar.R` + `python/99_verificar.py` (chequea
  existencia y no-vacuidad de todo el checklist de AGENTS §1, re-corre la
  comparación R↔Python de `comparacion_R_python.csv`, y además **byte-compara los
  10 `.md` de copia única** de `outputs/tables/` entre una corrida R y una
  Python), completar el cuerpo real de `run_all.ps1` y dejar
  `logs/corrida_<fecha>.txt`. `99_verificar` debe imprimir
  `TODAS LAS VERIFICACIONES PASARON`.

### Notas para la próxima sesión

- **Bug preexistente de `run_all.ps1`** (ya anotado en T9): bajo
  `Set-StrictMode -Version Latest`, `(Get-Command Rscript.exe -EA
  SilentlyContinue).Source` explota si el comando no existe. Guardar el resultado
  de `Get-Command` antes de `.Source`. Arreglar al completar `run_all` en T11.
- `99_verificar` puede reusar `comparar_archivo()` / `leer_csv()` de
  `98_comparacion` (misma lógica); mantener la regla de importar sólo
  `00_config` — copiar el helper, no cross-importar entre scripts numerados.
- El informe (T11) ya tiene todo el material: Acto 1 (T3–T6), Acto 2 (T7–T9,
  conclusión negativa robusta), y `comparacion_R_python.csv` para la sección de
  reproducibilidad. Var. expl. PC1 76/82 %; MNAR descartada con números en
  `analisis_descartados.md`.

---

## Sesión 12 — 2026-09-10 — T11 (cierre del proyecto)

### Decisiones tomadas (ninguna toca D1–D12)

- **Informe a mano, sin rmarkdown/pandoc.** No están en el toolchain (`rmarkdown`
  ni siquiera se instala en este build de R) y un motor intermedio nunca daría
  salida byte-comparable. `12_informe` arma `docs/informe.html` con la misma
  lógica de strings en R y Python: un mini Markdown→HTML propio (encabezados,
  párrafos, listas con anidado, tablas `| … |`, `**negrita**`/`` `código` ``/
  `*énfasis*`/`[t](u)`), un CSS inline y `<section>` por bloque.
- **Figuras incrustadas en base64** → informe **autocontenido**, pero como los
  PNG **no** son byte-idénticos ggplot2↔matplotlib, `informe.html` tampoco lo es.
  La paridad R/Python se chequea sobre `informe.textonly.html` (los `data:` →
  `src="[png]"`). base64 en R: implementación propia RFC 4648 (no hay paquete en
  renv), misma salida que `base64.b64encode`.
- **Sin marca de tiempo en el HTML** (rompería la paridad). La fecha va en
  `logs/corrida_<fecha>.txt`.
- **`run_all.ps1` en pasadas separadas por lenguaje** (modo `both`): PRIME
  python 00–11 → R 00–11/98/12 → python 00–11/98/12 → VERIFY R/99 + python/99.
  La pasada intercalada (R, luego Python, por paso) hacía que Python pisara
  siempre los `.md` de R antes de que `12_informe` R los viera → el byte-compare
  de `.md` entre lenguajes era trivial (Python vs Python). El PRIME (sin `98`)
  puebla `outputs/tables/python/` para que `R/98` tenga contraparte en frío.
  `-Only R|python` = una pasada, sin `98`/`99`.
- **El bloque de auditoría del informe excluye las filas `script == 12_informe`**
  de `procedencia.csv`/`verificaciones.csv` (`leer_csv_sin_este`): así queda
  idéntico R/Python e idempotente entre corridas (12 no se cuenta a sí mismo).
- **`docs/*` y `logs/*` a `.gitignore`** (misma lógica que `outputs/`: regenerable
  y una corrida real dejaría números inéditos). El formato del HTML es
  "listo para GitHub Pages", el archivo no se versiona.
- **`-FromSynthetic` cableado**: `run_all.ps1` exporta
  `MIA_LPS_FORZAR_SINTETICO=1` y `00_config` (R y Python) lo respeta en
  `ruta_datos()`/`fuente_datos()` (cambio aditivo, no toca salidas).
- **`run_all.ps1` fuerza UTF-8** (`PYTHONUTF8=1`, `[Console]::OutputEncoding`):
  con la salida redirigida, `00_config.py` caía en cp1252 al imprimir `♀`.
- **Bug preexistente de `run_all.ps1` arreglado**: `Get-Command` se guarda en una
  variable antes de leer `.Source` (bajo `Set-StrictMode` explotaba si faltaba).

### Qué se completó

- **`python/12_informe.py` + `R/12_informe.R`** (equivalentes). Leen los 11 `.md`
  de copia única de `outputs/tables/`, unos CSV de `outputs/tables/<lang>/` para
  los números del resumen, y los 13 PNG de `outputs/figures/`. Producen:
  - `docs/informe.html` — 8 secciones: 1 Resumen (números parseados de los CSV),
    2 Diseño y métodos (D1–D12 + cascada + `qc_reporte.md`), 3 Acto 1 (ELISA,
    cuantificación, modelos qPCR, pSTAT3), 4 Acto 2 (correlaciones, dispersión +
    test Δρ, simulación, sensibilidad), 5 Reproducibilidad
    (`comparacion_reporte.md`), 6 Análisis descartados (`analisis_descartados.md`),
    7 Limitaciones (D9, censura ELISA, il6 cerebro, sintético≠real, BRAIN_P1,
    baja potencia Acto 2), 8 Procedencia y verificaciones (tablas embebidas).
  - `docs/informe.pdf` — Edge headless (`--print-to-pdf`), best-effort. Motor:
    env `MIA_LPS_PDF_ENGINE` o búsqueda de `msedge.exe`/`chrome.exe`. En esta
    máquina: `C:\Program Files (x86)\…\msedge.exe` → PDF de ~4.2 MB, `pdf_status=ok`.
  - `outputs/intermediate/render/<lang>/` — snapshot para la paridad de `99`:
    `informe.html`, `informe.textonly.html`, `pdf_status.txt`, `tables/*.md`.
  - Filas en `procedencia.csv` (2) y `verificaciones.csv` (5, todas TRUE).
- **`python/99_verificar.py` + `R/99_verificar.R`** (equivalentes). Compuerta
  final. 77 chequeos duros: checklist de AGENTS §1 (raíz, 15 scripts × 2, 5
  sintéticos, 13 figuras, 11 reportes, 3 tablas de auditoría, `informe.html`),
  `informe.pdf` **blando** si `pdf_status != ok`; re-corre la concordancia
  numérica de los 31 CSV `outputs/tables/{R,python}/` en proceso (helpers copiados
  de `98`); byte-compara `render/R/` vs `render/python/` (`informe.textonly.html`
  + `tables/*.md` + disco==snapshot); `verificaciones.csv` todas TRUE. Escribe
  `logs/corrida_<fecha>.txt` (una sección por lenguaje, marcadores
  `<!-- R:inicio -->` / `<!-- python:inicio -->`). Imprime
  `TODAS LAS VERIFICACIONES PASARON` y sale 0; si algo duro falla, imprime los
  fallos y sale 1.
- **`run_all.ps1`** reescrito (cuerpo real, 4 pasadas, UTF-8, bug de
  `Get-Command`, flags `-SkipVerify`).
- **`00_config.{R,py}`**: `_forzar_sintetico()` / `.forzar_sintetico()` +
  `-FromSynthetic`.
- **`.gitignore`**: `docs/*`, `logs/*`.
- **AGENTS.md / CLAUDE.md** (md5 idéntico): checklist §1 todo a `[x]`; §6
  **T11 → HECHO (2026-09-10)** + nota de la estructura de `run_all`.
- **README.md**: §2 sin pandoc; §4 pasadas separadas + qué chequea `99`; §5 PDF.

### Resultados (datos reales; NO se versionan)

- **`.\run_all.ps1` de punta a punta: ~2.9 min, ambos `99_verificar` imprimen
  `TODAS LAS VERIFICACIONES PASARON`.** 77/77 chequeos duros; 31/31 CSV
  byte-idénticos R↔Python; 0 celdas fuera de tolerancia; paridad de render OK
  (`informe.textonly.html` + los 11 `.md` byte-idénticos entre la pasada R y la
  Python); `informe.pdf` presente (~4.2 MB).
- El informe recoge los números de la corrida real: pSTAT3 SEXO×TTO
  p = 1.68e-05 (post hoc ♀ Holm p = 1.27e-07); qPCR 7/7 interacciones en cerebro
  (0/7 en placenta), las 7 sobreviven BH; Acto 2 mínimo `p_bw` = 0.079
  (`score_compuesto`), simulación 2/20 FUERA (solo CONTROL), eigengene
  `p_bw` = 0.109, 0/11 fetos extremos cambian el veredicto.

### Paridad y robustez verificadas

- `12_informe`: `informe.textonly.html` y los 11 `.md` snapshot **byte-idénticos**
  R↔Python (probado con la secuencia python→R→python que replica `run_all`).
- `99_verificar`: mismo resultado en R y Python (77/77, log con ambas secciones).
- Idempotente: `run_all` corrido dos veces → mismos `verificaciones.csv` /
  `procedencia.csv` / `comparacion_R_python.csv` (el bloque de auditoría del
  informe excluye a `12_informe`, así que no crece entre corridas).
- Pendiente de confirmar en esta sesión: `.\run_all.ps1 -FromSynthetic`
  (corriendo; los números cambian, la estructura y las verificaciones no).

### Entorno

- Sin dependencias nuevas. `12_informe`: sólo stdlib de Python / base de R
  (`base64` propio en R). `99_verificar`: stdlib / base. `renv.lock` y
  `requirements.txt` sin cambios.

### Cierre de T11 — PROYECTO TERMINADO

- Checklist de AGENTS §1 completo y verificado por `.\run_all.ps1`.
- **T11 cerrado. El proyecto está terminado.** No queda tarea pendiente en la
  tabla de la Sección 6.

### Notas para una eventual sesión futura

- Si se agrega/renombra una figura o un `.md` de reporte, actualizar las listas
  `SECCIONES` (12) y `FIG_ACTO1`/`FIG_ACTO2`/`REPORTES_MD` (99) en **los dos**
  lenguajes.
- Si `12_informe` toca el texto fijo (D1–D12, limitaciones, CSS), el cambio va
  idéntico en `.py` y `.R` o rompe la paridad de `informe.textonly.html`.
- El mini Markdown→HTML de `12_informe` cubre lo que hoy usan los `.md`; si un
  reporte futuro usa sintaxis nueva (blockquotes, código con fence, listas
  con `+`), hay que extender el conversor en ambos lenguajes a la vez.
- `docs/` y `logs/` están en `.gitignore`: para publicar el informe en GitHub
  Pages hay que forzar el add (`git add -f docs/informe.html`) y aceptar que
  lleva números de datos reales.

---

## Sesión 13 — 2026-09-11 — pedido post-cierre: T7 estratificado por sexo

> El proyecto ya estaba **terminado** (T11, sesión 12). Esta sesión aplica un
> cambio de alcance pedido explícitamente por el usuario después del cierre,
> vía `pedidos/cambios_acto2_correlaciones_por_sexo.md` ("pegar como mensaje al
> agente"). No reabre D1–D12. El otro archivo en `pedidos/`
> (`boxplots_acto1_base_R.R`, reescritura de 07 en R base) **queda sin tocar**:
> es una tarea aparte, marcada NO PROBADO por su autor, con 2 de 5 nombres de
> columna ya verificados como incorrectos contra
> `qpcr_modelos_clasificacion.csv` real (`p_interaccion`→`p_SEXOxTTO`,
> `metodo`→`rama_cascada`) — pendiente para una sesión futura si se pide.

### Qué se completó

- **`R/00_config.R` + `python/00_config.py`**: dos listas nuevas,
  `GENES_SPLOM_PLACENTA` (7 transportadores + `il6` + `gp130`, 9) y
  `GENES_SPLOM_BRAIN` (7 transportadores + `gp130`, 8) — `il6` es cuantificable
  en placenta pero no en cerebro (D7); conjunto distinto por tejido, pedido
  explícito. `stopifnot`/`assert` de tamaño y contenido.
- **`R/08_acto2_correlaciones.R` + `python/08_acto2_correlaciones.py`**
  (equivalentes), reescritos sobre la base de T7 (sesión 8):
  - **Estratos**: de 3 (`GLOBAL`, `CONTROL`, `LPS`) a 7, agregando las 4 celdas
    `SEXO×TTO` (`HEMBRA_CONTROL`, `HEMBRA_LPS`, `MACHO_CONTROL`, `MACHO_LPS`).
    `SEXO` se lee directo de `qpcr_cuantificacion_long.tsv` (ya normalizado por
    `02_ingesta_qc`, no se derivó de `GRUPO`). **Revierte** la decisión anterior
    de T7 que limitaba a 3 estratos "porque las celdas SEXO×TTO darían IC
    inútiles" — documentado en `analisis_descartados.md` como corrección
    explícita, con la limitación (`n≤9`) declarada en el propio texto.
  - **`il6R` excluido de todo el Acto 2** (tablas y las 6 figuras SPLOM):
    detección insuficiente en cerebro (3/7/3/3) que, al estratificar, deja casi
    todas las celdas bajo el piso de 5. Se conserva en Acto 1 (07, 05) — la
    exclusión es solo para correlaciones. `ITEMS` pasa de 10 a 9 (8 genes +
    `score_compuesto`).
  - **9 figuras nuevas** `acto2_corr_placenta_cerebro_<item>.png`: paneles
    `Females`/`Males`, Control y LPS superpuestos (color+forma), ajuste lineal +
    banda 95% solo si `n≥5` por celda, leyenda al pie con
    `rho`/`p`/`n` de las 4 celdas, ejes en `2^(-ΔΔCt)` (log2) — salvo
    `score_compuesto`: es un z-score (puede ser negativo), eje lineal,
    desviación documentada.
  - **4 SPLOM nuevos** por tejido×sexo (`_HEMBRA`/`_MACHO`), más los 2
    existentes actualizados: conjunto de variables pasa a ser
    tejido-específico (9 en placenta, 8 en cerebro, antes 7 fijos en ambos);
    etiqueta del triángulo superior `rho:` en vez de `Corr:` (ambigüedad con
    Pearson); `n` visible en cada panel; respeta el piso de 5.
  - `tabla_coexpresion` gana columna `SEXO` (`AMBOS`/`HEMBRA`/`MACHO`).
  - Verificaciones nuevas: `acto2_estratos` (7 estratos exactos),
    `acto2_estratos_particion` (n hembra+macho == n agregado, por item),
    `acto2_paneles_n_leyenda` (n≤9 por panel), `acto2_sin_pearson` (ninguna
    **tabla** —no la prosa metodológica, que sí dice "no Pearson" a
    propósito— contiene valores de Pearson), `acto2_il6R_fuera`,
    `acto2_splom_variables`. 12 verificaciones de este script, todas TRUE.

### Bugs propios encontrados y corregidos en esta sesión

- **Rutas absolutas en `procedencia.csv`**: las 16 figuras nuevas se registraron
  primero con la ruta absoluta de `RUTA_FIGURAS`/`cfg.RUTA_FIGURAS` en vez de
  `outputs/figures/<archivo>` (convención del proyecto). `98_comparacion`
  las marcaba "sin fila". Corregido con un helper `_rel_fig()` en ambos
  lenguajes.
- **Paridad de texto R/Python en `analisis_descartados.md`**: dos frases de la
  sección `08_acto2_correlaciones` quedaron redactadas con una diferencia
  menor entre `DESCARTES` (R) y `_DESCARTES` (Python) — `99_verificar` lo
  detectó como falla de "paridad render" (`informe.textonly.html` +
  `analisis_descartados.md` no byte-idénticos). Corregido igualando el texto.
- **Falso positivo en `acto2_sin_pearson`**: la primera versión buscaba la
  palabra "Pearson" en todo el texto generado, incluida la prosa que
  **explica que no se usa** ("Spearman rho (no Pearson)"). Corregido para que
  el chequeo mire solo las columnas/valores de las tablas (`COLS_CORR`,
  `COLS_COEXP`, `corr`, `coexp`), que es lo que pide el punto 3 del pedido.
- **Layout de la figura por gen en matplotlib**: la primera versión usaba
  `fig.text` a coordenadas fijas + `tight_layout()`, que no reserva espacio
  para texto puesto a mano → la etiqueta del eje X se superponía con la
  leyenda al pie. Corregido con `subplots_adjust` explícito (sin
  `tight_layout`) y las cuatro zonas (título, subtítulo, ejes, leyenda+nota)
  en posiciones fijas con margen suficiente.

### Verificado

- **R↔Python byte-idénticos**: `acto2_correlaciones.csv` (63 filas, antes 30),
  `acto2_coexpresion_transportadores.csv` (192 filas, antes 42),
  `acto2_correlaciones_reporte.md`, `analisis_descartados.md` (sección de este
  script), `procedencia.csv`, `verificaciones.csv` — sobre datos reales.
  Cruza-verificación R vs `cor.test(exact=FALSE)`: peor `|dif|` = 8.9e-16.
- **`.\run_all.ps1` completo (PRIME→R→Python→VERIFY): `TODAS LAS
  VERIFICACIONES PASARON` en R y en Python, 77/77 chequeos duros, 31/31 CSV
  byte-idénticos, paridad de render OK.** ~4.3 min.
- Las figuras nuevas se revisaron visualmente (una por gen, un SPLOM por sexo)
  en R y en Python: paneles, colores, KDE de la diagonal (Python ahora
  también separa Control/LPS en la diagonal del SPLOM, antes era un
  histograma único — mejora de equivalencia visual con R, no solo para las
  figuras nuevas), leyendas y piso de 5 se comportan como pide el texto.

### Resultados (datos reales; NO se versionan)

- GLOBAL sin cambios respecto de T7 (mismos 8 genes, antes 9 con il6R
  incluido): `slc38a1` rho=0.47 p=6.4e-3 sigue siendo el único con p<.05;
  `score_compuesto` rho=0.24 p=0.165.
  Por sexo (nuevo, descriptivo — **no se compara**, prohibición 4): varios
  genes muestran rho más alto en hembras que en machos en algunas celdas
  (p. ej. `slc38a1` ♀LPS rho=0.67 p=0.050 n=9 vs ♂LPS rho=0.26 p=0.531 n=8),
  con IC anchos por el n≤9 — el test formal de esa diferencia sigue siendo T8
  (no se re-corrió T8 por sexo en esta sesión; el pedido no lo pidió).

### Pendiente / siguiente paso concreto

- **`pedidos/boxplots_acto1_base_R.R`** sigue sin aplicar: reescritura de los
  boxplots del Acto 1 en R base con cascada D11 ampliada (bracket de efecto
  principal si la interacción no es significativa). Antes de tocar código:
  confirmar con el usuario si se aplica tal cual (con las 2 correcciones de
  nombre de columna ya detectadas) y si corresponde una versión Python
  equivalente (el pedido no trae una).
- No se actualizó `docs/informe.html` con las 9+4 figuras nuevas del Acto 2
  más allá de lo que `12_informe` ya arma automáticamente (que solo referencia
  las figuras que el propio `12_informe` lista por nombre). El pedido no pidió
  explícitamente incorporarlas al cuerpo del informe; si se quiere, es un
  cambio aparte en `12_informe.{R,py}` (dos listas `FIG_ACTO2` a extender).
- Falta decidir si esta sesión ameritaba un número de tarea nuevo en la tabla
  de la Sección 6 de AGENTS/CLAUDE.md (el proyecto ya estaba cerrado en T11);
  se dejó sin tocar esa tabla para no sugerir que el checklist original quedó
  incompleto.

---

## Sesión 14 — 2026-09-11 — pedido post-cierre: boxplots Acto 1 + D11 ampliada

> Segundo pedido pendiente de la sesión 13, ahora aplicado: `pedidos/
> boxplots_acto1_base_R.R`, tratado como **especificación de estilo y de
> lógica de anotación, no como código a copiar** (instrucción explícita del
> usuario). Cambia D11 (protegida por la prohibición 9) por **pedido explícito
> del usuario** — no es una reapertura por cuenta propia.

### Qué se completó

- **AGENTS.md / CLAUDE.md** (vueltos a copiar byte a byte): fila **D11**
  reescrita (cascada de 3 ramas) + nueva **sección 4.2** ("Cascada de
  anotación de brackets, D11 ampliada") con la tabla de las 3 situaciones, la
  aclaración de que (b)/(c) no son excluyentes, y la nota de por qué 07 pasa a
  R base para los boxplots de expresión.
- **`R/07_figuras_acto1.R`**: reescrito.
  - **D11 ampliada** (`d11_texto` + `d11_brackets_especificacion`, UNA función
    para las 3 ramas, compartida por los paneles de expresión y por pSTAT3):
    (a) interacción significativa → brackets por par (post hoc D6, sin
    cambios); (b) sin interacción y `p_TTO` <0.1 → bracket único 0–3,
    `Control vs LPS`; (c) sin interacción y `p_SEXO` <0.1 → bracket 0.5–2.5,
    `♀ vs ♂`. (b)/(c) se apilan si ambos aplican.
  - **2 nombres de columna corregidos** contra `qpcr_modelos_clasificacion.csv`
    real (detectados en la sesión 13, aplicados ahora): `p_interaccion`
    (no existe) → **`p_SEXOxTTO`** (solo para mostrarla en el subtítulo del
    panel; el gate sigue siendo `interaccion_significativa`); `metodo` (no
    existe en esa tabla) → **`rama_cascada`**.
  - **Boxplots de expresión pasan de ggplot2 a R base** (`boxplot()`), único
    cambio de motor gráfico pedido explícitamente: el estilo necesita bigote
    punteado más claro y tope del bigote sólido con trazo distinto del borde
    de la caja, algo que `geom_boxplot` no expone (`whisklty/whiskcol`,
    `staplelty/staplecol`, `border` sí lo permiten en `boxplot()`/`bxp()`).
    **El Acto 2 (08) sigue en ggplot2/GGally**; pSTAT3 (mismo script) también
    sigue en ggplot2 — el pedido de estilo aplicaba solo a las figuras de
    expresión.
  - **Paleta y layout nuevos** (estilo pedido): filas del panel agrupadas por
    vía metabólica (lipídos 3 / glucosa 2 / aminoácidos 2 / IL-6 hasta 3, no
    el orden crudo de `GENES`), color Control por sexo (celeste) y LPS por
    sexo × vía metabólica (una paleta por vía). `il6R` sale del panel de
    cerebro por prolijidad visual (se conserva en la tabla de modelos y en el
    panel de placenta). Subtítulo por panel: `rama_cascada · p SEXOxTTO · n`.
    Caja solo si el grupo tiene ≥3 detectados (sin cambios respecto de antes).
  - pSTAT3 (ggplot2) sin cambios de estilo, solo cambia la cascada que decide
    los brackets (en los datos reales cae igual en la rama (a), sin cambios
    visibles: interacción ya era significativa).
- **`python/07_figuras_acto1.py`**: espejado funcionalmente.
  - matplotlib **ya** permite estilar caja/bigote/tope por separado
    (`boxprops`/`whiskerprops`/`capprops`): este lenguaje **no** tuvo que
    cambiar de librería, solo replica el mismo estilo (paleta por vía,
    bigote punteado, tope sólido, filas agrupadas por vía metabólica).
  - Misma cascada D11 (`d11_texto` + `d11_brackets_especificacion`, texto
    idéntico a R carácter por carácter donde corresponde).

### Bugs propios encontrados y corregidos en esta sesión

- **Bug real preexistente en `leer_tabla`/`split_keep` (R), desde T6
  (2026-09-09), nunca disparado hasta ahora**: al retipear la función se
  reemplazó sin querer un carácter de control invisible (`\001`, centinela
  usado para evitar que `strsplit` descarte la última columna) por una cadena
  vacía `""`. Con `""` la función **siempre** pierde la última columna de
  cualquier TSV leído con `leer_tabla(...,"\t")`. Nunca se notó antes porque
  la única vez que la última columna importa es `PSTAT3` en
  `pstat3_long.tsv` — y ahí crasheaba (`data.frame(...)`: "argumentos implican
  un número diferente de filas: 1, 0"). Diagnosticado comparando byte a byte
  contra `git show HEAD:R/07_figuras_acto1.R`. Corregido restaurando el
  centinela como `"\001"` (escape octal, en vez de un byte crudo invisible en
  el fuente). **Se recomienda revisar si este patrón (`leer_tabla`/
  `split_keep`) se reutiliza en algún otro script** — por ahora solo lo usa
  07, y con el fix ya verificado extremo a extremo.
- **Layout de matplotlib**: los brackets apilados con espaciado multiplicativo
  fijo (`ymax * 1.30**k`) colisionaban con el título/subtítulo del panel en
  paneles de rango de datos angosto (el margen queda enorme) o con más de un
  nivel de bracket en paneles de rango ancho (el margen entre niveles queda
  minúsculo) — el mismo salto absoluto en escala log es una fracción muy
  distinta del alto visible según cuánto rango cubran los datos. Corregido
  con `_brackets_geometria()`: techo del eje y posición de cada nivel como
  **fracciones fijas del alto total del panel en log10** (los datos ocupan
  una fracción fija, el resto se reparte en franjas iguales por nivel + un
  margen), consistente sin importar el rango de cada gen; título/subtítulo
  quedan afuera de los ejes vía `annotate` con offset en puntos, y los
  brackets con `clip_on=True` para que nunca puedan dibujar sobre ese margen.

### Verificado

- **R↔Python**: mismos 21 brackets de expresión + 2 de pSTAT3, con el mismo
  texto exacto (`Control vs LPS *`, `♀ vs ♂ **`, etc.) en ambos lenguajes —
  comparado por stdout, no solo por CSV. `analisis_descartados.md` (sección
  `07_figuras_acto1`) byte-idéntica R↔Python.
  Figuras revisadas visualmente (ambos lenguajes, PLACENTA_E15 y BRAIN_E15):
  filas por vía metabólica, paleta, bigote/tope diferenciados, brackets
  nuevos en placenta (antes 0 brackets ahí, ahora 13) sin superposiciones.
- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 31/31 CSV byte-idénticos, paridad de render
  OK, ~4.3 min.**

### Resultados (datos reales; NO se versionan)

- **Placenta E15 pasa de 0 a 13 brackets** (0/10 genes marcados → 9/10, todos
  vía efecto principal: ninguna interacción SEXO×TTO es significativa en
  placenta, pero varios genes tienen `p_TTO` y/o `p_SEXO` <0.1). Ejemplos:
  `il6` `Control vs LPS ***` (`p_TTO`=6.0e-06); `glut1` con **ambos**
  brackets apilados (`Control vs LPS *` y `♀ vs ♂ **`); `il6R` sigue sin
  marca (ni interacción ni efectos principales cruzan 0.1).
- **Cerebro E15**: los 7 genes con interacción significativa (T5) mantienen
  sus brackets por par sin cambios; `glut3` (interacción no significativa,
  `p_TTO`=0.0178) gana un bracket nuevo `Control vs LPS *` que antes no
  existía.
- pSTAT3: sin cambios (interacción ya significativa, cae en la rama (a)).

### Pendiente / siguiente paso concreto

- Ninguna tarea pendiente conocida de los dos pedidos post-cierre (sesión 13
  y esta). Si aparece un pedido nuevo, seguir el mismo patrón: leer AGENTS.md
  + ESTADO.md, avisar en dos líneas antes de escribir código, adaptar nombres
  de columna contra los CSV reales (no asumir los de la especificación), y
  cerrar con `.\run_all.ps1` completo antes de dar por terminado.

---

## Sesión 15 — 2026-09-14 — pedido post-cierre: brackets/eje, grilla centrada, paleta pSTAT3

> Tres cambios visuales pedidos explícitamente sobre `07_figuras_acto1.{R,py}`,
> uno por commit: (1) brackets sin perpendiculares + eje Y fijado por los
> datos; (2) paneles de 2 genes centrados (grilla de 6 columnas); (3) paleta
> de pSTAT3 alineada a la de los boxplots de expresión, colores tomados de
> `00_config.{R,py}`. No reabre D11 (la cascada de 3 ramas no cambia, solo su
> geometría de dibujo).

### Cambio 1/3 — brackets sin perpendiculares + eje Y fijado por los datos

**Bug reportado**: en `acto1_expresion_PLACENTA_E15.png` el eje llegaba a 10⁴
con datos que no pasan de 10; en `il6R` (sin un solo bracket) llegaba a 10⁵.
Causa: el techo del eje (`ylim`) se calculaba como `ymax * 1.30^(n_niveles+1)`
con `n_niveles = max(length(especs), 1)` — **siempre** con un piso de 1 nivel,
incluso en paneles sin ningún bracket, y creciendo multiplicativamente con la
cantidad de niveles apilados.

- **`bracket_base`/`_bracket`**: pasan a dibujar una **línea horizontal
  simple** (sin las dos perpendiculares en los extremos) con el texto encima.
  Más correcto: un bracket de efecto principal (ramas (b)/(c) de D11) marca un
  efecto que abarca los 4 grupos o los dos sexos, no una comparación puntual
  entre dos extremos. La versión de pSTAT3 (ggplot2/matplotlib) **no cambia**
  — sigue con el rectángulo de extremos (fuera del alcance del pedido, que
  aplicaba a "los boxplots de expresión"); en Python se separó en
  `_bracket_pstat3` para no romper esa figura al simplificar `_bracket`.
- **`rango_eje_paneles`/`_rango_eje_paneles`** (función nueva): el rango del
  eje sale **solo de los datos** (`ymin*0.6` .. `ymax*1.15`, mismo padding de
  siempre). **Solo si hay brackets** (`n_niveles = length(especs)`, sin piso
  artificial) se reserva una **fracción fija** del alto total del panel en
  log10 (`frac_reservada = 0.25`, constante — no crece con la cantidad de
  niveles), repartida en franjas iguales por nivel. Sin brackets, el eje no
  reserva nada: el boxplot es el protagonista.
- Verificado visualmente (R y Python, PLACENTA_E15 y BRAIN_E15): paneles sin
  bracket (`il6R`) con eje ajustado a los datos; paneles con 2 brackets
  apilados (`fatp4`, `glut1`, `slc38a2`) sin colisión con el título/subtítulo.

### Bug propio encontrado y corregido

- **Paridad de texto R/Python** en la sección nueva de `analisis_descartados.md`:
  la primera redacción citaba literalmente la fórmula del bug (`ymax * 1.30^(n_niveles+1)`
  en R vs `ymax * 1.30**k` en Python) y los nombres de función (`rango_eje_paneles()`
  vs `_rango_eje_paneles()`), que difieren entre lenguajes por convención de
  nombres. `99_verificar` lo detectó como falla de paridad de render
  (`informe.textonly.html` + `analisis_descartados.md` no byte-idénticos).
  Corregido redactando en prosa neutra (sin literal de código ni nombre de
  función específico de un lenguaje).

### Verificado

- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 31/31 CSV byte-idénticos, paridad de render
  OK, ~4.6 min.**

### Cambio 2/3 — grilla de 6 columnas para centrar las filas de 2 genes

**Bug reportado**: las filas de 2 genes (`SLC38A1`/`SLC38A2`, `gp130`/`il6`)
quedaban corridas hacia un costado en vez de centradas. Causa: `layout()` (R)
y `subplots` (Python) usaban una grilla de 3 columnas real, con la fila de 2
desplazada 1 columna (pegada al borde derecho, no centrada).

- Grilla de **6 columnas**, cada panel ocupa 2: la fila de 3 usa las columnas
  1-2, 3-4 y 5-6; la de 2 usa 2-3 y 4-5 (centrada dentro del ancho de la fila
  de 3). R: `layout()` con una matriz de 6 columnas donde cada panel repite su
  índice en las 2 columnas que ocupa (fusiona la región). Python:
  `matplotlib.gridspec.GridSpec(nrow, 6)` + `fig.add_subplot(gs[i, c0:c0+2])`.
- Verificado visualmente (R y Python, PLACENTA_E15 y BRAIN_E15): las 2 filas
  de 2 genes por tejido quedan centradas, sin cambios en las filas de 3.

### Verificado

- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 31/31 CSV byte-idénticos, paridad de render
  OK, ~4.5 min.**

### Cambio 3/3 — paleta de pSTAT3 alineada a los boxplots de expresión

**Pedido**: que pSTAT3 siga la misma lógica de color que los boxplots de
expresión (celeste `COL_CTRL` para Control, macho más oscuro que hembra) y
que el LPS use el violeta de la vía IL-6 (`COL_LPS$IL6`, `#C5A3E0` hembra /
`#7B4EA8` macho) — coherente porque pSTAT3 es señalización de IL-6. Colores
tomados de `00_config`, no reescritos a mano en el script de figuras.

- **`COL_CTRL` y `COL_LPS` se mueven de `07_figuras_acto1.{R,py}` a
  `00_config.{R,py}`**: única fuente de estos colores. `07_figuras_acto1` los
  referencia (`cfg.COL_CTRL`/`cfg.COL_LPS` en Python; en R quedan en el mismo
  entorno vía `source()`).
- **`figura_pstat3`** deja de usar la paleta Okabe-Ito (`COL_TTO`, azul/naranja
  genérica, eliminada del script) y pasa a colorear por grupo SEXO×TTO (4
  colores, no 2) vía el helper nuevo `color_pstat3(sexo, tto)` — misma lógica
  que `color_for()` pero sin el argumento `gen` (siempre usa la vía IL6). En R
  esto requirió agregar la columna `sexo`/`grupo` al data.frame y pasar de
  `scale_colour_manual(values = COL_TTO)` (2 niveles) a una paleta de 4
  (`pal4`, una entrada por `SEXO TTO`).
- Verificado visualmente (R y Python): Control celeste (macho más oscuro),
  LPS violeta IL-6 (macho más oscuro que hembra), en las 4 cajas.

### Verificado

- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 31/31 CSV byte-idénticos, paridad de render
  OK, ~4.3 min.**

### Pendiente / siguiente paso concreto

- Ninguna tarea pendiente conocida de esta sesión (los 3 cambios pedidos
  quedaron aplicados y verificados, un commit por cambio). Si aparece un
  pedido nuevo: leer AGENTS.md + ESTADO.md, avisar en dos líneas antes de
  escribir código, y cerrar con `.\run_all.ps1` completo.

---

## Sesión 16 — 2026-09-14 — bugfix: `12_informe` no incrustaba 13 figuras nuevas

> Bug reportado por el usuario: las 9 figuras `acto2_corr_placenta_cerebro_*.png`
> (T7/sesión de correlaciones por sexo) y los 4 SPLOM por sexo existen en
> `outputs/figures/` y tienen fila en `procedencia.csv`, pero no aparecían en
> `docs/informe.html` — la sección "Acto 2.1-2.2" seguía mostrando solo
> `acto2_dispersion_placenta_cerebro.png` (la vista global vieja, sin separar
> por sexo). Causa: `SECCIONES` en `12_informe.{R,py}` traía la lista de PNG
> de cada sección **escrita a mano**, y nunca se actualizó cuando `08_acto2_correlaciones`
> agregó esas 13 figuras en una sesión anterior.

### Qué se completó

- **`SECCIONES`** deja de traer una lista de PNG a mano en su 3er elemento y
  pasa a un **criterio** (`script` de `procedencia.csv` + `patron` opcional de
  nombre): las figuras de cada sección se derivan solas filtrando
  `procedencia.csv` en `figuras_de_seccion()` (nueva). El `patron` solo hace
  falta para separar los dos usos de `07_figuras_acto1` (boxplots de expresión
  vs pSTAT3, mismo `script`, filas distintas); para el resto de las secciones
  alcanza con el `script`.
- Efecto directo: la sección "Acto 2.1-2.2" ahora muestra las **9 figuras
  `acto2_corr_placenta_cerebro_<item>.png`** (por gen + score compuesto), los
  **6 SPLOM** (2 globales + 4 por sexo) y la vista global vieja
  (`acto2_dispersion_placenta_cerebro.png`, que sigue siendo una figura
  legítima del mismo script — no se borra, solo deja de ser la única).
  **Figuras incrustadas: 13 → 26.**
- **Verificación nueva `informe_figuras_procedencia_embebidas`**: compara el
  set de figuras con fila en `procedencia.csv` contra las efectivamente
  incrustadas (`figuras_embebidas()`, la unión de `figuras_de_seccion()` sobre
  todas las `SECCIONES`) y falla si sobra alguna sin embeber. **Esta es la
  verificación que habría detectado el bug** — antes solo existía la
  verificación inversa (`informe_figuras_incrustadas`: que las figuras
  *listadas a mano* existan en disco, que no detecta figuras *ausentes de la
  lista*).
- Mismo patrón en R (`figuras_procedencia`/`figuras_de_seccion`/`figuras_embebidas`,
  usando `.col()`/`leer_csv_sin_este()` ya existentes) y Python (ídem con
  `_col()`/`leer_csv_sin_este()`).

### Verificado

- Recuento de `figcaption` en `informe.textonly.html`: las 26 figuras de
  `procedencia.csv` aparecen, en el orden esperado (SPLOM, luego correlación
  por gen, luego la vista global, dentro de la sección de Acto 2.1-2.2).
- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 31/31 CSV byte-idénticos, **91/91**
  verificaciones en TRUE (antes 90 — la nueva se suma), paridad de render OK,
  ~4.3 min.**

### Pendiente / siguiente paso concreto

- Ninguna. Si se agregan figuras nuevas a algún script existente, aparecerán
  solas en su sección (mismo `script` en `procedencia.csv`); si se agrega una
  sección nueva en `SECCIONES`, sigue siendo manual (título + `.md` + criterio),
  como corresponde.

---

## Sesión 17 — 2026-09-14 — pedido: `cambios_acto2_dispersion_por_sexo.md`

> Pedido completo en `pedidos/cambios_acto2_dispersion_por_sexo.md`: dispersión
> y Δρ separados por sexo + test de interacción SEXO×TTO sobre dispersión +
> simulación extendida a los estratos por sexo. **Puntos 1–4 confirmados**;
> punto 5 (test de pendientes) queda **pendiente, no se implementa** (el
> pedido exige preguntar antes, y así se hizo). Corrección del usuario sobre
> el punto 3: el test de interacción sobre dispersión **no aplica la cascada
> D5** — ANOVA tipo III directo, sin selección de rama (ver ese punto).

### Punto 1/4 — Δρ por sexo (sección 2.4)

**Problema**: la tabla de Δρ (Fisher z, prohibición 4) agrupaba los sexos
(n=16–18), mientras que desde T7 las figuras de correlación ya están
separadas por sexo. El caso que lo deja claro es `fatcd36`: ♀Control
rho=0.86 (n=8), ♂Control rho=−0.43 (n=8); agrupados dan 0.47 (verificado:
0.8571/−0.4286/0.4706 en los datos reales).

- **`pares()`/`_pares()`** ganan un filtro opcional de sexo (`sexo_filtro`);
  `cargar()` carga `SEXO` por feto (ya estaba en `qpcr_cuantificacion_long.tsv`).
- **`COLS_TEST`** gana la columna `ESTRATO` (`AMBOS_SEXOS`/`HEMBRA`/`MACHO`).
  `AMBOS_SEXOS` = filas previas, sin cambios (aditivo); `HEMBRA`/`MACHO` son
  filas nuevas, mismo Fisher z (SE Bonett-Wright primario + clásico), sin
  cambiar el método.
  - **BH (D12) dentro de cada estrato**, no a través de los 27 `p` juntos
    (pedido explícito: la potencia difiere demasiado entre estratos).
- **Potencia declarada en el cuerpo del reporte** (no nota al pie): con
  `n<=9` por celda de sexo, el Fisher z tiene poca potencia.
- `figura_test_delta_rho` sigue mostrando solo `AMBOS_SEXOS` (un forest por
  sexo saturaría el panel); los estratos se leen en la tabla completa.
- Verificación nueva `acto2_estrato_particion_test`: `n(HEMBRA)+n(MACHO) ==
  n(AMBOS_SEXOS)` por item y por grupo — verificado 10/10 items.

### Verificado

- R↔Python: mismos valores (`fatcd36` HEMBRA rho=0.8571429, MACHO
  rho=−0.4285714, AMBOS_SEXOS rho=0.4705882, coincide con el ejemplo del
  pedido).
- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 31/31 CSV byte-idénticos, **93/93**
  verificaciones en TRUE (antes 91), paridad de render OK, ~4.4 min.**

### Punto 2/4 — Dispersión por sexo (sección 2.3)

Misma estratificación que el punto 1, aplicada a la tabla de dispersión
(`acto2_dispersion.csv`): `COLS_DISP` gana `ESTRATO` (insertada antes de
`TEJIDO`); el Levene Brown-Forsythe Control vs LPS por lado se corre también
dentro de `HEMBRA` y `MACHO`, además de `AMBOS_SEXOS` (sin cambios,
aditivo). `figura_dispersion_sd` sigue mostrando solo `AMBOS_SEXOS` (mismo
criterio que la figura de Δρ). Verificación nueva
`acto2_estrato_particion_dispersion`: partición por item×tejido (20 celdas),
verificada 20/20.

### Verificado

- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 31/31 CSV byte-idénticos, **94/94**
  verificaciones en TRUE (antes 93), paridad de render OK, ~4.3 min.**

### Punto 3/4 — NUEVO: test de interacción SEXO×TTO sobre la dispersión

Extensión factorial del Levene Brown-Forsythe (2.3): en vez de un factor
(`TTO`), ANOVA III de 2 factores (`SEXO`, `TTO`) sobre `z = |x − mediana de
su celda SEXO×TTO|`, por gen×tejido (universo = T5, `via == "modelo"`,
18/20 — se lee `qpcr_modelos_clasificacion.csv` del propio idioma). Tabla
nueva `acto2_dispersion_interaccion.csv`.

**Corrección del usuario aplicada**: esta es la única parte del proyecto
donde **D5 no se aplica** — ANOVA tipo III directo con contrastes suma-cero,
sin cascada. Razón (declarada en el reporte y en `analisis_descartados.md`):
el propio Levene Brown-Forsythe ya es un ANOVA sobre desvíos absolutos; los
desvíos absolutos son positivos y sesgados por construcción (Shapiro fallaría
casi siempre → ART sin motivo real), y evaluar homocedasticidad sobre una
variable que ya es una medida de dispersión es circular. Núcleo ANOVA III
(`resolver`/`ajustar`/diseño suma-cero/`anova3_terminos`) portado de
`05_qpcr_modelos` (mismo diseño de 4 columnas `[1, s, t, s*t]`).

**Verificación 7.2 del pedido** (`acto2_interaccion_vs_levene`): el test
colapsado a un solo factor (`TTO`, ignorando `SEXO`) tiene que reproducir
exactamente el Levene Brown-Forsythe. Verificado con `fatcd36@PLACENTA_E15`
(fijo, determinista): **F 2.009666 vs 2.009666, p 1.653996e-01 vs
1.653996e-01** — coincide a 6 decimales, tolerancia 1e-8 superada.

BH (D12) por término, dentro de cada tejido. Resultado (datos reales): 8/18
gen×tejido con `p_SEXOxTTO < .05` (incluye `fatcd36`, `fatp4`, `gp130` en
cerebro — 3 de los 4 genes que el pedido señala como "interacción
significativa en T5 sin post hoc que sobreviva Holm").

### Verificado

- R↔Python: mismos 8/18 genes significativos, mismos F/p exactos.
- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, **32/32** CSV byte-idénticos (nueva tabla
  incluida), **98/98** verificaciones en TRUE (antes 94), paridad de render
  OK, ~4.3 min.**
- Bug propio encontrado y corregido en el camino: la primera redacción de
  `analisis_descartados.md` citaba `fila_diseno`/`_fila_diseno` y
  `05_qpcr_modelos.R`/`.py` (nombres específicos de cada lenguaje) → falla de
  paridad de render; corregido con prosa neutra.

### Punto 4/4 — Simulación cubre los estratos HEMBRA/MACHO

`10_acto2_simulacion` gana columna `ESTRATO` (`AMBOS_SEXOS`/`HEMBRA`/`MACHO`):
mismo diseño (normal bivariada en `-ddCt`, `r` verdadera común como rango
GLOBAL/CONTROL), con las SD marginales = SD **observadas dentro de ese
estrato de sexo**. `figura_simulacion` sigue mostrando solo `AMBOS_SEXOS`
(mismo criterio que las demás figuras de 09/10). Verificación de partición
`acto2_sim_estrato_particion`: 20/20 celdas item×escenario.

**Resultado (limitación declarada en el reporte, no oculta)**: con `n<=9`
por celda de sexo, el veredicto `DENTRO` es casi automático — `MACHO` dio
18/18 DENTRO; `HEMBRA` 13/18 DENTRO, 5/18 FUERA. Con este `n` no se puede
distinguir cambio de coordinación de cambio de dispersión dentro de cada
sexo por separado; la lectura confiable sigue siendo `AMBOS_SEXOS`.

**Bug propio encontrado y corregido**: la fila `sin_test` (celdas bajo el
piso `n_par>=5`) tenía un campo de menos que columnas la tabla — bug
preexistente en R (nunca disparado porque `AMBOS_SEXOS` siempre superaba el
piso con estos datos), expuesto ahora porque `HEMBRA`/`MACHO` sí caen bajo
el piso para algunos items. Corregido completando los 20 campos de
`COLS_SIM` explícitamente (`rep("", 11L)` en R, `[""] * 11` en Python).

### Verificado

- R↔Python: mismos resultados exactos (`HEMBRA` 18 celdas/13 DENTRO/5 FUERA,
  `MACHO` 18/18/0, partición 20/20 OK).
- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 32/32 CSV byte-idénticos, **100/100**
  verificaciones en TRUE (antes 98), paridad de render OK, ~4.6 min.**

### Pendiente / siguiente paso concreto

- **Puntos 1–4 del pedido completos y verificados** (`pedidos/cambios_acto2_dispersion_por_sexo.md`),
  cada uno con su commit. **Punto 5 (test de pendientes)**: quedó pendiente en
  esta sesión (exigía preguntar antes de implementar) y se **descartó en la
  sesión 19** tras consultar al usuario — ver esa sesión y
  `analisis_descartados.md` (`09_acto2_dispersion`). No es tarea pendiente.

---

## Sesión 18 — 2026-09-21 — pedido: `cambios_informe_conclusiones.md`

> Pedido completo en `pedidos/cambios_informe_conclusiones.md`: agregar
> "Conclusión de la sección" al final de cada sección del informe, una
> síntesis del eje madre→placenta→cerebro, una conclusión revisada (Acto 1
> frente a Acto 2) y reescribir el Resumen. **Alcance: solo redacción** — no
> se corrió ningún test nuevo ni se tocó ningún script de análisis (02–11).
> Único archivo tocado: `12_informe.{R,py}` (más el propio pedido).
> **Puntos 1–4 aplicados y verificados; punto 5 (opcional, correcciones de
> figuras) queda sin aplicar**, tal como pide el pedido ("aplicalas solo si
> confirmo").

### Qué se completó

- **Verificación previa de cada número citado en el borrador del pedido**
  contra los CSV reales en `outputs/tables/R/` (fork de solo lectura, sin
  tocar el repo): **todos los números coincidieron** con los CSV — Fisher MS
  p=4.995e-03, Peto-Peto LA ♀p=0.302/♂p=0.333, 5 genes con `p_TTO<.05` en
  placenta, 7 genes con interacción SEXO×TTO en cerebro (glut1/slc38a2/fatp1
  con efecto en media; fatcd36/fatp4/gp130/slc38a1 con Holm mínimo=0.074),
  pSTAT3 ♀1.976→4.805 (Holm p=1.271e-7) / ♂2.695→2.752 (p=0.883), Δρ 0/28
  significativos (mínimo fatcd36 ♀ p=0.075), interacción sobre dispersión 5/18
  en cerebro sobreviven BH (gp130/fatcd36/fatp4/fatp1/slc38a2; slc38a1
  tendencia) y 0/10 en placenta, eigengene 76 %/82 % (placenta/cerebro,
  81.55 % redondea a 82), simulación fatcd36/score_compuesto en el límite del
  IC (estrato HEMBRA, escenario GLOBAL: 0.757 vs q975=0.755; 0.75 vs
  q975=0.733), exclusión de feto fatcd36 0.107→0.211 / score 0.079→0.175
  (0/11 cambian veredicto). **Ningún número tuvo que corregirse** — la única
  precisión: "CD36" del borrador es `fatcd36` en los datos (se usó el nombre
  de dato, no el informal).
- **`numeros_conclusiones()`** (nueva, una por lenguaje): mismo principio que
  `resumen_numeros()` — todo número y toda lista de genes citados en la nueva
  prosa se leen/derivan de los CSV en tiempo de generación del informe, nunca
  se escriben a mano (si los datos cambian en una corrida futura, el texto
  cambia solo). Incluye clasificación programática de los 7 genes de cerebro
  en "efecto en media" (algún contraste post hoc con `p_holm<.05`) vs "efecto
  en dispersión" (ninguno), y detección de los ítems en el límite del IC de
  la simulación (`ESTRATO=HEMBRA, ESCENARIO=GLOBAL, veredicto=FUERA`) — no se
  hardcodeó qué genes/ítems caen en cada grupo.
- **`conclusion_seccion(k, n)`**: bloque "Conclusión de la sección" (`<h4>`)
  insertado al final de cada sección de `SECCIONES` excepto Acto 1.2
  (cuantificación, de método — sin conclusión, como pide el punto 1). Prosa
  fija idéntica R/Python, con los números/listas interpolados.
- **Dos secciones nuevas**: "4. Síntesis del eje madre → placenta → cerebro"
  (entre Acto 1 y Acto 2) y "6. Conclusión revisada (Acto 1 frente a Acto 2)"
  (tabla de 3 filas + párrafo, entre Acto 2 y Reproducibilidad). Esto
  renumeró las secciones existentes: Reproducibilidad 5→7, Descartados 6→8,
  Limitaciones 7→9, Auditoría 8→10 (`items_toc`, títulos `<h2>` y el chequeo
  `n_sec == 10L`/`== 10` en `main()`, antes 8).
- **Resumen**: el bullet "Coordinación placenta↔cerebro" ahora agrega, tras
  la conclusión negativa de Δρ/simulación, una frase sobre el resultado
  positivo del Acto 2 (interacción SEXO×TTO sobre la dispersión en cerebro),
  para que el Resumen no quede desactualizado respecto de la sesión 17.
- **Redondeo determinista `.round_fmt()`/`_round_fmt()`** (floor manual, no
  `round()`): las medias/porcentajes de la nueva prosa se redondean igual
  bit a bit en R y Python sobre el mismo double, necesario porque este texto
  se byte-compara entre lenguajes (`informe.textonly.html`). Los p-valores se
  muestran tal cual vienen del CSV (mismo criterio que el Resumen preexistente),
  sin redondear.
- **Bug propio encontrado y corregido en el camino**: la tabla de
  "Conclusión revisada" se armó primero con `c(...)` en R (tres fragmentos
  de la fila de encabezado como elementos separados) mientras que Python los
  concatenaba en un único string — daba paridad de render rota (una línea de
  más en R al hacer `paste(..., collapse="\n")`). Corregido usando `paste0()`
  para que ambos lenguajes produzcan una sola línea. También se cambiaron
  `&female;`/`&male;` (no son entidades HTML5 válidas) por los caracteres
  Unicode literales ♀/♂ (`♀`/`♂`, idénticos en ambos lenguajes).

### Verificado

- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 32/32 CSV byte-idénticos, paridad de
  render OK, 100/100 verificaciones en TRUE (mismo total que sesión 17: no
  se agregó ninguna verificación nueva, solo se regeneró el informe), PDF
  ok, ~5.5 min.**
- `docs/informe.html`: 10 secciones (antes 8), 26 figuras incrustadas (sin
  cambios), `informe.textonly.html` R↔Python byte-idéntico (`diff` manual
  sobre ambos snapshots, exit 0).

### Pendiente / siguiente paso concreto

- **Punto 5 del pedido principal (test de pendientes, Acto 2.5.1)**: sigue
  pendiente al cerrar esta sesión — **descartado en la sesión 19** (ver esa
  sesión). No es parte de este pedido.
- **Sección 5 de `cambios_informe_conclusiones.md` (OPCIONAL, correcciones de
  consistencia)**: puntos (a) y (b) **aplicados en la sesión 19** (confirmados
  por el usuario). Sigue pendiente **(c)**: el triángulo superior de los SPLOM
  por sexo muestra un solo ρ (Control+LPS juntos) en vez de un ρ por
  tratamiento — no se aplicó, sin confirmación del usuario.

---

## Sesión 19 — 2026-09-21 — sección 5 (a)(b) de `cambios_informe_conclusiones.md` + descarte del test de pendientes

> Pedido: aplicar solo los puntos (a) y (b) de la sección 5 (opcional) del
> pedido de la sesión 18. (a) conservar el panel agrupado y agregar HEMBRA/
> MACHO al lado en las 3 figuras del Acto 2.3-2.5, de modo que la figura
> muestre el contraste que describe el texto; ajustar la conclusión si hace
> falta. (b) `il6R` seguía apareciendo en tablas/figuras del Acto 2 pese a
> estar excluido — corregir. El punto (c) (SPLOM) queda sin tocar. Además:
> descartar formalmente el "test de pendientes" (punto 5 de
> `cambios_acto2_dispersion_por_sexo.md`, sesión 17) en `analisis_descartados.md`,
> con las razones que dio el usuario, y sacarlo de "pendiente" en `ESTADO.md`.

### Qué se completó

- **Bug (b) encontrado y corregido**: `08_acto2_correlaciones` excluye `il6R`
  de todo el Acto 2 (`GEN_EXCLUIDO_CORR`, detección insuficiente en cerebro),
  pero `09_acto2_dispersion`, `10_acto2_simulacion` y `11_sensibilidad`
  definían `GENES_CORR` cada uno por su cuenta sin esa exclusión — `il6R`
  seguía en sus tablas y figuras (9 genes en vez de 8). Mismo
  `GEN_EXCLUIDO_CORR <- "il6R"` agregado a los tres scripts (R y Python).
  **No toca** la Sección 3 (test de interacción SEXO×TTO sobre dispersión,
  sesión 17): ese universo sale de T5 (`via == "modelo"`), no de
  `GENES_CORR` — es una pregunta distinta (por tejido, no del par
  placenta-cerebro) y `il6R@PLACENTA_E15` se sigue modelando ahí.
  **Efecto colateral encontrado y corregido**: esto bajaba `ITEMS_LOO` de 11
  a 10 en `11_sensibilidad`, y una verificación (`sens_excl_extremo_items`)
  tenía el `11` hardcodeado (`length(ITEMS_LOO) == 11L`) — habría quedado en
  FALSE sin este segundo fix.
- **(a) Las 3 figuras ganan grilla de 3 estratos** (`AMBOS_SEXOS` conservado +
  `HEMBRA` + `MACHO`), en vez de mostrar solo el agrupado:
  - `acto2_dispersion_sd.png` (09): `facet_grid(estrato ~ item)` en R /
    grilla `subplots(3, len(ITEMS))` en Python — filas = estrato, columnas =
    item. Se ve directamente el patrón cruzado (LPS baja SD en HEMBRA, la
    sube en MACHO).
  - `acto2_test_delta_rho.png` (09): `facet_grid(~ estrato)` en R / 3 ejes
    lado a lado en Python — mismo eje y (items) y mismo rango x en los 3
    paneles. Se extendió el `xlim` para que la anotación de texto quede
    dentro de cada panel (con `facet_grid` ya no sirve el truco de
    `clip="off"` + margen grande de un solo panel).
  - `acto2_simulacion_delta_rho.png` (10): mismo patrón que 09 (a),
    `facet_grid(estrato ~ item)`.
  - Verificado visualmente (R): las 3 figuras muestran el contraste
    correctamente (ver figuras en `outputs/figures/`).
- **`12_informe`**: la frase de la conclusión de 2.3-2.4 que decía "la figura
  de SD agrupa los sexos y por eso no muestra este efecto" ya no era cierta
  — reescrita para decir que la figura ahora muestra los 3 estratos y que el
  patrón se ve comparando las filas HEMBRA/MACHO (texto idéntico R/Python,
  verificado por la paridad de render).
- **Descarte del test de pendientes** (`analisis_descartados.md`, sección
  `09_acto2_dispersion`, nueva): tres razones dadas por el usuario — (1)
  contesta una pregunta distinta (co-expresión entre pares de genes dentro
  de un tejido, no coordinación placenta↔cerebro); (2) sumaría más de 100
  tests con `n<=9` por celda; (3) la pregunta de coordinación ya está
  respondida por Δρ (2.4) + simulación (2.5). Ya hay un control de
  co-expresión entre transportadores en 2.6 (eigengene) que cubre lo que el
  test de pendientes buscaría, sin necesidad de un test por par.
- Comentarios y textos de "9 genes"/conteos de items actualizados a "8
  genes"/nuevos totales en `09_acto2_dispersion`, `10_acto2_simulacion` y
  `11_sensibilidad` (R y Python), incluida la sección "Alcance y piso" de
  `analisis_descartados.md` (texto byte-idéntico R/Python, con un párrafo
  nuevo "BUGFIX" explicando la corrección de `il6R`).

### Verificado

- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 32/32 CSV byte-idénticos, paridad de
  render OK, 100/100 verificaciones en TRUE, ~4.9 min.**
- `il6R` ausente de `acto2_dispersion.csv`, `acto2_test_correlaciones.csv`,
  `acto2_simulacion.csv` y `acto2_sensibilidad_excl_extremo.csv` (grep, 0
  coincidencias en las 4 tablas).
- Números del informe se recalculan solos con el nuevo universo de 8 genes
  (nada hardcodeado): Δρ pasa de "0 de 28" a "0 de 27" (se pierde 1 test
  válido de `il6R`), `sens_excl_extremo_items` pasa a "10 items; 10 con
  Fisher z" (antes 11).
- Figuras inspeccionadas visualmente (ver arriba): las 3 muestran el
  contraste HEMBRA/MACHO correctamente.

### Pendiente / siguiente paso concreto

- **Punto (c) de la sección 5** (SPLOM por sexo, un solo ρ en el triángulo
  superior en vez de uno por tratamiento): sigue sin aplicarse, sin
  confirmación del usuario.
- Ningún otro pendiente conocido de esta sesión.

---

## Sesión 20 — 2026-09-21 — `pedidos/cambios_revision_codex.md`: puntos 1, 3 y 2

> El repositorio se revisó con un agente externo (Codex) sin contexto previo,
> en dos pruebas sobre clones limpios: `revisiones/AUDITORIA.md` (auditoría
> escéptica) y `revisiones/PRUEBA_REPRODUCCION.md` (prueba de reproducción).
> Pedido de correcciones en `pedidos/cambios_revision_codex.md`, con orden
> explícito: puntos 1+3 primero (cambian el informe), después 2
> (verificaciones), después 4 (reproducibilidad), después 5+6. **No se agregó
> ningún análisis estadístico nuevo** — D1–D13, los datos y los resultados no
> cambiaron. Todo lo de esta sesión es documentación, calidad de las
> verificaciones y derivación del texto desde las tablas.

### Punto 1 — D13 y moderar lenguaje (commit `a42e792`)

- **D13** (sin término de camada) agregada a `AGENTS.md`/`CLAUDE.md` como
  decisión ya tomada pero nunca escrita; no reabre D1–D12 (la regla "no
  reabrir" ahora cubre D1–D12, D13 es la nueva "fresca").
- Entrada en `analisis_descartados.md` (vía `05_qpcr_modelos`) con marcador
  `<<< COMPLETAR: evidencia de analisis previos >>>` para que el usuario
  complete el dato concreto — **sigue pendiente de completar**.
- Limitación nueva en el informe: independencia asumida entre fetos.
- Moderadas 5 frases (2 señaladas explícitamente por la auditoría + 3 más
  con el mismo patrón: "no son efectos independientes", "no dependen de")
  que atribuían una causa única donde los tests solo muestran compatibilidad.

### Punto 3 — texto derivado de tablas (commit `a42e792`)

- "Todas robustas a BH" (C1) ahora cuenta `p_SEXOxTTO_BH` real: "7 de 7
  sobreviven a BH" en vez de una afirmación sin cálculo.
- "Ningún gen mostró interacción en placenta" (C3) ahora se deriva del
  conteo real (`pla_int_n`); si placenta tuviera alguna interacción
  significativa en una corrida futura, la frase cambia sola.
- **3.4 — el informe sabe con qué datos se generó**: con fuente sintética,
  las 9 secciones interpretativas (7 conclusiones de sección + síntesis +
  conclusión revisada) nunca llaman a las funciones que arman esa prosa —
  se reemplazan por un aviso fijo. Verificación nueva
  `informe_sintetico_sin_interpretacion` (recalculo) cuenta el aviso en el
  HTML final (9/9 con `MIA_LPS_FORZAR_SINTETICO=1`, 0/0 con datos reales,
  probado en R y Python).
- Bug propio encontrado y corregido en el camino: el chequeo de esa
  verificación usaba `identical(gregexpr(...), -1L)`, que nunca es `TRUE`
  porque `gregexpr` sin match devuelve un objeto con atributos — exactamente
  el tipo de "verificación que no puede fallar" que señala la auditoría.

### Punto 2 — dejar de contar lo que no verifica (commit siguiente)

- **Columna `tipo`** agregada a `verificaciones.csv` (posición 2, después de
  `id`) en los 13 scripts que escriben ahí (02–11, 98, 12; R y Python),
  clasificando las **108** filas actuales: **50 `recalculo`** (recalculan o
  contrastan un resultado de forma independiente), **9 `existencia`**
  (confirman que un archivo existe y no está vacío), **49 `declaracion`**
  (registran una decisión de método o una constante, no verifican su
  aplicación caso por caso). Mecánica: header + índice de columna `script`
  en el `clave_orden` de cada `registrar_verificaciones` (script→posición 7
  en vez de 6 en R, `f[6]` en vez de `f[5]` en Python).
- **12_informe** (§10 Auditoría) y **98_comparacion** (`comparacion_reporte.md`
  y consola) ya no dicen "X/Y en TRUE" sin desglose: reportan
  recalculo/existencia/declaración por separado, más un aviso `NO_EJECUTADA`
  cuando corresponda (ver abajo). Ninguna otra mención de "100/100
  verificaciones" en el repo (grep confirmado: no está en README).
- **Estado `NO_EJECUTADA`** agregado a la columna `ok` (además de
  `TRUE`/`FALSE`): una verificación que no corrió en esta corrida (PDF sin
  motor, o comparación R↔Python en `-Only R`, punto 4) no cuenta como pase
  ni como fallo. `99_verificar` (chequeo duro #4) y `98_comparacion`
  (`audit_verificaciones_ok`) tratan `NO_EJECUTADA` como advertencia, no
  como fallo — "TODAS LAS VERIFICACIONES PASARON" sigue siendo válido con
  filas `NO_EJECUTADA`, pero quedan listadas aparte.
- **Los 4 checks rotos, arreglados** (no solo reclasificados):
  - `informe_pdf`: antes aprobaba también con `fallo:*`. Ahora `ok` → TRUE,
    `sin_motor` → NO_EJECUTADA, `fallo:*` → FALSE.
  - `pstat3_rama_cascada` (06_pstat3): antes exigía `rama == "anova3"` fijo.
    Ahora recalcula la rama esperada desde los p de Shapiro/Levene
    registrados y la compara contra la elegida (cualquiera sea).
  - `modelos_rama_conteo` (05_qpcr_modelos): antes solo sumaba que las 3
    ramas dieran 18. Ahora valida, **gen por gen**, que la rama de cada uno
    corresponde a sus propios diagnósticos.
  - `acto2_test_items` (09_acto2_dispersion): antes solo reportaba
    `n_test`/`length(ITEMS)`. Ahora exige `n_test == length(ITEMS)`.
- **Bugs propios encontrados en el camino** (independientes del pedido):
  - `acto2_estratos` (08_acto2_correlaciones) estaba mal clasificado como
    `recalculo` en mi primer pase: `ESTRATOS` ahí es una constante
    hardcodeada (igual que en 09/10), no algo derivado de datos — corregido
    a `declaracion` antes de commitear.
  - Al borrar `verificaciones.csv` para poder migrar el schema (el guard de
    encabezado incompatible lo exige), la primera corrida completa mostró
    una paridad R/Python rota (94 vs 101 filas en el snapshot de
    `comparacion_reporte.md`) — es un artefacto transitorio de arrancar de
    un archivo vacío (98_comparacion en la pasada R ve el archivo antes de
    que 12_informe haya corrido nunca; en la pasada Python ya lo ve con las
    filas de 12_informe de la pasada R anterior). Se autocorrige solo en la
    segunda corrida completa (el archivo persiste entre corridas en uso
    normal). No se tocó código por esto.

### Verificado

- **`.\run_all.ps1` completo: `TODAS LAS VERIFICACIONES PASARON` en R y
  Python — 77/77 chequeos duros, 32/32 CSV byte-idénticos, paridad de
  render OK, ~4.8 min (segunda corrida; la primera mostró el artefacto
  transitorio de arriba).**
- `verificaciones.csv`: 108 filas, 50/50 recalculo + 9/9 existencia + 49/49
  declaracion en TRUE (verificado con `csv.DictReader` de Python, no con
  parseo naive de comas).
- `informe_pdf` en TRUE real (antes hubiera sido TRUE incluso con
  `fallo:*`); confirmado leyendo la fila de `verificaciones.csv`.

### Pendiente / siguiente paso concreto

- **`<<< COMPLETAR: evidencia de analisis previos >>>`** en
  `analisis_descartados.md` (sección `05_qpcr_modelos`, D13): el usuario
  tiene que dar el dato concreto de qué análisis previo mostró que `MADRE`
  no modificaba los resultados.
- **Punto 4 (reproducibilidad)** sin empezar: detectar Python en PATH y
  documentarlo en el README; `run_all.ps1 -Only R` no debe exigir Python;
  `99_verificar` con una sola implementación (usar el `NO_EJECUTADA` ya
  agregado en este punto para la comparación R↔Python); ruta de R buscada
  en PATH/instalación estándar, no hardcodeada; documentar
  `-ExecutionPolicy Bypass` y la caché de red de `renv` en el README.
- **Puntos 5 y 6** sin empezar: `revisiones/conciliacion_auditoria.md` (las
  5 afirmaciones de la auditoría contra la fila exacta del CSV) y
  `revisiones/RESPUESTA.md` (una fila por hallazgo de los dos informes).
- **Corrida final pedida**: `run_all.ps1` completo + `run_all.ps1 -Only R
  -FromSynthetic` (para probar que la corrida de una sola implementación
  con datos sintéticos funciona y no muestra conclusiones biológicas) —
  esto último depende de que el punto 4 esté resuelto (`-Only R` hoy sigue
  exigiendo Python).
- Punto (c) de la sección 5 de la sesión 19 (SPLOM) sigue sin aplicar.
