# Reanálisis reproducible — MIA-LPS en placenta E15 y cerebro fetal E15

Reanálisis bioestadístico del modelo murino de activación inmune materna (LPS 100 µg/kg
i.p., día 15 de gestación, colecta a las 6 h). Caracteriza la respuesta de la placenta y
del cerebro fetal a la activación inmune materna y evalúa si depende del sexo del feto.

El repositorio corre **de punta a punta sobre datos sintéticos**. Los datos crudos son
inéditos y **no** están versionados (ver "Datos" abajo).

> ⚠️ **Las conclusiones biológicas solo son válidas con los datos reales.** Sobre datos
> sintéticos, las salidas demuestran que el pipeline es completo y reproducible, nada más.
> Los efectos simulados son arbitrarios.

---

## 1. Estructura

Ver `AGENTS.md` §3. En resumen: implementación **gemela** en `R/` y `python/` (mismos
números y nombres de script), salidas en `outputs/`, informe en `docs/`.

**Antes de tocar nada, leé `AGENTS.md` y `ESTADO.md` completos.**

---

## 2. Instalación del entorno

Hacen falta **dos** intérpretes: R y Python. `run_all.ps1` exige solo el que el modo
elegido va a usar — `-Only R` no necesita Python, y `-Only python` no necesita R.

### 2.0 Cómo encuentra el pipeline a los intérpretes

`run_all.ps1` **no tiene rutas fijas**. Las rutas concretas que aparecen más abajo son
ejemplos de esta máquina, no requisitos. El script busca, en orden, y valida cada
candidato ejecutándolo:

| | `Rscript` | `python` |
|---|---|---|
| 1 | `Rscript.exe` en el `PATH` | `.venv\Scripts\python.exe` en la raíz del repo |
| 2 | registro `HKLM\SOFTWARE\R-core\R` → `InstallPath\bin\Rscript.exe` (lo escribe el instalador oficial) | `python.exe` / `python3.exe` en el `PATH` |
| 3 | `%ProgramFiles%\R\R-*\bin\Rscript.exe` y `%LOCALAPPDATA%\Programs\R\R-*\bin\Rscript.exe`, versión más nueva primero | el lanzador `py -3`, preguntándole `sys.executable` |
| 4 | — | `%LOCALAPPDATA%\Programs\Python\Python3*\python.exe`, `%ProgramFiles%\Python3*\python.exe`, `C:\Python3*\python.exe` |

Al arrancar, el script imprime las rutas que resolvió (`Rscript:` / `Python:`), o
`(no se usa en este modo)`. Si no encuentra uno que sí necesita, aborta con un
mensaje que remite a esta sección. Los candidatos se validan corriéndolos, así que
el alias de la Microsoft Store (`python3.exe` en `WindowsApps`, que solo abre la
tienda cuando Python no está instalado) no se confunde con un intérprete.

### 2.1 R (≥ 4.4; probado con 4.6.1)

Instalar R desde CRAN (<https://cran.r-project.org/bin/windows/base/>) con el
instalador oficial para Windows; deja el registro `R-core` que usa el punto 2 de la
tabla, así que **no hace falta agregar nada al `PATH`**. En esta máquina quedó en
`C:\Program Files\R\R-4.6.1` (ejemplo: cualquier versión ≥ 4.4 sirve).

```powershell
# Desde la raíz del repo. Reemplazar la ruta por la de tu instalación,
# o usar `Rscript` a secas si está en el PATH.
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "install.packages('renv'); renv::restore()"
```

**`renv::restore()` necesita red.** La primera vez baja ~80 paquetes de CRAN: en un
entorno aislado sin salida a internet, la restauración falla (le pasó a la prueba de
reproducción externa). Además `renv` guarda su caché **fuera del repositorio**, en
`%LOCALAPPDATA%\R\cache\R\renv`; el entorno donde se corra tiene que poder escribir
ahí. Una vez poblada la caché, restauraciones posteriores funcionan sin red.

Si todavía no hay `renv.lock`, instalar a mano:

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "install.packages(c('car','ARTool','emmeans','NADA','survival','sandwich','lmtest','nlme','rmarkdown','ggplot2','GGally','readxl','writexl','jsonlite'))"
```

**No hace falta pandoc ni `rmarkdown`.** `12_informe` arma `docs/informe.html` a
mano (misma lógica de strings en R y Python → salida comparable byte a byte). Para
`docs/informe.pdf` alcanza con tener Edge o Chrome instalado (impresión headless);
si no hay ninguno, el pipeline deja el HTML y sigue.

### 2.2 Python (≥ 3.11; probado con 3.13)

Instalar CPython desde <https://www.python.org/downloads/windows/> con el instalador
oficial (**no** el alias de la Microsoft Store), marcando *Add python.exe to PATH*.
En esta máquina es **Python 3.13 per-user**, en
`%LOCALAPPDATA%\Programs\Python\Python313\python.exe`, y está en el `PATH` — que es
como lo encuentra el pipeline (punto 2 de la tabla de 2.0). El instalador deja además
el lanzador `py`, que sirve de respaldo (punto 3).

Comprobación rápida:

```powershell
python --version    # Python 3.13.x
```

Las dependencias, en un entorno virtual dentro del repo (el pipeline lo prefiere
sobre el `PATH`, así que no hace falta activarlo para correr `run_all.ps1`):

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2.3 Política de ejecución de PowerShell

Si `.\run_all.ps1` no arranca y PowerShell contesta *"no se puede cargar porque la
ejecución de scripts está deshabilitada en este sistema"*, correrlo con la política
relajada **solo para ese proceso**, sin tocar la política de la máquina ni del
usuario:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_all.ps1
```

El `-ExecutionPolicy` de la línea de comandos vale únicamente para ese proceso; al
cerrarlo no queda nada cambiado.

---

## 3. Datos

Los tres archivos crudos **no están en el repo** (`data/raw/` en `.gitignore`). Para
correr con datos reales, copiar a `data/raw/` con estos nombres exactos:

| Archivo en `data/raw/` | Hoja |
|---|---|
| `Raw data CTs.xlsx` | `Sheet1` |
| `ELISA IL6 2026 Dosis 100.xlsx` | `Sueros y LA` |
| `pstat3 placenta.xlsx` | `Sheet1` |

**Sin los crudos, el pipeline usa `data/synthetic/`** (generado por `01_generar_sinteticos`),
que replica estructura, nombres de columna, tipos y patologías de los reales.

---

## 4. Ejecución completa desde cero

```powershell
.\run_all.ps1
```

Corre `00_config` → `01_generar_sinteticos` → `02_ingesta_qc` → `03_elisa` →
`04_qpcr_cuantificacion` → `05_qpcr_modelos` → `06_pstat3` → `07_figuras_acto1` →
`08_acto2_correlaciones` → `09_acto2_dispersion` → `10_acto2_simulacion` →
`11_sensibilidad` → `98_comparacion` → `12_informe` → `99_verificar`.

En modo `both` (el default) lo hace en **pasadas separadas por lenguaje**: primero
Python 00–11 (prime), después R 00–11/98/12, después Python 00–11/98/12, y al final
`99_verificar` en R y en Python. Es la única forma de que `12_informe` deje en disco
los `.md` y el `informe.html` de cada lenguaje para que `99_verificar` los
byte-compare. Toma ~3 min.

`run_all.ps1` acepta flags (ver cabecera del script): `-Only R|python` (una sola
pasada, sin `98_comparacion` pero **con** `99_verificar`), `-FromSynthetic` (fuerza
`data/synthetic/` aunque haya crudos), `-SkipReport`, `-SkipVerify`.

`99_verificar` imprime `TODAS LAS VERIFICACIONES PASARON` cuando (a) existe y no
está vacío cada ítem del checklist de `AGENTS.md` §1, (b) los CSV de
`outputs/tables/{R,python}/` concuerdan celda a celda, (c) los `.md` y el
`informe.html` sin figuras son byte-idénticos entre la corrida R y la Python,
(d) ninguna fila de `verificaciones.csv` quedó distinta de `TRUE`, **y** (e) nada
quedó sin ejecutar.

### Verificar con una sola implementación

Quien tenga solo R (o solo Python) puede verificar su mitad:

```powershell
.\run_all.ps1 -Only R -FromSynthetic
```

En ese modo no corre `98_comparacion`, así que los chequeos que cruzan R contra
Python —concordancia de los CSV, paridad de los `.md` y del `informe.html`,
`comparacion_R_python.csv`, `comparacion_reporte.md`— **no se pueden ejecutar**.
No cuentan como fallo ni como aprobación: quedan en un tercer estado,
`NO_EJECUTADA`, listado en la salida y en `logs/corrida_<fecha>.txt`. La corrida
termina con `VERIFICACION PARCIAL: 0 fallos, N chequeo(s) y M fila(s)
NO_EJECUTADA` y con código de salida 0. `TODAS LAS VERIFICACIONES PASARON` exige
cero fallos **y** cero `NO_EJECUTADA`, es decir la corrida completa con las dos
implementaciones.

### PDF del informe

`12_informe` genera `docs/informe.pdf` imprimiendo el HTML con Chrome/Edge headless
(busca `msedge.exe` / `chrome.exe` en las rutas habituales; se puede fijar otro con
la variable de entorno `MIA_LPS_PDF_ENGINE`). Equivale a:

```powershell
& "$env:ProgramFiles (x86)\Microsoft\Edge\Application\msedge.exe" --headless --disable-gpu --print-to-pdf="docs\informe.pdf" --no-pdf-header-footer "docs\informe.html"
```

El pipeline **no falla** si el PDF no se puede generar: deja el HTML y avisa. En ese
caso `99_verificar` trata la ausencia del PDF como aviso, no como falla.

---

## 5. Salidas esperadas

| Ruta | Contenido |
|---|---|
| `outputs/figures/` | Boxplots de expresión (gen × tejido), pSTAT3, ELISA (Acto 1); correlaciones placenta–cerebro por gen, pair plots, dispersión, simulación (Acto 2) |
| `outputs/tables/R/`, `outputs/tables/python/` | Resultados numéricos, nombres idénticos entre implementaciones |
| `outputs/tables/procedencia.csv` | Una fila por figura y por tabla (script, función, entrada, n, origen del código, decisión, alternativa descartada, por qué) |
| `outputs/tables/verificaciones.csv` | Una fila por resultado principal (resultado, verificación, cómo, qué mostró, pasó) |
| `outputs/tables/analisis_descartados.md` | Qué se probó, por qué no funcionó, qué se hizo en su lugar |
| `outputs/tables/comparacion_R_python.csv` | Concordancia numérica R vs Python con tolerancia declarada |
| `docs/informe.html`, `docs/informe.pdf` | Informe completo autocontenido (Acto 1 + Acto 2 + reproducibilidad + análisis descartados + limitaciones + procedencia/verificaciones); figuras en base64. No versionado (`.gitignore`), regenerable |
| `logs/corrida_<fecha>.txt` | Fecha, SO, versiones de R/Python y paquetes, semilla, fuente, concordancia R↔Python, paridad de render, verificaciones pasadas/totales — una sección por lenguaje |

---

## 6. Decisiones metodológicas

Fijadas en `AGENTS.md` §4 (D1–D12) y no se reabren. Puntos salientes: análisis sobre
**−ΔΔCt** (no fold-change); **sin imputación** de no detectados; calibrador = ♀Control por
gen × tejido sobre detectados; cascada ANOVA-III / HC3 / ART según diagnóstico de
residuos; `il6` en cerebro solo como proporción de detección; ELISA con censura a
izquierda tratada con métodos para datos censurados; score compuesto de transportadores
(**nunca** "TONE").
