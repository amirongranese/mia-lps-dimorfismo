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

### 2.1 R (≥ 4.4; probado con 4.6.1)

```powershell
# Desde la raíz del repo
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "install.packages('renv'); renv::restore()"
```

Si todavía no hay `renv.lock` (se genera al cerrar T1), instalar a mano:

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "install.packages(c('car','ARTool','emmeans','NADA','survival','sandwich','lmtest','nlme','rmarkdown','ggplot2','GGally','readxl','writexl','jsonlite'))"
```

Pandoc (para renderizar el informe):

```powershell
& "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "install.packages('installr'); installr::install.pandoc()"
```

### 2.2 Python (≥ 3.11; probado con 3.13)

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

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

Corre, en orden: `00_config` → `01_generar_sinteticos` → `02_ingesta_qc` → `03_elisa` →
`04_qpcr_cuantificacion` → `05_qpcr_modelos` → `06_pstat3` → `07_figuras_acto1` →
`08_acto2_correlaciones` → `09_acto2_dispersion` → `10_acto2_simulacion` →
`11_sensibilidad` → `12_informe` → `99_verificar`, en R y en Python.

`run_all.ps1` acepta flags (ver cabecera del script): `-Only R|python`, `-FromSynthetic`,
`-SkipReport`.

### PDF del informe sin LaTeX

Si `12_informe` no encuentra LaTeX, genera `docs/informe.pdf` imprimiendo el HTML con
Chrome/Edge headless:

```powershell
& "$env:ProgramFiles\Microsoft\Edge\Application\msedge.exe" --headless --disable-gpu --print-to-pdf="docs\informe.pdf" --no-pdf-header-footer "docs\informe.html"
```

El pipeline **no falla** si el PDF no se puede generar: deja el HTML y avisa.

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
| `docs/informe.html`, `docs/informe.pdf` | Informe completo autocontenido (Acto 1 + Acto 2 + conclusión comparada) |
| `logs/corrida_<fecha>.txt` | Fecha, SO, versiones de R/Python y paquetes, semilla, tiempo, verificaciones pasadas/totales |

`99_verificar` imprime `TODAS LAS VERIFICACIONES PASARON` cuando el proyecto está completo.

---

## 6. Decisiones metodológicas

Fijadas en `AGENTS.md` §4 (D1–D12) y no se reabren. Puntos salientes: análisis sobre
**−ΔΔCt** (no fold-change); **sin imputación** de no detectados; calibrador = ♀Control por
gen × tejido sobre detectados; cascada ANOVA-III / HC3 / ART según diagnóstico de
residuos; `il6` en cerebro solo como proporción de detección; ELISA con censura a
izquierda tratada con métodos para datos censurados; score compuesto de transportadores
(**nunca** "TONE").
