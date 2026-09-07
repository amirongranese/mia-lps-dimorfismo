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
