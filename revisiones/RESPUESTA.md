# Respuesta a la revisión externa

Punto 6 de `pedidos/cambios_revision_codex.md`. Fecha: 2026-09-28.

Una fila por cada hallazgo de los dos informes de la revisión —
`revisiones/AUDITORIA.md` (auditoría escéptica) y `revisiones/PRUEBA_REPRODUCCION.md`
(prueba de reproducción)— con qué se hizo y en qué commit. **Los hallazgos que no se
resolvieron también están, con el motivo.**

Ningún cambio de esta tanda tocó los datos, los análisis ni los resultados: D1–D12 no
se reabrieron, no se agregó ningún test nuevo y los números del informe son los mismos.
Lo que cambió es documentación, calidad de las verificaciones, derivación del texto
desde las tablas y reproducibilidad.

## Commits

| Commit | Contenido |
|---|---|
| `a42e792` | Punto 1 (D13 + moderar lenguaje causal) y punto 3 (texto derivado de las tablas) |
| `8ffbec2` | Punto 2 (columna `tipo` en `verificaciones.csv` + 4 checks arreglados) |
| `281891f` | Punto 4 (reproducibilidad: intérpretes, `-Only`, `NO_EJECUTADA`, README) |
| `48e89df` | Punto 5 (`revisiones/conciliacion_auditoria.md` + script de conciliación) |
| `527ee01` | Residuo del punto 2: las 3 menciones de «100/100 verificaciones» que quedaban en `ESTADO.md` |
| este commit | Punto 6 (este archivo) |

## Estados usados

- **Resuelto** — el hallazgo ya no aplica: el código o el texto cambió.
- **Resuelto parcialmente** — se corrigió lo que el pedido pidió corregir, y queda
  declarado qué parte no.
- **Declarado como limitación** — no se cambia el análisis; la limitación se escribe
  donde corresponde (AGENTS, informe, `analisis_descartados.md`).
- **No aplica** — no era un defecto del repositorio, o es un artefacto del entorno de
  prueba.
- **Pendiente** — reconocido, sin resolver todavía, con el motivo.

---

## A. `AUDITORIA.md` — auditoría escéptica

| # | Hallazgo | Qué se hizo | Commit |
|---|---|---|---|
| A1 | **Dictamen**: el clon no contiene `docs/informe.html`, `informe.pdf`, `verificaciones.csv` ni CSV en `outputs/tables/{R,python}/` (están en `.gitignore`), así que no se pudo contrastar ningún número. | **Pendiente.** Es correcto: las salidas no se versionan. Qué se publica depende de una decisión del usuario que sigue abierta; el pedido lo dice explícitamente («No se versiona todavía nada de `docs/` ni de `outputs/`»). Mientras tanto, el punto 5 concilia las cinco afirmaciones contra las tablas de esta máquina y deja el script para repetirlo. | `48e89df` |
| A2 | **Las cinco afirmaciones**, todas marcadas `NO VERIFICABLE` por falta de tablas. | **Resuelto.** `revisiones/conciliacion_auditoria.md`: las cinco conciliadas contra la fila exacta del CSV, incluido el redondeo. Las cinco **coinciden**; ninguna requirió corregir el informe. | `48e89df` |
| A2.4 | Nota sobre la afirmación 4: «el denominador incluye los tres estratos». | **Resuelto** (confirmado). 27 = 9 ítems × 3 estratos (`AMBOS_SEXOS`, `HEMBRA`, `MACHO`), y el informe lo dice en la misma frase («ni agrupando sexos ni dentro de cada sexo»). | `48e89df` |
| A2.5 | Nota sobre la afirmación 5: «el informe multiplica por 100 y redondea». | **Resuelto** (confirmado). `.round_fmt` redondea medio hacia arriba: 81.551649 → 82 (truncar daría 81). El informe dice 82. | `48e89df` |
| A3 | «Que el informe vuelva a leer los CSV aporta trazabilidad, no validación independiente: un error que escribió un CSV se replica en el HTML.» | **Declarado como limitación.** Es cierto y no se puede resolver leyendo tablas. Lo que sí hay: 50 verificaciones de tipo `recalculo` y la concordancia celda a celda entre dos implementaciones escritas por separado (`98_comparacion` / `99_verificar`, 32 CSV byte-idénticos). La limitación queda escrita en `revisiones/conciliacion_auditoria.md` §«Lo que esta conciliación no demuestra». | `48e89df` |
| A4 | **Verificaciones no discriminantes**: buena parte de `verificaciones.csv` son `TRUE` literales que no pueden fallar (lista por script: `elisa_la_metodo`, `cuant_sd_metodo`, `modelos_posthoc_holm`, `modelos_art_core_vs_artool`, `modelos_bh_d12`, `pstat3_modelo_d9`, `pstat3_limitacion_d9`, `pstat3_core_vs_libs`, las 6 de `figuras_*`, `acto2_coeficiente`, `acto2_no_compara_grupos`, `acto2_pareo_por_feto`, las 7 de `09`, las 5 de `10`, las 4 de `11`). | **Resuelto.** Columna `tipo` en `verificaciones.csv` con tres valores: `recalculo` / `existencia` / `declaracion`. Las 108 filas actuales quedaron **50 `recalculo` + 9 `existencia` + 49 `declaracion`**. El pedido no exigió convertirlas en recálculos («alcanza con que estén rotuladas como lo que son»); las que sí estaban mal se arreglaron (A6–A9). | `8ffbec2` |
| A5 | El informe dice «100/100 verificaciones» sin desglose, lo que sobrestima la evidencia. | **Resuelto.** La expresión sin desglose desapareció del repo. `12_informe` §10 y `98_comparacion` reportan recalculo / existencia / declaración por separado, más el aviso de `NO_EJECUTADA` cuando corresponde. **Residuo encontrado al escribir esta respuesta:** la sesión 20 verificó el grep sólo contra el README y había dejado tres menciones en sus propias entradas de bitácora de `ESTADO.md`; reformuladas (no se puede desglosar por tipo una corrida anterior a que existiera la columna `tipo`, así que se dice «100 filas en TRUE, sin desglose por tipo» en vez de inventar uno). | `8ffbec2`, `527ee01` |
| A6 | `informe_pdf` devuelve `TRUE` para `ok`, `sin_motor` y cualquier estado `fallo:*`: un fallo de PDF figura como aprobado. | **Resuelto.** `ok` → `TRUE`; `sin_motor` → `NO_EJECUTADA`; `fallo:*` → `FALSE`. | `8ffbec2` |
| A7 | **C2** — `pstat3_rama_cascada` aprueba sólo si la rama es `anova3`: puede reprobar una elección correcta distinta de la histórica. | **Resuelto.** Recalcula la rama esperada desde los p de Shapiro-Wilk y Levene registrados y la compara con la elegida, cualquiera sea. | `8ffbec2` |
| A8 | **C2** (segunda parte) — `05_qpcr_modelos` sólo comprueba que los conteos de ramas sumen 18: aprueba una clasificación errónea cuyo total siga siendo 18. | **Resuelto.** `modelos_rama_conteo` valida **gen por gen** que la rama usada corresponde a sus propios diagnósticos. | `8ffbec2` |
| A9 | `acto2_test_items` no exige el número esperado de ítems, sólo lo reporta. | **Resuelto.** Ahora exige `n_test == length(ITEMS)`. | `8ffbec2` |
| A10 | `acto2_interaccion_vs_levene` es `TRUE` literal; el `assert` previo (dif. < 1e-8) es útil, la fila no agrega evidencia. | **Resuelto parcialmente.** Reclasificada como `declaracion`, que es exactamente lo que es. El `assert` se mantiene (aborta el script si falla, así que sí protege). No se convirtió en `recalculo`: la evidencia real está en el `assert`, y rotularla `recalculo` sería el problema que el hallazgo señala. | `8ffbec2` |
| A11 | Nueve filas verifican sólo existencia o tamaño, no contenido: `elisa_figura_ms`, `elisa_figura_la`, `figuras_acto1_generadas`, `acto2_figuras`, `acto2_dispersion_figuras`, `acto2_sim_figura`, `sens_figuras`, `informe_html_generado`, `informe_snapshot_paridad`. | **Resuelto.** Son las 9 filas rotuladas `existencia`. Se informan aparte en el desglose del informe. | `8ffbec2` |
| A12 | Ocho filas comparan una constante con el valor declarado por el mismo módulo (`elisa_lod`, `modelos_alfa`, `modelos_piso_celda`, `pstat3_alfa`, `acto2_test_piso`, `acto2_estratos`, `acto2_sim_escenarios`, `acto2_sim_B_semilla`): son metadatos, no pruebas de aplicación. | **Resuelto.** Las ocho quedaron rotuladas `declaracion`. (`acto2_estratos` estaba mal clasificada como `recalculo` en el primer pase y se corrigió antes de commitear: `ESTRATOS` ahí es una constante del módulo, no algo derivado de los datos.) | `8ffbec2` |
| A13 | **C1** — El informe declara que las interacciones de cerebro son «todas robustas a la corrección BH» sin leer `p_SEXOxTTO_BH` ni informar numerador/denominador, como pide D12. | **Resuelto.** La frase se calcula desde la columna: «**7 de 7** sobreviven a BH». BH sigue siendo suplementario (D12): el conteo se informa, no redefine qué es significativo. Verificado en el punto 5: los 7 tienen `p_SEXOxTTO_BH = 2.992194e-02`. | `a42e792` |
| A14 | **C3** — «Ningún gen mostró interacción SEXO×TTO» en placenta es texto fijo; el cálculo que alimenta esa sección sólo filtra `p_TTO < .05`. Una corrida futura con interacción placentaria conservaría una frase falsa. | **Resuelto.** La frase se deriva del conteo real de interacciones en placenta (`pla_int_n`); si en una corrida futura hubiera alguna, la frase cambia sola. | `a42e792` |
| A15 | **C4** — Los modelos no tienen término `MADRE` ni de camada, aunque el diseño declara 18 madres y 36 fetos, y el informe no declara posible dependencia intra-camada. | **Declarado como limitación.** **No se agrega `MADRE` al modelo** (decisión del usuario, ya tomada, que no estaba escrita). Ahora es **D13** en `AGENTS.md`/`CLAUDE.md`, tiene entrada en `analisis_descartados.md` y una limitación nueva en el informe: el LPS se administra a la madre, cada camada aporta un feto de cada sexo, el análisis asume independencia, y si hubiera variación entre camadas afectaría sobre todo a la precisión de los efectos principales de tratamiento. | `a42e792` |
| A16 | **C4** (segunda parte) — El lenguaje interpretativo excede los tests: «cada individuo responde distinto», «no es un segundo hallazgo, sino la misma reducción de dispersión». | **Resuelto.** Las dos frases se moderaron según el pedido, y se aplicó el mismo criterio a **3 frases más** con el mismo patrón («no son efectos independientes», «no dependen de»): 5 en total. Donde los tests sólo muestran compatibilidad, el texto dice compatibilidad. | `a42e792` |
| A17 | **Evidencia mínima 1** — Entregar por commit HTML, CSV R/Python, `verificaciones.csv`, `procedencia.csv`, comparación y log, con hashes SHA-256. | **Pendiente.** Misma razón que A1: depende de la decisión de publicación, que el pedido deja explícitamente fuera de esta tanda. | — |
| A18 | **Evidencia mínima 2** — Conciliar las cinco afirmaciones contra filas exactas, incluido redondeo. | **Resuelto** (= A2). | `48e89df` |
| A19 | **Evidencia mínima 3** — Reemplazar literales `TRUE` por pruebas de salidas o recálculos desde TSV. | **Resuelto parcialmente, por decisión explícita del usuario.** El pedido acota el alcance: «No hace falta reescribir todas las `declaracion` como recálculos: alcanza con que estén rotuladas como lo que son». Se rotularon las 108 filas y se arreglaron de verdad las 4 que estaban mal (A6–A9). Las 49 `declaracion` siguen siendo declaraciones, ahora dichas como tales. | `8ffbec2` |
| A20 | **Evidencia mínima 4** — Resolver C1 sin redefinir la conclusión primaria con BH. | **Resuelto** (= A13). | `a42e792` |
| A21 | **Evidencia mínima 5** — Declarar la independencia fetal asumida o moderar el lenguaje causal. | **Resuelto** — se hizo lo uno **y** lo otro (= A15 + A16). | `a42e792` |

---

## B. `PRUEBA_REPRODUCCION.md` — prueba de reproducción

| # | Hallazgo | Qué se hizo | Commit |
|---|---|---|---|
| B1 | No había intérprete de Python: `python`, `py`, `python3` y `winget` no estaban disponibles. El README exige Python ≥ 3.11 pero no dice cómo instalarlo ni de qué fuente. | **Resuelto.** Averiguado y documentado: en esta máquina es **CPython 3.13 per-user de python.org**, en `%LOCALAPPDATA%\Programs\Python\Python313\python.exe`, presente en el `PATH`. README §2.2 explica versión, fuente e instalación (con la advertencia de no usar el alias de la Microsoft Store), y §2.0 documenta los cuatro lugares donde el pipeline lo busca. | `281891f` |
| B2 | `./run_all.ps1` bloqueado por la política de ejecución de PowerShell; el README no contempla esa situación. | **Resuelto.** README §2.3: `powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_all.ps1`, con la aclaración de que vale sólo para ese proceso y no cambia la política del sistema. | `281891f` |
| B3 | `run_all.ps1 -Only R` aborta con «No se encontro python en PATH», aunque ese modo no ejecuta Python: la comprobación se hace antes de interpretar `-Only`. | **Resuelto.** Cada intérprete se resuelve y se exige **sólo si el modo lo va a usar**. Con `-Only R` el script imprime `Python: (no se usa en este modo)` y no busca Python. | `281891f` |
| B4 | Reproducción completa, `98_comparacion` y `99_verificar` no ejecutables sin Python: no se pueden crear las salidas gemelas ni la comparación. Quien tiene una sola implementación no puede verificar nada. | **Resuelto.** `-Only R|python` ahora **corre `99_verificar`**. Los chequeos que cruzan R contra Python (concordancia de los CSV, paridad de `.md`/`informe.html`, `comparacion_R_python.csv`, `comparacion_reporte.md`, CSV de la contraparte) pasan a un tercer estado `NO_EJECUTADA`: ni pase ni fallo, listados en consola y en `logs/corrida_<fecha>.txt`. La salida final es `VERIFICACION PARCIAL: 0 fallos, N chequeo(s) y M fila(s) NO_EJECUTADA`, con código de salida 0. `TODAS LAS VERIFICACIONES PASARON` exige cero fallos **y** cero `NO_EJECUTADA`. | `281891f` |
| B5 | `renv::restore()` falló sin acceso de red a CRAN, y la caché de `renv` vive fuera del repositorio (`AppData\Local\R\cache\R\renv`); el README no lo aclaraba. | **Resuelto.** README §2.1: la primera restauración necesita red a CRAN (~80 paquetes) y la caché está en `%LOCALAPPDATA%\R\cache\R\renv`, fuera del repo, con permiso de escritura. | `281891f` |
| B6 | El README usa una ruta fija de R (`C:\Program Files\R\R-4.6.1\...`). La prueba corrió en esta misma máquina, así que no pudo detectarlo; en otra computadora fallaría. (Hallazgo señalado por el usuario en el pedido, no por la prueba.) | **Resuelto.** `run_all.ps1` resuelve `Rscript` por `PATH` → registro `R-core` → instalaciones estándar, y `python` por `.venv` del repo → `PATH` → lanzador `py -3` → instalaciones estándar; cada candidato se valida ejecutándolo. Las rutas del README quedaron como ejemplos, y `AGENTS.md` §8 lo aclara. **Verificado en esta máquina: `Rscript` no está en el `PATH` y se resuelve igual** (por el registro). | `281891f` |
| B7 | «Ejecutar los scripts R uno por uno fue una inferencia operativa para poder probar la parte R.» | **Resuelto** (= B3 + B4): ya no hace falta, `-Only R -FromSynthetic` corre y verifica la parte R de punta a punta. Documentado en README §4 «Verificar con una sola implementación». | `281891f` |
| B8 | No se generaron `outputs/tables/python/`, `comparacion_R_python.csv` ni `logs/corrida_<fecha>.txt`; `verificaciones.csv` y `procedencia.csv` se escribieron pero no pasaron por la compuerta final. | **Resuelto** como consecuencia de B4: con una sola implementación ahora se escribe el log y se pasa por `99_verificar`, y los tres artefactos que dependen de la pasada gemela quedan visibles como `NO_EJECUTADA` en vez de faltar sin explicación. | `281891f` |
| B9 | Avisos de locale (`LC_COLLATE`, `LC_CTYPE`, `LC_MONETARY`, `LC_TIME`) en R; no detuvieron ningún script ni afectaron las salidas. | **No aplica.** Es un aviso del entorno de prueba, no del repositorio. Los ordenamientos que importan usan `method = "radix"` (independiente del locale) precisamente para que el locale no cambie las salidas, y la paridad byte a byte R↔Python lo confirma en cada corrida. | — |
| B10 | «La frase de corrida exitosa en `ESTADO.md` no aporta evidencia forense disponible en este checkout.» | **Resuelto parcialmente.** Correcto como crítica de la evidencia. Lo que cambió: `ESTADO.md` ya no dice «100/100 verificaciones» sin desglose, y el log de corrida informa modo, desglose por tipo y `NO_EJECUTADA`. Lo que falta es el paquete de artefactos con hashes (A1/A17), pendiente de la decisión de publicación. | `8ffbec2`, `527ee01`, `281891f` |

---

## Pendientes conocidos, fuera del alcance de esta tanda

1. **Paquete inmutable de artefactos con hashes SHA-256** por commit (A1, A17, B10):
   depende de la decisión de qué se publica de `docs/` y `outputs/`, que el pedido deja
   explícitamente para después.
2. **`<<< COMPLETAR: evidencia de analisis previos >>>`** en
   `outputs/tables/analisis_descartados.md` (sección `05_qpcr_modelos`, D13): el marcador
   está puesto a propósito y lo tiene que completar el usuario con el dato concreto de
   qué análisis previo mostró que `MADRE` no modificaba los resultados. No se inventó
   ningún número.
3. **Las 49 verificaciones de tipo `declaracion`** siguen siendo declaraciones (A19):
   rotuladas como tales por decisión explícita del pedido, no convertidas en recálculos.
