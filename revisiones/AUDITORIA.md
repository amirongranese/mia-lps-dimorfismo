# Auditoría escéptica del repositorio

Fecha: 2026-09-21. Alcance: sólo lectura; no ejecuté el pipeline ni modifiqué salidas.

## Dictamen

Este clon no permite confiar todavía en resultados numéricos auditados. No contiene `docs/informe.html`, `docs/informe.pdf`, `outputs/tables/verificaciones.csv`, ni CSV bajo `outputs/tables/R` o `outputs/tables/python`: las carpetas sólo contienen `.gitkeep` y las salidas están ignoradas por Git. Dado que se pidió no correr nada que pueda sobrescribir, no pude contrastar números de informe y tablas.

Esto no demuestra que la corrida histórica sea falsa. Demuestra que la frase de corrida exitosa en `ESTADO.md` no aporta evidencia forense disponible en este checkout. Hace falta el paquete inmutable de artefactos de esa corrida, o una ejecución aislada.

## Cinco afirmaciones: trazabilidad, pero sin conciliación posible

| # | Afirmación generada por `12_informe.py` | Cálculo fuente | Tabla y campos fuente | ¿Coincide? |
|---|---|---|---|---|
| 1 | “IL-6 detectable en LPS `x/n` frente a Control `x/n` (Fisher p=`p`).” | `03_elisa.py`: `tabla_fisher()` (l. 358), `fisher_2x2()`; informe: `numeros_conclusiones()` l. 553–561. | `elisa_fisher_deteccion.csv`; fila `bloque=MS`, campos `control_detectado`, `control_n`, `lps_detectado`, `lps_n`, `p_valor`. | **NO VERIFICABLE**. |
| 2 | “En líquido amniótico ningún contraste alcanzó significancia (Peto–Peto ♀/♂).” | `03_elisa.py`: `tabla_petopeto_la()` (l. 296); informe l. 563–568 y 731–743. | `elisa_petopeto_la.csv`; `ESTRATO=HEMBRA/MACHO`, `p_valor`. | **NO VERIFICABLE**. |
| 3 | “`N` genes de cerebro mostraron interacción SEXO×TTO.” | `05_qpcr_modelos.py`: `clasificar()` y `_ajustar_cascada()`; informe l. 587–628, 746–758. | `qpcr_modelos_clasificacion.csv`: `TEJIDO`, `GEN`, `via`, `interaccion_significativa`; post hoc: `qpcr_modelos_posthoc.csv`, `p_holm`. | **NO VERIFICABLE**. |
| 4 | “Ningún Δrho significativo (`x/y`).” | `09_acto2_dispersion.py`: `tabla_test()` (l. 445); informe l. 648–660. | `acto2_test_correlaciones.csv`: `ITEM`, `ESTRATO`, `p_bw`. | **NO VERIFICABLE**; el denominador incluye los tres estratos. |
| 5 | “PC1 explica `x%` en placenta e `y%` en cerebro.” | `11_sensibilidad.py`, PCA; informe l. 691–704. | `acto2_sensibilidad_pca_varianza.csv`: `TEJIDO`, `PC`, `prop_var`; el informe multiplica por 100 y redondea. | **NO VERIFICABLE**. |

Que el informe vuelva a leer los CSV aporta trazabilidad, no validación independiente: un error que escribió un CSV se replica en el HTML.

## Verificaciones que no pueden refutar el comportamiento anunciado

No pude abrir `verificaciones.csv`; no existe. Estas son las filas que generaría el código Python (R tiene equivalentes).

| Script | IDs no discriminantes | Razón |
|---|---|---|
| `03_elisa.py` | `elisa_la_metodo` | Literal `TRUE`; no prueba que se usó Peto–Peto estratificado. |
| `04_qpcr_cuantificacion.py` | `cuant_sd_metodo` | Literal `TRUE`; no recalcula z-scores. |
| `05_qpcr_modelos.py` | `modelos_posthoc_holm`, `modelos_art_core_vs_artool`, `modelos_bh_d12` | El primero prueba sólo `n_posthoc >= 0`; los otros son `TRUE` literal, incluso en Python sin ARTool. |
| `06_pstat3.py` | `pstat3_modelo_d9`, `pstat3_limitacion_d9`, `pstat3_core_vs_libs` | Literales; no buscan la limitación ni validan el cruce declarado. |
| `07_figuras_acto1.py` | `figuras_d11_una_funcion`, `figuras_d11_cascada`, `figuras_fc_log`, `figuras_il6_brain_deteccion`, `figuras_pstat3_membrana`, `figuras_expresion_r_base` | `TRUE` por construcción; no inspeccionan PNG ni brackets. |
| `08_acto2_correlaciones.py` | `acto2_coeficiente`, `acto2_no_compara_grupos`, `acto2_pareo_por_feto` | Literales; no recalculan rho ni pareo. |
| `09_acto2_dispersion.py` | `acto2_test_metodo`, `acto2_test_items`, `acto2_test_bh_suplementario`, `acto2_dispersion_pares`, `acto2_dispersion_levene`, `acto2_dispersion_no_concluye`, `acto2_interaccion_sin_cascada` | Descripciones aprobadas por construcción; `acto2_test_items` no exige el conteo esperado. |
| `10_acto2_simulacion.py` | `acto2_sim_proposito`, `acto2_sim_inversion_rho_r`, `acto2_sim_SD_por_grupo`, `acto2_sim_veredicto`, `acto2_sim_tasa_fisher` | Literales; no comprueban simulaciones ni sus SD. |
| `11_sensibilidad.py` | `sens_eigengene_correlacion_global`, `sens_excl_extremo_criterio`, `sens_excl_extremo_veredicto`, `sens_excl_extremo_no_reestima_pc1` | Literales; el veredicto sólo informa un conteo. |

`acto2_interaccion_vs_levene` también se guarda como `TRUE` literal. Hay un `assert` anterior de diferencia menor que `1e-8`; si falla, el proceso aborta y esa fila no se escribe. El `assert` es útil; la fila no agrega evidencia.

Los siguientes sólo verifican existencia o tamaño, no contenido: `elisa_figura_ms`, `elisa_figura_la`, `figuras_acto1_generadas`, `acto2_figuras`, `acto2_dispersion_figuras`, `acto2_sim_figura`, `sens_figuras`, `informe_html_generado` e `informe_snapshot_paridad`. `informe_pdf` devuelve `TRUE` para `ok`, `sin_motor` y cualquier estado que empiece por `fallo`: un fallo de PDF figura como aprobado.

`elisa_lod`, `modelos_alfa`, `modelos_piso_celda`, `pstat3_alfa`, `acto2_test_piso`, `acto2_estratos`, `acto2_sim_escenarios` y `acto2_sim_B_semilla` comparan una constante con el valor declarado por el mismo módulo. Son metadatos, no pruebas de aplicación.

## Conflictos entre AGENTS.md y código

### C1 — El informe afirma robustez BH sin calcularla para esa frase

D12 dice que BH es suplementario, no cambia conclusiones y que el informe debe decir cuántos resultados sobreviven. En `12_informe.py` l. 606–608, `bra_int_genes` se forma sólo con `interaccion_significativa == "TRUE"` de `qpcr_modelos_clasificacion.csv` (p primario). En l. 746–758 el texto declara que esas interacciones son “**todas robustas a la corrección BH**”. Esa función no lee `p_SEXOxTTO_BH` ni cuenta supervivientes BH.

Es una contradicción entre afirmación y mecanismo. Aunque los datos concretos hicieran verdadera la frase, no está demostrada programáticamente ni informa el numerador/denominador pedido por D12. Presenta resultados primarios como robustos a multiplicidad.

### C2 — La verificación de pSTAT3 fija una rama, no la cascada D5

D5 exige elegir ANOVA III, HC3 o ART según Shapiro y Levene. Sin embargo, `06_pstat3.py` l. 951–954 aprueba `pstat3_rama_cascada` sólo si `rec["rama"] == "anova3"`. Con datos donde la cascada eligiera correctamente HC3 o ART, esa fila sería `FALSE`. A la vez, `05_qpcr_modelos.py` l. 1046–1048 sólo comprueba que los conteos de ramas sumen 18, no que cada gen siguió sus diagnósticos.

La cascada sí está implementada; no es una violación demostrada de la corrida histórica. Es un conflicto entre D5 como regla general y el sentido de estas verificaciones: pueden reprobar una elección correcta distinta de la histórica y aprobar una clasificación errónea cuyo total siga siendo 18.

### C3 — La frase placentaria no se deriva de la tabla

`conclusion_seccion()` (`12_informe.py` l. 746–748) escribe “Ningún gen mostró interacción SEXO×TTO” para placenta. El cálculo que alimenta esa sección, `pla_tto_genes` (l. 601–605), filtra sólo `p_TTO < .05`; no calcula ni lee el número de interacciones placentarias. Una corrida futura con interacción placentaria podría conservar una frase falsa. Sin las tablas no puedo decidir si ya es falsa en la corrida histórica.

### C4 — El lenguaje interpretativo excede los tests

La conclusión revisada sostiene que la pérdida de acoplamiento es “la misma reducción de dispersión” y que en machos “cada individuo responde distinto”. El diseño declara 18 madres y 36 fetos, pero `05_qpcr_modelos.py` ajusta sobre fetos sin término `MADRE` o camada. La simulación evalúa restricción de rango, no identifica causalmente una única explicación de correlación. Esto no viola literalmente D4/D5, pero el informe es más fuerte que la inferencia cubierta y no declara posible dependencia intra-camada.

## Evidencia mínima antes de confiar

1. Entregar por commit HTML, CSV R/Python, `verificaciones.csv`, `procedencia.csv`, comparación R/Python y log, con hashes SHA-256.
2. Conciliar las cinco afirmaciones anteriores contra filas exactas, incluido redondeo.
3. Reemplazar literales `TRUE` por pruebas de salidas o recálculos desde TSV.
4. Resolver C1: leer `p_SEXOxTTO_BH`, informar cuántos sobreviven y no redefinir la conclusión primaria con BH si D12 lo prohíbe.
5. Declarar la independencia fetal asumida o moderar el lenguaje causal.

## Anexo de alcance

Inspeccioné `AGENTS.md`, `ESTADO.md` y scripts Python 03–12, con referencias R cuando hizo falta. No inferí valores numéricos ausentes ni traté frases históricas de `ESTADO.md` como evidencia verificable.
