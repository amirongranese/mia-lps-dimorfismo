# Reauditoría de reproducción — 2026-09-28

Commit examinado: `c2beacf`. Continuación de la prueba interrumpida en este clon.
Se leyeron AGENTS.md y ESTADO.md completos, los dos informes originales y RESPUESTA.md.
No se modificaron los scripts ni las decisiones metodológicas. No hay datos crudos
en `data/raw/`; las comprobaciones de esta sesión corresponden a datos sintéticos.

## Dictamen

**No está demostrada la reproducción completa desde un clon limpio.** Se confirmó
concordancia byte a byte de los 32 CSV R/Python disponibles, pero el informe Python
falla con los datos sintéticos y no existe su snapshot. Ambos verificadores finales
informan verificación parcial, no `TODAS LAS VERIFICACIONES PASARON`.

Al retomar ya existían las tablas de ambos lenguajes, el informe R y su PDF; no había
log de verificación final ni snapshot Python. No se atribuye a esta sesión la
ejecución completa de los pasos 00–11. Se ejecutaron directamente Python/12 y los
dos verificadores 99 sobre ese estado conservado.

## Hallazgos

### H1 — Bloqueante: el informe Python falla sobre sintéticos

Comando ejecutado, desde la raíz y con `MIA_LPS_FORZAR_SINTETICO=1` y `PYTHONUTF8=1`:

```powershell
& .\.venv\Scripts\python.exe python/12_informe.py
```

Resultado: exit 1, `OverflowError: cannot convert float infinity to integer`.
Traza: `main()` (línea 1356) → `numeros_conclusiones()` (673) → `_round_fmt()` (189).

`disp_holm_min` se inicializa con `math.inf` (661). Las tablas sintéticas actuales
no tienen genes de cerebro con interacción significativa, de modo que el bucle
no cambia ese valor. La línea 673 intenta redondearlo usando `math.floor`.
`main()` calcula estos números **antes** de decidir reemplazar las conclusiones
por el aviso de datos sintéticos. Por eso el aviso no evita la excepción.

R mantiene el mismo cálculo (`R/12_informe.R:529–542`), pero su redondeo de infinito
no produce la misma excepción. Tener éxito con `-Only R -FromSynthetic` no prueba
el éxito de Python ni del modo completo.

Corrección sugerida: no calcular conclusiones biológicas cuando la fuente es
sintética y representar explícitamente los mínimos no estimables para fuentes
reales con conjuntos vacíos. Verificar ambos lenguajes y la paridad del texto.

### H2 — Alta: evidencia de dependencia del estado previo en la paridad

El snapshot R de `comparacion_reporte.md` contiene **94** filas en la auditoría de
verificaciones; el archivo actual contiene **101**. Los dos dicen que todas sus
filas están en TRUE. Esta diferencia está presente en los archivos de este clon.

`auditar_verificaciones()` excluye sólo las filas del propio `98_comparacion`
(`python/98_comparacion.py:261–263`; equivalente R), pero incluye las de
`12_informe`. En el orden de arranque, la primera pasada de 98 precede al primer
informe; la segunda sucede después. ESTADO, sesión 20, ya describe este defecto
y que una segunda corrida lo disimula.

La diferencia observada es suficiente para impedir la paridad de esos reportes
si el snapshot Python copiara el archivo actual. **No se completó esa comparación
de snapshots**, porque H1 impidió crear el de Python. No se presenta como una
segunda ejecución completa fallida realizada en esta sesión.

Corrección sugerida: definir una población de filas estable entre pasadas y probar
el primer arranque sin salidas previas. Repetir sobre una carpeta ya poblada no
demuestra el requisito de AGENTS §1.

### H3 — Media: snapshot ausente se trata como omisión, con motivo incorrecto

Se ejecutaron `python/99_verificar.py` y `R/99_verificar.R`, con fuente sintética
y sin modo de implementación única. Ambos dieron:

```text
modo = R + Python
checklist: 76/77 chequeos duros OK; 1 NO_EJECUTADA
CSV R<->Python: 32; byte-identicos 32; fuera de tol 0
paridad render: NO_EJECUTADA
VERIFICACION PARCIAL: 0 fallos, 1 chequeo(s) y 0 fila(s) NO_EJECUTADA
```

Ambos terminaron con exit 0. El motivo impreso fue «no hay salidas de la contraparte
en outputs/tables/ para comparar», aunque **sí existen las 32 tablas de cada lado**.
Lo ausente es `outputs/intermediate/render/python/`.

Ubicación: `python/99_verificar.py:330–354`, `R/99_verificar.R:294–322`.
El verificador no afirma que todo pasó: eso está corregido. Sin embargo, ejecutado
directamente no distingue este informe fallido de una comparación legítimamente
no ejecutable. El lanzador sí abortaría ante el exit 1 de H1; no se afirma que
`run_all.ps1` oculte ese fallo.

Corrección sugerida: distinguir modo único de artefactos faltantes en modo completo
y dar un motivo específico para el snapshot ausente.

### H4 — Media: el desglose por tipo sigue incompleto

Las 108 filas actuales sí están rotuladas: 50 `recalculo`, 9 `existencia` y 49
`declaracion`. Pero `comparacion_reporte.md` §4 sigue informando sólo totales
agregados (94 o 101 en TRUE), sin desglose. Lo generan
`python/98_comparacion.py:497` y `R/98_comparacion.R:457`.
La salida de consola de ambos 99 tampoco muestra el desglose de esas 108 filas.
Esto contradice el alcance de cierre documentado en ESTADO/RESPUESTA y deja
incompleto el punto 2.2 del pedido. No implica que todas esas filas sean pruebas
independientes ni invalida la clasificación que sí se agregó.

## Evidencia y límites

- Logs nuevos: `logs/auditoria_2026-09-28_python99.txt`,
  `logs/auditoria_2026-09-28_R99.txt`, `logs/corrida_2026-09-28.txt`.
- Ambos 99 confirman 32 CSV byte-idénticos, cero celdas fuera de tolerancia.
- El PDF existente tiene estado R `ok`; no se efectuó una auditoría visual del PDF.
- El intérprete de `.venv` requirió ejecución fuera del aislamiento por acceso
  denegado al Python instalado; con autorización se obtuvo la excepción H1.
- R emitió avisos de locale, sin impedir la ejecución del verificador.
- No se conciliaron resultados biológicos reales: esos datos no están en este clon.
- No se versionaron `docs/`, `outputs/`, logs ni datos. Los seis TSV sintéticos
  ya figuraban modificados al retomar. Este informe se entrega como archivo,
  sin commit solicitado. Se retiró la entrada que se había agregado a ESTADO.md.

Hashes SHA-256 de evidencia local al auditar:

| Archivo | SHA-256 |
|---|---|
| `python/12_informe.py` | `AD025A3ADA27D60A2DA9AD429ED0B2E4FE9A92E8E68ED2CA80638F6CE2C10FF1` |
| `outputs/tables/python/qpcr_modelos_clasificacion.csv` | `1FA49DF7183F70FD14C69A36851A05CC6A48F828956D2E6F066630A3F71F0475` |
| `outputs/tables/python/qpcr_modelos_posthoc.csv` | `6E04461C92E9F425430D650784046284B990A4D2150EDE4D82A8472621371A59` |
| `outputs/intermediate/render/R/tables/comparacion_reporte.md` | `9D5EA99EC1F3D2B728BDDC066083FFEE7C3DD706A047D6D446090D428E384436` |
| `outputs/tables/comparacion_reporte.md` | `F605E767DDDF92C19CAA4ECEB42E8429130DBB712C95D53D61DC5F023CE77A4D` |

Los hashes identifican archivos de esta prueba; no sustituyen un paquete de
artefactos accesible a terceros. La decisión de publicación sigue pendiente.

## Siguiente paso

Corregir H1 y H2 en ambas implementaciones, completar H3/H4 y repetir el arranque
sintético desde salidas vacías en un directorio aislado. Conservar por separado la
prueba de `-Only R`. El cierre histórico sobre datos reales no sustituye esa prueba.
