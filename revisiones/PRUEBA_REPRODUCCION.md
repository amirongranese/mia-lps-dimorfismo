# Prueba de reproducción — 2026-09-21

## Resultado

La reproducción completa R + Python **no pudo finalizar** en este clon porque no hay
un intérprete de Python disponible. Sí se reprodujo íntegramente la implementación R
con datos sintéticos: 32 tablas numéricas, 26 figuras, `docs/informe.html` y
`docs/informe.pdf`.

No se modificó código, documentación ni configuración del repositorio. Los únicos
archivos modificados fuera de este informe son los que regeneró el propio pipeline
(datos sintéticos y salidas).

## Pasos seguidos

1. Leí `AGENTS.md`, `ESTADO.md` y `README.md` antes de ejecutar el pipeline.
2. Comprobé las herramientas:
   - R 4.6.1 estaba instalado en la ruta indicada en README.
   - `python`, `py`, `python3` y `winget` no estaban disponibles.
3. Ejecuté el comando de instalación de R documentado:

   ```powershell
   & "C:\Program Files\R\R-4.6.1\bin\Rscript.exe" -e "install.packages('renv'); renv::restore()"
   ```

   Se instaló `renv` 1.2.4 y se restauraron los 80 paquetes requeridos por el
   `renv.lock`.
4. Intenté el comando de arranque literal, `./run_all.ps1 -Only R -FromSynthetic`.
   PowerShell lo bloqueó por la política de ejecución.
5. Lo reintenté sin cambiar la política permanente, mediante
   `powershell.exe -NoProfile -ExecutionPolicy Bypass -File .\run_all.ps1 -Only R -FromSynthetic`.
   El script se detuvo antes de ejecutar R porque exige encontrar Python incluso con
   `-Only R`.
6. Para ejecutar la parte posible, corrí directamente la secuencia R
   `00_config` a `12_informe` con `MIA_LPS_FORZAR_SINTETICO=1`. `renv` necesitó usar
   su caché habitual de usuario para ello. Los scripts 10--12 se ejecutaron también
   de forma individual para aislar la finalización de Acto 2, sensibilidad e informe.

## Inferencias necesarias que no están documentadas

- El README exige Python >=3.11 (probado con 3.13), pero no explica cómo instalarlo
  ni proporciona un intérprete, instalador, enlace de descarga o alternativa si
  falta. No instalé uno de una fuente externa porque no estaría indicado "usando solo
  lo que está escrito en el repo".
- Para sortear la política local que bloqueaba `./run_all.ps1`, usé
  `-ExecutionPolicy Bypass` solo en el proceso invocado. El README no contempla esa
  situación.
- `run_all.ps1 -Only R` está documentado como una pasada de una sola implementación,
  pero realiza la comprobación de Python antes de interpretar `-Only`. Ejecutar los
  scripts R uno por uno fue una inferencia operativa para poder probar la parte R.
- La caché de `renv` se ubica fuera del repositorio (`AppData\Local\R\cache\R\renv`);
  el README no aclara que el entorno aislado debe permitir esa escritura.

## Fallos y causa observada

| Punto | Resultado | Causa |
|---|---|---|
| Primera restauración de R | Falló | El entorno aislado no tenía acceso de red a CRAN. Con acceso al repositorio oficial, la restauración terminó correctamente. |
| `./run_all.ps1` | Falló | La política de ejecución local de PowerShell prohíbe cargar scripts. |
| `run_all.ps1 -Only R -FromSynthetic` con bypass temporal | Falló | `run_all.ps1` aborta con `No se encontro python en PATH`, aunque el modo `-Only R` no debería ejecutar Python. |
| Reproducción completa / `98_comparacion` / `99_verificar` | No ejecutable | No existe intérprete Python; por eso no se pueden crear las salidas gemelas ni hacer la comparación R↔Python. |

También aparecieron avisos de locale (`LC_COLLATE`, `LC_CTYPE`, `LC_MONETARY`,
`LC_TIME`); no detuvieron ningún script R ni afectaron las salidas generadas.

## Salidas generadas y verificadas

- `data/synthetic/`: regenerado por `01_generar_sinteticos.R` con la semilla del
  proyecto.
- `data/processed/`: intermedios de ingesta, cuantificación y modelos.
- `outputs/tables/R/`: 32 CSV no vacíos, incluidos QC, ELISA, cuantificación qPCR,
  modelos, pSTAT3, correlaciones, dispersión, simulación y sensibilidad.
- `outputs/figures/`: 26 PNG no vacíos (Acto 1 y Acto 2).
- `docs/informe.html`: generado, 4,737,771 bytes.
- `docs/informe.pdf`: generado, 8,550,645 bytes.

No se generaron `outputs/tables/python/` (más allá de `.gitkeep`),
`outputs/tables/comparacion_R_python.csv` ni `logs/corrida_<fecha>.txt`, porque
requieren la pasada Python y la verificación final que no fue posible realizar.
`outputs/tables/verificaciones.csv` y `procedencia.csv` sí fueron escritos durante
la pasada R, pero no se sometieron a la compuerta final `99_verificar`.
