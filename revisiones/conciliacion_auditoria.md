# Conciliación de las cinco afirmaciones de `AUDITORIA.md`

Punto 5 de `pedidos/cambios_revision_codex.md`.

`revisiones/AUDITORIA.md` rastreó cinco afirmaciones del informe hasta el cálculo y la
tabla que las produce, pero las marcó todas como **NO VERIFICABLE**: las salidas no
están versionadas, así que en el clon no había CSV con los que contrastar los números.
Acá se hace ese contraste.

- **Fecha de la conciliación:** 2026-09-28.
- **Corrida conciliada:** datos **reales** (`fuente = real`), la que dejó
  `docs/informe.html` y `outputs/tables/R/` en esta máquina.
- **Cómo se hizo:** cada número se recalculó desde el CSV (`csv.DictReader`, no parseo
  de comas) y se buscó la frase **literal** resultante dentro del texto del informe
  (`outputs/intermediate/render/R/informe.textonly.html`, sin etiquetas, con las
  entidades HTML resueltas). Si el número del informe fuera distinto del de la tabla,
  la frase no aparecería y el verdicto sería `NO COINCIDE`.
- **Reproducible:** `python revisiones/conciliar_afirmaciones.py` desde la raíz del
  repo, después de una corrida. Es una herramienta de auditoría de un solo uso, no
  parte del pipeline.

## Resultado

| # | Afirmación del informe | Tabla y fila | Valor del informe | Valor de la tabla | Redondeo | ¿Coincide? |
|---|---|---|---|---|---|---|
| 1 | «la IL-6 fue detectable en 9/9 madres LPS frente a 1/5 control (Fisher p = 4.995005e-03)» | `elisa_fisher_deteccion.csv`, fila `bloque=MS`, `estrato=(sin estrato)` | 9/9, 1/5, p = 4.995005e-03 | `control_detectado=1`, `control_n=5`, `lps_detectado=9`, `lps_n=9`, `p_valor=4.995005e-03` | ninguno (el informe copia el string `%.6e` del CSV) | **SÍ** |
| 2 | «En líquido amniótico ningún contraste alcanzó significancia (Peto-Peto ♀ p = 3.017712e-01; ♂ p = 3.329216e-01)» | `elisa_petopeto_la.csv`, filas `ESTRATO=HEMBRA/MACHO` | ♀ 3.017712e-01; ♂ 3.329216e-01 | `p_valor` = 3.017712e-01 (♀), 3.329216e-01 (♂) | ninguno | **SÍ** |
| 3 | «De 18 gen × tejido modelados, la interacción SEXO×TTO es significativa en 7 de cerebro E15 y 0 de placenta E15; en cerebro las 7 sobreviven la corrección BH» / «Placenta. Ningún gen mostró interacción SEXO×TTO» | `qpcr_modelos_clasificacion.csv`, columnas `via`, `TEJIDO`, `interaccion_significativa`, `p_SEXOxTTO_BH` | 18 modelados, 7 cerebro, 0 placenta, 7/7 BH | `via=modelo`: 18 filas (10 placenta + 8 cerebro). `interaccion_significativa=TRUE`: 7 en cerebro, 0 en placenta. `p_SEXOxTTO_BH = 2.992194e-02` en las 7 → 7 < 0.05 | conteos, sin redondeo | **SÍ** |
| 4 | «Ningún Δρ Control vs LPS es significativo, ni agrupando sexos ni dentro de cada sexo (0 de 27; mínimo fatcd36 en hembras, p = 7.504615e-02)» | `acto2_test_correlaciones.csv`, columna `p_bw` | 0 de 27; mínimo fatcd36 ♀ p = 7.504615e-02 | 27 filas = 9 ítems × 3 estratos; `p_bw < 0.05` en 0; mínimo `p_bw=7.504615e-02` en `ITEM=fatcd36`, `ESTRATO=HEMBRA` | ninguno | **SÍ** |
| 5 | «El eigengene (PC1) explica el 76 % de la varianza de los transportadores en placenta y el 82 % en cerebro» | `acto2_sensibilidad_pca_varianza.csv`, filas `PC=1` | 76 % y 82 % | `prop_var` = 0.7632848468 (placenta) y 0.8155164853 (cerebro) | ×100 y redondeo a 0 decimales, medio hacia arriba (`.round_fmt`): 76.328485 → **76**; 81.551649 → **82** | **SÍ** |

Las cinco coinciden. Ninguna requirió corregir el informe.

## Detalle por afirmación

### 1 — Fisher de detección de IL-6 en suero materno

Fila exacta de `outputs/tables/R/elisa_fisher_deteccion.csv`:

```
bloque,estrato,control_detectado,control_n,lps_detectado,lps_n,or_haldane,p_valor
MS,(sin estrato),1,5,9,9,0.01754385965,4.995005e-03
```

El informe no redondea: imprime el campo `p_valor` tal como está en el CSV, que ya
viene formateado `%.6e` por `03_elisa` (convención del repo para los p que pasan por
funciones trascendentes). Las proporciones `9/9` y `1/5` son los cuatro conteos de la
fila, en el orden LPS-primero que usa la frase.

### 2 — Líquido amniótico: ningún contraste alcanza significancia

Filas exactas de `outputs/tables/R/elisa_petopeto_la.csv`:

```
estrato,n_control,n_lps,cens_control,cens_lps,chisq,df,p_valor,metodo
HEMBRA,5,9,3,4,1.066350711,1,3.017712e-01,Peto-Peto G-rho=1 (survdiff) sobre Conc reflejada
MACHO,5,9,0,3,0.9375,1,3.329216e-01,Peto-Peto G-rho=1 (survdiff) sobre Conc reflejada
```

La auditoría anotó que el denominador de la afirmación no está explícito. Los dos p
citados son los del contraste primario (Peto-Peto por sexo). Para dejarlo cerrado se
comprobaron **los cuatro** contrastes que el pipeline corre sobre LA: los dos
Peto-Peto de arriba más el Fisher de detección por sexo
(`elisa_fisher_deteccion.csv`, `bloque=LA`: ♀ p = 1.000000e+00, ♂ p = 2.582418e-01).
**0 de 4 con p < 0.05**, así que «ningún contraste alcanzó significancia» vale también
con el denominador amplio.

### 3 — Interacción SEXO×TTO por tejido

Derivado de `outputs/tables/R/qpcr_modelos_clasificacion.csv`:

- `via = modelo` en 18 filas: 10 de `PLACENTA_E15` y 8 de `BRAIN_E15`. Las dos que
  faltan para 20 son `BRAIN_E15/il6` (`via = D7_deteccion`, 0 detectados) y
  `BRAIN_E15/il6R` (`via = descriptivo_n_bajo`, celda con < 5 detectados) — los dos
  casos que D7 y el piso de celda sacan del modelo factorial.
- `interaccion_significativa = TRUE`: **7 en cerebro** (fatcd36, fatp1, fatp4, glut1,
  gp130, slc38a1, slc38a2) y **0 en placenta**.
- El único gen de cerebro modelado sin interacción es `glut3`, con
  `p_SEXOxTTO = 5.216446e-02` — justo del lado no significativo de α = 0.05. Es el
  caso que hace que la frase valga la pena verificar y no dar por obvia.
- BH entre genes (D12, suplementario): los 7 tienen
  `p_SEXOxTTO_BH = 2.992194e-02` < 0.05 → **7 de 7 sobreviven**. Es el numerador y el
  denominador que pedía C1 de la auditoría, y desde el punto 3.1 del pedido se cuenta
  desde la columna, no se afirma.

### 4 — Δρ Control vs LPS

`outputs/tables/R/acto2_test_correlaciones.csv` tiene 27 filas = 9 ítems (8 genes +
`score_compuesto`) × 3 estratos (`AMBOS_SEXOS`, `HEMBRA`, `MACHO`). La auditoría
señaló que «el denominador incluye los tres estratos»: es correcto, y el informe lo
dice en la misma frase («ni agrupando sexos ni dentro de cada sexo»). Ninguna fila
tiene `p_bw < 0.05`. El mínimo es:

```
ITEM=fatcd36  ESTRATO=HEMBRA  delta_rho=0.7571428571  p_bw=7.504615e-02
```

que es el que el informe cita como «mínimo fatcd36 en hembras».

### 5 — Varianza explicada por PC1

Filas `PC=1` de `outputs/tables/R/acto2_sensibilidad_pca_varianza.csv`:

```
TEJIDO,PC,autovalor,prop_var,prop_var_acum
PLACENTA_E15,1,5.342993928,0.7632848468,0.7632848468
BRAIN_E15,1,5.708615397,0.8155164853,0.8155164853
```

El informe multiplica por 100 y redondea a 0 decimales con `.round_fmt`
(`floor(|x|·10^n + 0.5)`, es decir medio hacia arriba; ninguno de los dos valores cae
en un empate, así que acá da lo mismo que «medio al par»):

| Tejido | `prop_var` | ×100 | Redondeado | En el informe |
|---|---|---|---|---|
| `PLACENTA_E15` | 0.7632848468 | 76.328485 | 76 | 76 % |
| `BRAIN_E15` | 0.8155164853 | 81.551649 | 82 | 82 % |

El caso de cerebro es el que hace falta mirar: 81.55 **truncado** daría 81 y
**redondeado** da 82. El informe dice 82, consistente con el redondeo que
efectivamente aplica.

## Lo que esta conciliación no demuestra

Que el informe y las tablas dicen lo mismo, no que el cálculo que produjo las tablas
sea correcto: como observa la auditoría, «un error que escribió un CSV se replica en
el HTML». Eso lo cubren, hasta donde pueden, las verificaciones de tipo `recalculo` de
`outputs/tables/verificaciones.csv` y la concordancia independiente R↔Python de
`98_comparacion` y `99_verificar`. Tampoco reemplaza el punto 1 de «Evidencia mínima»
de la auditoría (entregar los artefactos con hashes por commit), que sigue pendiente
de la decisión de qué se publica.
