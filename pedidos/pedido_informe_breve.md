# PEDIDO — Informe breve del proyecto (máximo 5 páginas)

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

---

## 0. Qué es y qué no es

Un informe **de máximo 5 páginas** que cuente el proyecto: el problema, cómo se
trabajó, qué se encontró y qué límites tiene. Es un documento para leer, distinto de
tres cosas que ya existen y que **no se tocan**:

- `docs/informe.html`: el informe técnico largo, con procedencia y verificaciones.
- `docs/index.html` y `docs/presentacion.pdf`: la presentación.
- El repositorio y todos los análisis.

**No se corre ningún análisis nuevo.** Todos los números salen de
`outputs/tables/`, leídos al generar el documento, nunca escritos a mano. Si un
número que pido no está en ninguna tabla, avisame en vez de calcularlo.

**Script nuevo:** `R/14_informe_breve.R`, solo en R, misma excepción ya documentada
para `13_presentacion.R`. Reusá la maquinaria de render que ya existe.

**Dos versiones, igual que la presentación:**

| Versión | Datos | Destino | ¿Va a GitHub? |
|---|---|---|---|
| Pública | sintéticos | `docs/informe_breve.pdf` (y su HTML) | Sí |
| Real | reales | `outputs/informe_breve_real/` | No, ignorada por git |

Antes de escribir nada, verificá que `.gitignore` excluya
`outputs/informe_breve_real/`. Y acordate de la falla de la sesión pasada: **no
alcanza con bloquear el texto, las figuras también tienen que venir de la corrida
sintética.**

---

## 1. Estructura y reparto

### Página 1 — Problema, modelo experimental y método de trabajo

- Media página de contexto: inflamación materna y neurodesarrollo de la descendencia;
  por qué importa el sexo del feto; la placenta como interfaz. Podés apoyarte en el
  texto de la diapositiva 2 de la presentación y en `docs/referencias.md`.
- Media página de diseño experimental, con `assets/modelo-experimental.png`: LPS
  100 µg/kg i.p. en E15, colecta a las 6 h, 18 camadas, un feto de cada sexo por
  camada, 36 fetos, placenta y cerebro del mismo individuo, tres mediciones, cuatro
  grupos de n = 9.
- Un párrafo sobre cómo se organizó el trabajo: decisiones metodológicas fijadas por
  escrito antes de correr ningún análisis, trabajo dividido en tareas acotadas,
  repositorio con registro de procedencia y verificaciones, implementación duplicada
  en R y Python.

### Página 2 — Métodos

En prosa, sin figuras. Corto y preciso:

- Cuantificación relativa (D1) y por qué se analiza en escala log2 y se grafica como
  fold-change con eje logarítmico (D2).
- No detectados como NA, sin imputación, y por qué se descartó imputar (D3).
- La cascada de supuestos (D5) y el post hoc sobre las cuatro comparaciones de
  interés con corrección de Holm (D6).
- El caso de `il6` en cerebro, analizado solo como proporción de detección (D7).
- La censura a izquierda del ELISA de líquido amniótico (D10).
- Una línea sobre el modelo de pSTAT3 con membrana como bloque (D9).

### Página 3 — Placenta

- Una línea de validación del modelo: IL-6 detectable en suero materno, con el
  conteo y el p de Fisher leídos de la tabla.
- Los transportadores responden al LPS y lo hacen **igual en ambos sexos**: ningún
  gen con interacción; efecto principal de tratamiento en los genes que
  correspondan, con sus p. **Figura**: el subconjunto de boxplots de placenta que ya
  generaste para la presentación.
- pSTAT3 responde **solo en hembras**, con los valores y el p del post hoc.
  **Figura**: `acto1_pstat3_presentacion.png`.
- Cierra con la disociación: los transportadores placentarios cambian en ambos sexos
  aunque solo las hembras activan STAT3, así que esos cambios no dependen de esa vía
  o los machos llegan por otro camino. Presentalo como observación, no como
  conclusión demostrada.

### Página 4 — Cerebro fetal

El resultado principal. Dos formas de dimorfismo:

- **Desplazamiento del nivel en hembras**: los genes con post hoc significativo, con
  el patrón ♀Control ≠ ♀LPS y ♀LPS ≠ ♂LPS, sin cambios en machos. **Figura**: el
  subconjunto de boxplots de cerebro de la presentación.
- **Cambio de dispersión**: los genes con interacción significativa donde ninguna
  comparación de medias sobreviviera a Holm, y el test de interacción SEXO×TTO sobre
  la dispersión, con cuántos sobreviven a BH. El LPS compacta la expresión en hembras
  y la dispersa en machos. **Figura nueva**, ver sección 2.
- El argumento que une las dos partes: el boxplot mostraba que algo dependía del
  sexo pero no dónde estaba; el análisis de dispersión lo explicó.

### Página 5 — Eje entre tejidos, límites y conclusión

- Correlación placenta–cerebro: ninguna de las comparaciones fue significativa
  (conteo leído de la tabla). La caída de correlación observada en hembras se
  explica por la compactación de la expresión, según la simulación de restricción de
  rango. Sin figura: se cuenta en texto.
- Los diagramas triangulares, en dos o tres líneas, como exploración que sugiere
  estructura compartida entre transportadores, con el respaldo cuantitativo del PC1
  (porcentaje de varianza explicada en cada tejido, leído de la tabla). Decir
  explícitamente que las diferencias de rho entre grupos no se interpretan.
- Controles de sensibilidad en dos líneas: eigengene y exclusión del individuo
  extremo.
- **Limitaciones**, sin adornos: n = 9 por grupo; se asume independencia entre fetos
  aunque el tratamiento se aplica a la madre (D13); pSTAT3 normalizado a proteína
  total sin STAT3 total, así que refleja abundancia y no fracción fosforilada (D9);
  el efecto de dispersión es un patrón compartido por los transportadores, no
  hallazgos independientes por gen.
- **Conclusión**: el sexo del feto influye en la respuesta al LPS, pero de manera
  distinta en cada tejido. En placenta afecta la señalización y no el transporte; en
  cerebro afecta a los transportadores, desplazando el nivel en unos genes y la
  dispersión en otros.

---

## 2. Figura nueva: densidades de expresión por grupo

Reemplaza a `acto2_dispersion_sd.png` en este informe, porque muestra el mismo
resultado mejor: se ve la distribución completa en vez de un número resumen.

- Densidades de `-ΔΔCt` de **Control y LPS superpuestas**, con el mismo estilo que la
  diagonal de los diagramas triangulares (curvas suavizadas, rellenas, semitransparentes).
- **Cerebro E15**, dos filas: hembras arriba, machos abajo.
- **Tres o cuatro genes** en columnas, elegidos entre los que tienen interacción
  significativa sobre la dispersión. Elegilos leyendo `acto2_dispersion_interaccion.csv`
  y decime cuáles usaste.
- Lo que tiene que verse: en hembras el LPS es un pico angosto frente a un control
  ancho; en machos, al revés.
- No recalcula nada: usa los mismos `-ΔΔCt` que el resto del pipeline. Fila propia en
  `procedencia.csv`.

Si al generarla el contraste no se ve como lo describo, **no la fuerces**: avisame y
lo revisamos.

---

## 3. Reglas de formato

- Máximo 5 páginas en el PDF. Si se pasa, recortá texto, no figuras.
- Cinco figuras en total: modelo experimental, boxplots de placenta, pSTAT3, boxplots
  de cerebro, densidades. Cada una con su epígrafe breve.
- Prosa continua, sin viñetas salvo donde sea inevitable. Es un informe, no una
  presentación.
- Nada de lenguaje causal que los tests no respalden: "compatible con", "sugiere",
  "no se puede descartar", según corresponda.
- En la versión pública, las afirmaciones biológicas se reemplazan por el aviso de
  datos simulados, igual que en el informe técnico y la presentación.

---

## 4. Verificaciones

1. El PDF público tiene 5 páginas o menos.
2. La versión pública no contiene ningún número real ni frase interpretativa, y sus
   figuras provienen de la corrida sintética (comprobalo comparando estructura, no
   solo texto, como hiciste con pSTAT3).
3. `git check-ignore` excluye `outputs/informe_breve_real/`.
4. Todas las figuras referenciadas existen.

---

## 5. Cierre

Generá las dos versiones, decime ruta, tamaño y número de páginas de cada una, y qué
genes elegiste para la figura de densidades. Un commit por bloque. Actualizá
`ESTADO.md` y `procedencia.csv`.

Antes de escribir código decime en dos líneas qué vas a hacer, y avisame si alguna
sección no se puede armar con lo que hay en el repo.
