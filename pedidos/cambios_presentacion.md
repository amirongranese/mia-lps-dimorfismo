# CAMBIO — Página de presentación (HTML pública + PDF para exponer)

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

---

## 0. Qué hay que producir y por qué son dos cosas

La consigna del curso pide, dentro del repositorio, **una página HTML y un PDF**
desde los cuales se presenta. La página HTML va a quedar enlazada desde una página
pública con los proyectos de todo el curso.

Los resultados biológicos de este proyecto son **inéditos y no se publican**. Por eso
el mismo material se genera dos veces, cambiando solo la fuente de datos:

| Versión | Datos | Destino | ¿Va a GitHub? |
|---|---|---|---|
| Pública | sintéticos | `docs/index.html` y `docs/presentacion.pdf` | **Sí** |
| Para exponer | reales | `outputs/presentacion_real/` | **No**, ignorada por git |

Las dos salen del mismo código. La diferencia la da el mismo mecanismo que ya usa
`12_informe`: con datos sintéticos, ninguna afirmación biológica aparece; se muestra
el aviso de que los efectos son simulados. Con datos reales, la presentación
completa.

**Antes de escribir nada, verificá que `.gitignore` excluya
`outputs/presentacion_real/`.** Es el punto donde un error publica resultados
inéditos.

---

## 1. Script nuevo: `13_presentacion`

**Solo en R.** Es capa de presentación, no análisis: no hay resultado numérico nuevo
ni salida que comparar. Documentá esta excepción a la regla de scripts gemelos en
`AGENTS.md` y en `procedencia.csv`; no crees el gemelo en Python.

Usá la misma maquinaria de render que ya genera `docs/informe.html` y su PDF, sin
agregar dependencias nuevas. Si hace falta alguna, decímelo antes de instalarla.

Todas las figuras se toman de `outputs/figures/`; ninguna se regenera ni se dibuja
dentro del script. Los números que aparezcan en texto se leen de `outputs/tables/`,
igual que en el informe: nada escrito a mano.

---

## 2. Estructura: 10 minutos, 11 diapositivas

Es una exposición oral corta. Poco texto por diapositiva, una idea cada una, la
figura como protagonista.

**Contexto (3 diapositivas, ~2 min)**

1. **Título.** Título del proyecto, autora, curso, fecha.
2. **El problema.** Tres elementos en una sola diapositiva:
   - A la izquierda, `assets/ilustracion-mia.png` (embarazada). **Dibujá los rayos
     en SVG** al lado del vientre, con el color de acento de la presentación: no
     están en el PNG. Tres o cuatro trazos en zigzag, sobrios.
   - Al centro, una flecha y un recuadro con el texto **"Trastornos del
     neurodesarrollo"**, en HTML/SVG, no como imagen.
   - A la derecha, el gráfico de razones de prevalencia (sección 2bis).

   La idea que transmite: la inflamación materna durante la gestación es un factor
   de riesgo para trastornos del neurodesarrollo, y esas patologías afectan de forma
   distinta a varones y mujeres.
3. **El experimento.** `assets/modelo-experimental.png`, **sin modificar**, con el
   texto mínimo alrededor: LPS 100 µg/kg i.p. en E15, colecta a las 6 h; 18 camadas,
   un feto de cada sexo por camada, 36 fetos; placenta y cerebro del mismo individuo.
   Tres mediciones: IL-6 por ELISA, 10 genes por RT-qPCR, pSTAT3 por Western blot.
   Cuatro grupos, n = 9.

**Lo que produjo el agente (5 diapositivas, ~5 min)**

4. **El punto de partida fue un prompt.** Mostrar un fragmento **real** de las
   decisiones fijas de `AGENTS.md` (D1, D2 y D5 alcanzan) para que se vea el nivel de
   especificación que hizo falta: qué calibrador, en qué escala se analiza y en cuál
   se grafica, y la cascada según supuestos. La idea que transmite: no fue "analizá
   mis datos".
5. **Lo que salió de ahí.** Conteos leídos del repo, no escritos a mano: cantidad de
   scripts, de figuras, de tablas, de filas de `procedencia.csv` y de
   `verificaciones.csv` por tipo (`recalculo` / `existencia` / `declaracion`), y que
   R y Python producen CSV byte-idénticos.
6. **Placenta.** `acto1_pstat3.png`. La idea: la vía IL-6/STAT3 se activa solo en
   placentas de fetos hembra.
7. **Cerebro fetal.** El panel de expresión de cerebro, o un recorte con los genes
   con efecto específico de hembras. La idea: la respuesta depende del sexo del feto.
8. **Lo que agregó explorar con el agente.** La figura de dispersión agrupada al lado
   de la separada por sexo (`acto2_dispersion_sd.png` ya muestra los tres estratos).
   La idea: agrupando los sexos no se ve nada, porque los efectos opuestos se
   cancelan; separando aparece el cruce. **Esta es la diapositiva más importante de
   la presentación**: dedicale el espacio que haga falta.

**Errores y cómo se detectaron (3 diapositivas, ~3 min)**

9. **El agente toma decisiones y las documenta como propias del usuario.** Fijó los
   estratos del Acto 2 en GLOBAL/CONTROL/LPS, sin separar por sexo, y lo justificó en
   su propio archivo de descartes. Cómo se detectó: leyendo el informe y notando que
   la figura mostraba una cosa y el test medía otra. Los detalles están en
   `analisis_descartados.md`.
10. **Las verificaciones no verificaban.** Dos hechos concretos, tomados de
    `revisiones/AUDITORIA.md` y `revisiones/RESPUESTA.md`: las verificaciones pasaban
    en verde mientras el informe mostraba figuras viejas, porque comprobaban que el
    archivo existiera y no que estuviera embebido; y "todas las verificaciones
    pasaron" se imprimía aunque seis chequeos no se hubieran ejecutado. Cómo se
    detectó: con un agente externo, de otra empresa, sobre un clon limpio del repo.
    Qué se hizo: clasificar las verificaciones por tipo y agregar el estado
    `NO_EJECUTADA`.
11. **Cierre.** El agente hizo en días lo que llevaría semanas, y ninguna de sus
    verificaciones detectó sus propios errores: aparecieron al leer el informe y al
    auditarlo con otro agente. Qué haría falta para confiar en un análisis hecho así.

---

## 2bis. El gráfico de razones de prevalencia (diapositiva 2)

Hay un PNG en `assets/epidemiologia.png`, pero **no se usa**: se regenera para que
lleve el estilo de la presentación, la etiqueta de eje correcta y los nombres
completos. Es el único gráfico de la presentación que no sale de `outputs/figures/`,
porque no es un resultado del proyecto sino contexto bibliográfico. Anotalo así en
`procedencia.csv`.

**Qué representa:** razón de prevalencia entre sexos. La barra indica cuántas veces
más frecuente es la patología en el sexo más afectado respecto del otro. Etiqueta del
eje X: `Razón de prevalencia entre sexos`. Marca de referencia en 1 (sin diferencia).

**Datos** (de la revisión citada abajo; los rangos se muestran como rango, no como un
valor único):

| Patología | Razón | Sexo más afectado |
|---|---|---|
| Trastorno del espectro autista | 4:1 | varones |
| Enfermedad de Parkinson | 3,5:1 | varones |
| Trastorno por déficit de atención e hiperactividad | 3:1 | varones |
| Esclerosis lateral amiotrófica | 1,6:1 | varones |
| Esquizofrenia | 1,4:1 | varones |
| Esclerosis múltiple | 2–3:1 | mujeres |
| Enfermedad de Alzheimer | 1,6–3:1 | mujeres |
| Depresión y trastornos de ansiedad | 2:1 | mujeres |

**Cómo mostrarlo:** dos bloques, "Predominio en varones" y "Predominio en mujeres",
cada uno ordenado de mayor a menor razón. Nombres completos en el eje, sin siglas, de
modo que no haga falta ninguna leyenda de abreviaturas. Un color por bloque, tomados
de la paleta de la presentación.

**Fuente:** al pie del gráfico, una línea corta del tipo
`Razones de prevalencia compiladas de 10 estudios epidemiológicos (ver referencias.md)`.
Creá `docs/referencias.md` con las diez citas completas que están en
`pedidos/referencias_epidemiologia.md`. No inventes ninguna cita ni ningún número: si
algo no está en esa lista, no va.

---

## 3. Reglas de contenido

- **Nada escrito a mano que dependa de los datos.** Todo número o lista de genes se
  lee de las tablas, como en `12_informe`.
- **En la versión sintética**, las diapositivas 6, 7 y 8 muestran las figuras con el
  aviso de datos simulados, y el texto interpretativo se reemplaza por ese aviso.
  Las diapositivas 1 a 5 y 9 a 11 son iguales en las dos versiones: describen el
  proceso, no los resultados.
- Cada afirmación sobre errores del agente tiene que poder rastrearse a un archivo
  del repo (`analisis_descartados.md`, `revisiones/`). Si alguna no se sostiene con
  lo que hay escrito, no la pongas y avisame.
- Sin logos, sin texto de relleno. La única cita bibliográfica es la del gráfico de
  la diapositiva 2.

**Imágenes.** Las tres de `assets/` se versionan y son de autoría propia, así que van
en las dos versiones. `modelo-experimental.png` no se modifica ni se recorta.
`ilustracion-mia.png` se usa tal cual y los rayos se dibujan aparte, en SVG.
`epidemiologia.png` no se usa: queda en `assets/` como referencia de lo que había.

**Estilo.** Definí la paleta, la tipografía y el espaciado en un único bloque de
variables CSS al inicio, para poder cambiarlos en un solo lugar. Tomá como base la
paleta que ya usan las figuras del proyecto (`00_config.R`): azules para control,
y los colores por vía metabólica como acentos. Fondo claro, mucho aire, texto grande:
se proyecta en una sala. Nada de animaciones.

---

## 4. Verificaciones nuevas

1. La versión pública no contiene ninguna frase interpretativa (mismo chequeo que ya
   existe para el informe sintético).
2. `docs/index.html` y `docs/presentacion.pdf` existen, no están vacíos y se
   generaron con fuente de datos sintética.
3. `outputs/presentacion_real/` está cubierto por `.gitignore`: `git check-ignore` lo
   confirma.
4. Toda figura referenciada por la presentación existe en `outputs/figures/`.

---

## 5. Cierre

Generá las dos versiones y decime el tamaño y la ruta de cada una. Un commit para el
script y otro para las salidas públicas. Actualizá `ESTADO.md` y `procedencia.csv`.

Antes de escribir código decime en dos líneas qué vas a hacer, y avisame si alguna
diapositiva no se puede armar con lo que hay en el repo.
