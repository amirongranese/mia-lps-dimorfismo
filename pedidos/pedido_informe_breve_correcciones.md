# PEDIDO — Informe breve, correcciones

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.
Todo esto es sobre `R/14_informe_breve.R` y sus dos versiones. **No se corre ningún
análisis nuevo**; los números siguen leyéndose de `outputs/tables/`.

---

## 1. Correcciones generales, en todo el documento

**1.1 Tildes.** El documento está escrito sin tildes: "activacion", "metodo",
"senalizacion", "expresion", "analisis", "deteccion". Corregir la acentuación en
todo el texto, incluidos títulos, epígrafes y notas al pie. Revisá que la
codificación del archivo no esté perdiendo los caracteres acentuados y la ñ.

**1.2 Formato de los p-valores.** Aparecen en notación científica cruda
(`p = 4.995005e-03`, `p = 1.000000e+00`, `p = 6.032341e-06`). Reformatearlos como se
reportan en un texto científico:

- Tres decimales: `p = 0,005`.
- Por debajo de 0,001: `p < 0,001`.
- Igual a 1: `p = 1,00`.
- Coma decimal, no punto, por ser texto en español.

Aplicar a todos los p del documento, incluidos los de Fisher, los del post hoc y los
de las listas de genes. Escribí una función de formato y usala en todos lados, en vez
de formatear en cada lugar.

**1.3 Códigos internos del repositorio.** El texto menciona `(D13)`, `(D7)` y
`ver docs/referencias.md`. Para un lector que no conoce el proyecto no significan
nada. Reemplazar por la explicación en palabras, y las referencias por citas
bibliográficas normales.

**1.4 Numerar las figuras.** Cada epígrafe empieza con `Figura 1.`, `Figura 2.`, y
así, y el texto las llama por su número donde corresponda.

---

## 2. Página 1

**2.1** Reemplazar el primer párrafo por:

> La inflamación materna durante la gestación es considerada un factor de riesgo para
> trastornos del neurodesarrollo en la descendencia, y la evidencia epidemiológica
> muestra que esas patologías afectan de manera distinta a varones y mujeres. En este
> contexto, la placenta constituye la interfaz entre la madre y el ambiente fetal, y
> su correcto funcionamiento resulta esencial para el desarrollo del feto gestante.
> Este proyecto emplea un modelo animal de activación inmune materna para
> caracterizar la respuesta placentaria y del cerebro fetal, en términos de la
> expresión de transportadores de nutrientes, frente a un estímulo inflamatorio
> materno inducido con lipopolisacárido (LPS), componente de la membrana externa de
> las bacterias gram negativas.

Las citas epidemiológicas van como referencias numeradas al pie o al final, tomadas
de `docs/referencias.md`, no como una remisión a un archivo del repositorio.

**2.2** Reemplazar el párrafo "Cómo se organizó el trabajo" por:

> Las decisiones metodológicas centrales se fijaron por escrito antes de ejecutar
> cualquier análisis y no se reabrieron una vez conocidos los resultados: la escala de
> análisis, el tratamiento de los valores no detectados, el modelo a ajustar y el
> criterio para realizar comparaciones post hoc. El trabajo se dividió en tareas
> acotadas, cada una con su propio cierre. El repositorio conserva el registro de
> procedencia de cada tabla y cada figura, y el de las verificaciones aplicadas a cada
> resultado. La implementación existe por duplicado, en R y en Python, y ambas se
> comparan celda a celda.

---

## 3. Página de placenta

Reemplazar el párrafo de pSTAT3 por:

> pSTAT3 aumentó en placentas de fetos hembra y no se modificó en machos. El
> dimorfismo placentario se manifiesta, entonces, en la señalización de esta vía a
> nivel proteico, y no en la expresión génica de los transportadores de nutrientes.

Los p del post hoc se incluyen con el formato de 1.2.

---

## 4. Página de cerebro fetal

**4.1** Reemplazar la frase final del párrafo de resultados por:

> En fatcd36, fatp4, gp130 y slc38a1 la interacción resultó significativa sin que
> ninguna de las comparaciones entre medias diera cuenta de ella. Esa falta de
> explicación por las medias motivó la exploración de los datos que se describe en la
> sección siguiente.

**4.2 Composición de las figuras.** Hoy la figura de proporción de detección de il6 y
los boxplots de cerebro están en la misma fila pero con tamaños muy distintos: la de
detección ocupa mucho más y los boxplots quedan ilegibles. Ponerlas **en una misma
fila y del mismo alto**, con los boxplots ocupando el ancho que necesiten para
leerse. Si con el mismo alto los boxplots siguen sin leerse, ponelas en dos filas y
avisame.

---

## 5. Página de exploración: figura nueva

Agregar `acto2_corr_placenta_cerebro_fatp4.png` en la subsección "Correlación entre
tejidos", en **tamaño moderado**, de modo que no empuje el documento más allá de las
5 páginas. **Si al agregarla el documento pasa de 5 páginas, no la incluyas** y
avisame.

**El epígrafe es crítico.** En esa figura, el panel de hembras muestra una
correlación aparente en el grupo control y ninguna en el grupo LPS. Un lector puede
concluir que el LPS rompió la correlación, que es justamente lo que el texto sostiene
que no puede afirmarse. El epígrafe tiene que decir, en una línea, que la diferencia
entre grupos se evaluó con un test formal y no resultó significativa, y que el
patrón es compatible con la reducción de variabilidad descripta en la sección.

---

## 6. Cierre

Generá las dos versiones y decime ruta, tamaño y número de páginas de cada una, y si
la figura nueva entró o no. Confirmame que ninguna sección de la versión pública
quedó vacía y que `git check-ignore` sigue excluyendo `outputs/informe_breve_real/`.
Un commit por bloque. Actualizá `ESTADO.md` y `procedencia.csv`.

Antes de escribir código decime en dos líneas qué vas a hacer.
