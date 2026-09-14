# CAMBIO EN T8 — Dispersión y Δρ separados por sexo, + test de interacción sobre dispersión

Pegar como mensaje al agente. Aplica a `09_acto2_dispersion` y `10_acto2_simulacion`
(R y Python). Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

---

## 0. Por qué

Hoy todo el Acto 2.3–2.4 agrupa los sexos: los n de la tabla de Δρ son 16–18, o sea
hembras y machos juntos. Desde T7 las figuras de correlación están separadas por
sexo, así que el informe **muestra una cosa y testea otra**.

El caso que lo deja claro es `fatcd36`: ♀Control tiene rho = 0.86 (n=8) y ♂Control
tiene rho = −0.43 (n=8); agrupados dan 0.47, que no describe a ninguno de los dos.

El cambio es **aditivo**: los estratos agrupados se conservan tal como están. La
intención declarada es ver primero el comportamiento conjunto y después si al
separar por sexo aparecen efectos o tendencias que el promedio tapaba.

---

## 1. Δρ por sexo (sección 2.4)

Hoy la tabla compara `rho_control` vs `rho_lps` sobre los sexos juntos. Agregar las
dos comparaciones dentro de cada sexo:

- ♀Control vs ♀LPS
- ♂Control vs ♂LPS

Misma maquinaria que ya está: Fisher z sobre rho de Spearman, SE primario
Bonett-Wright, SE clásico como columna al lado, p a dos colas. **No cambiar el
método**, solo el subconjunto de pares.

La tabla gana una columna `ESTRATO` con los valores `AMBOS_SEXOS`, `HEMBRA`,
`MACHO`. Las filas actuales pasan a ser `AMBOS_SEXOS`.

**Potencia.** Con n ≤ 9 por celda, el Fisher z tiene muy poca potencia: un Δρ chico
es indetectable y uno grande puede no alcanzar significancia. Declaralo en el
reporte de la sección, no en una nota al pie. Sigue valiendo el piso de n ≥ 5 pares.

La columna BH suplementaria (D12) se calcula **dentro de cada estrato**, no a través
de los tres.

---

## 2. Dispersión por sexo (sección 2.3)

Hoy el Levene Brown-Forsythe compara Control vs LPS con los sexos juntos. Agregar el
mismo test **dentro de cada sexo**, por gen × tejido, con la misma columna `ESTRATO`.

Conservar las filas agrupadas actuales.

---

## 3. NUEVO — Test de interacción SEXO×TTO sobre la dispersión

Esto estaba en la especificación del Acto 2.3 y no se implementó. Es el análisis que
más falta hace.

Procedimiento, por gen × tejido, sobre los 36 fetos (no sobre pares):

1. Para cada feto, calcular `z = |x − mediana de su celda SEXO×TTO|`, donde `x` es el
   `-ΔΔCt`. Son las cuatro celdas de siempre, con la mediana de cada una.
2. Ajustar `z ~ SEXO * TTO` y reportar el ANOVA tipo III: F y p de `SEXO`, de `TTO` y
   de la **interacción**.
3. Tabla nueva `acto2_dispersion_interaccion.csv` con una fila por gen × tejido:
   gen, tejido, n por celda, mediana y desvío absoluto medio de cada celda, F y p de
   los tres términos, y columna BH suplementaria entre genes dentro de cada tejido.

Es la extensión factorial del test de Levene: el término de interacción pregunta si
el efecto del LPS **sobre la variabilidad** difiere entre sexos.

**Por qué importa acá.** En cerebro E15 el Acto 1 encontró interacción SEXO×TTO
significativa en siete genes. En tres de ellos el post hoc explica la interacción con
un patrón limpio y consistente (`glut1`, `slc38a2`, `fatp1`: sobreviven a Holm
♀Control vs ♀LPS y ♀LPS vs ♂LPS, o sea una respuesta específica de hembras).

En los otros cuatro —`fatcd36`, `fatp4`, `gp130`, `slc38a1`— la interacción es
significativa pero **ninguna comparación puntual sobrevive a Holm** (los mínimos van
de 0.074 a 0.267). Ese patrón es compatible con un efecto que está en la dispersión y
no en la media: el modelo detecta que algo difiere entre celdas, pero ningún
contraste de medias lo captura. Este test lo distingue, y si da significativo en
alguno de esos cuatro es un resultado del Acto 2 que el Acto 1 no podía ver.

Se aplica a todos los gen × tejido que se modelaron en T5, con el mismo piso de 5
detectados por celda.

---

## 4. La simulación tiene que cubrir los estratos nuevos

`10_acto2_simulacion` corre hoy sobre los Δρ agrupados. Cada Δρ nuevo por sexo
necesita su propio control de restricción de rango, o queda sin la verificación que
lo hace interpretable.

Extender la simulación a los estratos `HEMBRA` y `MACHO`, usando las SD observadas de
cada celda. Mismo diseño: una sola r verdadera común a los dos grupos, y lo único que
difiere entre ellos es la dispersión. Mismo veredicto DENTRO / FUERA.

Con n ≤ 9 el intervalo simulado va a ser ancho, lo cual hace que el veredicto DENTRO
sea casi automático. **Eso también es un resultado y hay que decirlo**: con este n no
se puede distinguir cambio de coordinación de cambio de dispersión.

---

## 5. Test de pendientes (confirmar antes de hacer)

El reporte actual dice que se eligió Fisher z "entre Fisher z / interacción de
pendientes / permutación". En la especificación original eran análisis distintos, no
alternativas: el Acto 2.5 pedía el test de pendientes (regresión entre pares de genes
con interacción X×TTO) **además** del test de correlaciones.

Antes de implementarlo, preguntá. Si se implementa: `Y ~ X * TTO` dentro de cada sexo
y tejido, reportando el p de la interacción, sobre los pares de transportadores.

---

## 6. Lo que NO cambia

- Prohibición 4: no se compara rho entre estratos sin test formal. Eso incluye los
  estratos nuevos: **no leer "significativo en hembras y no en machos" como
  dimorfismo**. Si esa comparación hace falta, es otro Fisher z, no una lectura.
- Prohibición 5: un cambio de rho no se interpreta como coordinación sin descartar
  restricción de rango.
- Emparejamiento por `FETO`, `il6R` excluido, solo Spearman, nada de Pearson.
- D1–D12 siguen fijas.

---

## 7. Verificaciones a agregar

1. Los estratos por sexo particionan al agrupado: para cada ítem y tratamiento,
   `n(HEMBRA) + n(MACHO) = n(AMBOS_SEXOS)`. Si no cierra, hay error de filtrado.
2. El test de interacción sobre dispersión reproduce el Levene Brown-Forsythe cuando
   se corre sobre un solo factor: verificalo con un gen × tejido y dejá la
   comparación registrada.
3. Toda fila nueva de las tablas tiene su `ESTRATO` poblado; ninguna queda en blanco.
4. Ninguna figura ni tabla nueva contiene Pearson.

---

## 8. Informe y procedencia

Filas nuevas en `procedencia.csv` para cada tabla y figura nueva. En la columna de
decisión, para los estratos por sexo: *"estratificación pedida explícitamente para
ver si el promedio entre sexos tapaba efectos; alternativa descartada: dejar solo el
estrato agrupado, descartada porque promedia correlaciones de signo opuesto (ver
fatcd36). Limitación declarada: n ≤ 9, potencia baja."*

Entrada nueva en `analisis_descartados.md` con el problema de agrupar sexos y por qué
se corrige.

En el informe, la sección 2.3–2.4 tiene que mostrar el estrato agrupado y los dos por
sexo en la misma tabla, y el texto tiene que decir explícitamente que agrupar
`fatcd36` promedia correlaciones de signo opuesto.

Regenerá el informe HTML y el PDF. Un commit por punto.

Antes de escribir código decime en dos líneas qué vas a hacer.
