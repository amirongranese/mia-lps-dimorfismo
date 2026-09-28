# PEDIDO — Informe breve del proyecto (máximo 5 páginas)

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.
Este pedido reemplaza a cualquier versión anterior del informe breve.

---

## 0. Marco

Script `R/14_informe_breve.R`, solo en R (misma excepción ya documentada para
`13_presentacion.R`). Dos versiones:

| Versión | Datos | Destino | ¿Va a GitHub? |
|---|---|---|---|
| Pública | sintéticos | `docs/informe_breve.pdf` y su HTML | Sí |
| Real | reales | `outputs/informe_breve_real/` | No, ignorada por git |

No se corre ningún análisis nuevo. Todos los números se leen de `outputs/tables/` al
generar el documento; ninguno se escribe a mano. Si un número que pide el texto no
está en ninguna tabla, avisame en vez de calcularlo.

Las figuras de la versión pública tienen que provenir de una corrida con datos
sintéticos, no de `outputs/figures/` con datos reales.

**El PDF tiene que entrar en 5 páginas.** Si se pasa, recortá texto de las secciones
de métodos y de exploración, nunca figuras, y avisame qué sacaste. No achiques las
figuras para ganar espacio.

**En la versión pública ninguna sección queda vacía.** Se conservan la narrativa, los
métodos, la descripción de qué análisis se hicieron, las figuras y las limitaciones.
Solo se reemplaza por el aviso de datos simulados aquello que afirme **qué dio** un
análisis sobre los datos reales.

En todo el documento, nada de `--` literales: puntuación normal.

---

## 1. Página 1 — Problema, modelo experimental y método de trabajo

Sin cambios respecto de la versión actual. Está bien como está.

---

## 2. Página 2 — Métodos

Texto, sin figuras:

> Se entregaron al agente las tablas de datos crudos, sobre las que aplicó los
> cálculos estandarizados para obtener la cuantificación relativa de cada target. Las
> comparaciones entre los grupos experimentales de interés se definieron de antemano.
> Cuando el grupo calibrador no pudo cuantificarse por limitaciones del método
> experimental, ese gen no se analizó como expresión relativa sino como proporción de
> detección, mediante el test exacto de Fisher.
>
> La placenta y el cerebro analizados provienen del mismo individuo, de modo que las
> observaciones de ambos tejidos están pareadas por feto.
>
> El análisis se realizó sobre −ΔΔCt, en escala logarítmica de base 2, donde las
> diferencias son simétricas y aditivas; el fold-change se reservó para la
> representación gráfica, con eje logarítmico. Los valores no detectados se
> mantuvieron como faltantes y no se imputaron. Para cada gen y tejido se ajustó el
> modelo −ΔΔCt ~ sexo × tratamiento, eligiendo el método según el diagnóstico de los
> residuos: ANOVA de tipo III cuando se cumplían normalidad y homocedasticidad,
> errores robustos cuando fallaba la homocedasticidad, y una transformación por
> rangos alineados cuando fallaba la normalidad. Las comparaciones post hoc se
> realizaron únicamente cuando la interacción entre sexo y tratamiento resultó
> significativa, sobre cuatro contrastes definidos de antemano y con corrección de
> Holm.
>
> En el ELISA de líquido amniótico, los valores por debajo del límite de detección se
> trataron como censurados a izquierda y se analizaron con métodos específicos para
> datos censurados.

---

## 3. Página 3 — Placenta

> El agente devolvió el análisis por tejido y por gen, separando por sexo. En placenta
> se evaluó la expresión de diez genes: tres componentes de la vía de IL-6, una
> citoquina proinflamatoria fuertemente vinculada a patologías del neurodesarrollo, y
> siete transportadores de nutrientes esenciales para el desarrollo fetal,
> correspondientes al transporte de glucosa, aminoácidos y lípidos.

Completar a continuación, con los valores leídos de las tablas:

- **Validación del modelo**, en texto y sin figura: IL-6 detectable en suero materno,
  con el conteo por grupo y el p de Fisher. En líquido amniótico, que ningún
  contraste alcanzó significancia, con el n por grupo.
- **Resultado principal**: cuántos genes mostraron efecto de tratamiento, cuáles, y
  que ninguno mostró interacción entre sexo y tratamiento. **Verificá la dirección
  del cambio** de cada gen contra las medianas por grupo antes de escribirla. Si no
  todos se mueven en el mismo sentido, escribí "modificando la expresión" y detallá
  la dirección solo de los que puedas verificar.
- **pSTAT3, en una sola frase y sin figura**: aumenta en placentas de fetos hembra y
  no se modifica en machos, con el p del post hoc. Cerrar con la idea de que el
  dimorfismo placentario aparece en la señalización y no en el transporte.

**Figura:** boxplots de placenta de **il6, fatp1, slc38a2 y glut1**, en una sola
fila. Reusá la función de los subconjuntos cambiando la lista de genes.
**Sin figura de pSTAT3.**

---

## 4. Página 4 — Cerebro fetal

> El mismo análisis se aplicó al cerebro fetal. Mientras la placenta respondió al
> estímulo inflamatorio de manera equivalente en ambos sexos, en el cerebro fetal la
> respuesta de los transportadores de nutrientes resultó dimórfica.

Completar, verificado contra las tablas:

- Cuántos genes mostraron interacción entre sexo y tratamiento.
- En cuáles el post hoc ubica el efecto en las hembras, describiendo el patrón: el
  grupo control de hembras difiere del grupo LPS de hembras, y este último difiere
  del grupo LPS de machos, mientras que los machos no se modifican.
- Que en el resto la interacción fue significativa pero ninguna comparación entre
  medias la explicaba, lo que motivó el análisis de dispersión de la sección
  siguiente.

**Figuras:** la figura de proporción de detección de `il6` en cerebro (el gen no es
cuantificable por el método) y los boxplots de cerebro de **fatp1, slc38a2 y glut1**.
Si no entran las dos a igual tamaño, la de detección va más chica al costado.

---

## 5. Página 5 — Exploración, límites y conclusión

> Con el análisis convencional cerrado, el agente propuso profundizar por dos vías:
> correlacionar la expresión de cada gen entre la placenta y el cerebro de un mismo
> feto, y explorar la co-expresión entre genes dentro de cada tejido.

Después, en este orden:

**Correlación entre tejidos.** En las hembras, las correlaciones altas observadas en
el grupo control caían bajo LPS, lo que sugería una pérdida de acoplamiento entre
placenta y cerebro. El test formal no mostró diferencias significativas en ninguna de
las comparaciones (dar el conteo leído de la tabla), y una simulación mostró que una
reducción de la variabilidad alcanza por sí sola para producir esa caída, sin que la
relación subyacente cambie. La aparente pérdida de acoplamiento resultó entonces un
efecto de la menor dispersión bajo LPS.

**Co-expresión dentro de cada tejido.** Los diagramas triangulares son exploratorios y
describen cómo se acompañan los genes entre sí. **Sobre ellos no se testeó ninguna
diferencia entre grupos y no se interpretan diferencias de rho.** Lo que sí queda
establecido es que los transportadores varían mayormente juntos: el primer componente
principal explica una proporción alta de la varianza en ambos tejidos (dar los
porcentajes leídos de la tabla).

**Figura:** recorte del diagrama triangular de **cerebro, hembras**, mostrando la
diagonal con las densidades de Control y LPS superpuestas, con cuatro o cinco genes.
No la matriz completa, que no se lee impresa. En el epígrafe: la distribución del
grupo LPS es marcadamente más concentrada que la del control, y esa reducción de
variabilidad es el mecanismo que explica la caída de correlación.

**Limitaciones**, en dos o tres oraciones corridas, no como lista: n = 9 por grupo; se
asume independencia entre fetos aunque el tratamiento se administra a la madre (D13);
el efecto sobre la dispersión constituye un patrón compartido por los transportadores
y no hallazgos independientes por gen.

**Conclusión:** el sexo del feto influye en la respuesta al LPS, pero de manera
distinta en cada tejido. En la placenta, los transportadores responden de forma
equivalente en ambos sexos y el dimorfismo aparece en la señalización. En el cerebro
fetal, el dimorfismo se expresa en los transportadores, desplazando el nivel de
expresión en unos genes y la variabilidad entre individuos en otros.

---

## 6. Cierre

Cuatro figuras además de la imagen del modelo experimental de la página 1: boxplots
de placenta, detección de il6 en cerebro, boxplots de cerebro y recorte del diagrama
triangular. Cada una con epígrafe breve.

Prosa continua, sin viñetas. Nada de lenguaje causal que los tests no respalden:
"compatible con", "sugiere", "no permite descartar", según corresponda.

Generá las dos versiones y decime ruta, tamaño y número de páginas de cada una.
Confirmame que ninguna sección de la versión pública quedó vacía y que
`git check-ignore` sigue excluyendo `outputs/informe_breve_real/`. Un commit por
bloque. Actualizá `ESTADO.md` y `procedencia.csv`.

Antes de escribir código decime en dos líneas qué vas a hacer, y avisame si alguna
figura no se puede armar con lo que hay en el repo.
