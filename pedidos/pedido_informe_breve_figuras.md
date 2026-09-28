# PEDIDO — Informe breve, ajustes de figuras y bibliografía

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

**Alcance acotado: solo lo que está en este pedido.** No modifiques el texto, los
números, la estructura ni ninguna otra figura. No se corre ningún análisis nuevo.

**Todas las figuras se generan desde el código**, como el resto. No incorpores
imágenes prearmadas: una imagen fija se usaría también en la versión pública y
expondría datos reales.

---

## 1. Figura 3 (cerebro): que se vea como la Figura 2

La Figura 2 (placenta) tiene cuatro paneles en una fila, todos del mismo tamaño y
proporción, y se lee bien. La Figura 3 tiene el gráfico de barras de il6 y los tres
boxplots con proporciones distintas, y los boxplots quedan estirados en vertical y
comprimidos en horizontal.

Componé la Figura 3 con **los cuatro paneles en una sola fila, del mismo ancho y el
mismo alto que los de la Figura 2**: el gráfico de detección de il6 y los boxplots de
glut1, slc38a2 y fatp1. La figura completa tiene que tener las mismas dimensiones y
la misma relación de aspecto que la Figura 2.

---

## 2. Figura 5: cambiar el recorte del diagrama triangular

Reemplazar la figura actual de densidades por un **recorte del diagrama triangular de
cerebro, hembras**, que conserve la estructura del diagrama y no solo la diagonal.
Es decir, un bloque de la matriz con:

- la **diagonal** con las densidades de Control y LPS superpuestas,
- el **triángulo inferior** con los puntos individuales,
- el **triángulo superior** con el valor de rho y el n.

Con los primeros cuatro o cinco genes alcanza; no la matriz completa. Tiene que verse
el título de cada gen arriba de su columna, como en el diagrama original.

El sentido del cambio: se ve la reducción de variabilidad en la diagonal **y**, al
mismo tiempo, el contexto de co-expresión del que habla el párrafo. El epígrafe
actual se mantiene, adaptado a que ahora se ve la estructura completa del recorte.

---

## 3. Figura 4 (correlación de fatp4): más grande, con la leyenda al costado

Hoy la leyenda de Spearman va debajo del gráfico, lo que obliga a achicar la figura
entera. Regenerá esa figura para el informe con la **leyenda a la derecha**, en una
columna: las cuatro líneas de rho, p y n, más la nota de los ejes.

Con la leyenda al costado, el área de los paneles puede crecer. Ampliá la figura lo
que permita el espacio disponible sin pasar de 5 páginas.

Esto es una variante de la figura solo para el informe, generada desde el mismo
código: no modifiques la figura original que usan el informe técnico y la
presentación. Fila propia en `procedencia.csv`.

---

## 4. Bibliografía al final

Las diez referencias están hoy al final de la primera página y ocupan mucho espacio.
Moverlas al **final del documento**, después de la conclusión, con **tipografía más
chica** que el cuerpo del texto. Las llamadas numeradas en el texto se mantienen.

---

## 5. Cierre

Generá las dos versiones y decime ruta, tamaño y número de páginas de cada una.
Confirmame que las figuras de la versión pública siguen saliendo de la corrida
sintética y que `git check-ignore` sigue excluyendo `outputs/informe_breve_real/`.
Un commit por bloque. Actualizá `ESTADO.md` y `procedencia.csv`.

Antes de escribir código decime en dos líneas qué vas a hacer.
