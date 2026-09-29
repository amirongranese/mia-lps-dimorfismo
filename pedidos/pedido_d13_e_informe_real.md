# PEDIDO — Evidencia de D13, informe breve con datos reales en docs/

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

**No se corre ningún análisis nuevo.** Los tres puntos son documentación y
configuración de salida.

---

## 1. Completar la evidencia de D13 (camada)

En `analisis_descartados.md` está el marcador
`<<< COMPLETAR: evidencia de análisis previos >>>` en la entrada de la decisión de no
incluir MADRE en el modelo. Reemplazalo por:

> En análisis previos sobre este mismo modelo experimental, el componente de varianza
> atribuido a la madre resultó igual a cero. Con varianza de camada nula, el modelo
> con madre como efecto aleatorio colapsa al modelo sin ella: las estimaciones y los
> errores estándar coinciden, de modo que incluirla no modifica ningún resultado.

Revisá que no quede ningún otro marcador `COMPLETAR` en archivos versionados, y
avisame si encontrás alguno.

---

## 2. El informe breve publicado pasa a ser el de datos reales

Cambio de criterio, decidido por la autora: el repositorio sigue funcionando
íntegramente con datos sintéticos (crudos excluidos, generador sintético, informe
técnico y presentación en sus versiones sintéticas), pero **el informe breve de
`docs/` pasa a ser el generado con los datos experimentales**.

**2.1** `docs/informe_breve.pdf` y `docs/informe_breve.html` se generan con datos
reales y se versionan.

**2.2 Nota en la primera página**, en tipografía menor, al pie o debajo del
subtítulo:

> Este informe fue generado a partir de los datos experimentales del proyecto, que
> por tratarse de resultados inéditos no se incluyen en el repositorio. El resto del
> repositorio funciona con datos sintéticos de idéntica estructura, de modo que el
> análisis completo puede reproducirse sin acceso a los datos originales; los valores
> obtenidos en esa reproducción no coinciden con los de este informe.

**2.3 Corregir la remisión final.** Hoy el informe termina con "Detalle completo,
procedencia y verificaciones: `docs/informe.html`". Como ese informe técnico es la
versión sintética, sus números no coinciden con los de este documento y un lector que
vaya a verificar se encuentra con otra cosa. Reformular esa línea para que remita a
la **metodología, el registro de procedencia y las verificaciones**, aclarando que
esos documentos corresponden a la corrida con datos sintéticos.

**2.4 Protección contra sobrescritura.** Hasta ahora `docs/` se llenaba con la
corrida sintética. Con este cambio, una corrida sintética completa podría pisar el
informe breve real sin que nadie lo note. Resolvelo de forma explícita: que el paso
que escribe `docs/informe_breve.*` solo se ejecute en la corrida con datos reales, y
que la corrida sintética no toque esos dos archivos. Dejalo documentado en
`AGENTS.md` y en el README, porque contradice la regla general de que `docs/` es
sintético.

**2.5** La versión sintética del informe breve se sigue generando, en
`outputs/informe_breve_sintetico/`, ignorada por git, para poder compararlas.

---

## 3. Verificaciones

1. `docs/informe_breve.pdf` se generó con datos reales y contiene la nota de 2.2.
2. Una corrida sintética completa **no modifica** `docs/informe_breve.pdf` ni su
   HTML: comprobalo por hash antes y después.
3. `docs/informe.html`, `docs/index.html` y `docs/presentacion.pdf` siguen siendo las
   versiones sintéticas, sin ningún número real.
4. `git status` no muestra nada bajo `data/`, `outputs/informe_breve_real/` ni
   `outputs/informe_breve_sintetico/`.
5. No queda ningún marcador `COMPLETAR` en archivos versionados.

---

## 4. Cierre

Antes de commitear, listame **exactamente qué archivos con datos reales quedarían
versionados** después de este cambio. Quiero revisarlo antes de que se suba.

Un commit por punto. Actualizá `ESTADO.md`, `AGENTS.md` y `procedencia.csv`.

Antes de escribir código decime en dos líneas qué vas a hacer.
