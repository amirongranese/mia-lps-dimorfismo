# CAMBIO — Correcciones a partir de la revisión con un agente externo

Pegar como mensaje al agente. Leé `AGENTS.md` y `ESTADO.md` antes de tocar nada.

---

## 0. Contexto

El repositorio se revisó con un agente externo (Codex) que nunca vio el proyecto,
en dos pruebas sobre clones limpios de GitHub:

- **Prueba de reproducción** → `revisiones/PRUEBA_REPRODUCCION.md`
- **Auditoría escéptica** → `revisiones/AUDITORIA.md`

Los dos informes los voy a copiar en `revisiones/` (creá la carpeta si no existe).
Leelos completos antes de empezar: este pedido los resume, pero el detalle está ahí.

**Qué NO cambia:** D1–D12, los análisis, los datos, los resultados. No se agrega
ningún análisis estadístico nuevo. Todo lo de abajo es documentación, calidad de las
verificaciones, derivación del texto desde las tablas y reproducibilidad.

---

## 1. Decisión sobre la camada (MADRE) — documentar como D13

Codex señala (C4) que los modelos no tienen término de camada aunque el diseño tiene
18 madres con dos fetos cada una. **No se agrega MADRE al modelo.** La decisión ya
estaba tomada y es mía, pero no está escrita en ningún lado, y eso es lo que hay que
corregir.

**1.1** Agregar a `AGENTS.md` una decisión fija nueva:

> **D13 — Sin término de camada.** Los modelos no incluyen MADRE, ni como efecto
> fijo ni aleatorio. En análisis previos de este mismo modelo experimental se evaluó
> incluirla y no modificaba los resultados. Se asume independencia entre los fetos
> para el análisis; la limitación se declara en el informe.

**1.2** Entrada en `analisis_descartados.md`: se evaluó incluir MADRE como efecto
aleatorio; no se incorporó porque en análisis previos no modificaba los resultados.
Dejá un lugar marcado `<<< COMPLETAR: evidencia de análisis previos >>>` para que yo
agregue el dato concreto. No inventes ningún número.

**1.3** Agregar a Limitaciones del informe: el LPS se administra a la madre y cada
camada aporta un feto de cada sexo, así que los dos fetos de una camada no son
estrictamente independientes. El análisis asume independencia (D13). Si existiera
variación entre camadas, afectaría sobre todo a la precisión de los efectos
principales de tratamiento.

**1.4** Moderar dos frases de la conclusión revisada que la auditoría marca como más
fuertes que la inferencia:

- "cada individuo responde distinto" → "aumenta la variabilidad entre fetos".
- "no es un segundo hallazgo, sino la misma reducción de dispersión vista a través de
  la correlación" → "es compatible con la reducción de dispersión; la simulación
  muestra que esta alcanza para explicarla, aunque no permite descartar un cambio de
  coordinación".

Aplicar el mismo criterio a cualquier otra frase del informe que atribuya una causa
única donde los tests solo muestran compatibilidad.

---

## 2. Verificaciones: dejar de contar lo que no verifica

La auditoría encontró que buena parte de `verificaciones.csv` no puede fallar:
`TRUE` literales, chequeos de existencia de archivos, y constantes comparadas
consigo mismas. El informe dice "100/100 verificaciones", y eso sobrestima la
evidencia.

**2.1 Clasificar.** Agregar a `verificaciones.csv` una columna `tipo` con uno de
estos valores:

- `recalculo`: recalcula o contrasta un resultado de forma independiente.
- `existencia`: confirma que un archivo existe o no está vacío.
- `declaracion`: registra una decisión de método o una constante; no verifica su
  aplicación.

Usá la clasificación de `AUDITORIA.md` como punto de partida y revisá fila por fila.

**2.2 Informar por tipo.** El informe y la salida de `99_verificar` reportan cuántas
verificaciones hay de cada tipo y cuántas pasaron, por separado. La expresión
"100/100 verificaciones" sin desglose desaparece de todo el repo (informe, README,
ESTADO).

**2.3 Arreglar las que están mal**, no solo reclasificarlas:

- `informe_pdf`: hoy devuelve `TRUE` también cuando la generación del PDF falló. Un
  fallo tiene que dar `FALSE`. Si no hay motor de PDF, marcarlo como `NO_EJECUTADA`,
  no como aprobada.
- `pstat3_rama_cascada` (C2): hoy exige que la rama sea `anova3`. Tiene que
  comprobar que la rama elegida es la que corresponde a los p de Shapiro-Wilk y
  Levene registrados, cualquiera sea.
- En `05_qpcr_modelos`: reemplazar el chequeo de que las ramas sumen 18 por uno que
  compruebe, **gen por gen**, que la rama elegida corresponde a sus diagnósticos.
- `acto2_test_items`: tiene que exigir el número esperado de ítems, no solo
  reportarlo.

No hace falta reescribir todas las `declaracion` como recálculos: alcanza con que
estén rotuladas como lo que son.

---

## 3. El texto del informe se deriva de las tablas (C1, C3)

Hay afirmaciones interpretativas escritas como texto fijo en las conclusiones. En los
datos reales son ciertas, pero no se calculan: con otros datos podrían quedar falsas
sin que nadie lo note.

**3.1** "todas robustas a la corrección BH" (C1): calcularlo desde
`p_SEXOxTTO_BH` e informar numerador y denominador ("7 de 7 sobreviven a BH"). BH
sigue siendo suplementario (D12): el conteo se informa, no redefine qué es
significativo.

**3.2** "Ningún gen mostró interacción SEXO×TTO" en placenta (C3): derivarlo del
conteo real de interacciones en placenta.

**3.3** Revisá **todas** las secciones "Conclusión de la sección", la Síntesis y la
Conclusión revisada buscando el mismo patrón: toda afirmación que dependa de los
datos (qué genes, cuántos, en qué sexo, en qué dirección) se calcula de las tablas.
Si una frase no se puede derivar, se condiciona: solo aparece si los datos la
sostienen.

**3.4 El informe sabe con qué datos se generó.** Cuando la fuente es sintética, las
secciones interpretativas (conclusiones por sección, Síntesis, Conclusión revisada)
se reemplazan por un aviso: *"Informe generado con datos sintéticos. Los efectos son
simulados y arbitrarios; las conclusiones biológicas corresponden a los datos reales,
que no se incluyen en este repositorio."* Los números, tablas y figuras se muestran
igual. Agregar una verificación de tipo `recalculo` que compruebe que con datos
sintéticos no aparece ninguna frase interpretativa.

---

## 4. Reproducibilidad (de la prueba de reproducción)

**4.1 Python.** Averiguá dónde está el Python que usás en cada corrida: la prueba
externa no encontró `python`, `py` ni `python3` en el PATH de esta misma máquina.
Después documentá en el README cómo instalar Python (versión, fuente) y cómo lo
encuentra el pipeline. No dependas de una ruta que solo vos conocés.

**4.2 `run_all.ps1 -Only R`** no tiene que exigir Python. Chequeá cada intérprete
solo si se va a usar.

**4.3 `99_verificar` con una sola implementación.** Quien solo tenga R tiene que
poder verificar la parte R. La comparación R↔Python queda como `NO_EJECUTADA`,
visible en la salida, no como fallo ni como aprobada. La salida final no puede decir
"TODAS LAS VERIFICACIONES PASARON" si alguna quedó sin ejecutar.

**4.4 Ruta de R.** El README usa una ruta fija (`C:\Program Files\R\R-4.6.1\...`).
Que el pipeline busque `Rscript` en el PATH o en la instalación estándar, y que la
ruta fija quede solo como ejemplo. La prueba externa corrió en esta misma máquina,
así que no pudo detectar este problema; en otra computadora fallaría.

**4.5 Política de ejecución de PowerShell.** Documentar en el README cómo correr
`run_all.ps1` cuando la política lo bloquea (con `-ExecutionPolicy Bypass` solo para
ese proceso, sin cambiar la política del sistema).

**4.6** Mencionar en el README que `renv` necesita acceso de red a CRAN la primera
vez y usa su caché fuera del repositorio.

---

## 5. Conciliar las cinco afirmaciones

La auditoría rastreó cinco afirmaciones del informe hasta su tabla, pero no pudo
comprobar los números porque las tablas no están en el clon. Vos sí las tenés. Para
cada una, compará el valor del informe con la fila exacta del CSV (incluido el
redondeo) y dejá el resultado en `revisiones/conciliacion_auditoria.md`.

---

## 6. Respuesta a la revisión

Crear `revisiones/RESPUESTA.md` con una fila por cada hallazgo de los dos informes:
el hallazgo, qué se hizo (resuelto / declarado como limitación / no aplica) y en qué
commit. Los hallazgos que no se resuelvan también van, con el motivo.

No se versiona todavía nada de `docs/` ni de `outputs/`: qué se publica depende de
una decisión pendiente.

---

## 7. Orden y cierre

1. Punto 1 (D13 y lenguaje) y punto 3 (texto derivado): son los que cambian el
   informe.
2. Punto 2 (verificaciones).
3. Punto 4 (reproducibilidad).
4. Puntos 5 y 6.

Un commit por punto. Corré `run_all.ps1` completo al final y también con
`-Only R -FromSynthetic`, para comprobar que la corrida con una sola implementación y
datos sintéticos funciona y no muestra conclusiones biológicas. Actualizá `ESTADO.md`.

Esta es una tanda grande: si el contexto se te llena, cerrá la sesión al terminar el
punto 2 y dejá escrito en `ESTADO.md` dónde seguir.

Antes de escribir código decime en dos líneas qué vas a hacer.
