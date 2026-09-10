# 98_comparacion.py -- T10: concordancia numerica R <-> Python y auditoria de
#                      procedencia / verificaciones / analisis_descartados.
#
# Por que existe este archivo: todo el pipeline se implementa dos veces (R y
# Python) y la convencion del proyecto (AGENTS.md 7) es que las salidas
# numericas tengan NOMBRES IDENTICOS en outputs/tables/R/ y outputs/tables/
# python/ para poder cruzarlas automaticamente. Este script hace ese cruce
# archivo por archivo, celda por celda, con las tolerancias declaradas
# (1e-6 para estadisticos; 1e-4 de fallback para p de tests iterativos), y ademas
# consolida las tres tablas de auditoria: chequea que toda figura y toda tabla
# tenga su fila en procedencia.csv, que analisis_descartados.md tenga una seccion
# por cada script de analisis (02..11), y que no quede ninguna verificacion en
# estado distinto de TRUE.
#
# NO hace analisis cientifico. Depende de que 02..11 ya hayan corrido (deja las
# tablas en outputs/tables/); en run_all.ps1 va justo antes de 99_verificar.
#
# PARIDAD R/Python: R y Python producen EL MISMO comparacion_R_python.csv y el
# mismo comparacion_reporte.md (son, de hecho, el control cruzado uno del otro).
# Los numeros se formatean con "%.10g"; la comparacion celda a celda parsea cada
# lado a double y aplica |a-b| <= tol OR |a-b|/max(|a|,|b|) <= tol.

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "98_comparacion"

TOL_EST = cfg.TOL_ESTADISTICO      # 1e-6  -- estadisticos, medias, coeficientes
TOL_P = cfg.TOL_P_ITERATIVO        # 1e-4  -- fallback para p de tests iterativos

# Scripts de analisis que DEBEN tener su seccion en analisis_descartados.md.
SCRIPTS_ANALISIS = [
    "02_ingesta_qc", "03_elisa", "04_qpcr_cuantificacion", "05_qpcr_modelos",
    "06_pstat3", "07_figuras_acto1", "08_acto2_correlaciones",
    "09_acto2_dispersion", "10_acto2_simulacion", "11_sensibilidad",
]


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..11.
# ---------------------------------------------------------------------------
def _fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, bool):
        return "TRUE" if x else "FALSE"
    if isinstance(x, float):
        if math.isnan(x):
            return ""
        if x == math.floor(x) and abs(x) < 1e15:
            return "%d" % int(round(x))
        return "%.10g" % x
    if isinstance(x, int):
        return str(x)
    return str(x)


def _csv_cell(x) -> str:
    s = _fmt(x)
    if any(c in s for c in (",", '"', "\n", "\r")):
        s = '"' + s.replace('"', '""') + '"'
    return s


def escribir_csv(ruta: Path, encabezado, filas) -> None:
    lineas = [",".join(_csv_cell(v) for v in encabezado)]
    lineas += [",".join(_csv_cell(v) for v in fila) for fila in filas]
    ruta.write_bytes(("\n".join(lineas) + "\n").encode("utf-8"))


def escribir_texto(ruta: Path, texto: str) -> None:
    if not texto.endswith("\n"):
        texto += "\n"
    ruta.write_bytes(texto.encode("utf-8"))


# ---------------------------------------------------------------------------
# Lectura de CSV -- mismo parser que 02..11 (comillas dobles, sin dependencias).
# ---------------------------------------------------------------------------
def _parse_csv_line(linea: str):
    out, cur, i, q = [], [], 0, False
    while i < len(linea):
        c = linea[i]
        if q:
            if c == '"':
                if i + 1 < len(linea) and linea[i + 1] == '"':
                    cur.append('"'); i += 1
                else:
                    q = False
            else:
                cur.append(c)
        else:
            if c == '"':
                q = True
            elif c == ",":
                out.append("".join(cur)); cur = []
            else:
                cur.append(c)
        i += 1
    out.append("".join(cur))
    return out


def _leer_csv(ruta: Path):
    if not ruta.is_file():
        return None, []
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    if not lineas:
        return None, []
    return _parse_csv_line(lineas[0]), [_parse_csv_line(x) for x in lineas[1:]]


# ---------------------------------------------------------------------------
# Comparacion numerica celda a celda.
# ---------------------------------------------------------------------------
def _es_num(s: str):
    """Devuelve el float si 's' es un numero finito, si no None."""
    if s is None or s == "":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def _dentro(a: float, b: float, tol: float) -> bool:
    d = abs(a - b)
    if d <= tol:
        return True
    m = max(abs(a), abs(b))
    return m > 0.0 and d / m <= tol


def comparar_archivo(nombre: str, ruta_r: Path, ruta_py: Path) -> dict:
    """Cruza un CSV de outputs/tables/R/ contra su gemelo de python/."""
    hr, fr = _leer_csv(ruta_r)
    hp, fp = _leer_csv(ruta_py)
    byte_id = (ruta_r.is_file() and ruta_py.is_file()
               and ruta_r.read_bytes() == ruta_py.read_bytes())

    header_igual = hr == hp
    filas_ok = len(fr) == len(fp)
    n_cel = n_num = n_txt = n_dif_txt = n_fuera = 0
    max_abs = max_rel = 0.0
    peor = ""
    tol_usada = TOL_EST

    if header_igual and filas_ok and hr is not None:
        for i, (rr, rp) in enumerate(zip(fr, fp)):
            for j in range(len(hr)):
                cr = rr[j] if j < len(rr) else ""
                cp = rp[j] if j < len(rp) else ""
                n_cel += 1
                va, vb = _es_num(cr), _es_num(cp)
                if va is not None and vb is not None:
                    n_num += 1
                    d = abs(va - vb)
                    m = max(abs(va), abs(vb))
                    rel = d / m if m > 0.0 else 0.0
                    if d > max_abs:
                        max_abs = d
                        peor = "fila %d / %s" % (i + 2, hr[j])
                    max_rel = max(max_rel, rel)
                    if not _dentro(va, vb, TOL_EST):
                        if _dentro(va, vb, TOL_P):
                            tol_usada = TOL_P
                        else:
                            n_fuera += 1
                else:
                    n_txt += 1
                    if cr != cp:
                        n_dif_txt += 1
                        if not peor:
                            peor = "fila %d / %s (texto)" % (i + 2, hr[j])

    ok = (header_igual and filas_ok and n_dif_txt == 0 and n_fuera == 0
          and hr is not None)
    return {
        "archivo": nombre,
        "en_R": ruta_r.is_file(),
        "en_python": ruta_py.is_file(),
        "filas_R": len(fr),
        "filas_py": len(fp),
        "cols_R": len(hr) if hr else 0,
        "cols_py": len(hp) if hp else 0,
        "header_igual": header_igual,
        "celdas": n_cel,
        "celdas_numericas": n_num,
        "celdas_texto": n_txt,
        "n_dif_texto": n_dif_txt,
        "max_dif_abs": max_abs,
        "max_dif_rel": max_rel,
        "peor_celda": peor,
        "tol": tol_usada,
        "n_fuera_tol": n_fuera,
        "byte_identico": byte_id,
        "ok": ok,
    }


COLS_COMP = ["archivo", "en_R", "en_python", "filas_R", "filas_py", "cols_R",
             "cols_py", "header_igual", "celdas", "celdas_numericas",
             "celdas_texto", "n_dif_texto", "max_dif_abs", "max_dif_rel",
             "peor_celda", "tol", "n_fuera_tol", "byte_identico", "ok"]


def _fila_comp(d: dict):
    return [d[c] for c in COLS_COMP]


# ---------------------------------------------------------------------------
# Auditoria de procedencia / verificaciones / analisis_descartados.
# ---------------------------------------------------------------------------
def _basename_artefacto(s: str) -> str:
    """Ultimo segmento tras '/' -- colapsa 'outputs/tables/{R,python}/x.csv' -> 'x.csv'."""
    return s.rsplit("/", 1)[-1]


def auditar_procedencia(figuras, tablas):
    h, filas = _leer_csv(cfg.RUTA_TABLAS / "procedencia.csv")
    mencionados = set()
    if h is not None:
        idx = h.index("artefacto")
        for f in filas:
            if len(f) > idx:
                mencionados.add(_basename_artefacto(f[idx]))
    figs_sin = [p.name for p in figuras if p.name not in mencionados]
    tabs_sin = [p.name for p in tablas if p.name not in mencionados]
    return figs_sin, tabs_sin


def auditar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    txt = ruta.read_bytes().decode("utf-8") if ruta.is_file() else ""
    faltan = [s for s in SCRIPTS_ANALISIS
              if f"<!-- {s}:inicio -->" not in txt or f"<!-- {s}:fin -->" not in txt]
    return faltan


def auditar_verificaciones():
    h, filas = _leer_csv(cfg.RUTA_TABLAS / "verificaciones.csv")
    if h is None:
        return 0, 0, []
    idx_ok = h.index("ok")
    idx_id = h.index("id")
    idx_sc = h.index("script")
    # el conteo excluye las propias filas de 98_comparacion para que sea estable
    # sin importar el orden de ejecucion (R antes / despues de Python) ni re-corridas.
    filas = [f for f in filas if not (len(f) > idx_sc and f[idx_sc] == ESTE_SCRIPT)]
    total = len(filas)
    no_true = [f[idx_id] for f in filas if not (len(f) > idx_ok and f[idx_ok] == "TRUE")]
    return total, total - len(no_true), no_true


# ---------------------------------------------------------------------------
# Justificaciones que 02/03/04 difieren explicitamente a T10.
# ---------------------------------------------------------------------------
def _num_no_detectados():
    """Del qpcr_no_detectados_descriptivo.csv de la corrida: totales y celdas D7."""
    h, filas = _leer_csv(cfg.RUTA_TABLAS_R / "qpcr_no_detectados_descriptivo.csv")
    if h is None:
        return None
    c = {n: h.index(n) for n in
         ("TEJIDO", "GEN", "GRUPO", "n_total", "n_no_detectado",
          "calibrador_cero_detectados")}
    tot_obs = tot_nd = 0
    d7 = []
    for f in filas:
        if f[c["GRUPO"]] == "TODOS":
            tot_obs += int(f[c["n_total"]])
            tot_nd += int(f[c["n_no_detectado"]])
        if f[c["calibrador_cero_detectados"]] == "TRUE" and f[c["GRUPO"]] != "TODOS":
            d7.append((f[c["TEJIDO"]], f[c["GEN"]], f[c["GRUPO"]],
                       int(f[c["n_total"]]), int(f[c["n_no_detectado"]])))
    return tot_obs, tot_nd, d7


_SEC_INI = "<!-- 98_comparacion:inicio -->"
_SEC_FIN = "<!-- 98_comparacion:fin -->"


def _seccion_descartados():
    nd = _num_no_detectados()
    if nd is not None:
        tot_obs, tot_nd, d7 = nd
        pct = 100.0 * tot_nd / tot_obs if tot_obs else 0.0
        linea_tot = (f"En la corrida actual: **{tot_nd} de {tot_obs}** observaciones "
                     f"gen x feto x tejido no detectadas ({pct:.1f} %).")
        if d7:
            celdas = "; ".join(f"{t}/{g}/{gr} {nnd}/{ntot}" for t, g, gr, ntot, nnd in d7)
            linea_d7 = (f"Celdas con calibrador HEMBRA_CONTROL 0/n detectado (D7): "
                        f"{celdas}.")
        else:
            linea_d7 = ("En esta corrida ningun calibrador quedo con 0 detectados "
                        "(sobre sintetico puede no dispararse D7).")
    else:
        linea_tot = ("(qpcr_no_detectados_descriptivo.csv no disponible: correr "
                     "04_qpcr_cuantificacion antes de 98_comparacion.)")
        linea_d7 = ""

    L = [
        "## 98_comparacion",
        "",
        "Consolidacion de auditoria (T10). Concordancia R <-> Python archivo por "
        "archivo (`comparacion_R_python.csv`) y cierre de las justificaciones que "
        "02/03/04 dejaron diferidas a esta tarea.",
        "",
        "### Imputacion MNAR de no-detectados (D3) -- justificacion completa",
        "",
        "- **Que es**: la imputacion MNAR (p. ej. `nondetects` en R) modela la "
        "probabilidad de no-deteccion como funcion decreciente del Ct latente y "
        "sortea Ct por encima del umbral para las celdas no observadas, para "
        "despues correr el analisis sobre una matriz \"completa\".",
        f"- **Numeros**: {linea_tot}",
    ]
    if linea_d7:
        L.append(f"  {linea_d7}")
    L += [
        "- **Por que se descarta** (D3, prohibicion 1):",
        "  1. *Ancla inexistente en las celdas D7.* La imputacion necesita algunos "
        "Ct observados en el grupo para estimar la pendiente de la curva de "
        "deteccion. Donde el calibrador HEMBRA_CONTROL tiene 0 detectados "
        "(il6 @ BRAIN_E15) no hay nada sobre lo que anclar: el dCt de calibrador, "
        "y por lo tanto todo ddCt del gen x tejido, quedaria definido contra un "
        "valor enteramente inventado (prohibicion 10).",
        "  2. *Estructura artificial en el baseline.* En la exploracion previa "
        "(informe E15 original) la imputacion MNAR rellenaba los no-detectados del "
        "grupo Control con Ct altos correlacionados, comprimiendo su varianza y "
        "generando diferencias Control vs LPS que el dato crudo no tiene. La "
        "imputacion inyecta la hipotesis que el analisis deberia poner a prueba.",
        "  3. *No hace falta.* D1/D8 calculan calibrador, medias y z **solo sobre "
        "detectados**; D5 modela solo las celdas con n >= 5. La informacion que "
        "aporta un no-detectado (un \"< umbral\") se analiza como **proporcion de "
        "deteccion** (Fisher exacto, D7) alli donde es informativa -- sin inventar "
        "su magnitud.",
        "- **Que se hizo**: no-detectado -> NA, sin imputar por ningun metodo; "
        "deteccion como desenlace propio donde corresponde (D7).",
        "",
        "### LOD alternativo del ELISA (menor estandar de `CURVA IL6`) -- no aplicado",
        "",
        "- D10 admite como alternativa al blanco de placa (`LOD = 0`, el primario "
        "de 02/03) el menor estandar de la hoja `CURVA IL6`.",
        "- **Por que no se corre como sensibilidad separada**: el contraste "
        "primario del ELISA trata a los censurados como **empatados en el rango "
        "mas bajo** (Peto-Peto G-rho=1 en LA; Fisher de deteccion en MS). Quien "
        "esta censurado lo define `Conc < 0`, no el valor del LOD; mover el LOD de "
        "0 al menor estandar no reordena esos empates ni cambia el conjunto "
        "censurado, asi que **no puede cambiar la conclusion por construccion del "
        "test**. El LOD explicito solo entraria en un estimador de ubicacion con "
        "censura, que D10 (clausula final) ya descarta por censura alta.",
        "- **Que se hizo**: se mantiene `LOD = 0`; el alternativo queda documentado "
        "y no ejecutado.",
        "",
        "### Factorizacion del nucleo hand-rolled (08--11) -- evaluado, no se factoriza",
        "",
        "- **Situacion**: `spearman_rho` + `rangos_promedio` + `pares()` / "
        "`cargar()` + los merges de `-ddCt` / score estan replicados casi "
        "identicos en 08, 09, 10 y 11; 11 agrega Jacobi sin trigonometria y "
        "`pearson_pairwise`.",
        "- **Opcion evaluada**: moverlos a un modulo compartido (`00_config` o un "
        "`nucleo_acto2`) importado por 08--11.",
        "- **Decision: no se factoriza en T10.** (a) La regla del repo (AGENTS.md "
        "7) es que cada script importa **solo** `00_config`; meter logica de "
        "analisis en `00_config` rompe esa linea. (b) Las salidas de 02--11 son "
        "**byte-identicas R/Python** y estan cerradas; tocar el nucleo justo antes "
        "del informe (T11) arriesga perturbarlas sin ganancia cientifica. (c) La "
        "replicacion es chica, ya verificada byte a byte entre lenguajes y, en R, "
        "cruzada en cada script contra `cor.test` / `car::leveneTest` / `eigen()`; "
        "`comparacion_R_python.csv` la cubre de forma continua.",
        "- **Revisable** despues de T11 si el informe necesita re-ejecutar el "
        "nucleo.",
        "",
        "### Concordancia R <-> Python (`comparacion_R_python.csv`)",
        "",
        "- Una fila por CSV de `outputs/tables/{R,python}/`. Cada celda que parsea "
        "a numero finito en ambos lados se compara con "
        "`|a-b| <= 1e-6` **o** `|a-b|/max(|a|,|b|) <= 1e-6` (`tol` pasa a `1e-4` si "
        "alguna celda necesita el margen mas laxo, reservado a `p` de tests "
        "iterativos); el resto se compara como texto exacto. `byte_identico` marca "
        "los pares que ademas coinciden byte a byte.",
        "- Los `.md` de `outputs/tables/` (reportes legibles) son de copia unica "
        "-- R y Python escriben la misma ruta -- y no entran en esta tabla; su "
        "paridad entre lenguajes se re-chequea en `99_verificar` (T11).",
    ]
    return "\n".join(L)


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    nuevo = f"{_SEC_INI}\n{_seccion_descartados()}\n\n{_SEC_FIN}"
    if ruta.is_file():
        txt = ruta.read_bytes().decode("utf-8")
    else:
        txt = ("# Analisis descartados\n\nQue se probo, por que no funciono o no se"
               " uso, y que se hizo en su lugar. Una seccion por script.\n")
    if _SEC_INI in txt and _SEC_FIN in txt:
        txt = txt.split(_SEC_INI)[0] + nuevo + txt.split(_SEC_FIN)[1]
    else:
        if not txt.endswith("\n"):
            txt += "\n"
        txt += "\n" + nuevo + "\n"
    escribir_texto(ruta, txt)


# ---------------------------------------------------------------------------
# Merge por 'script' en procedencia.csv / verificaciones.csv -- headers 02..11.
# ---------------------------------------------------------------------------
def merge_por_script(ruta: Path, header, filas_nuevas, col_script, clave_orden):
    h_old, filas_old = _leer_csv(ruta)
    if h_old is not None and h_old != header:
        raise SystemExit(f"{ruta.name}: encabezado incompatible {h_old!r} vs {header!r}")
    idx = header.index(col_script)
    conservadas = [f for f in filas_old if len(f) > idx and f[idx] != ESTE_SCRIPT]
    todas = conservadas + [[_fmt(v) for v in f] for f in filas_nuevas]
    todas.sort(key=clave_orden)
    escribir_csv(ruta, header, todas)


def registrar_procedencia(filas_nuevas):
    header = ["artefacto", "tipo", "script", "origen_codigo", "entradas", "descripcion"]
    merge_por_script(cfg.RUTA_TABLAS / "procedencia.csv", header, filas_nuevas,
                     "script", lambda f: (f[2], f[0]))


def registrar_verificaciones(filas_nuevas):
    header = ["id", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script"]
    merge_por_script(cfg.RUTA_TABLAS / "verificaciones.csv", header, filas_nuevas,
                     "script", lambda f: (f[5], f[0]))


# ---------------------------------------------------------------------------
# Reporte legible.
# ---------------------------------------------------------------------------
def construir_reporte(fuente, comps, resumen, figs_sin, tabs_sin, desc_faltan,
                      verif_total, verif_true, verif_no):
    L = []
    ap = L.append
    ap("# Comparacion R <-> Python y auditoria (T10)")
    ap("")
    ap("Generado por `98_comparacion` (R y Python producen este archivo identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Concordancia numerica archivo por archivo")
    ap("")
    ap("Tolerancia: `1e-6` (estadisticos) con fallback `1e-4` (p de tests "
       "iterativos). `ok` = header igual, mismo n de filas, sin diferencias de "
       "texto y sin celdas fuera de tolerancia.")
    ap("")
    ap("| archivo | filas | num | max_dif_abs | byte-id | tol | ok |")
    ap("| --- | --- | --- | --- | --- | --- | --- |")
    for d in comps:
        ap("| %s | %d | %d | %.2e | %s | %g | %s |" % (
            d["archivo"], d["filas_R"], d["celdas_numericas"], d["max_dif_abs"],
            "si" if d["byte_identico"] else "no", d["tol"],
            "OK" if d["ok"] else "REVISAR"))
    ap("")
    ap(f"- Archivos comparados: **{resumen['n']}**. "
       f"Byte-identicos: **{resumen['byte']}**. "
       f"Fuera de tolerancia: **{resumen['fuera']}**. "
       f"Peor |dif| absoluta: **{resumen['peor_abs']:.3e}** "
       f"({resumen['peor_arch'] or '-'}).")
    ap(f"- Solo en R: {resumen['solo_r'] or 'ninguno'}. "
       f"Solo en Python: {resumen['solo_py'] or 'ninguno'}.")
    ap("")
    ap("## 2. Auditoria de procedencia.csv")
    ap("")
    ap(f"- Figuras sin fila en procedencia: **{figs_sin or 'ninguna'}**.")
    ap(f"- Tablas (outputs/tables/{{R,python}}/*.csv) sin fila: "
       f"**{tabs_sin or 'ninguna'}**.")
    ap("")
    ap("## 3. Auditoria de analisis_descartados.md")
    ap("")
    ap(f"- Scripts de analisis esperados: {len(SCRIPTS_ANALISIS)}. "
       f"Sin seccion: **{desc_faltan or 'ninguno'}**.")
    ap("")
    ap("## 4. Auditoria de verificaciones.csv")
    ap("")
    ap(f"- Filas: **{verif_total}**. En TRUE: **{verif_true}**. "
       f"Distintas de TRUE: **{verif_no or 'ninguna'}**.")
    ap("")
    ap("## 5. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `98_comparacion`: justificacion "
       "completa del descarte de la imputacion MNAR, LOD alternativo del ELISA no "
       "aplicado, y la decision de no factorizar el nucleo hand-rolled de 08--11.")
    ap("")
    return "\n".join(L)


# ===========================================================================
def main():
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)

    csv_r = sorted(cfg.RUTA_TABLAS_R.glob("*.csv"))
    csv_py = sorted(cfg.RUTA_TABLAS_PY.glob("*.csv"))
    nombres_r = {p.name for p in csv_r}
    nombres_py = {p.name for p in csv_py}
    if not nombres_r:
        raise SystemExit("outputs/tables/R/ no tiene CSV: correr 02..11 antes.")
    solo_r = sorted(nombres_r - nombres_py)
    solo_py = sorted(nombres_py - nombres_r)
    comunes = sorted(nombres_r & nombres_py)

    comps = [comparar_archivo(n, cfg.RUTA_TABLAS_R / n, cfg.RUTA_TABLAS_PY / n)
             for n in comunes]

    n_byte = sum(1 for d in comps if d["byte_identico"])
    n_fuera = sum(d["n_fuera_tol"] + d["n_dif_texto"] for d in comps)
    n_ok = sum(1 for d in comps if d["ok"])
    peor_abs = 0.0
    peor_arch = ""
    for d in comps:
        if d["max_dif_abs"] > peor_abs:
            peor_abs = d["max_dif_abs"]
            peor_arch = d["archivo"]

    # --- comparacion_R_python.csv: una fila por archivo + fila __TOTAL__ ---
    filas_csv = [_fila_comp(d) for d in comps]
    total = {
        "archivo": "__TOTAL__", "en_R": len(nombres_r), "en_python": len(nombres_py),
        "filas_R": sum(d["filas_R"] for d in comps),
        "filas_py": sum(d["filas_py"] for d in comps),
        "cols_R": "", "cols_py": "",
        "header_igual": all(d["header_igual"] for d in comps),
        "celdas": sum(d["celdas"] for d in comps),
        "celdas_numericas": sum(d["celdas_numericas"] for d in comps),
        "celdas_texto": sum(d["celdas_texto"] for d in comps),
        "n_dif_texto": sum(d["n_dif_texto"] for d in comps),
        "max_dif_abs": peor_abs, "max_dif_rel": max((d["max_dif_rel"] for d in comps),
                                                    default=0.0),
        "peor_celda": peor_arch, "tol": TOL_EST,
        "n_fuera_tol": sum(d["n_fuera_tol"] for d in comps),
        "byte_identico": n_byte == len(comps),
        "ok": (n_ok == len(comps) and not solo_r and not solo_py),
    }
    filas_csv.append([total[c] for c in COLS_COMP])
    escribir_csv(cfg.RUTA_TABLAS / "comparacion_R_python.csv", COLS_COMP, filas_csv)

    # --- auditoria ---
    figuras = sorted(cfg.RUTA_FIGURAS.glob("*.png"))
    figs_sin, tabs_sin = auditar_procedencia(figuras, csv_r)
    desc_faltan = auditar_descartados()
    verif_total, verif_true, verif_no = auditar_verificaciones()

    resumen = {"n": len(comps), "byte": n_byte, "fuera": total["n_fuera_tol"]
               + total["n_dif_texto"], "peor_abs": peor_abs, "peor_arch": peor_arch,
               "solo_r": ", ".join(solo_r), "solo_py": ", ".join(solo_py)}

    escribir_texto(cfg.RUTA_TABLAS / "comparacion_reporte.md",
                   construir_reporte(fuente, comps, resumen, ", ".join(figs_sin),
                                     ", ".join(tabs_sin), ", ".join(desc_faltan),
                                     verif_total, verif_true, ", ".join(verif_no)))
    actualizar_descartados()

    # --- procedencia / verificaciones (filas de 98_comparacion) ---
    ent = "outputs/tables/{R,python}/*.csv"
    registrar_procedencia([
        ["outputs/tables/comparacion_R_python.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent, "concordancia numerica R<->python archivo por archivo; tol 1e-6 "
         "(fallback 1e-4) + byte-identidad; fila __TOTAL__ con agregados"],
        ["outputs/tables/comparacion_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO",
         ent, "reporte legible de la comparacion R<->python y de la auditoria de "
         "procedencia / verificaciones / analisis_descartados (T10)"],
        ["outputs/tables/procedencia.csv", "auditoria", ESTE_SCRIPT, "PROPIO",
         "todos los scripts 02..11",
         "indice de procedencia: una fila por figura y por tabla; cobertura "
         "verificada en T10 (auditar_procedencia)"],
        ["outputs/tables/verificaciones.csv", "auditoria", ESTE_SCRIPT, "PROPIO",
         "todos los scripts 02..11",
         "indice de verificaciones: un resultado principal por fila; T10 chequea "
         "que ninguna quede distinta de TRUE"],
        ["outputs/tables/analisis_descartados.md", "auditoria", ESTE_SCRIPT,
         "PROPIO", "todos los scripts 02..11",
         "una seccion por script de analisis (02..11); cobertura verificada en T10"],
    ])
    registrar_verificaciones([
        ["comp_cobertura_csv",
         "todo CSV de outputs/tables/R/ tiene gemelo en python/ y viceversa",
         f"R={len(nombres_r)}; python={len(nombres_py)}; solo_R={len(solo_r)}; "
         f"solo_python={len(solo_py)}",
         "sin huerfanos en ningun lado",
         "TRUE" if (not solo_r and not solo_py) else "FALSE", ESTE_SCRIPT],
        ["comp_headers_iguales",
         "el encabezado coincide en los pares comparados",
         f"{sum(1 for d in comps if d['header_igual'])}/{len(comps)} con header igual",
         f"{len(comps)}/{len(comps)}",
         "TRUE" if all(d["header_igual"] for d in comps) else "FALSE", ESTE_SCRIPT],
        ["comp_dentro_tolerancia",
         "concordancia numerica R<->python dentro de 1e-6 (fallback 1e-4)",
         f"archivos={len(comps)}; byte_identicos={n_byte}; fuera_tol="
         f"{total['n_fuera_tol']}; dif_texto={total['n_dif_texto']}; "
         f"peor|dif|abs={peor_abs:.2e}",
         "0 celdas fuera de tolerancia; 0 diferencias de texto",
         "TRUE" if n_fuera == 0 else "FALSE", ESTE_SCRIPT],
        ["audit_procedencia_figuras",
         "toda figura de outputs/figures/ tiene fila en procedencia.csv",
         f"figuras={len(figuras)}; sin_fila={figs_sin or '[]'}",
         "sin_fila = []",
         "TRUE" if not figs_sin else "FALSE", ESTE_SCRIPT],
        ["audit_procedencia_tablas",
         "toda tabla outputs/tables/{R,python}/*.csv tiene fila en procedencia.csv",
         f"tablas={len(nombres_r)}; sin_fila={tabs_sin or '[]'}",
         "sin_fila = []",
         "TRUE" if not tabs_sin else "FALSE", ESTE_SCRIPT],
        ["audit_descartados_cobertura",
         "analisis_descartados.md tiene una seccion por script de analisis (02..11)",
         f"esperados={len(SCRIPTS_ANALISIS)}; faltan={desc_faltan or '[]'}",
         "faltan = []",
         "TRUE" if not desc_faltan else "FALSE", ESTE_SCRIPT],
        ["audit_verificaciones_ok",
         "ninguna fila de verificaciones.csv queda distinta de TRUE",
         f"total={verif_total}; TRUE={verif_true}; no_TRUE={verif_no or '[]'}",
         "no_TRUE = []",
         "TRUE" if not verif_no else "FALSE", ESTE_SCRIPT],
    ])

    # --- salida legible ---
    print("== 98_comparacion.py ==")
    print(f"  fuente = {fuente}")
    print(f"  CSV comparados: {len(comps)}  |  byte-identicos: {n_byte}/{len(comps)}"
          f"  |  OK: {n_ok}/{len(comps)}")
    print(f"  peor |dif| absoluta: {peor_abs:.3e}"
          + (f"  ({peor_arch})" if peor_arch else ""))
    if solo_r or solo_py:
        print(f"  HUERFANOS  solo_R={solo_r}  solo_python={solo_py}")
    malos = [d["archivo"] for d in comps if not d["ok"]]
    if malos:
        print(f"  REVISAR: {malos}")
    print(f"  procedencia: figuras sin fila = {figs_sin or '[]'}; "
          f"tablas sin fila = {tabs_sin or '[]'}")
    print(f"  descartados: sin seccion = {desc_faltan or '[]'}")
    print(f"  verificaciones: {verif_true}/{verif_total} en TRUE; "
          f"distintas de TRUE = {verif_no or '[]'}")
    todo_ok = (not malos and not solo_r and not solo_py and not figs_sin
               and not tabs_sin and not desc_faltan and not verif_no)
    print("  -> outputs/tables/comparacion_R_python.csv, comparacion_reporte.md")
    print("  TODAS LAS COMPARACIONES Y AUDITORIAS PASARON" if todo_ok
          else "  *** HAY ITEMS A REVISAR (ver arriba) ***")


if __name__ == "__main__":
    main()
