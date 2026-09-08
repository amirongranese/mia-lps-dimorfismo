# 04_qpcr_cuantificacion.py -- Cuantificacion relativa de la qPCR (D1, D2, D8).
#
# Por que existe este archivo: convierte los CT crudos de 02_ingesta_qc en las
# medidas sobre las que modela T5. Aplica -- sin imputar NADA -- las decisiones:
#
#   * D1: dCt = CT_gen - CT_rsp29 por muestra. Calibrador = HEMBRA_CONTROL, por
#     gen x tejido, promediando SOLO valores detectados. ddCt = dCt - dCt_calib.
#   * D2: la medida de analisis es -ddCt (log2, simetrica y aditiva). El
#     fold-change FC = 2^(-ddCt) es solo para graficar (T6/T7): NO se guarda aca
#     para no arrastrar el redondeo de 2^x entre lenguajes.
#   * D7: si el calibrador HEMBRA_CONTROL de un gen x tejido tiene 0 detectados,
#     no hay dCt de calibrador -> ese gen x tejido NO es cuantificable (ddCt/-ddCt/z
#     quedan NA, cuantificable=FALSE). Se detecta programaticamente. Hoy: solo
#     il6 @ BRAIN_E15.
#   * D8: z-score de cada gen por separado, dentro de cada tejido, sobre los 36
#     fetos, usando SOLO los detectados; desvio estandar MUESTRAL (n-1). Score
#     compuesto = promedio de los 7 z de transportadores por feto x tejido
#     (los z disponibles si hay <7). PROHIBIDO llamarlo "TONE".
#
# D3: los no detectados quedan NA y no se imputan por ningun metodo. Este script
# ademas emite qpcr_no_detectados_descriptivo.csv -- conteo y % de NA por
# gen x tejido x grupo + marca de calibrador sin detectados -- que es el sustento
# numerico de D3 y D7 a la vez (la redaccion del descarte MNAR se hace en T10
# apoyandose en esa tabla).
#
# Este script NO modela y NO grafica.

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "04_qpcr_cuantificacion"

CAL_GRUPO = "HEMBRA_CONTROL"          # D1: calibrador
GRUPOS_4 = ["HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
ORDEN_GRUPO = {g: i for i, g in enumerate(GRUPOS_4)}
ORDEN_TEJIDO = {t: i for i, t in enumerate(cfg.TEJIDOS_E15)}
ORDEN_GEN = {g: i for i, g in enumerate(cfg.GENES)}
SET_TRANSP = set(cfg.GENES_TRANSPORTADORES)

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02/03 (paridad R/Python).
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


def escribir_tsv(ruta: Path, encabezado, filas) -> None:
    lineas = ["\t".join(encabezado)]
    lineas += ["\t".join(_fmt(v) for v in fila) for fila in filas]
    ruta.write_bytes(("\n".join(lineas) + "\n").encode("utf-8"))


def escribir_csv(ruta: Path, encabezado, filas) -> None:
    lineas = [",".join(_csv_cell(v) for v in encabezado)]
    lineas += [",".join(_csv_cell(v) for v in fila) for fila in filas]
    ruta.write_bytes(("\n".join(lineas) + "\n").encode("utf-8"))


def escribir_texto(ruta: Path, texto: str) -> None:
    if not texto.endswith("\n"):
        texto += "\n"
    ruta.write_bytes(texto.encode("utf-8"))


# ---------------------------------------------------------------------------
# Aritmetica de promedios/desvios -- acumulador double explicito, mismo orden
# de suma en R y Python (R sum() usa long double -> se evita a proposito).
# ---------------------------------------------------------------------------
def suma(xs):
    a = 0.0
    for v in xs:
        a += v
    return a


def promedio(xs):
    return suma(xs) / len(xs)


def desvio_muestral(xs, media):
    if len(xs) < 2:
        return None
    s = 0.0
    for v in xs:
        s += (v - media) * (v - media)
    return math.sqrt(s / (len(xs) - 1))


def mediana(xs):
    s = sorted(xs)
    k = len(s)
    if k == 0:
        return None
    if k % 2 == 1:
        return s[k // 2]
    return (s[k // 2 - 1] + s[k // 2]) / 2.0


# ===========================================================================
# 1. Lectura de qpcr_e15_long.tsv
# ===========================================================================
def cargar():
    ruta = cfg.RUTA_DATOS_PROC / "qpcr_e15_long.tsv"
    if not ruta.is_file():
        raise SystemExit(f"04: falta {ruta} (correr 02_ingesta_qc primero)")
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    h = lineas[0].split("\t")
    return [dict(zip(h, ln.split("\t"))) for ln in lineas[1:]]


# ===========================================================================
# 2. Estructuras base: metadata de feto, CT, rsp29
# ===========================================================================
def indexar(rows):
    meta = {}
    for r in rows:
        meta.setdefault(r["FETO"],
                        (r["MADRE_ID"], r["FETO"], r["SEXO"], r["TTO"], r["GRUPO"]))
    fetos = sorted(meta.values(), key=lambda m: (m[0], m[1]))

    ct = {}       # (feto, tejido, gen) -> float | None
    no_det = {}   # (feto, tejido, gen) -> bool
    rsp = {}      # (feto, tejido) -> float
    for r in rows:
        k = (r["FETO"], r["TEJIDO"], r["GEN"])
        no_det[k] = r["no_detectado"] == "TRUE"
        ct[k] = None if r["CT"] == "" else float(r["CT"])
        rsp[(r["FETO"], r["TEJIDO"])] = float(r["rsp29"])
    return fetos, ct, no_det, rsp


def dct_de(ct, rsp, feto, tej, gen):
    c = ct[(feto, tej, gen)]
    return None if c is None else c - rsp[(feto, tej)]


# ===========================================================================
# 3. Calibrador (D1) + regla D7
# ===========================================================================
def calcular_calibradores(fetos, ct, rsp):
    cal_fetos = [m[1] for m in fetos if m[4] == CAL_GRUPO]
    calib = {}  # (tej, gen) -> dict(n_total, n_detectado, dct_calibrador, cuantificable)
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            vals = [dct_de(ct, rsp, f, tej, gen) for f in cal_fetos]
            det = [v for v in vals if v is not None]
            dcal = promedio(det) if det else None
            calib[(tej, gen)] = {
                "n_total": len(vals), "n_detectado": len(det),
                "dct_calibrador": dcal, "cuantificable": len(det) > 0,
            }
    return calib


# ===========================================================================
# 4. dCt / ddCt / -ddCt / z (D2, D8) y score compuesto
# ===========================================================================
def cuantificar(fetos, ct, no_det, rsp, calib):
    long_rows = []
    negddct = {}  # (feto, tej, gen) -> float | None
    for m in fetos:
        for tej in cfg.TEJIDOS_E15:
            for gen in cfg.GENES:
                k = (m[1], tej, gen)
                c = ct[k]
                d = None if c is None else c - rsp[(m[1], tej)]
                cal = calib[(tej, gen)]
                dcal = cal["dct_calibrador"]
                dd = None if (d is None or dcal is None) else d - dcal
                nd = None if dd is None else -dd
                negddct[k] = nd
                long_rows.append({
                    "MADRE_ID": m[0], "FETO": m[1], "SEXO": m[2], "TTO": m[3],
                    "GRUPO": m[4], "TEJIDO": tej, "GEN": gen,
                    "es_transportador": gen in SET_TRANSP,
                    "no_detectado": no_det[k], "cuantificable": cal["cuantificable"],
                    "rsp29": rsp[(m[1], tej)], "CT": c,
                    "dCt": d, "dCt_calibrador": dcal, "ddCt": dd, "neg_ddCt": nd,
                })

    # z-score por gen x tejido sobre los 36 fetos (solo detectados, sd n-1)
    z = {}
    zstats = {}  # (tej, gen) -> (n_det, media, sd)
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            pares = [(m[1], negddct[(m[1], tej, gen)]) for m in fetos]
            det = [v for (_f, v) in pares if v is not None]
            if len(det) >= 2:
                mu = promedio(det)
                sd = desvio_muestral(det, mu)
                for (f, v) in pares:
                    z[(f, tej, gen)] = None if (v is None or sd is None or sd == 0) \
                        else (v - mu) / sd
                zstats[(tej, gen)] = (len(det), mu, sd)
            else:
                for (f, _v) in pares:
                    z[(f, tej, gen)] = None
                zstats[(tej, gen)] = (len(det), None, None)

    for row in long_rows:
        row["z"] = z[(row["FETO"], row["TEJIDO"], row["GEN"])]

    # score compuesto por feto x tejido: promedio de los z de transportadores
    score_rows = []
    for m in fetos:
        for tej in cfg.TEJIDOS_E15:
            zs = [z[(m[1], tej, g)] for g in cfg.GENES_TRANSPORTADORES]
            zs = [v for v in zs if v is not None]
            score_rows.append({
                "MADRE_ID": m[0], "FETO": m[1], "SEXO": m[2], "TTO": m[3],
                "GRUPO": m[4], "TEJIDO": tej, "n_z_disponibles": len(zs),
                "score_compuesto": promedio(zs) if zs else None,
            })
    return long_rows, score_rows, zstats


# ===========================================================================
# 5. Tabla descriptiva de no detectados (sustento numerico de D3 y D7)
# ===========================================================================
def tabla_no_detectados(fetos, no_det, calib):
    header = ["TEJIDO", "GEN", "GRUPO", "n_total", "n_no_detectado",
              "pct_no_detectado", "calibrador_cero_detectados"]
    grupo_de = {m[1]: m[4] for m in fetos}
    filas = []
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            cal_cero = not calib[(tej, gen)]["cuantificable"]
            for grp in GRUPOS_4 + ["TODOS"]:
                fs = [m[1] for m in fetos if grp == "TODOS" or m[4] == grp]
                n = len(fs)
                nnd = sum(1 for f in fs if no_det[(f, tej, gen)])
                filas.append([tej, gen, grp, n, nnd,
                              (100.0 * nnd / n) if n else None, cal_cero])
    return header, filas


# ===========================================================================
# 6. Resumenes
# ===========================================================================
def tabla_calibradores(calib):
    header = ["TEJIDO", "GEN", "es_transportador", "n_calibrador",
              "n_calibrador_detectado", "dCt_calibrador", "cuantificable"]
    filas = []
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            c = calib[(tej, gen)]
            filas.append([tej, gen, gen in SET_TRANSP, c["n_total"],
                          c["n_detectado"], c["dct_calibrador"], c["cuantificable"]])
    return header, filas


def tabla_resumen_cuant(long_rows, calib, zstats):
    header = ["TEJIDO", "GEN", "es_transportador", "cuantificable",
              "n_detectado_de_36", "mean_neg_ddCt", "sd_neg_ddCt",
              "median_neg_ddCt"]
    filas = []
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            vals = [r["neg_ddCt"] for r in long_rows
                    if r["TEJIDO"] == tej and r["GEN"] == gen
                    and r["neg_ddCt"] is not None]
            n = len(vals)
            mu = promedio(vals) if n else None
            sd = desvio_muestral(vals, mu) if n >= 2 else None
            md = mediana(vals) if n else None
            filas.append([tej, gen, gen in SET_TRANSP,
                          calib[(tej, gen)]["cuantificable"], n, mu, sd, md])
    return header, filas


def tabla_resumen_score(score_rows):
    header = ["TEJIDO", "GRUPO", "n_fetos", "mean_score_compuesto",
              "sd_score_compuesto"]
    filas = []
    for tej in cfg.TEJIDOS_E15:
        for grp in GRUPOS_4:
            vals = [r["score_compuesto"] for r in score_rows
                    if r["TEJIDO"] == tej and r["GRUPO"] == grp
                    and r["score_compuesto"] is not None]
            n = len(vals)
            mu = promedio(vals) if n else None
            sd = desvio_muestral(vals, mu) if n >= 2 else None
            filas.append([tej, grp, n, mu, sd])
    return header, filas


# ===========================================================================
# 7. Reporte legible
# ===========================================================================
def _md_tabla(header, filas):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for f in filas:
        out.append("| " + " | ".join(_fmt(v) for v in f) + " |")
    return "\n".join(out)


def construir_reporte(fuente, cal_h, cal, nd_h, nd, res_h, res, sc_h, sc,
                      no_cuant):
    L = []
    ap = L.append
    ap("# Reporte de cuantificacion qPCR (T4)")
    ap("")
    ap("Generado por `04_qpcr_cuantificacion` (R y Python producen este archivo "
       "identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Metodo (D1, D2, D8) -- sin imputacion (D3)")
    ap("")
    ap("- `dCt = CT_gen - CT_rsp29` por muestra. Calibrador = promedio de `dCt` "
       "en **HEMBRA_CONTROL**, por gen x tejido, **solo detectados** (D1).")
    ap("- `ddCt = dCt - dCt_calibrador`; la medida de analisis es `-ddCt` (D2). "
       "`FC = 2^(-ddCt)` se calcula al graficar, no se guarda.")
    ap("- z-score por gen dentro de tejido sobre los 36 fetos, solo detectados, "
       "**desvio muestral (n-1)** (D8).")
    ap("- Score compuesto = promedio de los z de los 7 transportadores por "
       "feto x tejido (los disponibles si hay <7) (D8).")
    ap("")
    ap("## 2. Calibradores por gen x tejido (D1) y regla D7")
    ap("")
    if no_cuant:
        ap(f"**No cuantificable(s) por D7** (calibrador HEMBRA_CONTROL con 0 "
           f"detectados): {', '.join(no_cuant)}.")
    else:
        ap("Todos los gen x tejido tienen >=1 detectado en el calibrador.")
    ap("")
    ap(_md_tabla(cal_h, cal))
    ap("")
    ap("## 3. No detectados por gen x tejido x grupo (sustento numerico de D3 y D7)")
    ap("")
    ap("Columna `calibrador_cero_detectados` = TRUE marca el gen x tejido donde "
       "el calibrador HEMBRA_CONTROL no tiene ningun detectado: ahi la imputacion "
       "MNAR no tendria nada sobre lo que anclar (se argumenta en T10). Tabla "
       "completa en `qpcr_no_detectados_descriptivo.csv`; aca solo la fila `TODOS` "
       "por gen x tejido.")
    ap("")
    ap(_md_tabla(nd_h, [f for f in nd if f[2] == "TODOS"]))
    ap("")
    ap("## 4. Resumen de -ddCt por gen x tejido (36 fetos)")
    ap("")
    ap(_md_tabla(res_h, res))
    ap("")
    ap("## 5. Score compuesto por grupo x tejido")
    ap("")
    ap(_md_tabla(sc_h, sc))
    ap("")
    ap("## 6. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `04_qpcr_cuantificacion`: "
       "il6 @ BRAIN_E15 fuera de la cuantificacion (D7), y el recordatorio de que "
       "no se imputa (D3) con el puntero a la tabla de no detectados.")
    ap("")
    return "\n".join(L)


# ===========================================================================
# 8. Artefactos compartidos (merge por 'script') -- headers identicos a 02/03.
# ===========================================================================
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


BLOQUE_DESCARTES = """\
## 04_qpcr_cuantificacion

### il6 @ BRAIN_E15: fuera de la cuantificacion (D7)

- El calibrador HEMBRA_CONTROL de il6 @ BRAIN_E15 tiene 0/9 detectados (ver
  `qpcr_no_detectados_descriptivo.csv`: TEJIDO=BRAIN_E15, GEN=il6,
  GRUPO=HEMBRA_CONTROL, n_no_detectado = n_total). Sin dCt de calibrador no hay
  ddCt: se marca `cuantificable = FALSE`, sus `ddCt` / `-ddCt` / `z` quedan NA,
  y el gen entra al analisis solo como proporcion de deteccion en T5 (D7).
- Es el unico gen x tejido donde ocurre. Se detecta programaticamente:
  gen x tejido con HEMBRA_CONTROL `n_detectado == 0 & n_total > 0`.

### No se imputan los no detectados (D3) -- sustento numerico

- Los CT no detectados (`CT_CRUDO == 40` o celda vacia, unificados en
  02_ingesta_qc) quedan NA y NO se imputan por ningun metodo. El conteo y % de
  NA por gen x tejido x grupo esta en `qpcr_no_detectados_descriptivo.csv`.
- El calibrador y todos los promedios/desvios (dCt de calibrador, media y sd
  del z-score) se calculan **solo sobre valores detectados** (D1/D8). Ningun
  `dCt` se rellena: la cantidad de `dCt` NA es exactamente la de no detectados.
- La justificacion completa del descarte de la imputacion MNAR (`nondetects`) se
  redacta en T10 apoyandose en esta tabla: la imputacion MNAR estima lo no
  observado a partir de un modelo de la censura, y en las celdas donde el
  calibrador no tiene detectados (D7) no hay nada sobre lo que anclar la
  estimacion, asi que el fold-change quedaria definido contra un valor inventado.
"""


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 04_qpcr_cuantificacion:inicio -->"
    marca_fin = "<!-- 04_qpcr_cuantificacion:fin -->"
    nuevo = f"{marca_ini}\n{BLOQUE_DESCARTES}\n{marca_fin}"
    if ruta.is_file():
        txt = ruta.read_bytes().decode("utf-8")
    else:
        txt = ("# Analisis descartados\n\nQue se probo, por que no funciono o no se"
               " uso, y que se hizo en su lugar. Una seccion por script.\n")
    if marca_ini in txt and marca_fin in txt:
        pre = txt.split(marca_ini)[0]
        post = txt.split(marca_fin)[1]
        txt = pre + nuevo + post
    else:
        if not txt.endswith("\n"):
            txt += "\n"
        txt += "\n" + nuevo + "\n"
    escribir_texto(ruta, txt)


# ===========================================================================
def main():
    rows = cargar()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)
    fetos, ct, no_det, rsp = indexar(rows)
    calib = calcular_calibradores(fetos, ct, rsp)
    long_rows, score_rows, zstats = cuantificar(fetos, ct, no_det, rsp, calib)

    no_cuant = [f"{g}@{t}" for t in cfg.TEJIDOS_E15 for g in cfg.GENES
                if not calib[(t, g)]["cuantificable"]]

    # --- intermedios largos (data/processed/) ---------------------------
    cols_long = ["MADRE_ID", "FETO", "SEXO", "TTO", "GRUPO", "TEJIDO", "GEN",
                 "es_transportador", "no_detectado", "cuantificable", "rsp29",
                 "CT", "dCt", "dCt_calibrador", "ddCt", "neg_ddCt", "z"]
    escribir_tsv(cfg.RUTA_DATOS_PROC / "qpcr_cuantificacion_long.tsv", cols_long,
                 [[r[c] for c in cols_long] for r in long_rows])
    cols_sc = ["MADRE_ID", "FETO", "SEXO", "TTO", "GRUPO", "TEJIDO",
               "n_z_disponibles", "score_compuesto"]
    escribir_tsv(cfg.RUTA_DATOS_PROC / "qpcr_score_compuesto_long.tsv", cols_sc,
                 [[r[c] for c in cols_sc] for r in score_rows])

    # --- tablas de salida (una copia por implementacion) ---------------
    cal_h, cal = tabla_calibradores(calib)
    nd_h, nd = tabla_no_detectados(fetos, no_det, calib)
    res_h, res = tabla_resumen_cuant(long_rows, calib, zstats)
    sc_h, sc = tabla_resumen_score(score_rows)
    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "qpcr_calibradores.csv", cal_h, cal)
        escribir_csv(base / "qpcr_no_detectados_descriptivo.csv", nd_h, nd)
        escribir_csv(base / "qpcr_cuantificacion_resumen.csv", res_h, res)
        escribir_csv(base / "qpcr_score_compuesto_resumen.csv", sc_h, sc)

    escribir_texto(cfg.RUTA_TABLAS / "qpcr_cuantificacion_reporte.md",
                   construir_reporte(fuente, cal_h, cal, nd_h, nd, res_h, res,
                                     sc_h, sc, no_cuant))
    actualizar_descartados()

    # --- verificaciones internas (sin imputar) ------------------------
    n_dct_na = sum(1 for r in long_rows if r["dCt"] is None)
    n_nodet = sum(1 for r in long_rows if r["no_detectado"])

    ent = f"data/processed/qpcr_e15_long.tsv (de data/{fuente}/{cfg.ARCHIVO_QPCR})"
    registrar_procedencia([
        ["data/processed/qpcr_cuantificacion_long.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "dCt, calibrador, ddCt, -ddCt, z por feto x tejido x gen (D1/D2/D8)"],
        ["data/processed/qpcr_score_compuesto_long.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "score compuesto (promedio z de 7 transportadores) por feto x tejido (D8)"],
        ["outputs/tables/{R,python}/qpcr_calibradores.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "dCt de calibrador HEMBRA_CONTROL por gen x tejido + flag cuantificable (D1/D7)"],
        ["outputs/tables/{R,python}/qpcr_no_detectados_descriptivo.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "n y % de no detectados por gen x tejido x grupo + calibrador_cero_detectados (sustento D3/D7)"],
        ["outputs/tables/{R,python}/qpcr_cuantificacion_resumen.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "n detectados, media/sd/mediana de -ddCt por gen x tejido"],
        ["outputs/tables/{R,python}/qpcr_score_compuesto_resumen.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "media/sd del score compuesto por grupo x tejido"],
        ["outputs/tables/qpcr_cuantificacion_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible de la cuantificacion qPCR"],
    ])
    registrar_verificaciones([
        ["cuant_filas_long", "qpcr_cuantificacion_long tiene 36x2x10 = 720 filas",
         str(len(long_rows)), "720", "TRUE" if len(long_rows) == 720 else "FALSE",
         ESTE_SCRIPT],
        ["cuant_score_filas", "score compuesto: 36x2 = 72 filas",
         str(len(score_rows)), "72", "TRUE" if len(score_rows) == 72 else "FALSE",
         ESTE_SCRIPT],
        ["cuant_no_cuantificables", "gen x tejido no cuantificables por D7",
         ";".join(no_cuant) if no_cuant else "(ninguno)", "il6@BRAIN_E15",
         "TRUE" if no_cuant == ["il6@BRAIN_E15"] else "FALSE", ESTE_SCRIPT],
        ["cuant_sin_imputacion", "cantidad de dCt NA == cantidad de no detectados (no se imputa)",
         f"{n_dct_na}=={n_nodet}", "iguales", "TRUE" if n_dct_na == n_nodet else "FALSE",
         ESTE_SCRIPT],
        ["cuant_sd_metodo", "z-score con desvio estandar muestral (n-1)",
         "muestral_n-1", "muestral_n-1", "TRUE", ESTE_SCRIPT],
        ["cuant_transportadores_score", "score compuesto sobre 7 transportadores",
         str(len(cfg.GENES_TRANSPORTADORES)), "7",
         "TRUE" if len(cfg.GENES_TRANSPORTADORES) == 7 else "FALSE", ESTE_SCRIPT],
        ["cuant_calibrador_grupo", "calibrador de la cuantificacion relativa (D1)",
         CAL_GRUPO, "HEMBRA_CONTROL", "TRUE" if CAL_GRUPO == "HEMBRA_CONTROL" else "FALSE",
         ESTE_SCRIPT],
    ])

    print("== 04_qpcr_cuantificacion.py ==")
    print(f"  fuente qPCR = {fuente}")
    print(f"  filas long = {len(long_rows)}   score compuesto filas = {len(score_rows)}")
    print(f"  no cuantificable (D7): {', '.join(no_cuant) if no_cuant else '(ninguno)'}")
    print(f"  dCt NA = {n_dct_na}  (== no detectados = {n_nodet}) -> sin imputar")
    for f in sc:
        print(f"    score {f[0]:12} {f[1]:14} n={f[2]:2}  media={_fmt(f[3])}  sd={_fmt(f[4])}")
    print("  -> data/processed/qpcr_{cuantificacion,score_compuesto}_long.tsv,")
    print("     outputs/tables/{R,python}/qpcr_*.csv, outputs/tables/qpcr_cuantificacion_reporte.md")


if __name__ == "__main__":
    main()
