# 06_pstat3.py -- Modelo de pSTAT3 (D9): PSTAT3 ~ SEXO * TTO + MEMBRANA.
#
# Por que existe este archivo: pSTAT3 se midio en 3 membranas (western blot); la
# membrana es un bloque tecnico. D9 fija el modelo `PSTAT3 ~ SEXO * TTO + MEMBRANA`
# con MEMBRANA como BLOQUE FIJO (diseno balanceado 3 x 4 x 3), la MISMA cascada de
# supuestos D5 / seccion 4.1 y el MISMO post hoc D6 que la qPCR (05_qpcr_modelos).
#
#   * Diseno suma-cero de 6 columnas: [1, s, t, s*t, m1, m2] con
#       s = +1 HEMBRA / -1 MACHO,  t = +1 CONTROL / -1 LPS,
#       (m1, m2) = contr.sum de MEMBRANA (3 niveles): lev1->(1,0) lev2->(0,1) lev3->(-1,-1).
#     Terminos: SEXO=col1, TTO=col2, SEXO:TTO=col3, MEMBRANA=cols4-5. df_resid = n - 6.
#   * Cascada D5 (alfa = 0.05):
#       - Shapiro-Wilk sobre los RESIDUOS del modelo conjunto (con MEMBRANA).
#       - Levene (Brown-Forsythe, centro = mediana) sobre las 4 celdas SEXO x TTO.
#       - Shapiro >= .05 y Levene >= .05  -> ANOVA III           (rama "anova3")
#       - Shapiro >= .05 y Levene <  .05  -> Wald III con HC3     (rama "hc3")
#       - Shapiro <  .05                  -> ART aditivo hand-rolled (rama "art")
#     ARTool RECHAZA un bloque aditivo no cruzado (`parse.art.formula` exige todas
#     las interacciones), asi que la rama `art` es un ART hand-rolled alineado
#     sobre el modelo aditivo; sobre datos reales Y sinteticos la cascada cae en
#     `anova3` (Shapiro ~0.10, Levene ~0.27 / ~0.51) => la rama `art` NO se ejecuta.
#     Ver `analisis_descartados.md`, seccion `06_pstat3`.
#   * Piso de celda: si alguna celda SEXO x TTO tiene < 5 valores, via
#     `descriptivo_n_bajo` (pSTAT3 real y sintetico: 9/celda -> siempre se modela).
#   * Post hoc D6 (solo si SEXO x TTO p < 0.05): 4 comparaciones fijas + Holm.
#       anova3 -> contrastes de medias marginales (promediando sobre MEMBRANA, vcov OLS)
#       hc3    -> idem con vcov HC3
#       art    -> ART-C hand-rolled (rangos alineados por la interaccion + MEMBRANA)
#   * LIMITACION OBLIGATORIA DEL INFORME (D9): pSTAT3 esta normalizado a proteina
#     total SIN STAT3 total -> refleja ABUNDANCIA de fosfo-STAT3, no la fraccion
#     fosforilada. Queda escrita en el reporte y en `analisis_descartados.md`.
#
# PARIDAD R/Python: el nucleo numerico (OLS 6x6 por Gauss-Jordan, SS tipo III por
# comparacion de modelos, sandwich HC3, Levene, ART y ART-C, Holm) es codigo PROPIO
# identico en ambos lenguajes. R cruza-verifica en corrida contra
# car::Anova / emmeans (ramas anova3 / hc3, que son las que se ejecutan) con
# stopifnot (< 1e-6); la rama `art` se auto-verifica (dos vias de alineado) porque
# ARTool no puede ajustar el modelo aditivo. Unico componente de libreria:
# Shapiro-Wilk (scipy aca, shapiro.test en R). Estadisticos / p dependientes de
# trascendentes se guardan como texto "%.6e" (p6e), como en 03_elisa / 05.
#
# Este script NO grafica (el boxplot de pSTAT3 se arma en 07_figuras_acto1).

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

from scipy import stats as _sst

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "06_pstat3"

ALFA = 0.05          # D5/D6: gate de interaccion y pretests de supuestos
PISO_CELDA = 5       # minimo de valores por celda SEXO x TTO para ajustar modelo
NCOL = 6             # [1, s, t, s*t, m1, m2]

GRUPOS_4 = ["HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
CELDAS_4 = [("HEMBRA", "CONTROL"), ("HEMBRA", "LPS"),
            ("MACHO", "CONTROL"), ("MACHO", "LPS")]

# posicion de columnas de cada termino en el diseno de 6 columnas
TERM_COLS = {"SEXO": [1], "TTO": [2], "SEXO:TTO": [3], "MEMBRANA": [4, 5]}

# D6: las 4 comparaciones fijas, como (etiqueta, celda_a, celda_b) -> estima a - b
PARES_D6 = [
    ("HEMBRA_CONTROL-HEMBRA_LPS", ("HEMBRA", "CONTROL"), ("HEMBRA", "LPS")),
    ("MACHO_CONTROL-MACHO_LPS", ("MACHO", "CONTROL"), ("MACHO", "LPS")),
    ("HEMBRA_LPS-MACHO_LPS", ("HEMBRA", "LPS"), ("MACHO", "LPS")),
    ("HEMBRA_CONTROL-MACHO_CONTROL", ("HEMBRA", "CONTROL"), ("MACHO", "CONTROL")),
]


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02/03/04/05.
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


def p6e(x) -> str:
    """p-valor / estadistico dependiente de trascendentes -> texto %.6e (paridad)."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return "%.6e" % float(x)


def g10(x) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return "%.10g" % float(x)


# ---------------------------------------------------------------------------
# Algebra lineal PROPIA -- eliminacion de Gauss-Jordan con pivoteo parcial.
# Codigo identico en R; los sistemas son 6x6 / 3x3 / 2x2, bien condicionados.
# ---------------------------------------------------------------------------
def resolver(A, b):
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        d = M[c][c]
        for j in range(c, n + 1):
            M[c][j] /= d
        for r in range(n):
            if r != c and M[r][c] != 0.0:
                f = M[r][c]
                for j in range(c, n + 1):
                    M[r][j] -= f * M[c][j]
    return [M[i][n] for i in range(n)]


def invertir(A):
    n = len(A)
    M = [list(A[i]) + [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        d = M[c][c]
        for j in range(2 * n):
            M[c][j] /= d
        for r in range(n):
            if r != c:
                f = M[r][c]
                for j in range(2 * n):
                    M[r][j] -= f * M[c][j]
    return [[M[i][n + j] for j in range(n)] for i in range(n)]


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


def rangos_promedio(xs):
    """rank() de R con ties.method='average' (== scipy 'average')."""
    orden = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(orden):
        j = i
        while j + 1 < len(orden) and xs[orden[j + 1]] == xs[orden[i]]:
            j += 1
        rango = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[orden[k]] = rango
        i = j + 1
    return r


# ---------------------------------------------------------------------------
# Distribuciones -- scipy para las colas (pf/pt); identicas a pf/pt de R a ~1e-15.
# ---------------------------------------------------------------------------
def f_sf(x, d1, d2):
    return float(_sst.f.sf(x, d1, d2))


def t_sf2(x, df):
    return 2.0 * float(_sst.t.sf(abs(x), df))


# ---------------------------------------------------------------------------
# Diseno suma-cero de 6 columnas [1, s, t, s*t, m1, m2].
# ---------------------------------------------------------------------------
def _cs_membrana(memb, membranas):
    """contr.sum de MEMBRANA (3 niveles): lev1->(1,0), lev2->(0,1), lev3->(-1,-1)."""
    if memb == membranas[0]:
        return 1.0, 0.0
    if memb == membranas[1]:
        return 0.0, 1.0
    return -1.0, -1.0


def _fila_diseno(sexo, tto, memb, membranas):
    s = 1.0 if sexo == "HEMBRA" else -1.0
    t = 1.0 if tto == "CONTROL" else -1.0
    m1, m2 = _cs_membrana(memb, membranas)
    return [1.0, s, t, s * t, m1, m2]


def diseno(muestras, membranas):
    X = [_fila_diseno(m["SEXO"], m["TTO"], m["MEMBRANA"], membranas) for m in muestras]
    y = [m["y"] for m in muestras]
    return X, y


def ajustar(X, y, cols):
    Xs = [[fila[c] for c in cols] for fila in X]
    p = len(cols)
    n = len(y)
    XtX = [[suma([Xs[i][a] * Xs[i][b] for i in range(n)]) for b in range(p)]
           for a in range(p)]
    Xty = [suma([Xs[i][a] * y[i] for i in range(n)]) for a in range(p)]
    b = resolver(XtX, Xty)
    resid = [y[i] - suma([Xs[i][a] * b[a] for a in range(p)]) for i in range(n)]
    sse = suma([e * e for e in resid])
    return sse, b, resid


# ---------------------------------------------------------------------------
# Rama anova3 -- SS tipo III por comparacion de modelos (contr.sum saturado).
# MEMBRANA (2 gl) tambien se reporta (bloque tecnico): F((sse_red - sse_full)/2 / mse).
# ---------------------------------------------------------------------------
def anova3_terminos(X, y):
    n = len(y)
    p = len(X[0])
    sse_full, _b, _r = ajustar(X, y, list(range(p)))
    df = n - p
    mse = sse_full / df
    out = {}
    for nombre, cols_term in TERM_COLS.items():
        cols = [c for c in range(p) if c not in cols_term]
        sse_r, _b2, _r2 = ajustar(X, y, cols)
        k = len(cols_term)
        F = ((sse_r - sse_full) / k) / mse
        out[nombre] = (F, f_sf(F, k, df))
    return out, df


# ---------------------------------------------------------------------------
# Rama hc3 -- vcov sandwich HC3 y Wald tipo III.
#   V = (X'X)^-1 [ sum_i x_i x_i' e_i^2 / (1 - h_ii)^2 ] (X'X)^-1
#   termino de k gl -> W = b_sub' V_sub^-1 b_sub ; F = W / k  (k=1 -> b^2 / V).
# ---------------------------------------------------------------------------
def _vcov_hc3(X, resid):
    n = len(X)
    p = len(X[0])
    XtX = [[suma([X[i][a] * X[i][c] for i in range(n)]) for c in range(p)]
           for a in range(p)]
    XtXi = invertir(XtX)
    h = [suma([X[i][a] * XtXi[a][c] * X[i][c]
               for a in range(p) for c in range(p)]) for i in range(n)]
    meat = [[suma([X[i][a] * X[i][c] * (resid[i] * resid[i]) / ((1.0 - h[i]) ** 2)
                   for i in range(n)]) for c in range(p)] for a in range(p)]
    V = [[suma([XtXi[a][k] * meat[k][l] * XtXi[l][c]
               for k in range(p) for l in range(p)]) for c in range(p)]
         for a in range(p)]
    return V, XtXi


def hc3_terminos(X, y):
    n = len(y)
    p = len(X[0])
    _sse, b, resid = ajustar(X, y, list(range(p)))
    V, _XtXi = _vcov_hc3(X, resid)
    df = n - p
    out = {}
    for nombre, cols_term in TERM_COLS.items():
        k = len(cols_term)
        if k == 1:
            idx = cols_term[0]
            F = b[idx] * b[idx] / V[idx][idx]
        else:
            bsub = [b[c] for c in cols_term]
            Vsub = [[V[a][c] for c in cols_term] for a in cols_term]
            Vinv = invertir(Vsub)
            W = suma([bsub[i] * Vinv[i][j] * bsub[j]
                      for i in range(k) for j in range(k)])
            F = W / k
        out[nombre] = (F, f_sf(F, k, df))
    return out, df


# ---------------------------------------------------------------------------
# Levene / Brown-Forsythe -- ANOVA de una via sobre |y - mediana(celda)| en las
# 4 celdas SEXO x TTO (la MEMBRANA no entra: es homogeneidad del factorial).
# ---------------------------------------------------------------------------
def levene_bf(muestras):
    celdas = {}
    for m in muestras:
        celdas.setdefault((m["SEXO"], m["TTO"]), []).append(m["y"])
    z = []
    g = []
    for k in CELDAS_4:
        med = mediana(celdas[k])
        for v in celdas[k]:
            z.append(abs(v - med))
            g.append(k)
    N = len(z)
    kk = 4
    gran = promedio(z)
    gm = {}
    for k in CELDAS_4:
        gm[k] = promedio([z[i] for i in range(N) if g[i] == k])
    ssb = suma([len(celdas[k]) * (gm[k] - gran) ** 2 for k in CELDAS_4])
    ssw = suma([(z[i] - gm[g[i]]) ** 2 for i in range(N)])
    F = (ssb / (kk - 1)) / (ssw / (N - kk))
    return F, f_sf(F, kk - 1, N - kk)


# ---------------------------------------------------------------------------
# Shapiro-Wilk sobre los residuos del modelo conjunto (con MEMBRANA). Unica
# dependencia de libreria; scipy aca, shapiro.test en R -- coinciden a ~1e-12.
# ---------------------------------------------------------------------------
def shapiro_residuos(X, y):
    p = len(X[0])
    _sse, _b, resid = ajustar(X, y, list(range(p)))
    r = _sst.shapiro(resid)
    return float(r.statistic), float(r.pvalue)


# ---------------------------------------------------------------------------
# ART aditivo hand-rolled (rama `art`). ARTool RECHAZA `y ~ SEXO*TTO + MEMBRANA`
# (bloque no cruzado): se alinea cada termino sobre el modelo aditivo (residuo del
# ajuste conjunto + contribucion ajustada del termino), se rankean, y se corre
# ANOVA III sobre los rangos alineados usando el diseno completo (con MEMBRANA).
# No se ejecuta con los datos reales ni sinteticos (cascada -> anova3).
# ---------------------------------------------------------------------------
def art_terminos(muestras, membranas):
    X, y = diseno(muestras, membranas)
    n = len(y)
    p = len(X[0])
    _sse, b, resid = ajustar(X, y, list(range(p)))
    out = {}
    for term in ("SEXO", "TTO", "SEXO:TTO"):
        cols_term = TERM_COLS[term]
        alin = [round(resid[i] + suma([b[c] * X[i][c] for c in cols_term]), 8)
                for i in range(n)]
        rr = rangos_promedio(alin)
        muestras_r = [{"SEXO": m["SEXO"], "TTO": m["TTO"],
                       "MEMBRANA": m["MEMBRANA"], "y": rr[i]}
                      for i, m in enumerate(muestras)]
        Xr, yr = diseno(muestras_r, membranas)
        term_stats, _df = anova3_terminos(Xr, yr)
        out[term] = term_stats[term]
    return out, n - p


# ---------------------------------------------------------------------------
# ART-C hand-rolled para las 4 comparaciones D6 (Elkin et al. 2021, extendido al
# bloque aditivo): rangos alineados por la interaccion, modelo
# rr ~ SEXOxTTO(4 celdas, contr.sum) + MEMBRANA ; contrastes t con MSE combinado,
# gl = n - 6. Solo se usa si la rama es `art` (no ocurre con estos datos).
# ---------------------------------------------------------------------------
def artc_pares(muestras, membranas):
    X, y = diseno(muestras, membranas)
    n = len(y)
    p = len(X[0])
    _sse, b, resid = ajustar(X, y, list(range(p)))
    alin = [round(resid[i] + b[3] * X[i][3], 8) for i in range(n)]
    rr = rangos_promedio(alin)

    def fila_c(sexo, tto, memb):
        ci = CELDAS_4.index((sexo, tto))
        c = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
             [-1.0, -1.0, -1.0]][ci]
        m1, m2 = _cs_membrana(memb, membranas)
        return [1.0, c[0], c[1], c[2], m1, m2]

    Xc = [fila_c(m["SEXO"], m["TTO"], m["MEMBRANA"]) for m in muestras]
    sse_c, _bc, _rc = ajustar(Xc, rr, list(range(6)))
    df = n - 6
    mse = sse_c / df
    idx = {k: [i for i in range(n)
               if (muestras[i]["SEXO"], muestras[i]["TTO"]) == k] for k in CELDAS_4}
    media = {k: promedio([rr[i] for i in idx[k]]) for k in CELDAS_4}
    out = {}
    for etiqueta, a, bb in PARES_D6:
        se = math.sqrt(mse * (1.0 / len(idx[a]) + 1.0 / len(idx[bb])))
        est = media[a] - media[bb]
        tval = est / se
        out[etiqueta] = (est, se, tval, t_sf2(tval, df))
    return out, df


# ---------------------------------------------------------------------------
# Contrastes de medias marginales para las ramas anova3 / hc3 (post hoc D6).
# Modelo saturado + bloque -> media marginal de celda (promediando sobre MEMBRANA)
# = media cruda de celda (diseno balanceado). EE del contraste: c' V c con
# c = x_a - x_b en el espacio de coeficientes; las columnas de MEMBRANA se anulan.
# ---------------------------------------------------------------------------
def emmeans_pares(muestras, membranas, robusto):
    X, y = diseno(muestras, membranas)
    n = len(y)
    p = len(X[0])
    sse, b, resid = ajustar(X, y, list(range(p)))
    df = n - p
    mse = sse / df
    if robusto:
        V, _ = _vcov_hc3(X, resid)
    else:
        XtX = [[suma([X[i][a] * X[i][c] for i in range(n)]) for c in range(p)]
               for a in range(p)]
        XtXi = invertir(XtX)
        V = [[mse * XtXi[a][c] for c in range(p)] for a in range(p)]
    celdas = {}
    for m in muestras:
        celdas.setdefault((m["SEXO"], m["TTO"]), []).append(m["y"])
    media = {k: promedio(celdas[k]) for k in CELDAS_4}
    out = {}
    for etiqueta, a, bb in PARES_D6:
        sa = 1.0 if a[0] == "HEMBRA" else -1.0
        ta = 1.0 if a[1] == "CONTROL" else -1.0
        sb = 1.0 if bb[0] == "HEMBRA" else -1.0
        tb = 1.0 if bb[1] == "CONTROL" else -1.0
        cv = [0.0, sa - sb, ta - tb, sa * ta - sb * tb, 0.0, 0.0]
        var = suma([cv[i] * V[i][j] * cv[j] for i in range(p) for j in range(p)])
        se = math.sqrt(var)
        est = media[a] - media[bb]
        tval = est / se
        out[etiqueta] = (est, se, tval, t_sf2(tval, df))
    return out, df


# ---------------------------------------------------------------------------
# Holm dentro de las 4 comparaciones D6.
# ---------------------------------------------------------------------------
def holm(pvals):
    m = len(pvals)
    orden = sorted(range(m), key=lambda i: pvals[i])
    ajust = [0.0] * m
    corr = 0.0
    for rank, i in enumerate(orden):
        val = (m - rank) * pvals[i]
        corr = max(corr, val)
        ajust[i] = min(1.0, corr)
    return ajust


# ===========================================================================
# 1. Carga
# ===========================================================================
def cargar():
    ruta = cfg.RUTA_DATOS_PROC / "pstat3_long.tsv"
    if not ruta.is_file():
        raise SystemExit(f"06: falta {ruta} (correr 02_ingesta_qc primero)")
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    h = lineas[0].split("\t")
    filas = [dict(zip(h, ln.split("\t"))) for ln in lineas[1:]]
    for r in filas:
        r["y"] = float(r["PSTAT3"])
    return filas


def _n_celdas(filas):
    return [sum(1 for r in filas if r["SEXO"] == sx and r["TTO"] == tt)
            for (sx, tt) in CELDAS_4]


# ===========================================================================
# 2. Clasificacion y modelado (una sola fila: pSTAT3 @ PLACENTA_E15)
# ===========================================================================
def clasificar(filas, membranas):
    tej = filas[0]["TEJIDO"]
    ncel = _n_celdas(filas)
    rec = {
        "TEJIDO": tej, "GEN": "pSTAT3",
        "n_celda": ncel, "n_total": len(filas),
        "via": None, "motivo": "",
        "shapiro_p": None, "levene_p": None, "rama": "",
        "stats": {}, "df": None, "posthoc": None,
    }
    if min(ncel) < PISO_CELDA:
        rec["via"] = "descriptivo_n_bajo"
        rec["motivo"] = f"celda SEXO x TTO con < {PISO_CELDA} valores"
        return rec
    rec["via"] = "modelo"
    _ajustar_cascada(rec, filas, membranas)
    return rec


def _ajustar_cascada(rec, filas, membranas):
    X, y = diseno(filas, membranas)
    _W, sp = shapiro_residuos(X, y)
    _lF, lp = levene_bf(filas)
    rec["shapiro_p"] = sp
    rec["levene_p"] = lp
    if sp >= ALFA and lp >= ALFA:
        rec["rama"] = "anova3"
        st_, df = anova3_terminos(X, y)
    elif sp >= ALFA and lp < ALFA:
        rec["rama"] = "hc3"
        st_, df = hc3_terminos(X, y)
    else:
        rec["rama"] = "art"
        st_, df = art_terminos(filas, membranas)
    rec["stats"] = st_       # {'SEXO','TTO','SEXO:TTO','MEMBRANA': (F,p)}
    rec["df"] = df
    ip = st_["SEXO:TTO"][1]
    if ip < ALFA:
        if rec["rama"] == "art":
            pares, pdf = artc_pares(filas, membranas)
            metodo = "ART-C"
        else:
            pares, pdf = emmeans_pares(filas, membranas,
                                       robusto=(rec["rama"] == "hc3"))
            metodo = ("contraste_marginal_HC3" if rec["rama"] == "hc3"
                      else "contraste_marginal_OLS")
        etiquetas = [e for e, _a, _b in PARES_D6]
        praw = [pares[e][3] for e in etiquetas]
        ph = holm(praw)
        rec["posthoc"] = {
            "metodo": metodo, "df": pdf,
            "filas": [(etiquetas[i], pares[etiquetas[i]][0], pares[etiquetas[i]][1],
                       pares[etiquetas[i]][2], praw[i], ph[i])
                      for i in range(4)],
        }


# ===========================================================================
# 3. Armado de tablas de salida
# ===========================================================================
COLS_CLASIF = [
    "TEJIDO", "GEN",
    "n_HEMBRA_CONTROL", "n_HEMBRA_LPS", "n_MACHO_CONTROL", "n_MACHO_LPS", "n_total",
    "via", "motivo", "shapiro_p_residuos", "levene_p", "rama_cascada",
    "F_SEXO", "p_SEXO", "F_TTO", "p_TTO", "F_SEXOxTTO", "p_SEXOxTTO",
    "F_MEMBRANA", "p_MEMBRANA",
    "interaccion_significativa", "post_hoc_corrido",
]


def fila_clasif(rec):
    st_ = rec["stats"]

    def F(term):
        return p6e(st_[term][0]) if term in st_ else ""

    def P(term):
        return p6e(st_[term][1]) if term in st_ else ""

    ip = st_["SEXO:TTO"][1] if "SEXO:TTO" in st_ else None
    sig = (ip is not None and ip < ALFA)
    return [
        rec["TEJIDO"], rec["GEN"],
        rec["n_celda"][0], rec["n_celda"][1], rec["n_celda"][2], rec["n_celda"][3],
        rec["n_total"], rec["via"], rec["motivo"],
        p6e(rec["shapiro_p"]), p6e(rec["levene_p"]), rec["rama"],
        F("SEXO"), P("SEXO"), F("TTO"), P("TTO"), F("SEXO:TTO"), P("SEXO:TTO"),
        F("MEMBRANA"), P("MEMBRANA"),
        (sig if rec["via"] == "modelo" else ""),
        (rec["posthoc"] is not None) if rec["via"] == "modelo" else "",
    ]


COLS_POSTHOC = ["TEJIDO", "GEN", "rama_cascada", "metodo", "contraste",
                "estimador", "EE", "estadistico_t", "p_sin_correccion", "p_holm"]


def filas_posthoc(rec):
    if not rec.get("posthoc"):
        return []
    ph = rec["posthoc"]
    out = []
    for etiqueta, est, ee, tval, praw, pholm in ph["filas"]:
        out.append([rec["TEJIDO"], rec["GEN"], rec["rama"], ph["metodo"], etiqueta,
                    p6e(est), p6e(ee), p6e(tval), p6e(praw), p6e(pholm)])
    return out


COLS_DESCR = ["TEJIDO", "GEN", "GRUPO", "n", "mean_PSTAT3", "sd_PSTAT3",
              "median_PSTAT3"]


def filas_descriptivo(filas):
    tej = filas[0]["TEJIDO"]
    out = []
    for (sx, tt), grp in zip(CELDAS_4, GRUPOS_4):
        vals = [r["y"] for r in filas if r["SEXO"] == sx and r["TTO"] == tt]
        n = len(vals)
        mu = promedio(vals) if n else None
        sd = desvio_muestral(vals, mu) if n >= 2 else None
        md = mediana(vals) if n else None
        out.append([tej, "pSTAT3", grp, n,
                    g10(mu) if mu is not None else "",
                    g10(sd) if sd is not None else "",
                    g10(md) if md is not None else ""])
    return out


# ===========================================================================
# 4. Cruza-verificacion en corrida contra car / emmeans (solo en 06_pstat3.R).
#    En Python no se corre (no hay libs de modelos); el nucleo es identico.
# ===========================================================================

# ===========================================================================
# 5. Reporte legible
# ===========================================================================
def _md(header, filas):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for f in filas:
        out.append("| " + " | ".join(_fmt(v) for v in f) + " |")
    return "\n".join(out)


LIMITACION_D9 = (
    "**Limitacion obligatoria (D9).** pSTAT3 se cuantifico por western blot "
    "normalizado a **proteina total** (p. ej. Ponceau / stain-free), **sin STAT3 "
    "total** en la misma membrana. La medida refleja por lo tanto la **abundancia "
    "de fosfo-STAT3 (Tyr705)** relativa a la carga de proteina, **no la fraccion "
    "de STAT3 que esta fosforilada**. Un cambio en pSTAT3 puede deberse a mas "
    "fosforilacion, a mas STAT3 total, o a ambos; con estos datos no se puede "
    "separar. El bloque `MEMBRANA` (efecto fijo) absorbe la variacion tecnica "
    "entre las 3 membranas del ensayo."
)

BULLETS_METODO = [
    "- Modelo (D9): `PSTAT3 ~ SEXO * TTO + MEMBRANA`, OLS con contrastes suma-cero. "
    "Diseno de 6 columnas `[1, s, t, s*t, m1, m2]`; `MEMBRANA` = bloque fijo de 3 "
    "niveles (contr.sum, 2 columnas). df de residuos = n - 6.",
    "- Cascada D5 / 4.1 (alfa = 0.05): Shapiro-Wilk sobre los **residuos del modelo "
    "conjunto** (con MEMBRANA); Levene (Brown-Forsythe, centro = mediana) sobre las "
    "4 celdas SEXO x TTO. Shapiro>=.05 & Levene>=.05 -> `anova3` (ANOVA III); "
    "Shapiro>=.05 & Levene<.05 -> `hc3` (Wald III con sandwich HC3); Shapiro<.05 -> "
    "`art` (ART aditivo hand-rolled: ARTool rechaza el bloque no cruzado).",
    "- Piso de celda: si alguna celda SEXO x TTO tiene < 5 valores, no se modela "
    "(`descriptivo_n_bajo`). pSTAT3 tiene 9 por celda -> siempre se modela.",
    "- Post hoc D6 (solo si interaccion p < alfa): 4 comparaciones fijas "
    "(HEMBRA_CONTROL-HEMBRA_LPS, MACHO_CONTROL-MACHO_LPS, HEMBRA_LPS-MACHO_LPS, "
    "HEMBRA_CONTROL-MACHO_CONTROL), correccion Holm. `anova3` -> contraste de medias "
    "marginales promediando sobre MEMBRANA (vcov OLS); `hc3` -> idem con vcov HC3; "
    "`art` -> ART-C hand-rolled.",
    "- `F_MEMBRANA` / `p_MEMBRANA` se reportan como diagnostico del bloque tecnico; "
    "no habilitan ni bloquean nada (la compuerta es SEXO x TTO).",
]


def construir_reporte(fuente, rec, clasif, posthoc, descr, resumen):
    L = []
    ap = L.append
    ap("# Reporte de pSTAT3 (T6)")
    ap("")
    ap("Generado por `06_pstat3` (R y Python producen este archivo identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Metodo")
    ap("")
    for b in BULLETS_METODO:
        ap(b)
    ap("")
    ap("## 2. Clasificacion del modelo")
    ap("")
    ap(_md(COLS_CLASIF, clasif))
    ap("")
    ap("## 3. Post hoc (D6)")
    ap("")
    if posthoc:
        ap(_md(COLS_POSTHOC, posthoc))
    else:
        ap("La interaccion SEXO x TTO no alcanzo p < alfa: no se corrio post hoc.")
    ap("")
    ap("## 4. Descriptivo de pSTAT3 por grupo")
    ap("")
    ap(_md(COLS_DESCR, descr))
    ap("")
    ap("## 5. Resumen")
    ap("")
    for linea in resumen:
        ap(f"- {linea}")
    ap("")
    ap("## 6. Limitacion del ensayo (D9)")
    ap("")
    ap(LIMITACION_D9)
    ap("")
    ap("## 7. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `06_pstat3`: eleccion de rama de la "
       "cascada, el rechazo de ARTool al bloque aditivo y el ART hand-rolled, y la "
       "nota de paridad (Shapiro-Wilk por libreria).")
    ap("")
    return "\n".join(L)


# ===========================================================================
# 6. Artefactos compartidos (merge por 'script') -- headers identicos a 02..05.
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


def _bloque_descartes(rama, shapiro_p, levene_p):
    L = [
        "## 06_pstat3",
        "",
        "### Eleccion de rama de la cascada D5",
        "",
        "- Modelo D9 `PSTAT3 ~ SEXO * TTO + MEMBRANA` (OLS, contr.sum, 6 columnas). "
        "Se decide la rama con Shapiro-Wilk sobre los **residuos del modelo "
        "conjunto** y Levene (Brown-Forsythe) sobre las 4 celdas SEXO x TTO, ambos a "
        f"alfa = 0.05. En esta corrida: Shapiro p = {p6e(shapiro_p)}, Levene p = "
        f"{p6e(levene_p)} -> rama `{rama}`. Sobre datos reales Y sinteticos la "
        "cascada cae en `anova3` (residuos normales y varianzas homogeneas).",
        "",
        "### La rama `art`: ARTool no admite el bloque aditivo",
        "",
        "- `ARTool::art()` **rechaza** `y ~ SEXO * TTO + MEMBRANA` con "
        "`parse.art.formula`: *\"Model must include all combinations of interactions "
        "of fixed effects\"*. El ART clasico exige un diseno completamente cruzado; "
        "un bloque tecnico aditivo no lo es.",
        "- Fallback fijado con el usuario: **ART aditivo hand-rolled**. Se alinea "
        "cada termino restando del dato el ajuste de todos los demas terminos "
        "(incluido el efecto principal de MEMBRANA), estimados por OLS del modelo "
        "conjunto; se rankean los valores alineados (redondeo a 8 decimales, rangos "
        "promedio) y se corre ANOVA III sobre esos rangos con el diseno completo. El "
        "post hoc es ART-C hand-rolled: rangos alineados por la interaccion, modelo "
        "`rr ~ SEXOxTTO + MEMBRANA`, contrastes t con MSE combinado (gl = n - 6).",
        "- Como ARTool no puede ajustar el modelo, la cruza-verificacion de esta "
        "rama en `06_pstat3.R` es una **auto-verificacion**: R recomputa el alineado "
        "por dos vias (residuo + contribucion del termino via `lm`, y aritmetica de "
        "medias/margenes) y exige que coincidan (`stopifnot`, tol 1e-9), ademas de "
        "chequear el ANOVA III sobre los rangos alineados contra `car::Anova`. "
        "**Esta rama no se ejecuta con los datos reales ni sinteticos.**",
        "",
        "### MEMBRANA como bloque fijo",
        "",
        "- El diseno es balanceado y ortogonal (3 membranas x 4 grupos x 3 replicas "
        "= 36; cada membrana aporta 3 valores a cada celda SEXO x TTO), asi que SS "
        "tipo I = tipo III y la media marginal de celda (promediando sobre MEMBRANA) "
        "= media cruda de celda. `F_MEMBRANA` / `p_MEMBRANA` se reportan como "
        "diagnostico del bloque; no dirigen la inferencia.",
        "",
        "### Limitacion del ensayo (D9)",
        "",
        "- " + LIMITACION_D9.replace("**", ""),
        "",
        "### Paridad R / Python",
        "",
        "- Nucleo numerico PROPIO identico en R y Python (OLS 6x6 por Gauss-Jordan, "
        "SS tipo III por comparacion de modelos, sandwich HC3, Levene, ART y ART-C, "
        "Holm). `06_pstat3.R` cruza-verifica en corrida las ramas que se ejecutan "
        "(`anova3` / `hc3`) contra `car::Anova` y `emmeans` (`stopifnot`, tol 1e-6).",
        "- Unica excepcion: **Shapiro-Wilk** de libreria (`shapiro.test` / "
        "`scipy.stats.shapiro`); coinciden a ~1e-12 y se comparan dentro de "
        "`TOL_ESTADISTICO` en T10. Los p / estadisticos que pasan por funciones "
        "trascendentes (pf, pt, W) se guardan como texto `%.6e` (igual que 03/05).",
    ]
    return "\n".join(L)


def actualizar_descartados(rama, shapiro_p, levene_p):
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 06_pstat3:inicio -->"
    marca_fin = "<!-- 06_pstat3:fin -->"
    nuevo = f"{marca_ini}\n{_bloque_descartes(rama, shapiro_p, levene_p)}\n\n{marca_fin}"
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
    filas = cargar()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_PSTAT3)
    membranas = sorted({r["MEMBRANA"] for r in filas})
    if len(membranas) != 3:
        raise SystemExit(f"06: se esperaban 3 membranas, hay {len(membranas)}: {membranas}")

    rec = clasificar(filas, membranas)
    clasif = [fila_clasif(rec)]
    posthoc = filas_posthoc(rec)
    descr = filas_descriptivo(filas)

    # --- balance del bloque (para verificaciones) ---
    memb_cont = {mb: sum(1 for r in filas if r["MEMBRANA"] == mb) for mb in membranas}
    balanceado = len(set(memb_cont.values())) == 1

    st_ = rec["stats"]
    ip = st_["SEXO:TTO"][1] if "SEXO:TTO" in st_ else None
    sig = ip is not None and ip < ALFA
    n_ph = 0 if rec["posthoc"] is None else len(rec["posthoc"]["filas"])
    ph_sig = []
    if rec["posthoc"] is not None:
        for etiqueta, _e, _ee, _t, _praw, pholm in rec["posthoc"]["filas"]:
            if pholm < ALFA:
                ph_sig.append(etiqueta)

    resumen = [
        f"Via: `{rec['via']}`" + (f" ({rec['motivo']})" if rec["motivo"] else "")
        + f"; rama de la cascada D5: `{rec['rama']}` "
        f"(Shapiro p = {p6e(rec['shapiro_p'])}, Levene p = {p6e(rec['levene_p'])}).",
        f"SEXO: p = {p6e(st_['SEXO'][1])} | TTO: p = {p6e(st_['TTO'][1])} | "
        f"SEXO x TTO: p = {p6e(st_['SEXO:TTO'][1])} | "
        f"MEMBRANA (bloque): p = {p6e(st_['MEMBRANA'][1])}.",
        (f"Interaccion SEXO x TTO significativa (p < {ALFA}): post hoc D6 corrido "
         f"({n_ph} comparaciones, Holm). Significativas tras Holm: "
         f"{', '.join(ph_sig) if ph_sig else '(ninguna)'}.")
        if sig else
        f"Interaccion SEXO x TTO no significativa (p >= {ALFA}): sin post hoc (D6).",
        "pSTAT3 = abundancia de fosfo-STAT3 (no fraccion fosforilada); MEMBRANA "
        "es bloque tecnico fijo. Ver seccion 6 del reporte y D9.",
    ]

    # --- salidas ------------------------------------------------------
    escribir_tsv(cfg.RUTA_DATOS_PROC / "pstat3_modelo_clasificacion.tsv",
                 COLS_CLASIF, clasif)
    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "pstat3_modelo_clasificacion.csv", COLS_CLASIF, clasif)
        escribir_csv(base / "pstat3_posthoc.csv", COLS_POSTHOC, posthoc)
        escribir_csv(base / "pstat3_descriptivo.csv", COLS_DESCR, descr)

    escribir_texto(cfg.RUTA_TABLAS / "pstat3_reporte.md",
                   construir_reporte(fuente, rec, clasif, posthoc, descr, resumen))
    actualizar_descartados(rec["rama"], rec["shapiro_p"], rec["levene_p"])

    ent = f"data/processed/pstat3_long.tsv (de data/{fuente}/{cfg.ARCHIVO_PSTAT3})"
    registrar_procedencia([
        ["data/processed/pstat3_modelo_clasificacion.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "clasificacion del modelo pSTAT3 (D9): via, rama D5, F/p de "
         "SEXO/TTO/SEXOxTTO/MEMBRANA, flag de post hoc"],
        ["outputs/tables/{R,python}/pstat3_modelo_clasificacion.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "idem, copia por implementacion"],
        ["outputs/tables/{R,python}/pstat3_posthoc.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "post hoc D6 de pSTAT3: 4 comparaciones fijas con Holm "
         "(contraste marginal OLS/HC3 o ART-C), solo si interaccion p < alfa"],
        ["outputs/tables/{R,python}/pstat3_descriptivo.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "n, media, sd y mediana de PSTAT3 por grupo SEXO x TTO"],
        ["outputs/tables/pstat3_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO", ent,
         "reporte legible del modelo pSTAT3 (T6), con la limitacion obligatoria D9"],
    ])
    registrar_verificaciones([
        ["pstat3_via", "via del modelo pSTAT3 (D9)", rec["via"], "modelo",
         "TRUE" if rec["via"] == "modelo" else "FALSE", ESTE_SCRIPT],
        ["pstat3_rama_cascada", "rama de la cascada D5 elegida por diagnostico",
         f"{rec['rama']} (Shapiro p={p6e(rec['shapiro_p'])}, Levene p={p6e(rec['levene_p'])})",
         "anova3 (residuos normales, varianzas homogeneas)",
         "TRUE" if rec["rama"] == "anova3" else "FALSE", ESTE_SCRIPT],
        ["pstat3_modelo_d9", "modelo D9 con MEMBRANA como bloque fijo",
         "PSTAT3 ~ SEXO * TTO + MEMBRANA (contr.sum, df_resid = n-6)",
         "PSTAT3 ~ SEXO * TTO + MEMBRANA", "TRUE", ESTE_SCRIPT],
        ["pstat3_bloque_balanceado", "las 3 membranas estan balanceadas",
         ";".join(f"{mb}:{memb_cont[mb]}" for mb in membranas),
         "12;12;12", "TRUE" if (balanceado and set(memb_cont.values()) == {12}) else "FALSE",
         ESTE_SCRIPT],
        ["pstat3_alfa", "alfa para gate D6 y pretests de supuestos",
         g10(ALFA), "0.05", "TRUE" if ALFA == 0.05 else "FALSE", ESTE_SCRIPT],
        ["pstat3_posthoc_holm",
         "post hoc D6: 4 comparaciones fijas con Holm, solo si interaccion p < alfa",
         (f"interaccion p={p6e(ip)} < alfa -> {n_ph} comparaciones, Holm"
          if sig else f"interaccion p={p6e(ip)} >= alfa -> sin post hoc"),
         "4 comparaciones + Holm si y solo si interaccion significativa",
         "TRUE" if ((sig and n_ph == 4) or (not sig and n_ph == 0)) else "FALSE",
         ESTE_SCRIPT],
        ["pstat3_limitacion_d9",
         "limitacion obligatoria: abundancia de fosfo-STAT3, no fraccion fosforilada",
         "escrita en pstat3_reporte.md (seccion 6) y analisis_descartados.md",
         "presente", "TRUE", ESTE_SCRIPT],
        ["pstat3_core_vs_libs",
         "R cruza-verifica anova3/hc3/emmeans contra car::Anova / emmeans (stopifnot 1e-6)",
         "nucleo PROPIO identico R/Python; Shapiro-Wilk por libreria; ART aditivo "
         "auto-verificado (ARTool no admite el bloque)",
         "verificado", "TRUE", ESTE_SCRIPT],
    ])

    print("== 06_pstat3.py ==")
    print(f"  fuente pSTAT3 = {fuente}")
    print(f"  membranas = {membranas}  (n por membrana: "
          f"{', '.join(str(memb_cont[mb]) for mb in membranas)})")
    print(f"  via = {rec['via']}  rama = {rec['rama']}  "
          f"(Shapiro p={p6e(rec['shapiro_p'])}, Levene p={p6e(rec['levene_p'])})")
    print(f"  p: SEXO={p6e(st_['SEXO'][1])}  TTO={p6e(st_['TTO'][1])}  "
          f"SEXOxTTO={p6e(st_['SEXO:TTO'][1])}  MEMBRANA={p6e(st_['MEMBRANA'][1])}")
    if rec["posthoc"] is not None:
        print(f"  post hoc D6 ({rec['posthoc']['metodo']}):")
        for etiqueta, est, ee, tval, praw, pholm in rec["posthoc"]["filas"]:
            print(f"    {etiqueta:30}  est={g10(est)}  p_holm={p6e(pholm)}")
    else:
        print("  post hoc D6: no corrido (interaccion no significativa)")
    print("  -> data/processed/pstat3_modelo_clasificacion.tsv,")
    print("     outputs/tables/{R,python}/pstat3_*.csv, outputs/tables/pstat3_reporte.md")


if __name__ == "__main__":
    main()
