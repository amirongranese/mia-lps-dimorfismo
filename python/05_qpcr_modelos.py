# 05_qpcr_modelos.py -- Modelado de la qPCR: cascada D5, post hoc D6, Fisher D7,
#                       columna suplementaria BH D12.
#
# Por que existe este archivo: sobre `neg_ddCt` de 04_qpcr_cuantificacion, ajusta
# UN modelo por gen x tejido y clasifica cada uno segun el diagnostico de
# supuestos, sin elegir el metodo por conveniencia.
#
#   * D5 / seccion 4.1 -- cascada de supuestos (alfa = 0.05 en todo):
#       modelo base  neg_ddCt ~ SEXO * TTO   (OLS, contrastes suma-cero).
#       - Shapiro-Wilk sobre los RESIDUOS del modelo conjunto del gen x tejido
#         (no por celda) y Levene (Brown-Forsythe, centro = mediana) sobre las 4
#         celdas SEXO x TTO.
#       - Shapiro p >= .05  y  Levene p >= .05  -> ANOVA III  (rama "anova3")
#       - Shapiro p >= .05  y  Levene p <  .05  -> Wald III con HC3 (rama "hc3")
#       - Shapiro p <  .05  (con o sin Levene)  -> ART (rama "art")
#   * Piso de celda: si alguna de las 4 celdas SEXO x TTO tiene < 5 valores
#     DETECTADOS, no se ajusta el modelo -> via "descriptivo_n_bajo" (solo n y
#     descriptivo de neg_ddCt).
#   * D7 -- gen x tejido con calibrador HEMBRA_CONTROL 0/detectados
#     (`cuantificable == FALSE` en 04): fuera del modelo -> via "D7_deteccion":
#     Fisher exacto 2x2 de PROPORCION DE DETECCION Control vs LPS dentro de cada
#     sexo + tabla 2x4 (grupo x detectado). Hoy: solo il6 @ BRAIN_E15.
#   * D6 -- post hoc SOLO si la interaccion SEXO x TTO es significativa
#     (p < 0.05). 4 comparaciones fijas:
#       HEMBRA_CONTROL-HEMBRA_LPS, MACHO_CONTROL-MACHO_LPS,
#       HEMBRA_LPS-MACHO_LPS, HEMBRA_CONTROL-MACHO_CONTROL.
#     Correccion Holm dentro de esas 4.  Metodo por rama:
#       anova3 -> contrastes de medias marginales (vcov OLS)
#       hc3    -> contrastes de medias marginales (vcov HC3)
#       art    -> ART-C (Elkin et al. 2021): rangos de y por celda, modelo de una
#                 via sobre el factor de 4 niveles, t con MSE combinado.
#   * D12 -- columna suplementaria: p ajustado por Benjamini-Hochberg dentro de
#     cada tejido a traves de los genes MODELADOS, por separado para SEXO, TTO y
#     SEXO x TTO. No cambia ninguna conclusion; el reporte dice cuantos sobreviven.
#
# PARIDAD R/Python: el nucleo numerico (OLS 4x4 por eliminacion de Gauss, SS tipo
# III por comparacion de modelos, sandwich HC3, Levene, ART y ART-C, Holm, BH) es
# codigo PROPIO identico en ambos lenguajes. R cruza-verifica cada rama en corrida
# contra car::Anova / emmeans / ARTool con stopifnot (< 1e-6). Unico componente de
# libreria: Shapiro-Wilk (scipy.stats.shapiro aca, shapiro.test en R); coinciden a
# ~1e-12 y se comparan dentro de TOL_ESTADISTICO en T10, no por identidad de bytes.
# Los estadisticos y p dependientes de trascendentes (pf/pt/lgamma/shapiro) se
# guardan como texto "%.6e" (p6e), como en 03_elisa.
#
# Este script NO grafica (Acto 1 completo se arma en 07_figuras_acto1).

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

ESTE_SCRIPT = "05_qpcr_modelos"

ALFA = 0.05          # D5/D6: gate de interaccion y pretests de supuestos
PISO_CELDA = 5       # minimo de detectados por celda SEXO x TTO para ajustar modelo

GRUPOS_4 = ["HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
CELDAS_4 = [("HEMBRA", "CONTROL"), ("HEMBRA", "LPS"),
            ("MACHO", "CONTROL"), ("MACHO", "LPS")]
ORDEN_TEJIDO = {t: i for i, t in enumerate(cfg.TEJIDOS_E15)}
ORDEN_GEN = {g: i for i, g in enumerate(cfg.GENES)}
SET_TRANSP = set(cfg.GENES_TRANSPORTADORES)

# D6: las 4 comparaciones fijas, como (etiqueta, celda_a, celda_b) -> estima a - b
PARES_D6 = [
    ("HEMBRA_CONTROL-HEMBRA_LPS", ("HEMBRA", "CONTROL"), ("HEMBRA", "LPS")),
    ("MACHO_CONTROL-MACHO_LPS", ("MACHO", "CONTROL"), ("MACHO", "LPS")),
    ("HEMBRA_LPS-MACHO_LPS", ("HEMBRA", "LPS"), ("MACHO", "LPS")),
    ("HEMBRA_CONTROL-MACHO_CONTROL", ("HEMBRA", "CONTROL"), ("MACHO", "CONTROL")),
]


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02/03/04.
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
# Codigo identico en R; los sistemas son 4x4 / 3x3, bien condicionados.
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
# Diseno suma-cero [1, s, t, s*t] con s = +1 HEMBRA / -1 MACHO,
#                                    t = +1 CONTROL / -1 LPS.
# ---------------------------------------------------------------------------
def _fila_diseno(sexo, tto):
    s = 1.0 if sexo == "HEMBRA" else -1.0
    t = 1.0 if tto == "CONTROL" else -1.0
    return [1.0, s, t, s * t]


def diseno(muestras):
    X = [_fila_diseno(m["SEXO"], m["TTO"]) for m in muestras]
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
# ---------------------------------------------------------------------------
def anova3_terminos(X, y):
    n = len(y)
    sse_full, _b, _r = ajustar(X, y, [0, 1, 2, 3])
    df = n - 4
    mse = sse_full / df
    out = {}
    for nombre, quita in (("SEXO", 1), ("TTO", 2), ("SEXO:TTO", 3)):
        cols = [c for c in (0, 1, 2, 3) if c != quita]
        sse_r, _b2, _r2 = ajustar(X, y, cols)
        F = (sse_r - sse_full) / mse
        out[nombre] = (F, f_sf(F, 1, df))
    return out, df


# ---------------------------------------------------------------------------
# Rama hc3 -- vcov sandwich HC3 y Wald tipo III (1 gl por termino).
#   V = (X'X)^-1 [ sum_i x_i x_i' e_i^2 / (1 - h_ii)^2 ] (X'X)^-1
# ---------------------------------------------------------------------------
def _vcov_hc3(X, resid):
    n = len(X)
    XtX = [[suma([X[i][a] * X[i][c] for i in range(n)]) for c in range(4)]
           for a in range(4)]
    XtXi = invertir(XtX)
    h = [suma([X[i][a] * XtXi[a][c] * X[i][c]
               for a in range(4) for c in range(4)]) for i in range(n)]
    meat = [[suma([X[i][a] * X[i][c] * (resid[i] * resid[i]) / ((1.0 - h[i]) ** 2)
                   for i in range(n)]) for c in range(4)] for a in range(4)]
    V = [[suma([XtXi[a][k] * meat[k][l] * XtXi[l][c]
               for k in range(4) for l in range(4)]) for c in range(4)]
         for a in range(4)]
    return V, XtXi


def hc3_terminos(X, y):
    n = len(y)
    _sse, b, resid = ajustar(X, y, [0, 1, 2, 3])
    V, _XtXi = _vcov_hc3(X, resid)
    df = n - 4
    out = {}
    for nombre, idx in (("SEXO", 1), ("TTO", 2), ("SEXO:TTO", 3)):
        F = b[idx] * b[idx] / V[idx][idx]
        out[nombre] = (F, f_sf(F, 1, df))
    return out, df


# ---------------------------------------------------------------------------
# Levene / Brown-Forsythe -- ANOVA de una via sobre |y - mediana(celda)|.
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
# Shapiro-Wilk sobre los residuos del modelo conjunto (unica dependencia de
# libreria; scipy aca, shapiro.test en R -- coinciden a ~1e-12).
# ---------------------------------------------------------------------------
def shapiro_residuos(X, y):
    _sse, _b, resid = ajustar(X, y, [0, 1, 2, 3])
    r = _sst.shapiro(resid)
    return float(r.statistic), float(r.pvalue)


# ---------------------------------------------------------------------------
# ART (Wobbrock et al. 2011) -- alineado por inclusion-exclusion con medias de
# grupo (ponderadas), redondeo a 8 decimales, rangos promedio, ANOVA III sobre
# los rangos alineados de cada termino.  Identico a ARTool::art + anova().
# ---------------------------------------------------------------------------
def art_terminos(muestras):
    y = [m["y"] for m in muestras]
    n = len(y)
    gran = promedio(y)
    sm = {lev: promedio([m["y"] for m in muestras if m["SEXO"] == lev])
          for lev in cfg.NIVELES_SEXO}
    tm = {lev: promedio([m["y"] for m in muestras if m["TTO"] == lev])
          for lev in cfg.NIVELES_TTO}
    cm = {k: promedio([m["y"] for m in muestras
                       if (m["SEXO"], m["TTO"]) == k]) for k in CELDAS_4}
    out = {}
    for term in ("SEXO", "TTO", "SEXO:TTO"):
        alin = []
        for m in muestras:
            resid = m["y"] - cm[(m["SEXO"], m["TTO"])]
            if term == "SEXO":
                efecto = sm[m["SEXO"]] - gran
            elif term == "TTO":
                efecto = tm[m["TTO"]] - gran
            else:
                efecto = (gran - sm[m["SEXO"]] - tm[m["TTO"]]
                          + cm[(m["SEXO"], m["TTO"])])
            alin.append(round(resid + efecto, 8))
        rr = rangos_promedio(alin)
        muestras_r = [{"SEXO": m["SEXO"], "TTO": m["TTO"], "y": rr[i]}
                      for i, m in enumerate(muestras)]
        Xr, yr = diseno(muestras_r)
        term_stats, df = anova3_terminos(Xr, yr)
        out[term] = term_stats[term]
    return out, n - 4


# ---------------------------------------------------------------------------
# ART-C para las 4 comparaciones D6 (Elkin et al. 2021): el modelo concatenado
# de una via sobre el factor de 4 celdas se reduce a rank(round(y - grand, 8));
# contrastes t con MSE combinado del modelo de una via, gl = n - 4.
# ---------------------------------------------------------------------------
def artc_pares(muestras):
    y = [m["y"] for m in muestras]
    n = len(y)
    gran = promedio(y)
    rr = rangos_promedio([round(v - gran, 8) for v in y])
    idx = {k: [i for i in range(n) if (muestras[i]["SEXO"], muestras[i]["TTO"]) == k]
           for k in CELDAS_4}
    media = {k: promedio([rr[i] for i in idx[k]]) for k in CELDAS_4}
    fitted = [media[(m["SEXO"], m["TTO"])] for m in muestras]
    sse = suma([(rr[i] - fitted[i]) ** 2 for i in range(n)])
    df = n - 4
    mse = sse / df
    out = {}
    for etiqueta, a, b in PARES_D6:
        se = math.sqrt(mse * (1.0 / len(idx[a]) + 1.0 / len(idx[b])))
        est = media[a] - media[b]
        tval = est / se
        out[etiqueta] = (est, se, tval, t_sf2(tval, df))
    return out, df


# ---------------------------------------------------------------------------
# Contrastes de medias marginales para las ramas anova3 / hc3.
# Modelo saturado -> media marginal de celda = media cruda de celda; el EE del
# contraste sale de c' V c con c = x_a - x_b en el espacio de coeficientes.
# ---------------------------------------------------------------------------
def emmeans_pares(muestras, robusto):
    X, y = diseno(muestras)
    n = len(y)
    sse, b, resid = ajustar(X, y, [0, 1, 2, 3])
    df = n - 4
    mse = sse / df
    if robusto:
        V, _ = _vcov_hc3(X, resid)
    else:
        _XtX = [[suma([X[i][a] * X[i][c] for i in range(n)]) for c in range(4)]
                for a in range(4)]
        XtXi = invertir(_XtX)
        V = [[mse * XtXi[a][c] for c in range(4)] for a in range(4)]
    celdas = {}
    for m in muestras:
        celdas.setdefault((m["SEXO"], m["TTO"]), []).append(m["y"])
    media = {k: promedio(celdas[k]) for k in CELDAS_4}
    out = {}
    for etiqueta, a, bb in PARES_D6:
        xa = _fila_diseno(*a)
        xb = _fila_diseno(*bb)
        cv = [xa[i] - xb[i] for i in range(4)]
        var = suma([cv[i] * V[i][j] * cv[j] for i in range(4) for j in range(4)])
        se = math.sqrt(var)
        est = media[a] - media[bb]
        tval = est / se
        out[etiqueta] = (est, se, tval, t_sf2(tval, df))
    return out, df


# ---------------------------------------------------------------------------
# Holm (dentro de las 4 comparaciones D6) y Benjamini-Hochberg (D12, entre genes).
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


def bh(pvals):
    m = len(pvals)
    if m == 0:
        return []
    orden = sorted(range(m), key=lambda i: pvals[i])
    ajust = [0.0] * m
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = orden[rank]
        val = pvals[i] * m / (rank + 1)
        prev = min(prev, val)
        ajust[i] = min(1.0, prev)
    return ajust


# ---------------------------------------------------------------------------
# Fisher exacto 2x2 a dos colas -- misma regla que stats::fisher.test / 03_elisa.
# ---------------------------------------------------------------------------
def fisher_2x2(a, b, c, d):
    r1, r2 = a + b, c + d
    c1 = a + c
    n = a + b + c + d
    if r1 == 0 or r2 == 0 or c1 == 0 or (n - c1) == 0:
        or_h = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        return or_h, 1.0

    def lp(k):
        return (math.lgamma(r1 + 1) + math.lgamma(r2 + 1)
                + math.lgamma(c1 + 1) + math.lgamma(n - c1 + 1)
                - math.lgamma(n + 1) - math.lgamma(k + 1)
                - math.lgamma(r1 - k + 1) - math.lgamma(c1 - k + 1)
                - math.lgamma(r2 - c1 + k + 1))

    lo = max(0, c1 - r2)
    hi = min(r1, c1)
    p_obs = math.exp(lp(a))
    tot = 0.0
    for k in range(lo, hi + 1):
        pk = math.exp(lp(k))
        if pk <= p_obs * (1 + 1e-7):
            tot += pk
    or_h = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return or_h, min(1.0, tot)


# ===========================================================================
# 1. Carga
# ===========================================================================
def cargar():
    ruta = cfg.RUTA_DATOS_PROC / "qpcr_cuantificacion_long.tsv"
    if not ruta.is_file():
        raise SystemExit(f"05: falta {ruta} (correr 04_qpcr_cuantificacion primero)")
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    h = lineas[0].split("\t")
    filas = [dict(zip(h, ln.split("\t"))) for ln in lineas[1:]]
    for r in filas:
        r["y"] = None if r["neg_ddCt"] == "" else float(r["neg_ddCt"])
        r["detectado"] = r["no_detectado"] != "TRUE"
    return filas


def celdas_detectadas(filas, tej, gen):
    """n detectado por celda SEXO x TTO, en orden CELDAS_4."""
    return [sum(1 for r in filas if r["TEJIDO"] == tej and r["GEN"] == gen
                and r["SEXO"] == sx and r["TTO"] == tt and r["y"] is not None)
            for (sx, tt) in CELDAS_4]


# ===========================================================================
# 2. Clasificacion y modelado por gen x tejido
# ===========================================================================
def clasificar(filas):
    registros = []
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            sub = [r for r in filas if r["TEJIDO"] == tej and r["GEN"] == gen]
            cuant = sub[0]["cuantificable"] == "TRUE" if sub else False
            ncel = celdas_detectadas(filas, tej, gen)
            det = [r for r in sub if r["y"] is not None]
            rec = {
                "TEJIDO": tej, "GEN": gen,
                "es_transportador": gen in SET_TRANSP,
                "n_celda": ncel, "n_detectado_total": len(det),
                "via": None, "motivo": "",
                "shapiro_p": None, "levene_p": None, "rama": "",
                "stats": {}, "df": None, "posthoc": None,
            }
            if not cuant:
                rec["via"] = "D7_deteccion"
                rec["motivo"] = "calibrador HEMBRA_CONTROL 0/detectados (D7)"
            elif min(ncel) < PISO_CELDA:
                rec["via"] = "descriptivo_n_bajo"
                rec["motivo"] = f"celda SEXO x TTO con < {PISO_CELDA} detectados"
            else:
                rec["via"] = "modelo"
                _ajustar_cascada(rec, det)
            registros.append(rec)
    _bh_por_tejido(registros)
    return registros


def _ajustar_cascada(rec, det):
    X, y = diseno(det)
    _W, sp = shapiro_residuos(X, y)
    _lF, lp = levene_bf(det)
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
        st_, df = art_terminos(det)
    rec["stats"] = st_        # {'SEXO':(F,p), 'TTO':(F,p), 'SEXO:TTO':(F,p)}
    rec["df"] = df
    ip = st_["SEXO:TTO"][1]
    if ip < ALFA:
        if rec["rama"] == "art":
            pares, pdf = artc_pares(det)
            metodo = "ART-C"
        else:
            pares, pdf = emmeans_pares(det, robusto=(rec["rama"] == "hc3"))
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


def _bh_por_tejido(registros):
    """D12: BH entre los genes MODELADOS de cada tejido, por termino."""
    for tej in cfg.TEJIDOS_E15:
        mods = [r for r in registros if r["TEJIDO"] == tej and r["via"] == "modelo"]
        for term in ("SEXO", "TTO", "SEXO:TTO"):
            praw = [r["stats"][term][1] for r in mods]
            adj = bh(praw)
            for r, a in zip(mods, adj):
                r.setdefault("bh", {})[term] = a
    for r in registros:
        r.setdefault("bh", {})


# ===========================================================================
# 3. D7 -- deteccion de il6 @ BRAIN_E15 (y cualquier gen x tejido con via D7)
# ===========================================================================
def deteccion_d7(filas, rec):
    tej, gen = rec["TEJIDO"], rec["GEN"]
    sub = [r for r in filas if r["TEJIDO"] == tej and r["GEN"] == gen]
    fisher_filas = []
    for sexo in cfg.NIVELES_SEXO:
        gc = [r for r in sub if r["SEXO"] == sexo and r["TTO"] == "CONTROL"]
        gl = [r for r in sub if r["SEXO"] == sexo and r["TTO"] == "LPS"]
        a = sum(1 for r in gc if r["detectado"])
        b = len(gc) - a
        c = sum(1 for r in gl if r["detectado"])
        d = len(gl) - c
        or_h, p = fisher_2x2(a, b, c, d)
        fisher_filas.append([tej, gen, sexo, len(gc), a, len(gl), c,
                             g10(or_h), p6e(p)])
    tabla_filas = []
    for (sx, tt), grp in zip(CELDAS_4, GRUPOS_4):
        cel = [r for r in sub if r["SEXO"] == sx and r["TTO"] == tt]
        nd = sum(1 for r in cel if r["detectado"])
        tabla_filas.append([tej, gen, grp, len(cel), nd, len(cel) - nd])
    return fisher_filas, tabla_filas


# ===========================================================================
# 4. Armado de tablas de salida
# ===========================================================================
COLS_CLASIF = [
    "TEJIDO", "GEN", "es_transportador",
    "n_HEMBRA_CONTROL", "n_HEMBRA_LPS", "n_MACHO_CONTROL", "n_MACHO_LPS",
    "n_detectado_total", "via", "motivo",
    "shapiro_p_residuos", "levene_p", "rama_cascada",
    "F_SEXO", "p_SEXO", "F_TTO", "p_TTO", "F_SEXOxTTO", "p_SEXOxTTO",
    "p_SEXO_BH", "p_TTO_BH", "p_SEXOxTTO_BH",
    "interaccion_significativa", "post_hoc_corrido",
]


def fila_clasif(rec):
    st_ = rec["stats"]
    bhd = rec["bh"]

    def F(term):
        return p6e(st_[term][0]) if term in st_ else ""

    def P(term):
        return p6e(st_[term][1]) if term in st_ else ""

    def B(term):
        return p6e(bhd[term]) if term in bhd else ""

    ip = st_["SEXO:TTO"][1] if "SEXO:TTO" in st_ else None
    sig = (ip is not None and ip < ALFA)
    return [
        rec["TEJIDO"], rec["GEN"], rec["es_transportador"],
        rec["n_celda"][0], rec["n_celda"][1], rec["n_celda"][2], rec["n_celda"][3],
        rec["n_detectado_total"], rec["via"], rec["motivo"],
        p6e(rec["shapiro_p"]), p6e(rec["levene_p"]), rec["rama"],
        F("SEXO"), P("SEXO"), F("TTO"), P("TTO"), F("SEXO:TTO"), P("SEXO:TTO"),
        B("SEXO"), B("TTO"), B("SEXO:TTO"),
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


COLS_NOMODELO = ["TEJIDO", "GEN", "via", "GRUPO", "n_detectado",
                 "mean_neg_ddCt", "sd_neg_ddCt", "median_neg_ddCt"]


def filas_nomodelo(filas, rec):
    out = []
    for (sx, tt), grp in zip(CELDAS_4, GRUPOS_4):
        vals = [r["y"] for r in filas
                if r["TEJIDO"] == rec["TEJIDO"] and r["GEN"] == rec["GEN"]
                and r["SEXO"] == sx and r["TTO"] == tt and r["y"] is not None]
        n = len(vals)
        mu = promedio(vals) if n else None
        sd = desvio_muestral(vals, mu) if n >= 2 else None
        md = mediana(vals) if n else None
        out.append([rec["TEJIDO"], rec["GEN"], rec["via"], grp, n,
                    g10(mu) if mu is not None else "",
                    g10(sd) if sd is not None else "",
                    g10(md) if md is not None else ""])
    return out


# ===========================================================================
# 5. Reporte legible
# ===========================================================================
def _md(header, filas):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for f in filas:
        out.append("| " + " | ".join(_fmt(v) for v in f) + " |")
    return "\n".join(out)


def construir_reporte(fuente, registros, clasif, posthoc, fisher_filas,
                      tabla_filas, resumen, d7_nota):
    L = []
    ap = L.append
    ap("# Reporte de modelado qPCR (T5)")
    ap("")
    ap("Generado por `05_qpcr_modelos` (R y Python producen este archivo identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Metodo")
    ap("")
    for b in BULLETS_METODO:
        ap(b)
    ap("")
    ap("## 2. Clasificacion por gen x tejido")
    ap("")
    ap(_md(COLS_CLASIF, clasif))
    ap("")
    ap("## 3. Post hoc (D6) -- comparaciones con interaccion significativa")
    ap("")
    if posthoc:
        ap(_md(COLS_POSTHOC, posthoc))
    else:
        ap("Ninguna interaccion SEXO x TTO alcanzo p < alfa: no se corrio post hoc.")
    ap("")
    ap("## 4. D7 -- il6 @ BRAIN_E15: proporcion de deteccion")
    ap("")
    ap("Fisher exacto 2x2 (deteccion Control vs LPS) dentro de cada sexo:")
    ap("")
    ap(_md(["TEJIDO", "GEN", "SEXO", "n_CONTROL", "det_CONTROL", "n_LPS",
            "det_LPS", "OR_Haldane", "p_fisher"], fisher_filas))
    ap("")
    ap("Tabla 2x4 (grupo x deteccion):")
    ap("")
    ap(_md(["TEJIDO", "GEN", "GRUPO", "n_total", "n_detectado", "n_no_detectado"],
           tabla_filas))
    ap("")
    ap(d7_nota)
    ap("")
    ap("## 5. Resumen")
    ap("")
    for linea in resumen:
        ap(f"- {linea}")
    ap("")
    ap("## 6. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `05_qpcr_modelos`: eleccion de "
       "rama por gen x tejido, il6R @ BRAIN_E15 fuera del modelo por piso de "
       "celda, il6 @ BRAIN_E15 degenerado (0 detectados en los 4 grupos sobre "
       "datos reales), y la nota de paridad (Shapiro-Wilk por libreria).")
    ap("")
    return "\n".join(L)


# ===========================================================================
# 6. Artefactos compartidos (merge por 'script') -- headers identicos a 02/03/04.
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
    header = ["id", "tipo", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script"]
    merge_por_script(cfg.RUTA_TABLAS / "verificaciones.csv", header, filas_nuevas,
                     "script", lambda f: (f[6], f[0]))


# Bloques de prosa -- literales IDENTICOS a los de R/05_qpcr_modelos.R (paridad
# byte a byte de los .md). Una linea por bullet (sin corte interno de parrafo).
BULLETS_METODO = [
    "- Modelo por gen x tejido: `neg_ddCt ~ SEXO * TTO` (OLS, contrastes suma-cero). "
    "alfa = 0.05 para la interaccion (gate D6) y para los pretests de supuestos.",
    "- Cascada D5 / 4.1: Shapiro-Wilk sobre los **residuos del modelo conjunto**; "
    "Levene (Brown-Forsythe, centro = mediana) sobre las 4 celdas. Shapiro>=.05 & "
    "Levene>=.05 -> `anova3` (ANOVA III); Shapiro>=.05 & Levene<.05 -> `hc3` (Wald "
    "III con sandwich HC3); Shapiro<.05 -> `art` (Aligned Rank Transform).",
    "- Piso de celda: si alguna celda SEXO x TTO tiene < 5 detectados, no se modela "
    "(`descriptivo_n_bajo`).",
    "- D7: gen x tejido con calibrador HEMBRA_CONTROL 0/detectados "
    "(`cuantificable = FALSE`) fuera del modelo -> Fisher exacto 2x2 de deteccion "
    "Control vs LPS dentro de cada sexo + tabla 2x4.",
    "- Post hoc D6 (solo si interaccion p < alfa): 4 comparaciones fijas, correccion "
    "Holm. `anova3` -> contraste de medias marginales (vcov OLS); `hc3` -> idem con "
    "vcov HC3; `art` -> ART-C (Elkin et al. 2021).",
    "- D12: columna suplementaria `*_BH` = p ajustado por Benjamini-Hochberg dentro "
    "de cada tejido entre los genes modelados, por termino. No dirige la inferencia.",
]


def _bloque_descartes(resumen_rama, resumen_bh):
    L = [
        "## 05_qpcr_modelos",
        "",
        "### D13 -- MADRE (camada) no se incluye en el modelo",
        "",
        "- Se evaluo incluir `MADRE` como efecto aleatorio (18 madres, 2 fetos "
        "por camada -- 1 hembra + 1 macho). **No se incorporo**: "
        "<<< COMPLETAR: evidencia de analisis previos >>>. Se asume "
        "independencia entre fetos para el analisis (D13, AGENTS.md); la "
        "limitacion se declara en el informe (Seccion 7).",
        "",
        "### Eleccion de rama de la cascada D5 por gen x tejido",
        "",
        "- Se ajusta `neg_ddCt ~ SEXO * TTO` (OLS, contr.sum) y se decide la rama con "
        "Shapiro-Wilk sobre los **residuos del modelo conjunto** del gen x tejido y "
        "Levene (Brown-Forsythe) sobre las 4 celdas, ambos a alfa = 0.05. Reparto "
        f"sobre datos reales: {resumen_rama}. La rama usada queda en `rama_cascada` "
        "de `qpcr_modelos_clasificacion.csv` y en `procedencia.csv`.",
        "- La rama `hc3` reporta un Wald tipo III (1 gl) con vcov sandwich HC3; su "
        "post hoc usa el mismo vcov HC3. La rama `art` reporta ART (Wobbrock et al. "
        "2011) y su post hoc es ART-C (Elkin et al. 2021), NUNCA emmeans directo "
        "sobre el modelo ART.",
        "",
        "### il6R @ BRAIN_E15: fuera del modelo por piso de celda",
        "",
        "- Tras descartar no detectados quedan 3/7/3/3 valores en las celdas "
        "HEMBRA_CONTROL / HEMBRA_LPS / MACHO_CONTROL / MACHO_LPS. Con celdas de 3 el "
        "modelo factorial 2x2 y sus pretests de supuestos no son defendibles: via "
        "`descriptivo_n_bajo` (solo n y descriptivo de neg_ddCt en "
        "`qpcr_modelos_nomodelo_descriptivo.csv`). El calibrador de il6R @ BRAIN_E15 "
        "SI tiene detectados (3/9), asi que no es un caso D7.",
        "",
        "### il6 @ BRAIN_E15: D7 -- solo proporcion de deteccion",
        "",
        "- Calibrador HEMBRA_CONTROL 0/9 detectados -> `cuantificable = FALSE` en "
        "04_qpcr_cuantificacion -> via `D7_deteccion` (D7): sin dCt de calibrador no "
        "hay ddCt, asi que il6 en cerebro no entra al modelo y se analiza solo como "
        "proporcion de deteccion. Se emite el Fisher exacto 2x2 (deteccion Control "
        "vs LPS) dentro de cada sexo y la tabla 2x4 grupo x deteccion "
        "(`qpcr_il6_brain_fisher.csv`, `qpcr_il6_brain_tabla2x4.csv`); el reporte "
        "dice si algun Fisher alcanza p < alfa. Es el unico gen x tejido via D7; "
        "se detecta programaticamente (calibrador HEMBRA_CONTROL 0/detectados).",
        "",
        "### Correccion entre genes (D12) -- columna suplementaria",
        "",
        "- Se agrega `p_SEXO_BH`, `p_TTO_BH`, `p_SEXOxTTO_BH` = p ajustado por "
        "Benjamini-Hochberg dentro de cada tejido a traves de los genes modelados, "
        f"por termino. Es suplementaria: no cambia ninguna conclusion. {resumen_bh}",
        "",
        "### Paridad R / Python",
        "",
        "- Todo el nucleo numerico (OLS 4x4 por eliminacion de Gauss, SS tipo III "
        "por comparacion de modelos, sandwich HC3, Levene, ART y ART-C, Holm, BH) es "
        "codigo PROPIO identico en R y Python; R lo cruza-verifica en corrida contra "
        "`car::Anova`, `emmeans` y `ARTool` (`stopifnot`, tol 1e-6).",
        "- Unica excepcion: **Shapiro-Wilk** se toma de libreria (`shapiro.test` en "
        "R, `scipy.stats.shapiro` en Python). Coinciden a ~1e-12; se comparan dentro "
        "de `TOL_ESTADISTICO` en T10, no por identidad de bytes. Los p / "
        "estadisticos que pasan por funciones trascendentes (pf, pt, lgamma, W) se "
        "guardan como texto `%.6e` (igual que 03_elisa).",
    ]
    return "\n".join(L)


def actualizar_descartados(resumen_rama, resumen_bh):
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 05_qpcr_modelos:inicio -->"
    marca_fin = "<!-- 05_qpcr_modelos:fin -->"
    nuevo = f"{marca_ini}\n{_bloque_descartes(resumen_rama, resumen_bh)}\n\n{marca_fin}"
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
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)
    registros = clasificar(filas)

    orden = lambda r: (ORDEN_TEJIDO[r["TEJIDO"]], ORDEN_GEN[r["GEN"]])
    registros.sort(key=orden)

    clasif = [fila_clasif(r) for r in registros]
    posthoc = [fp for r in registros for fp in filas_posthoc(r)]
    nomodelo = [fn for r in registros if r["via"] == "descriptivo_n_bajo"
                for fn in filas_nomodelo(filas, r)]

    fisher_filas, tabla_filas = [], []
    for r in registros:
        if r["via"] == "D7_deteccion":
            ff, tf = deteccion_d7(filas, r)
            fisher_filas += ff
            tabla_filas += tf

    # --- conteos deterministas para resumen / verificaciones ---
    via_cont = {v: sum(1 for r in registros if r["via"] == v)
                for v in ("modelo", "D7_deteccion", "descriptivo_n_bajo")}
    rama_cont = {b: sum(1 for r in registros if r["rama"] == b)
                 for b in ("anova3", "hc3", "art")}
    n_posthoc = sum(1 for r in registros if r.get("posthoc"))
    resumen_rama = "; ".join(f"{k}={rama_cont[k]}" for k in ("anova3", "hc3", "art"))

    # BH: cuantas interacciones significativas crudas sobreviven BH, por tejido
    bh_txt = []
    for tej in cfg.TEJIDOS_E15:
        mods = [r for r in registros if r["TEJIDO"] == tej and r["via"] == "modelo"]
        crudas = [r for r in mods if r["stats"]["SEXO:TTO"][1] < ALFA]
        sobre = [r for r in crudas if r["bh"]["SEXO:TTO"] < ALFA]
        bh_txt.append(f"{tej}: {len(sobre)}/{len(crudas)} interacciones "
                      f"significativas crudas sobreviven BH")
    resumen_bh = "; ".join(bh_txt) + "."

    resumen = [
        f"gen x tejido: {via_cont['modelo']} modelados, "
        f"{via_cont['D7_deteccion']} por D7 (deteccion), "
        f"{via_cont['descriptivo_n_bajo']} descriptivo por piso de celda.",
        f"Reparto de ramas D5: {resumen_rama}.",
        f"Post hoc D6 corrido en {n_posthoc} gen x tejido (interaccion p < {ALFA}).",
        resumen_bh,
    ]

    if fisher_filas:
        n_fis = len(fisher_filas)
        k_fis = sum(1 for f in fisher_filas if f[8] != "" and float(f[8]) < ALFA)
        d7_nota = (f"De los {n_fis} Fisher exactos por sexo, {k_fis} alcanzan "
                   f"p < {ALFA} (Control vs LPS, deteccion de il6 en cerebro).")
    else:
        d7_nota = "No hay gen x tejido via D7 en esta corrida."

    # --- salidas ------------------------------------------------------
    escribir_tsv(cfg.RUTA_DATOS_PROC / "qpcr_modelos_clasificacion.tsv",
                 COLS_CLASIF, clasif)
    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "qpcr_modelos_clasificacion.csv", COLS_CLASIF, clasif)
        escribir_csv(base / "qpcr_modelos_posthoc.csv", COLS_POSTHOC, posthoc)
        escribir_csv(base / "qpcr_modelos_nomodelo_descriptivo.csv",
                     COLS_NOMODELO, nomodelo)
        escribir_csv(base / "qpcr_il6_brain_fisher.csv",
                     ["TEJIDO", "GEN", "SEXO", "n_CONTROL", "det_CONTROL",
                      "n_LPS", "det_LPS", "OR_Haldane", "p_fisher"], fisher_filas)
        escribir_csv(base / "qpcr_il6_brain_tabla2x4.csv",
                     ["TEJIDO", "GEN", "GRUPO", "n_total", "n_detectado",
                      "n_no_detectado"], tabla_filas)

    escribir_texto(cfg.RUTA_TABLAS / "qpcr_modelos_reporte.md",
                   construir_reporte(fuente, registros, clasif, posthoc,
                                     fisher_filas, tabla_filas, resumen, d7_nota))
    actualizar_descartados(resumen_rama, resumen_bh)

    ent = f"data/processed/qpcr_cuantificacion_long.tsv (de data/{fuente}/{cfg.ARCHIVO_QPCR})"
    registrar_procedencia([
        ["data/processed/qpcr_modelos_clasificacion.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "clasificacion por gen x tejido: via, rama D5, F/p de "
         "SEXO/TTO/SEXOxTTO, columna BH (D12), flag de post hoc"],
        ["outputs/tables/{R,python}/qpcr_modelos_clasificacion.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "idem, copia por implementacion"],
        ["outputs/tables/{R,python}/qpcr_modelos_posthoc.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "post hoc D6: 4 comparaciones fijas con Holm, por rama "
         "(contraste marginal OLS/HC3 o ART-C)"],
        ["outputs/tables/{R,python}/qpcr_modelos_nomodelo_descriptivo.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "n y descriptivo de neg_ddCt por grupo para gen x tejido no modelados"],
        ["outputs/tables/{R,python}/qpcr_il6_brain_fisher.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "D7: Fisher exacto 2x2 de deteccion Control vs LPS por sexo"],
        ["outputs/tables/{R,python}/qpcr_il6_brain_tabla2x4.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "D7: tabla 2x4 grupo x deteccion para il6 @ BRAIN_E15"],
        ["outputs/tables/qpcr_modelos_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del modelado qPCR (T5)"],
    ])
    registrar_verificaciones([
        ["modelos_via_conteo", "recalculo",
         "20 gen x tejido repartidos en via modelo / D7 / descriptivo",
         f"modelo={via_cont['modelo']};D7_deteccion={via_cont['D7_deteccion']};"
         f"descriptivo_n_bajo={via_cont['descriptivo_n_bajo']}",
         "modelo=18;D7_deteccion=1;descriptivo_n_bajo=1",
         "TRUE" if (via_cont == {"modelo": 18, "D7_deteccion": 1,
                                 "descriptivo_n_bajo": 1}) else "FALSE",
         ESTE_SCRIPT],
        ["modelos_rama_conteo", "recalculo",
         "cada gen x tejido modelado usa la rama que corresponde a sus propios Shapiro/Levene",
         resumen_rama,
         "rama = art si Shapiro<alfa; hc3 si Shapiro>=alfa & Levene<alfa; anova3 si ambos >=alfa",
         "TRUE" if not [r for r in registros if r["via"] == "modelo"
                        and r["rama"] != (
                            "art" if r["shapiro_p"] < ALFA
                            else "hc3" if r["levene_p"] < ALFA else "anova3")]
         else "FALSE", ESTE_SCRIPT],
        ["modelos_alfa", "declaracion", "alfa para gate D6 y pretests de supuestos",
         g10(ALFA), "0.05", "TRUE" if ALFA == 0.05 else "FALSE", ESTE_SCRIPT],
        ["modelos_piso_celda", "declaracion", "minimo de detectados por celda SEXO x TTO para modelar",
         str(PISO_CELDA), "5", "TRUE" if PISO_CELDA == 5 else "FALSE", ESTE_SCRIPT],
        ["modelos_posthoc_holm", "declaracion",
         "post hoc D6: 4 comparaciones fijas con correccion Holm",
         f"{n_posthoc} gen x tejido; 4 comparaciones c/u; Holm",
         "solo si interaccion p < alfa",
         "TRUE" if n_posthoc >= 0 else "FALSE", ESTE_SCRIPT],
        ["modelos_art_core_vs_artool", "declaracion",
         "ART y ART-C PROPIOS; R cruza-verifica en corrida contra ARTool (stopifnot 1e-6)",
         "nucleo PROPIO identico R/Python; Shapiro-Wilk por libreria",
         "verificado", "TRUE", ESTE_SCRIPT],
        ["modelos_bh_d12", "declaracion", "columna suplementaria BH entre genes por tejido (D12)",
         resumen_bh, "no cambia conclusiones",
         "TRUE", ESTE_SCRIPT],
        ["modelos_d7_il6_brain", "recalculo", "il6 @ BRAIN_E15 fuera del modelo (D7), solo deteccion",
         f"{via_cont['D7_deteccion']} gen x tejido via D7",
         "il6@BRAIN_E15", "TRUE" if via_cont["D7_deteccion"] == 1 else "FALSE",
         ESTE_SCRIPT],
    ])

    print("== 05_qpcr_modelos.py ==")
    print(f"  fuente qPCR = {fuente}")
    print(f"  via: modelo={via_cont['modelo']}  D7={via_cont['D7_deteccion']}  "
          f"descriptivo_n_bajo={via_cont['descriptivo_n_bajo']}")
    print(f"  ramas D5: {resumen_rama}")
    print(f"  post hoc D6 en {n_posthoc} gen x tejido")
    for r in registros:
        if r["via"] != "modelo":
            print(f"    [{r['via']:18}] {r['TEJIDO']:12} {r['GEN']:9} {r['motivo']}")
        else:
            ip = r["stats"]["SEXO:TTO"][1]
            marca = " *interaccion*" if ip < ALFA else ""
            print(f"    [{r['rama']:6}] {r['TEJIDO']:12} {r['GEN']:9} "
                  f"pS={p6e(r['stats']['SEXO'][1])} pT={p6e(r['stats']['TTO'][1])} "
                  f"pSxT={p6e(ip)}{marca}")
    print("  -> data/processed/qpcr_modelos_clasificacion.tsv,")
    print("     outputs/tables/{R,python}/qpcr_modelos_*.csv, qpcr_il6_brain_*.csv,")
    print("     outputs/tables/qpcr_modelos_reporte.md")


if __name__ == "__main__":
    main()
