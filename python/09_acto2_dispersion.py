# 09_acto2_dispersion.py -- ACTO 2.3-2.4: dispersion por grupo + test formal de
#                           diferencia de correlaciones placenta <-> cerebro.
#
# Por que existe este archivo: T7 (08) describio la correlacion placenta <->
# cerebro por feto (Spearman rho), global y por TTO, pero tiene PROHIBIDO
# compararlas entre grupos (prohibicion 4) e interpretarlas como coordinacion
# biologica sin descartar restriccion de rango (prohibicion 5). Este script hace
# lo primero de forma formal y prepara lo segundo:
#
#   * 2.3 -- DISPERSION (insumo de la prohibicion 5). Sobre los MISMOS pares por
#     feto que entran en la correlacion, la SD (n-1) de `-ddCt` de cada lado
#     (placenta y cerebro) en Control y en LPS, el cociente de varianzas
#     LPS/Control y un test de Levene (Brown-Forsythe, centro = mediana) por lado.
#     Una caida marcada de la SD bajo LPS es lo que haria sospechar que el cambio
#     de rho es un artefacto de rango -- eso lo dirime la simulacion de 10.
#
#   * 2.4 -- TEST REPORTADO (prohibicion 4; decision del usuario): **Fisher z
#     sobre rho de Spearman**. z = atanh(rho); estadistico
#     (z_control - z_lps) / sqrt(SE_control^2 + SE_lps^2); p a dos colas por la
#     normal. **SE primario = Bonett-Wright** `sqrt((1 + rho^2/2)/(n-3))` (el
#     mismo del IC de T7, coherencia dentro del Acto 2); **SE clasico**
#     `1/sqrt(n-3)` va como columna al lado y NO cambia conclusiones. Se agrega
#     ademas `p_bw_bh` (Benjamini-Hochberg entre los items testeados) como
#     columna suplementaria rotulada, en el espiritu de D12 -- tampoco dirige la
#     inferencia. Piso: se testea solo si Control y LPS tienen n_par >= 5.
#
# Items: los 8 genes con `-ddCt` en ambos tejidos (todos menos il6, D7, y il6R,
# excluido de todo el Acto 2 -- mismo criterio que 08_acto2_correlaciones,
# GEN_EXCLUIDO_CORR) + el score compuesto (D8). Piso: si algun estrato de sexo
# cae por debajo de n_par >= 5, se deja la fila con los n y sin estadistico.
#
# PARIDAD R/Python: rho, SD, cociente y sumas usan acumulador double explicito
# (mismo orden que R) -> texto "%.10g" bit-identico; todo lo que pasa por
# trascendentes (z por atanh, p por la normal, Levene por pf) -> texto "%.6e".
# Las figuras son PNG: equivalentes, no byte-identicas. En R (solo alli) se
# cruza-verifica cada rho contra cor.test(exact=FALSE) y cada F de Levene contra
# car::leveneTest(center=median), con stopifnot.

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from scipy import stats as _sst  # noqa: E402

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "09_acto2_dispersion"

COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}   # Okabe-Ito, igual que 03/07/08
Z975 = 1.959963984540054                             # qnorm(0.975), literal exacto
PISO_PAR = 5                                         # min fetos emparejados por grupo

GEN_SIN_CEREBRO = "il6"                              # D7: sin -ddCt en BRAIN_E15
GEN_EXCLUIDO_CORR = "il6R"   # pedido explicito (08_acto2_correlaciones): deteccion insuficiente en cerebro
GENES_CORR = [g for g in cfg.GENES
              if g not in (GEN_SIN_CEREBRO, GEN_EXCLUIDO_CORR)]   # 8 genes
ITEMS = GENES_CORR + ["score_compuesto"]
TEJIDOS = list(cfg.TEJIDOS_E15)                      # PLACENTA_E15, BRAIN_E15
DPI = 300

# Estratificacion por sexo (pedido explicito, pedidos/cambios_acto2_dispersion_por_sexo.md):
# AMBOS_SEXOS = comportamiento previo (agrupado); HEMBRA/MACHO = dentro de cada
# sexo. Aditivo: AMBOS_SEXOS se conserva tal cual estaba.
ESTRATOS = ["AMBOS_SEXOS", "HEMBRA", "MACHO"]


def _sexo_de_estrato(e):
    return None if e == "AMBOS_SEXOS" else e


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..08.
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


def p6e(x) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return "%.6e" % float(x)


def g10(x) -> str:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return "%.10g" % float(x)


# ---------------------------------------------------------------------------
# Nucleo PROPIO -- acumulador double explicito (mismo orden que R).
# ---------------------------------------------------------------------------
def suma(xs):
    a = 0.0
    for v in xs:
        a += v
    return a


def promedio(xs):
    return suma(xs) / len(xs)


def desvio(xs):
    """SD muestral (n-1), acumulador explicito. NA si n < 2."""
    n = len(xs)
    if n < 2:
        return None
    m = suma(xs) / n
    s = 0.0
    for v in xs:
        s += (v - m) * (v - m)
    return math.sqrt(s / (n - 1))


def mediana(xs):
    if len(xs) == 0:
        return None
    s = sorted(xs)
    k = len(s)
    if k % 2 == 1:
        return s[k // 2]
    return (s[k // 2 - 1] + s[k // 2]) / 2.0


def rangos_promedio(xs):
    """rank() de R con ties.method='average'."""
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
# ANOVA tipo III (contrastes suma-cero) -- nucleo PROPIO, portado de
# 05_qpcr_modelos.py (mismo Gauss-Jordan). Usado SOLO por la Seccion 3 (test
# de interaccion SEXO x TTO sobre la dispersion, NUEVO): a diferencia de D5
# (05/06), aca NO hay cascada de supuestos -- ver esa seccion para la
# justificacion (correccion explicita del usuario).
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


def _fila_diseno(sexo, tto):
    """Diseno suma-cero [1, s, t, s*t]; s=+1 HEMBRA/-1 MACHO, t=+1 CONTROL/-1 LPS."""
    s = 1.0 if sexo == "HEMBRA" else -1.0
    t = 1.0 if tto == "CONTROL" else -1.0
    return [1.0, s, t, s * t]


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


def anova3_terminos(X, y):
    """SS tipo III por comparacion de modelos (con vs sin cada termino)."""
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


def spearman_rho(x, y):
    """rho = Pearson sobre los rangos (corrige empates). Sumas explicitas."""
    n = len(x)
    rx = rangos_promedio(x)
    ry = rangos_promedio(y)
    mx = suma(rx) / n
    my = suma(ry) / n
    sxy = 0.0
    sxx = 0.0
    syy = 0.0
    for i in range(n):
        dx = rx[i] - mx
        dy = ry[i] - my
        sxy += dx * dy
        sxx += dx * dx
        syy += dy * dy
    if sxx <= 0.0 or syy <= 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def f_sf(x, d1, d2):
    return float(_sst.f.sf(x, d1, d2))


def norm_sf2(x):
    """p a dos colas por la normal estandar."""
    return 2.0 * float(_sst.norm.sf(abs(x)))


def levene_bf_2(a, b):
    """Brown-Forsythe de 2 grupos: ANOVA de una via sobre |y - mediana(grupo)|.
    df1 = 1, df2 = nA + nB - 2. Devuelve (F, p) o (None, None) si no aplica."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return None, None
    ma, mb = mediana(a), mediana(b)
    za = [abs(v - ma) for v in a]
    zb = [abs(v - mb) for v in b]
    N = na + nb
    gran = suma(za + zb) / N
    gma = suma(za) / na
    gmb = suma(zb) / nb
    ssb = na * (gma - gran) ** 2 + nb * (gmb - gran) ** 2
    ssw = suma([(v - gma) ** 2 for v in za]) + suma([(v - gmb) ** 2 for v in zb])
    if ssw <= 0.0:
        return None, None
    F = (ssb / 1.0) / (ssw / (N - 2))
    return F, f_sf(F, 1, N - 2)


def bh(pvals):
    """Benjamini-Hochberg -- identico a 05_qpcr_modelos."""
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


# ===========================================================================
# 1. Carga y emparejamiento por feto -- identico a 08.
# ===========================================================================
def _leer_tsv(ruta: Path):
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    h = lineas[0].split("\t")
    return [dict(zip(h, ln.split("\t"))) for ln in lineas[1:]]


def cargar():
    proc = cfg.RUTA_DATOS_PROC
    cuant = _leer_tsv(proc / "qpcr_cuantificacion_long.tsv")
    score = _leer_tsv(proc / "qpcr_score_compuesto_long.tsv")

    madre = {}
    tto = {}
    sexo = {}
    negdd = {}
    for r in cuant:
        f = r["FETO"]
        madre[f] = r["MADRE_ID"]
        tto[f] = r["TTO"]
        sexo[f] = r["SEXO"]
        v = None if r["neg_ddCt"] == "" else float(r["neg_ddCt"])
        negdd[(f, r["TEJIDO"], r["GEN"])] = v
    sc = {}
    for r in score:
        v = None if r["score_compuesto"] == "" else float(r["score_compuesto"])
        sc[(r["FETO"], r["TEJIDO"])] = v

    fetos = sorted(madre.keys(), key=lambda f: (madre[f], f))
    return dict(fetos=fetos, tto=tto, sexo=sexo, negdd=negdd, sc=sc)


def _pares(D, item, estrato, sexo_filtro=None):
    """(placenta[], cerebro[], tto[]) de los fetos con ambos lados detectados.
    Identico a 08: la entrada de T8 son los mismos pares por feto de T7.
    `sexo_filtro`: None = ambos sexos (comportamiento previo); 'HEMBRA'/'MACHO'
    restringe ademas por sexo -- usado por la estratificacion pedida (ver
    analisis_descartados 09_acto2_dispersion): agrupar sexos promedia
    correlaciones de signo opuesto (ver fatcd36 en el reporte)."""
    xs, ys, ts = [], [], []
    for f in D["fetos"]:
        if estrato != "GLOBAL" and D["tto"][f] != estrato:
            continue
        if sexo_filtro is not None and D["sexo"][f] != sexo_filtro:
            continue
        if item == "score_compuesto":
            xp = D["sc"].get((f, "PLACENTA_E15"))
            yb = D["sc"].get((f, "BRAIN_E15"))
        else:
            xp = D["negdd"].get((f, "PLACENTA_E15", item))
            yb = D["negdd"].get((f, "BRAIN_E15", item))
        if xp is None or yb is None:
            continue
        xs.append(xp)
        ys.append(yb)
        ts.append(D["tto"][f])
    return xs, ys, ts


def _lado(xs, ys, tej):
    return xs if tej == "PLACENTA_E15" else ys


# ===========================================================================
# 2.3 -- Tabla de dispersion (insumo de la prohibicion 5). Estratificada por
# ESTRATO (AMBOS_SEXOS/HEMBRA/MACHO, pedido explicito, misma logica que 2.4):
# AMBOS_SEXOS = filas previas, sin cambios; HEMBRA/MACHO son aditivas.
# ===========================================================================
COLS_DISP = ["ITEM", "TIPO", "ESTRATO", "TEJIDO", "n_control", "sd_control",
             "n_lps", "sd_lps", "ratio_var_lps_control", "levene_bf_F",
             "levene_bf_p"]


def tabla_dispersion(D):
    filas = []
    for item in ITEMS:
        tipo = "score" if item == "score_compuesto" else "gen"
        for estrato in ESTRATOS:
            sx = _sexo_de_estrato(estrato)
            cx, cy, _ = _pares(D, item, "CONTROL", sx)
            lx, ly, _ = _pares(D, item, "LPS", sx)
            for tej in TEJIDOS:
                vc = _lado(cx, cy, tej)
                vl = _lado(lx, ly, tej)
                nc, nl = len(vc), len(vl)
                sdc = desvio(vc)
                sdl = desvio(vl)
                ratio = (sdl * sdl) / (sdc * sdc) if (sdc and sdl and sdc > 0) else None
                if nc >= PISO_PAR and nl >= PISO_PAR:
                    F, p = levene_bf_2(vc, vl)
                else:
                    F, p = None, None
                filas.append([item, tipo, estrato, tej, nc, g10(sdc), nl, g10(sdl),
                              g10(ratio), p6e(F), p6e(p)])
    return filas


# ===========================================================================
# 2.4 -- Test reportado: Fisher z sobre rho de Spearman (prohibicion 4).
# Estratificado por ESTRATO (AMBOS_SEXOS/HEMBRA/MACHO, pedido explicito): el
# caso que lo justifica es fatcd36, donde HEMBRA Control rho=0.86 y MACHO
# Control rho=-0.43 -- agrupados dan 0.47, que no describe a ninguno de los
# dos (ver reporte y analisis_descartados). AMBOS_SEXOS = filas previas, sin
# cambios; HEMBRA/MACHO son aditivas.
# ===========================================================================
COLS_TEST = ["ITEM", "TIPO", "ESTRATO", "n_control", "rho_control", "n_lps",
             "rho_lps", "delta_rho", "z_control", "z_lps", "se_bw_control",
             "se_bw_lps", "stat_z_bw", "p_bw", "p_bw_bh", "se_clasico_control",
             "se_clasico_lps", "stat_z_clasico", "p_clasico"]


def _se_bw(rho, n):
    return math.sqrt((1.0 + rho * rho / 2.0) / (n - 3))


def _se_clasico(n):
    return 1.0 / math.sqrt(n - 3)


def _fisher_z(rc, nc, rl, nl, se_fn_c, se_fn_l):
    zc = math.atanh(rc)
    zl = math.atanh(rl)
    sec = se_fn_c
    sel = se_fn_l
    stat = (zc - zl) / math.sqrt(sec * sec + sel * sel)
    return zc, zl, sec, sel, stat, norm_sf2(stat)


# BH (D12) se calcula DENTRO de cada estrato (pedido explicito: la potencia y
# el n difieren mucho entre AMBOS_SEXOS y HEMBRA/MACHO por separado; mezclar
# los 27 p-values en un solo ajuste no tendria sentido).
def tabla_test(D):
    filas = []
    pend_bh = {e: [] for e in ESTRATOS}   # estrato -> [(indice_fila, p_bw), ...]
    for item in ITEMS:
        tipo = "score" if item == "score_compuesto" else "gen"
        for estrato in ESTRATOS:
            sx = _sexo_de_estrato(estrato)
            cx, cy, _ = _pares(D, item, "CONTROL", sx)
            lx, ly, _ = _pares(D, item, "LPS", sx)
            nc, nl = len(cx), len(lx)
            idx = len(filas)
            if nc < PISO_PAR or nl < PISO_PAR:
                filas.append([item, tipo, estrato, nc, "", nl, ""] + [""] * 12)
                continue
            rc = spearman_rho(cx, cy)
            rl = spearman_rho(lx, ly)
            drho = rc - rl
            zc, zl, sbc, sbl, stat_bw, p_bw = _fisher_z(
                rc, nc, rl, nl, _se_bw(rc, nc), _se_bw(rl, nl))
            _zc2, _zl2, scc, scl, stat_cl, p_cl = _fisher_z(
                rc, nc, rl, nl, _se_clasico(nc), _se_clasico(nl))
            pend_bh[estrato].append((idx, p_bw))
            filas.append([item, tipo, estrato, nc, g10(rc), nl, g10(rl), g10(drho),
                          p6e(zc), p6e(zl), p6e(sbc), p6e(sbl), p6e(stat_bw),
                          p6e(p_bw), "", p6e(scc), p6e(scl), p6e(stat_cl), p6e(p_cl)])
    for estrato in ESTRATOS:
        pend = pend_bh[estrato]
        if not pend:
            continue
        ajust = bh([p for _, p in pend])
        for (idx, _), pa in zip(pend, ajust):
            filas[idx][14] = p6e(pa)
    return filas


# ===========================================================================
# 2.5 -- NUEVO: test de interaccion SEXO x TTO sobre la dispersion. Extension
# factorial del Levene Brown-Forsythe de 2.3: en vez de un factor (TTO,
# 2 grupos), un ANOVA III de 2 factores (SEXO, TTO) sobre |x - mediana de su
# celda SEXO x TTO| -- pregunta si el efecto del LPS SOBRE LA VARIABILIDAD
# difiere entre sexos. Aplicado a los (TEJIDO, GEN) que T5 modelo (via ==
# "modelo" en qpcr_modelos_clasificacion.csv; hereda el piso de 5 detectados
# por celda que T5 ya aplico, no se recalcula aca).
#
# CORRECCION DEL USUARIO -- esta es la UNICA parte del proyecto donde D5 (la
# cascada de supuestos) NO se aplica: el propio Levene Brown-Forsythe (lo que
# D5 usaria si fallara Shapiro) YA ES un ANOVA ordinario sobre desvios
# absolutos. Los desvios absolutos son positivos y sesgados por construccion
# -> Shapiro fallaria casi siempre y mandaria todo a ART sin motivo real; y
# evaluar homocedasticidad (Levene) SOBRE una variable que ya es una medida
# de dispersion no tiene sentido logico (es circular). Se usa ANOVA tipo III
# directo con contrastes suma-cero, sin seleccion de rama.
# ===========================================================================
CELDAS_INTER = [("HEMBRA", "CONTROL"), ("HEMBRA", "LPS"),
                ("MACHO", "CONTROL"), ("MACHO", "LPS")]
COLS_INTER = ["GEN", "TEJIDO", "n_HC", "n_HL", "n_MC", "n_ML",
              "mediana_HC", "mad_HC", "mediana_HL", "mad_HL",
              "mediana_MC", "mad_MC", "mediana_ML", "mad_ML",
              "F_SEXO", "p_SEXO", "F_TTO", "p_TTO", "F_SEXOxTTO", "p_SEXOxTTO",
              "p_SEXO_BH", "p_TTO_BH", "p_SEXOxTTO_BH"]


def universo_interaccion():
    """Los (TEJIDO, GEN) con via == 'modelo' en la clasificacion de T5
    (05_qpcr_modelos) -- 18 de 20 (fuera: il6@BRAIN_E15 D7, il6R@BRAIN_E15
    descriptivo_n_bajo). Se lee del propio idioma (RUTA_TABLAS_PY aca)."""
    h, filas = _leer_csv(cfg.RUTA_TABLAS_PY / "qpcr_modelos_clasificacion.csv")
    j_tej, j_gen, j_via = h.index("TEJIDO"), h.index("GEN"), h.index("via")
    return [(f[j_tej], f[j_gen]) for f in filas if f[j_via] == "modelo"]


def _celda_valores(D, tej, gen, sexo, tto):
    """Valores DETECTADOS (-ddCt) de un gen x tejido, en una celda SEXO x TTO
    -- los mismos datos que 04/05 (D['negdd']), no los pares placenta-cerebro
    de 08/09."""
    ys = []
    for f in D["fetos"]:
        if D["sexo"][f] != sexo or D["tto"][f] != tto:
            continue
        v = D["negdd"].get((f, tej, gen))
        if v is None:
            continue
        ys.append(v)
    return ys


def tabla_interaccion(D):
    universo = universo_interaccion()
    filas = []
    pend = {tej: {"SEXO": [], "TTO": [], "SEXO:TTO": []} for tej in TEJIDOS}
    for tej, gen in universo:
        resumen = {}
        m_sexo, m_tto, m_y = [], [], []
        for sexo, tto in CELDAS_INTER:
            ys = _celda_valores(D, tej, gen, sexo, tto)
            med = mediana(ys)
            zs = [abs(v - med) for v in ys]
            resumen[(sexo, tto)] = (len(ys), med, promedio(zs) if zs else None)
            m_sexo += [sexo] * len(zs)
            m_tto += [tto] * len(zs)
            m_y += zs
        X = [_fila_diseno(s, t) for s, t in zip(m_sexo, m_tto)]
        stats, _df = anova3_terminos(X, m_y)
        idx = len(filas)
        for term in ("SEXO", "TTO", "SEXO:TTO"):
            pend[tej][term].append((idx, stats[term][1]))
        r = resumen
        filas.append([
            gen, tej,
            r[CELDAS_INTER[0]][0], r[CELDAS_INTER[1]][0],
            r[CELDAS_INTER[2]][0], r[CELDAS_INTER[3]][0],
            g10(r[CELDAS_INTER[0]][1]), g10(r[CELDAS_INTER[0]][2]),
            g10(r[CELDAS_INTER[1]][1]), g10(r[CELDAS_INTER[1]][2]),
            g10(r[CELDAS_INTER[2]][1]), g10(r[CELDAS_INTER[2]][2]),
            g10(r[CELDAS_INTER[3]][1]), g10(r[CELDAS_INTER[3]][2]),
            p6e(stats["SEXO"][0]), p6e(stats["SEXO"][1]),
            p6e(stats["TTO"][0]), p6e(stats["TTO"][1]),
            p6e(stats["SEXO:TTO"][0]), p6e(stats["SEXO:TTO"][1]),
            "", "", "",
        ])
    # BH dentro de cada tejido (D12, mismo criterio que T5), por termino.
    col = {"SEXO": 20, "TTO": 21, "SEXO:TTO": 22}
    for tej in TEJIDOS:
        for term in ("SEXO", "TTO", "SEXO:TTO"):
            pd = pend[tej][term]
            if not pd:
                continue
            ajust = bh([p for _, p in pd])
            for (idx, _), pa in zip(pd, ajust):
                filas[idx][col[term]] = p6e(pa)
    return filas


def verificar_interaccion_vs_levene(D):
    """Verificacion 7.2 del pedido: el test de interaccion reproduce el
    Levene Brown-Forsythe cuando se colapsa a UN solo factor (TTO, ignorando
    SEXO) -- es la misma matematica (ANOVA tipo III con 1 factor de 2
    niveles == ANOVA de una via == Brown-Forsythe sobre los mismos desvios
    absolutos). Un gen x tejido fijo (determinista), dejado registrado en el
    log."""
    tej, gen = "PLACENTA_E15", "fatcd36"
    vc = _celda_valores(D, tej, gen, "HEMBRA", "CONTROL") + \
        _celda_valores(D, tej, gen, "MACHO", "CONTROL")
    vl = _celda_valores(D, tej, gen, "HEMBRA", "LPS") + \
        _celda_valores(D, tej, gen, "MACHO", "LPS")
    F_lev, p_lev = levene_bf_2(vc, vl)
    mc, ml = mediana(vc), mediana(vl)
    zc = [abs(v - mc) for v in vc]
    zl = [abs(v - ml) for v in vl]
    y = zc + zl
    ttov = ["CONTROL"] * len(zc) + ["LPS"] * len(zl)
    X1 = [[1.0, 1.0 if t == "CONTROL" else -1.0] for t in ttov]
    sse_full, _b, _r = ajustar(X1, y, [0, 1])
    df = len(y) - 2
    mse = sse_full / df
    sse_r, _b2, _r2 = ajustar(X1, y, [0])
    F = (sse_r - sse_full) / mse
    return dict(F_generico=F, p_generico=f_sf(F, 1, df), F_levene=F_lev,
                p_levene=p_lev, tej=tej, gen=gen)


# ===========================================================================
# 3. Figuras.
# ===========================================================================
def figura_dispersion_sd(disp, ruta):
    """SD de -ddCt por item: placenta/cerebro x Control/LPS. Grilla de filas =
    estrato (AMBOS_SEXOS + HEMBRA + MACHO) x columnas = item, para que el
    patron cruzado que describe la seccion 2.3-2.4 (LPS baja la dispersion en
    hembras, la sube en machos) se vea directamente en la figura."""
    idx = {}
    for f in disp:
        idx[(f[0], f[3], f[2])] = f
    ncol = len(ITEMS)
    fig, axes = plt.subplots(len(ESTRATOS), ncol, figsize=(18.0, 9.5))
    for i, estrato in enumerate(ESTRATOS):
        for k, item in enumerate(ITEMS):
            ax = axes[i, k]
            etiquetas = ["pla\nCtrl", "pla\nLPS", "cer\nCtrl", "cer\nLPS"]
            vals, cols = [], []
            for tej in ("PLACENTA_E15", "BRAIN_E15"):
                row = idx[(item, tej, estrato)]
                sdc = float(row[5]) if row[5] != "" else float("nan")
                sdl = float(row[7]) if row[7] != "" else float("nan")
                vals += [sdc, sdl]
                cols += [COL_TTO["CONTROL"], COL_TTO["LPS"]]
            ax.bar(range(4), vals, color=cols, edgecolor="0.3", linewidth=0.5)
            ax.set_xticks(range(4))
            ax.set_xticklabels(etiquetas, fontsize=5.5)
            ax.tick_params(axis="y", labelsize=6.5)
            if i == 0:
                nom = "score compuesto" if item == "score_compuesto" else item
                ax.set_title(nom, fontsize=8,
                             style="normal" if item == "score_compuesto" else "italic")
            if k == 0:
                ax.set_ylabel(f"{estrato}\nSD  -ΔΔCt", fontsize=6.5)
    fig.suptitle("Dispersion de -ΔΔCt dentro del par por feto (placenta y cerebro) "
                 "por grupo\nFilas = estrato de sexo. Una caida marcada de la SD "
                 "bajo LPS es la sospecha de restriccion de rango que dirime la "
                 "simulacion (10, prohibicion 5)", fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


def figura_test_delta_rho(test, ruta):
    """Por item: rho_control vs rho_lps (puntos unidos) + Delta rho y p_bw.
    Los 3 estratos lado a lado (columnas): AMBOS_SEXOS (agrupado, como antes)
    + HEMBRA + MACHO, mismo eje y y mismo rango x en los tres paneles."""
    fig, axes = plt.subplots(1, len(ESTRATOS), figsize=(16.0, 6.5), sharey=True)
    ys = list(range(len(ITEMS)))[::-1]
    for col, estrato in enumerate(ESTRATOS):
        ax = axes[col]
        filas = [f for f in test if f[2] == estrato]
        for y, item, row in zip(ys, ITEMS, filas):
            nom = "score compuesto" if item == "score_compuesto" else item
            if row[4] == "" or row[6] == "":
                ax.text(0.0, y, f"  {nom}: n<{PISO_PAR} en algun grupo (sin test)",
                        va="center", fontsize=6.5, color="0.5")
                continue
            rc, rl = float(row[4]), float(row[6])
            drho, p_bw = float(row[7]), float(row[13])
            ax.plot([rc, rl], [y, y], "-", color="0.7", lw=1.2, zorder=1)
            ax.plot(rc, y, "o", ms=6, color=COL_TTO["CONTROL"], zorder=2,
                    label="Control" if (col == 0 and y == ys[0]) else None)
            ax.plot(rl, y, "o", ms=6, color=COL_TTO["LPS"], zorder=2,
                    label="LPS" if (col == 0 and y == ys[0]) else None)
            estrellas = ("***" if p_bw < 0.001 else "**" if p_bw < 0.01
                         else "*" if p_bw < 0.05 else "n.s.")
            ax.text(1.1, y, f"Δρ={drho:+.2f}  p={p_bw:.3f} {estrellas}",
                    va="center", fontsize=6)
        ax.axvline(0.0, color="0.4", lw=0.8, ls="--")
        ax.set_yticks(ys)
        ax.set_yticklabels(["score compuesto" if i == "score_compuesto" else i
                            for i in ITEMS], fontsize=7.5)
        ax.set_xlim(-1.05, 2.55)
        ax.set_title(estrato, fontsize=9)
        ax.set_xlabel("ρ de Spearman", fontsize=8)
    axes[0].legend(fontsize=8, loc="lower left", frameon=False)
    fig.suptitle("Test reportado (prohibicion 4): Fisher z sobre ρ de Spearman, "
                 "Control vs LPS\nColumnas = estrato de sexo. SE Bonett-Wright; "
                 "p = p_bw con BH entre items (suplementario)", fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


# ===========================================================================
# 4. Artefactos compartidos (merge por 'script') -- headers identicos a 02..08.
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


_DESCARTES = "\n".join([
    "## 09_acto2_dispersion",
    "",
    "### Que resuelve",
    "",
    "- **2.4 -- test formal de diferencia de correlaciones (prohibicion 4).** T7 "
    "tenia prohibido comparar la rho de Control con la de LPS; aca se hace con el "
    "test que **se reporta**: **Fisher z sobre rho de Spearman** (decision del "
    "usuario, entre Fisher z / interaccion de pendientes / permutacion). "
    "`z = atanh(rho)`; estadistico `(z_control - z_lps)/sqrt(SE_control^2 + "
    "SE_lps^2)`; `p` a dos colas por la normal.",
    "- **2.3 -- dispersion por grupo (insumo de la prohibicion 5).** Sobre los "
    "MISMOS pares por feto que la correlacion, la SD (n-1) de `-ddCt` de cada "
    "lado (placenta y cerebro) en Control y en LPS, el cociente de varianzas "
    "LPS/Control y Levene Brown-Forsythe (centro = mediana) por lado. NO alcanza "
    "para descartar restriccion de rango: eso lo hace la simulacion de "
    "`10_acto2_simulacion`.",
    "",
    "### Error estandar del Fisher z (decision del usuario: ambos, BW primario)",
    "",
    "- **Primario = Bonett-Wright**: `SE_i = sqrt((1 + rho_i^2/2)/(n_i - 3))` -- "
    "el mismo SE que ya usa el IC de T7, para que todo el Acto 2 sea coherente. "
    "Corrige levemente por la varianza extra de rho de Spearman respecto de "
    "Pearson.",
    "- **Clasico** `SE_i = 1/sqrt(n_i - 3)` (Fisher de libro) va como columna al "
    "lado (`stat_z_clasico`, `p_clasico`). **No cambia ninguna conclusion**; se "
    "reporta para que el lector vea que la eleccion de SE no mueve el resultado.",
    "- `p_bw_bh`: Benjamini-Hochberg de `p_bw` entre los items testeados. "
    "Columna suplementaria en el espiritu de D12 -- **no dirige la inferencia**; "
    "se menciona en el informe cuantos items sobreviven.",
    "",
    "### Estratificacion por sexo del Delta rho (pedido explicito, aditiva)",
    "",
    "- **Problema**: hasta esta sesion, la tabla de Delta rho (2.4) agrupaba "
    "los sexos (n=16-18 por celda), mientras que desde T7 las figuras de "
    "correlacion ya estan separadas por sexo -- el informe mostraba una cosa "
    "(correlaciones por sexo en las figuras) y testeaba otra (Delta rho "
    "agrupado). El caso que lo deja claro es `fatcd36`: `HEMBRA` Control rho "
    "= 0.86 (n=8) y `MACHO` Control rho = -0.43 (n=8); agrupados dan 0.47, "
    "que no describe a ninguno de los dos.",
    "- **Correccion** (`pedidos/cambios_acto2_dispersion_por_sexo.md`): "
    "columna `ESTRATO` con `AMBOS_SEXOS`/`HEMBRA`/`MACHO` en la tabla de "
    "Delta rho (2.4) y en la de dispersion (2.3, ver abajo). **Aditivo**: el "
    "estrato `AMBOS_SEXOS` es exactamente lo que habia antes, sin cambios; "
    "`HEMBRA`/`MACHO` son filas nuevas.",
    "- **Alternativa descartada**: dejar solo el estrato agrupado. Se "
    "descarta porque promedia correlaciones de signo opuesto (ver `fatcd36` "
    "arriba) -- no es una simplificacion neutral, esconde el patron.",
    "- **Limitacion declarada**: con `n <= 9` por celda de sexo, el Fisher z "
    "tiene potencia baja. Se dice en el cuerpo del reporte (Seccion 1), no en "
    "una nota al pie.",
    "- **BH dentro de cada estrato** (no a traves de los tres): los tres "
    "estratos tienen n y potencia muy distintos (`AMBOS_SEXOS` vs `HEMBRA`/"
    "`MACHO` por separado); mezclar sus 27 `p` en un solo ajuste Benjamini-"
    "Hochberg no tendria sentido estadistico.",
    "",
    "### Dispersion por sexo (2.3), misma estratificacion",
    "",
    "- La tabla de dispersion (`acto2_dispersion.csv`) gana la misma columna "
    "`ESTRATO`: el Levene Brown-Forsythe Control vs LPS por lado se corre "
    "tambien dentro de `HEMBRA` y dentro de `MACHO`, ademas del `AMBOS_SEXOS` "
    "agrupado (sin cambios). Es insumo directo del test de interaccion SEXO x "
    "TTO sobre la dispersion (Seccion 3, nuevo) y de la simulacion "
    "estratificada de `10_acto2_simulacion`.",
    "",
    "### DESCARTADO -- Test de interaccion de pendientes (punto 5 de "
    "`cambios_acto2_dispersion_por_sexo.md`)",
    "",
    "- **Pedido original**: dentro de cada sexo y tejido, `Y ~ X * TTO` sobre "
    "los pares de transportadores (un gen contra otro, ambos en el mismo "
    "tejido), reportando el `p` de la interaccion -- test de si la "
    "PENDIENTE de la relacion entre dos genes cambia con el tratamiento. "
    "Quedo pendiente en la sesion de ese pedido (exigia preguntar antes de "
    "implementar) y se descarta ahora, con el usuario consultado.",
    "- **Por que se descarta (tres razones)**: (1) **contesta una pregunta "
    "distinta a la de este proyecto** -- la relacion entre PARES DE GENES "
    "dentro de un mismo tejido (co-expresion) no es la coordinacion "
    "PLACENTA<->CEREBRO (mismo gen, dos tejidos) que es el objeto del Acto "
    "2; ya existe un control de co-expresion entre transportadores en el "
    "Acto 2.6 (`acto2_coexpresion_transportadores.csv`, eigengene) para eso, "
    "sin necesidad de una interaccion de pendientes por par. (2) **sumaria "
    "muchos tests con `n <= 9` por celda**: 8 genes x 7 pares (C(8,2)=28) x "
    "2 tejidos x 2 sexos = mas de 100 interacciones posibles, cada una con "
    "la misma potencia baja que ya limita al resto del Acto 2 estratificado "
    "por sexo -- inflaria el numero de comparaciones sin agregar potencia. "
    "(3) **la pregunta que si es del proyecto -- si la coordinacion "
    "placenta<->cerebro cambia entre Control y LPS -- ya esta respondida**: "
    "el test formal de Delta rho (2.4, Fisher z, 0/28 significativos) y la "
    "simulacion de restriccion de rango (2.5) cubren esa pregunta sin este "
    "test adicional.",
    "",
    "### NUEVO -- Test de interaccion SEXO x TTO sobre la dispersion (2.5)",
    "",
    "- **Por que hace falta**: estaba en la especificacion original del Acto "
    "2.3 y no se habia implementado. Extension factorial del Levene "
    "Brown-Forsythe: en cerebro E15, siete genes tienen interaccion SEXO x "
    "TTO significativa en el Acto 1 (T5); en tres (`glut1`, `slc38a2`, "
    "`fatp1`) el post hoc explica la interaccion con un patron limpio "
    "(respuesta especifica de hembras), pero en los otros cuatro (`fatcd36`, "
    "`fatp4`, `gp130`, `slc38a1`) la interaccion es significativa y NINGUNA "
    "comparacion puntual sobrevive a Holm -- patron compatible con un efecto "
    "en la DISPERSION, no en la media, que este test puede distinguir y que "
    "el Acto 1 no podia ver.",
    "- **Procedimiento**: por gen x tejido (universo = T5, `via == "
    "\"modelo\"`, 18/20), `z = |x - mediana de la celda SEXO x TTO de x|` "
    "sobre los 36 fetos (valores detectados); ANOVA tipo III de `z ~ SEXO * "
    "TTO`, F y p de `SEXO`, `TTO` y `SEXO:TTO`.",
    "- **CORRECCION DEL USUARIO -- unica parte del proyecto donde D5 NO se "
    "aplica.** Primera version intentaba correr la cascada D5 tambien aca; "
    "el usuario senalo que no corresponde: el Levene Brown-Forsythe (lo que "
    "D5 usaria si fallara Shapiro) YA ES un ANOVA ordinario sobre desvios "
    "absolutos. (a) Los desvios absolutos (`z`) son positivos y sesgados "
    "por construccion -> Shapiro fallaria casi siempre -> cascada mandaria "
    "todo a ART sin motivo real. (b) Evaluar homocedasticidad (Levene) "
    "SOBRE una variable que ya es una medida de dispersion es circular -- "
    "no tiene sentido logico. **Se usa ANOVA tipo III directo con "
    "contrastes suma-cero, sin seleccion de rama, en todos los gen x "
    "tejido de esta seccion.**",
    "- **Nucleo PROPIO** (`resolver`/`ajustar`/`anova3_terminos` y el "
    "diseno suma-cero, portados de `05_qpcr_modelos`, mismo diseno de "
    "4 columnas `[1, s, t, s*t]`). **Verificacion (pedido 7.2)**: "
    "el test colapsado a UN solo factor (`TTO`, ignorando `SEXO`) tiene que "
    "reproducir exactamente el Levene Brown-Forsythe -- son la misma "
    "matematica (ANOVA tipo III con 1 factor de 2 niveles == ANOVA de una "
    "via == Brown-Forsythe sobre los mismos desvios absolutos). Verificado "
    "en corrida con `fatcd36@PLACENTA_E15` (fijo, determinista), `stopifnot` "
    "tol 1e-8, resultado en el log.",
    "- BH (D12) de cada termino, DENTRO de cada tejido (mismo criterio que "
    "T5), suplementario.",
    "",
    "### Alcance y piso",
    "",
    "- Items: 8 genes con `-ddCt` en ambos tejidos (todos menos `il6`, D7, y "
    "`il6R`, excluido de todo el Acto 2 por deteccion insuficiente en "
    "cerebro -- mismo criterio que `08_acto2_correlaciones`, ver bugfix mas "
    "abajo) + score compuesto (D8).",
    "- Se testea un item solo si Control **y** LPS tienen `n_par >= 5`. Si "
    "algun estrato de sexo (`HEMBRA`/`MACHO`) cae por debajo, se deja la "
    "fila con los `n` y sin estadistico, no se omite en silencio.",
    "",
    "### BUGFIX -- `il6R` no se excluia de este script (pedido "
    "`cambios_informe_conclusiones.md`, punto 5b)",
    "",
    "- **Problema**: `08_acto2_correlaciones` excluye `il6R` de todo el Acto 2 "
    "(`GEN_EXCLUIDO_CORR`, deteccion insuficiente en cerebro), pero "
    "`09_acto2_dispersion`, `10_acto2_simulacion` y `11_sensibilidad` "
    "definian `GENES_CORR` cada uno por su cuenta sin esa exclusion -> "
    "`il6R` seguia apareciendo en sus tablas y figuras (Delta rho, "
    "dispersion, simulacion, sensibilidad), aunque el informe ya decia que "
    "estaba excluido.",
    "- **Correccion**: mismo `GEN_EXCLUIDO_CORR <- \"il6R\"` agregado a "
    "`GENES_CORR` en los tres scripts (antes 9 genes, ahora 8 -- consistente "
    "con `08_acto2_correlaciones`). No afecta a la Seccion 3 (test de "
    "interaccion sobre la dispersion, arriba): ese universo sale de la "
    "clasificacion de T5 (`via == \"modelo\"`), no de `GENES_CORR`, y ahi "
    "`il6R@PLACENTA_E15` se sigue modelando (es una pregunta distinta, por "
    "tejido, no del par placenta-cerebro).",
    "",
    "### Paridad R / Python",
    "",
    "- `rho`, `SD`, cociente de varianzas y todas las sumas usan acumulador "
    "`double` explicito (mismo orden) -> texto `%.10g` bit-identico. Lo que pasa "
    "por trascendentes (`z` por `atanh`, `p` por la normal, `F` de Levene por "
    "`pf`) -> texto `%.6e`. Las figuras son PNG: equivalentes, no byte-identicas.",
    "- Solo en R se cruza-verifica cada `rho` contra `cor.test(method=\"spearman\""
    ", exact=FALSE)` (tol 1e-9) y cada `F` de Levene contra "
    "`car::leveneTest(center=median)` (tol 1e-8), con `stopifnot`.",
])


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 09_acto2_dispersion:inicio -->"
    marca_fin = "<!-- 09_acto2_dispersion:fin -->"
    nuevo = f"{marca_ini}\n{_DESCARTES}\n\n{marca_fin}"
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
# 5. Reporte legible.
# ===========================================================================
def _md(header, filas):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for f in filas:
        out.append("| " + " | ".join(_fmt(v) for v in f) + " |")
    return "\n".join(out)


def construir_reporte(fuente, disp, test, inter):
    L = []
    ap = L.append
    ap("# Reporte de dispersion y test de correlaciones Acto 2.3-2.4 (T8)")
    ap("")
    ap("Generado por `09_acto2_dispersion` (R y Python producen este archivo "
       "identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Test reportado: Fisher z sobre rho de Spearman (prohibicion 4)")
    ap("")
    ap("- Compara `rho_control` con `rho_lps` de la correlacion placenta <-> cerebro "
       "por feto. `z = atanh(rho)`, estadistico "
       "`(z_control - z_lps)/sqrt(SE_control^2 + SE_lps^2)`, `p` normal a dos colas.")
    ap("- **SE primario = Bonett-Wright** (coherente con el IC de T7); **SE "
       "clasico** `1/sqrt(n-3)` como columna al lado, no cambia conclusiones. "
       "`p_bw_bh` = BH entre items, suplementario (D12), **calculado dentro "
       "de cada ESTRATO** (no a traves de los tres).")
    ap("- Se testea solo con `n_par >= 5` en Control **y** LPS.")
    ap("- **Estratificado por sexo** (pedido explicito, aditivo): columna "
       "`ESTRATO` = `AMBOS_SEXOS` (agrupado, como antes) / `HEMBRA` / `MACHO`. "
       "**Por que hace falta**: agrupar los sexos promedia correlaciones de "
       "signo opuesto. El caso que lo deja claro es `fatcd36`: dentro de "
       "`HEMBRA` Control la placenta y el cerebro coordinan con signo positivo "
       "fuerte, dentro de `MACHO` Control coordinan con signo NEGATIVO; el "
       "`AMBOS_SEXOS` agrupado da un rho intermedio que no describe a ninguno "
       "de los dos sexos por separado (ver filas `fatcd36` en la tabla, por "
       "`ESTRATO`).")
    ap("- **Potencia**: con `n <= 9` por celda (HEMBRA/MACHO), el Fisher z "
       "tiene muy poca potencia -- un `Delta rho` chico es indetectable y uno "
       "grande puede no alcanzar significancia. Esta limitacion es real y se "
       "declara aca, no en una nota al pie: **no leer \"no significativo en "
       "HEMBRA/MACHO\" como evidencia de que el efecto desaparece al "
       "estratificar** (prohibicion 4 sigue aplicando dentro de cada "
       "estrato).")
    ap("")
    ap(_md(COLS_TEST, test))
    ap("")
    ap("## 2. Dispersion por grupo (insumo de la prohibicion 5)")
    ap("")
    ap("- SD (n-1) de `-ddCt` de cada lado en los mismos pares por feto, cociente "
       "de varianzas LPS/Control y Levene Brown-Forsythe por lado. La lectura "
       "biologica del cambio de correlacion queda pendiente de "
       "`10_acto2_simulacion`.")
    ap("- **Estratificado por sexo** (pedido explicito, aditivo, misma columna "
       "`ESTRATO` que la Seccion 1): el mismo Levene Brown-Forsythe Control vs "
       "LPS, ahora tambien **dentro de cada sexo**. Insumo de la Seccion 3 (test "
       "de interaccion) y de la simulacion estratificada de "
       "`10_acto2_simulacion`.")
    ap("")
    ap(_md(COLS_DISP, disp))
    ap("")
    ap("## 3. Test de interaccion SEXO x TTO sobre la dispersion (NUEVO)")
    ap("")
    ap("- Extension factorial del Levene Brown-Forsythe de la Seccion 2: en vez de "
       "un factor (`TTO`, 2 grupos), un ANOVA tipo III de **2 factores** (`SEXO`, "
       "`TTO`) sobre `z = |x - mediana de su celda SEXO x TTO|` -- pregunta si el "
       "efecto del LPS **sobre la variabilidad** difiere entre sexos. Aplicado a "
       "los (`TEJIDO`, `GEN`) que T5 modelo (`via == \"modelo\"`, 18 de 20; hereda "
       "el piso de 5 detectados por celda que T5 ya aplico).")
    ap("- **Esta es la UNICA parte del proyecto donde D5 (la cascada de supuestos) "
       "NO se aplica** -- correccion explicita del usuario. El propio Levene "
       "Brown-Forsythe (lo que D5 usaria si fallara Shapiro) **ya es** un ANOVA "
       "ordinario sobre desvios absolutos: (a) los desvios absolutos son positivos "
       "y sesgados por construccion, Shapiro fallaria casi siempre y mandaria todo "
       "a ART sin motivo real; (b) evaluar homocedasticidad (Levene) **sobre una "
       "variable que ya es una medida de dispersion** es circular. Se usa **ANOVA "
       "tipo III directo con contrastes suma-cero, sin seleccion de rama**.")
    ap("- `p_*_BH`: Benjamini-Hochberg de cada termino (`SEXO`, `TTO`, "
       "`SEXO:TTO`) **dentro de cada tejido**, mismo criterio que T5 (D12) -- "
       "suplementario, no dirige la inferencia.")
    ap("")
    ap(_md(COLS_INTER, inter))
    ap("")
    ap("## 4. Figuras")
    ap("")
    ap("- `outputs/figures/acto2_dispersion_sd.png` -- SD de -ddCt por item "
       "(placenta/cerebro x Control/LPS).")
    ap("- `outputs/figures/acto2_test_delta_rho.png` -- rho_control vs rho_lps por "
       "item con Delta rho y `p_bw` anotados.")
    ap("")
    ap("## 5. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `09_acto2_dispersion`: eleccion del "
       "test y del SE, rol de la tabla de dispersion, piso de n y items de baja "
       "potencia, y la justificacion completa de por que la Seccion 3 no aplica "
       "D5.")
    ap("")
    return "\n".join(L)


# ===========================================================================
def main():
    D = cargar()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)

    disp = tabla_dispersion(D)
    test = tabla_test(D)

    inter = tabla_interaccion(D)
    chk = verificar_interaccion_vs_levene(D)
    print(f"  [interaccion colapsada a 1 factor vs Levene BF, {chk['gen']}@{chk['tej']}] "
          f"F: {chk['F_generico']:.6f} vs {chk['F_levene']:.6f}  |  "
          f"p: {chk['p_generico']:.6e} vs {chk['p_levene']:.6e}")
    assert abs(chk["F_generico"] - chk["F_levene"]) < 1e-8
    assert abs(chk["p_generico"] - chk["p_levene"]) < 1e-8

    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "acto2_dispersion.csv", COLS_DISP, disp)
        escribir_csv(base / "acto2_test_correlaciones.csv", COLS_TEST, test)
        escribir_csv(base / "acto2_dispersion_interaccion.csv", COLS_INTER, inter)

    fig_sd = cfg.RUTA_FIGURAS / "acto2_dispersion_sd.png"
    fig_dr = cfg.RUTA_FIGURAS / "acto2_test_delta_rho.png"
    figura_dispersion_sd(disp, fig_sd)
    figura_test_delta_rho(test, fig_dr)

    escribir_texto(cfg.RUTA_TABLAS / "acto2_dispersion_reporte.md",
                   construir_reporte(fuente, disp, test, inter))
    actualizar_descartados()

    # Conteos de resumen (print + verificaciones): AMBOS_SEXOS, para no romper
    # el denominador "de N items" de las verificaciones ya existentes;
    # HEMBRA/MACHO se resumen aparte en el print() final.
    def _por_estrato(e):
        return [f for f in test if f[2] == e]

    def _n_test(fs):
        return sum(1 for f in fs if f[4] != "")

    def _n_sig_bw(fs):
        return sum(1 for f in fs if f[13] != "" and float(f[13]) < 0.05)

    def _n_sig_bh(fs):
        return sum(1 for f in fs if f[14] != "" and float(f[14]) < 0.05)

    test_ambos = _por_estrato("AMBOS_SEXOS")
    n_test = _n_test(test_ambos)
    n_sig_bw = _n_sig_bw(test_ambos)
    n_sig_bh = _n_sig_bh(test_ambos)

    # Particion: n(HEMBRA) + n(MACHO) == n(AMBOS_SEXOS), por item y por grupo
    # (verificacion 7.1 del pedido) -- si no cierra, hay error de filtrado.
    n_por = {(f[0], f[2]): (int(f[3]), int(f[5])) for f in test}
    particion_ok = True
    particion_detalle = []
    for item in ITEMS:
        nc_a, nl_a = n_por[(item, "AMBOS_SEXOS")]
        nc_h, nl_h = n_por[(item, "HEMBRA")]
        nc_m, nl_m = n_por[(item, "MACHO")]
        if nc_h + nc_m != nc_a or nl_h + nl_m != nl_a:
            particion_ok = False
            particion_detalle.append(item)

    # Misma particion para la tabla de dispersion (2.3), por item x tejido.
    n_por_disp = {(f[0], f[2], f[3]): (int(f[4]), int(f[6])) for f in disp}
    particion_disp_ok = True
    particion_disp_detalle = []
    for item in ITEMS:
        for tej in TEJIDOS:
            nc_a, nl_a = n_por_disp[(item, "AMBOS_SEXOS", tej)]
            nc_h, nl_h = n_por_disp[(item, "HEMBRA", tej)]
            nc_m, nl_m = n_por_disp[(item, "MACHO", tej)]
            if nc_h + nc_m != nc_a or nl_h + nl_m != nl_a:
                particion_disp_ok = False
                particion_disp_detalle.append(f"{item} {tej}")

    ent = (f"data/processed/qpcr_cuantificacion_long.tsv + qpcr_score_compuesto_long.tsv "
           f"(de data/{fuente}/{cfg.ARCHIVO_QPCR})")
    registrar_procedencia([
        ["outputs/tables/{R,python}/acto2_dispersion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "SD (n-1) de -ddCt por lado (placenta/cerebro) x grupo "
         "sobre los pares por feto; cociente de varianzas LPS/Control; Levene "
         "Brown-Forsythe por lado; estratificado por ESTRATO "
         "(AMBOS_SEXOS/HEMBRA/MACHO, pedido explicito -- mismo motivo que "
         "acto2_test_correlaciones.csv)"],
        ["outputs/tables/{R,python}/acto2_test_correlaciones.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "test reportado (prohibicion 4): Fisher z "
         "sobre rho de Spearman Control vs LPS; SE Bonett-Wright primario + SE "
         "clasico + p_bh suplementario; estratificado por ESTRATO "
         "(AMBOS_SEXOS/HEMBRA/MACHO, pedido explicito -- estratificacion pedida "
         "para ver si el promedio entre sexos tapaba efectos; alternativa "
         "descartada: dejar solo el estrato agrupado, descartada porque promedia "
         "correlaciones de signo opuesto, ver fatcd36; limitacion declarada: "
         "n<=9, potencia baja; BH dentro de cada estrato, no a traves de los tres)"],
        ["outputs/tables/{R,python}/acto2_dispersion_interaccion.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "NUEVO (pedido explicito): ANOVA tipo III "
         "de z=|x-mediana celda SEXO x TTO| ~ SEXO*TTO por gen x tejido "
         "(universo = T5, via==modelo, 18/20); F y p de SEXO/TTO/SEXO:TTO + BH "
         "dentro de tejido; extension factorial del Levene de "
         "acto2_dispersion.csv. UNICA parte del proyecto sin cascada D5 "
         "(correccion del usuario: Levene BF ya es un ANOVA sobre desvios "
         "absolutos, positivos y sesgados por construccion -- Shapiro fallaria "
         "casi siempre y mandaria a ART sin motivo; evaluar homocedasticidad "
         "sobre una medida de dispersion es circular)"],
        ["outputs/figures/acto2_dispersion_sd.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent, "SD de -ddCt por item: placenta/cerebro x Control/LPS"],
        ["outputs/figures/acto2_test_delta_rho.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent, "rho_control vs rho_lps por item con Delta rho y p_bw"],
        ["outputs/tables/acto2_dispersion_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.3-2.4 (T8)"],
    ])
    registrar_verificaciones([
        ["acto2_test_metodo",
         "test reportado de diferencia de correlaciones (prohibicion 4)",
         "Fisher z sobre rho de Spearman; SE Bonett-Wright primario; SE clasico y "
         "p_bh suplementarios",
         "Fisher z (decision del usuario)", "TRUE", ESTE_SCRIPT],
        ["acto2_test_items",
         "items con test formal Control vs LPS (n_par >= 5 en ambos grupos)",
         f"{n_test} de {len(ITEMS)} items testeados",
         "8 genes + score, los que superan el piso", "TRUE", ESTE_SCRIPT],
        ["acto2_test_piso",
         "piso de n_par por grupo para el Fisher z", str(PISO_PAR), "5",
         "TRUE" if PISO_PAR == 5 else "FALSE", ESTE_SCRIPT],
        ["acto2_test_se_doble",
         "se reportan SE Bonett-Wright (primario) y SE clasico 1/sqrt(n-3)",
         "columnas se_bw_* y se_clasico_* + stat/p de cada uno", "ambos SE",
         "TRUE" if COLS_TEST[10] == "se_bw_control" and COLS_TEST[15] ==
         "se_clasico_control" else "FALSE", ESTE_SCRIPT],
        ["acto2_test_bh_suplementario",
         "p_bw_bh (BH DENTRO de cada ESTRATO) es suplementario y no cambia "
         "conclusiones (D12)",
         f"AMBOS_SEXOS: {n_sig_bw} items p_bw<.05; {n_sig_bh} items p_bw_bh<.05",
         "columna rotulada, no dirige", "TRUE", ESTE_SCRIPT],
        ["acto2_estratos_sexo",
         "ESTRATO = AMBOS_SEXOS/HEMBRA/MACHO en Delta rho (2.4) y dispersion (2.3)",
         ";".join(ESTRATOS), "AMBOS_SEXOS;HEMBRA;MACHO",
         "TRUE" if ESTRATOS == ["AMBOS_SEXOS", "HEMBRA", "MACHO"] else "FALSE",
         ESTE_SCRIPT],
        ["acto2_estrato_particion_test",
         "n(HEMBRA) + n(MACHO) = n(AMBOS_SEXOS) por item y por grupo (Delta rho, 2.4)",
         "particiona en %d/%d items%s" % (
             len(ITEMS) - len(particion_detalle), len(ITEMS),
             ("; falla en: " + ", ".join(particion_detalle)) if particion_detalle else ""),
         "particiona en %d/%d items" % (len(ITEMS), len(ITEMS)),
         "TRUE" if particion_ok else "FALSE", ESTE_SCRIPT],
        ["acto2_estrato_particion_dispersion",
         "n(HEMBRA) + n(MACHO) = n(AMBOS_SEXOS) por item x tejido (dispersion, 2.3)",
         "particiona en %d/%d celdas item x tejido%s" % (
             len(ITEMS) * len(TEJIDOS) - len(particion_disp_detalle),
             len(ITEMS) * len(TEJIDOS),
             ("; falla en: " + ", ".join(particion_disp_detalle))
             if particion_disp_detalle else ""),
         "particiona en %d/%d celdas" % (len(ITEMS) * len(TEJIDOS),
                                          len(ITEMS) * len(TEJIDOS)),
         "TRUE" if particion_disp_ok else "FALSE", ESTE_SCRIPT],
        ["acto2_dispersion_pares",
         "la dispersion se mide sobre los MISMOS pares por feto que la correlacion",
         "vector por lado = componente placenta/cerebro de pares(item, grupo)",
         "mismos pares que T7/2.4", "TRUE", ESTE_SCRIPT],
        ["acto2_dispersion_levene",
         "test de dispersion por lado = Levene Brown-Forsythe (centro = mediana)",
         "ANOVA de una via sobre |y - mediana(grupo)|, df1=1 df2=nC+nL-2",
         "Brown-Forsythe", "TRUE", ESTE_SCRIPT],
        ["acto2_dispersion_no_concluye",
         "la tabla de dispersion NO descarta restriccion de rango por si sola",
         "prohibicion 5 la dirime 10_acto2_simulacion; el reporte lo dice",
         "insumo, no conclusion", "TRUE", ESTE_SCRIPT],
        ["acto2_dispersion_figuras",
         "figuras Acto 2.3-2.4: SD por grupo + forest de Delta rho",
         f"sd={fig_sd.is_file()};delta_rho={fig_dr.is_file()}".replace(
             "True", "TRUE").replace("False", "FALSE"),
         "2 figuras existen",
         "TRUE" if (fig_sd.is_file() and fig_dr.is_file()) else "FALSE",
         ESTE_SCRIPT],
        ["acto2_interaccion_universo",
         "el test de interaccion sobre dispersion se aplica al universo modelado por T5",
         f"{len(inter)} filas (T5: {len(inter)} via=='modelo')",
         "18 gen x tejido (20 - il6@BRAIN D7 - il6R@BRAIN descriptivo_n_bajo)",
         "TRUE" if len(inter) == 18 else "FALSE", ESTE_SCRIPT],
        ["acto2_interaccion_sin_cascada",
         "el test de interaccion NO aplica la cascada D5 (correccion explicita "
         "del usuario) -- ANOVA tipo III directo, sin rama",
         "anova3_terminos() unico metodo, sin Shapiro/Levene previos, sin columna rama",
         "ANOVA III directo en las 18 filas", "TRUE", ESTE_SCRIPT],
        ["acto2_interaccion_vs_levene",
         "el test de interaccion colapsado a 1 factor (TTO) reproduce el Levene "
         "Brown-Forsythe (pedido 7.2)",
         "%s@%s: F %.6f vs %.6f; p %.6e vs %.6e; |dif F|=%.2e" % (
             chk["gen"], chk["tej"], chk["F_generico"], chk["F_levene"],
             chk["p_generico"], chk["p_levene"],
             abs(chk["F_generico"] - chk["F_levene"])),
         "|dif F| y |dif p| < 1e-8", "TRUE", ESTE_SCRIPT],
        ["acto2_interaccion_bh",
         "BH (D12) de SEXO/TTO/SEXO:TTO dentro de cada tejido, suplementario",
         "3 columnas p_*_BH, ajuste separado por TEJIDO", "BH por tejido, no dirige",
         "TRUE" if (COLS_INTER[20] == "p_SEXO_BH" and COLS_INTER[21] == "p_TTO_BH"
                    and COLS_INTER[22] == "p_SEXOxTTO_BH") else "FALSE", ESTE_SCRIPT],
    ])

    print("== 09_acto2_dispersion.py ==")
    print(f"  fuente = {fuente}")
    print("  items testeados (n_par>=5 ambos grupos), AMBOS_SEXOS: "
          f"{n_test} / {len(ITEMS)}")
    for f in test_ambos:
        if f[4] != "":
            print(f"    {f[0]:16} rhoC={f[4]:>8}  rhoL={f[6]:>8}  dRho={f[7]:>8}  "
                  f"p_bw={f[13]}  p_bh={f[14]}")
        else:
            print(f"    {f[0]:16} nC={f[3]} nL={f[5]}  (sin test, piso)")
    print(f"  significativos AMBOS_SEXOS: p_bw<.05 -> {n_sig_bw};  "
          f"p_bw_bh<.05 -> {n_sig_bh}")
    for estrato in ("HEMBRA", "MACHO"):
        fs = _por_estrato(estrato)
        print(f"  items testeados {estrato}: {_n_test(fs)} / {len(ITEMS)}  "
              f"(p_bw<.05 -> {_n_sig_bw(fs)}; p_bw_bh<.05 -> {_n_sig_bh(fs)})")
    print("  particion n(HEMBRA)+n(MACHO)=n(AMBOS_SEXOS): " + (
        f"OK en {len(ITEMS)}/{len(ITEMS)} items" if particion_ok
        else "FALLA en: " + ", ".join(particion_detalle)))
    print(f"  -> outputs/tables/{{R,python}}/acto2_dispersion.csv, "
          f"acto2_test_correlaciones.csv")
    print(f"  -> {fig_sd.name}, {fig_dr.name}")

    n_sig_inter = sum(1 for f in inter if f[19] != "" and float(f[19]) < 0.05)
    print(f"  interaccion SEXO x TTO sobre dispersion (2.5, sin cascada D5): "
          f"{n_sig_inter}/{len(inter)} gen x tejido con p_SEXOxTTO < .05")
    for f in inter:
        if f[19] != "" and float(f[19]) < 0.05:
            print(f"    {f[0]:10} {f[1]:14} F_inter={f[18]}  p_inter={f[19]}")
    print("  -> outputs/tables/{R,python}/acto2_dispersion_interaccion.csv")


if __name__ == "__main__":
    main()
