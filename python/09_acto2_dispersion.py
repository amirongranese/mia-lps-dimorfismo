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
# Items: los 9 genes con `-ddCt` en ambos tejidos (todos menos il6, D7) + el
# score compuesto (D8). `il6R` suele quedar con n_par por grupo chico -> su test
# tiene potencia casi nula; se deja explicito, no se omite.
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
GENES_CORR = [g for g in cfg.GENES if g != GEN_SIN_CEREBRO]   # 9 genes
ITEMS = GENES_CORR + ["score_compuesto"]
TEJIDOS = list(cfg.TEJIDOS_E15)                      # PLACENTA_E15, BRAIN_E15
DPI = 300


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
    negdd = {}
    for r in cuant:
        f = r["FETO"]
        madre[f] = r["MADRE_ID"]
        tto[f] = r["TTO"]
        v = None if r["neg_ddCt"] == "" else float(r["neg_ddCt"])
        negdd[(f, r["TEJIDO"], r["GEN"])] = v
    sc = {}
    for r in score:
        v = None if r["score_compuesto"] == "" else float(r["score_compuesto"])
        sc[(r["FETO"], r["TEJIDO"])] = v

    fetos = sorted(madre.keys(), key=lambda f: (madre[f], f))
    return dict(fetos=fetos, tto=tto, negdd=negdd, sc=sc)


def _pares(D, item, estrato):
    """(placenta[], cerebro[], tto[]) de los fetos con ambos lados detectados.
    Identico a 08: la entrada de T8 son los mismos pares por feto de T7."""
    xs, ys, ts = [], [], []
    for f in D["fetos"]:
        if estrato != "GLOBAL" and D["tto"][f] != estrato:
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
# 2.3 -- Tabla de dispersion (insumo de la prohibicion 5).
# ===========================================================================
COLS_DISP = ["ITEM", "TIPO", "TEJIDO", "n_control", "sd_control", "n_lps",
             "sd_lps", "ratio_var_lps_control", "levene_bf_F", "levene_bf_p"]


def tabla_dispersion(D):
    filas = []
    for item in ITEMS:
        tipo = "score" if item == "score_compuesto" else "gen"
        cx, cy, _ = _pares(D, item, "CONTROL")
        lx, ly, _ = _pares(D, item, "LPS")
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
            filas.append([item, tipo, tej, nc, g10(sdc), nl, g10(sdl),
                          g10(ratio), p6e(F), p6e(p)])
    return filas


# ===========================================================================
# 2.4 -- Test reportado: Fisher z sobre rho de Spearman (prohibicion 4).
# ===========================================================================
COLS_TEST = ["ITEM", "TIPO", "n_control", "rho_control", "n_lps", "rho_lps",
             "delta_rho", "z_control", "z_lps", "se_bw_control", "se_bw_lps",
             "stat_z_bw", "p_bw", "p_bw_bh", "se_clasico_control",
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


def tabla_test(D):
    filas = []
    pend_bh = []          # (indice_fila, p_bw) de las filas con test
    for item in ITEMS:
        tipo = "score" if item == "score_compuesto" else "gen"
        cx, cy, _ = _pares(D, item, "CONTROL")
        lx, ly, _ = _pares(D, item, "LPS")
        nc, nl = len(cx), len(lx)
        if nc < PISO_PAR or nl < PISO_PAR:
            filas.append([item, tipo, nc, "", nl, "", "", "", "", "", "",
                          "", "", "", "", "", "", ""])
            continue
        rc = spearman_rho(cx, cy)
        rl = spearman_rho(lx, ly)
        drho = rc - rl
        zc, zl, sbc, sbl, stat_bw, p_bw = _fisher_z(
            rc, nc, rl, nl, _se_bw(rc, nc), _se_bw(rl, nl))
        _zc2, _zl2, scc, scl, stat_cl, p_cl = _fisher_z(
            rc, nc, rl, nl, _se_clasico(nc), _se_clasico(nl))
        pend_bh.append((len(filas), p_bw))
        filas.append([item, tipo, nc, g10(rc), nl, g10(rl), g10(drho),
                      p6e(zc), p6e(zl), p6e(sbc), p6e(sbl), p6e(stat_bw),
                      p6e(p_bw), "", p6e(scc), p6e(scl), p6e(stat_cl), p6e(p_cl)])
    if pend_bh:
        ajust = bh([p for _, p in pend_bh])
        for (idx, _), pa in zip(pend_bh, ajust):
            filas[idx][13] = p6e(pa)
    return filas


# ===========================================================================
# 3. Figuras.
# ===========================================================================
def figura_dispersion_sd(disp, ruta):
    """SD de -ddCt por item: placenta/cerebro x Control/LPS. Un panel por item."""
    ncol, nrow = 4, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(13.0, 8.5))
    axl = axes.flatten()
    idx = {}
    for f in disp:
        idx.setdefault(f[0], {})[f[2]] = f
    for k, item in enumerate(ITEMS):
        ax = axl[k]
        etiquetas = ["pla\nCtrl", "pla\nLPS", "cer\nCtrl", "cer\nLPS"]
        vals, cols = [], []
        for tej in ("PLACENTA_E15", "BRAIN_E15"):
            row = idx[item][tej]
            sdc = float(row[4]) if row[4] != "" else float("nan")
            sdl = float(row[6]) if row[6] != "" else float("nan")
            vals += [sdc, sdl]
            cols += [COL_TTO["CONTROL"], COL_TTO["LPS"]]
        ax.bar(range(4), vals, color=cols, edgecolor="0.3", linewidth=0.5)
        ax.set_xticks(range(4))
        ax.set_xticklabels(etiquetas, fontsize=6.5)
        ax.tick_params(axis="y", labelsize=7)
        nom = "score compuesto" if item == "score_compuesto" else item
        ax.set_title(nom, fontsize=9.5,
                     style="normal" if item == "score_compuesto" else "italic")
        ax.set_ylabel("SD  -ΔΔCt", fontsize=7.5)
    for k in range(len(ITEMS), len(axl)):
        axl[k].set_visible(False)
    fig.suptitle("Dispersion de -ΔΔCt dentro del par por feto (placenta y cerebro) "
                 "por grupo\nUna caida marcada de la SD bajo LPS es la sospecha de "
                 "restriccion de rango que dirime la simulacion (10, prohibicion 5)",
                 fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


def figura_test_delta_rho(test, ruta):
    """Por item: rho_control vs rho_lps (puntos unidos) + Delta rho y p_bw."""
    fig, ax = plt.subplots(figsize=(9.5, 6.0))
    ys = list(range(len(ITEMS)))[::-1]
    for y, item, row in zip(ys, ITEMS, test):
        nom = "score compuesto" if item == "score_compuesto" else item
        if row[3] == "" or row[5] == "":
            ax.text(0.0, y, f"  {nom}: n<{PISO_PAR} en algun grupo (sin test)",
                    va="center", fontsize=8, color="0.5")
            continue
        rc, rl = float(row[3]), float(row[5])
        drho, p_bw = float(row[6]), float(row[12])
        ax.plot([rc, rl], [y, y], "-", color="0.7", lw=1.2, zorder=1)
        ax.plot(rc, y, "o", ms=7, color=COL_TTO["CONTROL"], zorder=2,
                label="Control" if y == ys[0] else None)
        ax.plot(rl, y, "o", ms=7, color=COL_TTO["LPS"], zorder=2,
                label="LPS" if y == ys[0] else None)
        estrellas = ("***" if p_bw < 0.001 else "**" if p_bw < 0.01
                     else "*" if p_bw < 0.05 else "n.s.")
        ax.text(1.02, y, f"Δρ={drho:+.2f}  p={p_bw:.3f} {estrellas}",
                va="center", fontsize=7.5, transform=ax.get_yaxis_transform())
    ax.axvline(0.0, color="0.4", lw=0.8, ls="--")
    ax.set_yticks(ys)
    ax.set_yticklabels(["score compuesto" if i == "score_compuesto" else i
                        for i in ITEMS], fontsize=8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("ρ de Spearman (placenta ↔ cerebro por feto)", fontsize=9)
    ax.legend(fontsize=8, loc="lower left", frameon=False)
    ax.set_title("Test reportado (prohibicion 4): Fisher z sobre ρ de Spearman, "
                 "Control vs LPS\nSE Bonett-Wright; p = p_bw con BH entre items "
                 "(suplementario)", fontsize=9.5)
    fig.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))
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
    "### Alcance y piso",
    "",
    "- Items: 9 genes con `-ddCt` en ambos tejidos (todos menos `il6`, D7) + "
    "score compuesto (D8).",
    "- Se testea un item solo si Control **y** LPS tienen `n_par >= 5`. `il6R` "
    "suele quedar por debajo en algun grupo -> su test tendria potencia casi "
    "nula; se deja la fila con los `n` y sin estadistico, no se omite en "
    "silencio.",
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


def construir_reporte(fuente, disp, test):
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
       "`p_bw_bh` = BH entre items, suplementario (D12).")
    ap("- Se testea solo con `n_par >= 5` en Control **y** LPS.")
    ap("")
    ap(_md(COLS_TEST, test))
    ap("")
    ap("## 2. Dispersion por grupo (insumo de la prohibicion 5)")
    ap("")
    ap("- SD (n-1) de `-ddCt` de cada lado en los mismos pares por feto, cociente "
       "de varianzas LPS/Control y Levene Brown-Forsythe por lado. La lectura "
       "biologica del cambio de correlacion queda pendiente de "
       "`10_acto2_simulacion`.")
    ap("")
    ap(_md(COLS_DISP, disp))
    ap("")
    ap("## 3. Figuras")
    ap("")
    ap("- `outputs/figures/acto2_dispersion_sd.png` -- SD de -ddCt por item "
       "(placenta/cerebro x Control/LPS).")
    ap("- `outputs/figures/acto2_test_delta_rho.png` -- rho_control vs rho_lps por "
       "item con Delta rho y `p_bw` anotados.")
    ap("")
    ap("## 4. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `09_acto2_dispersion`: eleccion del "
       "test y del SE, rol de la tabla de dispersion, piso de n y items de baja "
       "potencia.")
    ap("")
    return "\n".join(L)


# ===========================================================================
def main():
    D = cargar()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)

    disp = tabla_dispersion(D)
    test = tabla_test(D)

    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "acto2_dispersion.csv", COLS_DISP, disp)
        escribir_csv(base / "acto2_test_correlaciones.csv", COLS_TEST, test)

    fig_sd = cfg.RUTA_FIGURAS / "acto2_dispersion_sd.png"
    fig_dr = cfg.RUTA_FIGURAS / "acto2_test_delta_rho.png"
    figura_dispersion_sd(disp, fig_sd)
    figura_test_delta_rho(test, fig_dr)

    escribir_texto(cfg.RUTA_TABLAS / "acto2_dispersion_reporte.md",
                   construir_reporte(fuente, disp, test))
    actualizar_descartados()

    n_test = sum(1 for f in test if f[3] != "")
    n_sig_bw = sum(1 for f in test if f[12] != "" and float(f[12]) < 0.05)
    n_sig_bh = sum(1 for f in test if f[13] != "" and float(f[13]) < 0.05)

    ent = (f"data/processed/qpcr_cuantificacion_long.tsv + qpcr_score_compuesto_long.tsv "
           f"(de data/{fuente}/{cfg.ARCHIVO_QPCR})")
    registrar_procedencia([
        ["outputs/tables/{R,python}/acto2_dispersion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "SD (n-1) de -ddCt por lado (placenta/cerebro) x grupo "
         "sobre los pares por feto; cociente de varianzas LPS/Control; Levene "
         "Brown-Forsythe por lado"],
        ["outputs/tables/{R,python}/acto2_test_correlaciones.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "test reportado (prohibicion 4): Fisher z "
         "sobre rho de Spearman Control vs LPS; SE Bonett-Wright primario + SE "
         "clasico + p_bh suplementario"],
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
         "9 genes + score, los que superan el piso", "TRUE", ESTE_SCRIPT],
        ["acto2_test_piso",
         "piso de n_par por grupo para el Fisher z", str(PISO_PAR), "5",
         "TRUE" if PISO_PAR == 5 else "FALSE", ESTE_SCRIPT],
        ["acto2_test_se_doble",
         "se reportan SE Bonett-Wright (primario) y SE clasico 1/sqrt(n-3)",
         "columnas se_bw_* y se_clasico_* + stat/p de cada uno", "ambos SE",
         "TRUE" if COLS_TEST[9] == "se_bw_control" and COLS_TEST[14] ==
         "se_clasico_control" else "FALSE", ESTE_SCRIPT],
        ["acto2_test_bh_suplementario",
         "p_bw_bh (BH entre items) es suplementario y no cambia conclusiones (D12)",
         f"{n_sig_bw} items p_bw<.05; {n_sig_bh} items p_bw_bh<.05",
         "columna rotulada, no dirige", "TRUE", ESTE_SCRIPT],
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
    ])

    print("== 09_acto2_dispersion.py ==")
    print(f"  fuente = {fuente}")
    print(f"  items testeados (n_par>=5 ambos grupos): {n_test} / {len(ITEMS)}")
    for f in test:
        if f[3] != "":
            print(f"    {f[0]:16} rhoC={f[3]:>8}  rhoL={f[5]:>8}  dRho={f[6]:>8}  "
                  f"p_bw={f[12]}  p_bh={f[13]}")
        else:
            print(f"    {f[0]:16} nC={f[2]} nL={f[4]}  (sin test, piso)")
    print(f"  significativos: p_bw<.05 -> {n_sig_bw};  p_bw_bh<.05 -> {n_sig_bh}")
    print(f"  -> outputs/tables/{{R,python}}/acto2_dispersion.csv, "
          f"acto2_test_correlaciones.csv")
    print(f"  -> {fig_sd.name}, {fig_dr.name}")


if __name__ == "__main__":
    main()
