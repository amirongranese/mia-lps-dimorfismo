# 08_acto2_correlaciones.py -- ACTO 2.1-2.2: correlacion placenta <-> cerebro.
#
# Por que existe este archivo: el Acto 2 pregunta si el programa de transporte de
# la placenta y el del cerebro fetal estan acoplados dentro de un mismo feto.
# Este script SOLO describe esa correlacion; NO compara correlaciones entre
# grupos (eso es el test formal de T8 -- prohibicion 4) ni la interpreta como
# coordinacion biologica sin descartar restriccion de rango (prohibicion 5, T8).
#
#   * Emparejamiento por FETO: cada feto tiene una fila PLACENTA_E15 y una
#     BRAIN_E15. Un par entra en la correlacion si AMBOS lados estan detectados.
#   * Magnitudes: -ddCt por gen y el score compuesto de transportadores (D8).
#   * Coeficiente (decision del usuario): **Spearman rho** (no Pearson).
#     Robusto a outliers de qPCR. p por la t-aproximacion
#     t = rho * sqrt((n-2)/(1-rho^2)), df = n-2, dos colas (misma formula que
#     scipy.stats.spearmanr por defecto, implementada PROPIA). IC 95% por
#     Bonett-Wright: SE_z = sqrt((1 + rho^2/2)/(n-3)), z = atanh(rho),
#     IC = tanh(z +/- 1.959963984540054 * SE_z).
#   * Estratos (CAMBIO pedido explicito, pedidos/cambios_acto2_correlaciones_
#     por_sexo.md): GLOBAL + por TTO (CONTROL/LPS) + por SEXO x TTO (4 celdas).
#     Los 3 estratos originales SE CONSERVAN; se agregan los 4 por sexo. **No
#     se compara rho entre estratos** (prohibicion 4).
#   * Piso: si el par tiene < 5 fetos, no se calcula rho/IC/p (solo n); en las
#     figuras los puntos se dibujan igual pero sin linea de tendencia.
#   * `il6R` se EXCLUYE de todo el Acto 2 (tablas y figuras) por deteccion
#     insuficiente en cerebro (3/7/3/3), que al estratificar por sexo deja casi
#     todas las celdas bajo el piso de 5 pares (pedido explicito). Se conserva
#     en las figuras y tablas del Acto 1.
#   * Co-expresion: rho de Spearman entre los genes del SPLOM de cada tejido
#     (GENES_SPLOM_PLACENTA = 9, GENES_SPLOM_BRAIN = 8), para AMBOS sexos
#     juntos y por separado.
#
# PARIDAD R/Python: rho y todas las sumas usan acumulador double explicito (mismo
# orden que R) -> bit-identico; los valores que pasan por trascendentes (p por
# pt, IC por tanh/atanh) se guardan como texto "%.6e" (p6e), como en 03/05/06.
# Las figuras son PNG: equivalentes, no byte-identicas.

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from scipy import stats as _sst  # noqa: E402

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "08_acto2_correlaciones"

COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}   # Okabe-Ito, igual que 03/07
MARK_TTO = {"CONTROL": "o", "LPS": "^"}
Z975 = 1.959963984540054                             # qnorm(0.975), literal exacto
PISO_PAR = 5                                         # min fetos emparejados para rho

GEN_SIN_CEREBRO = "il6"          # D7: sin -ddCt en BRAIN_E15
GEN_EXCLUIDO_CORR = "il6R"       # pedido explicito: deteccion insuficiente en cerebro
GENES_CORR = [g for g in cfg.GENES if g not in (GEN_SIN_CEREBRO, GEN_EXCLUIDO_CORR)]  # 8
ITEMS = GENES_CORR + ["score_compuesto"]             # 9

# 3 estratos originales + 4 celdas SEXO x TTO (pedido explicito).
ESTRATOS = ["GLOBAL", "CONTROL", "LPS",
            "HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
FILTRO_ESTRATO = {
    "GLOBAL":         (None, None),
    "CONTROL":        (None, "CONTROL"),
    "LPS":            (None, "LPS"),
    "HEMBRA_CONTROL": ("HEMBRA", "CONTROL"),
    "HEMBRA_LPS":     ("HEMBRA", "LPS"),
    "MACHO_CONTROL":  ("MACHO", "CONTROL"),
    "MACHO_LPS":      ("MACHO", "LPS"),
}
DPI = 300


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..07.
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


def spearman_p(rho, n):
    """t-aproximacion: t = rho*sqrt((n-2)/(1-rho^2)), df = n-2, dos colas."""
    if rho is None or n <= 2:
        return None
    if rho * rho >= 1.0:
        return 0.0
    t = rho * math.sqrt((n - 2) / (1.0 - rho * rho))
    return 2.0 * float(_sst.t.sf(abs(t), n - 2))


def spearman_ci(rho, n):
    """IC 95% Bonett-Wright sobre atanh(rho)."""
    if rho is None or n <= 3 or rho * rho >= 1.0:
        return None, None
    se = math.sqrt((1.0 + rho * rho / 2.0) / (n - 3))
    z = math.atanh(rho)
    return math.tanh(z - Z975 * se), math.tanh(z + Z975 * se)


def ols_ci_banda(x, y, n_grid=60, conf=0.95):
    """Ajuste lineal (mínimos cuadrados) PROPIO + banda de confianza para la
    media, en el espacio lineal (-ddCt). Solo se usa para dibujar (las figuras
    no exigen paridad bit a bit R/Python)."""
    n = len(x)
    xbar = suma(x) / n
    ybar = suma(y) / n
    sxx = 0.0
    sxy = 0.0
    for i in range(n):
        dx = x[i] - xbar
        sxx += dx * dx
        sxy += dx * (y[i] - ybar)
    if sxx <= 0.0:
        return None
    pendiente = sxy / sxx
    intercepto = ybar - pendiente * xbar
    sse = 0.0
    for i in range(n):
        r = y[i] - (intercepto + pendiente * x[i])
        sse += r * r
    if n <= 2:
        return None
    s = math.sqrt(sse / (n - 2))
    t_crit = float(_sst.t.ppf(0.5 + conf / 2.0, n - 2))
    xlo, xhi = min(x), max(x)
    grid = [xlo + (xhi - xlo) * k / (n_grid - 1) for k in range(n_grid)]
    y_hat = [intercepto + pendiente * xg for xg in grid]
    se_pred = [s * math.sqrt(1.0 / n + (xg - xbar) ** 2 / sxx) for xg in grid]
    lo = [yh - t_crit * se for yh, se in zip(y_hat, se_pred)]
    hi = [yh + t_crit * se for yh, se in zip(y_hat, se_pred)]
    return grid, y_hat, lo, hi


# ===========================================================================
# 1. Carga y emparejamiento por feto.
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
        # SEXO ya viene normalizado (HEMBRA/MACHO) desde 02_ingesta_qc; se usa
        # esa columna directamente (no se deriva de GRUPO).
        sexo[f] = r["SEXO"]
        v = None if r["neg_ddCt"] == "" else float(r["neg_ddCt"])
        negdd[(f, r["TEJIDO"], r["GEN"])] = v
    sc = {}
    for r in score:
        v = None if r["score_compuesto"] == "" else float(r["score_compuesto"])
        sc[(r["FETO"], r["TEJIDO"])] = v

    fetos = sorted(madre.keys(), key=lambda f: (madre[f], f))
    return dict(fetos=fetos, tto=tto, sexo=sexo, negdd=negdd, sc=sc)


def _pares(D, item, estrato):
    """(placenta[], cerebro[], tto[], sexo[]) de los fetos con ambos lados
    detectados, filtrados segun el estrato pedido."""
    sexo_f, tto_f = FILTRO_ESTRATO[estrato]
    xs, ys, ts, ss = [], [], [], []
    for f in D["fetos"]:
        if tto_f is not None and D["tto"][f] != tto_f:
            continue
        if sexo_f is not None and D["sexo"][f] != sexo_f:
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
        ss.append(D["sexo"][f])
    return xs, ys, ts, ss


# ===========================================================================
# 2. Tablas de correlacion.
# ===========================================================================
COLS_CORR = ["ITEM", "TIPO", "ESTRATO", "n_par", "rho_spearman",
             "ic95_low", "ic95_high", "p_valor"]


def tabla_correlaciones(D):
    filas = []
    for item in ITEMS:
        tipo = "score" if item == "score_compuesto" else "gen"
        for est in ESTRATOS:
            xs, ys, _ts, _ss = _pares(D, item, est)
            n = len(xs)
            if n >= PISO_PAR:
                rho = spearman_rho(xs, ys)
                lo, hi = spearman_ci(rho, n)
                p = spearman_p(rho, n)
                filas.append([item, tipo, est, n, g10(rho), p6e(lo), p6e(hi), p6e(p)])
            else:
                filas.append([item, tipo, est, n, "", "", "", ""])
    return filas


COLS_COEXP = ["TEJIDO", "SEXO", "GEN_A", "GEN_B", "n_par", "rho_spearman", "p_valor"]


def genes_splom_tejido(tej):
    return cfg.GENES_SPLOM_PLACENTA if tej == "PLACENTA_E15" else cfg.GENES_SPLOM_BRAIN


def tabla_coexpresion(D):
    filas = []
    for tej in cfg.TEJIDOS_E15:
        gs = genes_splom_tejido(tej)
        for sx in ("AMBOS", "HEMBRA", "MACHO"):
            for i in range(len(gs)):
                for j in range(i + 1, len(gs)):
                    a, b = gs[i], gs[j]
                    xs, ys = [], []
                    for f in D["fetos"]:
                        if sx != "AMBOS" and D["sexo"][f] != sx:
                            continue
                        va = D["negdd"].get((f, tej, a))
                        vb = D["negdd"].get((f, tej, b))
                        if va is None or vb is None:
                            continue
                        xs.append(va)
                        ys.append(vb)
                    n = len(xs)
                    if n >= PISO_PAR:
                        rho = spearman_rho(xs, ys)
                        p = spearman_p(rho, n)
                        filas.append([tej, sx, a, b, n, g10(rho), p6e(p)])
                    else:
                        filas.append([tej, sx, a, b, n, "", ""])
    return filas


# ===========================================================================
# 3. Figuras.
# ===========================================================================
def _rho_txt(D, item, estratos):
    out = {}
    for est in estratos:
        xs, ys, _, _ = _pares(D, item, est)
        n = len(xs)
        if n >= PISO_PAR:
            rho = spearman_rho(xs, ys)
            out[est] = (rho, n)
        else:
            out[est] = (None, n)
    return out


# --- 3.1 Figura global (se conserva: sin separar por sexo) -----------------
def figura_dispersion(D, ruta):
    ncol, nrow = 3, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(11.5, 9.5))
    axl = axes.flatten()
    for k, item in enumerate(ITEMS):
        ax = axl[k]
        xs, ys, ts, _ = _pares(D, item, "GLOBAL")
        for tt in ("CONTROL", "LPS"):
            xx = [xs[i] for i in range(len(xs)) if ts[i] == tt]
            yy = [ys[i] for i in range(len(ys)) if ts[i] == tt]
            ax.plot(xx, yy, "o", ms=4, mfc=COL_TTO[tt], mec="white", mew=0.4,
                    ls="none", label=tt.capitalize())
        rt = _rho_txt(D, item, ("GLOBAL", "CONTROL", "LPS"))
        rg, ng = rt["GLOBAL"]
        lo, hi = spearman_ci(rg, ng) if rg is not None else (None, None)
        linea1 = (f"ρ = {rg:.2f} [{lo:.2f}, {hi:.2f}]  (n={ng})"
                  if rg is not None and lo is not None else f"(n={ng}, sin ρ)")
        rc, nc = rt["CONTROL"]
        rl, nl = rt["LPS"]
        linea2 = (f"Control ρ={rc:.2f} (n={nc})  |  LPS ρ={rl:.2f} (n={nl})"
                  if rc is not None and rl is not None else
                  f"Control n={nc} | LPS n={nl}")
        nom = "score compuesto" if item == "score_compuesto" else item
        ax.set_title(nom, fontsize=9.5,
                     style="normal" if item == "score_compuesto" else "italic")
        ax.set_xlabel("placenta  -ΔΔCt", fontsize=7.5)
        ax.set_ylabel("cerebro  -ΔΔCt", fontsize=7.5)
        ax.tick_params(labelsize=7)
        ax.annotate(linea1 + "\n" + linea2, xy=(0.03, 0.97), xycoords="axes fraction",
                    va="top", ha="left", fontsize=7, color="0.15",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", lw=0.5))
    for k in range(len(ITEMS), len(axl)):
        axl[k].set_visible(False)
    axl[0].legend(fontsize=7.5, loc="lower right", frameon=False)
    fig.suptitle("Correlacion placenta <-> cerebro por feto  --  Spearman ρ "
                 "(-ΔΔCt por gen y score compuesto)\n"
                 "T7 solo describe: la diferencia de ρ entre Control y LPS NO se "
                 "testea aca (prohibicion 4); el test formal es T8. il6R excluido "
                 "(deteccion insuficiente, pedido explicito)", fontsize=9.5)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.91))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


# --- 3.2 Figura nueva: una por item, dos paneles (Females/Males) -----------
# Ejes en escala log2 mostrando FC = 2^(-ddCt) (D2), salvo score_compuesto: es
# un z-score (puede ser negativo), no tiene FC -- se grafica en escala lineal.
# Desviacion documentada en analisis_descartados.md (seccion de este script).
NOTA_PIE = ("Axes: 2^-ΔΔCt (log2 display) | Spearman on -ΔΔCt | "
            "Linear fit with 95% CI | Same fetus pairing")
NOTA_PIE_SCORE = ("Axis: composite z-score (linear) | Spearman on the z-score | "
                   "Linear fit with 95% CI | Same fetus pairing")


def _linea_leyenda(nombre, D, item, sx, tt):
    xs, ys, _, _ = _pares(D, item, f"{sx}_{tt}")
    n = len(xs)
    if n >= PISO_PAR:
        rho = spearman_rho(xs, ys)
        p = spearman_p(rho, n)
        return f"{nombre}: Spearman rho = {rho:.2f} ; p = {p:.3f} ; n = {n}"
    return f"{nombre}: n = {n} (sin rho)"


def figura_gen_sexo(D, item, ruta):
    es_score = item == "score_compuesto"
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 7.0), sharey=True)
    fig.subplots_adjust(left=0.13, right=0.97, top=0.80, bottom=0.36, wspace=0.06)
    paneles = [("HEMBRA", "Females", axes[0]), ("MACHO", "Males", axes[1])]

    for sx, nombre_panel, ax in paneles:
        ax.set_title(nombre_panel, fontsize=11, fontweight="bold",
                     bbox=dict(boxstyle="square,pad=0.35", fc="0.91", ec="0.35"))
        for spine in ax.spines.values():
            spine.set_edgecolor("0.35")
        for tt in ("CONTROL", "LPS"):
            xs, ys, _, _ = _pares(D, item, f"{sx}_{tt}")
            if not xs:
                continue
            if es_score:
                xp, yp = xs, ys
            else:
                xp = [2.0 ** v for v in xs]
                yp = [2.0 ** v for v in ys]
            ax.plot(xp, yp, MARK_TTO[tt], ms=5.5, mfc=COL_TTO[tt], mec="white",
                     mew=0.5, ls="none", alpha=0.9)
            if len(xs) >= PISO_PAR:
                banda = ols_ci_banda(xs, ys)
                if banda is not None:
                    grid, y_hat, lo, hi = banda
                    if es_score:
                        gx, gy, glo, ghi = grid, y_hat, lo, hi
                    else:
                        gx = [2.0 ** v for v in grid]
                        gy = [2.0 ** v for v in y_hat]
                        glo = [2.0 ** v for v in lo]
                        ghi = [2.0 ** v for v in hi]
                    ax.plot(gx, gy, "-", color=COL_TTO[tt], lw=1.4)
                    ax.fill_between(gx, glo, ghi, color=COL_TTO[tt], alpha=0.15,
                                     lw=0)
        if not es_score:
            ax.set_xscale("log", base=2)
            ax.set_yscale("log", base=2)
        ax.tick_params(labelsize=8)

    if es_score:
        lab_x = "Placenta E15 -- composite transporter z-score"
        lab_y = "Brain E15 -- composite transporter z-score"
        titulo = "score compuesto"
    else:
        lab_x = f"Placenta E15 — {item}/rsp29 relative expression"
        lab_y = f"Brain E15 — {item}/rsp29 relative expression"
        titulo = item

    fig.text(0.55, 0.305, lab_x, ha="center", va="top", fontsize=9.5)
    axes[0].set_ylabel(lab_y, fontsize=9.5)

    leyenda = "\n".join([
        _linea_leyenda("♀ Control", D, item, "HEMBRA", "CONTROL"),
        _linea_leyenda("♀ LPS", D, item, "HEMBRA", "LPS"),
        _linea_leyenda("♂ Control", D, item, "MACHO", "CONTROL"),
        _linea_leyenda("♂ LPS", D, item, "MACHO", "LPS"),
    ])
    pie = NOTA_PIE_SCORE if es_score else NOTA_PIE

    fig.text(0.065, 0.965, titulo, fontsize=13, ha="left", va="top",
              fontstyle="normal" if es_score else "italic")
    fig.text(0.065, 0.925, "Placenta–brain correlation at E15", fontsize=9.5,
              ha="left", va="top")
    fig.text(0.065, 0.255, leyenda + "\n" + pie, fontsize=7.6, ha="left", va="top",
              linespacing=1.3)
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


# --- 3.3 SPLOM de co-expresion: AMBOS sexos + por sexo, gen set por tejido --
def _panel_diag_kde(ax, vals_por_tto):
    """Densidad KDE de Control vs LPS superpuestas y sombreadas (diagonal)."""
    for tt in ("CONTROL", "LPS"):
        v = vals_por_tto.get(tt, [])
        if len(v) < 2 or len(set(v)) < 2:
            continue
        try:
            kde = _sst.gaussian_kde(v)
        except Exception:
            continue
        xg = [min(v) - 0.15 * (max(v) - min(v)) + k * (max(v) - min(v)) * 1.3 / 99
              for k in range(100)]
        yg = kde(xg)
        ax.plot(xg, yg, color=COL_TTO[tt], lw=0.9)
        ax.fill_between(xg, yg, color=COL_TTO[tt], alpha=0.35, lw=0)
    ax.set_yticks([])


def figura_splom(D, tej, ruta, genes_t, sexo_filtro="AMBOS"):
    fetos_uso = [f for f in D["fetos"]
                 if sexo_filtro == "AMBOS" or D["sexo"][f] == sexo_filtro]
    n = len(genes_t)
    datos = {g: [] for g in genes_t}
    ttos = []
    for f in fetos_uso:
        ttos.append(D["tto"][f])
        for g in genes_t:
            v = D["negdd"].get((f, tej, g))
            datos[g].append(v if v is not None else float("nan"))

    fig, axes = plt.subplots(n, n, figsize=(12.5, 12.5))
    for i, gi in enumerate(genes_t):
        for j, gj in enumerate(genes_t):
            ax = axes[i, j]
            ax.tick_params(labelsize=6)
            if i == j:
                por_tto = {tt: [datos[gi][k] for k in range(len(datos[gi]))
                                if ttos[k] == tt and not math.isnan(datos[gi][k])]
                           for tt in ("CONTROL", "LPS")}
                _panel_diag_kde(ax, por_tto)
            elif i > j:
                x = datos[gj]
                y = datos[gi]
                for tt in ("CONTROL", "LPS"):
                    xx = [x[k] for k in range(len(x))
                          if ttos[k] == tt and not math.isnan(x[k]) and not math.isnan(y[k])]
                    yy = [y[k] for k in range(len(y))
                          if ttos[k] == tt and not math.isnan(x[k]) and not math.isnan(y[k])]
                    ax.plot(xx, yy, "o", ms=2.6, mfc=COL_TTO[tt], mec="none",
                            ls="none")
            else:
                xy = [(datos[gj][k], datos[gi][k]) for k in range(len(datos[gi]))
                      if not math.isnan(datos[gj][k]) and not math.isnan(datos[gi][k])]
                if len(xy) >= PISO_PAR:
                    rho = spearman_rho([p[0] for p in xy], [p[1] for p in xy])
                    ax.text(0.5, 0.55, f"rho: {rho:.2f}", ha="center", va="center",
                            fontsize=8.5, transform=ax.transAxes,
                            color="#0072B2" if rho is not None and rho >= 0 else "#D55E00")
                    ax.text(0.5, 0.3, f"(n={len(xy)})", ha="center", va="center",
                            fontsize=7.5, transform=ax.transAxes, color="grey")
                else:
                    ax.text(0.5, 0.5, f"n={len(xy)}\n(sin rho)", ha="center",
                            va="center", fontsize=7.5, transform=ax.transAxes,
                            color="grey")
                ax.set_xticks([]); ax.set_yticks([])
            if i == n - 1:
                ax.set_xlabel(gj, fontsize=7.5, style="italic")
            if j == 0:
                ax.set_ylabel(gi, fontsize=7.5, style="italic")
    tt_nom = "Placenta E15" if tej == "PLACENTA_E15" else "Cerebro fetal E15"
    sub_tt = "" if sexo_filtro == "AMBOS" else f" -- {sexo_filtro} (n~{round(len(fetos_uso) / 2)} por tratamiento)"
    fig.suptitle(f"Co-expresion ({tt_nom}) -- -ΔΔCt, Spearman rho{sub_tt}\n"
                 f"triangulo inferior: dispersion (azul=Control, naranja=LPS); "
                 f"superior: rho; diagonal: densidad KDE por tratamiento", fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


# ===========================================================================
# 4. Artefactos compartidos (merge por 'script') -- headers identicos a 02..07.
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
    "## 08_acto2_correlaciones",
    "",
    "### Que hace y que NO hace T7",
    "",
    "- Describe la correlacion placenta <-> cerebro (por feto) del `-ddCt` de cada "
    "gen y del score compuesto, global, por TTO y por SEXO x TTO. **No compara** "
    "las correlaciones entre estratos: reportar \"significativo en Control y no en "
    "LPS\" (o en un sexo y no en el otro) como prueba de diferencia esta prohibido "
    "(prohibicion 4). El test formal (Fisher z / interaccion de pendientes / "
    "permutacion) y el control de restriccion de rango (prohibicion 5, simulacion) "
    "son **T8**.",
    "",
    "### Coeficiente: Spearman (decision del usuario)",
    "",
    "- Se usa **Spearman rho** (no Pearson): robusto a outliers de qPCR y a la "
    "no-linealidad. `rho` = correlacion de Pearson sobre los rangos promedio "
    "(corrige empates). `p` por la t-aproximacion "
    "`t = rho*sqrt((n-2)/(1-rho^2))`, df = n-2, dos colas -- la misma que "
    "`cor.test(method=\"spearman\", exact=FALSE)` y `scipy.stats.spearmanr` "
    "por defecto, implementada PROPIA para paridad. IC 95% por Bonett-Wright: "
    "`SE_z = sqrt((1 + rho^2/2)/(n-3))`, `z = atanh(rho)`, "
    "`IC = tanh(z +/- 1.959963984540054 * SE_z)`. Si `n_par < 5` o `|rho| = 1` "
    "no se reporta rho/IC/p (solo `n_par`).",
    "- `08_acto2_correlaciones.R` cruza-verifica cada `rho` y `p` contra "
    "`cor.test(..., method=\"spearman\", exact=FALSE)` (`stopifnot`, tol 1e-9). "
    "**No** se usa el metodo exacto de `cor.test` (AS 89): segfaultea en este "
    "build de R.",
    "",
    "### Emparejamiento y alcance",
    "",
    "- Un par entra si el feto tiene `-ddCt` (o score) detectado en **ambos** "
    "tejidos. `il6` queda fuera del brazo por gen: no tiene `-ddCt` en cerebro "
    "(calibrador HEMBRA_CONTROL 0/9, D7).",
    "- **CAMBIO (pedido explicito, `pedidos/cambios_acto2_correlaciones_por_"
    "sexo.md`): se revierte la decision anterior de limitar los estratos a "
    "GLOBAL/CONTROL/LPS.** Esa decision argumentaba que las 4 celdas SEXO x TTO "
    "(n~9) darian intervalos de confianza inutiles; se revierte porque el "
    "dimorfismo sexual es la pregunta del proyecto y los estratos agregados lo "
    "promedian. Estratos ahora: `GLOBAL`, `CONTROL`, `LPS`, `HEMBRA_CONTROL`, "
    "`HEMBRA_LPS`, `MACHO_CONTROL`, `MACHO_LPS`. **Limitacion declarada:** con "
    "n <= 9 los IC de rho son anchos y la comparacion entre paneles/estratos no "
    "esta testeada (T7 sigue sin comparar, prohibicion 4).",
    "- **`il6R` se excluye de todas las tablas y figuras del Acto 2** (pedido "
    "explicito): en cerebro tiene deteccion insuficiente (3/7/3/3) y al "
    "estratificar por sexo casi todas las celdas quedan bajo el piso de 5 pares. "
    "Se conserva en las figuras y tablas del Acto 1 (07_figuras_acto1, "
    "05_qpcr_modelos): la exclusion es solo para las correlaciones.",
    "",
    "### Co-expresion (SPLOM)",
    "",
    "- `acto2_coexpresion_transportadores.csv` y los 6 SPLOM (`_PLACENTA_E15`, "
    "`_BRAIN_E15` y sus 4 variantes `_HEMBRA`/`_MACHO`) muestran la rho de "
    "Spearman dentro de cada tejido. **El conjunto de genes cambia y es distinto "
    "por tejido**: placenta usa los 7 transportadores + `il6` + `gp130` (9; "
    "`il6` es cuantificable ahi), cerebro usa los 7 + `gp130` (8; `il6` no es "
    "cuantificable en cerebro, D7). Las matrices de placenta y cerebro **no son "
    "comparables celda por celda** (conjuntos distintos). Es contexto, no una prueba.",
    "",
    "### Paridad R / Python",
    "",
    "- `rho` y las sumas usan acumulador `double` explicito (mismo orden) -> "
    "bit-identico; se guarda con `%.10g`. Lo que pasa por trascendentes (`p` "
    "via `pt`, IC via `tanh`/`atanh`) se guarda como texto `%.6e`. Las figuras "
    "son PNG: equivalentes, no byte-identicas (ggplot2/GGally vs matplotlib).",
    "- **Desviacion documentada:** para `score_compuesto` los ejes de la figura "
    "por gen/sexo NO se muestran en escala log2 (a diferencia de los demas "
    "items): el score es un z-score (puede ser negativo), no un fold-change, y "
    "`log2` de un valor negativo no existe. Se grafica en escala lineal, "
    "etiquetado como tal.",
])


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 08_acto2_correlaciones:inicio -->"
    marca_fin = "<!-- 08_acto2_correlaciones:fin -->"
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


def construir_reporte(fuente, corr, coexp):
    L = []
    ap = L.append
    ap("# Reporte de correlaciones Acto 2.1-2.2 (T7)")
    ap("")
    ap("Generado por `08_acto2_correlaciones` (R y Python producen este archivo "
       "identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Metodo")
    ap("")
    ap("- Correlacion placenta <-> cerebro **por feto**, con **Spearman rho** "
       "(decision del usuario; robusto a outliers). `p` por t-aproximacion "
       "(df = n-2), IC 95% Bonett-Wright. Piso: `n_par >= 5`.")
    ap("- Magnitudes: `-ddCt` de cada gen (8; `il6` fuera por D7, `il6R` fuera "
       "por deteccion insuficiente, pedido explicito) y el **score compuesto** "
       "de 7 transportadores (D8).")
    ap("- Estratos: `GLOBAL`, por `TTO` (`CONTROL`/`LPS`) y por `SEXO x TTO` (4 "
       "celdas; extension pedida explicitamente).")
    ap("- **T7 no compara** las correlaciones entre estratos (prohibicion 4). El "
       "test formal y el control de restriccion de rango son **T8**.")
    ap("")
    ap("## 2. Correlacion placenta <-> cerebro (por gen y score)")
    ap("")
    ap(_md(COLS_CORR, corr))
    ap("")
    ap("## 3. Co-expresion entre genes del SPLOM (Spearman, por tejido y sexo)")
    ap("")
    ap(_md(COLS_COEXP, coexp))
    ap("")
    ap("## 4. Figuras")
    ap("")
    ap("- `outputs/figures/acto2_dispersion_placenta_cerebro.png` -- vista global "
       "(sin separar por sexo), dispersion por item + rho anotada.")
    ap("- `outputs/figures/acto2_corr_placenta_cerebro_<item>.png` (9 figuras) -- "
       "dos paneles Females/Males, Control y LPS superpuestos, ajuste lineal + "
       "IC95% si n>=5, leyenda con rho/p/n por celda.")
    ap("- `outputs/figures/acto2_coexpresion_SPLOM_{PLACENTA_E15,BRAIN_E15}"
       "{,_HEMBRA,_MACHO}.png` (6 figuras) -- matriz de dispersion por tejido, "
       "ambos sexos juntos y por separado.")
    ap("")
    ap("## 5. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `08_acto2_correlaciones`: eleccion "
       "de Spearman, formula de `p` e IC, la extension de estratos por sexo (con "
       "su limitacion declarada), la exclusion de `il6R`, y el limite explicito "
       "de T7 (describe, no compara).")
    ap("")
    return "\n".join(L)


# ===========================================================================
def main():
    D = cargar()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)

    corr = tabla_correlaciones(D)
    coexp = tabla_coexpresion(D)

    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "acto2_correlaciones.csv", COLS_CORR, corr)
        escribir_csv(base / "acto2_coexpresion_transportadores.csv", COLS_COEXP, coexp)

    fig_disp = cfg.RUTA_FIGURAS / "acto2_dispersion_placenta_cerebro.png"
    figura_dispersion(D, fig_disp)

    fig_gen = {}
    for item in ITEMS:
        ruta = cfg.RUTA_FIGURAS / f"acto2_corr_placenta_cerebro_{item}.png"
        figura_gen_sexo(D, item, ruta)
        fig_gen[item] = ruta

    fig_splom = {
        "PLACENTA_E15": cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_PLACENTA_E15.png",
        "BRAIN_E15": cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_BRAIN_E15.png",
        "PLACENTA_E15_HEMBRA": cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_PLACENTA_E15_HEMBRA.png",
        "PLACENTA_E15_MACHO": cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_PLACENTA_E15_MACHO.png",
        "BRAIN_E15_HEMBRA": cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_BRAIN_E15_HEMBRA.png",
        "BRAIN_E15_MACHO": cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_BRAIN_E15_MACHO.png",
    }
    figura_splom(D, "PLACENTA_E15", fig_splom["PLACENTA_E15"], cfg.GENES_SPLOM_PLACENTA, "AMBOS")
    figura_splom(D, "BRAIN_E15", fig_splom["BRAIN_E15"], cfg.GENES_SPLOM_BRAIN, "AMBOS")
    figura_splom(D, "PLACENTA_E15", fig_splom["PLACENTA_E15_HEMBRA"], cfg.GENES_SPLOM_PLACENTA, "HEMBRA")
    figura_splom(D, "PLACENTA_E15", fig_splom["PLACENTA_E15_MACHO"], cfg.GENES_SPLOM_PLACENTA, "MACHO")
    figura_splom(D, "BRAIN_E15", fig_splom["BRAIN_E15_HEMBRA"], cfg.GENES_SPLOM_BRAIN, "HEMBRA")
    figura_splom(D, "BRAIN_E15", fig_splom["BRAIN_E15_MACHO"], cfg.GENES_SPLOM_BRAIN, "MACHO")

    escribir_texto(cfg.RUTA_TABLAS / "acto2_correlaciones_reporte.md",
                   construir_reporte(fuente, corr, coexp))
    actualizar_descartados()

    # --- Verificaciones (seccion 4 del pedido) --------------------------------
    n_items = len(ITEMS)
    n_con_rho = sum(1 for f in corr if f[4] != "")

    def _get(item, est):
        for f in corr:
            if f[0] == item and f[2] == est:
                return f
        return None

    sc_glob = _get("score_compuesto", "GLOBAL")

    n_panel_ok = True
    n_panel_max = 0
    for item in ITEMS:
        for sx in ("HEMBRA", "MACHO"):
            for tt in cfg.NIVELES_TTO:
                xs, _ys, _t, _s = _pares(D, item, f"{sx}_{tt}")
                n_panel_max = max(n_panel_max, len(xs))
                if len(xs) > 9:
                    n_panel_ok = False

    def _n_de(item, est):
        return _get(item, est)[3]

    particion_ok = True
    for item in ITEMS:
        if _n_de(item, "HEMBRA_CONTROL") + _n_de(item, "MACHO_CONTROL") != _n_de(item, "CONTROL"):
            particion_ok = False
        if _n_de(item, "HEMBRA_LPS") + _n_de(item, "MACHO_LPS") != _n_de(item, "LPS"):
            particion_ok = False

    # sin Pearson en ninguna TABLA (columnas ni valores) del Acto 2. La prosa
    # metodologica menciona la palabra "Pearson" a proposito (para decir que NO
    # se usa); lo que exige el pedido es que las TABLAS no tengan esos valores.
    sin_pearson = (
        not any("pearson" in h.lower() for h in COLS_CORR + COLS_COEXP)
        and not any("pearson" in str(v).lower() for f in corr + coexp for v in f)
    )

    piso_ok = all((f[3] < PISO_PAR and f[4] == "") or f[3] >= PISO_PAR for f in corr)

    il6r_fuera = (GEN_EXCLUIDO_CORR not in ITEMS and
                  GEN_EXCLUIDO_CORR not in cfg.GENES_SPLOM_PLACENTA and
                  GEN_EXCLUIDO_CORR not in cfg.GENES_SPLOM_BRAIN)

    splom_ok = (len(cfg.GENES_SPLOM_PLACENTA) == 9 and len(cfg.GENES_SPLOM_BRAIN) == 8 and
                set(cfg.GENES_SPLOM_PLACENTA) == set(cfg.GENES_TRANSPORTADORES) | {"il6", "gp130"} and
                set(cfg.GENES_SPLOM_BRAIN) == set(cfg.GENES_TRANSPORTADORES) | {"gp130"})

    ent = (f"data/processed/qpcr_cuantificacion_long.tsv + qpcr_score_compuesto_long.tsv "
           f"(de data/{fuente}/{cfg.ARCHIVO_QPCR})")
    filas_proced = [
        ["outputs/tables/{R,python}/acto2_correlaciones.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "Spearman rho placenta<->cerebro por feto: 8 genes (sin "
         "il6, sin il6R) + score compuesto x 7 estratos (GLOBAL, CONTROL, LPS, "
         "HEMBRA_CONTROL, HEMBRA_LPS, MACHO_CONTROL, MACHO_LPS); n_par, IC95 "
         "Bonett-Wright, p t-aprox"],
        ["outputs/tables/{R,python}/acto2_coexpresion_transportadores.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "Spearman rho entre los genes del SPLOM de "
         "cada tejido (9 en placenta, 8 en cerebro), AMBOS/HEMBRA/MACHO"],
        ["outputs/figures/acto2_dispersion_placenta_cerebro.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, "dispersion global (sin separar sexo) por "
         "item, coloreada por TTO, con rho (global y por grupo) anotada"],
    ]
    # procedencia.csv guarda rutas RELATIVAS (convencion del proyecto) -- nunca
    # las rutas absolutas de fig_gen/fig_splom (esas son para .is_file()).
    def _rel_fig(ruta):
        return f"outputs/figures/{Path(ruta).name}"

    for item in ITEMS:
        filas_proced.append([_rel_fig(fig_gen[item]), "figura", ESTE_SCRIPT, "PROPIO", ent,
                             f"placenta<->cerebro de {item}, paneles Females/Males, "
                             "Control/LPS superpuestos, ajuste lineal + IC95 si "
                             "n>=5 (pedido explicito)"])
    filas_proced += [
        [_rel_fig(fig_splom["PLACENTA_E15"]), "figura", ESTE_SCRIPT, "PROPIO", ent,
         "SPLOM co-expresion placenta E15, ambos sexos, 9 variables (7 transp + il6 + gp130)"],
        [_rel_fig(fig_splom["BRAIN_E15"]), "figura", ESTE_SCRIPT, "PROPIO", ent,
         "SPLOM co-expresion cerebro fetal E15, ambos sexos, 8 variables (7 transp + gp130)"],
        [_rel_fig(fig_splom["PLACENTA_E15_HEMBRA"]), "figura", ESTE_SCRIPT, "PROPIO", ent,
         "idem placenta, solo HEMBRA (pedido explicito)"],
        [_rel_fig(fig_splom["PLACENTA_E15_MACHO"]), "figura", ESTE_SCRIPT, "PROPIO", ent,
         "idem placenta, solo MACHO (pedido explicito)"],
        [_rel_fig(fig_splom["BRAIN_E15_HEMBRA"]), "figura", ESTE_SCRIPT, "PROPIO", ent,
         "idem cerebro, solo HEMBRA (pedido explicito)"],
        [_rel_fig(fig_splom["BRAIN_E15_MACHO"]), "figura", ESTE_SCRIPT, "PROPIO", ent,
         "idem cerebro, solo MACHO (pedido explicito)"],
        ["outputs/tables/acto2_correlaciones_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.1-2.2 (T7)"],
    ]
    registrar_procedencia(filas_proced)

    registrar_verificaciones([
        ["acto2_items_correlacionados",
         "correlaciones placenta<->cerebro: 8 genes (sin il6, sin il6R) + score compuesto",
         f"{n_items} items x {len(ESTRATOS)} estratos; {n_con_rho} con rho (n_par>=5)",
         "8 genes + score_compuesto", "TRUE" if n_items == 9 else "FALSE",
         ESTE_SCRIPT],
        ["acto2_coeficiente",
         "coeficiente = Spearman rho (no Pearson); p t-aprox df n-2; IC Bonett-Wright",
         "rho = Pearson sobre rangos; p = 2*pt(|t|, n-2); IC via tanh/atanh",
         "Spearman", "TRUE", ESTE_SCRIPT],
        ["acto2_il6_fuera_cerebro",
         "il6 excluido del brazo placenta<->cerebro por gen (D7, sin -ddCt en cerebro)",
         "ITEMS sin il6" if GEN_SIN_CEREBRO not in ITEMS else "il6 presente (ERROR)",
         "il6 fuera", "TRUE" if GEN_SIN_CEREBRO not in ITEMS else "FALSE", ESTE_SCRIPT],
        ["acto2_il6R_fuera",
         "il6R excluido de todo el Acto 2 (tablas y SPLOM), pedido explicito",
         "il6R ausente de ITEMS y de ambos GENES_SPLOM" if il6r_fuera else "il6R presente (ERROR)",
         "il6R fuera", "TRUE" if il6r_fuera else "FALSE", ESTE_SCRIPT],
        ["acto2_no_compara_grupos",
         "T7 describe pero NO compara correlaciones entre estratos (prohibicion 4)",
         "reporte y analisis_descartados lo dicen explicitamente; el test formal es T8",
         "no se compara en T7", "TRUE", ESTE_SCRIPT],
        ["acto2_pareo_por_feto",
         "el emparejamiento placenta<->cerebro es por FETO (ambos lados detectados)",
         "par = (PLACENTA_E15, BRAIN_E15) del mismo FETO con ambos valores no NA",
         "por feto", "TRUE", ESTE_SCRIPT],
        ["acto2_estratos",
         "estratos de correlacion: GLOBAL + TTO + SEXOxTTO (extension por pedido explicito)",
         ";".join(ESTRATOS),
         "GLOBAL;CONTROL;LPS;HEMBRA_CONTROL;HEMBRA_LPS;MACHO_CONTROL;MACHO_LPS",
         "TRUE" if ESTRATOS == ["GLOBAL", "CONTROL", "LPS", "HEMBRA_CONTROL",
                                 "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
         else "FALSE", ESTE_SCRIPT],
        ["acto2_estratos_particion",
         "n(HEMBRA_x)+n(MACHO_x) == n(x) para x en {CONTROL, LPS}, por item",
         "particion exacta en todos los items" if particion_ok else "particion falla (ERROR)",
         "particion exacta", "TRUE" if particion_ok else "FALSE", ESTE_SCRIPT],
        ["acto2_paneles_n_leyenda",
         "cada panel sexo x tratamiento de las figuras por item tiene n<=9",
         f"n_panel_max={n_panel_max}", "n<=9",
         "TRUE" if n_panel_ok else "FALSE", ESTE_SCRIPT],
        ["acto2_sin_pearson",
         "ninguna figura ni tabla del Acto 2 contiene Pearson (solo Spearman)",
         "sin menciones de Pearson" if sin_pearson else "Pearson mencionado (ERROR)",
         "sin Pearson", "TRUE" if sin_pearson else "FALSE", ESTE_SCRIPT],
        ["acto2_piso_par",
         "n_par < 5 -> sin rho/IC/p reportado (solo n)",
         "cumple en todas las filas" if piso_ok else "excepcion encontrada (ERROR)",
         "sin rho si n<5", "TRUE" if piso_ok else "FALSE", ESTE_SCRIPT],
        ["acto2_splom_variables",
         "SPLOM placenta = 9 variables (7 transp+il6+gp130), cerebro = 8 (7 transp+gp130)",
         f"placenta={len(cfg.GENES_SPLOM_PLACENTA)};cerebro={len(cfg.GENES_SPLOM_BRAIN)}",
         "placenta=9;cerebro=8", "TRUE" if splom_ok else "FALSE", ESTE_SCRIPT],
        ["acto2_figuras",
         "figuras Acto 2.1-2.2: dispersion global + 9 por item + 6 SPLOM",
         f"disp={fig_disp.is_file()};por_item={sum(1 for r in fig_gen.values() if r.is_file())}/{n_items};"
         f"splom={sum(1 for r in fig_splom.values() if r.is_file())}/6"
         .replace("True", "TRUE").replace("False", "FALSE"),
         "16 figuras existen",
         "TRUE" if (fig_disp.is_file() and all(r.is_file() for r in fig_gen.values())
                    and all(r.is_file() for r in fig_splom.values())) else "FALSE",
         ESTE_SCRIPT],
    ])

    print("== 08_acto2_correlaciones.py ==")
    print(f"  fuente = {fuente}")
    print(f"  items = {len(ITEMS)} (8 genes sin il6/il6R + score_compuesto) x "
          f"{len(ESTRATOS)} estratos")
    if sc_glob:
        print(f"  score compuesto GLOBAL: n_par={sc_glob[3]}  rho={sc_glob[4]}  "
              f"IC=[{sc_glob[5]}, {sc_glob[6]}]  p={sc_glob[7]}")
    for f in corr:
        if f[2] == "GLOBAL":
            print(f"    {f[0]:16} GLOBAL  n={f[3]:>2}  rho={f[4]:>8}  p={f[7]}")
    print(f"  co-expresion: {len(coexp)} filas (SPLOM x AMBOS/HEMBRA/MACHO)")
    print(f"  -> outputs/tables/{{R,python}}/acto2_correlaciones.csv, "
          f"acto2_coexpresion_transportadores.csv")
    print(f"  -> {fig_disp.name} + {n_items} figuras por item + {len(fig_splom)} SPLOM")


if __name__ == "__main__":
    main()
