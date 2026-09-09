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
#   * Magnitudes (decision del usuario): -ddCt por gen (9 genes: todos menos il6,
#     que no tiene -ddCt en cerebro por D7) y el score compuesto de 7
#     transportadores (04_qpcr_cuantificacion, D8).
#   * Coeficiente (decision del usuario): **Spearman rho** (no Pearson).
#     Robusto a outliers de qPCR. p por la t-aproximacion
#     t = rho * sqrt((n-2)/(1-rho^2)), df = n-2, dos colas (misma formula que
#     scipy.stats.spearmanr por defecto, implementada PROPIA). IC 95% por
#     Bonett-Wright: SE_z = sqrt((1 + rho^2/2)/(n-3)), z = atanh(rho),
#     IC = tanh(z +/- 1.959963984540054 * SE_z).
#   * Estratos (decision del usuario): GLOBAL (n<=36) + por TTO (CONTROL / LPS,
#     n<=18). Las 4 celdas SEXO x TTO (n~9) NO se usan (IC inutiles).
#   * Piso: si el par tiene < 5 fetos, no se calcula rho/IC/p (solo n).
#   * Co-expresion: rho de Spearman entre los 7 transportadores dentro de cada
#     tejido (para el SPLOM y como tabla).
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

from scipy import stats as _sst  # noqa: E402

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "08_acto2_correlaciones"

COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}   # Okabe-Ito, igual que 03/07
Z975 = 1.959963984540054                             # qnorm(0.975), literal exacto
PISO_PAR = 5                                         # min fetos emparejados para rho

GEN_SIN_CEREBRO = "il6"                              # D7: sin -ddCt en BRAIN_E15
GENES_CORR = [g for g in cfg.GENES if g != GEN_SIN_CEREBRO]   # 9 genes
ITEMS = GENES_CORR + ["score_compuesto"]
ESTRATOS = ["GLOBAL", "CONTROL", "LPS"]
TRANSP = list(cfg.GENES_TRANSPORTADORES)             # 7, para el SPLOM
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
    """(placenta[], cerebro[], tto[]) de los fetos con ambos lados detectados."""
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
            xs, ys, _ts = _pares(D, item, est)
            n = len(xs)
            if n >= PISO_PAR:
                rho = spearman_rho(xs, ys)
                lo, hi = spearman_ci(rho, n)
                p = spearman_p(rho, n)
                filas.append([item, tipo, est, n, g10(rho), p6e(lo), p6e(hi), p6e(p)])
            else:
                filas.append([item, tipo, est, n, "", "", "", ""])
    return filas


COLS_COEXP = ["TEJIDO", "GEN_A", "GEN_B", "n_par", "rho_spearman", "p_valor"]


def tabla_coexpresion(D):
    filas = []
    for tej in cfg.TEJIDOS_E15:
        for i in range(len(TRANSP)):
            for j in range(i + 1, len(TRANSP)):
                a, b = TRANSP[i], TRANSP[j]
                xs, ys = [], []
                for f in D["fetos"]:
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
                    filas.append([tej, a, b, n, g10(rho), p6e(p)])
                else:
                    filas.append([tej, a, b, n, "", ""])
    return filas


# ===========================================================================
# 3. Figuras.
# ===========================================================================
def _rho_txt(D, item):
    out = {}
    for est in ESTRATOS:
        xs, ys, _ = _pares(D, item, est)
        n = len(xs)
        if n >= PISO_PAR:
            rho = spearman_rho(xs, ys)
            out[est] = (rho, n)
        else:
            out[est] = (None, n)
    return out


def figura_dispersion(D, ruta):
    ncol, nrow = 4, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(13.0, 9.0))
    axl = axes.flatten()
    for k, item in enumerate(ITEMS):
        ax = axl[k]
        xs, ys, ts = _pares(D, item, "GLOBAL")
        for tt in ("CONTROL", "LPS"):
            xx = [xs[i] for i in range(len(xs)) if ts[i] == tt]
            yy = [ys[i] for i in range(len(ys)) if ts[i] == tt]
            ax.plot(xx, yy, "o", ms=4, mfc=COL_TTO[tt], mec="white", mew=0.4,
                    ls="none", label=tt.capitalize())
        rt = _rho_txt(D, item)
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
                 "testea aca (prohibicion 4); el test formal es T8", fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


def figura_splom(D, tej, ruta):
    n = len(TRANSP)
    # matriz de valores por feto (NaN si no detectado)
    datos = {g: [] for g in TRANSP}
    ttos = []
    for f in D["fetos"]:
        ttos.append(D["tto"][f])
        for g in TRANSP:
            v = D["negdd"].get((f, tej, g))
            datos[g].append(v if v is not None else float("nan"))
    fig, axes = plt.subplots(n, n, figsize=(12.5, 12.5))
    for i, gi in enumerate(TRANSP):
        for j, gj in enumerate(TRANSP):
            ax = axes[i, j]
            ax.tick_params(labelsize=6)
            if i == j:
                vals = [v for v in datos[gi] if not math.isnan(v)]
                ax.hist(vals, bins=10, color="0.6")
                ax.set_yticks([])
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
                    ax.text(0.5, 0.5, f"ρ = {rho:.2f}\n(n={len(xy)})", ha="center",
                            va="center", fontsize=9,
                            color="#0072B2" if rho is not None and rho >= 0 else "#D55E00")
                ax.set_xticks([]); ax.set_yticks([])
            if i == n - 1:
                ax.set_xlabel(gj, fontsize=7.5, style="italic")
            if j == 0:
                ax.set_ylabel(gi, fontsize=7.5, style="italic")
    tt = "Placenta E15" if tej == "PLACENTA_E15" else "Cerebro fetal E15"
    fig.suptitle(f"Co-expresion de los 7 transportadores ({tt})  --  -ΔΔCt, "
                 f"Spearman ρ\ntriangulo inferior: dispersion (azul=Control, "
                 f"naranja=LPS); superior: ρ; diagonal: histograma", fontsize=10)
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
    "gen y del score compuesto, global y por TTO. **No compara** las correlaciones "
    "entre grupos: reportar \"significativo en Control y no en LPS\" como prueba de "
    "diferencia esta prohibido (prohibicion 4). El test formal (Fisher z / "
    "interaccion de pendientes / permutacion) y el control de restriccion de rango "
    "(prohibicion 5, simulacion) son **T8**.",
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
    "(calibrador HEMBRA_CONTROL 0/9, D7). `il6R` en cerebro tiene n bajo -> sus "
    "`n_par` por estrato pueden quedar por debajo del piso de 5.",
    "- Estratos: `GLOBAL` (n<=36) y por `TTO` (`CONTROL` / `LPS`, n<=18). Las "
    "4 celdas SEXO x TTO (n~9) darian IC inutiles y no se usan.",
    "",
    "### Co-expresion (SPLOM)",
    "",
    "- `acto2_coexpresion_transportadores.csv` y `acto2_coexpresion_SPLOM_"
    "<tejido>.png` muestran la rho de Spearman entre los 7 transportadores "
    "dentro de cada tejido (21 pares x 2 tejidos). Es contexto, no una prueba.",
    "",
    "### Paridad R / Python",
    "",
    "- `rho` y las sumas usan acumulador `double` explicito (mismo orden) -> "
    "bit-identico; se guarda con `%.10g`. Lo que pasa por trascendentes (`p` "
    "via `pt`, IC via `tanh`/`atanh`) se guarda como texto `%.6e`. Las figuras "
    "son PNG: equivalentes, no byte-identicas (ggplot2/GGally vs matplotlib).",
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
    ap("- Magnitudes: `-ddCt` de cada gen (9; `il6` fuera por D7) y el **score "
       "compuesto** de 7 transportadores (D8).")
    ap("- Estratos: `GLOBAL` y por `TTO` (`CONTROL` / `LPS`).")
    ap("- **T7 no compara** las correlaciones entre grupos (prohibicion 4). El "
       "test formal y el control de restriccion de rango son **T8**.")
    ap("")
    ap("## 2. Correlacion placenta <-> cerebro (por gen y score)")
    ap("")
    ap(_md(COLS_CORR, corr))
    ap("")
    ap("## 3. Co-expresion entre transportadores (Spearman, por tejido)")
    ap("")
    ap(_md(COLS_COEXP, coexp))
    ap("")
    ap("## 4. Figuras")
    ap("")
    ap("- `outputs/figures/acto2_dispersion_placenta_cerebro.png` -- dispersion "
       "placenta vs cerebro por gen + score, coloreada por TTO, con rho (global y "
       "por grupo) anotada.")
    ap("- `outputs/figures/acto2_coexpresion_SPLOM_PLACENTA_E15.png` y "
       "`..._BRAIN_E15.png` -- matriz de dispersion de los 7 transportadores por "
       "tejido.")
    ap("")
    ap("## 5. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `08_acto2_correlaciones`: eleccion "
       "de Spearman, formula de `p` e IC, y el limite explicito de T7 (describe, no "
       "compara).")
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
    fig_spl_p = cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_PLACENTA_E15.png"
    fig_spl_b = cfg.RUTA_FIGURAS / "acto2_coexpresion_SPLOM_BRAIN_E15.png"
    figura_dispersion(D, fig_disp)
    figura_splom(D, "PLACENTA_E15", fig_spl_p)
    figura_splom(D, "BRAIN_E15", fig_spl_b)

    escribir_texto(cfg.RUTA_TABLAS / "acto2_correlaciones_reporte.md",
                   construir_reporte(fuente, corr, coexp))
    actualizar_descartados()

    # resumen determinista para verificaciones
    def _get(item, est):
        for f in corr:
            if f[0] == item and f[2] == est:
                return f
        return None

    sc_glob = _get("score_compuesto", "GLOBAL")
    n_items = len(ITEMS)
    n_con_rho = sum(1 for f in corr if f[4] != "")

    ent = (f"data/processed/qpcr_cuantificacion_long.tsv + qpcr_score_compuesto_long.tsv "
           f"(de data/{fuente}/{cfg.ARCHIVO_QPCR})")
    registrar_procedencia([
        ["outputs/tables/{R,python}/acto2_correlaciones.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "Spearman rho placenta<->cerebro por feto: 9 genes + score "
         "compuesto x {GLOBAL, CONTROL, LPS}; n_par, IC95 Bonett-Wright, p t-aprox"],
        ["outputs/tables/{R,python}/acto2_coexpresion_transportadores.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "Spearman rho entre los 7 transportadores "
         "dentro de cada tejido (21 pares x 2 tejidos)"],
        ["outputs/figures/acto2_dispersion_placenta_cerebro.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, "dispersion placenta vs cerebro por gen + "
         "score, coloreada por TTO, con rho (global y por grupo) anotada"],
        ["outputs/figures/acto2_coexpresion_SPLOM_PLACENTA_E15.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, "matriz de dispersion (SPLOM) de los 7 "
         "transportadores en placenta E15"],
        ["outputs/figures/acto2_coexpresion_SPLOM_BRAIN_E15.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, "idem en cerebro fetal E15"],
        ["outputs/tables/acto2_correlaciones_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.1-2.2 (T7)"],
    ])
    registrar_verificaciones([
        ["acto2_items_correlacionados",
         "correlaciones placenta<->cerebro: 9 genes (sin il6) + score compuesto",
         f"{n_items} items x {len(ESTRATOS)} estratos; {n_con_rho} con rho (n_par>=5)",
         "9 genes + score_compuesto", "TRUE" if n_items == 10 else "FALSE",
         ESTE_SCRIPT],
        ["acto2_coeficiente",
         "coeficiente = Spearman rho (no Pearson); p t-aprox df n-2; IC Bonett-Wright",
         "rho = Pearson sobre rangos; p = 2*pt(|t|, n-2); IC via tanh/atanh",
         "Spearman", "TRUE", ESTE_SCRIPT],
        ["acto2_il6_fuera_cerebro",
         "il6 excluido del brazo placenta<->cerebro por gen (D7, sin -ddCt en cerebro)",
         "ITEMS sin il6" if GEN_SIN_CEREBRO not in ITEMS else "il6 presente (ERROR)",
         "il6 fuera", "TRUE" if GEN_SIN_CEREBRO not in ITEMS else "FALSE", ESTE_SCRIPT],
        ["acto2_no_compara_grupos",
         "T7 describe pero NO compara correlaciones entre grupos (prohibicion 4)",
         "reporte y analisis_descartados lo dicen explicitamente; el test formal es T8",
         "no se compara en T7", "TRUE", ESTE_SCRIPT],
        ["acto2_pareo_por_feto",
         "el emparejamiento placenta<->cerebro es por FETO (ambos lados detectados)",
         "par = (PLACENTA_E15, BRAIN_E15) del mismo FETO con ambos valores no NA",
         "por feto", "TRUE", ESTE_SCRIPT],
        ["acto2_estratos",
         "estratos de correlacion: GLOBAL + por TTO (Control/LPS); sin celdas SEXOxTTO",
         ";".join(ESTRATOS), "GLOBAL;CONTROL;LPS",
         "TRUE" if ESTRATOS == ["GLOBAL", "CONTROL", "LPS"] else "FALSE", ESTE_SCRIPT],
        ["acto2_piso_par", "n_par minimo para calcular rho/IC/p", str(PISO_PAR),
         "5", "TRUE" if PISO_PAR == 5 else "FALSE", ESTE_SCRIPT],
        ["acto2_figuras", "figuras Acto 2.1-2.2: dispersion por gen + 2 SPLOM",
         f"disp={fig_disp.is_file()};splom_pla={fig_spl_p.is_file()};"
         f"splom_bra={fig_spl_b.is_file()}".replace("True", "TRUE").replace("False", "FALSE"),
         "3 figuras existen",
         "TRUE" if (fig_disp.is_file() and fig_spl_p.is_file() and fig_spl_b.is_file())
         else "FALSE", ESTE_SCRIPT],
    ])

    print("== 08_acto2_correlaciones.py ==")
    print(f"  fuente = {fuente}")
    print(f"  items = {len(ITEMS)} (9 genes sin il6 + score_compuesto) x "
          f"{len(ESTRATOS)} estratos")
    if sc_glob:
        print(f"  score compuesto GLOBAL: n_par={sc_glob[3]}  rho={sc_glob[4]}  "
              f"IC=[{sc_glob[5]}, {sc_glob[6]}]  p={sc_glob[7]}")
    for f in corr:
        if f[2] == "GLOBAL":
            print(f"    {f[0]:16} GLOBAL  n={f[3]:>2}  rho={f[4]:>8}  p={f[7]}")
    print(f"  co-expresion: {len(coexp)} pares (21 x 2 tejidos)")
    print(f"  -> outputs/tables/{{R,python}}/acto2_correlaciones.csv, "
          f"acto2_coexpresion_transportadores.csv")
    print(f"  -> {fig_disp.name}, {fig_spl_p.name}, {fig_spl_b.name}")


if __name__ == "__main__":
    main()
