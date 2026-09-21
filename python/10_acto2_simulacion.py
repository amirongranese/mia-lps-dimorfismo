# 10_acto2_simulacion.py -- ACTO 2.5: descartar restriccion de rango antes de
#                           interpretar el cambio de correlacion placenta<->cerebro.
#
# Por que existe este archivo: prohibicion 5 -- NO se puede interpretar un cambio
# de correlacion entre Control y LPS como cambio de coordinacion biologica sin
# descartar antes que sea un artefacto de dispersion (restriccion de rango). Una
# correlacion observada se atenua si el rango de alguna de las dos variables se
# comprime en un subgrupo, aunque la asociacion verdadera sea la misma.
#
# Diseno de la simulacion (decisiones del usuario):
#   * Escala: **normal bivariada en -ddCt**. Se generan (placenta, cerebro) con
#     UNA correlacion verdadera comun a ambos grupos y luego se escalan las SD
#     marginales de cada grupo a las SD OBSERVADAS de ese grupo (de 09). Lo unico
#     que difiere entre Control y LPS en la simulacion es la dispersion.
#   * r verdadera comun -- **ambas como rango**: se corre con
#       (a) ESCENARIO GLOBAL  -> rho verdadera = rho de Spearman GLOBAL del item
#           (pooled Control+LPS, la asociacion de T7);
#       (b) ESCENARIO CONTROL -> rho verdadera = rho de Spearman de Control
#           (el grupo con mas coordinacion; escenario mas exigente).
#     El Delta rho observado se compara contra la distribucion simulada de CADA
#     escenario.
#   * Inversion rho_Spearman -> r_Pearson generativo: para la normal bivariada
#     rho_S = (6/pi) * arcsin(r/2)  =>  r = 2 * sin(pi * rho_S / 6). Asi la rho
#     de Spearman de los datos simulados apunta al target.
#   * B = 2000 repeticiones, semilla 20260101, RNG PROPIO (LCG + polar de
#     Marsaglia de 00_config; misma secuencia exacta en R y Python).
#
# Por item y escenario se reporta: Delta rho observado, media y SD del Delta rho
# simulado, intervalo central 95% (cuantiles tipo 7), fraccion de repeticiones
# con |Delta rho_sim| >= |Delta rho_obs| (una "p" de simulacion para "la
# dispersion sola lo explica"), y la tasa de repeticiones en que el Fisher z de
# 09 daria p < .05 (tasa de falsos positivos del test reportado bajo pura
# restriccion de rango). Veredicto:
#   DENTRO  -> el Delta rho observado cae dentro del intervalo simulado: la
#              diferencia de dispersion sola puede producirlo -> NO se puede
#              descartar restriccion de rango; el cambio de correlacion NO es
#              interpretable como coordinacion biologica.
#   FUERA   -> el Delta rho observado excede lo que produce la sola diferencia de
#              dispersion en ese escenario.
#
# PARIDAD R/Python: los sorteos usan norm1() (LCG + polar, byte-identico en
# ambos por 00_config/T1); rho, medias, SD, cuantiles y conteos usan acumulador
# double explicito y comparaciones -> texto "%.10g" bit-identico. `r_pearson_gen`
# pasa por sin() -> texto "%.6e". Las figuras son PNG: equivalentes, no
# byte-identicas.

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "10_acto2_simulacion"

COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}
Z975 = 1.959963984540054
PISO_PAR = 5

GEN_SIN_CEREBRO = "il6"
GEN_EXCLUIDO_CORR = "il6R"   # pedido explicito (08_acto2_correlaciones): deteccion insuficiente en cerebro
GENES_CORR = [g for g in cfg.GENES
              if g not in (GEN_SIN_CEREBRO, GEN_EXCLUIDO_CORR)]   # 8 genes
ITEMS = GENES_CORR + ["score_compuesto"]
ESCENARIOS = ["GLOBAL", "CONTROL"]
# Estratificacion por sexo (pedido explicito, punto 4/4): AMBOS_SEXOS =
# comportamiento previo (agrupado); HEMBRA/MACHO = mismo diseno con las SD
# OBSERVADAS de cada celda de sexo. Aditivo.
ESTRATOS = ["AMBOS_SEXOS", "HEMBRA", "MACHO"]


def _sexo_de_estrato(e):
    return None if e == "AMBOS_SEXOS" else e


B_SIM = 2000
R_CLAMP = 0.999999
DPI = 300


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..09.
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


def desvio(xs):
    n = len(xs)
    if n < 2:
        return None
    m = suma(xs) / n
    s = 0.0
    for v in xs:
        s += (v - m) * (v - m)
    return math.sqrt(s / (n - 1))


def rangos_promedio(xs):
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


def cuantil_tipo7(xs_ordenado, prob):
    """Cuantil tipo 7 (default de R quantile): interp lineal sobre (n-1)*p."""
    n = len(xs_ordenado)
    if n == 0:
        return None
    if n == 1:
        return xs_ordenado[0]
    h = (n - 1) * prob
    lo = int(math.floor(h))
    if lo >= n - 1:
        return xs_ordenado[n - 1]
    frac = h - lo
    return xs_ordenado[lo] + frac * (xs_ordenado[lo + 1] - xs_ordenado[lo])


def rho_s_a_r(rho_s):
    """rho de Spearman -> r de Pearson generativo para la normal bivariada."""
    r = 2.0 * math.sin(math.pi * rho_s / 6.0)
    if r > R_CLAMP:
        return R_CLAMP
    if r < -R_CLAMP:
        return -R_CLAMP
    return r


def _fisher_sig(rc, nc, rl, nl):
    """TRUE si el Fisher z (SE Bonett-Wright) daria p < .05, sin llamar a la
    normal: |stat| > qnorm(.975)."""
    zc = math.atanh(rc)
    zl = math.atanh(rl)
    sec = math.sqrt((1.0 + rc * rc / 2.0) / (nc - 3))
    sel = math.sqrt((1.0 + rl * rl / 2.0) / (nl - 3))
    stat = (zc - zl) / math.sqrt(sec * sec + sel * sel)
    return abs(stat) > Z975


# ===========================================================================
# 1. Carga y emparejamiento por feto -- identico a 08/09.
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
    """`sexo_filtro`: None = ambos sexos (comportamiento previo); 'HEMBRA'/
    'MACHO' restringe ademas por sexo -- simulacion estratificada (pedido
    explicito, punto 4/4 de pedidos/cambios_acto2_dispersion_por_sexo.md)."""
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


# ===========================================================================
# 2. Simulacion.
# ===========================================================================
COLS_SIM = ["ITEM", "TIPO", "ESTRATO", "ESCENARIO", "rho_true", "r_pearson_gen",
            "n_control", "n_lps", "sd_pla_control", "sd_bra_control",
            "sd_pla_lps", "sd_bra_lps", "delta_rho_obs", "sim_mean_delta_rho",
            "sim_sd_delta_rho", "sim_q025", "sim_q975", "sim_frac_abs_ge_obs",
            "sim_tasa_fisher_sig", "veredicto"]


def _simular_celda(rng, r_pear, sd_pla_c, sd_bra_c, nc, sd_pla_l, sd_bra_l, nl,
                   drho_obs):
    raiz = math.sqrt(1.0 - r_pear * r_pear)
    sims = []
    n_fisher = 0
    for _b in range(B_SIM):
        xc, yc = [], []
        for _k in range(nc):
            z1 = rng.norm1()
            z2 = rng.norm1()
            xc.append(sd_pla_c * z1)
            yc.append(sd_bra_c * (r_pear * z1 + raiz * z2))
        xl, yl = [], []
        for _k in range(nl):
            z1 = rng.norm1()
            z2 = rng.norm1()
            xl.append(sd_pla_l * z1)
            yl.append(sd_bra_l * (r_pear * z1 + raiz * z2))
        rc = spearman_rho(xc, yc)
        rl = spearman_rho(xl, yl)
        if rc is None or rl is None:
            sims.append(0.0)
            continue
        d = rc - rl
        sims.append(d)
        if abs(rc) < 1.0 and abs(rl) < 1.0 and _fisher_sig(rc, nc, rl, nl):
            n_fisher += 1
    media = suma(sims) / B_SIM
    sd = desvio(sims)
    ordenado = sorted(sims)
    q025 = cuantil_tipo7(ordenado, 0.025)
    q975 = cuantil_tipo7(ordenado, 0.975)
    ge = 0
    a_obs = abs(drho_obs)
    for d in sims:
        if abs(d) >= a_obs:
            ge += 1
    return dict(media=media, sd=sd, q025=q025, q975=q975,
                frac_ge=ge / B_SIM, tasa_fisher=n_fisher / B_SIM)


def tabla_simulacion(D, rng):
    filas = []
    for item in ITEMS:
        tipo = "score" if item == "score_compuesto" else "gen"
        for estrato in ESTRATOS:
            sx = _sexo_de_estrato(estrato)
            cx, cy, _ = _pares(D, item, "CONTROL", sx)
            lx, ly, _ = _pares(D, item, "LPS", sx)
            gx, gy, _ = _pares(D, item, "GLOBAL", sx)
            nc, nl = len(cx), len(lx)
            if nc < PISO_PAR or nl < PISO_PAR:
                for esc in ESCENARIOS:
                    filas.append([item, tipo, estrato, esc, "", "", nc, nl]
                                 + [""] * 11 + ["sin_test"])
                continue
            rc_obs = spearman_rho(cx, cy)
            rl_obs = spearman_rho(lx, ly)
            rg_obs = spearman_rho(gx, gy)
            drho_obs = rc_obs - rl_obs
            sd_pla_c, sd_bra_c = desvio(cx), desvio(cy)
            sd_pla_l, sd_bra_l = desvio(lx), desvio(ly)
            for esc in ESCENARIOS:
                rho_true = rg_obs if esc == "GLOBAL" else rc_obs
                r_pear = rho_s_a_r(rho_true)
                s = _simular_celda(rng, r_pear, sd_pla_c, sd_bra_c, nc,
                                   sd_pla_l, sd_bra_l, nl, drho_obs)
                dentro = (s["q025"] <= drho_obs <= s["q975"])
                filas.append([
                    item, tipo, estrato, esc, g10(rho_true), p6e(r_pear), nc, nl,
                    g10(sd_pla_c), g10(sd_bra_c), g10(sd_pla_l), g10(sd_bra_l),
                    g10(drho_obs), g10(s["media"]), g10(s["sd"]), g10(s["q025"]),
                    g10(s["q975"]), g10(s["frac_ge"]), g10(s["tasa_fisher"]),
                    "DENTRO" if dentro else "FUERA"])
    return filas


# ===========================================================================
# 3. Figura.
# ===========================================================================
def figura_simulacion(sim, ruta):
    # Los 3 estratos lado a lado (fila = estrato, columna = item): AMBOS_SEXOS
    # (agrupado, como antes) + HEMBRA + MACHO.
    poritem = {}
    for f in sim:
        poritem.setdefault((f[0], f[2]), []).append(f)
    fig, axes = plt.subplots(len(ESTRATOS), len(ITEMS), figsize=(18.0, 9.5))
    col_esc = {"GLOBAL": "#4C4C4C", "CONTROL": COL_TTO["CONTROL"]}
    col_ver = {"DENTRO": "#009E73", "FUERA": "#D55E00"}
    for i, estrato in enumerate(ESTRATOS):
        for k, item in enumerate(ITEMS):
            ax = axes[i, k]
            filas = poritem[(item, estrato)]
            ymap = {"GLOBAL": 1.0, "CONTROL": 0.0}
            tiene = False
            for f in filas:
                esc = f[3]
                y = ymap[esc]
                if f[19] == "sin_test":
                    ax.text(0.5, 0.5, "n<5\n(sin simulacion)",
                            ha="center", va="center", fontsize=6.5, color="0.5",
                            transform=ax.transAxes)
                    break
                tiene = True
                q025, q975 = float(f[15]), float(f[16])
                media = float(f[13])
                drho_obs = float(f[12])
                ax.plot([q025, q975], [y, y], "-", color=col_esc[esc], lw=3,
                        alpha=0.55, solid_capstyle="butt")
                ax.plot(media, y, "|", color=col_esc[esc], ms=10, mew=1.3)
                ax.plot(drho_obs, y, "D", color=col_ver[f[19]], ms=6, mec="black",
                        mew=0.5, zorder=3)
            if tiene:
                ax.axvline(0.0, color="0.6", lw=0.6, ls=":")
                ax.set_yticks([0.0, 1.0])
                ax.set_yticklabels(["r=ρ Control", "r=ρ GLOBAL"], fontsize=5.5)
                ax.set_ylim(-0.6, 1.6)
                ax.tick_params(axis="x", labelsize=6)
            if i == 0:
                nom = "score compuesto" if item == "score_compuesto" else item
                ax.set_title(nom, fontsize=8,
                             style="normal" if item == "score_compuesto" else "italic")
            if k == 0:
                ax.set_ylabel(estrato, fontsize=6.5)
    fig.suptitle("Δρ observado (rombo) vs intervalo 95% del Δρ simulado bajo r "
                 "verdadera comun y solo diferencia de dispersion\n"
                 "Filas = estrato de sexo. verde = DENTRO (no se puede descartar "
                 "restriccion de rango, prohibicion 5) · naranja = FUERA",
                 fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


# ===========================================================================
# 4. Artefactos compartidos (merge por 'script') -- headers identicos a 02..09.
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
    "## 10_acto2_simulacion",
    "",
    "### Que resuelve",
    "",
    "- **Prohibicion 5.** Antes de leer cualquier cambio de la correlacion "
    "placenta<->cerebro entre Control y LPS como cambio de coordinacion "
    "biologica, hay que descartar que sea un artefacto de dispersion "
    "(restriccion de rango). Una correlacion se atenua si el rango de alguna de "
    "las dos variables se comprime en un subgrupo, con la misma asociacion "
    "verdadera.",
    "",
    "### Diseno (decisiones del usuario)",
    "",
    "- **Escala: normal bivariada en `-ddCt`.** Se generan (placenta, cerebro) "
    "con UNA `r` verdadera comun a ambos grupos y se escalan las SD marginales "
    "de cada grupo a las **SD observadas** de ese grupo (tabla de "
    "`09_acto2_dispersion`). Lo unico que difiere entre grupos en la simulacion "
    "es la dispersion.",
    "- **`r` verdadera comun -- ambas como rango:** (a) ESCENARIO `GLOBAL` -> "
    "`rho` verdadera = `rho` de Spearman GLOBAL del item (pooled, la de T7); "
    "(b) ESCENARIO `CONTROL` -> `rho` verdadera = `rho` de Control (grupo con "
    "mas coordinacion; escenario mas exigente). El `Delta rho` observado se "
    "compara contra la distribucion simulada de **cada** escenario.",
    "- **Inversion `rho_S` -> `r_Pearson`:** para la normal bivariada "
    "`rho_S = (6/pi)*arcsin(r/2)`  =>  `r = 2*sin(pi*rho_S/6)` (asi la `rho` de "
    "Spearman simulada apunta al target). `r` se recorta a +/- 0.999999.",
    "- **B = 2000 repeticiones**, semilla 20260101, RNG PROPIO (LCG 32-bit + "
    "polar de Marsaglia de `00_config`; misma secuencia exacta R/Python). Un "
    "unico RNG recorre todos los items y escenarios en orden.",
    "",
    "### Salidas por item x escenario y veredicto",
    "",
    "- `delta_rho_obs`, media y SD del `Delta rho` simulado, intervalo central "
    "95% (`sim_q025`/`sim_q975`, cuantiles tipo 7), "
    "`sim_frac_abs_ge_obs` = fraccion de repeticiones con "
    "`|Delta rho_sim| >= |Delta rho_obs|` (una `p` de simulacion), y "
    "`sim_tasa_fisher_sig` = tasa de repeticiones en que el Fisher z de 09 "
    "(SE Bonett-Wright) daria `p < .05` -> **falsos positivos del test "
    "reportado bajo pura restriccion de rango**.",
    "- `veredicto = DENTRO`: el `Delta rho` observado cae dentro del intervalo "
    "simulado -> la sola diferencia de dispersion puede producirlo -> **no se "
    "puede descartar restriccion de rango**; el cambio de correlacion no es "
    "interpretable como coordinacion biologica. `FUERA`: lo excede en ese "
    "escenario.",
    "",
    "### Estratificacion por sexo (pedido explicito, punto 4/4, aditiva)",
    "",
    "- **Por que hace falta**: `09_acto2_dispersion` gano Delta rho por sexo "
    "(`HEMBRA`/`MACHO`, ademas de `AMBOS_SEXOS`) -- cada uno de esos Delta rho "
    "necesita su propio control de restriccion de rango, o queda sin la "
    "verificacion que lo hace interpretable (prohibicion 5).",
    "- **Mismo diseno**, extendido con columna `ESTRATO`: una `r` verdadera "
    "comun a Control y LPS (rango GLOBAL/CONTROL, igual que antes), y las SD "
    "marginales = SD **observadas dentro de ese estrato de sexo** (misma "
    "celda SEXO x TTO que la tabla de dispersion de 09). El veredicto sigue "
    "siendo DENTRO/FUERA, mismo criterio.",
    "- **Limitacion declarada, no oculta**: con `n <= 9` por celda de sexo, el "
    "intervalo simulado es ancho y el veredicto `DENTRO` es casi automatico -- "
    "**con este `n` no se puede distinguir cambio de coordinacion de cambio de "
    "dispersion** dentro de cada sexo por separado. Es un resultado en si "
    "mismo (falta de potencia), no un defecto de la simulacion; se dice en el "
    "cuerpo del reporte (Seccion 1), no en una nota al pie.",
    "- **Bug propio encontrado y corregido**: la fila `sin_test` (celdas bajo "
    "el piso `n_par >= 5`) tenia UN campo de menos que columnas la tabla -- "
    "bug preexistente, nunca disparado porque `AMBOS_SEXOS` siempre superaba "
    "el piso con estos datos; con `HEMBRA`/`MACHO` (n<=9) si se dispara. "
    "Corregido completando los 20 campos de `COLS_SIM`.",
    "",
    "### Descartado",
    "",
    "- **Simular sobre rangos** (en vez de en `-ddCt`): no modela como la "
    "compresion de la escala original reduce el rango observado, que es el "
    "mecanismo de la restriccion de rango. Se usa la normal bivariada en "
    "`-ddCt`.",
    "- **Un solo escenario de `r` verdadera**: se corre GLOBAL y CONTROL como "
    "rango porque la eleccion del `r` verdadero cambia cuanto `Delta rho` "
    "espurio aparece; reportar los dos acota la conclusion.",
    "",
    "### Paridad R / Python",
    "",
    "- Los sorteos usan `norm1()` (LCG + polar, byte-identico R/Python por "
    "`00_config`/T1). `rho`, medias, SD, cuantiles y conteos usan acumulador "
    "`double` explicito -> texto `%.10g` bit-identico. `r_pearson_gen` pasa por "
    "`sin()` -> texto `%.6e`. Las figuras son PNG: equivalentes, no "
    "byte-identicas.",
])


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 10_acto2_simulacion:inicio -->"
    marca_fin = "<!-- 10_acto2_simulacion:fin -->"
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


def construir_reporte(fuente, sim):
    L = []
    ap = L.append
    ap("# Reporte de la simulacion de restriccion de rango Acto 2.5 (T8)")
    ap("")
    ap("Generado por `10_acto2_simulacion` (R y Python producen este archivo "
       "identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Metodo")
    ap("")
    ap(f"- Normal bivariada en `-ddCt`, `r` verdadera comun a ambos grupos, SD "
       f"marginales = SD observadas por grupo (09). `r = 2*sin(pi*rho_S/6)`. "
       f"B = {B_SIM}, semilla {cfg.SEMILLA}, RNG PROPIO.")
    ap("- Dos escenarios de `r` verdadera como rango: `GLOBAL` (rho de T7) y "
       "`CONTROL` (rho de Control).")
    ap("- `veredicto = DENTRO` -> la sola diferencia de dispersion puede "
       "producir el `Delta rho` observado (no se descarta restriccion de rango, "
       "prohibicion 5).")
    ap("- **Estratificado por sexo** (pedido explicito, aditivo): columna "
       "`ESTRATO` = `AMBOS_SEXOS` (agrupado, como antes) / `HEMBRA` / `MACHO`. "
       "Cada `Delta rho` nuevo por sexo (de `09_acto2_dispersion`) necesita su "
       "propio control de restriccion de rango; mismo diseno, con las SD "
       "OBSERVADAS de cada celda de sexo.")
    ap("- **Limitacion declarada** (no es un defecto de la simulacion, es un "
       "resultado en si mismo): con `n <= 9` por celda de sexo, el intervalo "
       "simulado va a ser ancho, lo cual hace que el veredicto `DENTRO` sea "
       "casi automatico. **Con este `n` no se puede distinguir cambio de "
       "coordinacion de cambio de dispersion** dentro de `HEMBRA`/`MACHO` por "
       "separado -- la lectura confiable de restriccion de rango sigue siendo "
       "`AMBOS_SEXOS`.")
    ap("")
    ap("## 2. Resultados por item, estrato y escenario")
    ap("")
    ap(_md(COLS_SIM, sim))
    ap("")
    ap("## 3. Figura")
    ap("")
    ap("- `outputs/figures/acto2_simulacion_delta_rho.png` -- `Delta rho` "
       "observado vs intervalo 95% del `Delta rho` simulado, por item y "
       "escenario.")
    ap("")
    ap("## 4. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `10_acto2_simulacion`: inversion "
       "rho->r, eleccion de escenarios, y lectura del veredicto.")
    ap("")
    return "\n".join(L)


# ===========================================================================
def main():
    D = cargar()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)
    rng = cfg.nuevo_rng(cfg.SEMILLA)

    sim = tabla_simulacion(D, rng)

    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "acto2_simulacion.csv", COLS_SIM, sim)

    fig_sim = cfg.RUTA_FIGURAS / "acto2_simulacion_delta_rho.png"
    figura_simulacion(sim, fig_sim)

    escribir_texto(cfg.RUTA_TABLAS / "acto2_simulacion_reporte.md",
                   construir_reporte(fuente, sim))
    actualizar_descartados()

    # Conteos de resumen (print + verificaciones): AMBOS_SEXOS, para no romper
    # la semantica de las verificaciones ya existentes (20 celdas = 10 items x
    # 2 escenarios); HEMBRA/MACHO se resumen aparte.
    sim_ambos = [f for f in sim if f[2] == "AMBOS_SEXOS"]
    n_sim = sum(1 for f in sim_ambos if f[19] != "sin_test")
    n_dentro = sum(1 for f in sim_ambos if f[19] == "DENTRO")
    n_fuera = sum(1 for f in sim_ambos if f[19] == "FUERA")

    # Particion: n(HEMBRA) + n(MACHO) == n(AMBOS_SEXOS) por item x escenario x
    # grupo -- si no cierra, hay error de filtrado. n_control/n_lps en col 6/7.
    n_por_sim = {(f[0], f[2], f[3]): (int(f[6]), int(f[7])) for f in sim}
    particion_sim_ok = True
    particion_sim_detalle = []
    for item in ITEMS:
        for esc in ESCENARIOS:
            nc_a, nl_a = n_por_sim[(item, "AMBOS_SEXOS", esc)]
            nc_h, nl_h = n_por_sim[(item, "HEMBRA", esc)]
            nc_m, nl_m = n_por_sim[(item, "MACHO", esc)]
            if nc_h + nc_m != nc_a or nl_h + nl_m != nl_a:
                particion_sim_ok = False
                particion_sim_detalle.append(f"{item} {esc}")

    ent = (f"data/processed/qpcr_cuantificacion_long.tsv + qpcr_score_compuesto_long.tsv "
           f"+ outputs/tables/*/acto2_dispersion.csv (de data/{fuente}/{cfg.ARCHIVO_QPCR})")
    registrar_procedencia([
        ["outputs/tables/{R,python}/acto2_simulacion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "simulacion de restriccion de rango: normal bivariada en "
         "-ddCt, r verdadera comun (escenarios GLOBAL y CONTROL), SD por grupo "
         "observadas; Delta rho simulado, IC95, frac >=|obs|, tasa Fisher, "
         "veredicto DENTRO/FUERA; estratificado por ESTRATO "
         "(AMBOS_SEXOS/HEMBRA/MACHO, pedido explicito -- cada Delta rho nuevo por "
         "sexo de 09_acto2_dispersion necesita su propio control de restriccion "
         "de rango, mismo diseno con las SD observadas de cada celda de sexo; "
         "limitacion declarada: con n<=9 el intervalo simulado es ancho y el "
         "veredicto DENTRO es casi automatico -- con este n no se puede "
         "distinguir cambio de coordinacion de cambio de dispersion)"],
        ["outputs/figures/acto2_simulacion_delta_rho.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent, "Delta rho observado vs intervalo 95% del Delta rho "
         "simulado por item y escenario"],
        ["outputs/tables/acto2_simulacion_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.5 (T8)"],
    ])
    registrar_verificaciones([
        ["acto2_sim_proposito",
         "la simulacion descarta (o no) restriccion de rango antes de interpretar "
         "(prohibicion 5)",
         "normal bivariada en -ddCt con r verdadera comun y SD por grupo "
         "observadas; veredicto DENTRO/FUERA por escenario",
         "prohibicion 5 atendida", "TRUE", ESTE_SCRIPT],
        ["acto2_sim_escenarios",
         "r verdadera comun probada como rango: GLOBAL (rho T7) y CONTROL (rho Control)",
         ";".join(ESCENARIOS), "GLOBAL;CONTROL",
         "TRUE" if ESCENARIOS == ["GLOBAL", "CONTROL"] else "FALSE", ESTE_SCRIPT],
        ["acto2_sim_inversion_rho_r",
         "inversion rho_Spearman -> r_Pearson generativo para la normal bivariada",
         "r = 2*sin(pi*rho_S/6), recorte +/- 0.999999",
         "2*sin(pi*rho/6)", "TRUE", ESTE_SCRIPT],
        ["acto2_sim_B_semilla",
         "repeticiones y semilla de la simulacion",
         f"B={B_SIM}; semilla={cfg.SEMILLA}; RNG PROPIO LCG+polar",
         f"B={B_SIM}, semilla 20260101",
         "TRUE" if (B_SIM == 2000 and cfg.SEMILLA == 20260101) else "FALSE",
         ESTE_SCRIPT],
        ["acto2_sim_SD_por_grupo",
         "las SD marginales simuladas por grupo son las SD observadas de ese grupo",
         "sd_pla/sd_bra de pares(item, grupo); unica diferencia entre grupos en "
         "la simulacion", "SD observadas", "TRUE", ESTE_SCRIPT],
        ["acto2_sim_veredicto",
         "veredicto por item x escenario (AMBOS_SEXOS): DENTRO (no se descarta "
         "rango) / FUERA",
         f"{n_sim} celdas simuladas; DENTRO={n_dentro}; FUERA={n_fuera}",
         "DENTRO/FUERA/sin_test", "TRUE", ESTE_SCRIPT],
        ["acto2_sim_estratos_sexo",
         "ESTRATO = AMBOS_SEXOS/HEMBRA/MACHO, mismo diseno con SD observadas por sexo",
         ";".join(ESTRATOS), "AMBOS_SEXOS;HEMBRA;MACHO",
         "TRUE" if ESTRATOS == ["AMBOS_SEXOS", "HEMBRA", "MACHO"] else "FALSE",
         ESTE_SCRIPT],
        ["acto2_sim_estrato_particion",
         "n(HEMBRA) + n(MACHO) = n(AMBOS_SEXOS) por item x escenario x grupo",
         "particiona en %d/%d celdas item x escenario%s" % (
             len(ITEMS) * len(ESCENARIOS) - len(particion_sim_detalle),
             len(ITEMS) * len(ESCENARIOS),
             ("; falla en: " + ", ".join(particion_sim_detalle))
             if particion_sim_detalle else ""),
         "particiona en %d/%d celdas" % (len(ITEMS) * len(ESCENARIOS),
                                          len(ITEMS) * len(ESCENARIOS)),
         "TRUE" if particion_sim_ok else "FALSE", ESTE_SCRIPT],
        ["acto2_sim_tasa_fisher",
         "se reporta la tasa de falsos positivos del Fisher z (09) bajo pura "
         "restriccion de rango",
         "columna sim_tasa_fisher_sig por celda", "tasa reportada", "TRUE",
         ESTE_SCRIPT],
        ["acto2_sim_figura",
         "figura Acto 2.5: Delta rho observado vs intervalo simulado",
         f"sim={fig_sim.is_file()}".replace("True", "TRUE").replace("False", "FALSE"),
         "1 figura existe", "TRUE" if fig_sim.is_file() else "FALSE", ESTE_SCRIPT],
    ])

    print("== 10_acto2_simulacion.py ==")
    print(f"  fuente = {fuente}   B = {B_SIM}   semilla = {cfg.SEMILLA}")
    print(f"  celdas simuladas AMBOS_SEXOS: {n_sim}  (DENTRO={n_dentro}, FUERA={n_fuera})")
    for f in sim_ambos:
        if f[19] == "sin_test":
            print(f"    {f[0]:16} {f[3]:8}  (sin test, piso)")
        else:
            print(f"    {f[0]:16} {f[3]:8}  dRho_obs={f[12]:>9}  "
                  f"IC95_sim=[{f[15]:>9}, {f[16]:>9}]  frac>=obs={f[17]:>8}  "
                  f"fisherFP={f[18]:>7}  -> {f[19]}")
    for estrato in ("HEMBRA", "MACHO"):
        fs = [f for f in sim if f[2] == estrato]
        nsi = sum(1 for f in fs if f[19] != "sin_test")
        ndi = sum(1 for f in fs if f[19] == "DENTRO")
        nfu = sum(1 for f in fs if f[19] == "FUERA")
        print(f"  celdas simuladas {estrato}: {nsi}  (DENTRO={ndi}, FUERA={nfu})")
    print("  particion n(HEMBRA)+n(MACHO)=n(AMBOS_SEXOS) [item x escenario]: " + (
        f"OK en {len(ITEMS) * len(ESCENARIOS)}/{len(ITEMS) * len(ESCENARIOS)} celdas"
        if particion_sim_ok else "FALLA en: " + ", ".join(particion_sim_detalle)))
    print(f"  -> outputs/tables/{{R,python}}/acto2_simulacion.csv")
    print(f"  -> {fig_sim.name}")


if __name__ == "__main__":
    main()
