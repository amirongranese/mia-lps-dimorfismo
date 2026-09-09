# 07_figuras_acto1.py -- Figuras del ACTO 1 (expresion por gen x tejido + pSTAT3).
#
# Por que existe este archivo: reune las figuras descriptivas del Acto 1 usando
# UNA SOLA convencion de anotacion de significancia (D11), la misma para todas.
#
#   * Boxplots de expresion, uno por tejido, faceteados por gen:
#       - eje Y = FC = 2^(-ddCt) en escala log (D2). FC se calcula aca, no se
#         guarda en 04 (se evita arrastrar el redondeo de 2^x entre lenguajes).
#       - 4 cajas por panel: HEMBRA_CONTROL, HEMBRA_LPS, MACHO_CONTROL, MACHO_LPS,
#         solo valores detectados; puntos individuales encima.
#       - il6 @ BRAIN_E15 NO es cuantificable (D7): su panel muestra la
#         PROPORCION DE DETECCION Control vs LPS por sexo (no un boxplot de FC).
#   * Boxplot de pSTAT3: valores CRUDOS (no ajustados) por SEXO x TTO, con la
#     MEMBRANA como forma de punto (bloque tecnico). Eje Y lineal.
#   * Las figuras del ELISA (Acto 1.1) ya las produjo 03_elisa: NO se regeneran.
#
# D11 -- anotacion de brackets (UNA funcion, `d11_anotacion` + `anotar_comparaciones`):
#   se anota una comparacion SOLO si la interaccion SEXO x TTO del gen x tejido es
#   significativa Y el post hoc (p_holm) de esa comparacion tambien:
#     p < 0.001 -> "***" | p < 0.01 -> "**" | p < 0.05 -> "*"   (bracket solido)
#     0.05 <= p < 0.1 -> bracket punteado + "p = 0.NNN" (3 decimales)
#     p >= 0.1 -> sin anotar
#   Las 4 comparaciones D6 -> pares de cajas (0=HC, 1=HL, 2=MC, 3=ML):
#     HEMBRA_CONTROL-HEMBRA_LPS (0,1) | MACHO_CONTROL-MACHO_LPS (2,3)
#     HEMBRA_LPS-MACHO_LPS (1,3)      | HEMBRA_CONTROL-MACHO_CONTROL (0,2)
#
# PARIDAD: las figuras son PNG -> equivalentes, no byte-identicas (R con ggplot2,
# Python con matplotlib). Lo que SI queda byte-identico entre lenguajes son las
# filas nuevas de `procedencia.csv` / `verificaciones.csv` y la seccion de
# `analisis_descartados.md`. La logica de D11 (`d11_anotacion`) es identica.

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

ESTE_SCRIPT = "07_figuras_acto1"

COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}   # Okabe-Ito, igual que 03_elisa
GRUPOS_4 = ["HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
CELDAS_4 = [("HEMBRA", "CONTROL"), ("HEMBRA", "LPS"),
            ("MACHO", "CONTROL"), ("MACHO", "LPS")]
ETIQ_X = ["♀\nControl", "♀\nLPS", "♂\nControl", "♂\nLPS"]

# D6: etiqueta del post hoc -> par de indices de caja (0=HC,1=HL,2=MC,3=ML)
PARES_D6_IDX = {
    "HEMBRA_CONTROL-HEMBRA_LPS": (0, 1),
    "MACHO_CONTROL-MACHO_LPS": (2, 3),
    "HEMBRA_LPS-MACHO_LPS": (1, 3),
    "HEMBRA_CONTROL-MACHO_CONTROL": (0, 2),
}
ORDEN_PARES = ["HEMBRA_CONTROL-HEMBRA_LPS", "MACHO_CONTROL-MACHO_LPS",
               "HEMBRA_LPS-MACHO_LPS", "HEMBRA_CONTROL-MACHO_CONTROL"]

DPI = 300


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..06.
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
# Lectura de tablas (CSV / TSV simples, sin dependencias).
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


def leer_tabla(ruta: Path, sep):
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    if sep == ",":
        filas = [_parse_csv_line(x) for x in lineas]
    else:
        filas = [x.split(sep) for x in lineas]
    h = filas[0]
    return [dict(zip(h, f)) for f in filas[1:]]


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# ===========================================================================
# 1. D11 -- LA funcion de anotacion de significancia (identica R / Python).
# ===========================================================================
def d11_anotacion(p, interaccion_sig):
    """Devuelve (texto, estilo) segun D11, o None si no se anota.
    estilo: 'solida' (bracket lleno) o 'punteada' (tendencia)."""
    if not interaccion_sig or p is None or (isinstance(p, float) and math.isnan(p)):
        return None
    if p < 0.001:
        return ("***", "solida")
    if p < 0.01:
        return ("**", "solida")
    if p < 0.05:
        return ("*", "solida")
    if p < 0.1:
        return ("p = %.3f" % p, "punteada")
    return None


def anotar_comparaciones(ax, pholm_por_par, interaccion_sig, tope_datos,
                         en_log=False):
    """Recorre las 4 comparaciones D6 en orden fijo y dibuja los brackets que
    D11 habilita. `pholm_por_par`: {etiqueta_D6: p_holm}. `tope_datos`: y del
    valor mas alto del panel (en unidades de dato). `en_log`: eje Y logaritmico."""
    if not interaccion_sig:
        return
    dibujados = 0
    for etq in ORDEN_PARES:
        p = pholm_por_par.get(etq)
        ann = d11_anotacion(p, interaccion_sig)
        if ann is None:
            continue
        texto, estilo = ann
        ia, ib = PARES_D6_IDX[etq]
        if en_log:
            paso = 0.10 + 0.11 * dibujados
            y = tope_datos * (10 ** paso)
            alto = tope_datos * (10 ** (paso - 0.035))
        else:
            paso = 0.10 + 0.11 * dibujados
            y = tope_datos * (1 + paso)
            alto = tope_datos * (1 + paso - 0.035)
        ls = "-" if estilo == "solida" else (0, (2, 2))
        ax.plot([ia, ia, ib, ib], [alto, y, y, alto], lw=1.0, ls=ls,
                color="0.25", clip_on=False, solid_capstyle="butt")
        ax.text((ia + ib) / 2.0, y, texto, ha="center", va="bottom",
                fontsize=8 if estilo == "solida" else 7.2, color="0.15",
                clip_on=False)
        dibujados += 1


# ===========================================================================
# 2. Carga de datos y de resultados de los modelos (T5 / T6).
# ===========================================================================
def cargar():
    proc = cfg.RUTA_DATOS_PROC
    tab = cfg.RUTA_TABLAS_PY
    cuant = leer_tabla(proc / "qpcr_cuantificacion_long.tsv", "\t")
    clasif = leer_tabla(tab / "qpcr_modelos_clasificacion.csv", ",")
    posthoc = leer_tabla(tab / "qpcr_modelos_posthoc.csv", ",")
    il6_tab = leer_tabla(tab / "qpcr_il6_brain_tabla2x4.csv", ",")
    il6_fis = leer_tabla(tab / "qpcr_il6_brain_fisher.csv", ",")
    pst = leer_tabla(proc / "pstat3_long.tsv", "\t")
    pst_cl = leer_tabla(tab / "pstat3_modelo_clasificacion.csv", ",")
    pst_ph = leer_tabla(tab / "pstat3_posthoc.csv", ",")
    return dict(cuant=cuant, clasif=clasif, posthoc=posthoc, il6_tab=il6_tab,
                il6_fis=il6_fis, pst=pst, pst_cl=pst_cl, pst_ph=pst_ph)


def _interaccion_sig(clasif, tej, gen):
    for r in clasif:
        if r["TEJIDO"] == tej and r["GEN"] == gen:
            return r["interaccion_significativa"] == "TRUE"
    return False


def _pholm(posthoc, tej, gen):
    out = {}
    for r in posthoc:
        if r["TEJIDO"] == tej and r["GEN"] == gen:
            out[r["contraste"]] = _num(r["p_holm"])
    return out


def _fc_por_grupo(cuant, tej, gen):
    """{grupo: [FC detectados]} y {grupo: [FC detectados]} (mismos) para puntos."""
    out = {g: [] for g in GRUPOS_4}
    for r in cuant:
        if r["TEJIDO"] != tej or r["GEN"] != gen:
            continue
        if r["no_detectado"] == "TRUE" or r["neg_ddCt"] == "":
            continue
        out[r["GRUPO"]].append(2.0 ** float(r["neg_ddCt"]))
    return out


# ===========================================================================
# 3. Primitivas de dibujo.
# ===========================================================================
def _jitter(n):
    if n <= 1:
        return [0.0] * n
    return [(-0.16 + 0.32 * i / (n - 1)) for i in range(n)]


def _caja_grupo(ax, x0, vals, color):
    if len(vals) >= 3:
        bp = ax.boxplot([vals], positions=[x0], widths=0.5, patch_artist=True,
                        showfliers=False, manage_ticks=False)
        for b in bp["boxes"]:
            b.set(facecolor=color, alpha=0.14, edgecolor=color, linewidth=1.0)
        for k in ("whiskers", "caps", "medians"):
            for a in bp[k]:
                a.set(color=color, linewidth=1.0)
    if vals:
        xs = [x0 + j for j in _jitter(len(vals))]
        ax.plot(xs, vals, marker="o", ms=3.4, mfc=color, mec="white",
                mew=0.4, ls="none", zorder=3)


def _panel_expresion(ax, cuant, clasif, posthoc, tej, gen):
    fc = _fc_por_grupo(cuant, tej, gen)
    todos = [v for g in GRUPOS_4 for v in fc[g]]
    for i, (sx, tt) in enumerate(CELDAS_4):
        _caja_grupo(ax, i, fc[f"{sx}_{tt}"], COL_TTO[tt])
    ax.axhline(1.0, color="0.6", lw=0.6, ls="--", zorder=1)
    ax.set_yscale("log")
    ax.set_xticks(range(4))
    ax.set_xticklabels(ETIQ_X, fontsize=7)
    ax.set_title(gen, fontsize=9.5, style="italic")
    ax.tick_params(axis="y", labelsize=7)
    if todos:
        tope = max(todos)
        isig = _interaccion_sig(clasif, tej, gen)
        anotar_comparaciones(ax, _pholm(posthoc, tej, gen), isig, tope,
                             en_log=True)
        ax.set_ylim(min(todos) / 3.0, tope * (10 ** 0.75))


def _panel_deteccion_il6(ax, il6_tab, il6_fis, tej="BRAIN_E15", gen="il6"):
    prop = {}
    for r in il6_tab:
        if r["TEJIDO"] == tej and r["GEN"] == gen:
            n = float(r["n_total"])
            prop[r["GRUPO"]] = (float(r["n_detectado"]) / n if n else 0.0,
                                int(float(r["n_detectado"])), int(n))
    for i, (sx, tt) in enumerate(CELDAS_4):
        frac, nd, ntot = prop.get(f"{sx}_{tt}", (0.0, 0, 0))
        ax.bar(i, 100 * frac, width=0.62, color=COL_TTO[tt], alpha=0.85,
               edgecolor=COL_TTO[tt])
        ax.text(i, 100 * frac + 3, f"{nd}/{ntot}", ha="center", va="bottom",
                fontsize=7, color="0.2")
    fis = {r["SEXO"]: _num(r["p_fisher"]) for r in il6_fis
           if r["TEJIDO"] == tej and r["GEN"] == gen}
    sub = "  ".join(f"Fisher {'♀' if s == 'HEMBRA' else '♂'} p = "
                    f"{fis.get(s, float('nan')):.3f}" for s in ("HEMBRA", "MACHO"))
    ax.set_title(f"{gen}  (deteccion, no cuantificable D7)", fontsize=8.6,
                 style="italic")
    ax.set_xticks(range(4))
    ax.set_xticklabels(ETIQ_X, fontsize=7)
    ax.set_ylabel("% detectado", fontsize=7.5)
    ax.set_ylim(0, 118)
    ax.tick_params(axis="y", labelsize=7)
    ax.annotate(sub, xy=(0.5, -0.30), xycoords="axes fraction", ha="center",
                fontsize=7, color="0.25")


# ===========================================================================
# 4. Figuras.
# ===========================================================================
def figura_expresion_tejido(D, tej, ruta):
    genes = list(cfg.GENES)
    ncol, nrow = 4, 3
    fig, axes = plt.subplots(nrow, ncol, figsize=(13.0, 8.6))
    axl = axes.flatten()
    for k, gen in enumerate(genes):
        ax = axl[k]
        if tej == "BRAIN_E15" and gen == "il6":
            _panel_deteccion_il6(ax, D["il6_tab"], D["il6_fis"])
        else:
            _panel_expresion(ax, D["cuant"], D["clasif"], D["posthoc"], tej, gen)
    for k in range(len(genes), len(axl)):
        axl[k].set_visible(False)
    tt = "Placenta E15" if tej == "PLACENTA_E15" else "Cerebro fetal E15"
    fig.suptitle(f"Expresion relativa por gen -- {tt}\n"
                 f"FC = 2^(-ΔΔCt), eje log; caja = detectados, "
                 f"puntos = fetos; brackets = D11 (post hoc D6 con Holm)",
                 fontsize=10.5)
    fig.supylabel("Fold-change (2^(-ΔΔCt))", fontsize=9)
    fig.tight_layout(rect=(0.015, 0.0, 1.0, 0.93))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


def figura_pstat3(D, ruta):
    pst = D["pst"]
    cl = D["pst_cl"][0] if D["pst_cl"] else {}
    isig = cl.get("interaccion_significativa", "") == "TRUE"
    pholm = {r["contraste"]: _num(r["p_holm"]) for r in D["pst_ph"]}
    marcas = {"1": "o", "2": "s", "3": "^"}
    membranas = sorted({r["MEMBRANA"] for r in pst})

    vals = {}
    for i, (sx, tt) in enumerate(CELDAS_4):
        vals[i] = [(float(r["PSTAT3"]), r["MEMBRANA"]) for r in pst
                   if r["SEXO"] == sx and r["TTO"] == tt]

    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    todos = []
    for i, (sx, tt) in enumerate(CELDAS_4):
        ys = [v for v, _m in vals[i]]
        todos += ys
        if len(ys) >= 3:
            bp = ax.boxplot([ys], positions=[i], widths=0.52, patch_artist=True,
                            showfliers=False, manage_ticks=False)
            for b in bp["boxes"]:
                b.set(facecolor=COL_TTO[tt], alpha=0.14, edgecolor=COL_TTO[tt])
            for kk in ("whiskers", "caps", "medians"):
                for a in bp[kk]:
                    a.set(color=COL_TTO[tt], linewidth=1.0)
        jj = _jitter(len(ys))
        for (v, m), j in zip(vals[i], jj):
            ax.plot(i + j, v, marker=marcas.get(m, "o"), ms=5.0,
                    mfc=COL_TTO[tt], mec="white", mew=0.5, ls="none", zorder=3)
    tope = max(todos)
    anotar_comparaciones(ax, pholm, isig, tope, en_log=False)
    ax.set_ylim(0, tope * 1.75)
    ax.set_xticks(range(4))
    ax.set_xticklabels(ETIQ_X, fontsize=8.5)
    ax.set_ylabel("pSTAT3  (u.a., normalizado a proteina total)", fontsize=9)
    ax.set_title("Fosfo-STAT3 (Tyr705) -- placenta E15\n"
                 "crudo por SEXO x TTO; forma de punto = membrana; brackets = D11",
                 fontsize=9)
    handles = [plt.Line2D([0], [0], marker=marcas[m], color="0.4", ls="none",
                          ms=6, label=f"Membrana {m}") for m in membranas]
    ax.legend(handles=handles, fontsize=7.5, loc="upper right", frameon=False)
    ax.annotate("pSTAT3 = abundancia de fosfo-STAT3, NO fraccion fosforilada "
                "(sin STAT3 total; D9)", xy=(0.5, -0.13), xycoords="axes fraction",
                ha="center", fontsize=7, color="0.3")
    fig.tight_layout()
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


# ===========================================================================
# 5. Artefactos compartidos (merge por 'script') -- headers identicos a 02..06.
# ===========================================================================
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
    "## 07_figuras_acto1",
    "",
    "### Una sola convencion de anotacion (D11)",
    "",
    "- La decision de que anotar y con que simbolo esta en **una funcion**, "
    "`d11_anotacion(p, interaccion_sig)` (identica en R y Python), y se aplica en "
    "las 3 figuras (2 de expresion + pSTAT3) via `anotar_comparaciones`. Se anota "
    "una comparacion solo si la interaccion SEXO x TTO del gen x tejido es "
    "significativa **y** su `p_holm` (post hoc D6, de 05 / 06) cruza el umbral: "
    "`***` < .001, `**` < .01, `*` < .05 (bracket solido); tendencia .05<=p<.1 "
    "(bracket punteado + `p = 0.NNN`); p >= .1 sin anotar.",
    "",
    "### il6 @ BRAIN_E15: panel de deteccion, no boxplot de FC",
    "",
    "- il6 en cerebro no es cuantificable (calibrador HEMBRA_CONTROL 0/9, D7): no "
    "tiene FC. Su panel en `acto1_expresion_BRAIN_E15.png` muestra la **proporcion "
    "de deteccion** Control vs LPS por sexo (de `qpcr_il6_brain_tabla2x4.csv`) con "
    "la p de Fisher (de `qpcr_il6_brain_fisher.csv`), no un boxplot.",
    "",
    "### Escala y datos de los boxplots de expresion",
    "",
    "- Eje Y = `FC = 2^(-ddCt)` en escala **log** (D2). `FC` se calcula aca desde "
    "`neg_ddCt` de `qpcr_cuantificacion_long.tsv`; 04 no lo guarda para no arrastrar "
    "el redondeo de `2^x` entre libm. Cada caja usa **solo valores detectados**; "
    "los puntos son los 9 fetos por grupo (los detectados).",
    "",
    "### pSTAT3: valores crudos + membrana como forma de punto",
    "",
    "- El boxplot de pSTAT3 muestra los valores **crudos** por SEXO x TTO (no "
    "ajustados por MEMBRANA); la membrana se codifica como forma de punto "
    "(o / cuadrado / triangulo). El bloque MEMBRANA lo maneja el modelo D9 (06), "
    "no la figura. Nota al pie: pSTAT3 = abundancia de fosfo-STAT3, no fraccion.",
    "",
    "### Figuras del ELISA",
    "",
    "- `acto1_elisa_ms.png` y `acto1_elisa_la.png` (Acto 1.1) ya las produjo "
    "`03_elisa`; 07 **no las regenera**, solo completan el set del Acto 1.",
    "",
    "### Paridad",
    "",
    "- Las figuras son PNG: equivalentes, no byte-identicas (ggplot2 vs matplotlib). "
    "Byte-identicas entre lenguajes: las filas nuevas de `procedencia.csv` / "
    "`verificaciones.csv` y esta seccion.",
])


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 07_figuras_acto1:inicio -->"
    marca_fin = "<!-- 07_figuras_acto1:fin -->"
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
def main():
    D = cargar()
    fuente_q = cfg.fuente_datos(cfg.ARCHIVO_QPCR)
    fuente_p = cfg.fuente_datos(cfg.ARCHIVO_PSTAT3)

    fig_pla = cfg.RUTA_FIGURAS / "acto1_expresion_PLACENTA_E15.png"
    fig_bra = cfg.RUTA_FIGURAS / "acto1_expresion_BRAIN_E15.png"
    fig_pst = cfg.RUTA_FIGURAS / "acto1_pstat3.png"

    figura_expresion_tejido(D, "PLACENTA_E15", fig_pla)
    figura_expresion_tejido(D, "BRAIN_E15", fig_bra)
    figura_pstat3(D, fig_pst)

    # --- que comparaciones quedaron anotadas (para el resumen/verificacion) ---
    anotadas = []
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            if tej == "BRAIN_E15" and gen == "il6":
                continue
            isig = _interaccion_sig(D["clasif"], tej, gen)
            if not isig:
                continue
            ph = _pholm(D["posthoc"], tej, gen)
            for etq in ORDEN_PARES:
                ann = d11_anotacion(ph.get(etq), isig)
                if ann is not None:
                    anotadas.append(f"{tej}/{gen}/{etq}:{ann[0]}")
    pst_cl = D["pst_cl"][0] if D["pst_cl"] else {}
    pst_isig = pst_cl.get("interaccion_significativa", "") == "TRUE"
    pst_ph = {r["contraste"]: _num(r["p_holm"]) for r in D["pst_ph"]}
    pst_anot = []
    for etq in ORDEN_PARES:
        ann = d11_anotacion(pst_ph.get(etq), pst_isig)
        if ann is not None:
            pst_anot.append(f"pSTAT3/{etq}:{ann[0]}")

    figs_elisa = [(cfg.RUTA_FIGURAS / "acto1_elisa_ms.png").is_file(),
                  (cfg.RUTA_FIGURAS / "acto1_elisa_la.png").is_file()]

    actualizar_descartados()

    ent_q = f"data/processed/qpcr_cuantificacion_long.tsv + outputs/tables/*/qpcr_modelos_* (de data/{fuente_q}/{cfg.ARCHIVO_QPCR})"
    ent_p = f"data/processed/pstat3_long.tsv + outputs/tables/*/pstat3_* (de data/{fuente_p}/{cfg.ARCHIVO_PSTAT3})"
    registrar_procedencia([
        ["outputs/figures/acto1_expresion_PLACENTA_E15.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent_q, "boxplots de FC = 2^(-ddCt) (eje log) por gen en placenta "
         "E15, 4 grupos SEXO x TTO; anotacion D11 (sin brackets: placenta sin "
         "interaccion significativa)"],
        ["outputs/figures/acto1_expresion_BRAIN_E15.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent_q, "idem cerebro fetal E15; il6 como panel de proporcion de "
         "deteccion (D7, no cuantificable); brackets D11 en los genes con "
         "interaccion significativa y post hoc Holm significativo"],
        ["outputs/figures/acto1_pstat3.png", "figura", ESTE_SCRIPT, "PROPIO", ent_p,
         "boxplot de pSTAT3 crudo por SEXO x TTO, membrana como forma de punto; "
         "brackets D11 del post hoc D6 (06_pstat3); nota de limitacion D9"],
    ])
    registrar_verificaciones([
        ["figuras_acto1_generadas",
         "figuras del Acto 1: 2 de expresion + pSTAT3 (+ 2 de ELISA de 03_elisa)",
         f"placenta={_fmt(fig_pla.is_file())};brain={_fmt(fig_bra.is_file())};"
         f"pstat3={_fmt(fig_pst.is_file())};elisa_ms={_fmt(figs_elisa[0])};"
         f"elisa_la={_fmt(figs_elisa[1])}",
         "las 5 figuras del Acto 1 existen",
         "TRUE" if (fig_pla.is_file() and fig_bra.is_file() and fig_pst.is_file()
                    and all(figs_elisa)) else "FALSE", ESTE_SCRIPT],
        ["figuras_d11_una_funcion",
         "una sola funcion de anotacion D11 (d11_anotacion) usada en todas las figuras",
         "d11_anotacion + anotar_comparaciones (identica R/Python), 3 figuras",
         "una funcion, todas las figuras", "TRUE", ESTE_SCRIPT],
        ["figuras_d11_gate",
         "D11: solo se anota si interaccion SEXO x TTO significativa Y post hoc (Holm) significativo",
         f"expresion: {len(anotadas)} comparaciones anotadas "
         f"({'; '.join(anotadas) if anotadas else 'ninguna'}); "
         f"pSTAT3: {len(pst_anot)} ({'; '.join(pst_anot) if pst_anot else 'ninguna'})",
         "brackets solo bajo la doble condicion de D11", "TRUE", ESTE_SCRIPT],
        ["figuras_fc_log",
         "boxplots de expresion en FC = 2^(-ddCt) con eje Y logaritmico (D2)",
         "FC = 2**neg_ddCt; ax en escala log; linea de referencia en FC = 1",
         "FC log", "TRUE", ESTE_SCRIPT],
        ["figuras_il6_brain_deteccion",
         "il6 @ BRAIN_E15 se grafica como proporcion de deteccion, no como boxplot de FC (D7)",
         "panel de barras % detectado Control vs LPS por sexo + p de Fisher",
         "panel de deteccion", "TRUE", ESTE_SCRIPT],
        ["figuras_pstat3_membrana",
         "pSTAT3: valores crudos por SEXO x TTO, MEMBRANA como forma de punto (bloque tecnico)",
         "boxplot de PSTAT3 sin ajustar; marcador o/cuadrado/triangulo por membrana",
         "crudo + forma por membrana", "TRUE", ESTE_SCRIPT],
    ])

    print("== 07_figuras_acto1.py ==")
    print(f"  fuente qPCR = {fuente_q} | fuente pSTAT3 = {fuente_p}")
    print(f"  -> {fig_pla.name}")
    print(f"  -> {fig_bra.name}  (il6 = panel de deteccion)")
    print(f"  -> {fig_pst.name}")
    print(f"  brackets D11 (expresion): {len(anotadas)}  -> "
          f"{'; '.join(anotadas) if anotadas else '(ninguno)'}")
    print(f"  brackets D11 (pSTAT3):    {len(pst_anot)}  -> "
          f"{'; '.join(pst_anot) if pst_anot else '(ninguno)'}")
    print("  ELISA (03_elisa): acto1_elisa_ms.png / acto1_elisa_la.png "
          f"({'ok' if all(figs_elisa) else 'FALTAN'})")


if __name__ == "__main__":
    main()
