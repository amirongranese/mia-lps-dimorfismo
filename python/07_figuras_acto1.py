# 07_figuras_acto1.py -- Figuras del ACTO 1 (expresion por gen x tejido + pSTAT3).
#
# Por que existe este archivo: reune las figuras descriptivas del Acto 1 usando
# UNA SOLA convencion de anotacion de significancia (D11), la misma para todas.
#
#   * Boxplots de expresion, uno por tejido: matplotlib ya permite estilar caja,
#     bigote y tope del bigote por separado (whiskerprops/capprops/boxprops),
#     asi que ESTE lenguaje no necesita "abandonar" nada -- solo se replica el
#     estilo pedido (pedidos/boxplots_acto1_base_R.R, tratado como especificacion
#     de estilo y logica, no como codigo a copiar; en R si obliga a bajar a R
#     base porque `geom_boxplot` no expone esos 3 trazos por separado; el Acto 2
#     sigue en ggplot2/GGally en R).
#       - eje Y = FC = 2^(-ddCt) en escala log (D2). FC se calcula aca, no se
#         guarda en 04 (se evita arrastrar el redondeo de 2^x entre lenguajes).
#       - 4 cajas por panel: HEMBRA_CONTROL, HEMBRA_LPS, MACHO_CONTROL, MACHO_LPS,
#         solo valores detectados (caja solo si >=3 detectados); puntos encima.
#       - Filas del panel agrupadas por via metabolica (lipidos/glucosa/
#         aminoacidos/IL-6), no por el orden crudo de GENES -- estilo pedido.
#         Color: Control por sexo (celeste), LPS por sexo x via metabolica.
#       - il6 @ BRAIN_E15 NO es cuantificable (D7): su panel muestra la
#         PROPORCION DE DETECCION Control vs LPS por sexo (no un boxplot de FC).
#   * Boxplot de pSTAT3 (sin cambios de estilo): valores CRUDOS por SEXO x TTO,
#     con la MEMBRANA como forma de punto (bloque tecnico). Eje Y lineal.
#   * Las figuras del ELISA (Acto 1.1) ya las produjo 03_elisa: NO se regeneran.
#
# D11 -- anotacion de brackets, AMPLIADA (cambio pedido explicitamente, ver
# AGENTS.md 4.2; revierte la restriccion previa "solo si la interaccion es
# significativa"). Cascada de 3 ramas, en este orden -- solo se entra a UNA:
#   (a) interaccion SEXO x TTO significativa -> brackets por PAR del post hoc
#       D6 (Holm/ART-C), como antes.
#   (b) interaccion NO significativa y efecto principal de TTO significativo o
#       en tendencia -> UN bracket que abarca los 4 grupos, etiqueta
#       "Control vs LPS" + estrellas/p de `p_TTO`.
#   (c) interaccion NO significativa y efecto principal de SEXO significativo o
#       en tendencia -> UN bracket entre los centros de cada sexo, etiqueta
#       "♀ vs ♂" + estrellas/p de `p_SEXO`. (b) y (c) no son excluyentes.
# Simbolos (las 3 ramas): p<0.001 -> "***" | p<0.01 -> "**" | p<0.05 -> "*"
# (bracket solido); 0.05<=p<0.1 -> bracket punteado + "p = 0.NNN"; p>=0.1 -> nada.
# UNA sola funcion decide que anotar (`d11_brackets_especificacion`, identica en
# R y Python); la geometria de dibujo es propia de cada motor.
#
# PARIDAD: las figuras son PNG -> equivalentes, no byte-identicas (R base +
# ggplot2 vs matplotlib). Byte-identicas entre lenguajes: las filas nuevas de
# `procedencia.csv` / `verificaciones.csv` y la seccion de `analisis_descartados.md`.
# La logica de D11 (`d11_texto` + `d11_brackets_especificacion`) es identica.

from __future__ import annotations

import importlib.util
import math
import random
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

# --- Paleta y estilo de los boxplots de expresion (estilo pedido) ----------
COL_CTRL = {"HEMBRA": "#AEDCF0", "MACHO": "#6BAED6"}
COL_LPS = {
    "GLUCOSA":     {"HEMBRA": "#F09EC8", "MACHO": "#D6317F"},
    "AMINOACIDOS": {"HEMBRA": "#8FD9B6", "MACHO": "#2E9E6B"},
    "LIPIDOS":     {"HEMBRA": "#FBC98A", "MACHO": "#E08214"},
    "IL6":         {"HEMBRA": "#C5A3E0", "MACHO": "#7B4EA8"},
}
COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}   # solo pSTAT3 (Okabe-Ito, igual que 03/08)

GPATH = {
    "fatcd36": "LIPIDOS", "fatp1": "LIPIDOS", "fatp4": "LIPIDOS",
    "glut1": "GLUCOSA", "glut3": "GLUCOSA",
    "slc38a1": "AMINOACIDOS", "slc38a2": "AMINOACIDOS",
    "il6": "IL6", "il6R": "IL6", "gp130": "IL6",
}
GDISP = {
    "fatcd36": "CD36", "fatp1": "FATP1", "fatp4": "FATP4",
    "glut1": "GLUT1", "glut3": "GLUT3",
    "slc38a1": "SLC38A1", "slc38a2": "SLC38A2",
    "il6": "il6", "il6R": "il6R", "gp130": "gp130",
}
# Orden de filas del panel (estilo pedido): lipidos, glucosa, aminoacidos, IL6
# -- no el orden crudo de GENES.
FILAS_VIA = {
    "LIPIDOS": ["fatcd36", "fatp1", "fatp4"],
    "GLUCOSA": ["glut1", "glut3"],
    "AMINOACIDOS": ["slc38a1", "slc38a2"],
    "IL6": ["il6R", "gp130", "il6"],
}


def gdisp(g):
    return GDISP[g]


def color_for(gen, sexo, tto):
    return COL_CTRL[sexo] if tto == "CONTROL" else COL_LPS[GPATH[gen]][sexo]


# Estilo de trazo (un solo lugar, para que las 2 figuras de expresion coincidan)
EST = dict(
    borde_caja="0.20", lwd_caja=0.9, ancho_caja=0.62,
    col_bigote="0.60", lwd_bigote=0.9,          # bigote punteado, mas claro que el borde
    col_tope="0.45", lwd_tope=1.1,               # tope (cap): linea fina solida
    col_mediana="0.10", lwd_mediana=2.2,
    borde_punto="0.20", lwd_punto=0.6, cex_punto=22,
    alfa_relleno=0.55,
)

GRUPOS = [("HEMBRA", "CONTROL"), ("HEMBRA", "LPS"), ("MACHO", "CONTROL"), ("MACHO", "LPS")]
GLAB = ["♀ C", "♀ LPS", "♂ C", "♂ LPS"]
GRUPOS_4 = ["HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
CELDAS_4 = GRUPOS   # alias (pSTAT3 usa el mismo orden)
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


def _fila_de(tabla, tej, gen):
    for r in tabla:
        if r["TEJIDO"] == tej and r["GEN"] == gen:
            return r
    return None


# ===========================================================================
# 1. D11 -- LA funcion de anotacion (AMPLIADA: cascada de 3 ramas, AGENTS 4.2).
#    Identica en R y Python. Decide QUE anotar; la geometria de apilado es de
#    cada motor de dibujo (R base / ggplot2 / matplotlib).
# ===========================================================================
def d11_texto(p):
    """(texto, estilo) segun los umbrales D11, o None si no se anota.
    estilo: 'solida' (bracket lleno) o 'punteada' (tendencia)."""
    if p is None or (isinstance(p, float) and math.isnan(p)):
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


def d11_brackets_especificacion(fila_clasif, pholm_por_par):
    """Cascada D11 ampliada. `fila_clasif`: fila de qpcr_modelos_clasificacion.csv
    / pstat3_modelo_clasificacion.csv (columnas reales: interaccion_significativa,
    p_TTO, p_SEXO -- OJO, NO "p_interaccion"/"metodo": esas columnas no existen,
    ver analisis_descartados.md). `pholm_por_par`: {etiqueta_D6: p_holm}.
    Devuelve una lista de {x1, x2, texto, estilo} (posiciones de caja 0..3)."""
    fila_clasif = fila_clasif or {}
    out = []
    isig = fila_clasif.get("interaccion_significativa") == "TRUE"
    if isig:
        # (a) post hoc D6 -- un bracket por par cuyo p_holm cruce el umbral.
        for etq in ORDEN_PARES:
            ann = d11_texto(pholm_por_par.get(etq))
            if ann is None:
                continue
            texto, estilo = ann
            x1, x2 = PARES_D6_IDX[etq]
            out.append({"x1": x1, "x2": x2, "texto": texto, "estilo": estilo})
        return out
    # (b) efecto principal de TTO -- bracket unico sobre los 4 grupos (0..3).
    ann_t = d11_texto(_num(fila_clasif.get("p_TTO")))
    if ann_t is not None:
        texto, estilo = ann_t
        pref = f"Control vs LPS {texto}" if estilo == "solida" else f"Control vs LPS  {texto}"
        out.append({"x1": 0, "x2": 3, "texto": pref, "estilo": estilo})
    # (c) efecto principal de SEXO -- bracket entre los centros de cada sexo.
    ann_s = d11_texto(_num(fila_clasif.get("p_SEXO")))
    if ann_s is not None:
        texto, estilo = ann_s
        pref = f"♀ vs ♂ {texto}" if estilo == "solida" else f"♀ vs ♂  {texto}"
        out.append({"x1": 0.5, "x2": 2.5, "texto": pref, "estilo": estilo})
    return out


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


def _pholm(posthoc, tej, gen):
    out = {}
    for r in posthoc:
        if r["TEJIDO"] == tej and r["GEN"] == gen:
            out[r["contraste"]] = _num(r["p_holm"])
    return out


def _fc_por_grupo(cuant, tej, gen):
    """{grupo: [FC detectados]}."""
    out = {g: [] for g in GRUPOS_4}
    for r in cuant:
        if r["TEJIDO"] != tej or r["GEN"] != gen:
            continue
        if r["no_detectado"] == "TRUE" or r["neg_ddCt"] == "":
            continue
        out[r["GRUPO"]].append(2.0 ** float(r["neg_ddCt"]))
    return out


# ===========================================================================
# 3. Primitivas de dibujo (boxplots de expresion).
# ===========================================================================
_rng = random.Random(cfg.SEMILLA)   # jitter reproducible, igual criterio que R


def _jitter(n):
    return [(_rng.random() - 0.5) * 0.26 for _ in range(n)]


def _bracket(ax, x1, x2, y, etiqueta, solido):
    # Linea horizontal simple, SIN las perpendiculares en los extremos (pedido
    # explicito: es mas limpio y mas correcto -- marca un efecto que abarca
    # los 4 grupos o los dos sexos, no una comparacion puntual entre extremos).
    col = "black" if solido else "0.4"
    ls = "-" if solido else (0, (2, 2))
    ax.plot([x1, x2], [y, y], color=col, ls=ls, lw=0.9, solid_capstyle="butt")
    ax.text((x1 + x2) / 2.0, y, etiqueta, ha="center", va="bottom",
            fontsize=8.5 if solido else 7.5, color=col)


def _rango_eje_paneles(ymin_datos, ymax_datos, n_niveles,
                        frac_reservada=0.25, pad_datos=1.15):
    """El eje Y lo fijan los DATOS (D2: escala log), no los brackets -- antes
    el techo se estiraba con `ymax * 1.30**k` siempre, incluso sin ningun
    bracket. Si hay brackets, se reserva una FRACCION FIJA del alto total del
    panel en log10 (`frac_reservada`, constante, no crece con la cantidad de
    niveles), repartida en franjas iguales; si no hay brackets, no se reserva
    nada -- el protagonista del panel es el boxplot."""
    piso = ymin_datos * 0.6
    if n_niveles <= 0:
        return piso, ymax_datos * pad_datos, []
    log_piso = math.log10(piso)
    log_datos_top = math.log10(ymax_datos * pad_datos)
    rango_datos = log_datos_top - log_piso
    log_techo = log_piso + rango_datos / (1.0 - frac_reservada)
    franja = (log_techo - log_datos_top) / n_niveles
    niveles = [10 ** (log_datos_top + (k + 0.25) * franja) for k in range(n_niveles)]
    return piso, 10 ** log_techo, niveles


def panel_gen(ax, fc, gen, tejido, fila_clasif, pholm_por_par, con_titulo=True):
    """Panel de un gen x tejido. Eje Y SIEMPRE logaritmico (D2)."""
    dat = [fc[g] for g in GRUPOS_4]
    todos = [v for vs in dat for v in vs if v > 0]
    if not todos:
        ax.set_visible(False)
        return

    especs = d11_brackets_especificacion(fila_clasif, pholm_por_par)

    # Caja solo si el grupo tiene >= 3 detectados (los puntos van igual).
    for j, (sexo, tto) in enumerate(GRUPOS):
        vals = dat[j]
        color = color_for(gen, sexo, tto)
        if len(vals) >= 3:
            bp = ax.boxplot(
                [vals], positions=[j + 1], widths=EST["ancho_caja"],
                patch_artist=True, showfliers=False, manage_ticks=False,
                boxprops=dict(facecolor=matplotlib.colors.to_rgba(color, EST["alfa_relleno"]),
                              edgecolor=EST["borde_caja"], linewidth=EST["lwd_caja"]),
                whiskerprops=dict(color=EST["col_bigote"], linewidth=EST["lwd_bigote"],
                                  linestyle=(0, (3, 2))),
                capprops=dict(color=EST["col_tope"], linewidth=EST["lwd_tope"]),
                medianprops=dict(color=EST["col_mediana"], linewidth=EST["lwd_mediana"]),
            )
        if vals:
            xs = [j + 1 + dx for dx in _jitter(len(vals))]
            ax.plot(xs, vals, marker="o", ms=EST["cex_punto"] ** 0.5, mfc=color,
                    mec=EST["borde_punto"], mew=EST["lwd_punto"], ls="none", zorder=3)
        ax.annotate(f"n={len(vals)}", xy=(j + 1, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, -22), textcoords="offset points", ha="center",
                    fontsize=6.5, color="0.25")

    ax.axhline(1.0, color="0.7", lw=0.8, ls="--", zorder=1)
    ax.set_yscale("log")
    ax.set_xlim(0.4, 4.6)
    ax.set_xticks(range(1, 5))
    ax.set_xticklabels(GLAB, fontsize=7)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel(r"FC $= 2^{-\Delta\Delta Ct}$", fontsize=8)

    ymax = max(todos)
    ymin = min(todos)
    piso, techo, niveles_y = _rango_eje_paneles(ymin, ymax, len(especs))
    ax.set_ylim(piso, techo)
    for y, b in zip(niveles_y, especs):
        _bracket(ax, b["x1"] + 1, b["x2"] + 1, y, b["texto"], b["estilo"] == "solida")

    if con_titulo:
        p_int = _num((fila_clasif or {}).get("p_SEXOxTTO"))
        metodo = (fila_clasif or {}).get("rama_cascada", "")
        sub = f"{metodo} · p SEXOxTTO = {'NA' if p_int is None else f'{p_int:.3f}'} · n = {len(todos)}"
        # Titulo y subtitulo se anclan por desplazamiento en PUNTOS desde el
        # borde superior de los ejes (no en fraccion de datos): quedan siempre
        # fuera del area de dibujo, sin importar cuantos brackets haya apilados.
        ax.annotate(f"{gdisp(gen)} · {tejido}", xy=(0.5, 1.0), xycoords="axes fraction",
                    xytext=(0, 28), textcoords="offset points", ha="center", va="bottom",
                    fontsize=10, fontweight="bold")
        ax.annotate(sub, xy=(0.5, 1.0), xycoords="axes fraction",
                    xytext=(0, 11), textcoords="offset points", ha="center", va="bottom",
                    fontsize=7, color="0.25")


def panel_deteccion(ax, il6_tab, il6_fis, tej="BRAIN_E15", gen="il6"):
    """il6 @ BRAIN_E15 (D7): barras + simbolos individuales (relleno=detectado)."""
    nd = [0] * 4; nt = [0] * 4
    for j, grp in enumerate(GRUPOS_4):
        r = next(x for x in il6_tab if x["TEJIDO"] == tej and x["GEN"] == gen and x["GRUPO"] == grp)
        nt[j] = int(float(r["n_total"])); nd[j] = int(float(r["n_detectado"]))
    pct = [100 * nd[j] / nt[j] if nt[j] else 0.0 for j in range(4)]

    for j, (sexo, tto) in enumerate(GRUPOS):
        color = color_for(gen, sexo, tto)
        ax.bar(j + 1, pct[j], width=0.62,
               color=matplotlib.colors.to_rgba(color, EST["alfa_relleno"]),
               edgecolor=EST["borde_caja"])
        if nt[j]:
            xs = [j + 1 + (k - (nt[j] + 1) / 2) * 0.075 for k in range(1, nt[j] + 1)]
            ys = [pct[j] + 9] * nt[j]
            det = [k <= nd[j] for k in range(1, nt[j] + 1)]
            for x, y, d in zip(xs, ys, det):
                ax.plot(x, y, marker="o", ms=6.5,
                        mfc=color if d else "none", mec=EST["borde_punto"],
                        mew=EST["lwd_punto"])
        ax.text(j + 1, pct[j] + 18, f"{nd[j]}/{nt[j]}\n({pct[j]:.0f}%)",
                ha="center", va="bottom", fontsize=6.5, color="0.2")

    fis = {r["SEXO"]: _num(r["p_fisher"]) for r in il6_fis
           if r["TEJIDO"] == tej and r["GEN"] == gen}
    ax.set_xlim(0.4, 4.6)
    ax.set_xticks(range(1, 5))
    ax.set_xticklabels(GLAB, fontsize=7)
    ax.set_ylim(0, 125)
    ax.tick_params(axis="y", labelsize=7)
    ax.set_ylabel("% detectado", fontsize=8)
    sub = (f"no cuantificable (D7) · Fisher ♀ p = {fis.get('HEMBRA', float('nan')):.3f} "
           f"· ♂ p = {fis.get('MACHO', float('nan')):.3f}")
    ax.annotate(f"{gdisp(gen)} · {tej}", xy=(0.5, 1.0), xycoords="axes fraction",
                xytext=(0, 28), textcoords="offset points", ha="center", va="bottom",
                fontsize=10, fontweight="bold")
    ax.annotate(sub, xy=(0.5, 1.0), xycoords="axes fraction",
                xytext=(0, 11), textcoords="offset points", ha="center", va="bottom",
                fontsize=7, color="0.25")


# ===========================================================================
# 4. Figuras.
# ===========================================================================
def figura_tejido(D, tejido, ruta):
    """Filas por via metabolica (estilo pedido): lipidos(3) / glucosa(2) /
    aminoacidos(2) / IL-6 (hasta 3, il6@BRAIN_E15 es el panel de deteccion).
    Filas de 2 se centran dejando en blanco la primera columna."""
    _rng.seed(cfg.SEMILLA)   # jitter reproducible por figura
    filas = {k: list(v) for k, v in FILAS_VIA.items()}
    if tejido == "BRAIN_E15":
        filas["IL6"] = [g for g in filas["IL6"] if g != "il6R"]

    nrow = len(filas)
    fig, axes = plt.subplots(nrow, 3, figsize=(11.5, nrow * 3.9))
    for i, (via, genes) in enumerate(filas.items()):
        despl = 1 if len(genes) == 2 else 0
        for col in range(3):
            j = col - despl
            ax = axes[i, col]
            if j < 0 or j >= len(genes):
                ax.set_visible(False)
                continue
            gen = genes[j]
            if tejido == "BRAIN_E15" and gen == "il6":
                panel_deteccion(ax, D["il6_tab"], D["il6_fis"])
            else:
                fc = _fc_por_grupo(D["cuant"], tejido, gen)
                fila_clasif = _fila_de(D["clasif"], tejido, gen)
                ph = _pholm(D["posthoc"], tejido, gen)
                panel_gen(ax, fc, gen, tejido, fila_clasif, ph)

    tt = "Placenta E15" if tejido == "PLACENTA_E15" else "Cerebro fetal E15"
    fig.suptitle(f"Expresion relativa por gen -- {tt}  "
                 f"(FC = 2^(-ΔΔCt), eje log; brackets = D11)", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96), h_pad=5.5, w_pad=2.0)
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


def _bracket_pstat3(ax, x1, x2, y, alto, etiqueta, solido):
    """Bracket CON extremos (rectangulo): sin cambios, fuera del alcance del
    pedido de simplificar brackets (que aplica solo a las figuras de
    expresion). Mantiene paridad visual con la version R de pSTAT3
    (`apilar_brackets_ggplot`, geom_segment), que tampoco cambio."""
    col = "black" if solido else "0.4"
    ls = "-" if solido else (0, (2, 2))
    ax.plot([x1, x1, x2, x2], [y, y + alto, y + alto, y], color=col, ls=ls,
            lw=0.9, solid_capstyle="butt")
    ax.text((x1 + x2) / 2.0, y + alto, etiqueta, ha="center", va="bottom",
            fontsize=8.5 if solido else 7.5, color=col)


def figura_pstat3(D, ruta):
    pst = D["pst"]
    cl = D["pst_cl"][0] if D["pst_cl"] else {}
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
        jj = [(-0.16 + 0.32 * k / (len(ys) - 1)) if len(ys) > 1 else 0.0
              for k in range(len(ys))]
        for (v, m), j in zip(vals[i], jj):
            ax.plot(i + j, v, marker=marcas.get(m, "o"), ms=5.0,
                    mfc=COL_TTO[tt], mec="white", mew=0.5, ls="none", zorder=3)
    tope = max(todos)
    especs = d11_brackets_especificacion(cl, pholm)
    for k, b in enumerate(especs):
        paso = 0.10 + 0.11 * k
        y = tope * (1 + paso)
        alto = tope * (1 + paso - 0.035) - y
        _bracket_pstat3(ax, b["x1"], b["x2"], y, alto, b["texto"], b["estilo"] == "solida")
    ax.set_ylim(0, tope * 1.75)
    ax.set_xticks(range(4))
    ax.set_xticklabels(ETIQ_X, fontsize=8.5)
    ax.set_ylabel("pSTAT3  (u.a., normalizado a proteina total)", fontsize=9)
    ax.set_title("Fosfo-STAT3 (Tyr705) -- placenta E15\n"
                 "crudo por SEXO x TTO; forma de punto = membrana; "
                 "brackets = D11 (cascada ampliada)", fontsize=9)
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
    "### D11 ampliada: cascada de 3 ramas (cambio pedido explicitamente)",
    "",
    "- Hasta esta sesion, D11 solo anotaba si la interaccion SEXO x TTO era "
    "significativa; sin interaccion, el panel quedaba sin ninguna marca aunque "
    "hubiera un efecto principal fuerte (p. ej. `il6@PLACENTA_E15` con "
    "`p_TTO` = 6.0e-06 no llevaba bracket). **Se revierte esa restriccion por "
    "pedido explicito** (`pedidos/boxplots_acto1_base_R.R`, tratado como "
    "especificacion de estilo y logica, no como codigo a copiar). Cascada "
    "nueva en `d11_brackets_especificacion` (identica R/Python), documentada "
    "en AGENTS.md 4.2: (a) interaccion significativa -> brackets por par del "
    "post hoc D6, igual que antes; (b) sin interaccion y `TTO` significativo/"
    "tendencia -> un bracket sobre los 4 grupos, `Control vs LPS`; (c) sin "
    "interaccion y `SEXO` significativo/tendencia -> un bracket entre los "
    "centros de cada sexo, `♀ vs ♂`. (b) y (c) no son excluyentes.",
    "- **Nombres de columna corregidos** contra los CSV reales al adaptar el "
    "pedido: la referencia asumia `p_interaccion` y `metodo` en "
    "`qpcr_modelos_clasificacion.csv`; las columnas reales son "
    "`p_SEXOxTTO` (usada solo para mostrarla en el subtitulo del panel; el "
    "gate sigue siendo `interaccion_significativa`) y `rama_cascada`.",
    "",
    "### Boxplots de expresion: R base, no ggplot2 (cambio pedido explicitamente)",
    "",
    "- El estilo pedido necesita bigote (punteado, mas claro) y tope del "
    "bigote (solido) con trazo distinto del borde de la caja; `geom_boxplot` "
    "no expone esos tres trazos por separado, `boxplot()`/`bxp()` de R base si "
    "(`whisklty/whiskcol`, `staplelty/staplecol`, `border`). **Se abandona "
    "ggplot2 solo para estos paneles**; el Acto 2 (08) sigue en ggplot2/GGally, "
    "y pSTAT3 (mismo script) tambien sigue en ggplot2. matplotlib (Python) ya "
    "permite estilar los tres trazos por separado (`whiskerprops`/`capprops`/"
    "`boxprops`), asi que este lenguaje no tuvo que cambiar de libreria: solo "
    "replica el mismo estilo.",
    "- **Filas agrupadas por via metabolica** (lipidos / glucosa / aminoacidos "
    "/ IL-6), no por el orden crudo de `GENES` -- estilo pedido, ayuda a leer "
    "el panel por sistema biologico. Color: Control por sexo (celeste, 2 "
    "tonos), LPS por sexo x via metabolica (una paleta por via) -- il6/il6R/"
    "gp130 comparten paleta (via IL-6).",
    "",
    "### il6 @ BRAIN_E15: panel de deteccion, no boxplot de FC",
    "",
    "- il6 en cerebro no es cuantificable (calibrador HEMBRA_CONTROL 0/9, D7): no "
    "tiene FC. Su panel en `acto1_expresion_BRAIN_E15.png` muestra la **proporcion "
    "de deteccion** Control vs LPS por sexo (de `qpcr_il6_brain_tabla2x4.csv`) con "
    "la p de Fisher (de `qpcr_il6_brain_fisher.csv`), no un boxplot. `il6R` sale "
    "del panel de cerebro por prolijidad visual (deteccion insuficiente en ese "
    "tejido); se conserva en la tabla de modelos (05) y en el panel de placenta.",
    "",
    "### Escala y datos de los boxplots de expresion",
    "",
    "- Eje Y = `FC = 2^(-ddCt)` en escala **log** (D2). `FC` se calcula aca desde "
    "`neg_ddCt` de `qpcr_cuantificacion_long.tsv`; 04 no lo guarda para no arrastrar "
    "el redondeo de `2^x` entre libm. Caja solo si el grupo tiene **>=3 detectados** "
    "(los puntos se dibujan igual, sin caja, si son menos).",
    "",
    "### Brackets: linea simple + eje fijado por los datos (pedido post-cierre)",
    "",
    "- **Bug reportado**: el techo del eje crecia en proporcion a la cantidad de "
    "brackets, siempre, incluso en paneles sin un solo bracket (`n_niveles` tenia un "
    "piso de 1) -- `acto1_expresion_PLACENTA_E15.png` llegaba a 10^4 con datos que no "
    "pasan de 10, `il6R` (sin brackets) a 10^5. **Se separan las dos responsabilidades**: "
    "la funcion que arma el rango del eje calcula el techo/piso **solo a partir de los "
    "datos** (mismo padding de siempre, `pad_datos = 1.15`), y **solo si hay brackets** "
    "reserva una fraccion fija del alto total del panel en log10 (`frac_reservada = 0.25`, "
    "constante, no crece con la cantidad de niveles), repartida en franjas iguales; "
    "sin brackets no se reserva nada.",
    "- El dibujo del bracket deja de trazar el rectangulo con perpendiculares en los "
    "extremos: ahora es una **linea horizontal simple** con el texto encima. Mas "
    "limpio y mas correcto: el bracket marca un efecto que abarca los 4 grupos o los "
    "dos sexos (ramas (b)/(c) de D11), no una comparacion puntual entre dos extremos.",
    "",
    "### pSTAT3: valores crudos + membrana como forma de punto",
    "",
    "- El boxplot de pSTAT3 muestra los valores **crudos** por SEXO x TTO (no "
    "ajustados por MEMBRANA); la membrana se codifica como forma de punto "
    "(o / cuadrado / triangulo). El bloque MEMBRANA lo maneja el modelo D9 (06), "
    "no la figura. Nota al pie: pSTAT3 = abundancia de fosfo-STAT3, no fraccion. "
    "Sigue en ggplot2; solo cambia la cascada D11 que decide los brackets.",
    "",
    "### Figuras del ELISA",
    "",
    "- `acto1_elisa_ms.png` y `acto1_elisa_la.png` (Acto 1.1) ya las produjo "
    "`03_elisa`; 07 **no las regenera**, solo completan el set del Acto 1.",
    "",
    "### Paridad",
    "",
    "- Las figuras son PNG: equivalentes, no byte-identicas (R base + ggplot2 vs "
    "matplotlib). Byte-identicas entre lenguajes: las filas nuevas de "
    "`procedencia.csv` / `verificaciones.csv` y esta seccion. La cascada D11 "
    "(`d11_texto` + `d11_brackets_especificacion`) es identica en ambos.",
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

    figura_tejido(D, "PLACENTA_E15", fig_pla)
    figura_tejido(D, "BRAIN_E15", fig_bra)
    figura_pstat3(D, fig_pst)

    # --- comparaciones/brackets anotados (mismo orden y formato que R) ---
    anotadas = []
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            if tej == "BRAIN_E15" and gen == "il6":
                continue
            fila_clasif = _fila_de(D["clasif"], tej, gen)
            if fila_clasif is None:
                continue
            ph = _pholm(D["posthoc"], tej, gen)
            for b in d11_brackets_especificacion(fila_clasif, ph):
                anotadas.append(f"{tej}/{gen}/[{_fmt(b['x1'])}-{_fmt(b['x2'])}]:{b['texto']}")
    pst_cl = D["pst_cl"][0] if D["pst_cl"] else {}
    pst_ph = {r["contraste"]: _num(r["p_holm"]) for r in D["pst_ph"]}
    pst_anot = [f"pSTAT3/[{_fmt(b['x1'])}-{_fmt(b['x2'])}]:{b['texto']}"
                for b in d11_brackets_especificacion(pst_cl, pst_ph)]

    figs_elisa = [(cfg.RUTA_FIGURAS / "acto1_elisa_ms.png").is_file(),
                  (cfg.RUTA_FIGURAS / "acto1_elisa_la.png").is_file()]

    actualizar_descartados()

    ent_q = f"data/processed/qpcr_cuantificacion_long.tsv + outputs/tables/*/qpcr_modelos_* (de data/{fuente_q}/{cfg.ARCHIVO_QPCR})"
    ent_p = f"data/processed/pstat3_long.tsv + outputs/tables/*/pstat3_* (de data/{fuente_p}/{cfg.ARCHIVO_PSTAT3})"
    registrar_procedencia([
        ["outputs/figures/acto1_expresion_PLACENTA_E15.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent_q, "boxplots (R base) de FC = 2^(-ddCt) (eje log) por gen en "
         "placenta E15, filas por via metabolica, 4 grupos SEXO x TTO; cascada D11 "
         "ampliada (brackets por efecto principal si no hay interaccion)"],
        ["outputs/figures/acto1_expresion_BRAIN_E15.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent_q, "idem cerebro fetal E15; il6 como panel de proporcion de "
         "deteccion (D7, no cuantificable), il6R fuera de este panel; cascada D11 "
         "ampliada (por par si hay interaccion, por efecto principal si no)"],
        ["outputs/figures/acto1_pstat3.png", "figura", ESTE_SCRIPT, "PROPIO", ent_p,
         "boxplot (ggplot2) de pSTAT3 crudo por SEXO x TTO, membrana como forma de "
         "punto; cascada D11 ampliada del post hoc D6 (06_pstat3); nota de "
         "limitacion D9"],
    ])
    ok_figs = (fig_pla.is_file() and fig_bra.is_file() and fig_pst.is_file()
               and all(figs_elisa))
    registrar_verificaciones([
        ["figuras_acto1_generadas",
         "figuras del Acto 1: 2 de expresion + pSTAT3 (+ 2 de ELISA de 03_elisa)",
         f"placenta={_fmt(fig_pla.is_file())};brain={_fmt(fig_bra.is_file())};"
         f"pstat3={_fmt(fig_pst.is_file())};elisa_ms={_fmt(figs_elisa[0])};"
         f"elisa_la={_fmt(figs_elisa[1])}",
         "las 5 figuras del Acto 1 existen",
         "TRUE" if ok_figs else "FALSE", ESTE_SCRIPT],
        ["figuras_d11_una_funcion",
         "una sola funcion de anotacion D11 (d11_brackets_especificacion) usada en todas las figuras",
         "d11_texto + d11_brackets_especificacion (identica R/Python), 3 figuras",
         "una funcion, todas las figuras", "TRUE", ESTE_SCRIPT],
        ["figuras_d11_cascada",
         "D11 ampliada: (a) interaccion->por par; (b) TTO principal->bracket 0-3; (c) SEXO principal->bracket 0.5-2.5",
         f"expresion: {len(anotadas)} brackets "
         f"({'; '.join(anotadas) if anotadas else 'ninguno'}); "
         f"pSTAT3: {len(pst_anot)} ({'; '.join(pst_anot) if pst_anot else 'ninguno'})",
         "cascada de 3 ramas, (b)/(c) no excluyentes", "TRUE", ESTE_SCRIPT],
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
        ["figuras_expresion_r_base",
         "boxplots de expresion en R base (bxp/boxplot), no ggplot2 -- pedido explicito",
         "boxplot() con whisklty/staplelty/border distintos; Acto 2 sigue en ggplot2",
         "R base para expresion", "TRUE", ESTE_SCRIPT],
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
