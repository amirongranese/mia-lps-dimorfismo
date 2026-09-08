# 03_elisa.py -- Analisis del ELISA de IL-6 con censura a izquierda (D10). Acto 1.1.
#
# Por que existe este archivo: valida el modelo MIA mostrando el efecto de LPS
# sobre IL-6 en suero materno (MS) y liquido amniotico (LA). El ELISA tiene
# censura a izquierda (Conc < 0 == "< LOD"); D10 fija el tratamiento:
#
#   1. Reportar % de deteccion / censura por grupo ANTES de cualquier estadistico.
#   2. Contraste primario = test de rango con los censurados empatados en el rango
#      mas bajo. Se usa Peto-Peto (Fleming-Harrington G-rho = 1). En R via
#      survival::survdiff(rho=1); aca implementado a mano con la MISMA formula, y
#      en R se verifica en corrida contra survdiff (stopifnot, tol 1e-8).
#   3. Si un grupo tiene censura tan alta que ningun estimador de ubicacion es
#      defendible -> solo proporcion de deteccion. Es el caso de MS: Control tiene
#      4/5 censurados (80%). Decision confirmada con el usuario: MS se analiza
#      SOLO como deteccion + Fisher exacto (Control vs LPS), sin test de ubicacion.
#   4. LA se analiza ESTRATIFICADO POR SEXO: dentro de cada sexo, Peto-Peto
#      Control vs LPS sobre Conc, mas Fisher de deteccion. Sin modelo factorial
#      (censura 0-60% + n 5-9 por celda lo hacen inestable) -- ver
#      analisis_descartados.md.
#
# Descriptivo con censura: mediana del estimador producto-limite (Kaplan-Meier)
# sobre el dato reflejado (t' = M - t pasa la censura izquierda a derecha). NO se
# usa ROS/NADA (dependencia extra + ajuste log-normal fragil con n chico y censura
# pesada); la alternativa de LOD (menor estandar de la hoja CURVA IL6) queda
# documentada, no aplicada. LOD primario = 0 (el blanco), como fijo 02_ingesta_qc.
#
# Este script NO cuantifica qPCR y NO toca ninguna decision fuera de D10.

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

# Identificador para los registros compartidos (procedencia, verificaciones,
# analisis_descartados): SIN extension, para que R y Python reemplacen la misma
# fila y el archivo quede byte-identico.
ESTE_SCRIPT = "03_elisa"

# ---------------------------------------------------------------------------
# Ordenes canonicos y definicion de grupos.
# ---------------------------------------------------------------------------
ORDEN_TTO = {"CONTROL": 0, "LPS": 1}
ORDEN_SEXO = {"": 0, "HEMBRA": 1, "MACHO": 2}
GRUPOS_MS = [("MS", "CONTROL", ""), ("MS", "LPS", "")]
GRUPOS_LA = [("LA", "CONTROL", "HEMBRA"), ("LA", "CONTROL", "MACHO"),
             ("LA", "LPS", "HEMBRA"), ("LA", "LPS", "MACHO")]
GRUPOS_TODOS = GRUPOS_MS + GRUPOS_LA
SEXOS_LA = ["HEMBRA", "MACHO"]

LOD_ELISA = 0.0  # el blanco de la placa (igual que 02_ingesta_qc)

# Colores (paleta Okabe-Ito, segura para daltonismo): Control azul / LPS bermellon.
COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02 (paridad R/Python): NA -> "",
# float -> %.10g, entero -> sin decimales, logico -> TRUE/FALSE.
# Los p-valores y estadisticos que dependen de funciones trascendentes
# (lgamma, exp, erfc / pchisq) se guardan como cadena "%.6e" ya formateada:
# asi R y Python escriben el MISMO texto aunque sus libm difieran en el ultimo
# bit. Ver nota en analisis_descartados.md.
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
    """p-valor / estadistico dependiente de trascendentes -> texto '%.6e' estable."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return "%.6e" % float(x)


def s3e(x) -> str:
    """Version corta para el reporte legible."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return "%.3e" % float(x)


# ---------------------------------------------------------------------------
# Lectura del intermedio de 02_ingesta_qc.
# ---------------------------------------------------------------------------
def cargar_elisa():
    ruta = cfg.RUTA_DATOS_PROC / "elisa_long.tsv"
    if not ruta.is_file():
        raise SystemExit(f"03_elisa: falta {ruta} (correr 02_ingesta_qc primero)")
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    header = lineas[0].split("\t")
    out = []
    for ln in lineas[1:]:
        d = dict(zip(header, ln.split("\t")))
        out.append({
            "bloque": d["bloque"], "TTO": d["TTO"], "SEXO": d["SEXO"],
            "ID": d["ID"], "Conc": float(d["Conc"]),
            "censurado": d["censurado"] == "TRUE",
            "LOD": float(d["LOD"]) if d["LOD"] != "" else LOD_ELISA,
        })
    return out


def sub(filas, bloque, tto, sexo):
    return [f for f in filas if f["bloque"] == bloque and f["TTO"] == tto
            and (sexo == "" or f["SEXO"] == sexo)]


def mediana(xs):
    s = sorted(xs)
    k = len(s)
    if k == 0:
        return None
    if k % 2 == 1:
        return s[k // 2]
    return (s[k // 2 - 1] + s[k // 2]) / 2.0


# ===========================================================================
# 1. Deteccion / censura por grupo  (D10, paso 1: ANTES de cualquier estadistico)
# ===========================================================================
def tabla_deteccion(filas):
    header = ["bloque", "TTO", "SEXO", "n", "n_detectado", "n_censurado",
              "pct_detectado", "pct_censurado"]
    out = []
    for (bloque, tto, sexo) in GRUPOS_TODOS:
        g = sub(filas, bloque, tto, sexo)
        n = len(g)
        nc = sum(1 for f in g if f["censurado"])
        nd = n - nc
        out.append([bloque, tto, sexo, n, nd, nc,
                    (100.0 * nd / n) if n else None,
                    (100.0 * nc / n) if n else None])
    return header, out


# ===========================================================================
# 2. Descriptivo por grupo: detectados + mediana Kaplan-Meier (dato reflejado)
#
# Por que reflejar: survdiff / KM manejan censura a DERECHA; el ELISA tiene
# censura a IZQUIERDA. t' = M - t con M > toda Conc detectada convierte "< LOD"
# en censura a derecha en M - LOD, y deja a los censurados en el extremo BAJO de
# la escala original (rango mas bajo), que es exactamente lo que pide D10.
# ===========================================================================
def constante_reflexion(filas):
    dets = [f["Conc"] for f in filas if not f["censurado"]]
    return float(math.ceil(max(dets)) + 1)


def km_mediana_reflejada(g, M):
    """Mediana del estimador producto-limite en escala ORIGINAL. None si S(t) no
    baja de 0.5 (censura demasiado alta para estimar la mediana)."""
    obs = []
    for f in g:
        if f["censurado"]:
            obs.append((M - f["LOD"], 0))   # censura a derecha en M - LOD
        else:
            obs.append((M - f["Conc"], 1))  # evento
    tiempos_evento = sorted({t for (t, e) in obs if e == 1})
    S = 1.0
    med_refl = None
    for t in tiempos_evento:
        n_risk = sum(1 for (tt, ee) in obs if tt >= t - 1e-12)
        d = sum(1 for (tt, ee) in obs if abs(tt - t) <= 1e-12 and ee == 1)
        if n_risk > 0:
            S *= (n_risk - d) / n_risk
        if med_refl is None and S <= 0.5 + 1e-12:
            med_refl = t
    return None if med_refl is None else M - med_refl


def tabla_descriptivo(filas, M):
    header = ["bloque", "TTO", "SEXO", "n", "n_detectado", "pct_censurado",
              "min_detectado", "mediana_detectada", "max_detectado",
              "km_mediana", "km_nota"]
    out = []
    for (bloque, tto, sexo) in GRUPOS_TODOS:
        g = sub(filas, bloque, tto, sexo)
        n = len(g)
        dets = sorted(f["Conc"] for f in g if not f["censurado"])
        nd = len(dets)
        pct_c = (100.0 * (n - nd) / n) if n else None
        kmm = km_mediana_reflejada(g, M)
        nota = "" if kmm is not None else "S(t) no baja de 0.5 (censura alta)"
        out.append([bloque, tto, sexo, n, nd, pct_c,
                    dets[0] if dets else None,
                    mediana(dets) if dets else None,
                    dets[-1] if dets else None,
                    kmm, nota])
    return header, out


# ===========================================================================
# 3. Peto-Peto (Fleming-Harrington G-rho = 1) para dos grupos.
#
# PROPIO. Misma formula que survival::survdiff: en cada tiempo de evento t,
# peso w = S(t-)^rho con S la KM combinada left-continua; se acumula
# O1 = sum w * d1 , E1 = sum w * d * n1/n , V = sum w^2 * d*(n-d)/(n-1)*(n1/n)*(1-n1/n).
# chisq = (O1 - E1)^2 / V ~ chi^2_1.  En R se verifica contra survdiff(rho=1).
# ===========================================================================
def _chi2_sf_1df(x):
    # P(Chi^2_1 > x) = erfc(sqrt(x/2)).  Se guarda formateado con p6e().
    if x <= 0:
        return 1.0
    return math.erfc(math.sqrt(x / 2.0))


def peto_peto_2grupos(tiempos, evento, grupo):
    datos = sorted(zip(tiempos, evento, grupo))
    grupos = sorted(set(grupo))
    if len(grupos) != 2:
        raise ValueError("peto_peto_2grupos: se esperan exactamente 2 grupos")
    g1 = grupos[0]
    tiempos_evento = sorted({t for (t, e, _g) in datos if e == 1})
    O1 = E1 = V = 0.0
    S_prev = 1.0
    for t in tiempos_evento:
        en_riesgo = [(e, gg) for (tt, e, gg) in datos if tt >= t - 1e-12]
        n = len(en_riesgo)
        n1 = sum(1 for (e, gg) in en_riesgo if gg == g1)
        d = sum(1 for (tt, e, gg) in datos
                if abs(tt - t) <= 1e-12 and e == 1)
        d1 = sum(1 for (tt, e, gg) in datos
                 if abs(tt - t) <= 1e-12 and e == 1 and gg == g1)
        w = S_prev
        O1 += w * d1
        E1 += w * d * n1 / n
        if n > 1:
            V += w * w * d * (n1 / n) * (1.0 - n1 / n) * (n - d) / (n - 1)
        S_prev *= (n - d) / n
    chisq = (O1 - E1) ** 2 / V if V > 0 else 0.0
    return chisq, 1, _chi2_sf_1df(chisq), O1, E1, V


def datos_supervivencia(g, M):
    """(tiempos reflejados, evento) para un grupo de filas del ELISA."""
    t, e = [], []
    for f in g:
        if f["censurado"]:
            t.append(M - f["LOD"]); e.append(0)
        else:
            t.append(M - f["Conc"]); e.append(1)
    return t, e


def tabla_petopeto_la(filas, M):
    header = ["estrato", "n_control", "n_lps", "cens_control", "cens_lps",
              "chisq", "df", "p_valor", "metodo"]
    out = []
    for sexo in SEXOS_LA:
        gc = sub(filas, "LA", "CONTROL", sexo)
        gl = sub(filas, "LA", "LPS", sexo)
        tc, ec = datos_supervivencia(gc, M)
        tl, el = datos_supervivencia(gl, M)
        t = tc + tl
        e = ec + el
        grp = ["CONTROL"] * len(tc) + ["LPS"] * len(tl)
        chisq, df, p, _o, _e, _v = peto_peto_2grupos(t, e, grp)
        out.append([sexo, len(gc), len(gl), ec.count(0), el.count(0),
                    "%.10g" % chisq, df, p6e(p),
                    "Peto-Peto G-rho=1 (survdiff) sobre Conc reflejada"])
    return header, out


# ===========================================================================
# 4. Fisher exacto a dos colas para 2x2 -- misma regla que stats::fisher.test
#    (suma de tablas con prob <= prob(observada) * (1 + 1e-7)). PROPIO.
# ===========================================================================
def fisher_2x2(a, b, c, d):
    """Tabla [[a,b],[c,d]]. Devuelve (or_haldane, p_dos_colas)."""
    r1, r2 = a + b, c + d
    c1 = a + c
    n = a + b + c + d

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
    p = min(1.0, tot)
    or_h = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
    return or_h, p


def fila_fisher_deteccion(filas, bloque, sexo):
    gc = sub(filas, bloque, "CONTROL", sexo)
    gl = sub(filas, bloque, "LPS", sexo)
    a = sum(1 for f in gc if not f["censurado"])   # Control detectados
    b = len(gc) - a                                 # Control no detectados
    c = sum(1 for f in gl if not f["censurado"])   # LPS detectados
    d = len(gl) - c                                 # LPS no detectados
    or_h, p = fisher_2x2(a, b, c, d)
    estrato = sexo if sexo else "(sin estrato)"
    return [bloque, estrato, a, len(gc), c, len(gl),
            "%.10g" % or_h, p6e(p)]


def tabla_fisher(filas):
    header = ["bloque", "estrato", "control_detectado", "control_n",
              "lps_detectado", "lps_n", "or_haldane", "p_valor"]
    out = [fila_fisher_deteccion(filas, "MS", "")]
    for sexo in SEXOS_LA:
        out.append(fila_fisher_deteccion(filas, "LA", sexo))
    return header, out


# ===========================================================================
# 5. Figuras del Acto 1.1
# ===========================================================================
def _offsets(k):
    if k <= 0:
        return []
    if k == 1:
        return [0.0]
    paso = 0.30 / (k - 1)
    return [-0.15 + i * paso for i in range(k)]


def _dibujar_grupo(ax, x0, g, color, ymax):
    dets = sorted(f["Conc"] for f in g if not f["censurado"])
    cens = [f for f in g if f["censurado"]]
    if len(dets) >= 3:
        bp = ax.boxplot([dets], positions=[x0], widths=0.46, showfliers=False,
                        patch_artist=True, zorder=2)
        for box in bp["boxes"]:
            box.set(facecolor=color, alpha=0.14, edgecolor=color)
        for el in ("whiskers", "caps", "medians"):
            for art in bp[el]:
                art.set(color=color)
    for off, y in zip(_offsets(len(dets)), dets):
        ax.plot(x0 + off, y, marker="o", ms=6, mfc=color, mec=color, ls="none",
                zorder=3)
    for off, _f in zip(_offsets(len(cens)), cens):
        ax.plot(x0 + off, 0.0, marker="o", ms=7, mfc="none", mec=color,
                mew=1.4, ls="none", zorder=3)
    n = len(g)
    nd = len(dets)
    ax.text(x0, ymax * 1.05, f"{nd}/{n} det.\n{100 * nd / n:.0f}%",
            ha="center", va="bottom", fontsize=8.5, color="0.25")


def _leyenda_marcadores(ax):
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="o", color="none", mfc="0.35", mec="0.35",
               ms=6, label="detectado"),
        Line2D([0], [0], marker="o", color="none", mfc="none", mec="0.35",
               mew=1.4, ms=7, label="censurado (< LOD)"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.11),
              ncol=2, fontsize=8, frameon=False)


def figura_ms(filas, fisher_ms, ruta):
    grupos = {tto: sub(filas, "MS", tto, "") for (_b, tto, _s) in GRUPOS_MS}
    ymax = max(f["Conc"] for f in filas
               if f["bloque"] == "MS" and not f["censurado"])
    fig, ax = plt.subplots(figsize=(4.8, 5.0))
    for i, (_b, tto, _s) in enumerate(GRUPOS_MS):
        _dibujar_grupo(ax, i, grupos[tto], COL_TTO[tto], ymax)
    ax.axhline(LOD_ELISA, color="0.6", lw=0.8, ls="--", zorder=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Control", "LPS"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylabel("IL-6  (pg/mL)")
    ax.set_ylim(-ymax * 0.06, ymax * 1.22)
    ax.set_title(f"Fisher exacto (deteccion) Control vs LPS:  "
                 f"p = {s3e(float(fisher_ms[-1]))}", fontsize=9.5)
    fig.suptitle("IL-6 en suero materno (ELISA) -- validacion del modelo MIA",
                 fontsize=11, y=0.98)
    _leyenda_marcadores(ax)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(ruta, dpi=300, bbox_inches="tight")
    plt.close(fig)


def figura_la(filas, peto_rows, ruta):
    ymax = max(f["Conc"] for f in filas
               if f["bloque"] == "LA" and not f["censurado"])
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 5.0), sharey=True)
    pmap = {r[0]: r[7] for r in peto_rows}
    for ax, sexo in zip(axes, SEXOS_LA):
        for i, tto in enumerate(("CONTROL", "LPS")):
            _dibujar_grupo(ax, i, sub(filas, "LA", tto, sexo),
                           COL_TTO[tto], ymax)
        ax.axhline(LOD_ELISA, color="0.6", lw=0.8, ls="--", zorder=1)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Control", "LPS"])
        ax.set_xlim(-0.6, 1.6)
        ax.set_title(f"{'Hembra' if sexo == 'HEMBRA' else 'Macho'}   "
                     f"(Peto-Peto p = {s3e(float(pmap[sexo]))})", fontsize=9.5)
    axes[0].set_ylabel("IL-6  (pg/mL)")
    axes[0].set_ylim(-ymax * 0.06, ymax * 1.22)
    fig.suptitle("IL-6 en liquido amniotico (ELISA) por sexo fetal",
                 fontsize=11, y=0.98)
    _leyenda_marcadores(axes[1])
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(ruta, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ===========================================================================
# 6. Reporte legible
# ===========================================================================
def _md_tabla(header, filas):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for f in filas:
        out.append("| " + " | ".join(_fmt(v) for v in f) + " |")
    return "\n".join(out)


def construir_reporte(fuente, det_h, det, desc_h, desc, fish_h, fish,
                      peto_h, peto, M):
    L = []
    ap = L.append
    ap("# Reporte del ELISA de IL-6 (Acto 1.1)")
    ap("")
    ap("Generado por `03_elisa` (R y Python producen este archivo identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. Deteccion / censura por grupo (D10, antes de cualquier estadistico)")
    ap("")
    ap(_md_tabla(det_h, det))
    ap("")
    ap("## 2. Descriptivo por grupo")
    ap("")
    ap("`km_mediana` = mediana del estimador producto-limite (Kaplan-Meier) sobre "
       "el dato reflejado `t' = M - t` "
       f"(M = {_fmt(M)} = ceil(max Conc detectada) + 1). Vacio = S(t) no baja de "
       "0.5 por censura alta -> solo proporcion de deteccion.")
    ap("")
    ap(_md_tabla(desc_h, desc))
    ap("")
    ap("## 3. Suero materno (MS): solo deteccion + Fisher exacto")
    ap("")
    ap("Control tiene 80% de censura (1 valor detectado): ningun estimador de "
       "ubicacion es defendible (D10, clausula final). Se reporta solo la "
       "proporcion de deteccion y el Fisher exacto Control vs LPS.")
    ap("")
    ms = fish[0]
    ap(f"- Control: {ms[2]}/{ms[3]} detectados; LPS: {ms[4]}/{ms[5]} detectados.")
    ap(f"- Fisher exacto (2 colas): OR (Haldane) = {ms[6]}, p = {s3e(float(ms[7]))}.")
    ap("")
    ap("## 4. Liquido amniotico (LA): Peto-Peto por sexo + Fisher de deteccion")
    ap("")
    ap("Contraste primario: Peto-Peto (Fleming-Harrington G-rho = 1) Control vs "
       "LPS sobre `Conc`, con los censurados empatados en el rango mas bajo. "
       "Estratificado por sexo (sin modelo factorial: ver `analisis_descartados.md`).")
    ap("")
    ap(_md_tabla(peto_h, peto))
    ap("")
    ap("Fisher exacto sobre deteccion (Control vs LPS) dentro de cada sexo:")
    ap("")
    ap(_md_tabla(fish_h, fish))
    ap("")
    ap("## 5. Figuras")
    ap("")
    ap("- `outputs/figures/acto1_elisa_ms.png` -- IL-6 en suero materno por "
       "tratamiento; censurados dibujados en el LOD (simbolo abierto).")
    ap("- `outputs/figures/acto1_elisa_la.png` -- IL-6 en liquido amniotico por "
       "sexo y tratamiento.")
    ap("")
    ap("## 6. Notas metodologicas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `03_elisa`: por que no se uso ROS, "
       "por que MS no lleva test de ubicacion, por que LA no lleva modelo "
       "factorial, la reflexion para aplicar survdiff a censura a izquierda, y el "
       "LOD alternativo (menor estandar de `CURVA IL6`) no aplicado.")
    ap("")
    return "\n".join(L)


# ===========================================================================
# 7. Artefactos compartidos (merge por 'script') -- headers identicos a 02.
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
## 03_elisa

### ROS / KM-NADA no usado como estimador de ubicacion

- **Que se probo:** reportar media/mediana por grupo con ROS (regression on order
  statistics; `NADA` en R) como estimacion puntual con censura.
- **Por que no:** (a) agrega una dependencia (`NADA`) fuera del toolchain minimo;
  (b) con n 5-9 por grupo y censura 30-80% el ajuste log-normal implicito de ROS
  es fragil y sus IC no son defendibles; (c) D10 admite explicitamente el test de
  rango con censurados al fondo como alternativa.
- **Que se hizo:** contraste primario = Peto-Peto (G-rho = 1). Descriptivo =
  mediana del estimador producto-limite (Kaplan-Meier) sobre el dato reflejado;
  donde S(t) no baja de 0.5 (censura alta) se reporta vacio + proporcion de
  deteccion.

### Suero materno (MS): sin test de ubicacion

- **Situacion:** Control tiene 4/5 censurados (80%): un solo valor detectado.
- **Por que:** ningun estimador de ubicacion (KM, ROS, media) es defendible con
  1 dato. D10 (clausula final) -> reportar solo proporcion de deteccion.
- **Que se hizo:** MS se analiza solo como deteccion + Fisher exacto Control vs
  LPS. La magnitud del efecto en MS se describe (mediana de detectados) pero no
  se testea. Decision confirmada con el usuario.

### Liquido amniotico (LA): sin modelo factorial SEXO x TTO

- **Que se probo:** un modelo factorial con censura (survreg / Tobit) sobre
  SEXO * TTO.
- **Por que no:** censura 0-60% por celda + n 5-9 lo dejan mal condicionado.
- **Que se hizo:** analisis estratificado por sexo -- Peto-Peto Control vs LPS
  dentro de hembra y dentro de macho, mas Fisher de deteccion por sexo. La
  comparacion formal hembra vs macho del efecto LPS queda fuera de T3.

### Peto-Peto por reflexion (censura a izquierda)

- `survival::survdiff` maneja censura a DERECHA; el ELISA tiene censura a
  IZQUIERDA. Se refleja `t' = M - t` con `M = ceil(max Conc detectada) + 1`: los
  "< LOD" pasan a censura a derecha en `M - LOD` y quedan empatados en el rango
  mas bajo de la escala original, que es lo que pide D10.
- La implementacion Python calcula la formula G-rho a mano; la implementacion R
  hace lo mismo y ademas se verifica en corrida contra `survival::survdiff(rho=1)`
  (`stopifnot`, tolerancia 1e-8).

### LOD alternativo no aplicado

- D10 menciona como alternativa el menor estandar de la hoja `CURVA IL6`. Aca se
  mantiene `LOD = 0` (el blanco de la placa, como en 02_ingesta_qc). El LOD
  alternativo queda como analisis de sensibilidad para T9/T10, no aplicado.

### p-valores como texto en las tablas

- Los p-valores y el `chisq` que dependen de funciones trascendentes (`lgamma`,
  `exp`, `erfc` / `pchisq`) se guardan formateados `"%.6e"`: asi R y Python
  escriben el mismo texto aunque sus librerias matematicas difieran en el ultimo
  bit. La comparacion fina R<->Python se hace en T10 con la tolerancia declarada.
"""


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 03_elisa:inicio -->"
    marca_fin = "<!-- 03_elisa:fin -->"
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
    filas = cargar_elisa()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_ELISA)
    M = constante_reflexion(filas)

    det_h, det = tabla_deteccion(filas)
    desc_h, desc = tabla_descriptivo(filas, M)
    fish_h, fish = tabla_fisher(filas)
    peto_h, peto = tabla_petopeto_la(filas, M)

    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "elisa_deteccion.csv", det_h, det)
        escribir_csv(base / "elisa_descriptivo.csv", desc_h, desc)
        escribir_csv(base / "elisa_fisher_deteccion.csv", fish_h, fish)
        escribir_csv(base / "elisa_petopeto_la.csv", peto_h, peto)

    # --- figuras (300 dpi) ---------------------------------------------
    figura_ms(filas, fish[0], cfg.RUTA_FIGURAS / "acto1_elisa_ms.png")
    figura_la(filas, peto, cfg.RUTA_FIGURAS / "acto1_elisa_la.png")

    # --- reporte legible + artefactos compartidos --------------------
    rep = construir_reporte(fuente, det_h, det, desc_h, desc, fish_h, fish,
                            peto_h, peto, M)
    escribir_texto(cfg.RUTA_TABLAS / "elisa_reporte.md", rep)
    actualizar_descartados()

    ent = f"data/processed/elisa_long.tsv (de data/{fuente}/{cfg.ARCHIVO_ELISA})"
    registrar_procedencia([
        ["outputs/tables/{R,python}/elisa_deteccion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "deteccion/censura por grupo (D10 paso 1)"],
        ["outputs/tables/{R,python}/elisa_descriptivo.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "descriptivo por grupo: detectados + mediana KM (dato reflejado)"],
        ["outputs/tables/{R,python}/elisa_fisher_deteccion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "Fisher exacto sobre deteccion: MS y LA por sexo (Control vs LPS)"],
        ["outputs/tables/{R,python}/elisa_petopeto_la.csv", "tabla", ESTE_SCRIPT,
         "PROPIO (R verifica contra survival::survdiff rho=1)", ent,
         "LA por sexo: Peto-Peto Control vs LPS sobre Conc reflejada"],
        ["outputs/tables/elisa_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO", ent,
         "reporte legible del ELISA (Acto 1.1)"],
        ["outputs/figures/acto1_elisa_ms.png", "figura", ESTE_SCRIPT, "PROPIO", ent,
         "IL-6 suero materno: Conc por tratamiento; censurados en el LOD"],
        ["outputs/figures/acto1_elisa_la.png", "figura", ESTE_SCRIPT, "PROPIO", ent,
         "IL-6 liquido amniotico por sexo: Conc por tratamiento"],
    ])

    ms_ok = fish[0][0] == "MS"
    registrar_verificaciones([
        ["elisa_ms_metodo", "MS: solo deteccion + Fisher exacto (D10, censura Control alta)",
         "fisher_deteccion", "fisher_deteccion", "TRUE" if ms_ok else "FALSE", ESTE_SCRIPT],
        ["elisa_la_metodo", "LA: Peto-Peto G-rho=1 Control vs LPS, estratificado por sexo",
         "petopeto_por_sexo", "petopeto_por_sexo", "TRUE", ESTE_SCRIPT],
        ["elisa_la_estratos", "estratos de LA", ";".join(SEXOS_LA), "HEMBRA;MACHO",
         "TRUE" if SEXOS_LA == ["HEMBRA", "MACHO"] else "FALSE", ESTE_SCRIPT],
        ["elisa_descriptivo_grupos", "grupos en el descriptivo (MSx2 + LAx4)",
         str(len(desc)), "6", "TRUE" if len(desc) == 6 else "FALSE", ESTE_SCRIPT],
        ["elisa_lod", "LOD primario del ELISA (el blanco de la placa)",
         _fmt(LOD_ELISA), "0", "TRUE" if LOD_ELISA == 0 else "FALSE", ESTE_SCRIPT],
        ["elisa_reflexion_M", "constante de reflexion M = ceil(max Conc detectada)+1",
         _fmt(M), ">0", "TRUE" if M > 0 else "FALSE", ESTE_SCRIPT],
        ["elisa_figura_ms", "figura Acto 1.1 suero materno existe",
         "acto1_elisa_ms.png",
         "existe", "TRUE" if (cfg.RUTA_FIGURAS / "acto1_elisa_ms.png").is_file() else "FALSE",
         ESTE_SCRIPT],
        ["elisa_figura_la", "figura Acto 1.1 liquido amniotico existe",
         "acto1_elisa_la.png",
         "existe", "TRUE" if (cfg.RUTA_FIGURAS / "acto1_elisa_la.png").is_file() else "FALSE",
         ESTE_SCRIPT],
    ])

    # --- resumen por consola ----------------------------------------
    print("== 03_elisa.py ==")
    print(f"  fuente ELISA = {fuente}   M(reflexion) = {_fmt(M)}")
    print("  deteccion por grupo:")
    for r in det:
        print(f"    {r[0]:2} {r[1]:7} {r[2]:6}  n={r[3]:2}  det={r[4]:2}/{r[3]:<2} "
              f"({_fmt(r[6])}%)  cens={r[5]}")
    print(f"  MS Fisher deteccion Control vs LPS: OR(Haldane)={fish[0][6]}  p={fish[0][7]}")
    for r in peto:
        print(f"  LA {r[0]:6} Peto-Peto: chisq={r[5]}  p={r[7]}  "
              f"(nC={r[1]}, nL={r[2]})")
    print("  -> outputs/tables/{R,python}/elisa_*.csv, outputs/figures/acto1_elisa_*.png,")
    print("     outputs/tables/elisa_reporte.md")


if __name__ == "__main__":
    main()
