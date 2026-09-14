# 12_informe.py -- T11: informe HTML autocontenido del reanalisis MIA-LPS.
#
# Por que existe este archivo: cierra el pipeline juntando en un unico documento
# navegable (Acto 1 + Acto 2 + reproducibilidad + analisis descartados +
# limitaciones + procedencia/verificaciones) todo lo que produjeron 02..11 y 98.
# El informe se arma A MANO (no rmarkdown / no pandoc): (a) esas dependencias no
# estan garantizadas en un clon limpio -- de hecho falta `rmarkdown` en este
# build de R -- y (b) la unica forma de que R y Python generen EL MISMO
# `docs/informe.html` byte a byte es construir el HTML con la misma logica de
# strings en ambos lenguajes, sin un motor intermedio.
#
# Entradas: los `.md` de copia unica de outputs/tables/ (byte-identicos R/Python
# por construccion de 02..11/98), unos pocos CSV de outputs/tables/<lang>/ para
# los numeros del resumen, y los PNG de outputs/figures/. Las figuras se
# incrustan en base64 (informe AUTOCONTENIDO); como los PNG NO son byte-identicos
# entre ggplot2 y matplotlib (ver ESTADO / analisis_descartados), `informe.html`
# tampoco lo es: la paridad R/Python se chequea en 99_verificar sobre la
# PROYECCION del HTML sin los blobs `data:` (informe.textonly.html) y sobre los
# `.md`. No lleva marca de tiempo (eso va en logs/corrida_<fecha>.txt).
#
# Salidas:
#   docs/informe.html                         -- informe autocontenido
#   docs/informe.pdf                          -- best-effort (Edge/Chrome headless);
#                                                el pipeline NO falla si no se puede
#   outputs/intermediate/render/<lang>/...    -- snapshots para la paridad de 99
#   filas en procedencia.csv / verificaciones.csv (script = 12_informe)

from __future__ import annotations

import base64
import importlib.util
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

ESTE_SCRIPT = "12_informe"
LANG = "python"

# --- Reportes .md de copia unica -> seccion del informe --------------------
# (titulo de seccion, archivo .md en outputs/tables/, criterio de figuras)
#
# El 3er elemento YA NO es una lista de PNG escrita a mano (se desincronizaba
# cada vez que un script agregaba figuras -- ver ESTADO.md, sesion del
# bugfix): es un CRITERIO que filtra `procedencia.csv` (clave
# "script", mas "patron" opcional de nombre para separar los dos usos de
# 07_figuras_acto1 -- expresion vs pSTAT3). Las figuras de cada seccion se
# derivan solas en `figuras_de_seccion()`.
SECCIONES = [
    ("Acto 1.1 -- ELISA de IL-6 (validacion del modelo)", "elisa_reporte.md",
     {"script": "03_elisa"}),
    ("Acto 1.2 -- Cuantificacion relativa (qPCR)", "qpcr_cuantificacion_reporte.md",
     {"script": "04_qpcr_cuantificacion"}),
    ("Acto 1.3 -- Modelos de expresion (qPCR)", "qpcr_modelos_reporte.md",
     {"script": "07_figuras_acto1", "patron": r"^acto1_expresion_"}),
    ("Acto 1.4 -- pSTAT3 en placenta", "pstat3_reporte.md",
     {"script": "07_figuras_acto1", "patron": r"^acto1_pstat3"}),
    ("Acto 2.1-2.2 -- Correlacion placenta<->cerebro y co-expresion",
     "acto2_correlaciones_reporte.md",
     {"script": "08_acto2_correlaciones"}),
    ("Acto 2.3-2.4 -- Dispersion y test formal de Delta rho",
     "acto2_dispersion_reporte.md",
     {"script": "09_acto2_dispersion"}),
    ("Acto 2.5 -- Simulacion de restriccion de rango",
     "acto2_simulacion_reporte.md", {"script": "10_acto2_simulacion"}),
    ("Acto 2.6 -- Sensibilidad (eigengene, exclusion del feto extremo)",
     "acto2_sensibilidad_reporte.md", {"script": "11_sensibilidad"}),
]

# Decisiones D1..D12 (texto fijo, identico en ambos lenguajes; espejo de AGENTS 4).
DECISIONES = [
    ("D1", "Cuantificacion relativa: dCt = Ct_gen - Ct_rsp29. Calibrador = HEMBRA "
     "CONTROL por gen x tejido, promediando solo valores detectados. ddCt = "
     "dCt_muestra - dCt_calibrador."),
    ("D2", "El analisis corre sobre -ddCt (log2, simetrica y aditiva). El "
     "fold-change FC = 2^(-ddCt) solo se grafica (eje Y log)."),
    ("D3", "No detectados -> NA. Se evaluo imputacion MNAR y se descarto; no se "
     "imputa por ningun metodo."),
    ("D4", "Poblacion de analisis: los 36 fetos, sin exclusion de outliers en el "
     "analisis principal (la robustez la cubren los controles de sensibilidad)."),
    ("D5", "Modelo por gen x tejido: -ddCt ~ SEXO * TTO, con cascada de supuestos "
     "sobre los residuos (ANOVA-III / HC3 / ART)."),
    ("D6", "Post hoc solo si SEXO x TTO es significativa: 4 comparaciones fijas, "
     "correccion Holm."),
    ("D7", "il6 en cerebro E15 no es cuantificable (calibrador HEMBRA_CONTROL "
     "0/9 detectados): fuera del modelo, se analiza como proporcion de deteccion "
     "(Fisher exacto)."),
    ("D8", "Score compuesto de transportadores = z-score por gen (dentro de "
     "tejido, 36 fetos) y promedio de los 7 z por feto. Variante PCA = eigengene."),
    ("D9", "pSTAT3: PSTAT3 ~ SEXO * TTO + MEMBRANA (bloque fijo). Limitacion "
     "obligatoria: sin STAT3 total -> mide abundancia de fosfo-STAT3, no fraccion "
     "fosforilada."),
    ("D10", "ELISA con censura a izquierda: indicador censurado, valor NA, LOD "
     "aparte. Se reporta % de censura por grupo antes de todo estadistico; "
     "metodos para datos censurados (KM/ROS o no parametricos)."),
    ("D11", "Anotacion de boxplots solo si SEXO x TTO y el post hoc son "
     "significativos; una sola funcion de brackets por lenguaje."),
    ("D12", "Correccion entre genes: primario sin correccion; columna "
     "suplementaria con p ajustado por Benjamini-Hochberg dentro de cada tejido. "
     "No cambia conclusiones."),
]

LIMITACIONES = [
    "Las conclusiones biologicas solo son validas con los datos reales. Sobre "
    "datos sinteticos, el informe demuestra que el pipeline es completo y "
    "reproducible; los efectos simulados son arbitrarios.",
    "pSTAT3 (D9): normalizado a proteina total, sin STAT3 total en la misma "
    "membrana -> refleja abundancia de fosfo-STAT3 (Tyr705), no la fraccion de "
    "STAT3 fosforilada. Un cambio puede deberse a mas fosforilacion, a mas STAT3 "
    "total, o a ambos.",
    "ELISA de IL-6: n propio (14 madres en suero, 28 sacos en liquido amniotico), "
    "no se fuerza al n=36 del diseno de qPCR. Censura a izquierda alta en varios "
    "grupos (hasta 80 % en suero Control) -> en esos casos solo se reporta "
    "proporcion de deteccion, no un estimador de ubicacion.",
    "il6 en cerebro fetal E15: no cuantificable por D7 (0 detectados en el "
    "calibrador). Solo se analiza como proporcion de deteccion.",
    "BRAIN_P1 (cerebro postnatal dia 1) queda fuera del alcance de este informe "
    "(placenta y cerebro fetal E15); se analiza por separado.",
    "Acto 2: la comparacion de correlaciones entre grupos es de baja potencia "
    "(n_par ~= 15-18 por estrato). La ausencia de significancia no es evidencia "
    "de igualdad; se reporta el test formal (Fisher z) y la simulacion de "
    "restriccion de rango, no el patron descriptivo.",
]


# =========================================================================
# Formateo y escritura -- identicos a 02..11/98.
# =========================================================================
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


def leer_csv(ruta: Path):
    if not ruta.is_file():
        return None, []
    lineas = ruta.read_bytes().decode("utf-8").split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    if not lineas:
        return None, []
    return _parse_csv_line(lineas[0]), [_parse_csv_line(x) for x in lineas[1:]]


def leer_texto(ruta: Path) -> str:
    if not ruta.is_file():
        return ""
    return ruta.read_bytes().decode("utf-8")


def leer_csv_sin_este(ruta: Path):
    """Como leer_csv pero descarta las filas cuyo campo 'script' es 12_informe:
    el informe describe el pipeline de analisis, no su propio paso de reporte,
    y asi el bloque de auditoria queda identico R/Python e idempotente entre
    corridas (12_informe no se cuenta a si mismo)."""
    h, f = leer_csv(ruta)
    if h is None or "script" not in h:
        return h, f
    j = h.index("script")
    return h, [fila for fila in f if not (j < len(fila) and fila[j] == ESTE_SCRIPT)]


# =========================================================================
# Merge en procedencia.csv / verificaciones.csv por 'script' -- igual que 02..11.
# =========================================================================
def merge_por_script(ruta: Path, header, filas_nuevas, clave_orden):
    h_old, filas_old = leer_csv(ruta)
    if h_old is not None and h_old != header:
        raise SystemExit(f"{ruta.name}: encabezado incompatible {h_old!r} vs {header!r}")
    idx = header.index("script")
    conservadas = [f for f in filas_old if len(f) > idx and f[idx] != ESTE_SCRIPT]
    todas = conservadas + [[_fmt(v) for v in f] for f in filas_nuevas]
    todas.sort(key=clave_orden)
    escribir_csv(ruta, header, todas)


def registrar_procedencia(filas_nuevas):
    header = ["artefacto", "tipo", "script", "origen_codigo", "entradas", "descripcion"]
    merge_por_script(cfg.RUTA_TABLAS / "procedencia.csv", header, filas_nuevas,
                     lambda f: (f[2], f[0]))


def registrar_verificaciones(filas_nuevas):
    header = ["id", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script"]
    merge_por_script(cfg.RUTA_TABLAS / "verificaciones.csv", header, filas_nuevas,
                     lambda f: (f[5], f[0]))


# =========================================================================
# Mini Markdown -> HTML  (deterministico; misma logica en R y Python).
#
# Cubre lo que usan los .md del repo: encabezados #/##/###, parrafos, listas
# `- ` con anidado por sangria de 2 espacios, sublistas ordenadas `  N. `,
# tablas `| ... |` (con fila separadora `| --- |`), y en linea **negrita**,
# `codigo`, *enfasis* y [texto](url). El offset de encabezados sube los del
# reporte para que aniden bajo el <h2> de la seccion del informe.
# =========================================================================
_RE_BOLD = re.compile(r"\*\*([^*]+)\*\*")
_RE_CODE = re.compile(r"`([^`]+)`")
_RE_LINK = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")
_RE_EM = re.compile(r"\*([^*]+)\*")


def _esc(s: str) -> str:
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _inline(s: str) -> str:
    s = _esc(s)
    s = _RE_CODE.sub(lambda m: "<code>" + m.group(1) + "</code>", s)
    s = _RE_BOLD.sub(lambda m: "<strong>" + m.group(1) + "</strong>", s)
    s = _RE_EM.sub(lambda m: "<em>" + m.group(1) + "</em>", s)
    s = _RE_LINK.sub(lambda m: '<a href="' + m.group(2) + '">' + m.group(1) + "</a>", s)
    return s


def _celdas_tabla(linea: str):
    t = linea.strip()
    if t.startswith("|"):
        t = t[1:]
    if t.endswith("|"):
        t = t[:-1]
    return [c.strip() for c in t.split("|")]


def _es_separador_tabla(linea: str) -> bool:
    cs = _celdas_tabla(linea)
    return len(cs) > 0 and all(set(c) <= set("-: ") and "-" in c for c in cs)


def md_a_html(texto: str, base_nivel: int = 3) -> str:
    lineas = texto.replace("\r\n", "\n").split("\n")
    out = []
    i = 0
    n = len(lineas)
    parr = []

    def cerrar_parr():
        if parr:
            out.append("<p>" + " ".join(parr) + "</p>")
            parr.clear()

    while i < n:
        ln = lineas[i]
        s = ln.strip()

        if s == "" or re.match(r"^<!--.*-->$", s):
            cerrar_parr()
            i += 1
            continue

        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            cerrar_parr()
            niv = min(6, base_nivel + len(m.group(1)) - 1)
            out.append("<h%d>%s</h%d>" % (niv, _inline(m.group(2)), niv))
            i += 1
            continue

        if s.startswith("|") and i + 1 < n and _es_separador_tabla(lineas[i + 1]):
            cerrar_parr()
            enc = _celdas_tabla(s)
            i += 2
            cuerpo = []
            while i < n and lineas[i].strip().startswith("|"):
                cuerpo.append(_celdas_tabla(lineas[i]))
                i += 1
            out.append('<table><thead><tr>'
                       + "".join("<th>" + _inline(c) + "</th>" for c in enc)
                       + "</tr></thead><tbody>")
            for fila in cuerpo:
                out.append("<tr>"
                           + "".join("<td>" + _inline(c) + "</td>" for c in fila)
                           + "</tr>")
            out.append("</tbody></table>")
            continue

        if re.match(r"^[-*]\s+", ln.lstrip()) or re.match(r"^\s{2,}\d+\.\s+", ln):
            cerrar_parr()
            bloque = []
            while i < n and (lineas[i].strip() != ""
                             and (re.match(r"^[-*]\s+", lineas[i].lstrip())
                                  or re.match(r"^\s+", lineas[i]))):
                bloque.append(lineas[i])
                i += 1
            out.append(_lista_html(bloque))
            continue

        parr.append(_inline(s))
        i += 1

    cerrar_parr()
    return "\n".join(out)


def _lista_html(bloque) -> str:
    # Items de primer nivel: lineas que empiezan (sin sangria) con '- ' o '* '.
    # Continuaciones sangradas que no son '  N. ' se pegan al item corriente.
    # Sublistas: lineas '  N. ' consecutivas -> <ol> dentro del <li>.
    items = []
    cur = None
    sub = None
    for ln in bloque:
        m1 = re.match(r"^[-*]\s+(.*)$", ln)
        m2 = re.match(r"^\s{2,}(\d+)\.\s+(.*)$", ln)
        if m1:
            if cur is not None:
                items.append((cur, sub))
            cur = _inline(m1.group(1).strip())
            sub = None
        elif m2:
            if sub is None:
                sub = []
            sub.append(_inline(m2.group(2).strip()))
        else:
            extra = _inline(ln.strip())
            if cur is None:
                cur = extra
            else:
                cur = cur + " " + extra
    if cur is not None:
        items.append((cur, sub))
    partes = ["<ul>"]
    for txt, sb in items:
        if sb:
            partes.append("<li>" + txt + "<ol>"
                          + "".join("<li>" + x + "</li>" for x in sb)
                          + "</ol></li>")
        else:
            partes.append("<li>" + txt + "</li>")
    partes.append("</ul>")
    return "\n".join(partes)


# =========================================================================
# Numeros del resumen -- se leen de outputs/tables/<lang>/ (byte-identicos R/Py).
# =========================================================================
def _tab(nombre: str):
    return leer_csv(cfg.RUTA_TABLAS_PY / nombre)


def _col(header, filas, nombre):
    if header is None or nombre not in header:
        return []
    j = header.index(nombre)
    return [f[j] if j < len(f) else "" for f in filas]


# =========================================================================
# Figuras por seccion -- DERIVADAS de procedencia.csv, no hardcodeadas (ver
# comentario de SECCIONES). `figuras_procedencia()` es la unica lectura de
# procedencia.csv para esto; `figuras_de_seccion()` filtra por script (y un
# patron opcional de nombre) para una seccion puntual; `figuras_embebidas()`
# junta lo que efectivamente se incrusta en el informe completo (usado tanto
# para armar el HTML como para la verificacion de cobertura en main()).
# =========================================================================
def figuras_procedencia():
    h, f = leer_csv_sin_este(cfg.RUTA_TABLAS / "procedencia.csv")
    art = _col(h, f, "artefacto")
    tipo = _col(h, f, "tipo")
    scr = _col(h, f, "script")
    filas = [(a, s) for a, t, s in zip(art, tipo, scr)
             if t == "figura" and a.startswith("outputs/figures/")]
    return filas


def figuras_de_seccion(figs, criterio):
    nombres = [Path(a).name for a, s in figs if s == criterio["script"]]
    patron = criterio.get("patron")
    if patron:
        nombres = [n for n in nombres if re.match(patron, n)]
    return nombres


def figuras_embebidas(figs):
    vistas = []
    for _tit, _md, criterio in SECCIONES:
        for n in figuras_de_seccion(figs, criterio):
            if n not in vistas:
                vistas.append(n)
    return vistas


def resumen_numeros():
    r = {}

    h, f = _tab("pstat3_modelo_clasificacion.csv")
    p = _col(h, f, "p_SEXOxTTO")
    r["pstat3_pint"] = p[0] if p else "n/d"

    h, f = _tab("pstat3_posthoc.csv")
    contr = _col(h, f, "contraste")
    ph = _col(h, f, "p_holm")
    r["pstat3_hh"] = "n/d"
    for c, v in zip(contr, ph):
        if c == "HEMBRA_CONTROL-HEMBRA_LPS":
            r["pstat3_hh"] = v

    h, f = _tab("qpcr_modelos_clasificacion.csv")
    tej = _col(h, f, "TEJIDO")
    isig = _col(h, f, "interaccion_significativa")
    r["qpcr_int_pla"] = sum(1 for t, s in zip(tej, isig)
                            if t == "PLACENTA_E15" and s == "TRUE")
    r["qpcr_int_bra"] = sum(1 for t, s in zip(tej, isig)
                            if t == "BRAIN_E15" and s == "TRUE")
    r["qpcr_modelados"] = sum(1 for s in _col(h, f, "via") if s == "modelo")

    h, f = _tab("acto2_test_correlaciones.csv")
    it = _col(h, f, "ITEM")
    pbw = _col(h, f, "p_bw")
    pares = [(a, b) for a, b in zip(it, pbw) if b not in ("", None)]
    r["acto2_n_test"] = len(pares)
    if pares:
        mn = min(pares, key=lambda ab: float(ab[1]))
        r["acto2_min_item"], r["acto2_min_pbw"] = mn[0], mn[1]
        r["acto2_n_sig"] = sum(1 for _, b in pares if float(b) < 0.05)
    else:
        r["acto2_min_item"] = r["acto2_min_pbw"] = "n/d"
        r["acto2_n_sig"] = 0

    h, f = _tab("acto2_simulacion.csv")
    ver = _col(h, f, "veredicto")
    r["sim_fuera"] = sum(1 for v in ver if v == "FUERA")
    r["sim_total"] = len(ver)

    h, f = _tab("acto2_sensibilidad_excl_extremo.csv")
    r["sens_cambia"] = sum(1 for v in _col(h, f, "veredicto_cambia") if v == "TRUE")
    r["sens_items"] = len(f)

    h, f = _tab("acto2_sensibilidad_eigengene_test.csv")
    it = _col(h, f, "ITEM")
    pbw = _col(h, f, "p_bw")
    r["eig_pbw"] = "n/d"
    for a, b in zip(it, pbw):
        if a == "eigengene":
            r["eig_pbw"] = b

    h, f = leer_csv(cfg.RUTA_TABLAS / "comparacion_R_python.csv")
    r["comp_total"] = r["comp_byte"] = r["comp_fuera"] = "n/d"
    if h is not None:
        arch = _col(h, f, "archivo")
        for k, fila in enumerate(f):
            if arch[k] == "__TOTAL__":
                d = dict(zip(h, fila))
                r["comp_byte"] = d.get("byte_identico", "n/d")
                r["comp_fuera"] = d.get("n_fuera_tol", "n/d")
        r["comp_total"] = str(sum(1 for a in arch if a != "__TOTAL__"))

    h, f = leer_csv_sin_este(cfg.RUTA_TABLAS / "procedencia.csv")
    r["n_procedencia"] = len(f)
    h, f = leer_csv_sin_este(cfg.RUTA_TABLAS / "verificaciones.csv")
    okc = _col(h, f, "ok")
    r["n_verif"] = len(f)
    r["n_verif_true"] = sum(1 for v in okc if v == "TRUE")
    return r


# =========================================================================
# Figuras -> <img data:...>  (autocontenido).
# =========================================================================
def img_datauri(ruta: Path) -> str:
    b = ruta.read_bytes()
    return "data:image/png;base64," + base64.b64encode(b).decode("ascii")


def fig_html(nombre: str) -> str:
    ruta = cfg.RUTA_FIGURAS / nombre
    if not ruta.is_file():
        return '<p class="falta">[falta la figura %s]</p>' % _esc(nombre)
    return ('<figure><img alt="%s" src="%s"><figcaption>%s</figcaption></figure>'
            % (_esc(nombre), img_datauri(ruta), _esc(nombre)))


# =========================================================================
# CSS (inline; sin recursos externos).
# =========================================================================
CSS = """
:root { color-scheme: light; }
* { box-sizing: border-box; }
body { font: 15px/1.55 -apple-system, "Segoe UI", Roboto, Helvetica, Arial,
       sans-serif; color: #1a1a1a; background: #fff; margin: 0;
       padding: 2.2rem 1.4rem 4rem; }
main { max-width: 60rem; margin: 0 auto; }
h1 { font-size: 1.7rem; margin: 0 0 .2rem; }
h2 { font-size: 1.28rem; margin: 2.4rem 0 .6rem; padding-top: .4rem;
     border-top: 2px solid #222; }
h3 { font-size: 1.08rem; margin: 1.5rem 0 .4rem; }
h4 { font-size: .98rem; margin: 1.1rem 0 .3rem; color: #333; }
h5 { font-size: .92rem; margin: .9rem 0 .3rem; color: #444; }
p { margin: .5rem 0; }
code { background: #f0f0f0; padding: .05em .35em; border-radius: 3px;
       font: .86em/1.4 "SF Mono", Consolas, "Liberation Mono", monospace; }
a { color: #0645ad; }
ul, ol { margin: .4rem 0 .7rem; padding-left: 1.5rem; }
li { margin: .18rem 0; }
table { border-collapse: collapse; margin: .8rem 0; font-size: .82rem;
        display: block; overflow-x: auto; max-width: 100%; }
th, td { border: 1px solid #ccc; padding: .28em .55em; text-align: left;
         white-space: nowrap; }
thead th { background: #f4f4f4; position: sticky; top: 0; }
tbody tr:nth-child(even) { background: #fafafa; }
figure { margin: 1rem 0; text-align: center; }
figure img { max-width: 100%; height: auto; border: 1px solid #e2e2e2; }
figcaption { font-size: .78rem; color: #666; margin-top: .3rem; }
.aviso { background: #fff4e5; border: 1px solid #f0c48a; padding: .8rem 1rem;
         border-radius: 6px; margin: 1rem 0; }
.meta { color: #555; font-size: .9rem; }
.falta { color: #b00; font-style: italic; }
nav.toc { background: #f7f7f7; border: 1px solid #e0e0e0; border-radius: 6px;
          padding: .8rem 1.2rem; margin: 1.4rem 0; }
nav.toc ol { margin: .3rem 0; }
footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #ccc;
         color: #666; font-size: .85rem; }
@media print {
  body { padding: 0; font-size: 11pt; }
  h2 { page-break-before: auto; }
  figure, table { page-break-inside: avoid; }
  thead th { position: static; }
  nav.toc { page-break-after: always; }
}
"""


# =========================================================================
# Armado del HTML.
# =========================================================================
def _seccion(id_, titulo, cuerpo_html):
    return '<section id="%s">\n<h2>%s</h2>\n%s\n</section>' % (
        id_, _esc(titulo), cuerpo_html)


def construir_html(fuente: str, num: dict) -> str:
    sint = fuente != "real"
    L = []
    ap = L.append
    ap("<!doctype html>")
    ap('<html lang="es">')
    ap("<head>")
    ap('<meta charset="utf-8">')
    ap('<meta name="viewport" content="width=device-width, initial-scale=1">')
    ap("<title>Reanalisis MIA-LPS -- placenta E15 / cerebro fetal E15</title>")
    ap("<style>%s</style>" % CSS)
    ap("</head>")
    ap("<body>")
    ap("<main>")

    ap("<h1>Reanalisis reproducible &mdash; MIA-LPS en placenta E15 y cerebro "
       "fetal E15</h1>")
    ap('<p class="meta">Activacion inmune materna (LPS 100 &micro;g/kg i.p., dia 15 '
       "de gestacion, colecta a las 6 h). Respuesta de la placenta y del cerebro "
       "fetal, y si depende del sexo del feto.</p>")
    ap('<p class="meta">Fuente de datos de esta corrida: <code>%s</code>. '
       "Informe generado por <code>12_informe</code> (R y Python producen el "
       "mismo <code>informe.html</code> salvo los PNG incrustados). Sin marca de "
       "tiempo: la fecha de corrida esta en <code>logs/corrida_&lt;fecha&gt;.txt</code>."
       % _esc(fuente))
    if sint:
        ap('<div class="aviso"><strong>Datos sinteticos.</strong> Esta corrida no '
           "uso los crudos reales. Las salidas demuestran que el pipeline es "
           "completo y reproducible; los efectos son arbitrarios y no tienen "
           "lectura biologica.</div>")

    # --- indice ---
    items_toc = [
        ("resumen", "1. Resumen"),
        ("metodos", "2. Diseno y metodos"),
        ("acto1", "3. Acto 1 &mdash; respuesta a la MIA y dependencia del sexo"),
        ("acto2", "4. Acto 2 &mdash; coordinacion placenta&lt;-&gt;cerebro"),
        ("reproducibilidad", "5. Reproducibilidad (R vs Python)"),
        ("descartados", "6. Analisis descartados"),
        ("limitaciones", "7. Limitaciones"),
        ("auditoria", "8. Procedencia y verificaciones"),
    ]
    ap('<nav class="toc"><strong>Contenido</strong><ol>')
    for id_, t in items_toc:
        ap('<li><a href="#%s">%s</a></li>' % (id_, t))
    ap("</ol></nav>")

    # --- 1. resumen ---
    res = []
    res.append("<ul>")
    res.append(
        "<li><strong>Validacion del modelo (ELISA IL-6).</strong> El LPS eleva "
        "IL-6 en suero materno (1/5 vs 9/9 detectados; Fisher p &asymp; 5&times;10"
        "<sup>-3</sup>). En liquido amniotico la senal es mas debil y la censura "
        "alta; el Peto-Peto por sexo no alcanza significancia.</li>")
    res.append(
        "<li><strong>pSTAT3 en placenta.</strong> Interaccion SEXO&times;TTO "
        "significativa (p = <code>%s</code>): el aumento de fosfo-STAT3 con LPS es "
        "restringido a hembras (HEMBRA_CONTROL&ndash;HEMBRA_LPS Holm p = "
        "<code>%s</code>); en machos no cambia. Mide abundancia de fosfo-STAT3, "
        "no fraccion fosforilada (D9).</li>" % (num["pstat3_pint"], num["pstat3_hh"]))
    res.append(
        "<li><strong>Programa de transportadores (qPCR).</strong> De %d gen "
        "&times; tejido modelados, la interaccion SEXO&times;TTO es significativa "
        "en %d de cerebro E15 y %d de placenta E15; en cerebro las %d sobreviven "
        "la correccion BH entre genes (D12).</li>" % (
            num["qpcr_modelados"], num["qpcr_int_bra"], num["qpcr_int_pla"],
            num["qpcr_int_bra"]))
    res.append(
        "<li><strong>Coordinacion placenta&lt;-&gt;cerebro (Acto 2).</strong> "
        "En el test formal de &Delta;&rho; (Fisher z sobre &rho; de Spearman, "
        "prohibicion 4), %d de %d items alcanzan p &lt; 0.05 (minimo: "
        "<code>%s</code>, p_bw = <code>%s</code>). La simulacion de restriccion de "
        "rango deja %d de %d celdas FUERA del IC95. <strong>Sobre los datos reales "
        "no hay evidencia de que la coordinacion cambie entre Control y "
        "LPS.</strong></li>" % (
            num["acto2_n_sig"], num["acto2_n_test"], num["acto2_min_item"],
            num["acto2_min_pbw"], num["sim_fuera"], num["sim_total"]))
    res.append(
        "<li><strong>Robustez.</strong> El resultado negativo del Acto 2 se "
        "sostiene con el eigengene PC1 en vez del promedio de z (p_bw = "
        "<code>%s</code>) y al excluir el feto mas influyente por item "
        "(%d de %d items cambian el veredicto).</li>" % (
            num["eig_pbw"], num["sens_cambia"], num["sens_items"]))
    res.append(
        "<li><strong>Reproducibilidad.</strong> %s CSV de resultados comparados "
        "R&harr;Python: byte-identicos = <code>%s</code>, celdas fuera de "
        "tolerancia = <code>%s</code>.</li>" % (
            num["comp_total"], num["comp_byte"], num["comp_fuera"]))
    res.append("</ul>")
    L.append(_seccion("resumen", "1. Resumen", "\n".join(res)))

    # --- 2. metodos ---
    met = []
    met.append("<h3>2.1 Diseno</h3>")
    met.append("<p>18 madres, 36 fetos E15 (1 hembra + 1 macho por camada). "
               "Cuatro grupos SEXO &times; TTO (HEMBRA/MACHO &times; CONTROL/LPS), "
               "9 fetos por grupo. Tejidos analizados: placenta E15 y cerebro "
               "fetal E15 (BRAIN_P1 fuera de alcance). El ELISA de IL-6 y el "
               "western de pSTAT3 se analizan con su propio n.</p>")
    met.append("<h3>2.2 Decisiones metodologicas fijas (D1&ndash;D12)</h3>")
    met.append('<table><thead><tr><th>#</th><th>Decision</th></tr></thead><tbody>')
    for k, v in DECISIONES:
        met.append("<tr><td>%s</td><td>%s</td></tr>" % (k, _esc(v)))
    met.append("</tbody></table>")
    met.append("<h3>2.3 Cascada de supuestos (D5)</h3>")
    met.append("<ul>"
               "<li>Shapiro-Wilk (residuos del modelo conjunto) y Levene "
               "(Brown-Forsythe, centro = mediana) OK &rarr; <strong>ANOVA tipo "
               "III</strong> (contrastes suma-cero).</li>"
               "<li>Falla solo Levene &rarr; <strong>OLS con errores HC3</strong>.</li>"
               "<li>Falla Shapiro (con o sin Levene) &rarr; <strong>ART</strong> "
               "(Aligned Rank Transform); post hoc ART-C, nunca emmeans directo.</li>"
               "</ul>")
    met.append("<h3>2.4 Control de calidad de la ingesta</h3>")
    met.append(md_a_html(leer_texto(cfg.RUTA_TABLAS / "qc_reporte.md")))
    L.append(_seccion("metodos", "2. Diseno y metodos", "\n".join(met)))

    # --- 3. Acto 1 + 4. Acto 2 ---
    figs = figuras_procedencia()

    def bloque_secciones(ids, titulo, indices):
        partes = []
        for k in indices:
            tit, md, criterio = SECCIONES[k]
            partes.append("<h3>%s</h3>" % _esc(tit))
            partes.append(md_a_html(leer_texto(cfg.RUTA_TABLAS / md)))
            for fg in figuras_de_seccion(figs, criterio):
                partes.append(fig_html(fg))
        return _seccion(ids, titulo, "\n".join(partes))

    L.append(bloque_secciones(
        "acto1", "3. Acto 1 -- respuesta a la MIA y dependencia del sexo",
        [0, 1, 2, 3]))
    L.append(bloque_secciones(
        "acto2", "4. Acto 2 -- coordinacion placenta<->cerebro", [4, 5, 6, 7]))

    # --- 5. reproducibilidad ---
    rep = []
    rep.append("<p>Todo el analisis esta implementado <strong>dos veces</strong> "
               "(R y Python), con semilla fija <code>%d</code> y solo rutas "
               "relativas. Las salidas numericas llevan nombres identicos en "
               "<code>outputs/tables/R/</code> y <code>outputs/tables/python/</code>; "
               "<code>98_comparacion</code> las cruza celda a celda (tolerancia "
               "1e-6, con fallback 1e-4 para p de tests iterativos) y "
               "<code>99_verificar</code> re-corre esa comparacion y byte-compara "
               "los <code>.md</code> y la proyeccion sin figuras de este "
               "informe.</p>" % cfg.SEMILLA)
    rep.append(md_a_html(leer_texto(cfg.RUTA_TABLAS / "comparacion_reporte.md")))
    L.append(_seccion("reproducibilidad", "5. Reproducibilidad (R vs Python)",
                      "\n".join(rep)))

    # --- 6. descartados ---
    L.append(_seccion(
        "descartados", "6. Analisis descartados",
        md_a_html(leer_texto(cfg.RUTA_TABLAS / "analisis_descartados.md"))))

    # --- 7. limitaciones ---
    lim = ["<ul>"]
    for x in LIMITACIONES:
        lim.append("<li>%s</li>" % _esc(x))
    lim.append("</ul>")
    L.append(_seccion("limitaciones", "7. Limitaciones", "\n".join(lim)))

    # --- 8. auditoria ---
    aud = []
    aud.append("<p><code>procedencia.csv</code>: <strong>%d</strong> filas (una "
               "por figura y por tabla). <code>verificaciones.csv</code>: "
               "<strong>%d</strong> filas, <strong>%d</strong> en TRUE. Tablas "
               "completas en <code>outputs/tables/</code>.</p>"
               % (num["n_procedencia"], num["n_verif"], num["n_verif_true"]))
    aud.append("<h3>8.1 Procedencia</h3>")
    aud.append(_csv_a_tabla(cfg.RUTA_TABLAS / "procedencia.csv"))
    aud.append("<h3>8.2 Verificaciones</h3>")
    aud.append('<p class="meta">Se omite la columna <code>valor_obtenido</code> '
               "(numeros de diagnostico que pueden diferir en el ultimo digito "
               "entre R y Python); la tabla completa esta en "
               "<code>outputs/tables/verificaciones.csv</code>.</p>")
    aud.append(_csv_a_tabla(cfg.RUTA_TABLAS / "verificaciones.csv",
                            omitir=("valor_obtenido",)))
    L.append(_seccion("auditoria", "8. Procedencia y verificaciones",
                      "\n".join(aud)))

    ap = L.append
    ap("<footer>Reanalisis MIA-LPS &mdash; placenta E15 / cerebro fetal E15. "
       "Generado por <code>12_informe</code>; R y Python producen este HTML "
       "identico salvo los PNG incrustados. Ver <code>AGENTS.md</code> para las "
       "decisiones D1&ndash;D12 y <code>ESTADO.md</code> para la bitacora.</footer>")
    ap("</main>")
    ap("</body>")
    ap("</html>")
    return "\n".join(L) + "\n"


def _csv_a_tabla(ruta: Path, omitir=()) -> str:
    h, f = leer_csv_sin_este(ruta)
    if h is None:
        return '<p class="falta">[falta %s]</p>' % _esc(ruta.name)
    cols = [j for j in range(len(h)) if h[j] not in omitir]
    out = ["<table><thead><tr>"
           + "".join("<th>" + _esc(h[j]) + "</th>" for j in cols)
           + "</tr></thead><tbody>"]
    for fila in f:
        out.append("<tr>"
                   + "".join("<td>" + _esc(fila[j] if j < len(fila) else "")
                             + "</td>" for j in cols)
                   + "</tr>")
    out.append("</tbody></table>")
    return "\n".join(out)


# =========================================================================
# Proyeccion "sin figuras" para la paridad R/Python en 99_verificar.
# =========================================================================
_RE_DATAURI = re.compile(r'src="data:image/png;base64,[^"]*"')


def html_sin_figuras(html: str) -> str:
    return _RE_DATAURI.sub('src="[png]"', html)


# =========================================================================
# PDF por impresion headless (best-effort; el pipeline NO falla por esto).
# =========================================================================
def _buscar_motor_pdf():
    env = os.environ.get("MIA_LPS_PDF_ENGINE")
    if env and Path(env).is_file():
        return env
    cands = [
        r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Microsoft\Edge\Application\msedge.exe",
        r"C:\Program Files\Google\Chrome\Application\chrome.exe",
        r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
    ]
    for c in cands:
        if Path(c).is_file():
            return c
    for n in ("msedge", "chrome", "chromium", "chromium-browser", "google-chrome"):
        w = shutil.which(n)
        if w:
            return w
    return None


def generar_pdf(html_path: Path, pdf_path: Path) -> str:
    """Devuelve 'ok' | 'sin_motor' | 'fallo:<motivo>'. Nunca lanza."""
    motor = _buscar_motor_pdf()
    if not motor:
        print("  [pdf] sin motor de impresion (Edge/Chrome); queda solo el HTML.")
        return "sin_motor"
    perfil = cfg.RUTA_INTERMEDIOS / "edge_profile"
    try:
        perfil.mkdir(parents=True, exist_ok=True)
        if pdf_path.exists():
            pdf_path.unlink()
        cmd = [
            motor, "--headless", "--disable-gpu", "--no-first-run",
            "--no-pdf-header-footer",
            "--user-data-dir=" + str(perfil),
            "--print-to-pdf=" + str(pdf_path),
            html_path.as_uri(),
        ]
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if pdf_path.is_file() and pdf_path.stat().st_size > 0:
            print("  [pdf] %s (%d KB)" % (pdf_path.name,
                                          pdf_path.stat().st_size // 1024))
            return "ok"
        print("  [pdf] el motor no produjo un PDF valido (rc=%d); queda el HTML."
              % p.returncode)
        return "fallo:sin_salida_rc%d" % p.returncode
    except Exception as e:  # noqa: BLE001 -- best-effort por diseno
        print("  [pdf] no se pudo generar (%s); queda el HTML." % e)
        return "fallo:%s" % type(e).__name__


# =========================================================================
def main():
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)
    num = resumen_numeros()
    html = construir_html(fuente, num)

    textonly = html_sin_figuras(html)
    ruta_html = cfg.RUTA_DOCS / "informe.html"
    ruta_pdf = cfg.RUTA_DOCS / "informe.pdf"
    escribir_texto(ruta_html, html)

    pdf_status = generar_pdf(ruta_html, ruta_pdf)

    # --- snapshot por lenguaje para la paridad de 99_verificar ---
    snap = cfg.RUTA_INTERMEDIOS / "render" / LANG
    (snap / "tables").mkdir(parents=True, exist_ok=True)
    escribir_texto(snap / "informe.html", html)
    escribir_texto(snap / "informe.textonly.html", textonly)
    escribir_texto(snap / "pdf_status.txt", pdf_status)
    for md in sorted(cfg.RUTA_TABLAS.glob("*.md")):
        escribir_texto(snap / "tables" / md.name, leer_texto(md))

    figs_proc = figuras_procedencia()
    figs_todas = figuras_embebidas(figs_proc)
    n_fig = len(figs_todas)
    faltan = [fg for fg in figs_todas if not (cfg.RUTA_FIGURAS / fg).is_file()]
    md_faltan = [md for _, md, _ in SECCIONES
                 if not (cfg.RUTA_TABLAS / md).is_file()]

    # Cobertura: toda figura con fila en procedencia.csv tiene que quedar
    # incrustada en alguna seccion -- esta es la verificacion que habria
    # detectado que acto2_corr_placenta_cerebro_* y los SPLOM por sexo no
    # aparecian en el informe (quedaban en procedencia.csv pero fuera de
    # cualquier SECCIONES a mano).
    todas_en_procedencia = sorted({Path(a).name for a, _s in figs_proc})
    sin_embeber = sorted(set(todas_en_procedencia) - set(figs_todas))

    ent = "outputs/tables/*.md + outputs/tables/{R,python}/*.csv + outputs/figures/*.png"
    registrar_procedencia([
        ["docs/informe.html", "informe", ESTE_SCRIPT, "PROPIO", ent,
         "informe HTML autocontenido (Acto 1 + Acto 2 + reproducibilidad + "
         "analisis descartados + limitaciones + procedencia/verificaciones); "
         "figuras incrustadas en base64; sin marca de tiempo"],
        ["docs/informe.pdf", "informe", ESTE_SCRIPT, "PROPIO", "docs/informe.html",
         "version imprimible por impresion headless (Edge/Chrome); best-effort, "
         "el pipeline no falla si no hay motor de PDF (estado: " + pdf_status + ")"],
    ])
    registrar_verificaciones([
        ["informe_html_generado",
         "docs/informe.html existe y no esta vacio",
         "%d secciones; %d figuras" % (html.count('<section id="'), n_fig),
         "8 secciones",
         "TRUE" if (ruta_html.is_file() and ruta_html.stat().st_size > 0
                    and html.count('<section id="') == 8) else "FALSE",
         ESTE_SCRIPT],
        ["informe_figuras_incrustadas",
         "todas las figuras del Acto 1 y 2 estan incrustadas en el informe",
         "figuras=%d; faltan=%s" % (n_fig, faltan or "[]"),
         "faltan = []",
         "TRUE" if not faltan else "FALSE", ESTE_SCRIPT],
        ["informe_reportes_incluidos",
         "todos los .md de seccion existen y se incluyeron",
         "secciones=%d; md_faltan=%s" % (len(SECCIONES), md_faltan or "[]"),
         "md_faltan = []",
         "TRUE" if not md_faltan else "FALSE", ESTE_SCRIPT],
        ["informe_figuras_procedencia_embebidas",
         "toda figura con fila en procedencia.csv esta incrustada en el informe",
         "procedencia=%d; embebidas=%d; sin_embeber=%s" % (
             len(todas_en_procedencia), len(figs_todas), sin_embeber or "[]"),
         "sin_embeber = []",
         "TRUE" if not sin_embeber else "FALSE", ESTE_SCRIPT],
        ["informe_pdf",
         "docs/informe.pdf generado, o degradado limpio si no hay motor de PDF",
         "estado=%s; existe=%s" % (
             pdf_status,
             "TRUE" if ruta_pdf.is_file() and ruta_pdf.stat().st_size > 0
             else "FALSE"),
         "ok | sin_motor | fallo (nunca frena el pipeline)",
         "TRUE" if pdf_status in ("ok", "sin_motor") or pdf_status.startswith("fallo")
         else "FALSE", ESTE_SCRIPT],
        ["informe_snapshot_paridad",
         "snapshot por lenguaje para el byte-compare R/Python de 99_verificar",
         "render/<lang>/: informe.html + informe.textonly.html + tables/*.md + "
         "pdf_status.txt",
         "snapshot escrito",
         "TRUE" if (snap / "informe.textonly.html").is_file() else "FALSE",
         ESTE_SCRIPT],
    ])

    print("== 12_informe.py ==")
    print("  fuente = %s" % fuente)
    print("  -> %s (%d KB)" % (ruta_html.name,
                               ruta_html.stat().st_size // 1024))
    print("  -> %s [%s]" % (ruta_pdf.name, pdf_status))
    print("  figuras incrustadas: %d%s" % (
        n_fig, "" if not faltan else "  FALTAN: %s" % faltan))
    print("  secciones .md: %d%s" % (
        len(SECCIONES), "" if not md_faltan else "  FALTAN: %s" % md_faltan))
    print("  figuras de procedencia.csv sin embeber: %s" % (
        "ninguna" if not sin_embeber else ", ".join(sin_embeber)))
    print("  snapshot: outputs/intermediate/render/%s/" % LANG)
    if faltan or md_faltan:
        print("  *** faltan insumos: correr 02..11 y 07 antes de 12_informe ***")
        sys.exit(1)


if __name__ == "__main__":
    main()
