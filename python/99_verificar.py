# 99_verificar.py -- T11: verificacion final del reanalisis reproducible MIA-LPS.
#
# Por que existe este archivo: es la compuerta del proyecto. Partiendo de un
# clon limpio, `.\run_all.ps1` corre todo y este script comprueba, sin acceso a
# los datos crudos, que:
#   1. existe y no esta vacio cada item del checklist de AGENTS.md 1;
#   2. las salidas numericas R y Python concuerdan celda a celda dentro de la
#      tolerancia declarada (re-corre la comparacion de 98_comparacion en
#      proceso, sobre outputs/tables/{R,python}/*.csv);
#   3. los `.md` de copia unica de outputs/tables/ y la proyeccion sin figuras
#      de docs/informe.html son BYTE-IDENTICOS entre la corrida R y la Python
#      (via los snapshots que deja 12_informe en outputs/intermediate/render/);
#   4. ninguna fila de verificaciones.csv quedo distinta de TRUE.
# Escribe logs/corrida_<fecha>.txt (fecha, entorno, versiones, semilla, fuente,
# concordancia, cuantas verificaciones pasaron) e imprime
# `TODAS LAS VERIFICACIONES PASARON` si y solo si todo lo duro pasa.
#
# NO hace analisis. Solo lee salidas. `docs/informe.pdf` es blando: si 12_informe
# no encontro motor de PDF, su ausencia no frena la verificacion (AGENTS 8).
#
# PARIDAD: los helpers de lectura/comparacion son copia de 98_comparacion (la
# convencion del repo es que cada script importe solo 00_config, no que se
# cross-importen los numerados).

from __future__ import annotations

import datetime
import importlib.util
import math
import platform
import sys
from pathlib import Path

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

LANG = "python"
TOL_EST = cfg.TOL_ESTADISTICO
TOL_P = cfg.TOL_P_ITERATIVO

# --- Checklist de AGENTS.md 1 (existencia + no-vacuidad) ------------------
SCRIPTS_NUM = [
    "00_config", "01_generar_sinteticos", "02_ingesta_qc", "03_elisa",
    "04_qpcr_cuantificacion", "05_qpcr_modelos", "06_pstat3", "07_figuras_acto1",
    "08_acto2_correlaciones", "09_acto2_dispersion", "10_acto2_simulacion",
    "11_sensibilidad", "98_comparacion", "12_informe", "99_verificar",
]
FIG_ACTO1 = [
    "acto1_expresion_PLACENTA_E15.png", "acto1_expresion_BRAIN_E15.png",
    "acto1_pstat3.png", "acto1_elisa_ms.png", "acto1_elisa_la.png",
]
FIG_ACTO2 = [
    "acto2_dispersion_placenta_cerebro.png",
    "acto2_coexpresion_SPLOM_PLACENTA_E15.png",
    "acto2_coexpresion_SPLOM_BRAIN_E15.png",
    "acto2_dispersion_sd.png", "acto2_test_delta_rho.png",
    "acto2_simulacion_delta_rho.png",
    "acto2_sensibilidad_eigengene.png", "acto2_sensibilidad_excl_extremo.png",
]
SINTETICOS = [
    "Raw data CTs.tsv", "pstat3 placenta.tsv",
    "ELISA IL6 2026 Dosis 100__Sueros y LA.tsv",
    "ELISA IL6 2026 Dosis 100__CURVA IL6.tsv", "MANIFEST.tsv",
]
REPORTES_MD = [
    "qc_reporte.md", "elisa_reporte.md", "qpcr_cuantificacion_reporte.md",
    "qpcr_modelos_reporte.md", "pstat3_reporte.md",
    "acto2_correlaciones_reporte.md", "acto2_dispersion_reporte.md",
    "acto2_simulacion_reporte.md", "acto2_sensibilidad_reporte.md",
    "comparacion_reporte.md", "analisis_descartados.md",
]

PY_PAQUETES = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "seaborn",
               "sklearn", "pingouin", "openpyxl"]


# =========================================================================
# Lectura de CSV -- copia de 98_comparacion (sin dependencias).
# =========================================================================
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


# =========================================================================
# Comparacion numerica celda a celda -- copia de 98_comparacion.
# =========================================================================
def _es_num(s):
    if s is None or s == "":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def _dentro(a, b, tol):
    d = abs(a - b)
    if d <= tol:
        return True
    m = max(abs(a), abs(b))
    return m > 0.0 and d / m <= tol


def comparar_archivo(nombre: str, ruta_r: Path, ruta_py: Path) -> dict:
    hr, fr = _leer_csv(ruta_r)
    hp, fp = _leer_csv(ruta_py)
    byte_id = (ruta_r.is_file() and ruta_py.is_file()
               and ruta_r.read_bytes() == ruta_py.read_bytes())
    header_igual = hr == hp
    filas_ok = len(fr) == len(fp)
    n_num = n_dif_txt = n_fuera = 0
    max_abs = 0.0
    peor = ""
    if header_igual and filas_ok and hr is not None:
        for i, (rr, rp) in enumerate(zip(fr, fp)):
            for j in range(len(hr)):
                cr = rr[j] if j < len(rr) else ""
                cp = rp[j] if j < len(rp) else ""
                va, vb = _es_num(cr), _es_num(cp)
                if va is not None and vb is not None:
                    n_num += 1
                    d = abs(va - vb)
                    if d > max_abs:
                        max_abs = d
                        peor = "fila %d / %s" % (i + 2, hr[j])
                    if not _dentro(va, vb, TOL_EST) and not _dentro(va, vb, TOL_P):
                        n_fuera += 1
                elif cr != cp:
                    n_dif_txt += 1
                    if not peor:
                        peor = "fila %d / %s (texto)" % (i + 2, hr[j])
    ok = (header_igual and filas_ok and n_dif_txt == 0 and n_fuera == 0
          and hr is not None)
    return {"archivo": nombre, "ok": ok, "header_igual": header_igual,
            "filas_ok": filas_ok, "n_num": n_num, "n_dif_texto": n_dif_txt,
            "n_fuera_tol": n_fuera, "max_dif_abs": max_abs, "peor": peor,
            "byte_identico": byte_id}


# =========================================================================
# Utilidades.
# =========================================================================
def no_vacio(ruta: Path) -> bool:
    return ruta.is_file() and ruta.stat().st_size > 0


def _ver_paquetes():
    try:
        import importlib.metadata as _md
    except ImportError:  # pragma: no cover
        return {}
    out = {}
    for p in PY_PAQUETES:
        try:
            out[p] = _md.version(p if p != "sklearn" else "scikit-learn")
        except Exception:  # noqa: BLE001
            out[p] = "?"
    return out


# =========================================================================
def main():
    dur = []          # lineas de fallo duro
    warn = []         # avisos que NO frenan
    chk_ok = 0
    chk_tot = 0

    def chequear(desc: str, cond: bool, blando: bool = False):
        nonlocal chk_ok, chk_tot
        chk_tot += 1
        if cond:
            chk_ok += 1
        elif blando:
            warn.append(desc)
        else:
            dur.append(desc)
        return cond

    raiz = cfg.RAIZ_REPO
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)

    # --- 1. Checklist de AGENTS 1 -----------------------------------------
    for f in ("AGENTS.md", "CLAUDE.md", "README.md", "ESTADO.md",
              "run_all.ps1", "requirements.txt", "renv.lock"):
        chequear("existe/no vacio: %s" % f, no_vacio(raiz / f))

    for s in SCRIPTS_NUM:
        chequear("script R/%s.R" % s, no_vacio(raiz / "R" / ("%s.R" % s)))
        chequear("script python/%s.py" % s, no_vacio(raiz / "python" / ("%s.py" % s)))

    for t in SINTETICOS:
        chequear("dato sintetico: %s" % t, no_vacio(cfg.RUTA_DATOS_SINT / t))

    for fg in FIG_ACTO1 + FIG_ACTO2:
        chequear("figura: %s" % fg, no_vacio(cfg.RUTA_FIGURAS / fg))

    for m in REPORTES_MD:
        chequear("reporte: outputs/tables/%s" % m, no_vacio(cfg.RUTA_TABLAS / m))

    for t in ("procedencia.csv", "verificaciones.csv", "comparacion_R_python.csv"):
        chequear("auditoria: outputs/tables/%s" % t, no_vacio(cfg.RUTA_TABLAS / t))

    csv_r = sorted(p.name for p in cfg.RUTA_TABLAS_R.glob("*.csv"))
    csv_py = sorted(p.name for p in cfg.RUTA_TABLAS_PY.glob("*.csv"))
    chequear("outputs/tables/R/ tiene CSV", len(csv_r) > 0)
    chequear("outputs/tables/python/ tiene CSV", len(csv_py) > 0)

    chequear("informe: docs/informe.html", no_vacio(cfg.RUTA_DOCS / "informe.html"))

    # --- PDF: blando si 12_informe no tenia motor ------------------------
    pdf_status = "?"
    for lg in ("python", "R"):
        p = cfg.RUTA_INTERMEDIOS / "render" / lg / "pdf_status.txt"
        if p.is_file():
            pdf_status = p.read_text(encoding="utf-8").strip()
    pdf_ok = no_vacio(cfg.RUTA_DOCS / "informe.pdf")
    chequear("informe: docs/informe.pdf (blando si no hay motor de PDF)",
             pdf_ok or pdf_status != "ok", blando=not pdf_ok)

    # --- 2. Concordancia numerica R <-> Python --------------------------
    comunes = sorted(set(csv_r) & set(csv_py))
    solo_r = sorted(set(csv_r) - set(csv_py))
    solo_py = sorted(set(csv_py) - set(csv_r))
    chequear("sin CSV huerfanos entre R/ y python/", not solo_r and not solo_py)
    comps = [comparar_archivo(n, cfg.RUTA_TABLAS_R / n, cfg.RUTA_TABLAS_PY / n)
             for n in comunes]
    n_byte = sum(1 for d in comps if d["byte_identico"])
    n_fuera = sum(d["n_fuera_tol"] for d in comps)
    n_dif_txt = sum(d["n_dif_texto"] for d in comps)
    peor_abs = max((d["max_dif_abs"] for d in comps), default=0.0)
    malos = [d["archivo"] for d in comps if not d["ok"]]
    chequear("concordancia numerica R<->Python (%d CSV, tol 1e-6/1e-4)" % len(comps),
             not malos and n_fuera == 0 and n_dif_txt == 0)

    # --- 3. Paridad de render R <-> Python (.md + informe sin figuras) --
    ren_r = cfg.RUTA_INTERMEDIOS / "render" / "R"
    ren_py = cfg.RUTA_INTERMEDIOS / "render" / "python"
    paridad_render = "n/d"
    if ren_r.is_dir() and ren_py.is_dir():
        difs = []
        a = ren_r / "informe.textonly.html"
        b = ren_py / "informe.textonly.html"
        if not (a.is_file() and b.is_file() and a.read_bytes() == b.read_bytes()):
            difs.append("informe.textonly.html")
        for md in sorted((ren_r / "tables").glob("*.md")):
            bb = ren_py / "tables" / md.name
            if not (bb.is_file() and md.read_bytes() == bb.read_bytes()):
                difs.append("tables/" + md.name)
        # los .md snapshot deben coincidir con los .md en disco (ultima corrida)
        for md in sorted((ren_py / "tables").glob("*.md")):
            live = cfg.RUTA_TABLAS / md.name
            if not (live.is_file() and live.read_bytes() == md.read_bytes()):
                difs.append("disco != snapshot: " + md.name)
        paridad_render = "OK" if not difs else ("DIFIEREN: " + ", ".join(difs))
        chequear("paridad R<->Python de .md e informe (sin figuras)", not difs)
    else:
        paridad_render = "sin contraparte (corrida de un solo lenguaje)"
        warn.append("paridad de render no evaluada: falta render/R o render/python "
                    "(corrida de un solo lenguaje)")

    # --- 4. verificaciones.csv todas en TRUE ---------------------------
    hv, fv = _leer_csv(cfg.RUTA_TABLAS / "verificaciones.csv")
    v_no_true = []
    if hv is not None and "ok" in hv:
        j = hv.index("ok")
        jd = hv.index("id") if "id" in hv else 0
        v_no_true = [f[jd] for f in fv if not (j < len(f) and f[j] == "TRUE")]
    chequear("verificaciones.csv: %d filas, todas TRUE" % (len(fv)),
             hv is not None and not v_no_true)

    # --- 5. log de corrida -------------------------------------------
    fecha = datetime.date.today().isoformat()
    ver_pkgs = _ver_paquetes()
    seccion = []
    seccion.append("## Corrida %s -- %s" % (LANG, fecha))
    seccion.append("")
    seccion.append("- fecha            : %s" % fecha)
    seccion.append("- plataforma       : %s" % platform.platform())
    seccion.append("- lenguaje         : Python %s" % platform.python_version())
    seccion.append("- paquetes         : "
                   + "; ".join("%s %s" % (k, v) for k, v in ver_pkgs.items()))
    seccion.append("- semilla          : %d" % cfg.SEMILLA)
    seccion.append("- fuente de datos  : %s" % fuente)
    seccion.append("- CSV R<->Python   : %d comparados; byte-identicos %d; "
                   "fuera de tol %d; dif texto %d; peor |dif| %.2e"
                   % (len(comps), n_byte, n_fuera, n_dif_txt, peor_abs))
    seccion.append("- paridad render   : %s" % paridad_render)
    seccion.append("- informe.pdf      : %s (estado 12_informe: %s)"
                   % ("presente" if pdf_ok else "ausente", pdf_status))
    seccion.append("- verificaciones   : %d/%d chequeos duros OK; avisos: %d"
                   % (chk_ok, chk_tot, len(warn)))
    seccion.append("- resultado        : %s"
                   % ("TODAS LAS VERIFICACIONES PASARON" if not dur
                      else "VERIFICACION FALLIDA (%d)" % len(dur)))
    if dur:
        seccion.append("")
        seccion.append("### Fallos")
        seccion += ["  - " + d for d in dur]
    if warn:
        seccion.append("")
        seccion.append("### Avisos")
        seccion += ["  - " + w for w in warn]
    seccion.append("")
    _escribir_log(fecha, "\n".join(seccion))

    # --- salida legible ---------------------------------------------
    print("== 99_verificar.py ==")
    print("  fuente = %s" % fuente)
    print("  checklist: %d/%d chequeos duros OK" % (chk_ok, chk_tot))
    print("  CSV R<->Python: %d; byte-identicos %d; fuera de tol %d; peor |dif| %.2e"
          % (len(comps), n_byte, n_fuera, peor_abs))
    print("  paridad render: %s" % paridad_render)
    print("  informe.pdf: %s (%s)" % ("ok" if pdf_ok else "ausente", pdf_status))
    print("  log -> logs/corrida_%s.txt" % fecha)
    for w in warn:
        print("  aviso: %s" % w)
    if dur:
        print("")
        for d in dur:
            print("  FALLA: %s" % d)
        print("")
        print("  *** VERIFICACION FALLIDA: %d chequeos duros no pasaron ***" % len(dur))
        sys.exit(1)
    print("")
    print("  TODAS LAS VERIFICACIONES PASARON")


def _escribir_log(fecha: str, seccion_lang: str):
    ruta = cfg.RUTA_LOGS / ("corrida_%s.txt" % fecha)
    cab = "# Corrida de verificacion -- reanalisis MIA-LPS\n"
    marca_ini = "<!-- %s:inicio -->" % LANG
    marca_fin = "<!-- %s:fin -->" % LANG
    bloque = "%s\n%s\n%s" % (marca_ini, seccion_lang, marca_fin)
    if ruta.is_file():
        txt = ruta.read_text(encoding="utf-8")
    else:
        txt = cab + "\n"
    if marca_ini in txt and marca_fin in txt:
        txt = txt.split(marca_ini)[0] + bloque + txt.split(marca_fin)[1]
    else:
        if not txt.endswith("\n"):
            txt += "\n"
        txt += "\n" + bloque + "\n"
    ruta.write_bytes(txt.encode("utf-8"))


if __name__ == "__main__":
    main()
