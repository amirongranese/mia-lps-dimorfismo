# 02_ingesta_qc.py -- Ingesta y control de calidad de los tres archivos de entrada.
#
# Por que existe este archivo: es la unica puerta de entrada de los datos al
# pipeline. Lee los tres Excel (crudo real si esta en data/raw/, si no el
# sintetico versionado), los deja en formato largo y limpio, y aplica -- sin
# tocar ninguna otra decision -- las tres reglas de saneamiento de la Seccion 2
# del brief:
#
#   * qPCR: CT_CRUDO == 40 y celda vacia son la MISMA cosa -> "no detectado" ->
#     NA (D3). Nunca se imputa. El tejido BRAIN_P1 se excluye del alcance E15
#     (se analiza en otro informe) y queda registrado en analisis_descartados.md.
#   * ELISA: la hoja mezcla suero materno y liquido amniotico; se parte en dos
#     bloques por la columna TEJIDO. Conc < 0 es censura a izquierda (D10): se
#     marca censurado=TRUE, el valor se guarda como NA y el LOD (= 0, el blanco)
#     se registra aparte. La columna "IL-6" del crudo NO se usa (ver
#     analisis_descartados.md).
#   * pSTAT3: se normaliza SEXO/TTO y se verifica el balanceo de las 3 membranas.
#
# Entregable: el "reporte de QC" -- tablas de n real por archivo x grupo x sexo,
# de no-detectados por gen x tejido x grupo, de censura del ELISA por bloque x
# grupo, y de que madres/fetos faltan en cada bloque respecto del diseno qPCR.
# Este script NO cuantifica, NO modela y NO grafica.

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import openpyxl

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)

# Identificador para los registros compartidos (procedencia, verificaciones,
# analisis_descartados): SIN extension, para que R y Python -- que producen los
# mismos artefactos -- reemplacen la misma fila y el archivo quede byte-identico.
ESTE_SCRIPT = "02_ingesta_qc"

# ---------------------------------------------------------------------------
# Ordenes canonicos: todas las salidas se ordenan por estas claves para que la
# comparacion byte a byte con la implementacion R sea posible.
# ---------------------------------------------------------------------------
ORDEN_TEJIDO = {t: i for i, t in enumerate(cfg.TEJIDOS_E15)}
ORDEN_GEN = {g: i for i, g in enumerate(cfg.GENES)}
ORDEN_TTO = {"CONTROL": 0, "LPS": 1}
ORDEN_SEXO = {"": 0, "HEMBRA": 1, "MACHO": 2}
GRUPOS_4 = ["HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS"]
ORDEN_GRUPO = {g: i for i, g in enumerate(GRUPOS_4)}


def grupo_norm(sexo: str, tto: str) -> str:
    return f"{sexo}_{tto}"


# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a los de 01_generar_sinteticos (paridad
# R/Python): NA -> "", float -> %.10g, entero -> sin decimales.
# ---------------------------------------------------------------------------
def _fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, bool):
        return "TRUE" if x else "FALSE"
    if isinstance(x, float):
        if math.isnan(x):
            return ""
        return "%.10g" % x
    if isinstance(x, int):
        return str(x)
    return str(x)


def _csv_cell(x) -> str:
    s = _fmt(x)
    if any(c in s for c in (",", '"', "\n", "\r")):
        s = '"' + s.replace('"', '""') + '"'
    return s


def escribir_tsv(ruta: Path, encabezado: list[str], filas: list[list]) -> None:
    lineas = ["\t".join(encabezado)]
    lineas += ["\t".join(_fmt(v) for v in fila) for fila in filas]
    ruta.write_bytes(("\n".join(lineas) + "\n").encode("utf-8"))


def escribir_csv(ruta: Path, encabezado: list[str], filas: list[list]) -> None:
    lineas = [",".join(_csv_cell(v) for v in encabezado)]
    lineas += [",".join(_csv_cell(v) for v in fila) for fila in filas]
    ruta.write_bytes(("\n".join(lineas) + "\n").encode("utf-8"))


def escribir_texto(ruta: Path, texto: str) -> None:
    if not texto.endswith("\n"):
        texto += "\n"
    ruta.write_bytes(texto.encode("utf-8"))


# ---------------------------------------------------------------------------
def to_num(v):
    """Cualquier celda -> float o None. No lanza: una celda no numerica es None
    (se contabiliza aparte). Acepta coma decimal por si el crudo la trae."""
    if v is None or isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        f = float(v)
        return None if math.isnan(f) else f
    s = str(v).strip()
    if s == "":
        return None
    for cand in (s, s.replace(",", ".")):
        try:
            return float(cand)
        except ValueError:
            pass
    return None


def s_txt(v) -> str:
    return "" if v is None else str(v).strip()


def madre_id_de(feto: str) -> str:
    """Parte-madre de un ID de feto: quita el sufijo '.NN' final."""
    import re

    return re.sub(r"\.[0-9]+$", "", feto)


def _leer_grid_xlsx(ruta: Path, hoja: str) -> list[tuple]:
    wb = openpyxl.load_workbook(ruta, data_only=True, read_only=True)
    ws = wb[hoja]
    filas = [tuple(r) for r in ws.iter_rows(values_only=True)]
    wb.close()
    return filas


def _leer_grid_tsv(ruta: Path) -> list[tuple]:
    txt = ruta.read_bytes().decode("utf-8")
    lineas = txt.split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    return [tuple(None if c == "" else c for c in ln.split("\t")) for ln in lineas]


def leer_grid(nombre_archivo: str, hoja: str) -> list[tuple]:
    """Grilla cruda (lista de tuplas por fila, sin encabezado asumido).

    - Datos reales: se lee el .xlsx de data/raw/.
    - Datos sinteticos: se lee el mirror .tsv canonico de data/synthetic/ (forma
      versionada, byte-identica R/Python; el .xlsx sintetico no se versiona). Si
      no hubiera .tsv pero si el .xlsx, se cae a este ultimo.
    """
    ruta = cfg.ruta_datos(nombre_archivo)
    if cfg.fuente_datos(nombre_archivo) == "real":
        return _leer_grid_xlsx(ruta, hoja)
    tallo = nombre_archivo.rsplit(".", 1)[0]
    for tsv in (cfg.RUTA_DATOS_SINT / f"{tallo}.tsv",
                cfg.RUTA_DATOS_SINT / f"{tallo}__{hoja}.tsv"):
        if tsv.is_file():
            return _leer_grid_tsv(tsv)
    if ruta.is_file():
        return _leer_grid_xlsx(ruta, hoja)
    raise SystemExit(f"ingesta: no hay fuente para {nombre_archivo!r} (hoja {hoja!r})")


def _ordkey(*vals):
    return tuple(vals)


# ===========================================================================
# 1. qPCR: lectura, 40/vacio -> NA, exclusion de BRAIN_P1
#
# Por que aca: el resto del pipeline (cuantificacion, modelos) parte de esta
# tabla larga. Se conserva CT_CRUDO tal cual venia (para auditar) y se agrega
# no_detectado + CT limpio. BRAIN_P1 se separa: no entra al analisis E15.
# ===========================================================================
COLS_QPCR = ["MADRE", "NOMINACION", "FETO", "GRUPO", "SEXO", "TTO", "TEJIDO",
             "rsp29", "GEN", "CT_CRUDO"]


def ingesta_qpcr():
    grid = leer_grid(cfg.ARCHIVO_QPCR, cfg.HOJA_QPCR)
    encab = [s_txt(c) for c in grid[0]]
    if encab[:10] != COLS_QPCR:
        raise SystemExit(f"qPCR: encabezado inesperado {encab[:10]!r}")

    e15, bp1 = [], []
    n_ct40 = n_ctvacio = n_ct_no_num = 0
    for r in grid[1:]:
        if all(c is None or s_txt(c) == "" for c in r):
            continue
        tejido = s_txt(r[6])
        gen = s_txt(r[8])
        crudo = r[9]
        x = to_num(crudo)
        es_vacio = (crudo is None) or (s_txt(crudo) == "")
        es_40 = (x is not None) and (x == 40.0)
        if (not es_vacio) and x is None:
            n_ct_no_num += 1
        no_det = es_vacio or es_40
        if es_40:
            n_ct40 += 1
        if es_vacio:
            n_ctvacio += 1

        sexo = cfg.normalizar_sexo(r[4])
        tto = cfg.normalizar_tto(r[5])
        feto = s_txt(r[2])
        fila = {
            "MADRE": s_txt(r[0]), "MADRE_ID": madre_id_de(feto) if feto else "",
            "NOMINACION": s_txt(r[1]), "FETO": feto,
            "GRUPO_RAW": s_txt(r[3]), "SEXO_RAW": s_txt(r[4]), "TTO_RAW": s_txt(r[5]),
            "SEXO": sexo, "TTO": tto, "GRUPO": grupo_norm(sexo, tto),
            "TEJIDO": tejido, "rsp29": to_num(r[7]), "GEN": gen,
            "CT_CRUDO": x, "no_detectado": "TRUE" if no_det else "FALSE",
            "CT": None if no_det else x,
        }
        if tejido == cfg.TEJIDO_EXCLUIDO:
            bp1.append(fila)
        else:
            e15.append(fila)

    # --- invariantes de la Seccion 2.1 -----------------------------------
    rsp = [f["rsp29"] for f in e15 + bp1]
    assert all(v is not None and v != 40.0 for v in rsp), "rsp29 con NA o == 40"
    tej_e15 = sorted({f["TEJIDO"] for f in e15})
    assert tej_e15 == sorted(cfg.TEJIDOS_E15), tej_e15
    fetos_e15 = {f["FETO"] for f in e15}
    assert len(fetos_e15) == 36, f"E15 tiene {len(fetos_e15)} fetos, no 36"
    assert all(f["FETO"] for f in e15), "hay FETO vacio en E15"
    assert len(e15) == 36 * len(cfg.TEJIDOS_E15) * len(cfg.GENES), len(e15)

    diag = dict(n_ct40=n_ct40, n_ctvacio=n_ctvacio, n_ct_no_num=n_ct_no_num,
                n_bp1_filas=len(bp1),
                madres_bp1=sorted({f["MADRE"] for f in bp1}))
    return e15, bp1, diag


# ===========================================================================
# 2. ELISA: parte suero materno / liquido amniotico, censura a izquierda (D10)
#
# Por que aca: el crudo pone las dos matrices en una sola hoja con una fila de
# titulo y encabezados con espacios finales. Se localiza el encabezado por la
# celda "TEJIDO" y se filtran solo las filas cuya columna TEJIDO es una de las
# dos etiquetas conocidas (asi se ignoran las celdas de scratch a la derecha).
# ===========================================================================
LOD_ELISA = 0.0  # el blanco de la placa; alternativa (menor estandar de CURVA IL6) en D10
TEJIDO_MS = "Suero materno"
TEJIDO_LA = "Líquido amniótico"


def ingesta_elisa():
    grid = leer_grid(cfg.ARCHIVO_ELISA, cfg.HOJA_ELISA)
    fila_encab = None
    for i, r in enumerate(grid):
        if len(r) > 5 and s_txt(r[5]) == "TEJIDO":
            fila_encab = i
            break
    if fila_encab is None:
        raise SystemExit("ELISA: no se encontro la fila de encabezado (celda 'TEJIDO')")
    encab = [s_txt(c) for c in grid[fila_encab]]
    esperado = ["TRATAMIENTO", "Abs 450", "Conc", "IL-6", "IDMADRE", "TEJIDO"]
    for k, pref in enumerate(esperado):
        if not encab[k].startswith(pref):
            raise SystemExit(f"ELISA: columna {k} = {encab[k]!r}, se esperaba ~{pref!r}")

    filas = []
    n = dict(MS=0, LA=0, MS_cens=0, LA_cens=0)
    for r in grid[fila_encab + 1:]:
        tejido = s_txt(r[5])
        if tejido not in (TEJIDO_MS, TEJIDO_LA):
            continue
        bloque = "MS" if tejido == TEJIDO_MS else "LA"
        trat_raw = s_txt(r[0])
        tto = cfg.normalizar_tto(trat_raw.split()[0])
        if bloque == "LA":
            sexo = ("HEMBRA" if cfg.EMOJI_HEMBRA in trat_raw
                    else "MACHO" if cfg.EMOJI_MACHO in trat_raw else "")
            if sexo == "":
                raise SystemExit(f"ELISA-LA: no se pudo leer el sexo de {trat_raw!r}")
        else:
            sexo = ""  # suero materno: no hay sexo fetal
        idm = s_txt(r[4])
        conc = to_num(r[2])
        if conc is None:
            raise SystemExit(f"ELISA: Conc no numerica en fila {r!r}")
        censurado = conc < 0.0
        n[bloque] += 1
        if censurado:
            n[bloque + "_cens"] += 1
        filas.append({
            "bloque": bloque, "TRATAMIENTO_RAW": trat_raw,
            "TTO": tto, "SEXO": sexo,
            "ID": idm, "MADRE_ID": idm if bloque == "MS" else madre_id_de(idm),
            "TEJIDO": tejido,
            "Abs450": to_num(r[1]), "Conc": conc, "IL6_col_crudo": to_num(r[3]),
            "censurado": "TRUE" if censurado else "FALSE",
            "IL6_pgml": None if censurado else conc, "LOD": LOD_ELISA,
        })

    assert n["MS"] > 0 and n["LA"] > 0, n
    return filas, n


# ===========================================================================
# 3. pSTAT3: normaliza SEXO/TTO, verifica balanceo de membranas (bloque D9)
# ===========================================================================
COLS_PSTAT3 = ["MEMBRANA", "MADRE", "NOMINACION", "FETO", "GRUPO", "SEXO", "TTO",
               "TEJIDO", "PSTAT3"]


def ingesta_pstat3():
    grid = leer_grid(cfg.ARCHIVO_PSTAT3, cfg.HOJA_PSTAT3)
    encab = [s_txt(c) for c in grid[0]]
    if encab[:9] != COLS_PSTAT3:
        raise SystemExit(f"pSTAT3: encabezado inesperado {encab[:9]!r}")
    filas = []
    for r in grid[1:]:
        if all(c is None or s_txt(c) == "" for c in r):
            continue
        sexo = cfg.normalizar_sexo(r[5])
        tto = cfg.normalizar_tto(r[6])
        feto = s_txt(r[3])
        memb = to_num(r[0])
        pstat = to_num(r[8])
        assert memb is not None and float(memb).is_integer(), r
        assert pstat is not None, f"PSTAT3 NA en {r!r}"
        filas.append({
            "MEMBRANA": int(memb), "MADRE": s_txt(r[1]),
            "MADRE_ID": madre_id_de(feto), "NOMINACION": s_txt(r[2]), "FETO": feto,
            "GRUPO_RAW": s_txt(r[4]), "SEXO_RAW": s_txt(r[5]), "TTO_RAW": s_txt(r[6]),
            "SEXO": sexo, "TTO": tto, "GRUPO": grupo_norm(sexo, tto),
            "TEJIDO": s_txt(r[7]), "PSTAT3": pstat,
        })
    membs = sorted({f["MEMBRANA"] for f in filas})
    conteo = {m: sum(1 for f in filas if f["MEMBRANA"] == m) for m in membs}
    balinstr = "|".join(f"{m}:{conteo[m]}" for m in membs)
    balanceado = len(set(conteo.values())) == 1 and len(filas) % 4 == 0
    return filas, dict(membs=membs, conteo=conteo, balinstr=balinstr,
                       balanceado=balanceado, n=len(filas))


# ===========================================================================
# 4. Tabla de n real por archivo x grupo x sexo
# ===========================================================================
def tabla_n(qpcr_e15, elisa, pstat3):
    filas = []
    # qPCR: unidad = feto, por tejido x grupo
    for tej in cfg.TEJIDOS_E15:
        sub = [f for f in qpcr_e15 if f["TEJIDO"] == tej]
        for g in GRUPOS_4:
            sg = [f for f in sub if f["GRUPO"] == g]
            sx, tt = g.split("_")
            filas.append([cfg.ARCHIVO_QPCR, tej, "feto", sx, tt, g,
                          len({f["FETO"] for f in sg}), len(sg)])
    # pSTAT3: unidad = feto, por grupo
    for g in GRUPOS_4:
        sg = [f for f in pstat3 if f["GRUPO"] == g]
        sx, tt = g.split("_")
        filas.append([cfg.ARCHIVO_PSTAT3, "", "feto", sx, tt, g,
                      len({f["FETO"] for f in sg}), len(sg)])
    # ELISA-MS: unidad = madre, por tratamiento (sin sexo fetal)
    ms = [f for f in elisa if f["bloque"] == "MS"]
    for tt in ("CONTROL", "LPS"):
        sg = [f for f in ms if f["TTO"] == tt]
        filas.append([cfg.ARCHIVO_ELISA, "MS (Suero materno)", "madre", "", tt, "",
                      len({f["ID"] for f in sg}), len(sg)])
    # ELISA-LA: unidad = saco, por tratamiento x sexo
    la = [f for f in elisa if f["bloque"] == "LA"]
    for tt in ("CONTROL", "LPS"):
        for sx in ("HEMBRA", "MACHO"):
            sg = [f for f in la if f["TTO"] == tt and f["SEXO"] == sx]
            filas.append([cfg.ARCHIVO_ELISA, "LA (Liquido amniotico)", "saco", sx, tt,
                          grupo_norm(sx, tt), len({f["ID"] for f in sg}), len(sg)])

    orden_arch = {cfg.ARCHIVO_QPCR: 0, cfg.ARCHIVO_PSTAT3: 1, cfg.ARCHIVO_ELISA: 2}
    orden_bloque = {"PLACENTA_E15": 0, "BRAIN_E15": 1, "": 2,
                    "MS (Suero materno)": 3, "LA (Liquido amniotico)": 4}
    filas.sort(key=lambda x: (orden_arch[x[0]], orden_bloque[x[1]], ORDEN_TTO[x[4]],
                              ORDEN_SEXO[x[3]]))
    return ["archivo", "bloque", "unidad", "SEXO", "TTO", "GRUPO",
            "n_unidades", "n_registros"], filas


# ===========================================================================
# 5. No-detectados de qPCR por gen x tejido x grupo (consecuencia de 40 -> NA)
#
# Por que importa: ademas del QC, esta tabla es la que hace visible el caso D7
# (un gen x tejido con 0 detectados en el calibrador HEMBRA_CONTROL no es
# cuantificable). El script solo lo REPORTA; la regla se aplica en 05.
# ===========================================================================
def tabla_no_detectados(qpcr_e15):
    filas = []
    cero_calib = []
    for tej in cfg.TEJIDOS_E15:
        for gen in cfg.GENES:
            for g in GRUPOS_4:
                sub = [f for f in qpcr_e15
                       if f["TEJIDO"] == tej and f["GEN"] == gen and f["GRUPO"] == g]
                nt = len(sub)
                nd = sum(1 for f in sub if f["no_detectado"] == "FALSE")
                filas.append([tej, gen, g, nt, nd, nt - nd,
                              (nd / nt) if nt else None])
                if g == "HEMBRA_CONTROL" and nt > 0 and nd == 0:
                    cero_calib.append(f"{gen}@{tej}")
    filas.sort(key=lambda x: (ORDEN_TEJIDO[x[0]], ORDEN_GEN[x[1]], ORDEN_GRUPO[x[2]]))
    header = ["TEJIDO", "GEN", "GRUPO", "n_total", "n_detectado", "n_no_detectado",
              "prop_detectado"]
    return header, filas, cero_calib


# ===========================================================================
# 6. Censura del ELISA por bloque x grupo (D10: % de censura ANTES de estadistica)
# ===========================================================================
def tabla_censura(elisa):
    filas = []
    ms = [f for f in elisa if f["bloque"] == "MS"]
    for tt in ("CONTROL", "LPS"):
        sg = [f for f in ms if f["TTO"] == tt]
        nc = sum(1 for f in sg if f["censurado"] == "TRUE")
        filas.append(["MS", tt, "", len(sg), nc, len(sg) - nc,
                      (100.0 * nc / len(sg)) if sg else None])
    la = [f for f in elisa if f["bloque"] == "LA"]
    for tt in ("CONTROL", "LPS"):
        for sx in ("HEMBRA", "MACHO"):
            sg = [f for f in la if f["TTO"] == tt and f["SEXO"] == sx]
            nc = sum(1 for f in sg if f["censurado"] == "TRUE")
            filas.append(["LA", tt, sx, len(sg), nc, len(sg) - nc,
                          (100.0 * nc / len(sg)) if sg else None])
    filas.sort(key=lambda x: (0 if x[0] == "MS" else 1, ORDEN_TTO[x[1]],
                              ORDEN_SEXO[x[2]]))
    return ["bloque", "TTO", "SEXO", "n", "n_censurado", "n_detectado",
            "pct_censurado"], filas


# ===========================================================================
# 7. Que madres/fetos faltan en cada bloque respecto del diseno qPCR (Seccion 2.2)
#
# Por que asi: el ELISA se colecto con su propio n; no se completa ni se imputa.
# Se toma como universo de referencia los MADRE_ID y FETO de qPCR E15 y se
# reporta, en ambas direcciones, que IDs de la referencia no tienen dato y que
# IDs del dato no matchean la referencia.
# ===========================================================================
def tabla_faltantes(qpcr_e15, elisa, pstat3):
    ref_madres = sorted({f["MADRE_ID"] for f in qpcr_e15})
    ref_fetos = sorted({f["FETO"] for f in qpcr_e15})
    ms_ids = sorted({f["ID"] for f in elisa if f["bloque"] == "MS"})
    la_ids = sorted({f["ID"] for f in elisa if f["bloque"] == "LA"})
    la_madres = sorted({f["MADRE_ID"] for f in elisa if f["bloque"] == "LA"})
    ps_fetos = sorted({f["FETO"] for f in pstat3})

    def fila(bloque, referencia, ref, dato):
        rs, ds = set(ref), set(dato)
        return [bloque, referencia, len(rs), len(ds), len(rs & ds),
                ";".join(sorted(rs - ds)), ";".join(sorted(ds - rs))]

    filas = [
        fila("ELISA_MS", "qpcr_E15_madre_id", ref_madres, ms_ids),
        fila("ELISA_LA", "qpcr_E15_feto", ref_fetos, la_ids),
        fila("ELISA_LA", "qpcr_E15_madre_id", ref_madres, la_madres),
        fila("pSTAT3", "qpcr_E15_feto", ref_fetos, ps_fetos),
    ]
    return ["bloque", "referencia", "n_ref", "n_dato", "n_match",
            "ids_ref_sin_dato", "ids_dato_sin_ref"], filas


# ===========================================================================
# 8. Resumen maquina-legible + reporte legible
# ===========================================================================
def construir_resumen(qd, ed, pd_, cero_calib, cens):
    def row(chk, val, esp, ok):
        return [chk, _fmt(val), _fmt(esp), "TRUE" if ok else "FALSE"]

    n_ms = ed["MS"]
    n_la = ed["LA"]
    filas = [
        row("fuente_qpcr", cfg.fuente_datos(cfg.ARCHIVO_QPCR), "", True),
        row("fuente_elisa", cfg.fuente_datos(cfg.ARCHIVO_ELISA), "", True),
        row("fuente_pstat3", cfg.fuente_datos(cfg.ARCHIVO_PSTAT3), "", True),
        row("qpcr_fetos_e15", qd["n_fetos_e15"], 36, qd["n_fetos_e15"] == 36),
        row("qpcr_filas_e15", qd["n_filas_e15"], 720, qd["n_filas_e15"] == 720),
        row("qpcr_tejidos_e15", ";".join(cfg.TEJIDOS_E15), ";".join(cfg.TEJIDOS_E15), True),
        row("qpcr_brain_p1_filas_excluidas", qd["n_bp1_filas"], ">0", qd["n_bp1_filas"] > 0),
        row("qpcr_ct_igual_40", qd["n_ct40"], ">=0", True),
        row("qpcr_ct_vacio", qd["n_ctvacio"], ">=0", True),
        row("qpcr_ct_no_detectado", qd["n_ct40"] + qd["n_ctvacio"], ">=0", True),
        row("qpcr_ct_no_numerico", qd["n_ct_no_num"], 0, qd["n_ct_no_num"] == 0),
        row("qpcr_genes_cero_det_en_calibrador",
            ";".join(cero_calib) if cero_calib else "(ninguno)", "", True),
        row("elisa_ms_n", n_ms, 14, n_ms == 14),
        row("elisa_ms_censurado", ed["MS_cens"], ">=0", True),
        row("elisa_la_n", n_la, 28, n_la == 28),
        row("elisa_la_censurado", ed["LA_cens"], ">=0", True),
        row("elisa_columna_il6_usada", "no (ver analisis_descartados.md)", "no", True),
        row("pstat3_filas", pd_["n"], 36, pd_["n"] == 36),
        row("pstat3_membranas", pd_["balinstr"], "1:12|2:12|3:12", pd_["balanceado"]),
        row("pstat3_balanceado", pd_["balanceado"], True, pd_["balanceado"]),
    ]
    header = ["check", "valor", "esperado", "ok"]
    return header, filas


def construir_reporte_md(ntab, ndet_head, ndet, cens_head, cens, falt_head, falt,
                         res_head, res, qd, ed, pd_):
    L = []
    ap = L.append
    ap("# Reporte de QC -- ingesta MIA-LPS")
    ap("")
    ap("Generado por `02_ingesta_qc` (R y Python producen este archivo identico).")
    ap("Fuentes en uso: "
       f"qPCR=`{cfg.fuente_datos(cfg.ARCHIVO_QPCR)}`, "
       f"ELISA=`{cfg.fuente_datos(cfg.ARCHIVO_ELISA)}`, "
       f"pSTAT3=`{cfg.fuente_datos(cfg.ARCHIVO_PSTAT3)}`.")
    ap("")
    ap("## 1. Saneamiento aplicado")
    ap("")
    ap(f"- qPCR `CT_CRUDO`: {qd['n_ct40']} celdas `== 40` + {qd['n_ctvacio']} vacias "
       f"= {qd['n_ct40'] + qd['n_ctvacio']} no-detectados -> `NA` (D3, sin imputar).")
    ap(f"- qPCR `BRAIN_P1`: {qd['n_bp1_filas']} filas excluidas del alcance E15 "
       f"(madres {', '.join(qd['madres_bp1'])}); ver `analisis_descartados.md`.")
    ap(f"- ELISA: hoja partida en MS (n={ed['MS']}) y LA (n={ed['LA']}); "
       f"`Conc < 0` -> censura a izquierda (D10), valor `NA`, `LOD = {LOD_ELISA:g}`.")
    ap("- ELISA columna `IL-6` del crudo: NO usada (ver `analisis_descartados.md`).")
    ap(f"- pSTAT3: {pd_['n']} filas, membranas `{pd_['balinstr']}` "
       f"({'balanceado' if pd_['balanceado'] else 'DESBALANCEADO'}).")
    ap("")
    ap("## 2. n real por archivo x grupo x sexo")
    ap("")
    ap(_md_tabla(ntab[0], ntab[1]))
    ap("")
    ap("## 3. Censura del ELISA por bloque x grupo (D10)")
    ap("")
    ap(_md_tabla(cens_head, cens))
    ap("")
    ap("## 4. No-detectados de qPCR: casos con 0 detectados en el calibrador")
    ap("")
    cero = [f for f in ndet if f[2] == "HEMBRA_CONTROL" and f[3] > 0 and f[4] == 0]
    if cero:
        ap("Genes x tejido no cuantificables por D7 (se excluyen del modelo en 05):")
        ap("")
        ap(_md_tabla(ndet_head, cero))
    else:
        ap("Ninguno: todos los gen x tejido tienen >=1 deteccion en HEMBRA_CONTROL.")
    ap("")
    ap("(Tabla completa gen x tejido x grupo en `qc_no_detectados_qpcr.csv`.)")
    ap("")
    ap("## 5. Madres / fetos faltantes por bloque (respecto del diseno qPCR E15)")
    ap("")
    ap(_md_tabla(falt_head, [f[:5] + [_trunc(f[5]), _trunc(f[6])] for f in falt]))
    ap("")
    ap("(Listas completas de IDs en `qc_faltantes.csv`.)")
    ap("")
    ap("## 6. Chequeos")
    ap("")
    ap(_md_tabla(res_head, res))
    ap("")
    return "\n".join(L)


def _trunc(s: str, n: int = 60) -> str:
    return s if len(s) <= n else s[: n - 1] + "…"


def _md_tabla(header, filas) -> str:
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for f in filas:
        out.append("| " + " | ".join(_fmt(v) for v in f) + " |")
    return "\n".join(out)


# ===========================================================================
# 9. Artefactos compartidos: analisis_descartados / procedencia / verificaciones
#
# Por que con "merge por script": estos tres archivos los van escribiendo todas
# las tareas. Cada corrida reemplaza SOLO sus propias filas (identificadas por la
# columna 'script') y reescribe el archivo ordenado -> R y Python convergen al
# mismo contenido.
# ===========================================================================
def _leer_csv(ruta: Path):
    if not ruta.is_file():
        return None, []
    txt = ruta.read_bytes().decode("utf-8")
    lineas = txt.split("\n")
    if lineas and lineas[-1] == "":
        lineas.pop()
    if not lineas:
        return None, []
    return _parse_csv_line(lineas[0]), [_parse_csv_line(x) for x in lineas[1:]]


def _parse_csv_line(linea: str):
    out, cur, i, q = [], [], 0, False
    while i < len(linea):
        c = linea[i]
        if q:
            if c == '"':
                if i + 1 < len(linea) and linea[i + 1] == '"':
                    cur.append('"')
                    i += 1
                else:
                    q = False
            else:
                cur.append(c)
        else:
            if c == '"':
                q = True
            elif c == ",":
                out.append("".join(cur))
                cur = []
            else:
                cur.append(c)
        i += 1
    out.append("".join(cur))
    return out


def merge_por_script(ruta: Path, header: list[str], filas_nuevas: list[list],
                     col_script: str, clave_orden):
    h_old, filas_old = _leer_csv(ruta)
    if h_old is not None and h_old != header:
        raise SystemExit(f"{ruta.name}: encabezado incompatible {h_old!r} vs {header!r}")
    idx = header.index(col_script)
    conservadas = [f for f in filas_old if len(f) > idx and f[idx] != ESTE_SCRIPT]
    todas = conservadas + [[_fmt(v) for v in f] for f in filas_nuevas]
    todas.sort(key=clave_orden)  # orden por codepoint (== method="radix" en R)
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
## 02_ingesta_qc

### BRAIN_P1 excluido de la ingesta

- **Que se probo / que trae el crudo:** el archivo `Raw data CTs.xlsx` incluye
  filas con `TEJIDO == BRAIN_P1` (cerebro, dia postnatal 1) ademas de
  `PLACENTA_E15` y `BRAIN_E15`.
- **Por que no se usa aca:** queda fuera del alcance de la Seccion 3 del brief
  (placenta y cerebro fetal E15). El P1 se analiza en un informe aparte.
- **Que se hizo en su lugar:** `02_ingesta_qc` filtra esas filas antes de
  cualquier calculo, deja constancia de cuantas eran y de que madres provienen
  (ver `qc_resumen.csv`, fila `qpcr_brain_p1_filas_excluidas`), y el resto del
  pipeline opera solo sobre `PLACENTA_E15` + `BRAIN_E15`.

### Columna `IL-6` del ELISA: no se usa

- **Que se probo / que trae el crudo:** la hoja `Sueros y LA` tiene una columna
  `IL-6` (indice 3) ademas de `Conc` (indice 2). En casi todas las filas es
  `max(Conc, 0)`; en el bloque de suero materno / LPS es `Conc x 4` (un factor
  de dilucion no documentado en la hoja).
- **Por que no sirve:** mezcla dos transformaciones (piso en 0 y factor de
  dilucion), pisa la censura a izquierda y no es una medida homogenea entre
  filas.
- **Que se hizo en su lugar:** se ignora la columna `IL-6`. El analisis del
  ELISA (03) usa `Conc` (indice 2) con el tratamiento de censura de D10:
  `Conc < 0` -> `censurado = TRUE`, valor `NA`, `LOD` registrado aparte.

### Imputacion de no-detectados (recordatorio)

- D3 ya fija que los no-detectados de qPCR (`CT_CRUDO == 40` o celda vacia) van a
  `NA` y **no se imputan por ningun metodo**. El detalle de por que se descarto la
  imputacion MNAR (`nondetects`) se documenta con numeros en T4/T10.
"""


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 02_ingesta_qc:inicio -->"
    marca_fin = "<!-- 02_ingesta_qc:fin -->"
    nuevo = f"{marca_ini}\n{BLOQUE_DESCARTES}\n{marca_fin}"
    if ruta.is_file():
        txt = ruta.read_bytes().decode("utf-8")
    else:
        txt = ("# Analisis descartados\n\n"
               "Que se probo, por que no funciono o no se uso, y que se hizo en su"
               " lugar. Una seccion por script.\n")
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
    qpcr_e15, qpcr_bp1, qd = ingesta_qpcr()
    elisa, ed = ingesta_elisa()
    pstat3, pd_ = ingesta_pstat3()

    qd["n_fetos_e15"] = len({f["FETO"] for f in qpcr_e15})
    qd["n_filas_e15"] = len(qpcr_e15)
    pd_["balinstr"] = pd_["balinstr"]

    # --- intermedios largos (data/processed/, regenerables) --------------
    def ordena_qpcr(f):
        return (f["MADRE_ID"], f["FETO"], ORDEN_TEJIDO.get(f["TEJIDO"], 9),
                ORDEN_GEN[f["GEN"]])

    qpcr_e15_s = sorted(qpcr_e15, key=ordena_qpcr)
    cols_e15 = ["MADRE", "MADRE_ID", "NOMINACION", "FETO", "GRUPO_RAW", "SEXO_RAW",
                "TTO_RAW", "SEXO", "TTO", "GRUPO", "TEJIDO", "rsp29", "GEN",
                "CT_CRUDO", "no_detectado", "CT"]
    escribir_tsv(cfg.RUTA_DATOS_PROC / "qpcr_e15_long.tsv", cols_e15,
                 [[f[c] for c in cols_e15] for f in qpcr_e15_s])

    bp1_s = sorted(qpcr_bp1, key=lambda f: (f["MADRE"], f["FETO"], f["GEN"]))
    cols_bp1 = ["MADRE", "NOMINACION", "FETO", "GRUPO_RAW", "SEXO", "TTO", "TEJIDO",
                "rsp29", "GEN", "CT_CRUDO", "no_detectado"]
    escribir_tsv(cfg.RUTA_DATOS_PROC / "qpcr_brain_p1_excluido.tsv", cols_bp1,
                 [[f[c] for c in cols_bp1] for f in bp1_s])

    el_ord = {"MS": 0, "LA": 1}
    elisa_s = sorted(elisa, key=lambda f: (el_ord[f["bloque"]], ORDEN_TTO[f["TTO"]],
                                           ORDEN_SEXO[f["SEXO"]], f["ID"]))
    cols_el = ["bloque", "TRATAMIENTO_RAW", "TTO", "SEXO", "ID", "MADRE_ID", "TEJIDO",
               "Abs450", "Conc", "IL6_col_crudo", "censurado", "IL6_pgml", "LOD"]
    escribir_tsv(cfg.RUTA_DATOS_PROC / "elisa_long.tsv", cols_el,
                 [[f[c] for c in cols_el] for f in elisa_s])

    ps_s = sorted(pstat3, key=lambda f: (f["MEMBRANA"], ORDEN_GRUPO[f["GRUPO"]], f["FETO"]))
    cols_ps = ["MEMBRANA", "MADRE", "MADRE_ID", "NOMINACION", "FETO", "GRUPO_RAW",
               "SEXO_RAW", "TTO_RAW", "SEXO", "TTO", "GRUPO", "TEJIDO", "PSTAT3"]
    escribir_tsv(cfg.RUTA_DATOS_PROC / "pstat3_long.tsv", cols_ps,
                 [[f[c] for c in cols_ps] for f in ps_s])

    # --- tablas de QC (una copia por implementacion) --------------------
    n_head, n_filas = tabla_n(qpcr_e15, elisa, pstat3)
    nd_head, nd_filas, cero_calib = tabla_no_detectados(qpcr_e15)
    c_head, c_filas = tabla_censura(elisa)
    f_head, f_filas = tabla_faltantes(qpcr_e15, elisa, pstat3)
    r_head, r_filas = construir_resumen(qd, ed, pd_, cero_calib, c_filas)

    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "qc_n_por_grupo.csv", n_head, n_filas)
        escribir_csv(base / "qc_no_detectados_qpcr.csv", nd_head, nd_filas)
        escribir_csv(base / "qc_censura_elisa.csv", c_head, c_filas)
        escribir_csv(base / "qc_faltantes.csv", f_head, f_filas)
        escribir_csv(base / "qc_resumen.csv", r_head, r_filas)

    # --- reporte legible + artefactos compartidos ----------------------
    md = construir_reporte_md((n_head, n_filas), nd_head, nd_filas, c_head, c_filas,
                              f_head, f_filas, r_head, r_filas, qd, ed, pd_)
    escribir_texto(cfg.RUTA_TABLAS / "qc_reporte.md", md)
    actualizar_descartados()

    ent_q = f"data/{cfg.fuente_datos(cfg.ARCHIVO_QPCR)}/{cfg.ARCHIVO_QPCR}"
    ent_e = f"data/{cfg.fuente_datos(cfg.ARCHIVO_ELISA)}/{cfg.ARCHIVO_ELISA}"
    ent_p = f"data/{cfg.fuente_datos(cfg.ARCHIVO_PSTAT3)}/{cfg.ARCHIVO_PSTAT3}"
    registrar_procedencia([
        ["data/processed/qpcr_e15_long.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_q, "qPCR largo E15, 40/vacio->NA (D3), BRAIN_P1 excluido"],
        ["data/processed/qpcr_brain_p1_excluido.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_q, "filas BRAIN_P1 apartadas del alcance E15"],
        ["data/processed/elisa_long.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_e, "ELISA largo MS+LA, censura a izquierda D10"],
        ["data/processed/pstat3_long.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_p, "pSTAT3 largo, SEXO/TTO normalizados"],
        ["outputs/tables/{R,python}/qc_n_por_grupo.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_q + "; " + ent_e + "; " + ent_p, "n real por archivo x grupo x sexo"],
        ["outputs/tables/{R,python}/qc_no_detectados_qpcr.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent_q, "no-detectados por gen x tejido x grupo; marca D7"],
        ["outputs/tables/{R,python}/qc_censura_elisa.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent_e, "% censura ELISA por bloque x grupo (D10)"],
        ["outputs/tables/{R,python}/qc_faltantes.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_q + "; " + ent_e + "; " + ent_p,
         "IDs de referencia sin dato y viceversa por bloque"],
        ["outputs/tables/{R,python}/qc_resumen.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_q + "; " + ent_e + "; " + ent_p, "chequeos maquina-legibles de la ingesta"],
        ["outputs/tables/qc_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO",
         ent_q + "; " + ent_e + "; " + ent_p, "reporte de QC legible"],
    ])
    registrar_verificaciones([
        ["ingesta_qpcr_fetos_e15", "qPCR E15 tiene 36 fetos",
         _fmt(qd["n_fetos_e15"]), "36", "TRUE" if qd["n_fetos_e15"] == 36 else "FALSE",
         ESTE_SCRIPT],
        ["ingesta_qpcr_filas_e15", "qPCR E15 tiene 36x2x10 = 720 filas",
         _fmt(qd["n_filas_e15"]), "720", "TRUE" if qd["n_filas_e15"] == 720 else "FALSE",
         ESTE_SCRIPT],
        ["ingesta_qpcr_tejidos", "tejidos E15 = PLACENTA_E15;BRAIN_E15 (sin BRAIN_P1)",
         ";".join(sorted({f["TEJIDO"] for f in qpcr_e15})),
         ";".join(sorted(cfg.TEJIDOS_E15)),
         "TRUE" if sorted({f["TEJIDO"] for f in qpcr_e15}) == sorted(cfg.TEJIDOS_E15)
         else "FALSE", ESTE_SCRIPT],
        ["ingesta_ct_no_numerico", "no hay CT_CRUDO no numerico fuera de vacio",
         _fmt(qd["n_ct_no_num"]), "0", "TRUE" if qd["n_ct_no_num"] == 0 else "FALSE",
         ESTE_SCRIPT],
        ["ingesta_pstat3_balanceado", "pSTAT3: 3 membranas con igual n",
         pd_["balinstr"], "1:12|2:12|3:12", "TRUE" if pd_["balanceado"] else "FALSE",
         ESTE_SCRIPT],
        ["ingesta_elisa_bloques", "ELISA parte en MS y LA con n>0",
         f"MS={ed['MS']};LA={ed['LA']}", "MS>0;LA>0",
         "TRUE" if ed["MS"] > 0 and ed["LA"] > 0 else "FALSE", ESTE_SCRIPT],
    ])

    # --- resumen por consola ------------------------------------------
    print("== 02_ingesta_qc.py ==")
    print(f"  qPCR   fuente={cfg.fuente_datos(cfg.ARCHIVO_QPCR):<9} "
          f"E15: {qd['n_fetos_e15']} fetos / {qd['n_filas_e15']} filas   "
          f"BRAIN_P1 excluido: {qd['n_bp1_filas']} filas")
    print(f"         CT no detectado: {qd['n_ct40']} (==40) + {qd['n_ctvacio']} (vacio)"
          f" = {qd['n_ct40'] + qd['n_ctvacio']}")
    print(f"         no cuantificable por D7: "
          f"{', '.join(cero_calib) if cero_calib else '(ninguno)'}")
    print(f"  ELISA  fuente={cfg.fuente_datos(cfg.ARCHIVO_ELISA):<9} "
          f"MS n={ed['MS']} (cens {ed['MS_cens']})   LA n={ed['LA']} (cens {ed['LA_cens']})")
    print(f"  pSTAT3 fuente={cfg.fuente_datos(cfg.ARCHIVO_PSTAT3):<9} "
          f"{pd_['n']} filas   membranas {pd_['balinstr']}   "
          f"{'balanceado' if pd_['balanceado'] else 'DESBALANCEADO'}")
    print(f"  -> data/processed/*.tsv, outputs/tables/{{R,python}}/qc_*.csv, "
          f"outputs/tables/qc_reporte.md")
    print("  invariantes OK")


if __name__ == "__main__":
    main()
