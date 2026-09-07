# 01_generar_sinteticos.py -- Generador de datos sinteticos.
#
# Por que existe: el objetivo del proyecto es que un tercero, sin acceso a los
# datos crudos (que son ineditos y no se versionan), pueda correr TODO el
# pipeline. Para eso data/synthetic/ contiene tres archivos que imitan la
# estructura EXACTA de los tres crudos -- mismos nombres de hoja y de columna,
# mismos tipos, y las mismas PATOLOGIAS que fuerzan las decisiones del brief:
#
#   * 18 camadas / 36 fetos E15 (una hembra y un macho por madre) + BRAIN_P1
#     como tejido extra que 02_ingesta_qc va a excluir.
#   * CT_CRUDO con no-detectados escritos de dos formas (el valor literal 40 y la
#     celda vacia), igual que en el crudo. Ambas = "no detectado" (D3).
#   * il6 en BRAIN_E15 con 0 detectados en el calibrador (hembra-control): dispara
#     la regla D7 (gen no cuantificable -> solo proporcion de deteccion).
#   * ELISA con la hoja partida en suero materno / liquido amniotico y con
#     valores de Conc negativos = censura a izquierda (D10).
#   * pSTAT3 con 3 membranas balanceadas (bloque fijo de D9).
#
# Los EFECTOS simulados (medias, magnitudes) son arbitrarios: este archivo no
# prueba ninguna hipotesis, solo fabrica un fixture reproducible. La secuencia de
# numeros aleatorios se consume en el mismo orden que en R/01_generar_sinteticos.R
# (mismo RNG propio, misma semilla) -> los .tsv canonicos salen byte a byte
# identicos entre lenguajes (se verifica en data/synthetic/MANIFEST.tsv).

from __future__ import annotations

import hashlib
import importlib.util
import math
from pathlib import Path

import pandas as pd
from openpyxl import Workbook

_cfg_spec = importlib.util.spec_from_file_location(
    "cfg00", Path(__file__).resolve().parent / "00_config.py"
)
cfg = importlib.util.module_from_spec(_cfg_spec)
_cfg_spec.loader.exec_module(cfg)


# ===========================================================================
# 1. Tabla canonica de fetos (deterministica, sin RNG)
#
# Por que sin RNG: los identificadores son estructura, no medicion. Fijarlos sin
# tocar el stream aleatorio deja toda la secuencia de numeros para los valores
# medidos, lo que hace trivial la paridad con R.
# ===========================================================================
_FECHAS = [
    "011025", "021025", "031025", "091025", "101025",
    "171025", "241025", "251025", "311025",
]  # 9 fechas: una por indice de madre (C{i} y L{i} comparten fecha, distinto prefijo)


def _tabla_fetos() -> list[dict]:
    fetos: list[dict] = []
    for prefijo, tto in (("C", "CONTROL"), ("L", "LPS")):
        for i in range(1, 10):  # C1..C9 / L1..L9
            madre = f"{prefijo}{i}"
            madre_id = f"{prefijo}_{_FECHAS[i - 1]}_1"
            for k, sexo in enumerate(cfg.NIVELES_SEXO):  # HEMBRA, luego MACHO
                nn = 3 + 2 * (i - 1) + k  # distinto dentro de la madre, 2 digitos
                emo = cfg.EMOJI_HEMBRA if sexo == "HEMBRA" else cfg.EMOJI_MACHO
                fetos.append(
                    dict(
                        MADRE=madre,
                        MADRE_ID=madre_id,
                        FETO=f"{madre_id}.{nn:02d}",
                        NOMINACION=f"{madre}{emo}{cfg._VS16}",
                        SEXO=sexo,
                        TTO=tto,
                        GRUPO=cfg.etiqueta_grupo(sexo, tto),
                    )
                )
    return fetos


FETOS = _tabla_fetos()  # 36 fetos, orden C1-H, C1-M, C2-H, ... L9-M
FETO_POR_MADRE_SEXO = {(f["MADRE"], f["SEXO"]): f for f in FETOS}

# Madres con tejido BRAIN_P1 en el crudo (subconjunto): C1..C7 y L1..L7.
# Tres de ellas (C6, C7, L7) tienen la celda FETO vacia -> 22 FETO unicos + 6
# vacios = 28 fetos * 10 genes = 280 filas, igual que el crudo.
_MADRES_P1 = [f"C{i}" for i in range(1, 8)] + [f"L{i}" for i in range(1, 8)]
_MADRES_P1_FETO_VACIO = {"C6", "C7", "L7"}

# Genes cuya celda vacia (NA) aparece en el crudo (6 celdas cada uno = 30 NA).
_GENES_CON_NA = ["fatcd36", "glut1", "glut3", "slc38a1", "slc38a2"]
_N_NA_POR_GEN = 6


# ===========================================================================
# 2. qPCR: tabla larga de CT
#
# Por que este modelo: se necesita un Ct plausible por gen x tejido, con los
# genes de la via IL-6 (il6, il6R) altos -> muchos no-detectados, y il6 en
# cerebro practicamente ausente. Los transportadores quedan en rango detectable
# con pocos no-detectados. Un efecto SEXO*TTO minusculo en un par de genes evita
# que todo sea ruido puro; su tamano es arbitrario.
# ===========================================================================
_CT_BASE = {
    "fatcd36": 27.0, "fatp1": 32.0, "fatp4": 28.0, "glut1": 26.0, "glut3": 31.5,
    "slc38a1": 33.0, "slc38a2": 29.5, "gp130": 25.5, "il6R": 36.5, "il6": 37.5,
}
_CT_SHIFT_TEJIDO = {"PLACENTA_E15": 0.0, "BRAIN_E15": 2.0, "BRAIN_P1": 1.5}


def _media_ct(gen: str, tejido: str, sexo: str, tto: str) -> float:
    mu = _CT_BASE[gen] + _CT_SHIFT_TEJIDO[tejido]
    if gen == "il6":
        mu = {"PLACENTA_E15": 35.5, "BRAIN_E15": 40.7, "BRAIN_P1": 41.0}[tejido]
    if gen == "il6R" and tejido != "PLACENTA_E15":
        mu += 2.5
    # efecto arbitrario chico: LPS baja ~0.6 Ct (mas expresion) en machos
    if gen in ("glut1", "slc38a2") and sexo == "MACHO" and tto == "LPS":
        mu -= 0.6
    return mu


def _sigma_ct(gen: str, tejido: str) -> float:
    s = 1.3 if gen in ("il6", "il6R") else 0.8
    if tejido != "PLACENTA_E15":
        s += 0.2
    return s


def _generar_qpcr(rng: cfg.RNG) -> pd.DataFrame:
    filas: list[dict] = []

    def _bloque(feto: dict, tejido: str, feto_field: str) -> None:
        rsp29 = min(max(rng.norm1(25.0, 2.3), 19.5), 31.9)
        for gen in cfg.GENES:  # orden canonico
            lat = _media_ct(gen, tejido, feto["SEXO"], feto["TTO"])
            lat += rng.norm1(0.0, _sigma_ct(gen, tejido))
            filas.append(
                dict(
                    MADRE=feto["MADRE"], NOMINACION=feto["NOMINACION"],
                    FETO=feto_field, GRUPO=feto["GRUPO"], SEXO=feto["SEXO"],
                    TTO=feto["TTO"], TEJIDO=tejido, rsp29=rsp29, GEN=gen,
                    _lat=lat,
                )
            )

    # E15: los 36 fetos, placenta y cerebro
    for feto in FETOS:
        for tejido in cfg.TEJIDOS_E15:
            _bloque(feto, tejido, feto["FETO"])
    # BRAIN_P1: subconjunto de madres; algunas con FETO vacio
    for feto in FETOS:
        if feto["MADRE"] not in _MADRES_P1:
            continue
        ff = "" if feto["MADRE"] in _MADRES_P1_FETO_VACIO else feto["FETO"]
        _bloque(feto, cfg.TEJIDO_EXCLUIDO, ff)

    df = pd.DataFrame(filas)

    # --- Deteccion: latente >= 40 -> no detectado --------------------------
    ct = []
    for lat in df["_lat"]:
        if lat >= 40.0:
            ct.append(40.0)  # no detectado escrito como valor literal 40
        else:
            ct.append(round(min(max(lat, 21.0), 39.9), 2))
    df["CT_CRUDO"] = ct

    # --- D7: il6 / BRAIN_E15 sin detectados en el calibrador --------------
    # En el crudo il6 en cerebro tiene unos pocos detectados fuera del calibrador
    # pero 0/9 en hembra-control. Se deja lo estocastico y se fuerza el calibrador
    # a no-detectado (es la condicion exacta que dispara D7).
    m_cal = (
        (df["GEN"] == "il6")
        & (df["TEJIDO"] == "BRAIN_E15")
        & (df["SEXO"] == "HEMBRA")
        & (df["TTO"] == "CONTROL")
    )
    df.loc[m_cal, "CT_CRUDO"] = 40.0

    # --- Celdas vacias: 6 por gen en 5 genes = 30 NA ---------------------
    # (igual que el crudo: pozos dejados en blanco, con o sin senal)
    for gen in _GENES_CON_NA:  # orden fijo
        cand = list(df.index[df["GEN"] == gen])
        elegidas = rng.elegir_sin_reemplazo(cand, _N_NA_POR_GEN)
        df.loc[elegidas, "CT_CRUDO"] = math.nan

    df = df.drop(columns=["_lat"])
    df = df[
        ["MADRE", "NOMINACION", "FETO", "GRUPO", "SEXO", "TTO", "TEJIDO",
         "rsp29", "GEN", "CT_CRUDO"]
    ]
    return df


# ===========================================================================
# 3. ELISA IL-6: hoja "Sueros y LA" (44 x 14, sin encabezado real) + "CURVA IL6"
#
# Por que esta forma: el crudo es una hoja de calculo con un rotulo suelto en la
# fila 0, los encabezados en la fila 1, dos bloques apilados (suero materno y
# liquido amniotico) y celdas de scratch a la derecha. La ingesta (02) la parte
# por la columna 5. Se reproduce la relacion real Abs<->Conc: los coeficientes
# 0.0016 y 0.1057 son las celdas de scratch M/OO del crudo, y
#   Conc = (Abs - 0.1057) / 0.0016 ,   IL-6 = max(Conc, 0)   [columna a NO usar]
# Los valores negativos de Conc son censura a izquierda (D10).
# ===========================================================================
_ELISA_SLOPE = 0.0016
_ELISA_BLANK_ABS = 0.1057


def _abs_conc_il6(conc_objetivo: float) -> tuple[float, float, float]:
    abs450 = round(conc_objetivo * _ELISA_SLOPE + _ELISA_BLANK_ABS, 3)
    conc = (abs450 - _ELISA_BLANK_ABS) / _ELISA_SLOPE
    return abs450, conc, max(conc, 0.0)


def _generar_elisa(rng: cfg.RNG) -> tuple[list[list], list[list]]:
    ms: list[list] = []  # filas suero materno (14): 5 Control (C1..C5) + 9 LPS (L1..L9)
    for j, madre in enumerate([f"C{i}" for i in range(1, 6)]):
        madre_id = FETO_POR_MADRE_SEXO[(madre, "HEMBRA")]["MADRE_ID"]
        if j < 4:  # 4 de 5 controles por debajo del blanco (censura)
            conc_obj = -(5.0 + abs(rng.norm1(0.0, 15.0)))
        else:
            conc_obj = 30.0 + abs(rng.norm1(0.0, 20.0))
        a, c, i6 = _abs_conc_il6(conc_obj)
        ms.append(["Control", a, c, i6, madre_id, "Suero materno"])
    for madre in [f"L{i}" for i in range(1, 10)]:
        madre_id = FETO_POR_MADRE_SEXO[(madre, "HEMBRA")]["MADRE_ID"]
        conc_obj = 150.0 + abs(rng.norm1(650.0, 380.0))
        a, c, i6 = _abs_conc_il6(conc_obj)
        ms.append(["LPS", a, c, i6, madre_id, "Suero materno"])

    # liquido amniotico (28): Control M/H (C1..C5) + LPS M/H (L1..L9)
    la: list[list] = []

    def _bloque_la(madres: list[str], sexo: str, tto: str, media_baja: float) -> None:
        emo = cfg.EMOJI_MACHO if sexo == "MACHO" else cfg.EMOJI_HEMBRA
        etq = ("Control" if tto == "CONTROL" else "LPS") + " " + emo
        for madre in madres:
            feto = FETO_POR_MADRE_SEXO[(madre, sexo)]
            if tto == "LPS":
                p = rng.unif1()
                if p < 0.25:  # cola alta ocasional (igual que el crudo)
                    conc_obj = 150.0 + abs(rng.norm1(300.0, 250.0))
                else:
                    conc_obj = rng.norm1(media_baja, 22.0)
            else:
                conc_obj = rng.norm1(media_baja, 10.0)
            a, c, i6 = _abs_conc_il6(conc_obj)
            la.append([etq, a, c, i6, feto["FETO"], "Líquido amniótico"])

    c_madres = [f"C{i}" for i in range(1, 6)]
    l_madres = [f"L{i}" for i in range(1, 10)]
    _bloque_la(c_madres, "MACHO", "CONTROL", 7.0)   # ~0 negativos
    _bloque_la(c_madres, "HEMBRA", "CONTROL", -9.0)  # ~3-4/5 negativos
    _bloque_la(l_madres, "MACHO", "LPS", 1.0)
    _bloque_la(l_madres, "HEMBRA", "LPS", 0.0)
    return ms, la


_CURVA_PGML = [0.0, 15.6, 31.25, 62.5, 125.0, 250.0, 500.0, 1000.0]


def _generar_curva(rng: cfg.RNG) -> list[list]:
    filas = []
    for pg in _CURVA_PGML:
        a1 = round(pg * _ELISA_SLOPE + 0.0695 + rng.norm1(0.0, 0.004), 4)
        a2 = round(pg * _ELISA_SLOPE + 0.0695 + rng.norm1(0.0, 0.004), 4)
        # media sin redondear: evita el caso "mitad exacta" en el 4o decimal, que
        # round() resuelve distinto en R y en Python. Es una celda de scratch.
        filas.append([a1, a2, (a1 + a2) / 2, pg])
    return filas


# ===========================================================================
# 4. pSTAT3: 36 filas, 3 membranas balanceadas (bloque fijo de D9)
# ===========================================================================
_MEMBRANAS = [
    ["C1", "C2", "C3", "L1", "L2", "L3"],
    ["C4", "C5", "C6", "L4", "L5", "L6"],
    ["C7", "C8", "C9", "L7", "L8", "L9"],
]


def _generar_pstat3(rng: cfg.RNG) -> pd.DataFrame:
    filas: list[dict] = []
    for m, madres_m in enumerate(_MEMBRANAS, start=1):
        controles = [x for x in madres_m if x.startswith("C")]
        lpss = [x for x in madres_m if x.startswith("L")]
        for sexo in cfg.NIVELES_SEXO:  # HEMBRA, MACHO
            for tto, grp_madres in (("CONTROL", controles), ("LPS", lpss)):
                mu = 2.2 if tto == "CONTROL" else 6.0
                for madre in grp_madres:
                    feto = FETO_POR_MADRE_SEXO[(madre, sexo)]
                    val = min(max(rng.norm1(mu, 0.65), 0.85), 7.2)
                    filas.append(
                        dict(
                            MEMBRANA=m, MADRE=madre, NOMINACION=feto["NOMINACION"],
                            FETO=feto["FETO"], GRUPO=cfg.etiqueta_grupo(sexo, tto),
                            SEXO=("Hembra" if sexo == "HEMBRA" else "Macho"),
                            TTO=("Control" if tto == "CONTROL" else "LPS"),
                            TEJIDO="PLACENTA_E15", PSTAT3=val,
                        )
                    )
    return pd.DataFrame(filas)


# ===========================================================================
# 5. Escritura: .xlsx (para el pipeline) + .tsv canonico (paridad R/Python)
# ===========================================================================
def _fmt(x) -> str:
    if x is None:
        return ""
    if isinstance(x, float):
        if math.isnan(x):
            return ""
        return "%.10g" % x
    if isinstance(x, int):
        return str(x)
    return str(x)


def _escribir_tsv(ruta: Path, filas: list[list], encabezado: list[str] | None) -> str:
    lineas = []
    if encabezado is not None:
        lineas.append("\t".join(encabezado))
    for fila in filas:
        lineas.append("\t".join(_fmt(v) for v in fila))
    texto = "\n".join(lineas) + "\n"
    ruta.write_bytes(texto.encode("utf-8"))
    return hashlib.sha256(texto.encode("utf-8")).hexdigest()


def _df_a_filas(df: pd.DataFrame) -> list[list]:
    out = []
    for row in df.itertuples(index=False, name=None):
        out.append([None if (isinstance(v, float) and math.isnan(v)) else v for v in row])
    return out


def _escribir_xlsx_grid(ruta: Path, hojas: dict[str, list[list]]) -> None:
    wb = Workbook()
    wb.remove(wb.active)
    for nombre, grid in hojas.items():
        ws = wb.create_sheet(title=nombre)
        for fila in grid:
            ws.append([None if v == "" else v for v in fila])
    wb.save(ruta)


def main() -> None:
    rng = cfg.nuevo_rng(cfg.SEMILLA)

    # --- orden de consumo del RNG: qpcr -> elisa -> curva -> pstat3 --------
    qpcr = _generar_qpcr(rng)
    ms, la = _generar_elisa(rng)
    curva = _generar_curva(rng)
    pstat3 = _generar_pstat3(rng)

    dst = cfg.RUTA_DATOS_SINT
    dst.mkdir(parents=True, exist_ok=True)

    # --- qPCR .xlsx + .tsv ---------------------------------------------
    qpcr.to_excel(dst / cfg.ARCHIVO_QPCR, sheet_name=cfg.HOJA_QPCR, index=False)
    h_qpcr = _escribir_tsv(
        dst / "Raw data CTs.tsv", _df_a_filas(qpcr), list(qpcr.columns)
    )

    # --- pSTAT3 .xlsx + .tsv ------------------------------------------
    pstat3.to_excel(dst / cfg.ARCHIVO_PSTAT3, sheet_name=cfg.HOJA_PSTAT3, index=False)
    h_pstat3 = _escribir_tsv(
        dst / "pstat3 placenta.tsv", _df_a_filas(pstat3), list(pstat3.columns)
    )

    # --- ELISA .xlsx (2 hojas tipo grilla, sin encabezado real) ---------
    ancho = 14
    g_elisa: list[list] = []
    g_elisa.append(["IL6"] + [""] * (ancho - 1))
    g_elisa.append(
        ["TRATAMIENTO y sexo", "Abs 450 MS", "Conc ", "IL-6 ", "IDMADRE/FETO",
         "TEJIDO", "", "", "", "", "", "Abs", "m", "oo"]
    )
    fila_scratch = [""] * ancho
    fila_scratch[12] = _ELISA_SLOPE
    fila_scratch[13] = _ELISA_BLANK_ABS
    for k, fila in enumerate(ms + la):
        g = list(fila) + [""] * (ancho - len(fila))
        if k == 0:
            g[12], g[13] = fila_scratch[12], fila_scratch[13]
        g_elisa.append(g)

    alto_curva = 22
    g_curva = [[""] * 10 for _ in range(alto_curva)]
    g_curva[13] = ["", "", "Abs 450", "pg/ml", "", "", "", "", "", ""]
    for r, fila in enumerate(curva, start=14):
        g_curva[r] = list(fila) + [""] * (10 - len(fila))

    # orden de hojas como en el crudo: CURVA IL6 primero
    _escribir_xlsx_grid(
        dst / cfg.ARCHIVO_ELISA, {cfg.HOJA_CURVA: g_curva, cfg.HOJA_ELISA: g_elisa}
    )
    h_elisa = _escribir_tsv(
        dst / "ELISA IL6 2026 Dosis 100__Sueros y LA.tsv", g_elisa, None
    )
    h_curva = _escribir_tsv(
        dst / "ELISA IL6 2026 Dosis 100__CURVA IL6.tsv", g_curva, None
    )

    # --- MANIFEST.tsv: derivado solo del contenido -> byte-identico R/Py --
    manifest = [
        [cfg.ARCHIVO_QPCR, cfg.HOJA_QPCR, len(qpcr), qpcr.shape[1], h_qpcr],
        [cfg.ARCHIVO_PSTAT3, cfg.HOJA_PSTAT3, len(pstat3), pstat3.shape[1], h_pstat3],
        [cfg.ARCHIVO_ELISA, cfg.HOJA_ELISA, len(g_elisa), ancho, h_elisa],
        [cfg.ARCHIVO_ELISA, cfg.HOJA_CURVA, alto_curva, 10, h_curva],
    ]
    _escribir_tsv(
        dst / "MANIFEST.tsv", manifest,
        ["archivo", "hoja", "n_filas", "n_columnas", "sha256"],
    )

    # --- PATOLOGIAS.tsv: contadores que 02+ dan por sentados -------------
    ct40 = int((qpcr["CT_CRUDO"] == 40.0).sum())
    ctna = int(qpcr["CT_CRUDO"].isna().sum())
    il6_be15 = qpcr[(qpcr["GEN"] == "il6") & (qpcr["TEJIDO"] == "BRAIN_E15")]
    il6_be15_det = int(((il6_be15["CT_CRUDO"] < 40) & il6_be15["CT_CRUDO"].notna()).sum())
    il6_be15_cal = il6_be15[il6_be15["GRUPO"] == cfg.etiqueta_grupo("HEMBRA", "CONTROL")]
    il6_be15_cal_det = int(
        ((il6_be15_cal["CT_CRUDO"] < 40) & il6_be15_cal["CT_CRUDO"].notna()).sum()
    )
    ms_neg = sum(1 for f in ms if f[2] < 0)
    la_neg = sum(1 for f in la if f[2] < 0)
    pat = [
        ["qpcr_n_filas", len(qpcr)],
        ["qpcr_n_madres", qpcr["MADRE"].nunique()],
        ["qpcr_n_fetos_e15", qpcr[qpcr["TEJIDO"].isin(cfg.TEJIDOS_E15)]["FETO"].nunique()],
        ["qpcr_n_filas_e15", int(qpcr["TEJIDO"].isin(cfg.TEJIDOS_E15).sum())],
        ["qpcr_ct_igual_40", ct40],
        ["qpcr_ct_vacio", ctna],
        ["qpcr_il6_brain_e15_detectados", il6_be15_det],
        ["qpcr_il6_brain_e15_hembracontrol_detectados", il6_be15_cal_det],
        ["elisa_ms_n", len(ms)],
        ["elisa_ms_conc_neg", ms_neg],
        ["elisa_la_n", len(la)],
        ["elisa_la_conc_neg", la_neg],
        ["pstat3_n_filas", len(pstat3)],
        ["pstat3_membranas", "|".join(
            f"{m}:{int((pstat3['MEMBRANA'] == m).sum())}" for m in (1, 2, 3)
        )],
    ]
    _escribir_tsv(dst / "PATOLOGIAS.tsv", pat, ["patologia", "valor"])

    # --- invariantes duros: si fallan, el fixture no sirve --------------
    assert len(qpcr) == 1000, len(qpcr)
    assert qpcr["MADRE"].nunique() == 18
    assert qpcr[qpcr["TEJIDO"].isin(cfg.TEJIDOS_E15)]["FETO"].nunique() == 36
    assert int(qpcr["TEJIDO"].isin(cfg.TEJIDOS_E15).sum()) == 720
    assert ctna == 30, ctna
    assert ct40 > 0
    assert il6_be15_cal_det == 0, "D7 no se dispara: hay detectados en el calibrador"
    assert il6_be15_det >= 1, il6_be15_det
    assert ms_neg > 0 and la_neg > 0
    assert len(pstat3) == 36
    assert all(int((pstat3["MEMBRANA"] == m).sum()) == 12 for m in (1, 2, 3))

    print("== 01_generar_sinteticos.py ==")
    print(f"destino: {dst}")
    print(f"  {cfg.ARCHIVO_QPCR:<32} {len(qpcr):>4} filas  sha256(tsv)={h_qpcr[:12]}")
    print(f"  {cfg.ARCHIVO_PSTAT3:<32} {len(pstat3):>4} filas  sha256(tsv)={h_pstat3[:12]}")
    print(f"  {cfg.ARCHIVO_ELISA:<32} {len(g_elisa):>4}+{alto_curva} filas  "
          f"sha256(tsv)={h_elisa[:12]} / {h_curva[:12]}")
    print(f"  CT==40: {ct40}   CT vacio: {ctna}   il6/BRAIN_E15 detectados: {il6_be15_det}")
    print(f"  ELISA Conc<0: MS={ms_neg}  LA={la_neg}")
    print("  invariantes OK")


if __name__ == "__main__":
    main()
