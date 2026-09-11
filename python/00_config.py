# 00_config.py -- Constantes, rutas y utilidades deterministas compartidas.
#
# Por que existe este archivo: el pipeline se implementa dos veces (R y Python) y
# todo script arranca importando esta configuracion. Aca se fijan (a) las rutas,
# siempre relativas a la raiz del repo resuelta desde la ubicacion del script,
# (b) la semilla global, (c) las listas de genes / grupos / tejidos, y (d) un
# generador pseudoaleatorio propio que produce EXACTAMENTE la misma secuencia en
# R y en Python, para que los datos sinteticos sean identicos entre lenguajes.
#
# Este archivo NO hace analisis. Su unico efecto colateral es crear las carpetas
# de salida si no existen (idempotente).

from __future__ import annotations

import math
from pathlib import Path

# ---------------------------------------------------------------------------
# Raiz del repo: se deduce del propio archivo (prohibido hardcodear rutas
# absolutas). python/00_config.py -> parents[1] == raiz del repo.
# ---------------------------------------------------------------------------
RAIZ_REPO = Path(__file__).resolve().parents[1]

# --- Rutas de datos y salidas (todas relativas a RAIZ_REPO) -----------------
RUTA_DATOS_RAW = RAIZ_REPO / "data" / "raw"
RUTA_DATOS_SINT = RAIZ_REPO / "data" / "synthetic"
RUTA_DATOS_PROC = RAIZ_REPO / "data" / "processed"
RUTA_OUT = RAIZ_REPO / "outputs"
RUTA_FIGURAS = RUTA_OUT / "figures"
RUTA_TABLAS = RUTA_OUT / "tables"
RUTA_TABLAS_R = RUTA_TABLAS / "R"
RUTA_TABLAS_PY = RUTA_TABLAS / "python"
RUTA_INTERMEDIOS = RUTA_OUT / "intermediate"
RUTA_DOCS = RAIZ_REPO / "docs"
RUTA_LOGS = RAIZ_REPO / "logs"

# Crear lo regenerable si falta (un clon limpio trae solo data/synthetic).
for _d in (
    RUTA_DATOS_SINT,
    RUTA_DATOS_PROC,
    RUTA_FIGURAS,
    RUTA_TABLAS_R,
    RUTA_TABLAS_PY,
    RUTA_INTERMEDIOS,
    RUTA_DOCS,
    RUTA_LOGS,
):
    _d.mkdir(parents=True, exist_ok=True)

# --- Archivos de entrada y hojas ------------------------------------------
ARCHIVO_QPCR = "Raw data CTs.xlsx"
HOJA_QPCR = "Sheet1"
ARCHIVO_ELISA = "ELISA IL6 2026 Dosis 100.xlsx"
HOJA_ELISA = "Sueros y LA"
HOJA_CURVA = "CURVA IL6"
ARCHIVO_PSTAT3 = "pstat3 placenta.xlsx"
HOJA_PSTAT3 = "Sheet1"


def _forzar_sintetico() -> bool:
    """run_all.ps1 -FromSynthetic exporta MIA_LPS_FORZAR_SINTETICO=1 para correr
    el pipeline sobre data/synthetic/ aunque existan los crudos (util para el
    chequeo de reproducibilidad 'corre de punta a punta sin datos reales')."""
    import os
    return os.environ.get("MIA_LPS_FORZAR_SINTETICO", "").strip() not in ("", "0")


def ruta_datos(nombre_archivo: str) -> Path:
    """Devuelve el archivo de datos a usar: el crudo real de data/raw/ si existe,
    y si no el sintetico versionado de data/synthetic/. Nunca escribe en data/raw/.
    Un clon limpio (sin datos crudos) cae automaticamente al sintetico."""
    real = RUTA_DATOS_RAW / nombre_archivo
    if real.is_file() and not _forzar_sintetico():
        return real
    return RUTA_DATOS_SINT / nombre_archivo


def fuente_datos(nombre_archivo: str) -> str:
    """'real' si se usaria el crudo, 'sintetico' si se usaria el generado."""
    if _forzar_sintetico():
        return "sintetico"
    return "real" if (RUTA_DATOS_RAW / nombre_archivo).is_file() else "sintetico"


# --- Semilla global -----------------------------------------------------
SEMILLA = 20260101

# --- Genes y su clasificacion (orden canonico del archivo crudo) -----------
GENES = [
    "fatcd36",
    "fatp1",
    "fatp4",
    "glut1",
    "glut3",
    "gp130",
    "il6",
    "il6R",
    "slc38a1",
    "slc38a2",
]
GENES_TRANSPORTADORES = [
    "fatcd36",
    "fatp1",
    "fatp4",
    "glut1",
    "glut3",
    "slc38a1",
    "slc38a2",
]  # 7 transportadores de nutrientes
GENES_VIA_IL6 = ["gp130", "il6", "il6R"]  # 3 de la via IL-6
GEN_HOUSEKEEPING = "rsp29"

# Conjuntos de variables para los SPLOM de co-expresion del Acto 2 (T7,
# extension pedida en pedidos/cambios_acto2_correlaciones_por_sexo.md): el
# conjunto es distinto por tejido porque il6 es cuantificable en placenta pero
# no en cerebro (calibrador HEMBRA_CONTROL 0/9, D7); il6R no entra en ningun
# SPLOM (deteccion insuficiente: 3/7/3/3 en cerebro).
GENES_SPLOM_PLACENTA = GENES_TRANSPORTADORES + ["il6", "gp130"]  # 9 variables
GENES_SPLOM_BRAIN = GENES_TRANSPORTADORES + ["gp130"]  # 8 variables

assert set(GENES) == set(GENES_TRANSPORTADORES) | set(GENES_VIA_IL6)
assert len(GENES_TRANSPORTADORES) == 7 and len(GENES_VIA_IL6) == 3
assert len(set(GENES)) == 10
assert len(GENES_SPLOM_PLACENTA) == 9 and len(GENES_SPLOM_BRAIN) == 8
assert set(GENES_SPLOM_PLACENTA) == set(GENES_TRANSPORTADORES) | {"il6", "gp130"}
assert set(GENES_SPLOM_BRAIN) == set(GENES_TRANSPORTADORES) | {"gp130"}

# --- Tejidos ----------------------------------------------------------
TEJIDOS_E15 = ["PLACENTA_E15", "BRAIN_E15"]
TEJIDO_EXCLUIDO = "BRAIN_P1"  # se excluye en 02_ingesta_qc (se analiza en otro informe)
TEJIDOS_TODOS = TEJIDOS_E15 + [TEJIDO_EXCLUIDO]

# --- Diseno de grupos -------------------------------------------------
NIVELES_SEXO = ["HEMBRA", "MACHO"]
NIVELES_TTO = ["CONTROL", "LPS"]
# D1: el calibrador de la cuantificacion relativa es HEMBRA CONTROL.
CALIBRADOR = {"SEXO": "HEMBRA", "TTO": "CONTROL"}

EMOJI_HEMBRA = "♀"
EMOJI_MACHO = "♂"
_VS16 = "️"  # variation selector que traen las NOMINACION de los crudos


def etiqueta_grupo(sexo: str, tto: str) -> str:
    """Etiqueta GRUPO tal como aparece en los crudos, reproduciendo el espaciado
    inconsistente real: 'HEMBRA'+'CONTROL' va pegado, el resto lleva un espacio."""
    s = normalizar_sexo(sexo)
    t = normalizar_tto(tto)
    emo = EMOJI_HEMBRA if s == "HEMBRA" else EMOJI_MACHO
    tt = "Control" if t == "CONTROL" else "LPS"
    sep = "" if (s == "HEMBRA" and t == "CONTROL") else " "
    return f"{emo}{sep}{tt}"


def normalizar_sexo(x: str) -> str:
    """Cualquier capitalizacion -> 'HEMBRA' / 'MACHO'."""
    u = str(x).strip().upper()
    if u.startswith("H") or u.startswith("F") or u == EMOJI_HEMBRA:
        return "HEMBRA"
    if u.startswith("M") or u == EMOJI_MACHO:
        return "MACHO"
    raise ValueError(f"SEXO no reconocido: {x!r}")


def normalizar_tto(x: str) -> str:
    """Cualquier capitalizacion -> 'CONTROL' / 'LPS'."""
    u = str(x).strip().upper()
    if u.startswith("C"):
        return "CONTROL"
    if u.startswith("L"):
        return "LPS"
    raise ValueError(f"TTO no reconocido: {x!r}")


# --- Tolerancias de comparacion R <-> Python (convencion del proyecto) -----
TOL_ESTADISTICO = 1e-6  # estadisticos, medias, coeficientes
TOL_P_ITERATIVO = 1e-4  # p-valores de tests iterativos / permutacion


# ---------------------------------------------------------------------------
# Generador pseudoaleatorio propio -- PARIDAD EXACTA R <-> Python
#
# PROPIO. No se usa para nada cientifico: unicamente para fabricar los datos
# sinteticos de 01_generar_sinteticos de forma identica en ambos lenguajes.
# Implementacion: LCG de 32 bits (constantes de Numerical Recipes) para el
# stream uniforme, y metodo polar de Marsaglia para los normales (solo usa
# sqrt y log, que estan correctamente redondeados en R y en Python -> los bits
# coinciden; se evita sin/cos justamente por eso).
# ---------------------------------------------------------------------------
_LCG_A = 1664525
_LCG_C = 1013904223
_LCG_M = 4294967296  # 2**32


class RNG:
    __slots__ = ("_estado", "_cache", "_tiene_cache")

    def __init__(self, semilla: int):
        self._estado = int(semilla) % _LCG_M
        self._cache = 0.0
        self._tiene_cache = False

    def unif1(self) -> float:
        """Un uniforme en [0, 1)."""
        self._estado = (_LCG_A * self._estado + _LCG_C) % _LCG_M
        return self._estado / _LCG_M

    def norm1(self, media: float = 0.0, sd: float = 1.0) -> float:
        """Un normal(media, sd) por el metodo polar de Marsaglia con cache."""
        if self._tiene_cache:
            self._tiene_cache = False
            return media + sd * self._cache
        while True:
            v1 = 2.0 * self.unif1() - 1.0
            v2 = 2.0 * self.unif1() - 1.0
            w = v1 * v1 + v2 * v2
            if 0.0 < w < 1.0:
                break
        f = math.sqrt(-2.0 * math.log(w) / w)
        self._cache = v2 * f
        self._tiene_cache = True
        return media + sd * (v1 * f)

    def entero(self, lo: int, hi: int) -> int:
        """Entero uniforme en [lo, hi] (ambos inclusive)."""
        return lo + int(self.unif1() * (hi - lo + 1))

    def elegir_sin_reemplazo(self, opciones: list, k: int) -> list:
        """k elementos distintos de 'opciones', por Fisher-Yates parcial."""
        idx = list(range(len(opciones)))
        n = len(idx)
        for i in range(k):
            j = i + int(self.unif1() * (n - i))
            idx[i], idx[j] = idx[j], idx[i]
        return [opciones[i] for i in idx[:k]]


def nuevo_rng(semilla: int = SEMILLA) -> RNG:
    return RNG(semilla)


# ---------------------------------------------------------------------------
def _resumen() -> None:
    print("== 00_config.py ==")
    print(f"RAIZ_REPO        : {RAIZ_REPO}")
    print(f"SEMILLA          : {SEMILLA}")
    print(f"GENES ({len(GENES)})       : {', '.join(GENES)}")
    print(f"  transportadores: {', '.join(GENES_TRANSPORTADORES)}")
    print(f"  via IL-6       : {', '.join(GENES_VIA_IL6)}")
    print(f"  housekeeping   : {GEN_HOUSEKEEPING}")
    print(f"TEJIDOS_E15      : {', '.join(TEJIDOS_E15)}   (excluido: {TEJIDO_EXCLUIDO})")
    for s in NIVELES_SEXO:
        for t in NIVELES_TTO:
            print(f"  grupo {s:<6} {t:<7} -> {etiqueta_grupo(s, t)!r}")
    for f in (ARCHIVO_QPCR, ARCHIVO_ELISA, ARCHIVO_PSTAT3):
        print(f"  {f:<32} fuente actual: {fuente_datos(f)}")
    # chequeo minimo del RNG (mismo valor esperado en R)
    r = nuevo_rng(SEMILLA)
    u = [round(r.unif1(), 10) for _ in range(3)]
    print(f"RNG check unif1 x3: {u}")


if __name__ == "__main__":
    _resumen()
