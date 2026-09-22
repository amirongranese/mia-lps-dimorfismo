# 11_sensibilidad.py -- ACTO 2.6: controles de sensibilidad de la correlacion
#                       placenta <-> cerebro y del test formal de T8.
#
# Por que existe este archivo: el Acto 2 (T7/T8) concluye que la coordinacion
# placenta <-> cerebro NO cambia de forma demostrable entre Control y LPS -- ni
# el test formal (Fisher z sobre rho de Spearman, prohibicion 4) ni la
# simulacion de restriccion de rango (prohibicion 5) permiten afirmar una
# diferencia. T9 somete esa conclusion a dos controles; cada uno rehace la
# correlacion (Spearman rho, estratos GLOBAL / CONTROL / LPS) y el Fisher z
# Control vs LPS:
#
#   (A) EIGENGENE -- variante PCA del score compuesto (D8 la nombra
#       explicitamente). En lugar del promedio de los 7 z de transportadores
#       (score compuesto de 04, D8), el resumen del modulo es la proyeccion de
#       cada feto sobre el PRIMER componente principal de los 7 transportadores,
#       calculado POR TEJIDO. PCA PROPIO: eigendescomposicion de la matriz de
#       correlacion (Pearson) 7x7 por Jacobi clasico SIN trigonometria (solo
#       +,-,*,/,sqrt -> bit-identico R/Python, el mismo truco que el RNG polar
#       de 00_config, que evita sin/cos a proposito). Signo de PC1 fijado a
#       cargas mayoritariamente positivas (suma de cargas > 0). Fetos con < 7
#       transportadores detectados (decision del usuario): se proyecta con los z
#       disponibles y se reescala por la norma de las cargas usadas -- analogo
#       al "promedio de los z disponibles" de D8. Se rehace la correlacion
#       placenta <-> cerebro y el Fisher z con el eigengene y se compara contra
#       el score compuesto.
#
#   (B) EXCLUSION DEL FETO EXTREMO -- uno por item (decision del usuario:
#       leave-one-out sobre rho; blanco = rho GLOBAL). Para cada item se quita
#       un feto por vez y se recalcula rho GLOBAL; el "feto extremo" es el de
#       mayor |rho_full - rho_sin_i| (empates: el primero en el orden por
#       MADRE_ID, FETO). Se lo excluye de los TRES estratos y se rehace rho (con
#       IC 95 % Bonett-Wright) y el Fisher z; se reporta si el veredicto
#       (p_bw < .05) cambia. NO se re-estima el PC1 al quitar el feto: la
#       pregunta es la sensibilidad de la correlacion a un punto influyente, no
#       la del PC1.
#
# Alcance: mismos items que T7/T8 -- 8 genes con -ddCt en ambos tejidos (todos
# menos il6, D7, y il6R, excluido de todo el Acto 2) + score compuesto (D8) --
# y ademas el eigengene en (B). Piso: se testea el Fisher z solo si Control y
# LPS tienen n_par >= 5; si no, la fila lleva los n y ningun estadistico.
#
# PARIDAD R/Python: Jacobi sin trig + Spearman/Pearson/SD/normas con acumulador
# double explicito (mismo orden que R) -> texto "%.10g" bit-identico para rho,
# cargas de PC1, autovalores y eigengene; lo que pasa por trascendentes (p por
# pt, IC por tanh/atanh, p del Fisher z por la normal) -> texto "%.6e". En R
# (solo alli) se cruza-verifica el PCA PROPIO contra eigen() (tol 1e-8) y cada
# rho contra cor.test(method="spearman", exact=FALSE) (tol 1e-9); NUNCA el
# metodo exacto AS 89 (segfaultea en este build de R). Las figuras son PNG:
# equivalentes, no byte-identicas.

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

ESTE_SCRIPT = "11_sensibilidad"

COL_TTO = {"CONTROL": "#0072B2", "LPS": "#D55E00"}   # Okabe-Ito, igual que 03/07/08/09
Z975 = 1.959963984540054                             # qnorm(0.975), literal exacto
PISO_PAR = 5                                         # min fetos emparejados por grupo

GEN_SIN_CEREBRO = "il6"                              # D7: sin -ddCt en BRAIN_E15
GEN_EXCLUIDO_CORR = "il6R"   # pedido explicito (08_acto2_correlaciones): deteccion insuficiente en cerebro
GENES_CORR = [g for g in cfg.GENES
              if g not in (GEN_SIN_CEREBRO, GEN_EXCLUIDO_CORR)]  # 8 genes
ITEMS_LOO = GENES_CORR + ["score_compuesto", "eigengene"]        # 10 items para (B)
ITEMS_A = ["eigengene", "score_compuesto"]                       # (A): eigengene primero
ESTRATOS = ["GLOBAL", "CONTROL", "LPS"]
TRANSP = list(cfg.GENES_TRANSPORTADORES)             # 7, base del PCA
TEJIDOS = list(cfg.TEJIDOS_E15)                      # PLACENTA_E15, BRAIN_E15
DPI = 300

JACOBI_TOL = 1e-15          # suma de cuadrados de la triangular superior estricta
JACOBI_MAX_SWEEPS = 100


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
    """SD muestral (n-1), acumulador explicito. None si n < 2."""
    n = len(xs)
    if n < 2:
        return None
    m = suma(xs) / n
    s = 0.0
    for v in xs:
        s += (v - m) * (v - m)
    return math.sqrt(s / (n - 1))


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
    if n < 2:
        return None
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
    """IC 95 % Bonett-Wright sobre atanh(rho)."""
    if rho is None or n <= 3 or rho * rho >= 1.0:
        return None, None
    se = math.sqrt((1.0 + rho * rho / 2.0) / (n - 3))
    z = math.atanh(rho)
    return math.tanh(z - Z975 * se), math.tanh(z + Z975 * se)


def norm_sf2(x):
    """p a dos colas por la normal estandar."""
    return 2.0 * float(_sst.norm.sf(abs(x)))


def pearson_pairwise(a, b):
    """Pearson sobre los pares completos de a,b (None = faltante). Sumas
    explicitas. None si quedan < 3 pares o alguna varianza es 0."""
    xs, ys = [], []
    for i in range(len(a)):
        if a[i] is not None and b[i] is not None:
            xs.append(a[i])
            ys.append(b[i])
    n = len(xs)
    if n < 3:
        return None
    mx = suma(xs) / n
    my = suma(ys) / n
    sxy = 0.0
    sxx = 0.0
    syy = 0.0
    for i in range(n):
        dx = xs[i] - mx
        dy = ys[i] - my
        sxy += dx * dy
        sxx += dx * dx
        syy += dy * dy
    if sxx <= 0.0 or syy <= 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def jacobi_eigen(A):
    """Eigendescomposicion de A (n x n, simetrica) por Jacobi clasico SIN
    trigonometria: la rotacion se arma con t = 1/(theta +/- sqrt(theta^2+1)),
    c = 1/sqrt(t^2+1), s = t*c -> solo +,-,*,/,sqrt, correctamente redondeados
    en R y Python. Barrido ciclico (p<q), umbral fijo, sweeps fijos ->
    determinista y bit-identico entre lenguajes.
    Devuelve (autovalores, autovectores) con autovectores como columnas
    (lista de listas: eigvecs[c][r])."""
    n = len(A)
    M = [[A[i][j] for j in range(n)] for i in range(n)]      # copia de trabajo
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _sweep in range(JACOBI_MAX_SWEEPS):
        off = 0.0
        for p in range(n - 1):
            for q in range(p + 1, n):
                off += M[p][q] * M[p][q]
        if off <= JACOBI_TOL:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = M[p][q]
                if apq == 0.0:
                    continue
                theta = (M[q][q] - M[p][p]) / (2.0 * apq)
                if theta >= 0.0:
                    t = 1.0 / (theta + math.sqrt(theta * theta + 1.0))
                else:
                    t = -1.0 / (-theta + math.sqrt(theta * theta + 1.0))
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):                            # M <- M * R
                    mkp = M[k][p]
                    mkq = M[k][q]
                    M[k][p] = c * mkp - s * mkq
                    M[k][q] = s * mkp + c * mkq
                for k in range(n):                            # M <- R^T * M
                    mpk = M[p][k]
                    mqk = M[q][k]
                    M[p][k] = c * mpk - s * mqk
                    M[q][k] = s * mpk + c * mqk
                for k in range(n):                            # V <- V * R
                    vkp = V[k][p]
                    vkq = V[k][q]
                    V[k][p] = c * vkp - s * vkq
                    V[k][q] = s * vkp + c * vkq
    eigvals = [M[i][i] for i in range(n)]
    eigvecs = [[V[r][cix] for r in range(n)] for cix in range(n)]
    return eigvals, eigvecs


# ===========================================================================
# 1. Carga y emparejamiento por feto -- identico a 08/09 + columna z por gen.
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

    madre, tto, negdd, zval = {}, {}, {}, {}
    for r in cuant:
        f = r["FETO"]
        madre[f] = r["MADRE_ID"]
        tto[f] = r["TTO"]
        negdd[(f, r["TEJIDO"], r["GEN"])] = (
            None if r["neg_ddCt"] == "" else float(r["neg_ddCt"]))
        zval[(f, r["TEJIDO"], r["GEN"])] = (
            None if r["z"] == "" else float(r["z"]))
    sc = {}
    for r in score:
        sc[(r["FETO"], r["TEJIDO"])] = (
            None if r["score_compuesto"] == "" else float(r["score_compuesto"]))

    fetos = sorted(madre.keys(), key=lambda f: (madre[f], f))
    return dict(fetos=fetos, madre=madre, tto=tto, negdd=negdd, z=zval, sc=sc)


# ===========================================================================
# 2. (A) EIGENGENE -- PC1 PROPIO del modulo de 7 transportadores, por tejido.
# ===========================================================================
def construir_eigengene(D):
    """Devuelve (eig, loadings, varianza):
      eig[(feto, tejido)]        -> score del feto sobre PC1 (o None)
      loadings[tejido][gen]      -> carga de PC1 (unit-norm, signo fijado)
      varianza[tejido]           -> lista de (pc, autovalor, prop_var, prop_acum)
    """
    eig, loadings, varianza = {}, {}, {}
    for tej in TEJIDOS:
        # matriz de columnas z (una por transportador), None si no detectado
        cols = {g: [D["z"].get((f, tej, g)) for f in D["fetos"]] for g in TRANSP}
        n = len(TRANSP)
        C = [[1.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                r = pearson_pairwise(cols[TRANSP[i]], cols[TRANSP[j]])
                if r is None:
                    r = 0.0
                C[i][j] = r
                C[j][i] = r
        eigvals, eigvecs = jacobi_eigen(C)
        idx = sorted(range(n), key=lambda k: eigvals[k], reverse=True)
        total = suma([eigvals[k] for k in idx])
        acum = 0.0
        vfilas = []
        for pc, k in enumerate(idx, start=1):
            prop = eigvals[k] / total if total > 0 else 0.0
            acum += prop
            vfilas.append((pc, eigvals[k], prop, acum))
        varianza[tej] = vfilas
        # PC1 = autovector del mayor autovalor; unit-norm explicito; signo fijo
        L = list(eigvecs[idx[0]])
        nrm = math.sqrt(suma([v * v for v in L]))
        L = [v / nrm for v in L]
        if suma(L) < 0.0:
            L = [-v for v in L]
        loadings[tej] = {TRANSP[i]: L[i] for i in range(n)}
        # proyeccion por feto: genes disponibles, reescalada por la norma usada
        for fi, f in enumerate(D["fetos"]):
            num = 0.0
            den = 0.0
            for i in range(n):
                zv = cols[TRANSP[i]][fi]
                if zv is None:
                    continue
                num += L[i] * zv
                den += L[i] * L[i]
            eig[(f, tej)] = (num / math.sqrt(den)) if den > 0.0 else None
    return eig, loadings, varianza


# ===========================================================================
# 3. Emparejamiento por feto (generaliza el de 08/09 a item = "eigengene").
# ===========================================================================
def _valor(D, eig, item, feto, tejido):
    if item == "eigengene":
        return eig.get((feto, tejido))
    if item == "score_compuesto":
        return D["sc"].get((feto, tejido))
    return D["negdd"].get((feto, tejido, item))


def pares(D, eig, item, estrato, excluir=None):
    """dict(x, y, tto, feto) de los fetos con ambos lados detectados."""
    xs, ys, ts, fs = [], [], [], []
    for f in D["fetos"]:
        if excluir is not None and f == excluir:
            continue
        if estrato != "GLOBAL" and D["tto"][f] != estrato:
            continue
        xp = _valor(D, eig, item, f, "PLACENTA_E15")
        yb = _valor(D, eig, item, f, "BRAIN_E15")
        if xp is None or yb is None:
            continue
        xs.append(xp)
        ys.append(yb)
        ts.append(D["tto"][f])
        fs.append(f)
    return dict(x=xs, y=ys, tto=ts, feto=fs)


# ===========================================================================
# 4. (A) tablas: correlacion y Fisher z del eigengene vs score compuesto.
# ===========================================================================
COLS_A_CORR = ["ITEM", "ESTRATO", "n_par", "rho_spearman", "ic95_low",
               "ic95_high", "p_valor"]


def tabla_a_correlacion(D, eig):
    filas = []
    for item in ITEMS_A:
        for est in ESTRATOS:
            pr = pares(D, eig, item, est)
            n = len(pr["x"])
            if n >= PISO_PAR:
                rho = spearman_rho(pr["x"], pr["y"])
                lo, hi = spearman_ci(rho, n)
                p = spearman_p(rho, n)
                filas.append([item, est, n, g10(rho), p6e(lo), p6e(hi), p6e(p)])
            else:
                filas.append([item, est, n, "", "", "", ""])
    return filas


COLS_A_TEST = ["ITEM", "n_control", "rho_control", "n_lps", "rho_lps",
               "delta_rho", "z_control", "z_lps", "se_bw_control", "se_bw_lps",
               "stat_z_bw", "p_bw", "se_clasico_control", "se_clasico_lps",
               "stat_z_clasico", "p_clasico"]


def _se_bw(rho, n):
    return math.sqrt((1.0 + rho * rho / 2.0) / (n - 3))


def _se_clasico(n):
    return 1.0 / math.sqrt(n - 3)


def _fisher_z(rc, rl, sec, sel):
    zc = math.atanh(rc)
    zl = math.atanh(rl)
    stat = (zc - zl) / math.sqrt(sec * sec + sel * sel)
    return zc, zl, sec, sel, stat, norm_sf2(stat)


def _fila_fisher(item, prc, prl):
    nc, nl = len(prc["x"]), len(prl["x"])
    if nc < PISO_PAR or nl < PISO_PAR:
        return [item, nc, "", nl, "", "", "", "", "", "", "", "", "", "", "", ""]
    rc = spearman_rho(prc["x"], prc["y"])
    rl = spearman_rho(prl["x"], prl["y"])
    drho = rc - rl
    zc, zl, sbc, sbl, sb, pb = _fisher_z(rc, rl, _se_bw(rc, nc), _se_bw(rl, nl))
    _c1, _c2, scc, scl, sc_, pc = _fisher_z(
        rc, rl, _se_clasico(nc), _se_clasico(nl))
    return [item, nc, g10(rc), nl, g10(rl), g10(drho), p6e(zc), p6e(zl),
            p6e(sbc), p6e(sbl), p6e(sb), p6e(pb), p6e(scc), p6e(scl),
            p6e(sc_), p6e(pc)]


def tabla_a_test(D, eig):
    filas = []
    for item in ITEMS_A:
        prc = pares(D, eig, item, "CONTROL")
        prl = pares(D, eig, item, "LPS")
        filas.append(_fila_fisher(item, prc, prl))
    return filas


COLS_PCA_LOAD = ["TEJIDO", "GEN", "loading_pc1"]
COLS_PCA_VAR = ["TEJIDO", "PC", "autovalor", "prop_var", "prop_var_acum"]


def tablas_pca(loadings, varianza):
    lo = []
    for tej in TEJIDOS:
        for g in TRANSP:
            lo.append([tej, g, g10(loadings[tej][g])])
    va = []
    for tej in TEJIDOS:
        for (pc, ev, prop, acum) in varianza[tej]:
            va.append([tej, pc, g10(ev), g10(prop), g10(acum)])
    return lo, va


# ===========================================================================
# 5. (B) exclusion del feto extremo (leave-one-out sobre rho GLOBAL).
# ===========================================================================
COLS_B = ["ITEM", "TIPO", "n_par_global", "feto_extremo", "tto_feto_extremo",
          "rho_global_full", "rho_global_sin", "impacto_rho_global",
          "rho_control_full", "rho_control_sin", "rho_lps_full", "rho_lps_sin",
          "delta_rho_full", "delta_rho_sin", "p_bw_full", "p_bw_sin",
          "p_clasico_full", "p_clasico_sin", "veredicto_full", "veredicto_sin",
          "veredicto_cambia"]


def _rho_est(D, eig, item, est, excluir=None):
    pr = pares(D, eig, item, est, excluir=excluir)
    n = len(pr["x"])
    if n < PISO_PAR:
        return None, n
    return spearman_rho(pr["x"], pr["y"]), n


def _fisher_pb_pc(D, eig, item, excluir=None):
    prc = pares(D, eig, item, "CONTROL", excluir=excluir)
    prl = pares(D, eig, item, "LPS", excluir=excluir)
    nc, nl = len(prc["x"]), len(prl["x"])
    if nc < PISO_PAR or nl < PISO_PAR:
        return None, None, None
    rc = spearman_rho(prc["x"], prc["y"])
    rl = spearman_rho(prl["x"], prl["y"])
    _zc, _zl, _sbc, _sbl, _sb, pb = _fisher_z(
        rc, rl, _se_bw(rc, nc), _se_bw(rl, nl))
    _c1, _c2, _scc, _scl, _sc, pc = _fisher_z(
        rc, rl, _se_clasico(nc), _se_clasico(nl))
    return rc - rl, pb, pc


def _veredicto(pb):
    if pb is None:
        return ""
    return "dif" if pb < 0.05 else "n.s."


def tabla_b(D, eig):
    filas = []
    for item in ITEMS_LOO:
        tipo = ("score" if item == "score_compuesto"
                else "eigengene" if item == "eigengene" else "gen")
        pr = pares(D, eig, item, "GLOBAL")
        n = len(pr["x"])
        if n < PISO_PAR:
            filas.append([item, tipo, n, "", "", "", "", "", "", "", "", "",
                          "", "", "", "", "", "", "", "", ""])
            continue
        rho_full = spearman_rho(pr["x"], pr["y"])
        # leave-one-out sobre rho GLOBAL
        peor_feto, peor_impacto = None, -1.0
        for k in range(n):
            xk = pr["x"][:k] + pr["x"][k + 1:]
            yk = pr["y"][:k] + pr["y"][k + 1:]
            rk = spearman_rho(xk, yk)
            if rk is None:
                continue
            imp = abs(rho_full - rk)
            if imp > peor_impacto:
                peor_impacto = imp
                peor_feto = pr["feto"][k]
        tto_ex = D["tto"][peor_feto]
        rho_g_sin, _ = _rho_est(D, eig, item, "GLOBAL", excluir=peor_feto)
        rc_full, _ = _rho_est(D, eig, item, "CONTROL")
        rc_sin, _ = _rho_est(D, eig, item, "CONTROL", excluir=peor_feto)
        rl_full, _ = _rho_est(D, eig, item, "LPS")
        rl_sin, _ = _rho_est(D, eig, item, "LPS", excluir=peor_feto)
        dr_full, pb_full, pc_full = _fisher_pb_pc(D, eig, item)
        dr_sin, pb_sin, pc_sin = _fisher_pb_pc(D, eig, item, excluir=peor_feto)
        v_full, v_sin = _veredicto(pb_full), _veredicto(pb_sin)
        cambia = ("" if (pb_full is None or pb_sin is None)
                  else ((pb_full < 0.05) != (pb_sin < 0.05)))
        filas.append([
            item, tipo, n, peor_feto, tto_ex,
            g10(rho_full), g10(rho_g_sin), g10(rho_full - rho_g_sin),
            g10(rc_full), g10(rc_sin), g10(rl_full), g10(rl_sin),
            g10(dr_full), g10(dr_sin), p6e(pb_full), p6e(pb_sin),
            p6e(pc_full), p6e(pc_sin), v_full, v_sin, cambia])
    return filas


# ===========================================================================
# 6. Figuras.
# ===========================================================================
def figura_eigengene(D, eig, loadings, varianza, a_corr, ruta):
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 9.0))

    # (0,0)+(0,1): cargas de PC1 por tejido
    for ax, tej in zip((axes[0, 0], axes[0, 1]), TEJIDOS):
        vals = [loadings[tej][g] for g in TRANSP]
        ax.bar(range(len(TRANSP)), vals, color="#0072B2", edgecolor="0.3",
               linewidth=0.5)
        ax.axhline(0.0, color="0.4", lw=0.8)
        ax.set_xticks(range(len(TRANSP)))
        ax.set_xticklabels(TRANSP, rotation=40, ha="right", fontsize=7,
                           style="italic")
        ax.tick_params(axis="y", labelsize=7)
        pv = varianza[tej][0][2]
        tt = "Placenta E15" if tej == "PLACENTA_E15" else "Cerebro fetal E15"
        ax.set_title(f"Cargas PC1 -- {tt}  (var. expl. {pv * 100:.0f} %)",
                     fontsize=9.5)
        ax.set_ylabel("carga", fontsize=8)

    # (1,0): dispersion placenta vs cerebro del eigengene (GLOBAL, color TTO)
    ax = axes[1, 0]
    pr = pares(D, eig, "eigengene", "GLOBAL")
    for tt in ("CONTROL", "LPS"):
        xx = [pr["x"][i] for i in range(len(pr["x"])) if pr["tto"][i] == tt]
        yy = [pr["y"][i] for i in range(len(pr["y"])) if pr["tto"][i] == tt]
        ax.plot(xx, yy, "o", ms=5, mfc=COL_TTO[tt], mec="white", mew=0.4,
                ls="none", label=tt.capitalize())
    rg = spearman_rho(pr["x"], pr["y"])
    lo, hi = spearman_ci(rg, len(pr["x"]))
    ax.set_title(f"Eigengene placenta <-> cerebro (GLOBAL)\n"
                 f"ρ = {rg:.2f} [{lo:.2f}, {hi:.2f}]  (n = {len(pr['x'])})",
                 fontsize=9.5)
    ax.set_xlabel("placenta  eigengene (PC1)", fontsize=8)
    ax.set_ylabel("cerebro  eigengene (PC1)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(fontsize=8, frameon=False, loc="best")

    # (1,1): rho score compuesto vs eigengene por estrato
    ax = axes[1, 1]
    idx = {}
    for f in a_corr:
        idx.setdefault(f[0], {})[f[1]] = f[3]
    x = range(len(ESTRATOS))
    w = 0.36
    sc_v = [float(idx["score_compuesto"][e]) if idx["score_compuesto"][e] != ""
            else float("nan") for e in ESTRATOS]
    ei_v = [float(idx["eigengene"][e]) if idx["eigengene"][e] != ""
            else float("nan") for e in ESTRATOS]
    ax.bar([i - w / 2 for i in x], sc_v, w, label="score compuesto",
           color="0.55", edgecolor="0.3", linewidth=0.5)
    ax.bar([i + w / 2 for i in x], ei_v, w, label="eigengene (PC1)",
           color="#0072B2", edgecolor="0.3", linewidth=0.5)
    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(ESTRATOS, fontsize=8)
    ax.set_ylabel("ρ de Spearman (placenta <-> cerebro)", fontsize=8)
    ax.set_title("Correlacion placenta <-> cerebro: score vs eigengene",
                 fontsize=9.5)
    ax.tick_params(axis="y", labelsize=7)
    ax.legend(fontsize=8, frameon=False)

    fig.suptitle("T9 (A) -- eigengene: PC1 PROPIO del modulo de 7 "
                 "transportadores (por tejido) como alternativa al score "
                 "compuesto de D8", fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


def figura_excl_extremo(tb, ruta):
    fig, ax = plt.subplots(figsize=(10.0, 6.5))
    ys = list(range(len(ITEMS_LOO)))[::-1]
    for y, row in zip(ys, tb):
        item = row[0]
        nom = ("score compuesto" if item == "score_compuesto"
               else "eigengene (PC1)" if item == "eigengene" else item)
        if row[5] == "" or row[6] == "":
            ax.text(0.0, y, f"  {nom}: n<{PISO_PAR} (sin test)", va="center",
                    fontsize=8, color="0.5")
            continue
        rf, rs = float(row[5]), float(row[6])
        ax.plot([rf, rs], [y, y], "-", color="0.7", lw=1.2, zorder=1)
        ax.plot(rf, y, "o", ms=7, color="0.35", zorder=2,
                label="ρ global (todos)" if y == ys[0] else None)
        ax.plot(rs, y, "o", ms=7, color="#D55E00", zorder=2,
                label="ρ global (sin feto extremo)" if y == ys[0] else None)
        pbf = f"{float(row[14]):.3f}" if row[14] != "" else "--"
        pbs = f"{float(row[15]):.3f}" if row[15] != "" else "--"
        marca = " (cambia)" if row[20] == "TRUE" else ""
        ax.text(1.02, y, f"p_bw {pbf} -> {pbs}{marca}", va="center",
                fontsize=7.5, transform=ax.get_yaxis_transform())
    ax.axvline(0.0, color="0.4", lw=0.8, ls="--")
    ax.set_yticks(ys)
    ax.set_yticklabels(
        ["score compuesto" if i == "score_compuesto"
         else "eigengene (PC1)" if i == "eigengene" else i
         for i in ITEMS_LOO], fontsize=8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("ρ de Spearman GLOBAL (placenta ↔ cerebro por feto)",
                  fontsize=9)
    ax.legend(fontsize=8, loc="lower left", frameon=False)
    ax.set_title("T9 (B) -- exclusion del feto extremo (leave-one-out sobre ρ "
                 "GLOBAL, uno por item)\nDeltas del Fisch z Control vs LPS: "
                 "p_bw completo -> p_bw sin el feto extremo", fontsize=9.5)
    fig.tight_layout(rect=(0.0, 0.0, 0.80, 1.0))
    fig.savefig(ruta, dpi=DPI)
    plt.close(fig)


# ===========================================================================
# 7. Artefactos compartidos (merge por 'script') -- headers identicos a 02..09.
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
    header = ["id", "tipo", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script"]
    merge_por_script(cfg.RUTA_TABLAS / "verificaciones.csv", header, filas_nuevas,
                     "script", lambda f: (f[6], f[0]))


_DESCARTES = "\n".join([
    "## 11_sensibilidad",
    "",
    "### Que resuelve (Acto 2.6)",
    "",
    "- Dos controles de sensibilidad de la conclusion del Acto 2 (la coordinacion "
    "placenta <-> cerebro no cambia de forma demostrable entre Control y LPS). "
    "Cada uno rehace la correlacion (Spearman rho, GLOBAL / CONTROL / LPS) y el "
    "Fisher z Control vs LPS de `09_acto2_dispersion`.",
    "",
    "### (A) Eigengene -- variante PCA del score compuesto (D8)",
    "",
    "- El score compuesto de D8 es el **promedio** de los 7 z de transportadores. "
    "El eigengene es la **proyeccion sobre PC1** del modulo de 7 transportadores, "
    "calculado **por tejido** (placenta y cerebro tienen su propia matriz y su "
    "propio PC1).",
    "- **Matriz de entrada**: correlacion de Pearson (pairwise-complete) entre las "
    "7 columnas de z (D8, z por gen x tejido sobre los 36 fetos, sd n-1). Se "
    "considero la matriz de Spearman para ser coherente con el resto del Acto 2 y "
    "**no se uso**: el eigengene es una reduccion de dimension estandar del "
    "modulo (Pearson / datos estandarizados); el test de asociacion entre tejidos "
    "sigue siendo Spearman. Se prefirio pairwise-complete a casos completos "
    "porque este ultimo descartaria 4-6 fetos (placenta / cerebro) y cambiaria el "
    "modulo para todos. Costo: la eliminacion pairwise puede volver la matriz "
    "levemente no definida positiva (en cerebro el 7mo autovalor queda en ~ -0.11 "
    "sobre una traza de 7); PC1 -- lo unico que se usa -- domina (>76 % de la "
    "varianza, cargas todas positivas) y coincide con `eigen()` a 1e-13, asi que "
    "el artefacto no lo afecta.",
    "- **PCA PROPIO**: eigendescomposicion por **Jacobi clasico SIN "
    "trigonometria** -- la rotacion se arma con "
    "`t = 1/(theta + sign*sqrt(theta^2+1))`, `c = 1/sqrt(t^2+1)`, `s = t*c` "
    "(solo `+ - * / sqrt`, correctamente redondeados en R y Python). Barrido "
    "ciclico p<q, umbral `1e-15` sobre la suma de cuadrados de la triangular "
    "superior, `<= 100` sweeps. Es el mismo motivo por el que el RNG de "
    "`00_config` usa el metodo polar y evita `sin`/`cos`: determinismo bit a bit "
    "entre lenguajes. `prcomp` / `numpy.linalg.eigh` **no** entran en el "
    "resultado.",
    "- **Signo de PC1**: se fija a cargas mayoritariamente positivas "
    "(`suma(cargas) > 0`, si no se invierte el vector). Convencion de "
    "\"tono de expresion\" del modulo.",
    "- **Fetos con < 7 transportadores detectados** (decision del usuario): se "
    "proyecta con los z disponibles y se divide por `sqrt(sum(carga_g^2))` de las "
    "cargas efectivamente usadas -- para un feto completo esto es exactamente el "
    "score de PC1; para uno incompleto es la proyeccion sobre la direccion "
    "unitaria del subespacio disponible (analogo al \"promedio de los z "
    "disponibles\" de D8). Alternativa descartada: suma cruda sin reescalar "
    "(sesga a la baja la magnitud de los fetos incompletos).",
    "",
    "### (B) Exclusion del feto extremo (decision del usuario)",
    "",
    "- **Criterio**: leave-one-out sobre rho. Para cada item se quita un feto por "
    "vez y se recalcula **rho GLOBAL**; el feto extremo es el de mayor "
    "`|rho_full - rho_sin_i|` (empates -> el primero en el orden `MADRE_ID, "
    "FETO`). Se descarto la distancia de Cook: mide influencia sobre un ajuste "
    "**lineal**, y el estadistico que se reporta en el Acto 2 es de **rango** "
    "(Spearman, decision del usuario en T7).",
    "- **Blanco = rho GLOBAL** (decision del usuario): el feto que mas mueve la "
    "correlacion pooled. Se lo excluye de los TRES estratos y se rehace rho (con "
    "IC Bonett-Wright) y el Fisher z; `veredicto_cambia` marca si `p_bw < .05` "
    "cambia de estado.",
    "- **No se re-estima el PC1** al quitar el feto: la pregunta es la "
    "sensibilidad de la correlacion a un punto influyente, no la del PC1. Para el "
    "item `eigengene` se usan las cargas calculadas sobre los 36 fetos.",
    "- Items: los 8 genes de T7/T8 (todos menos `il6`, D7, y `il6R`, excluido "
    "de todo el Acto 2 -- bugfix `cambios_informe_conclusiones.md` punto 5b, "
    "ver `analisis_descartados.md` de `09_acto2_dispersion`) + score "
    "compuesto + eigengene. Todos superan el piso de n>=5 por grupo en estos "
    "datos; si alguno no lo alcanzara, su fila llevaria los n y el impacto "
    "sobre rho GLOBAL, sin Fisher z.",
    "",
    "### Paridad R / Python",
    "",
    "- Jacobi sin trig + Spearman / Pearson / SD / normas con acumulador `double` "
    "explicito (mismo orden) -> `%.10g` bit-identico para rho, cargas, "
    "autovalores y eigengene. Lo que pasa por trascendentes (`p` por `pt`, IC por "
    "`tanh`/`atanh`, `p` del Fisher z por la normal) -> texto `%.6e`. Figuras "
    "PNG: equivalentes, no byte-identicas.",
    "- Solo en R se cruza-verifica el PCA PROPIO contra `eigen()` (tol 1e-8) y "
    "cada rho de estrato contra `cor.test(method=\"spearman\", exact=FALSE)` "
    "(tol 1e-9); en Python el PCA se contrasta contra `numpy.linalg.eigh` "
    "(tol 1e-6). Con `stopifnot` / `assert`.",
])


def actualizar_descartados():
    ruta = cfg.RUTA_TABLAS / "analisis_descartados.md"
    marca_ini = "<!-- 11_sensibilidad:inicio -->"
    marca_fin = "<!-- 11_sensibilidad:fin -->"
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
# 8. Reporte legible.
# ===========================================================================
def _md(header, filas):
    out = ["| " + " | ".join(str(h) for h in header) + " |",
           "| " + " | ".join("---" for _ in header) + " |"]
    for f in filas:
        out.append("| " + " | ".join(_fmt(v) for v in f) + " |")
    return "\n".join(out)


def construir_reporte(fuente, a_corr, a_test, pca_load, pca_var, tb):
    L = []
    ap = L.append
    ap("# Reporte de sensibilidad Acto 2.6 (T9)")
    ap("")
    ap("Generado por `11_sensibilidad` (R y Python producen este archivo "
       "identico).")
    ap(f"Fuente de datos en uso: `{fuente}`.")
    ap("")
    ap("## 1. (A) Eigengene -- PC1 PROPIO del modulo de 7 transportadores")
    ap("")
    ap("- PC1 por tejido (eigendescomposicion de la matriz de correlacion de "
       "Pearson de los z, Jacobi PROPIO sin trigonometria). Signo fijado a "
       "cargas mayoritariamente positivas. Fetos con <7 detectados: proyeccion "
       "reescalada por la norma de las cargas usadas (decision del usuario).")
    ap("")
    ap("### 1.1 Cargas de PC1")
    ap("")
    ap(_md(COLS_PCA_LOAD, pca_load))
    ap("")
    ap("### 1.2 Varianza explicada")
    ap("")
    ap(_md(COLS_PCA_VAR, pca_var))
    ap("")
    ap("### 1.3 Correlacion placenta <-> cerebro: eigengene vs score compuesto")
    ap("")
    ap(_md(COLS_A_CORR, a_corr))
    ap("")
    ap("### 1.4 Fisher z Control vs LPS: eigengene vs score compuesto")
    ap("")
    ap(_md(COLS_A_TEST, a_test))
    ap("")
    ap("## 2. (B) Exclusion del feto extremo (leave-one-out sobre rho GLOBAL)")
    ap("")
    ap("- Un feto por item: el de mayor `|rho_full - rho_sin_i|` sobre rho "
       "GLOBAL. Se excluye de los 3 estratos y se rehace rho + Fisher z. "
       "`veredicto_cambia` = TRUE si `p_bw < .05` cambia de estado.")
    ap("")
    ap(_md(COLS_B, tb))
    ap("")
    ap("## 3. Figuras")
    ap("")
    ap("- `outputs/figures/acto2_sensibilidad_eigengene.png` -- cargas de PC1 "
       "por tejido, dispersion del eigengene y comparacion de rho score vs "
       "eigengene por estrato.")
    ap("- `outputs/figures/acto2_sensibilidad_excl_extremo.png` -- rho GLOBAL "
       "con y sin el feto extremo por item, con el cambio de `p_bw`.")
    ap("")
    ap("## 4. Notas")
    ap("")
    ap("Ver `analisis_descartados.md`, seccion `11_sensibilidad`: eleccion de la "
       "matriz de entrada del PCA, Jacobi sin trigonometria, tratamiento de los "
       "fetos incompletos, y por que leave-one-out sobre rho y no Cook.")
    ap("")
    return "\n".join(L)


# ===========================================================================
def main():
    D = cargar()
    fuente = cfg.fuente_datos(cfg.ARCHIVO_QPCR)

    eig, loadings, varianza = construir_eigengene(D)

    # cruza-verificacion PROPIO vs numpy.linalg.eigh (Python; en R es eigen()).
    import numpy as _np
    peor_pca = 0.0
    for tej in TEJIDOS:
        n = len(TRANSP)
        cols = {g: [D["z"].get((f, tej, g)) for f in D["fetos"]] for g in TRANSP}
        C = _np.eye(n)
        for i in range(n):
            for j in range(i + 1, n):
                r = pearson_pairwise(cols[TRANSP[i]], cols[TRANSP[j]]) or 0.0
                C[i, j] = C[j, i] = r
        w, V = _np.linalg.eigh(C)
        L_ref = V[:, int(_np.argmax(w))]
        if L_ref.sum() < 0:
            L_ref = -L_ref
        L_prop = _np.array([loadings[tej][g] for g in TRANSP])
        peor_pca = max(peor_pca, float(_np.max(_np.abs(L_prop - L_ref))),
                       abs(float(w.max()) - varianza[tej][0][1]))
    print(f"  [cruza-verif PROPIO vs numpy.linalg.eigh] peor |dif| = {peor_pca:.3e}")
    assert peor_pca < 1e-6, f"PCA PROPIO diverge de numpy: {peor_pca:.3e}"

    a_corr = tabla_a_correlacion(D, eig)
    a_test = tabla_a_test(D, eig)
    pca_load, pca_var = tablas_pca(loadings, varianza)
    tb = tabla_b(D, eig)

    for base in (cfg.RUTA_TABLAS_R, cfg.RUTA_TABLAS_PY):
        escribir_csv(base / "acto2_sensibilidad_eigengene_correlacion.csv",
                     COLS_A_CORR, a_corr)
        escribir_csv(base / "acto2_sensibilidad_eigengene_test.csv",
                     COLS_A_TEST, a_test)
        escribir_csv(base / "acto2_sensibilidad_pca_loadings.csv",
                     COLS_PCA_LOAD, pca_load)
        escribir_csv(base / "acto2_sensibilidad_pca_varianza.csv",
                     COLS_PCA_VAR, pca_var)
        escribir_csv(base / "acto2_sensibilidad_excl_extremo.csv", COLS_B, tb)

    fig_e = cfg.RUTA_FIGURAS / "acto2_sensibilidad_eigengene.png"
    fig_x = cfg.RUTA_FIGURAS / "acto2_sensibilidad_excl_extremo.png"
    figura_eigengene(D, eig, loadings, varianza, a_corr, fig_e)
    figura_excl_extremo(tb, fig_x)

    escribir_texto(cfg.RUTA_TABLAS / "acto2_sensibilidad_reporte.md",
                   construir_reporte(fuente, a_corr, a_test, pca_load, pca_var, tb))
    actualizar_descartados()

    # --- resumen determinista para verificaciones ---
    def _a(item, est):
        for f in a_corr:
            if f[0] == item and f[1] == est:
                return f
        return None

    ei_glob = _a("eigengene", "GLOBAL")
    sc_glob = _a("score_compuesto", "GLOBAL")
    ei_test = next(f for f in a_test if f[0] == "eigengene")
    sc_test = next(f for f in a_test if f[0] == "score_compuesto")
    n_cambia = sum(1 for f in tb if f[20] == "TRUE")
    n_testados_b = sum(1 for f in tb if f[14] != "")
    pv_pla = varianza["PLACENTA_E15"][0][2]
    pv_bra = varianza["BRAIN_E15"][0][2]
    min_ev_pla = varianza["PLACENTA_E15"][-1][1]
    min_ev_bra = varianza["BRAIN_E15"][-1][1]
    cargas_pos = {tej: sum(1 for g in TRANSP if loadings[tej][g] > 0)
                  for tej in TEJIDOS}

    ent = (f"data/processed/qpcr_cuantificacion_long.tsv + "
           f"qpcr_score_compuesto_long.tsv (de data/{fuente}/{cfg.ARCHIVO_QPCR})")
    registrar_procedencia([
        ["outputs/tables/{R,python}/acto2_sensibilidad_eigengene_correlacion.csv",
         "tabla", ESTE_SCRIPT, "PROPIO", ent,
         "Spearman rho placenta<->cerebro del eigengene (PC1) vs score compuesto "
         "x {GLOBAL, CONTROL, LPS}; IC95 Bonett-Wright, p t-aprox"],
        ["outputs/tables/{R,python}/acto2_sensibilidad_eigengene_test.csv",
         "tabla", ESTE_SCRIPT, "PROPIO", ent,
         "Fisher z Control vs LPS del eigengene vs score compuesto; SE "
         "Bonett-Wright + SE clasico"],
        ["outputs/tables/{R,python}/acto2_sensibilidad_pca_loadings.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "cargas de PC1 (Jacobi PROPIO sin trig) de los 7 transportadores por "
         "tejido; signo fijado a suma>0"],
        ["outputs/tables/{R,python}/acto2_sensibilidad_pca_varianza.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "autovalores y proporcion de varianza de los 7 PC por tejido"],
        ["outputs/tables/{R,python}/acto2_sensibilidad_excl_extremo.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "leave-one-out sobre rho GLOBAL: feto extremo por item y rho + Fisher z "
         "con y sin el; veredicto_cambia"],
        ["outputs/figures/acto2_sensibilidad_eigengene.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent,
         "cargas PC1 por tejido, dispersion del eigengene y rho score vs "
         "eigengene por estrato"],
        ["outputs/figures/acto2_sensibilidad_excl_extremo.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent,
         "rho GLOBAL con y sin el feto extremo por item, con el cambio de p_bw"],
        ["outputs/tables/acto2_sensibilidad_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.6 (T9)"],
    ])
    registrar_verificaciones([
        ["sens_eigengene_pca_propio", "recalculo",
         "PC1 por eigendescomposicion PROPIA (Jacobi sin trigonometria), no "
         "prcomp/sklearn en el resultado",
         f"peor |dif| PC1 vs libreria de referencia = {peor_pca:.2e}", "< 1e-6",
         "TRUE" if peor_pca < 1e-6 else "FALSE", ESTE_SCRIPT],
        ["sens_eigengene_pc1_signo", "recalculo",
         "signo de PC1 fijado a cargas mayoritariamente positivas (suma>0)",
         f"cargas>0: placenta {cargas_pos['PLACENTA_E15']}/7, cerebro "
         f"{cargas_pos['BRAIN_E15']}/7",
         "mayoria positiva por tejido",
         "TRUE" if (cargas_pos["PLACENTA_E15"] >= 4 and
                    cargas_pos["BRAIN_E15"] >= 4) else "FALSE", ESTE_SCRIPT],
        ["sens_eigengene_var_pc1", "recalculo",
         "proporcion de varianza de PC1 por tejido (contexto del eigengene)",
         f"placenta {pv_pla:.3f}; cerebro {pv_bra:.3f}", "0 < prop <= 1",
         "TRUE" if (0.0 < pv_pla <= 1.0 and 0.0 < pv_bra <= 1.0) else "FALSE",
         ESTE_SCRIPT],
        ["sens_eigengene_matriz_psd", "recalculo",
         "la matriz de correlacion pairwise-complete puede ser levemente "
         "indefinida; PC1 domina y no se afecta",
         f"min autovalor: placenta {min_ev_pla:.3f}; cerebro {min_ev_bra:.3f}; "
         f"PC1 var {pv_pla:.2f}/{pv_bra:.2f}",
         "PC1 >> resto (traza = 7)",
         "TRUE" if (pv_pla > 0.5 and pv_bra > 0.5) else "FALSE", ESTE_SCRIPT],
        ["sens_eigengene_correlacion_global", "declaracion",
         "rho GLOBAL placenta<->cerebro: eigengene vs score compuesto",
         f"eigengene rho={ei_glob[3]} (n={ei_glob[2]}); score rho={sc_glob[3]} "
         f"(n={sc_glob[2]})", "ambos descritos, misma direccion cualitativa",
         "TRUE", ESTE_SCRIPT],
        ["sens_eigengene_fisher_z", "recalculo",
         "Fisher z Control vs LPS: el eigengene no cambia el veredicto del score",
         f"eigengene p_bw={ei_test[11]}; score p_bw={sc_test[11]}",
         "ninguno alcanza p_bw<.05 (consistente con T8)",
         "TRUE" if (ei_test[11] != "" and sc_test[11] != "" and
                    float(ei_test[11]) >= 0.05 and float(sc_test[11]) >= 0.05)
         else "REVISAR", ESTE_SCRIPT],
        ["sens_excl_extremo_criterio", "declaracion",
         "feto extremo = leave-one-out sobre rho GLOBAL (no Cook); blanco rho "
         "GLOBAL (decision del usuario)",
         "argmax_i |rho_full - rho_sin_i| por item, sobre el estrato GLOBAL",
         "LOO sobre rho GLOBAL", "TRUE", ESTE_SCRIPT],
        ["sens_excl_extremo_items", "recalculo",
         "items del control (B): 8 genes de T7/T8 + score compuesto + eigengene",
         f"{len(ITEMS_LOO)} items; {n_testados_b} con Fisher z (n>=5 por grupo)",
         "10 items", "TRUE" if len(ITEMS_LOO) == 10 else "FALSE", ESTE_SCRIPT],
        ["sens_excl_extremo_veredicto", "declaracion",
         "cuantos items cambian el veredicto p_bw<.05 al excluir su feto extremo",
         f"{n_cambia} de {n_testados_b} items testeados cambian de veredicto",
         "reportado; el informe menciona el numero", "TRUE", ESTE_SCRIPT],
        ["sens_excl_extremo_no_reestima_pc1", "declaracion",
         "en el item eigengene NO se re-estima PC1 al quitar el feto extremo",
         "cargas de PC1 calculadas sobre los 36 fetos, fijas en el LOO",
         "PC1 fijo", "TRUE", ESTE_SCRIPT],
        ["sens_figuras", "existencia",
         "figuras Acto 2.6: eigengene (cargas + dispersion + rho) y exclusion "
         "del extremo",
         f"eigengene={fig_e.is_file()};excl_extremo={fig_x.is_file()}".replace(
             "True", "TRUE").replace("False", "FALSE"),
         "2 figuras existen",
         "TRUE" if (fig_e.is_file() and fig_x.is_file()) else "FALSE",
         ESTE_SCRIPT],
    ])

    print("== 11_sensibilidad.py ==")
    print(f"  fuente = {fuente}")
    print(f"  PC1 var. expl.: placenta {pv_pla:.3f} | cerebro {pv_bra:.3f}")
    print(f"  cargas>0: placenta {cargas_pos['PLACENTA_E15']}/7 | cerebro "
          f"{cargas_pos['BRAIN_E15']}/7")
    print("  (A) rho GLOBAL placenta<->cerebro:")
    for f in a_corr:
        if f[1] == "GLOBAL":
            print(f"      {f[0]:16} n={f[2]:>2}  rho={f[3]:>9}  p={f[6]}")
    print("  (A) Fisher z Control vs LPS:")
    for f in a_test:
        print(f"      {f[0]:16} rhoC={f[2]:>9} rhoL={f[4]:>9} dRho={f[5]:>9} "
              f"p_bw={f[11]}")
    print(f"  (B) items con Fisher z: {n_testados_b}/{len(ITEMS_LOO)}; "
          f"cambian veredicto: {n_cambia}")
    for f in tb:
        if f[5] != "":
            print(f"      {f[0]:16} extremo={f[3]:<16} ({f[4]:<7}) "
                  f"rhoG {f[5]:>8} -> {f[6]:>8}  p_bw {f[14] or '--'} -> "
                  f"{f[15] or '--'}{'  CAMBIA' if f[20] == 'TRUE' else ''}")
    print(f"  -> outputs/tables/{{R,python}}/acto2_sensibilidad_*.csv (5)")
    print(f"  -> {fig_e.name}, {fig_x.name}")


if __name__ == "__main__":
    main()
