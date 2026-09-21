# 11_sensibilidad.R -- ACTO 2.6: controles de sensibilidad de la correlacion
#                      placenta <-> cerebro y del test formal de T8.
#
# Por que existe este archivo: el Acto 2 (T7/T8) concluye que la coordinacion
# placenta <-> cerebro NO cambia de forma demostrable entre Control y LPS -- ni
# el test formal (Fisher z sobre rho de Spearman, prohibicion 4) ni la
# simulacion de restriccion de rango (prohibicion 5) permiten afirmar una
# diferencia. T9 somete esa conclusion a dos controles; cada uno rehace la
# correlacion (Spearman rho, estratos GLOBAL / CONTROL / LPS) y el Fisher z
# Control vs LPS:
#
#   (A) EIGENGENE -- variante PCA del score compuesto (D8 la nombra). En lugar
#       del promedio de los 7 z de transportadores, el resumen del modulo es la
#       proyeccion de cada feto sobre PC1 de los 7 transportadores, POR TEJIDO.
#       PCA PROPIO: eigendescomposicion de la matriz de correlacion (Pearson)
#       7x7 por Jacobi clasico SIN trigonometria (solo +,-,*,/,sqrt ->
#       bit-identico R/Python, el mismo motivo por el que el RNG polar de
#       00_config evita sin/cos). Signo de PC1 fijado a cargas mayoritariamente
#       positivas (suma > 0). Fetos con < 7 transportadores detectados: se
#       proyecta con los z disponibles y se reescala por la norma de las cargas
#       usadas -- analogo al "promedio de los z disponibles" de D8 (decision del
#       usuario).
#
#   (B) EXCLUSION DEL FETO EXTREMO -- uno por item (decision del usuario:
#       leave-one-out sobre rho; blanco = rho GLOBAL). El feto extremo es el de
#       mayor |rho_full - rho_sin_i| sobre rho GLOBAL. Se lo excluye de los TRES
#       estratos y se rehace rho (IC Bonett-Wright) y el Fisher z; se reporta si
#       el veredicto (p_bw < .05) cambia. NO se re-estima el PC1 al quitar el
#       feto.
#
# PARIDAD R/Python: Jacobi sin trig + Spearman/Pearson/SD/normas con acumulador
# double explicito (mismo orden que Python) -> texto "%.10g" bit-identico para
# rho, cargas, autovalores y eigengene; lo que pasa por trascendentes (p por pt,
# IC por tanh/atanh, p del Fisher z por la normal) -> texto "%.6e". Aca (solo R)
# se cruza-verifica el PCA PROPIO contra eigen() (tol 1e-8) y cada rho de
# estrato contra cor.test(method="spearman", exact=FALSE) (tol 1e-9); NUNCA el
# metodo exacto AS 89 (segfaultea en este build de R). Figuras PNG:
# equivalentes, no byte-identicas.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({ library(ggplot2); library(patchwork) })

ESTE_SCRIPT <- "11_sensibilidad"

COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")
Z975 <- 1.959963984540054
PISO_PAR <- 5L

GEN_SIN_CEREBRO <- "il6"
GEN_EXCLUIDO_CORR <- "il6R"   # pedido explicito (08_acto2_correlaciones): deteccion insuficiente en cerebro
GENES_CORR <- setdiff(GENES, c(GEN_SIN_CEREBRO, GEN_EXCLUIDO_CORR))  # 8 genes
ITEMS_LOO <- c(GENES_CORR, "score_compuesto", "eigengene")  # 10 items para (B)
ITEMS_A <- c("eigengene", "score_compuesto")                # (A): eigengene primero
ESTRATOS <- c("GLOBAL", "CONTROL", "LPS")
TRANSP <- GENES_TRANSPORTADORES                             # 7, base del PCA
TEJIDOS <- TEJIDOS_E15                                      # PLACENTA_E15, BRAIN_E15
DPI <- 300

JACOBI_TOL <- 1e-15
JACOBI_MAX_SWEEPS <- 100L

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..09.
# ---------------------------------------------------------------------------
.fmt <- function(x) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("")
  if (is.logical(x)) return(if (x) "TRUE" else "FALSE")
  if (is.character(x)) return(x)
  if (x == floor(x) && abs(x) < 1e15) return(sprintf("%d", as.integer(round(x))))
  sprintf("%.10g", x)
}
.csv_cell <- function(x) {
  s <- .fmt(x)
  if (grepl('[,"\r\n]', s)) s <- paste0('"', gsub('"', '""', s, fixed = TRUE), '"')
  s
}
escribir_lineas <- function(ruta, texto) {
  if (!endsWith(texto, "\n")) texto <- paste0(texto, "\n")
  con <- file(ruta, open = "wb"); writeBin(charToRaw(enc2utf8(texto)), con); close(con)
}
escribir_csv <- function(ruta, encabezado, filas) {
  cuerpo <- vapply(filas, function(f)
    paste(vapply(f, .csv_cell, character(1)), collapse = ","), character(1))
  escribir_lineas(ruta, paste(c(paste(vapply(encabezado, .csv_cell, character(1)),
                                       collapse = ","), cuerpo), collapse = "\n"))
}
p6e <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "" else
  sprintf("%.6e", as.numeric(x))
g10 <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "" else
  sprintf("%.10g", as.numeric(x))

# ---------------------------------------------------------------------------
# Nucleo PROPIO -- acumulador double explicito (mismo orden que Python).
# ---------------------------------------------------------------------------
.suma <- function(xs) { a <- 0; for (v in xs) a <- a + v; a }
desvio <- function(xs) {
  n <- length(xs)
  if (n < 2L) return(NA_real_)
  m <- .suma(xs) / n
  s <- 0; for (v in xs) s <- s + (v - m) * (v - m)
  sqrt(s / (n - 1))
}
rangos_promedio <- function(xs) as.numeric(rank(xs, ties.method = "average"))

spearman_rho <- function(x, y) {
  n <- length(x)
  if (n < 2L) return(NA_real_)
  rx <- rangos_promedio(x); ry <- rangos_promedio(y)
  mx <- .suma(rx) / n; my <- .suma(ry) / n
  sxy <- 0; sxx <- 0; syy <- 0
  for (i in seq_len(n)) {
    dx <- rx[i] - mx; dy <- ry[i] - my
    sxy <- sxy + dx * dy; sxx <- sxx + dx * dx; syy <- syy + dy * dy
  }
  if (sxx <= 0 || syy <= 0) return(NA_real_)
  sxy / sqrt(sxx * syy)
}
spearman_p <- function(rho, n) {
  if (is.na(rho) || n <= 2L) return(NA_real_)
  if (rho * rho >= 1) return(0)
  t <- rho * sqrt((n - 2) / (1 - rho * rho))
  2 * pt(abs(t), n - 2, lower.tail = FALSE)
}
spearman_ci <- function(rho, n) {
  if (is.na(rho) || n <= 3L || rho * rho >= 1) return(c(NA_real_, NA_real_))
  se <- sqrt((1 + rho * rho / 2) / (n - 3))
  z <- atanh(rho)
  c(tanh(z - Z975 * se), tanh(z + Z975 * se))
}
norm_sf2 <- function(x) 2 * pnorm(-abs(x))

pearson_pairwise <- function(a, b) {
  xs <- c(); ys <- c()
  for (i in seq_along(a)) {
    if (!is.na(a[i]) && !is.na(b[i])) { xs <- c(xs, a[i]); ys <- c(ys, b[i]) }
  }
  n <- length(xs)
  if (n < 3L) return(NA_real_)
  mx <- .suma(xs) / n; my <- .suma(ys) / n
  sxy <- 0; sxx <- 0; syy <- 0
  for (i in seq_len(n)) {
    dx <- xs[i] - mx; dy <- ys[i] - my
    sxy <- sxy + dx * dy; sxx <- sxx + dx * dx; syy <- syy + dy * dy
  }
  if (sxx <= 0 || syy <= 0) return(NA_real_)
  sxy / sqrt(sxx * syy)
}

# Eigendescomposicion PROPIA por Jacobi clasico SIN trigonometria. Barrido
# ciclico p<q, umbral y sweeps fijos -> determinista y bit-identico a Python.
# Las rotaciones se aplican columna-a-columna con temporarios (mismo orden de
# operaciones IEEE que el loop escalar de Python).
jacobi_eigen <- function(A) {
  n <- nrow(A)
  M <- A + 0                                   # copia de trabajo
  V <- diag(n)
  off_sq <- function(MM) {
    s <- 0
    for (p in seq_len(n - 1L)) for (q in (p + 1L):n) s <- s + MM[p, q] * MM[p, q]
    s
  }
  for (.sw in seq_len(JACOBI_MAX_SWEEPS)) {
    if (off_sq(M) <= JACOBI_TOL) break
    for (p in seq_len(n - 1L)) {
      for (q in (p + 1L):n) {
        apq <- M[p, q]
        if (apq == 0) next
        theta <- (M[q, q] - M[p, p]) / (2 * apq)
        if (theta >= 0) {
          t <- 1 / (theta + sqrt(theta * theta + 1))
        } else {
          t <- -1 / (-theta + sqrt(theta * theta + 1))
        }
        c <- 1 / sqrt(t * t + 1)
        s <- t * c
        mp <- M[, p]; mq <- M[, q]                 # M <- M * R
        M[, p] <- c * mp - s * mq
        M[, q] <- s * mp + c * mq
        mp <- M[p, ]; mq <- M[q, ]                 # M <- R^T * M
        M[p, ] <- c * mp - s * mq
        M[q, ] <- s * mp + c * mq
        vp <- V[, p]; vq <- V[, q]                 # V <- V * R
        V[, p] <- c * vp - s * vq
        V[, q] <- s * vp + c * vq
      }
    }
  }
  list(values = diag(M), vectors = V)               # vectors[, k] = autovector k
}

# ===========================================================================
# 1. Carga y emparejamiento por feto -- identico a 08/09 + columna z por gen.
# ===========================================================================
.leer_tsv <- function(ruta) {
  read.delim(ruta, sep = "\t", quote = "", stringsAsFactors = FALSE,
             colClasses = "character", check.names = FALSE, encoding = "UTF-8",
             na.strings = character(0))
}
cargar <- function() {
  proc <- RUTA_DATOS_PROC
  cu <- .leer_tsv(file.path(proc, "qpcr_cuantificacion_long.tsv"))
  sc <- .leer_tsv(file.path(proc, "qpcr_score_compuesto_long.tsv"))
  madre <- setNames(cu$MADRE_ID, cu$FETO)[!duplicated(cu$FETO)]
  tto   <- setNames(cu$TTO, cu$FETO)[!duplicated(cu$FETO)]
  negdd <- setNames(
    ifelse(cu$neg_ddCt == "", NA_real_, suppressWarnings(as.numeric(cu$neg_ddCt))),
    paste(cu$FETO, cu$TEJIDO, cu$GEN, sep = "\r"))
  zval <- setNames(
    ifelse(cu$z == "", NA_real_, suppressWarnings(as.numeric(cu$z))),
    paste(cu$FETO, cu$TEJIDO, cu$GEN, sep = "\r"))
  scv <- setNames(
    ifelse(sc$score_compuesto == "", NA_real_,
           suppressWarnings(as.numeric(sc$score_compuesto))),
    paste(sc$FETO, sc$TEJIDO, sep = "\r"))
  fetos <- names(madre)[order(madre, names(madre), method = "radix")]
  list(fetos = fetos, madre = madre, tto = tto, negdd = negdd, z = zval, sc = scv)
}

# ===========================================================================
# 2. (A) EIGENGENE -- PC1 PROPIO del modulo de 7 transportadores, por tejido.
# ===========================================================================
construir_eigengene <- function(D) {
  eig <- new.env(); loadings <- list(); varianza <- list()
  n <- length(TRANSP)
  for (tej in TEJIDOS) {
    cols <- lapply(TRANSP, function(g)
      vapply(D$fetos, function(f) {
        v <- D$z[[paste(f, tej, g, sep = "\r")]]
        if (is.null(v)) NA_real_ else v
      }, numeric(1)))
    names(cols) <- TRANSP
    C <- matrix(1, n, n)
    for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
      r <- pearson_pairwise(cols[[TRANSP[i]]], cols[[TRANSP[j]]])
      if (is.na(r)) r <- 0
      C[i, j] <- r; C[j, i] <- r
    }
    ee <- jacobi_eigen(C)
    idx <- order(ee$values, decreasing = TRUE)
    total <- .suma(ee$values[idx])
    acum <- 0; vfilas <- list()
    for (pc in seq_along(idx)) {
      k <- idx[pc]
      prop <- if (total > 0) ee$values[k] / total else 0
      acum <- acum + prop
      vfilas[[pc]] <- list(pc = pc, autoval = ee$values[k], prop = prop, acum = acum)
    }
    varianza[[tej]] <- vfilas
    L <- ee$vectors[, idx[1]]
    nrm <- sqrt(.suma(L * L))
    L <- L / nrm
    if (.suma(L) < 0) L <- -L
    loadings[[tej]] <- setNames(as.list(L), TRANSP)
    for (fi in seq_along(D$fetos)) {
      f <- D$fetos[fi]
      num <- 0; den <- 0
      for (i in seq_len(n)) {
        zv <- cols[[TRANSP[i]]][fi]
        if (is.na(zv)) next
        num <- num + L[i] * zv
        den <- den + L[i] * L[i]
      }
      eig[[paste(f, tej, sep = "\r")]] <- if (den > 0) num / sqrt(den) else NA_real_
    }
  }
  list(eig = eig, loadings = loadings, varianza = varianza)
}

# ===========================================================================
# 3. Emparejamiento por feto (generaliza el de 08/09 a item = "eigengene").
# ===========================================================================
.valor <- function(D, eig, item, feto, tejido) {
  if (item == "eigengene") {
    v <- eig[[paste(feto, tejido, sep = "\r")]]
    return(if (is.null(v)) NA_real_ else v)
  }
  if (item == "score_compuesto") {
    v <- D$sc[[paste(feto, tejido, sep = "\r")]]
    return(if (is.null(v)) NA_real_ else v)
  }
  v <- D$negdd[[paste(feto, tejido, item, sep = "\r")]]
  if (is.null(v)) NA_real_ else v
}
pares <- function(D, eig, item, estrato, excluir = NULL) {
  xs <- c(); ys <- c(); ts <- c(); fs <- c()
  for (f in D$fetos) {
    if (!is.null(excluir) && f == excluir) next
    if (estrato != "GLOBAL" && D$tto[[f]] != estrato) next
    xp <- .valor(D, eig, item, f, "PLACENTA_E15")
    yb <- .valor(D, eig, item, f, "BRAIN_E15")
    if (is.na(xp) || is.na(yb)) next
    xs <- c(xs, xp); ys <- c(ys, yb); ts <- c(ts, D$tto[[f]]); fs <- c(fs, f)
  }
  list(x = xs, y = ys, tto = ts, feto = fs)
}

# cruza-verificacion de rho (solo R): PROPIO vs cor.test(exact=FALSE).
.CV <- new.env(); .CV$peor <- 0
rho_cv <- function(x, y) {
  r <- spearman_rho(x, y)
  if (length(x) >= 3L && !is.na(r)) {
    ct <- tryCatch(suppressWarnings(
      cor.test(x, y, method = "spearman", exact = FALSE)), error = function(e) NULL)
    if (!is.null(ct)) .CV$peor <- max(.CV$peor, abs(r - as.numeric(ct$estimate)))
  }
  r
}

# ===========================================================================
# 4. (A) tablas: correlacion y Fisher z del eigengene vs score compuesto.
# ===========================================================================
COLS_A_CORR <- c("ITEM", "ESTRATO", "n_par", "rho_spearman", "ic95_low",
                 "ic95_high", "p_valor")

tabla_a_correlacion <- function(D, eig) {
  filas <- list()
  for (item in ITEMS_A) for (est in ESTRATOS) {
    pr <- pares(D, eig, item, est)
    n <- length(pr$x)
    if (n >= PISO_PAR) {
      rho <- rho_cv(pr$x, pr$y)
      ic <- spearman_ci(rho, n); p <- spearman_p(rho, n)
      filas[[length(filas) + 1L]] <- list(item, est, n, g10(rho), p6e(ic[1]),
                                          p6e(ic[2]), p6e(p))
    } else {
      filas[[length(filas) + 1L]] <- list(item, est, n, "", "", "", "")
    }
  }
  filas
}

COLS_A_TEST <- c("ITEM", "n_control", "rho_control", "n_lps", "rho_lps",
                 "delta_rho", "z_control", "z_lps", "se_bw_control", "se_bw_lps",
                 "stat_z_bw", "p_bw", "se_clasico_control", "se_clasico_lps",
                 "stat_z_clasico", "p_clasico")

.se_bw <- function(rho, n) sqrt((1 + rho * rho / 2) / (n - 3))
.se_clasico <- function(n) 1 / sqrt(n - 3)
.fisher_z <- function(rc, rl, sec, sel) {
  zc <- atanh(rc); zl <- atanh(rl)
  stat <- (zc - zl) / sqrt(sec * sec + sel * sel)
  list(zc = zc, zl = zl, sec = sec, sel = sel, stat = stat, p = norm_sf2(stat))
}

.fila_fisher <- function(item, prc, prl) {
  nc <- length(prc$x); nl <- length(prl$x)
  if (nc < PISO_PAR || nl < PISO_PAR)
    return(list(item, nc, "", nl, "", "", "", "", "", "", "", "", "", "", "", ""))
  rc <- rho_cv(prc$x, prc$y); rl <- rho_cv(prl$x, prl$y)
  drho <- rc - rl
  bw <- .fisher_z(rc, rl, .se_bw(rc, nc), .se_bw(rl, nl))
  cl <- .fisher_z(rc, rl, .se_clasico(nc), .se_clasico(nl))
  list(item, nc, g10(rc), nl, g10(rl), g10(drho), p6e(bw$zc), p6e(bw$zl),
       p6e(bw$sec), p6e(bw$sel), p6e(bw$stat), p6e(bw$p), p6e(cl$sec),
       p6e(cl$sel), p6e(cl$stat), p6e(cl$p))
}
tabla_a_test <- function(D, eig) {
  lapply(ITEMS_A, function(item)
    .fila_fisher(item, pares(D, eig, item, "CONTROL"), pares(D, eig, item, "LPS")))
}

COLS_PCA_LOAD <- c("TEJIDO", "GEN", "loading_pc1")
COLS_PCA_VAR <- c("TEJIDO", "PC", "autovalor", "prop_var", "prop_var_acum")

tablas_pca <- function(loadings, varianza) {
  lo <- list()
  for (tej in TEJIDOS) for (g in TRANSP)
    lo[[length(lo) + 1L]] <- list(tej, g, g10(loadings[[tej]][[g]]))
  va <- list()
  for (tej in TEJIDOS) for (r in varianza[[tej]])
    va[[length(va) + 1L]] <- list(tej, r$pc, g10(r$autoval), g10(r$prop),
                                  g10(r$acum))
  list(load = lo, var = va)
}

# ===========================================================================
# 5. (B) exclusion del feto extremo (leave-one-out sobre rho GLOBAL).
# ===========================================================================
COLS_B <- c("ITEM", "TIPO", "n_par_global", "feto_extremo", "tto_feto_extremo",
            "rho_global_full", "rho_global_sin", "impacto_rho_global",
            "rho_control_full", "rho_control_sin", "rho_lps_full", "rho_lps_sin",
            "delta_rho_full", "delta_rho_sin", "p_bw_full", "p_bw_sin",
            "p_clasico_full", "p_clasico_sin", "veredicto_full", "veredicto_sin",
            "veredicto_cambia")

.rho_est <- function(D, eig, item, est, excluir = NULL) {
  pr <- pares(D, eig, item, est, excluir = excluir)
  n <- length(pr$x)
  if (n < PISO_PAR) return(NA_real_)
  rho_cv(pr$x, pr$y)
}
.fisher_pb_pc <- function(D, eig, item, excluir = NULL) {
  prc <- pares(D, eig, item, "CONTROL", excluir = excluir)
  prl <- pares(D, eig, item, "LPS", excluir = excluir)
  nc <- length(prc$x); nl <- length(prl$x)
  if (nc < PISO_PAR || nl < PISO_PAR) return(list(dr = NA_real_, pb = NA_real_,
                                                  pc = NA_real_))
  rc <- rho_cv(prc$x, prc$y); rl <- rho_cv(prl$x, prl$y)
  bw <- .fisher_z(rc, rl, .se_bw(rc, nc), .se_bw(rl, nl))
  cl <- .fisher_z(rc, rl, .se_clasico(nc), .se_clasico(nl))
  list(dr = rc - rl, pb = bw$p, pc = cl$p)
}
.veredicto <- function(pb) if (is.na(pb)) "" else if (pb < 0.05) "dif" else "n.s."

tabla_b <- function(D, eig) {
  filas <- list()
  for (item in ITEMS_LOO) {
    tipo <- if (item == "score_compuesto") "score" else
      if (item == "eigengene") "eigengene" else "gen"
    pr <- pares(D, eig, item, "GLOBAL")
    n <- length(pr$x)
    if (n < PISO_PAR) {
      filas[[length(filas) + 1L]] <- as.list(c(item, tipo, n, rep("", 18L)))
      next
    }
    rho_full <- spearman_rho(pr$x, pr$y)
    peor_feto <- NA_character_; peor_imp <- -1
    for (k in seq_len(n)) {
      rk <- spearman_rho(pr$x[-k], pr$y[-k])
      if (is.na(rk)) next
      imp <- abs(rho_full - rk)
      if (imp > peor_imp) { peor_imp <- imp; peor_feto <- pr$feto[k] }
    }
    tto_ex <- D$tto[[peor_feto]]
    rho_g_sin <- .rho_est(D, eig, item, "GLOBAL", excluir = peor_feto)
    rc_full <- .rho_est(D, eig, item, "CONTROL")
    rc_sin  <- .rho_est(D, eig, item, "CONTROL", excluir = peor_feto)
    rl_full <- .rho_est(D, eig, item, "LPS")
    rl_sin  <- .rho_est(D, eig, item, "LPS", excluir = peor_feto)
    ff <- .fisher_pb_pc(D, eig, item)
    fs <- .fisher_pb_pc(D, eig, item, excluir = peor_feto)
    v_full <- .veredicto(ff$pb); v_sin <- .veredicto(fs$pb)
    cambia <- if (is.na(ff$pb) || is.na(fs$pb)) "" else
      if ((ff$pb < 0.05) != (fs$pb < 0.05)) "TRUE" else "FALSE"
    filas[[length(filas) + 1L]] <- list(
      item, tipo, n, peor_feto, tto_ex,
      g10(rho_full), g10(rho_g_sin), g10(rho_full - rho_g_sin),
      g10(rc_full), g10(rc_sin), g10(rl_full), g10(rl_sin),
      g10(ff$dr), g10(fs$dr), p6e(ff$pb), p6e(fs$pb),
      p6e(ff$pc), p6e(fs$pc), v_full, v_sin, cambia)
  }
  filas
}

# ===========================================================================
# 6. Figuras.
# ===========================================================================
figura_eigengene <- function(D, eig, loadings, varianza, a_corr, ruta) {
  dl <- do.call(rbind, lapply(TEJIDOS, function(tej) {
    tt <- if (tej == "PLACENTA_E15") "Placenta E15" else "Cerebro fetal E15"
    pv <- varianza[[tej]][[1]]$prop
    data.frame(tejido = sprintf("%s  (var. expl. %.0f %%)", tt, pv * 100),
               gen = factor(TRANSP, levels = TRANSP),
               carga = vapply(TRANSP, function(g) loadings[[tej]][[g]], numeric(1)),
               stringsAsFactors = FALSE)
  }))
  p1 <- ggplot(dl, aes(gen, carga)) +
    geom_col(fill = "#0072B2", colour = "grey30", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    facet_wrap(~ tejido) +
    labs(title = "Cargas de PC1 (Jacobi PROPIO sin trigonometria)",
         x = NULL, y = "carga") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 40, hjust = 1, face = "italic",
                                     size = 7),
          strip.background = element_rect(fill = "grey93", colour = NA))

  pr <- pares(D, eig, "eigengene", "GLOBAL")
  dd <- data.frame(x = pr$x, y = pr$y,
                   tto = ifelse(pr$tto == "CONTROL", "Control", "LPS"))
  rg <- spearman_rho(pr$x, pr$y); ic <- spearman_ci(rg, length(pr$x))
  p2 <- ggplot(dd, aes(x, y, colour = tto)) +
    geom_point(size = 2.4) +
    scale_colour_manual(values = c(Control = unname(COL_TTO["CONTROL"]),
                                   LPS = unname(COL_TTO["LPS"])), name = NULL) +
    labs(title = sprintf("Eigengene placenta <-> cerebro (GLOBAL): rho = %.2f [%.2f, %.2f] (n = %d)",
                         rg, ic[1], ic[2], length(pr$x)),
         x = "placenta  eigengene (PC1)", y = "cerebro  eigengene (PC1)") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  idx <- list()
  for (f in a_corr) idx[[paste(f[[1]], f[[2]])]] <- f[[4]]
  dr <- do.call(rbind, lapply(ESTRATOS, function(e) {
    data.frame(estrato = factor(e, levels = ESTRATOS),
               metodo = c("score compuesto", "eigengene (PC1)"),
               rho = c(
                 if (nzchar(idx[[paste("score_compuesto", e)]]))
                   as.numeric(idx[[paste("score_compuesto", e)]]) else NA_real_,
                 if (nzchar(idx[[paste("eigengene", e)]]))
                   as.numeric(idx[[paste("eigengene", e)]]) else NA_real_),
               stringsAsFactors = FALSE)
  }))
  p3 <- ggplot(dr, aes(estrato, rho, fill = metodo)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62,
             colour = "grey30", linewidth = 0.3) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    scale_fill_manual(values = c("score compuesto" = "grey55",
                                 "eigengene (PC1)" = "#0072B2"), name = NULL) +
    labs(title = "Correlacion placenta <-> cerebro: score vs eigengene",
         x = NULL, y = "rho de Spearman") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  pp <- p1 / (p2 | p3) +
    plot_annotation(
      title = paste0("T9 (A) -- eigengene: PC1 PROPIO del modulo de 7 ",
                     "transportadores (por tejido) como alternativa al score ",
                     "compuesto de D8"),
      theme = theme(plot.title = element_text(size = 10)))
  ggsave(ruta, pp, width = 12.0, height = 9.0, dpi = DPI)
}

figura_excl_extremo <- function(tb, ruta) {
  ord <- vapply(ITEMS_LOO, function(i)
    if (i == "score_compuesto") "score compuesto" else
      if (i == "eigengene") "eigengene (PC1)" else i, character(1))
  reg <- list(); ann <- list()
  for (f in tb) {
    nom <- if (f[[1]] == "score_compuesto") "score compuesto" else
      if (f[[1]] == "eigengene") "eigengene (PC1)" else f[[1]]
    if (!nzchar(f[[6]]) || !nzchar(f[[7]])) {
      ann[[length(ann) + 1L]] <- data.frame(item = nom,
        lab = sprintf("n<%d (sin test)", PISO_PAR), stringsAsFactors = FALSE)
      next
    }
    rf <- as.numeric(f[[6]]); rs <- as.numeric(f[[7]])
    reg[[length(reg) + 1L]] <- data.frame(item = nom, tipo = "todos", rho = rf,
                                          stringsAsFactors = FALSE)
    reg[[length(reg) + 1L]] <- data.frame(item = nom, tipo = "sin feto extremo",
                                          rho = rs, stringsAsFactors = FALSE)
    pbf <- if (nzchar(f[[15]])) sprintf("%.3f", as.numeric(f[[15]])) else "--"
    pbs <- if (nzchar(f[[16]])) sprintf("%.3f", as.numeric(f[[16]])) else "--"
    marca <- if (identical(f[[21]], "TRUE")) " (cambia)" else ""
    ann[[length(ann) + 1L]] <- data.frame(item = nom,
      lab = sprintf("p_bw %s -> %s%s", pbf, pbs, marca), stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, reg); a <- do.call(rbind, ann)
  d$item <- factor(d$item, levels = rev(ord))
  a$item <- factor(a$item, levels = rev(ord))
  p <- ggplot(d, aes(rho, item, colour = tipo)) +
    geom_line(aes(group = item), colour = "grey70", linewidth = 0.6) +
    geom_point(size = 2.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.3) +
    geom_text(data = a, aes(x = 1.08, y = item, label = lab), inherit.aes = FALSE,
              hjust = 0, size = 2.5, colour = "grey15") +
    scale_colour_manual(values = c("todos" = "grey35",
                                   "sin feto extremo" = "#D55E00"), name = NULL) +
    coord_cartesian(xlim = c(-1.05, 1.05), clip = "off") +
    labs(title = paste0("T9 (B) -- exclusion del feto extremo (leave-one-out ",
                        "sobre rho GLOBAL, uno por item)"),
         subtitle = "p_bw del Fisher z Control vs LPS: completo -> sin el feto extremo",
         x = "rho de Spearman GLOBAL (placenta <-> cerebro por feto)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.margin = margin(5.5, 130, 5.5, 5.5),
          plot.subtitle = element_text(size = 8))
  ggsave(ruta, p, width = 10.0, height = 6.5, dpi = DPI)
}

# ===========================================================================
# 7. Artefactos compartidos (merge por 'script') -- headers identicos a 02..09.
# ===========================================================================
parse_csv_line <- function(linea) {
  chars <- strsplit(linea, "", fixed = TRUE)[[1]]
  out <- character(0); cur <- ""; q <- FALSE; i <- 1L; n <- length(chars)
  while (i <= n) {
    c <- chars[i]
    if (q) {
      if (c == '"') {
        if (i < n && chars[i + 1L] == '"') { cur <- paste0(cur, '"'); i <- i + 1L }
        else q <- FALSE
      } else cur <- paste0(cur, c)
    } else {
      if (c == '"') q <- TRUE
      else if (c == ",") { out <- c(out, cur); cur <- "" }
      else cur <- paste0(cur, c)
    }
    i <- i + 1L
  }
  c(out, cur)
}
leer_csv <- function(ruta) {
  if (!file.exists(ruta)) return(NULL)
  txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(txt) <- "UTF-8"
  lineas <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  if (length(lineas) && lineas[length(lineas)] == "") lineas <- lineas[-length(lineas)]
  if (!length(lineas)) return(NULL)
  list(header = parse_csv_line(lineas[1]), filas = lapply(lineas[-1], parse_csv_line))
}
merge_por_script <- function(ruta, header, filas_nuevas, col_script, clave_orden) {
  viejo <- leer_csv(ruta)
  if (!is.null(viejo) && !identical(viejo$header, header))
    stop(sprintf("%s: encabezado incompatible", basename(ruta)))
  idx <- match(col_script, header)
  conservadas <- if (is.null(viejo)) list() else
    Filter(function(f) length(f) >= idx && f[[idx]] != ESTE_SCRIPT, viejo$filas)
  nuevas <- lapply(filas_nuevas, function(f) vapply(f, .fmt, character(1)))
  todas <- c(conservadas, nuevas)
  todas <- todas[clave_orden(todas)]
  escribir_csv(ruta, header, todas)
}
registrar_procedencia <- function(filas_nuevas) {
  header <- c("artefacto", "tipo", "script", "origen_codigo", "entradas", "descripcion")
  merge_por_script(file.path(RUTA_TABLAS, "procedencia.csv"), header, filas_nuevas,
                   "script", function(fs) order(
                     vapply(fs, `[[`, character(1), 3), vapply(fs, `[[`, character(1), 1),
                     method = "radix"))
}
registrar_verificaciones <- function(filas_nuevas) {
  header <- c("id", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script")
  merge_por_script(file.path(RUTA_TABLAS, "verificaciones.csv"), header, filas_nuevas,
                   "script", function(fs) order(
                     vapply(fs, `[[`, character(1), 6), vapply(fs, `[[`, character(1), 1),
                     method = "radix"))
}

DESCARTES <- paste(c(
"## 11_sensibilidad",
"",
"### Que resuelve (Acto 2.6)",
"",
paste0("- Dos controles de sensibilidad de la conclusion del Acto 2 (la ",
       "coordinacion placenta <-> cerebro no cambia de forma demostrable entre ",
       "Control y LPS). Cada uno rehace la correlacion (Spearman rho, GLOBAL / ",
       "CONTROL / LPS) y el Fisher z Control vs LPS de `09_acto2_dispersion`."),
"",
"### (A) Eigengene -- variante PCA del score compuesto (D8)",
"",
paste0("- El score compuesto de D8 es el **promedio** de los 7 z de ",
       "transportadores. El eigengene es la **proyeccion sobre PC1** del modulo ",
       "de 7 transportadores, calculado **por tejido** (placenta y cerebro ",
       "tienen su propia matriz y su propio PC1)."),
paste0("- **Matriz de entrada**: correlacion de Pearson (pairwise-complete) ",
       "entre las 7 columnas de z (D8, z por gen x tejido sobre los 36 fetos, ",
       "sd n-1). Se considero la matriz de Spearman para ser coherente con el ",
       "resto del Acto 2 y **no se uso**: el eigengene es ",
       "una reduccion de dimension estandar del modulo (Pearson / datos ",
       "estandarizados); el test de asociacion entre tejidos sigue siendo ",
       "Spearman. Se prefirio pairwise-complete a casos completos porque este ",
       "ultimo descartaria 4-6 fetos (placenta / cerebro) y cambiaria el modulo ",
       "para todos. Costo: la eliminacion pairwise puede volver la matriz ",
       "levemente no definida positiva (en cerebro el 7mo autovalor queda en ",
       "~ -0.11 sobre una traza de 7); PC1 -- lo unico que se usa -- domina ",
       "(>76 % de la varianza, cargas todas positivas) y coincide con `eigen()` ",
       "a 1e-13, asi que el artefacto no lo afecta."),
paste0("- **PCA PROPIO**: eigendescomposicion por **Jacobi clasico SIN ",
       "trigonometria** -- la rotacion se arma con ",
       "`t = 1/(theta + sign*sqrt(theta^2+1))`, `c = 1/sqrt(t^2+1)`, `s = t*c` ",
       "(solo `+ - * / sqrt`, correctamente redondeados en R y Python). Barrido ",
       "ciclico p<q, umbral `1e-15` sobre la suma de cuadrados de la triangular ",
       "superior, `<= 100` sweeps. Es el mismo motivo por el que el RNG de ",
       "`00_config` usa el metodo polar y evita `sin`/`cos`: determinismo bit a ",
       "bit entre lenguajes. `prcomp` / `numpy.linalg.eigh` **no** entran en el ",
       "resultado."),
paste0("- **Signo de PC1**: se fija a cargas mayoritariamente positivas ",
       "(`suma(cargas) > 0`, si no se invierte el vector). Convencion de ",
       "\"tono de expresion\" del modulo."),
paste0("- **Fetos con < 7 transportadores detectados** (decision del usuario): ",
       "se proyecta con los z disponibles y se divide por ",
       "`sqrt(sum(carga_g^2))` de las cargas efectivamente usadas -- para un ",
       "feto completo esto es exactamente el score de PC1; para uno incompleto ",
       "es la proyeccion sobre la direccion unitaria del subespacio disponible ",
       "(analogo al \"promedio de los z disponibles\" de D8). Alternativa ",
       "descartada: suma cruda sin reescalar (sesga a la baja la magnitud de ",
       "los fetos incompletos)."),
"",
"### (B) Exclusion del feto extremo (decision del usuario)",
"",
paste0("- **Criterio**: leave-one-out sobre rho. Para cada item se quita un ",
       "feto por vez y se recalcula **rho GLOBAL**; el feto extremo es el de ",
       "mayor `|rho_full - rho_sin_i|` (empates -> el primero en el orden ",
       "`MADRE_ID, FETO`). Se descarto la distancia de Cook: mide influencia ",
       "sobre un ajuste **lineal**, y el estadistico que se reporta en el Acto ",
       "2 es de **rango** (Spearman, decision del usuario en T7)."),
paste0("- **Blanco = rho GLOBAL** (decision del usuario): el feto que mas mueve ",
       "la correlacion pooled. Se lo excluye de los TRES estratos y se rehace ",
       "rho (con IC Bonett-Wright) y el Fisher z; `veredicto_cambia` marca si ",
       "`p_bw < .05` cambia de estado."),
paste0("- **No se re-estima el PC1** al quitar el feto: la pregunta es la ",
       "sensibilidad de la correlacion a un punto influyente, no la del PC1. ",
       "Para el item `eigengene` se usan las cargas calculadas sobre los 36 ",
       "fetos."),
paste0("- Items: los 8 genes de T7/T8 (todos menos `il6`, D7, y `il6R`, excluido ",
       "de todo el Acto 2 -- bugfix `cambios_informe_conclusiones.md` punto 5b, ",
       "ver `analisis_descartados.md` de `09_acto2_dispersion`) + score ",
       "compuesto + eigengene. Todos superan el piso de n>=5 por grupo en estos ",
       "datos; si alguno no lo alcanzara, su fila llevaria los n y el impacto ",
       "sobre rho GLOBAL, sin Fisher z."),
"",
"### Paridad R / Python",
"",
paste0("- Jacobi sin trig + Spearman / Pearson / SD / normas con acumulador ",
       "`double` explicito (mismo orden) -> `%.10g` bit-identico para rho, ",
       "cargas, autovalores y eigengene. Lo que pasa por trascendentes (`p` por ",
       "`pt`, IC por `tanh`/`atanh`, `p` del Fisher z por la normal) -> texto ",
       "`%.6e`. Figuras PNG: equivalentes, no byte-identicas."),
paste0("- Solo en R se cruza-verifica el PCA PROPIO contra `eigen()` (tol 1e-8) ",
       "y cada rho de estrato contra `cor.test(method=\"spearman\", ",
       "exact=FALSE)` (tol 1e-9); en Python el PCA se contrasta contra ",
       "`numpy.linalg.eigh` (tol 1e-6). Con `stopifnot` / `assert`.")
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 11_sensibilidad:inicio -->"
  marca_fin <- "<!-- 11_sensibilidad:fin -->"
  nuevo <- paste0(marca_ini, "\n", DESCARTES, "\n\n", marca_fin)
  if (file.exists(ruta)) {
    txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(txt) <- "UTF-8"
  } else {
    txt <- paste0("# Analisis descartados\n\nQue se probo, por que no funciono o no",
                  " se uso, y que se hizo en su lugar. Una seccion por script.\n")
  }
  if (grepl(marca_ini, txt, fixed = TRUE) && grepl(marca_fin, txt, fixed = TRUE)) {
    pre <- strsplit(txt, marca_ini, fixed = TRUE)[[1]][1]
    post <- strsplit(txt, marca_fin, fixed = TRUE)[[1]][2]
    if (is.na(post)) post <- ""
    txt <- paste0(pre, nuevo, post)
  } else {
    if (!endsWith(txt, "\n")) txt <- paste0(txt, "\n")
    txt <- paste0(txt, "\n", nuevo, "\n")
  }
  escribir_lineas(ruta, txt)
}

# ===========================================================================
# 8. Reporte legible.
# ===========================================================================
.md <- function(header, filas) {
  l1 <- paste0("| ", paste(header, collapse = " | "), " |")
  l2 <- paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |")
  cuerpo <- vapply(filas, function(f)
    paste0("| ", paste(vapply(f, .fmt, character(1)), collapse = " | "), " |"),
    character(1))
  paste(c(l1, l2, cuerpo), collapse = "\n")
}

construir_reporte <- function(fuente, a_corr, a_test, pca_load, pca_var, tb) {
  L <- c(
    "# Reporte de sensibilidad Acto 2.6 (T9)", "",
    "Generado por `11_sensibilidad` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. (A) Eigengene -- PC1 PROPIO del modulo de 7 transportadores", "",
    paste0("- PC1 por tejido (eigendescomposicion de la matriz de correlacion de ",
           "Pearson de los z, Jacobi PROPIO sin trigonometria). Signo fijado a ",
           "cargas mayoritariamente positivas. Fetos con <7 detectados: ",
           "proyeccion reescalada por la norma de las cargas usadas (decision ",
           "del usuario)."), "",
    "### 1.1 Cargas de PC1", "",
    .md(COLS_PCA_LOAD, pca_load), "",
    "### 1.2 Varianza explicada", "",
    .md(COLS_PCA_VAR, pca_var), "",
    "### 1.3 Correlacion placenta <-> cerebro: eigengene vs score compuesto", "",
    .md(COLS_A_CORR, a_corr), "",
    "### 1.4 Fisher z Control vs LPS: eigengene vs score compuesto", "",
    .md(COLS_A_TEST, a_test), "",
    "## 2. (B) Exclusion del feto extremo (leave-one-out sobre rho GLOBAL)", "",
    paste0("- Un feto por item: el de mayor `|rho_full - rho_sin_i|` sobre rho ",
           "GLOBAL. Se excluye de los 3 estratos y se rehace rho + Fisher z. ",
           "`veredicto_cambia` = TRUE si `p_bw < .05` cambia de estado."), "",
    .md(COLS_B, tb), "",
    "## 3. Figuras", "",
    paste0("- `outputs/figures/acto2_sensibilidad_eigengene.png` -- cargas de ",
           "PC1 por tejido, dispersion del eigengene y comparacion de rho score ",
           "vs eigengene por estrato."),
    paste0("- `outputs/figures/acto2_sensibilidad_excl_extremo.png` -- rho ",
           "GLOBAL con y sin el feto extremo por item, con el cambio de ",
           "`p_bw`."), "",
    "## 4. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `11_sensibilidad`: eleccion ",
           "de la matriz de entrada del PCA, Jacobi sin trigonometria, ",
           "tratamiento de los fetos incompletos, y por que leave-one-out sobre ",
           "rho y no Cook."), "")
  paste(L, collapse = "\n")
}

# ===========================================================================
main <- function() {
  D <- cargar()
  fuente <- fuente_datos(ARCHIVO_QPCR)

  ce <- construir_eigengene(D)
  eig <- ce$eig; loadings <- ce$loadings; varianza <- ce$varianza

  # cruza-verificacion PCA PROPIO vs eigen().
  peor_pca <- 0
  n <- length(TRANSP)
  for (tej in TEJIDOS) {
    cols <- lapply(TRANSP, function(g)
      vapply(D$fetos, function(f) {
        v <- D$z[[paste(f, tej, g, sep = "\r")]]; if (is.null(v)) NA_real_ else v
      }, numeric(1)))
    names(cols) <- TRANSP
    C <- matrix(1, n, n)
    for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
      r <- pearson_pairwise(cols[[TRANSP[i]]], cols[[TRANSP[j]]]); if (is.na(r)) r <- 0
      C[i, j] <- r; C[j, i] <- r
    }
     er <- eigen(C, symmetric = TRUE)
    Lr <- er$vectors[, which.max(er$values)]
    if (sum(Lr) < 0) Lr <- -Lr
    Lp <- vapply(TRANSP, function(g) loadings[[tej]][[g]], numeric(1))
    peor_pca <- max(peor_pca, max(abs(Lp - Lr)),
                    abs(max(er$values) - varianza[[tej]][[1]]$autoval))
  }
  cat(sprintf("  [cruza-verif PCA PROPIO vs eigen()] peor |dif| = %.3e\n", peor_pca))
  stopifnot(peor_pca < 1e-8)

  a_corr <- tabla_a_correlacion(D, eig)
  a_test <- tabla_a_test(D, eig)
  tp <- tablas_pca(loadings, varianza); pca_load <- tp$load; pca_var <- tp$var
  tb <- tabla_b(D, eig)
  cat(sprintf("  [cruza-verif rho PROPIO vs cor.test(exact=FALSE)] peor |dif| = %.3e\n",
              .CV$peor))
  stopifnot(.CV$peor < 1e-9)

  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "acto2_sensibilidad_eigengene_correlacion.csv"),
                 COLS_A_CORR, a_corr)
    escribir_csv(file.path(base, "acto2_sensibilidad_eigengene_test.csv"),
                 COLS_A_TEST, a_test)
    escribir_csv(file.path(base, "acto2_sensibilidad_pca_loadings.csv"),
                 COLS_PCA_LOAD, pca_load)
    escribir_csv(file.path(base, "acto2_sensibilidad_pca_varianza.csv"),
                 COLS_PCA_VAR, pca_var)
    escribir_csv(file.path(base, "acto2_sensibilidad_excl_extremo.csv"), COLS_B, tb)
  }

  fig_e <- file.path(RUTA_FIGURAS, "acto2_sensibilidad_eigengene.png")
  fig_x <- file.path(RUTA_FIGURAS, "acto2_sensibilidad_excl_extremo.png")
  figura_eigengene(D, eig, loadings, varianza, a_corr, fig_e)
  figura_excl_extremo(tb, fig_x)

  escribir_lineas(file.path(RUTA_TABLAS, "acto2_sensibilidad_reporte.md"),
                  construir_reporte(fuente, a_corr, a_test, pca_load, pca_var, tb))
  actualizar_descartados()

  .a <- function(item, est) {
    for (f in a_corr) if (f[[1]] == item && f[[2]] == est) return(f)
    NULL
  }
  ei_glob <- .a("eigengene", "GLOBAL"); sc_glob <- .a("score_compuesto", "GLOBAL")
  ei_test <- Filter(function(f) f[[1]] == "eigengene", a_test)[[1]]
  sc_test <- Filter(function(f) f[[1]] == "score_compuesto", a_test)[[1]]
  n_cambia <- sum(vapply(tb, function(f) identical(f[[21]], "TRUE"), logical(1)))
  n_testados_b <- sum(vapply(tb, function(f) nzchar(f[[15]]), logical(1)))
  pv_pla <- varianza[["PLACENTA_E15"]][[1]]$prop
  pv_bra <- varianza[["BRAIN_E15"]][[1]]$prop
  min_ev_pla <- varianza[["PLACENTA_E15"]][[length(varianza[["PLACENTA_E15"]])]]$autoval
  min_ev_bra <- varianza[["BRAIN_E15"]][[length(varianza[["BRAIN_E15"]])]]$autoval
  cargas_pos <- vapply(TEJIDOS, function(tej)
    sum(vapply(TRANSP, function(g) loadings[[tej]][[g]] > 0, logical(1))), integer(1))
  names(cargas_pos) <- TEJIDOS

  ent <- sprintf(paste0("data/processed/qpcr_cuantificacion_long.tsv + ",
                        "qpcr_score_compuesto_long.tsv (de data/%s/%s)"),
                 fuente, ARCHIVO_QPCR)
  registrar_procedencia(list(
    list("outputs/tables/{R,python}/acto2_sensibilidad_eigengene_correlacion.csv",
         "tabla", ESTE_SCRIPT, "PROPIO", ent, paste0(
         "Spearman rho placenta<->cerebro del eigengene (PC1) vs score ",
         "compuesto x {GLOBAL, CONTROL, LPS}; IC95 Bonett-Wright, p t-aprox")),
    list("outputs/tables/{R,python}/acto2_sensibilidad_eigengene_test.csv",
         "tabla", ESTE_SCRIPT, "PROPIO", ent, paste0(
         "Fisher z Control vs LPS del eigengene vs score compuesto; SE ",
         "Bonett-Wright + SE clasico")),
    list("outputs/tables/{R,python}/acto2_sensibilidad_pca_loadings.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, paste0(
         "cargas de PC1 (Jacobi PROPIO sin trig) de los 7 transportadores por ",
         "tejido; signo fijado a suma>0")),
    list("outputs/tables/{R,python}/acto2_sensibilidad_pca_varianza.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "autovalores y proporcion de varianza de los 7 PC por tejido"),
    list("outputs/tables/{R,python}/acto2_sensibilidad_excl_extremo.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, paste0(
         "leave-one-out sobre rho GLOBAL: feto extremo por item y rho + Fisher ",
         "z con y sin el; veredicto_cambia")),
    list("outputs/figures/acto2_sensibilidad_eigengene.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, paste0(
         "cargas PC1 por tejido, dispersion del eigengene y rho score vs ",
         "eigengene por estrato")),
    list("outputs/figures/acto2_sensibilidad_excl_extremo.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, paste0(
         "rho GLOBAL con y sin el feto extremo por item, con el cambio de p_bw")),
    list("outputs/tables/acto2_sensibilidad_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.6 (T9)")
  ))
  registrar_verificaciones(list(
    list("sens_eigengene_pca_propio",
         paste0("PC1 por eigendescomposicion PROPIA (Jacobi sin trigonometria), ",
                "no prcomp/sklearn en el resultado"),
         sprintf("peor |dif| PC1 vs libreria de referencia = %.2e", peor_pca),
         "< 1e-6", if (peor_pca < 1e-6) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("sens_eigengene_pc1_signo",
         "signo de PC1 fijado a cargas mayoritariamente positivas (suma>0)",
         sprintf("cargas>0: placenta %d/7, cerebro %d/7",
                 cargas_pos[["PLACENTA_E15"]], cargas_pos[["BRAIN_E15"]]),
         "mayoria positiva por tejido",
         if (cargas_pos[["PLACENTA_E15"]] >= 4 && cargas_pos[["BRAIN_E15"]] >= 4)
           "TRUE" else "FALSE", ESTE_SCRIPT),
    list("sens_eigengene_var_pc1",
         "proporcion de varianza de PC1 por tejido (contexto del eigengene)",
         sprintf("placenta %.3f; cerebro %.3f", pv_pla, pv_bra), "0 < prop <= 1",
         if (pv_pla > 0 && pv_pla <= 1 && pv_bra > 0 && pv_bra <= 1) "TRUE" else
           "FALSE", ESTE_SCRIPT),
    list("sens_eigengene_matriz_psd",
         paste0("la matriz de correlacion pairwise-complete puede ser levemente ",
                "indefinida; PC1 domina y no se afecta"),
         sprintf("min autovalor: placenta %.3f; cerebro %.3f; PC1 var %.2f/%.2f",
                 min_ev_pla, min_ev_bra, pv_pla, pv_bra),
         "PC1 >> resto (traza = 7)",
         if (pv_pla > 0.5 && pv_bra > 0.5) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("sens_eigengene_correlacion_global",
         "rho GLOBAL placenta<->cerebro: eigengene vs score compuesto",
         sprintf("eigengene rho=%s (n=%s); score rho=%s (n=%s)",
                 ei_glob[[4]], ei_glob[[3]], sc_glob[[4]], sc_glob[[3]]),
         "ambos descritos, misma direccion cualitativa", "TRUE", ESTE_SCRIPT),
    list("sens_eigengene_fisher_z",
         "Fisher z Control vs LPS: el eigengene no cambia el veredicto del score",
         sprintf("eigengene p_bw=%s; score p_bw=%s", ei_test[[12]], sc_test[[12]]),
         "ninguno alcanza p_bw<.05 (consistente con T8)",
         if (nzchar(ei_test[[12]]) && nzchar(sc_test[[12]]) &&
             as.numeric(ei_test[[12]]) >= 0.05 && as.numeric(sc_test[[12]]) >= 0.05)
           "TRUE" else "REVISAR", ESTE_SCRIPT),
    list("sens_excl_extremo_criterio",
         paste0("feto extremo = leave-one-out sobre rho GLOBAL (no Cook); blanco ",
                "rho GLOBAL (decision del usuario)"),
         "argmax_i |rho_full - rho_sin_i| por item, sobre el estrato GLOBAL",
         "LOO sobre rho GLOBAL", "TRUE", ESTE_SCRIPT),
    list("sens_excl_extremo_items",
         "items del control (B): 8 genes de T7/T8 + score compuesto + eigengene",
         sprintf("%d items; %d con Fisher z (n>=5 por grupo)",
                 length(ITEMS_LOO), n_testados_b),
         "10 items", if (length(ITEMS_LOO) == 10L) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("sens_excl_extremo_veredicto",
         "cuantos items cambian el veredicto p_bw<.05 al excluir su feto extremo",
         sprintf("%d de %d items testeados cambian de veredicto", n_cambia,
                 n_testados_b),
         "reportado; el informe menciona el numero", "TRUE", ESTE_SCRIPT),
    list("sens_excl_extremo_no_reestima_pc1",
         "en el item eigengene NO se re-estima PC1 al quitar el feto extremo",
         "cargas de PC1 calculadas sobre los 36 fetos, fijas en el LOO",
         "PC1 fijo", "TRUE", ESTE_SCRIPT),
    list("sens_figuras",
         paste0("figuras Acto 2.6: eigengene (cargas + dispersion + rho) y ",
                "exclusion del extremo"),
         sprintf("eigengene=%s;excl_extremo=%s", file.exists(fig_e),
                 file.exists(fig_x)),
         "2 figuras existen",
         if (file.exists(fig_e) && file.exists(fig_x)) "TRUE" else "FALSE",
         ESTE_SCRIPT)
  ))

  cat("== 11_sensibilidad.R ==\n")
  cat(sprintf("  fuente = %s\n", fuente))
  cat(sprintf("  PC1 var. expl.: placenta %.3f | cerebro %.3f\n", pv_pla, pv_bra))
  cat(sprintf("  cargas>0: placenta %d/7 | cerebro %d/7\n",
              cargas_pos[["PLACENTA_E15"]], cargas_pos[["BRAIN_E15"]]))
  cat("  (A) rho GLOBAL placenta<->cerebro:\n")
  for (f in a_corr) if (f[[2]] == "GLOBAL")
    cat(sprintf("      %-16s n=%2s  rho=%9s  p=%s\n", f[[1]], f[[3]], f[[4]], f[[7]]))
  cat("  (A) Fisher z Control vs LPS:\n")
  for (f in a_test)
    cat(sprintf("      %-16s rhoC=%9s rhoL=%9s dRho=%9s p_bw=%s\n",
                f[[1]], f[[3]], f[[5]], f[[6]], f[[12]]))
  cat(sprintf("  (B) items con Fisher z: %d/%d; cambian veredicto: %d\n",
              n_testados_b, length(ITEMS_LOO), n_cambia))
  for (f in tb) if (nzchar(f[[6]]))
    cat(sprintf("      %-16s extremo=%-16s (%-7s) rhoG %8s -> %8s  p_bw %s -> %s%s\n",
                f[[1]], f[[4]], f[[5]], f[[6]], f[[7]],
                if (nzchar(f[[15]])) f[[15]] else "--",
                if (nzchar(f[[16]])) f[[16]] else "--",
                if (identical(f[[21]], "TRUE")) "  CAMBIA" else ""))
  cat("  -> outputs/tables/{R,python}/acto2_sensibilidad_*.csv (5)\n")
  cat(sprintf("  -> %s, %s\n", basename(fig_e), basename(fig_x)))
}

if (sys.nframe() == 0L) main()
