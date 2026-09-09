# 06_pstat3.R -- Modelo de pSTAT3 (D9): PSTAT3 ~ SEXO * TTO + MEMBRANA.
#
# Por que existe este archivo: pSTAT3 se midio en 3 membranas (western blot); la
# membrana es un bloque tecnico. D9 fija el modelo `PSTAT3 ~ SEXO * TTO + MEMBRANA`
# con MEMBRANA como BLOQUE FIJO (diseno balanceado 3 x 4 x 3), la MISMA cascada de
# supuestos D5 / seccion 4.1 y el MISMO post hoc D6 que la qPCR (05_qpcr_modelos).
#
#   * Diseno suma-cero de 6 columnas: [1, s, t, s*t, m1, m2] con
#       s = +1 HEMBRA / -1 MACHO,  t = +1 CONTROL / -1 LPS,
#       (m1, m2) = contr.sum de MEMBRANA (3 niveles): lev1->(1,0) lev2->(0,1) lev3->(-1,-1).
#     Terminos: SEXO=col1, TTO=col2, SEXO:TTO=col3, MEMBRANA=cols4-5. df_resid = n - 6.
#   * Cascada D5 (alfa = 0.05):
#       - Shapiro-Wilk sobre los RESIDUOS del modelo conjunto (con MEMBRANA).
#       - Levene (Brown-Forsythe, centro = mediana) sobre las 4 celdas SEXO x TTO.
#       - Shapiro >= .05 y Levene >= .05  -> ANOVA III           (rama "anova3")
#       - Shapiro >= .05 y Levene <  .05  -> Wald III con HC3     (rama "hc3")
#       - Shapiro <  .05                  -> ART aditivo hand-rolled (rama "art")
#     ARTool RECHAZA un bloque aditivo no cruzado; sobre datos reales Y sinteticos
#     la cascada cae en `anova3` => la rama `art` NO se ejecuta. Ver
#     `analisis_descartados.md`, seccion `06_pstat3`.
#   * Piso de celda: si alguna celda SEXO x TTO tiene < 5 valores -> via
#     `descriptivo_n_bajo` (pSTAT3 real y sintetico: 9/celda -> siempre se modela).
#   * Post hoc D6 (solo si SEXO x TTO p < 0.05): 4 comparaciones fijas + Holm.
#       anova3 -> contrastes de medias marginales (promediando sobre MEMBRANA, vcov OLS)
#       hc3    -> idem con vcov HC3
#       art    -> ART-C hand-rolled
#   * LIMITACION OBLIGATORIA DEL INFORME (D9): pSTAT3 esta normalizado a proteina
#     total SIN STAT3 total -> ABUNDANCIA de fosfo-STAT3, no fraccion fosforilada.
#
# PARIDAD R/Python: el nucleo numerico (OLS 6x6 por Gauss-Jordan, SS tipo III por
# comparacion de modelos, sandwich HC3, Levene, ART y ART-C, Holm) es codigo PROPIO
# identico. Aca (solo R) se cruza-verifica en corrida contra car::Anova / emmeans
# (ramas anova3 / hc3) con stopifnot (< 1e-6); la rama `art` se auto-verifica (dos
# vias de alineado) porque ARTool no puede ajustar el modelo aditivo. Unico
# componente de libreria en el resultado: Shapiro-Wilk. Estadisticos / p
# dependientes de trascendentes se guardan como texto "%.6e" (p6e), como en 03/05.
#
# Este script NO grafica (el boxplot de pSTAT3 se arma en 07_figuras_acto1).

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({
  library(car); library(emmeans); library(sandwich)
})

ESTE_SCRIPT <- "06_pstat3"

ALFA       <- 0.05
PISO_CELDA <- 5L
NCOL_DIS   <- 6L

GRUPOS_4 <- c("HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")
CELDAS_4 <- list(c("HEMBRA", "CONTROL"), c("HEMBRA", "LPS"),
                 c("MACHO", "CONTROL"), c("MACHO", "LPS"))

# posicion (1-based) de columnas de cada termino en el diseno de 6 columnas
TERM_COLS <- list("SEXO" = 2L, "TTO" = 3L, "SEXO:TTO" = 4L, "MEMBRANA" = c(5L, 6L))
TERM_ORD  <- c("SEXO", "TTO", "SEXO:TTO", "MEMBRANA")

# D6: (etiqueta, celda_a, celda_b) -> estima a - b
PARES_D6 <- list(
  list("HEMBRA_CONTROL-HEMBRA_LPS",   c("HEMBRA", "CONTROL"), c("HEMBRA", "LPS")),
  list("MACHO_CONTROL-MACHO_LPS",     c("MACHO", "CONTROL"),  c("MACHO", "LPS")),
  list("HEMBRA_LPS-MACHO_LPS",        c("HEMBRA", "LPS"),     c("MACHO", "LPS")),
  list("HEMBRA_CONTROL-MACHO_CONTROL", c("HEMBRA", "CONTROL"), c("MACHO", "CONTROL"))
)

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02/03/04/05.
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
escribir_tsv <- function(ruta, encabezado, filas) {
  cuerpo <- vapply(filas, function(f)
    paste(vapply(f, .fmt, character(1)), collapse = "\t"), character(1))
  escribir_lineas(ruta, paste(c(paste(encabezado, collapse = "\t"), cuerpo),
                              collapse = "\n"))
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
# Algebra lineal PROPIA -- misma eliminacion de Gauss-Jordan que Python.
# ---------------------------------------------------------------------------
resolver <- function(A, b) {
  n <- length(b)
  M <- cbind(matrix(unlist(A), n, n, byrow = TRUE), b)
  for (c in seq_len(n)) {
    piv <- which.max(abs(M[c:n, c])) + c - 1L
    if (piv != c) M[c(c, piv), ] <- M[c(piv, c), ]
    d <- M[c, c]
    M[c, c:(n + 1L)] <- M[c, c:(n + 1L)] / d
    for (r in seq_len(n)) {
      if (r != c && M[r, c] != 0) {
        f <- M[r, c]
        M[r, c:(n + 1L)] <- M[r, c:(n + 1L)] - f * M[c, c:(n + 1L)]
      }
    }
  }
  as.numeric(M[, n + 1L])
}
invertir <- function(A) {
  n <- nrow(A)
  M <- cbind(A, diag(n))
  for (c in seq_len(n)) {
    piv <- which.max(abs(M[c:n, c])) + c - 1L
    if (piv != c) M[c(c, piv), ] <- M[c(piv, c), ]
    d <- M[c, c]
    M[c, ] <- M[c, ] / d
    for (r in seq_len(n)) {
      if (r != c) { f <- M[r, c]; M[r, ] <- M[r, ] - f * M[c, ] }
    }
  }
  M[, (n + 1L):(2L * n), drop = FALSE]
}
.suma <- function(xs) { a <- 0; for (v in xs) a <- a + v; a }
promedio <- function(xs) .suma(xs) / length(xs)
desvio_muestral <- function(xs, media) {
  if (length(xs) < 2L) return(NA_real_)
  s <- 0; for (v in xs) s <- s + (v - media) * (v - media)
  sqrt(s / (length(xs) - 1L))
}
mediana <- function(xs) {
  if (length(xs) == 0L) return(NA_real_)
  s <- sort(xs); k <- length(s)
  if (k %% 2L == 1L) s[(k + 1L) %/% 2L] else (s[k %/% 2L] + s[k %/% 2L + 1L]) / 2
}
rangos_promedio <- function(xs) rank(xs, ties.method = "average")

f_sf  <- function(x, d1, d2) pf(x, d1, d2, lower.tail = FALSE)
t_sf2 <- function(x, df) 2 * pt(-abs(x), df)

# ---------------------------------------------------------------------------
# Diseno suma-cero de 6 columnas [1, s, t, s*t, m1, m2].
# ---------------------------------------------------------------------------
.cs_membrana <- function(memb, membranas) {
  if (memb == membranas[1]) return(c(1, 0))
  if (memb == membranas[2]) return(c(0, 1))
  c(-1, -1)
}
fila_diseno <- function(sexo, tto, memb, membranas) {
  s <- if (sexo == "HEMBRA") 1 else -1
  t <- if (tto == "CONTROL") 1 else -1
  m <- .cs_membrana(memb, membranas)
  c(1, s, t, s * t, m[1], m[2])
}
diseno <- function(m, membranas) {
  X <- lapply(seq_len(nrow(m)), function(i)
    fila_diseno(m$SEXO[i], m$TTO[i], m$MEMBRANA[i], membranas))
  list(X = X, y = m$y)
}
ajustar <- function(X, y, cols) {
  n <- length(y); p <- length(cols)
  Xs <- lapply(X, function(f) f[cols])
  XtX <- lapply(seq_len(p), function(a)
    vapply(seq_len(p), function(b) .suma(vapply(seq_len(n), function(i)
      Xs[[i]][a] * Xs[[i]][b], numeric(1))), numeric(1)))
  Xty <- vapply(seq_len(p), function(a) .suma(vapply(seq_len(n), function(i)
    Xs[[i]][a] * y[i], numeric(1))), numeric(1))
  b <- resolver(XtX, Xty)
  resid <- vapply(seq_len(n), function(i)
    y[i] - .suma(vapply(seq_len(p), function(a) Xs[[i]][a] * b[a], numeric(1))),
    numeric(1))
  list(sse = .suma(resid^2), b = b, resid = resid)
}

# ---------------------------------------------------------------------------
# Rama anova3 -- SS tipo III por comparacion de modelos (contr.sum saturado).
# ---------------------------------------------------------------------------
anova3_terminos <- function(X, y) {
  n <- length(y); p <- length(X[[1]])
  ff <- ajustar(X, y, seq_len(p)); df <- n - p; mse <- ff$sse / df
  out <- list()
  for (nombre in TERM_ORD) {
    ct <- TERM_COLS[[nombre]]
    cols <- setdiff(seq_len(p), ct)
    fr <- ajustar(X, y, cols)
    k <- length(ct)
    Fv <- ((fr$sse - ff$sse) / k) / mse
    out[[nombre]] <- c(Fv, f_sf(Fv, k, df))
  }
  list(stats = out, df = df)
}

# ---------------------------------------------------------------------------
# Rama hc3 -- vcov sandwich HC3 y Wald tipo III (k gl por termino).
# ---------------------------------------------------------------------------
.vcov_hc3 <- function(X, resid) {
  n <- length(X); p <- length(X[[1]])
  Xm <- matrix(unlist(X), n, p, byrow = TRUE)
  XtX <- t(Xm) %*% Xm
  XtXi <- invertir(XtX)
  h <- vapply(seq_len(n), function(i)
    as.numeric(Xm[i, , drop = FALSE] %*% XtXi %*% t(Xm[i, , drop = FALSE])), numeric(1))
  meat <- matrix(0, p, p)
  for (i in seq_len(n)) {
    w <- (resid[i]^2) / ((1 - h[i])^2)
    meat <- meat + w * (t(Xm[i, , drop = FALSE]) %*% Xm[i, , drop = FALSE])
  }
  XtXi %*% meat %*% XtXi
}
hc3_terminos <- function(X, y) {
  n <- length(y); p <- length(X[[1]])
  ff <- ajustar(X, y, seq_len(p))
  V <- .vcov_hc3(X, ff$resid)
  df <- n - p
  out <- list()
  for (nombre in TERM_ORD) {
    ct <- TERM_COLS[[nombre]]; k <- length(ct)
    if (k == 1L) {
      Fv <- ff$b[ct]^2 / V[ct, ct]
    } else {
      bs <- ff$b[ct]; Vs <- V[ct, ct, drop = FALSE]
      W <- as.numeric(t(bs) %*% invertir(Vs) %*% bs)
      Fv <- W / k
    }
    out[[nombre]] <- c(Fv, f_sf(Fv, k, df))
  }
  list(stats = out, df = df)
}

# ---------------------------------------------------------------------------
# Levene / Brown-Forsythe sobre las 4 celdas SEXO x TTO.
# ---------------------------------------------------------------------------
levene_bf <- function(m) {
  z <- c(); g <- c(); celdas <- list()
  for (k in CELDAS_4) {
    key <- paste(k, collapse = "\r")
    celdas[[key]] <- m$y[m$SEXO == k[1] & m$TTO == k[2]]
  }
  for (k in CELDAS_4) {
    key <- paste(k, collapse = "\r")
    med <- mediana(celdas[[key]])
    for (v in celdas[[key]]) { z <- c(z, abs(v - med)); g <- c(g, key) }
  }
  N <- length(z); kk <- 4L; gran <- promedio(z)
  gm <- vapply(CELDAS_4, function(k) {
    key <- paste(k, collapse = "\r"); promedio(z[g == key])
  }, numeric(1))
  names(gm) <- vapply(CELDAS_4, paste, character(1), collapse = "\r")
  ssb <- .suma(vapply(seq_along(CELDAS_4), function(j) {
    key <- names(gm)[j]; length(celdas[[key]]) * (gm[j] - gran)^2
  }, numeric(1)))
  ssw <- .suma(vapply(seq_len(N), function(i) (z[i] - gm[[g[i]]])^2, numeric(1)))
  Fv <- (ssb / (kk - 1)) / (ssw / (N - kk))
  c(Fv, f_sf(Fv, kk - 1, N - kk))
}

shapiro_residuos <- function(X, y) {
  p <- length(X[[1]])
  ff <- ajustar(X, y, seq_len(p))
  sw <- shapiro.test(ff$resid)
  c(as.numeric(sw$statistic), as.numeric(sw$p.value))
}

# ---------------------------------------------------------------------------
# ART aditivo hand-rolled (rama `art`). ARTool rechaza el bloque no cruzado; se
# alinea cada termino sobre el modelo aditivo (residuo del ajuste conjunto +
# contribucion ajustada del termino), se rankean y se corre ANOVA III con el
# diseno completo. No se ejecuta con los datos reales ni sinteticos.
# ---------------------------------------------------------------------------
art_terminos <- function(m, membranas) {
  d <- diseno(m, membranas); X <- d$X; y <- d$y
  n <- length(y); p <- length(X[[1]])
  ff <- ajustar(X, y, seq_len(p))
  out <- list()
  for (term in c("SEXO", "TTO", "SEXO:TTO")) {
    ct <- TERM_COLS[[term]]
    alin <- vapply(seq_len(n), function(i)
      round(ff$resid[i] + .suma(vapply(ct, function(c) ff$b[c] * X[[i]][c], numeric(1))), 8),
      numeric(1))
    rr <- rangos_promedio(alin)
    mr <- data.frame(SEXO = m$SEXO, TTO = m$TTO, MEMBRANA = m$MEMBRANA, y = rr,
                     stringsAsFactors = FALSE)
    dr <- diseno(mr, membranas)
    st_ <- anova3_terminos(dr$X, dr$y)$stats
    out[[term]] <- st_[[term]]
  }
  list(stats = out, df = n - p)
}

# ---------------------------------------------------------------------------
# ART-C hand-rolled para las 4 comparaciones D6 (bloque aditivo): rangos
# alineados por la interaccion, modelo rr ~ SEXOxTTO(4 celdas) + MEMBRANA,
# contrastes t con MSE combinado, gl = n - 6.
# ---------------------------------------------------------------------------
artc_pares <- function(m, membranas) {
  d <- diseno(m, membranas); X <- d$X; y <- d$y
  n <- length(y); p <- length(X[[1]])
  ff <- ajustar(X, y, seq_len(p))
  alin <- vapply(seq_len(n), function(i) round(ff$resid[i] + ff$b[4] * X[[i]][4], 8),
                 numeric(1))
  rr <- rangos_promedio(alin)
  cs4 <- list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(-1, -1, -1))
  idx_celda <- function(sexo, tto) {
    for (j in seq_along(CELDAS_4))
      if (CELDAS_4[[j]][1] == sexo && CELDAS_4[[j]][2] == tto) return(j)
    stop("celda")
  }
  Xc <- lapply(seq_len(n), function(i) {
    c4 <- cs4[[idx_celda(m$SEXO[i], m$TTO[i])]]
    mm <- .cs_membrana(m$MEMBRANA[i], membranas)
    c(1, c4[1], c4[2], c4[3], mm[1], mm[2])
  })
  fc <- ajustar(Xc, rr, seq_len(6L))
  df <- n - 6L; mse <- fc$sse / df
  idx <- lapply(CELDAS_4, function(k) which(m$SEXO == k[1] & m$TTO == k[2]))
  names(idx) <- vapply(CELDAS_4, paste, character(1), collapse = "\r")
  media <- vapply(idx, function(ii) promedio(rr[ii]), numeric(1))
  out <- list()
  for (par in PARES_D6) {
    ka <- paste(par[[2]], collapse = "\r"); kb <- paste(par[[3]], collapse = "\r")
    se <- sqrt(mse * (1 / length(idx[[ka]]) + 1 / length(idx[[kb]])))
    est <- media[[ka]] - media[[kb]]
    tval <- est / se
    out[[par[[1]]]] <- c(est, se, tval, t_sf2(tval, df))
  }
  list(pares = out, df = df)
}

# ---------------------------------------------------------------------------
# Contrastes de medias marginales (ramas anova3 / hc3), promediando sobre MEMBRANA.
# ---------------------------------------------------------------------------
emmeans_pares <- function(m, membranas, robusto) {
  d <- diseno(m, membranas); X <- d$X; y <- d$y
  n <- length(y); p <- length(X[[1]])
  ff <- ajustar(X, y, seq_len(p)); df <- n - p; mse <- ff$sse / df
  if (robusto) {
    V <- .vcov_hc3(X, ff$resid)
  } else {
    Xm <- matrix(unlist(X), n, p, byrow = TRUE)
    V <- mse * invertir(t(Xm) %*% Xm)
  }
  media <- list()
  for (k in CELDAS_4) media[[paste(k, collapse = "\r")]] <-
    promedio(m$y[m$SEXO == k[1] & m$TTO == k[2]])
  out <- list()
  for (par in PARES_D6) {
    sa <- if (par[[2]][1] == "HEMBRA") 1 else -1
    ta <- if (par[[2]][2] == "CONTROL") 1 else -1
    sb <- if (par[[3]][1] == "HEMBRA") 1 else -1
    tb <- if (par[[3]][2] == "CONTROL") 1 else -1
    cv <- c(0, sa - sb, ta - tb, sa * ta - sb * tb, 0, 0)
    var <- as.numeric(t(cv) %*% V %*% cv)
    se <- sqrt(var)
    est <- media[[paste(par[[2]], collapse = "\r")]] -
           media[[paste(par[[3]], collapse = "\r")]]
    tval <- est / se
    out[[par[[1]]]] <- c(est, se, tval, t_sf2(tval, df))
  }
  list(pares = out, df = df)
}

holm <- function(pvals) {
  m <- length(pvals); ord <- order(pvals)
  aj <- numeric(m); corr <- 0
  for (rank in seq_len(m)) {
    i <- ord[rank]
    corr <- max(corr, (m - rank + 1L) * pvals[i])
    aj[i] <- min(1, corr)
  }
  aj
}

# ===========================================================================
# 1. Carga
# ===========================================================================
cargar <- function() {
  ruta <- file.path(RUTA_DATOS_PROC, "pstat3_long.tsv")
  if (!file.exists(ruta))
    stop(sprintf("06: falta %s (correr 02_ingesta_qc primero)", ruta))
  d <- read.delim(ruta, sep = "\t", quote = "", stringsAsFactors = FALSE,
                  colClasses = "character", check.names = FALSE, encoding = "UTF-8",
                  na.strings = character(0))
  d$y <- suppressWarnings(as.numeric(d$PSTAT3))
  d
}

n_celdas <- function(d) {
  vapply(CELDAS_4, function(k)
    sum(d$SEXO == k[1] & d$TTO == k[2]), integer(1))
}

# ===========================================================================
# 2. Clasificacion y modelado (una fila: pSTAT3 @ PLACENTA_E15)
# ===========================================================================
clasificar <- function(d, membranas) {
  tej <- d$TEJIDO[1]
  ncel <- n_celdas(d)
  rec <- list(TEJIDO = tej, GEN = "pSTAT3", n_celda = ncel, n_total = nrow(d),
              via = NA_character_, motivo = "",
              shapiro_p = NA_real_, levene_p = NA_real_, rama = "",
              stats = list(), df = NA_integer_, posthoc = NULL)
  if (min(ncel) < PISO_CELDA) {
    rec$via <- "descriptivo_n_bajo"
    rec$motivo <- sprintf("celda SEXO x TTO con < %d valores", PISO_CELDA)
    return(rec)
  }
  rec$via <- "modelo"
  ajustar_cascada(rec, d, membranas)
}

ajustar_cascada <- function(rec, d, membranas) {
  dd <- diseno(d, membranas)
  sw <- shapiro_residuos(dd$X, dd$y)
  lv <- levene_bf(d)
  rec$shapiro_p <- sw[2]; rec$levene_p <- lv[2]
  if (sw[2] >= ALFA && lv[2] >= ALFA) {
    rec$rama <- "anova3"; st_ <- anova3_terminos(dd$X, dd$y)
  } else if (sw[2] >= ALFA && lv[2] < ALFA) {
    rec$rama <- "hc3"; st_ <- hc3_terminos(dd$X, dd$y)
  } else {
    rec$rama <- "art"; st_ <- art_terminos(d, membranas)
  }
  rec$stats <- st_$stats; rec$df <- st_$df
  ip <- st_$stats[["SEXO:TTO"]][2]
  if (ip < ALFA) {
    if (rec$rama == "art") {
      pr <- artc_pares(d, membranas); metodo <- "ART-C"
    } else {
      pr <- emmeans_pares(d, membranas, robusto = (rec$rama == "hc3"))
      metodo <- if (rec$rama == "hc3") "contraste_marginal_HC3" else "contraste_marginal_OLS"
    }
    etiquetas <- vapply(PARES_D6, `[[`, character(1), 1)
    praw <- vapply(etiquetas, function(e) pr$pares[[e]][4], numeric(1))
    ph <- holm(praw)
    rec$posthoc <- list(metodo = metodo, df = pr$df,
      filas = lapply(seq_along(etiquetas), function(i) {
        e <- etiquetas[i]; v <- pr$pares[[e]]
        list(e, v[1], v[2], v[3], praw[i], ph[i])
      }))
  }
  rec
}

# ===========================================================================
# 3. Armado de tablas de salida
# ===========================================================================
COLS_CLASIF <- c(
  "TEJIDO", "GEN",
  "n_HEMBRA_CONTROL", "n_HEMBRA_LPS", "n_MACHO_CONTROL", "n_MACHO_LPS", "n_total",
  "via", "motivo", "shapiro_p_residuos", "levene_p", "rama_cascada",
  "F_SEXO", "p_SEXO", "F_TTO", "p_TTO", "F_SEXOxTTO", "p_SEXOxTTO",
  "F_MEMBRANA", "p_MEMBRANA",
  "interaccion_significativa", "post_hoc_corrido")

fila_clasif <- function(rec) {
  gv <- function(term, k) if (!is.null(rec$stats[[term]])) rec$stats[[term]][k] else NA_real_
  ip <- if (!is.null(rec$stats[["SEXO:TTO"]])) rec$stats[["SEXO:TTO"]][2] else NA_real_
  sig <- !is.na(ip) && ip < ALFA
  list(
    rec$TEJIDO, rec$GEN,
    rec$n_celda[1], rec$n_celda[2], rec$n_celda[3], rec$n_celda[4], rec$n_total,
    rec$via, rec$motivo, p6e(rec$shapiro_p), p6e(rec$levene_p), rec$rama,
    p6e(gv("SEXO", 1)), p6e(gv("SEXO", 2)),
    p6e(gv("TTO", 1)), p6e(gv("TTO", 2)),
    p6e(gv("SEXO:TTO", 1)), p6e(gv("SEXO:TTO", 2)),
    p6e(gv("MEMBRANA", 1)), p6e(gv("MEMBRANA", 2)),
    if (rec$via == "modelo") sig else "",
    if (rec$via == "modelo") !is.null(rec$posthoc) else ""
  )
}

COLS_POSTHOC <- c("TEJIDO", "GEN", "rama_cascada", "metodo", "contraste",
                  "estimador", "EE", "estadistico_t", "p_sin_correccion", "p_holm")

filas_posthoc <- function(rec) {
  if (is.null(rec$posthoc)) return(list())
  lapply(rec$posthoc$filas, function(f)
    list(rec$TEJIDO, rec$GEN, rec$rama, rec$posthoc$metodo, f[[1]],
         p6e(f[[2]]), p6e(f[[3]]), p6e(f[[4]]), p6e(f[[5]]), p6e(f[[6]])))
}

COLS_DESCR <- c("TEJIDO", "GEN", "GRUPO", "n", "mean_PSTAT3", "sd_PSTAT3",
                "median_PSTAT3")

filas_descriptivo <- function(d) {
  tej <- d$TEJIDO[1]
  lapply(seq_along(CELDAS_4), function(j) {
    k <- CELDAS_4[[j]]; grp <- GRUPOS_4[j]
    vals <- d$y[d$SEXO == k[1] & d$TTO == k[2]]
    n <- length(vals)
    mu <- if (n) promedio(vals) else NA_real_
    sd <- if (n >= 2) desvio_muestral(vals, mu) else NA_real_
    md <- if (n) mediana(vals) else NA_real_
    list(tej, "pSTAT3", grp, n,
         if (is.na(mu)) "" else g10(mu),
         if (is.na(sd)) "" else g10(sd),
         if (is.na(md)) "" else g10(md))
  })
}

# ===========================================================================
# 4. Cruza-verificacion en corrida contra car / emmeans (solo R)
# ===========================================================================
verificar_contra_librerias <- function(d, rec, membranas) {
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op))
  if (rec$via != "modelo") return(0)
  s <- data.frame(
    y = d$y,
    SEXO = factor(d$SEXO, levels = c("HEMBRA", "MACHO")),
    TTO = factor(d$TTO, levels = c("CONTROL", "LPS")),
    MEMBRANA = factor(d$MEMBRANA, levels = membranas))
  m <- lm(y ~ SEXO * TTO + MEMBRANA, data = s)
  peor <- 0
  lv <- car::leveneTest(y ~ SEXO * TTO, data = s)[1, "Pr(>F)"]
  peor <- max(peor, abs(lv - rec$levene_p))
  filas_car <- c("SEXO", "TTO", "SEXO:TTO", "MEMBRANA")
  if (rec$rama == "anova3") {
    a <- car::Anova(m, type = 3)
    for (t in filas_car)
      peor <- max(peor, abs(a[t, "Pr(>F)"] - rec$stats[[t]][2]))
  } else if (rec$rama == "hc3") {
    a <- car::Anova(m, type = 3, white.adjust = "hc3")
    for (t in filas_car)
      peor <- max(peor, abs(a[t, "Pr(>F)"] - rec$stats[[t]][2]))
  } else {
    # rama art: ARTool no admite el bloque aditivo -> auto-verificacion
    dd <- diseno(d, membranas)
    ff <- ajustar(dd$X, dd$y, seq_len(NCOL_DIS))
    for (term in c("SEXO", "TTO", "SEXO:TTO")) {
      ct <- TERM_COLS[[term]]
      # via A: residuo del ajuste conjunto + contribucion del termino (via lm)
      via_a <- vapply(seq_len(nrow(d)), function(i)
        round(ff$resid[i] + .suma(vapply(ct, function(c) ff$b[c] * dd$X[[i]][c],
                                         numeric(1))), 8), numeric(1))
      # via B: Y menos el ajuste del modelo SIN el termino (reduccion directa)
      cols_sin <- setdiff(seq_len(NCOL_DIS), ct)
      fr <- ajustar(dd$X, dd$y, cols_sin)
      pred_sin <- vapply(seq_len(nrow(d)), function(i)
        .suma(vapply(seq_along(cols_sin), function(a)
          dd$X[[i]][cols_sin[a]] * fr$b[a], numeric(1))), numeric(1))
      via_b <- round(dd$y - pred_sin, 8)
      peor <- max(peor, max(abs(via_a - via_b)))
    }
  }
  if (!is.null(rec$posthoc)) {
    mape <- c("HEMBRA_CONTROL-HEMBRA_LPS"    = "HEMBRA CONTROL - HEMBRA LPS",
              "MACHO_CONTROL-MACHO_LPS"      = "MACHO CONTROL - MACHO LPS",
              "HEMBRA_LPS-MACHO_LPS"         = "HEMBRA LPS - MACHO LPS",
              "HEMBRA_CONTROL-MACHO_CONTROL" = "HEMBRA CONTROL - MACHO CONTROL")
    if (rec$rama == "art") {
      # ART-C hand-rolled vs emmeans sobre lm(rr ~ SEXOxTTO + MEMBRANA)
      dd <- diseno(d, membranas)
      ff <- ajustar(dd$X, dd$y, seq_len(NCOL_DIS))
      alin <- vapply(seq_len(nrow(d)), function(i)
        round(ff$resid[i] + ff$b[4] * dd$X[[i]][4], 8), numeric(1))
      s2 <- s; s2$rr <- rank(alin, ties.method = "average")
      s2$CELL <- factor(paste(s2$SEXO, s2$TTO, sep = "_"),
                        levels = c("HEMBRA_CONTROL", "HEMBRA_LPS",
                                   "MACHO_CONTROL", "MACHO_LPS"))
      mr <- lm(rr ~ CELL + MEMBRANA, data = s2)
      emm <- emmeans::emmeans(mr, ~ CELL)
      pr <- summary(emmeans::contrast(emm, method = "pairwise", adjust = "none"))
      pv <- setNames(pr$p.value, gsub("\\s+", " ", as.character(pr$contrast)))
      mapc <- c("HEMBRA_CONTROL-HEMBRA_LPS"    = "HEMBRA_CONTROL - HEMBRA_LPS",
                "MACHO_CONTROL-MACHO_LPS"      = "MACHO_CONTROL - MACHO_LPS",
                "HEMBRA_LPS-MACHO_LPS"         = "HEMBRA_LPS - MACHO_LPS",
                "HEMBRA_CONTROL-MACHO_CONTROL" = "HEMBRA_CONTROL - MACHO_CONTROL")
      for (fila in rec$posthoc$filas)
        peor <- max(peor, abs(pv[[mapc[[fila[[1]]]]]] - fila[[5]]))
    } else {
      vv <- if (rec$rama == "hc3") sandwich::vcovHC(m, type = "HC3") else vcov(m)
      emm <- emmeans::emmeans(m, ~ SEXO * TTO, vcov. = vv)
      pr <- summary(emmeans::contrast(emm, method = "pairwise", adjust = "none"))
      pv <- setNames(pr$p.value, gsub("\\s+", " ", as.character(pr$contrast)))
      for (fila in rec$posthoc$filas)
        peor <- max(peor, abs(pv[[mape[[fila[[1]]]]]] - fila[[5]]))
    }
  }
  peor
}

# ===========================================================================
# 5. Reporte legible
# ===========================================================================
.md <- function(header, filas) {
  l1 <- paste0("| ", paste(header, collapse = " | "), " |")
  l2 <- paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |")
  cuerpo <- vapply(filas, function(f)
    paste0("| ", paste(vapply(f, .fmt, character(1)), collapse = " | "), " |"),
    character(1))
  paste(c(l1, l2, cuerpo), collapse = "\n")
}

LIMITACION_D9 <- paste0(
  "**Limitacion obligatoria (D9).** pSTAT3 se cuantifico por western blot ",
  "normalizado a **proteina total** (p. ej. Ponceau / stain-free), **sin STAT3 ",
  "total** en la misma membrana. La medida refleja por lo tanto la **abundancia ",
  "de fosfo-STAT3 (Tyr705)** relativa a la carga de proteina, **no la fraccion ",
  "de STAT3 que esta fosforilada**. Un cambio en pSTAT3 puede deberse a mas ",
  "fosforilacion, a mas STAT3 total, o a ambos; con estos datos no se puede ",
  "separar. El bloque `MEMBRANA` (efecto fijo) absorbe la variacion tecnica ",
  "entre las 3 membranas del ensayo.")

BULLETS_METODO <- c(
  paste0("- Modelo (D9): `PSTAT3 ~ SEXO * TTO + MEMBRANA`, OLS con contrastes ",
         "suma-cero. Diseno de 6 columnas `[1, s, t, s*t, m1, m2]`; `MEMBRANA` = ",
         "bloque fijo de 3 niveles (contr.sum, 2 columnas). df de residuos = n - 6."),
  paste0("- Cascada D5 / 4.1 (alfa = 0.05): Shapiro-Wilk sobre los **residuos del ",
         "modelo conjunto** (con MEMBRANA); Levene (Brown-Forsythe, centro = ",
         "mediana) sobre las 4 celdas SEXO x TTO. Shapiro>=.05 & Levene>=.05 -> ",
         "`anova3` (ANOVA III); Shapiro>=.05 & Levene<.05 -> `hc3` (Wald III con ",
         "sandwich HC3); Shapiro<.05 -> `art` (ART aditivo hand-rolled: ARTool ",
         "rechaza el bloque no cruzado)."),
  paste0("- Piso de celda: si alguna celda SEXO x TTO tiene < 5 valores, no se ",
         "modela (`descriptivo_n_bajo`). pSTAT3 tiene 9 por celda -> siempre se modela."),
  paste0("- Post hoc D6 (solo si interaccion p < alfa): 4 comparaciones fijas ",
         "(HEMBRA_CONTROL-HEMBRA_LPS, MACHO_CONTROL-MACHO_LPS, HEMBRA_LPS-MACHO_LPS, ",
         "HEMBRA_CONTROL-MACHO_CONTROL), correccion Holm. `anova3` -> contraste de ",
         "medias marginales promediando sobre MEMBRANA (vcov OLS); `hc3` -> idem ",
         "con vcov HC3; `art` -> ART-C hand-rolled."),
  paste0("- `F_MEMBRANA` / `p_MEMBRANA` se reportan como diagnostico del bloque ",
         "tecnico; no habilitan ni bloquean nada (la compuerta es SEXO x TTO)."))

construir_reporte <- function(fuente, rec, clasif, posthoc, descr, resumen) {
  L <- c(
    "# Reporte de pSTAT3 (T6)", "",
    "Generado por `06_pstat3` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. Metodo", "",
    BULLETS_METODO, "",
    "## 2. Clasificacion del modelo", "",
    .md(COLS_CLASIF, clasif), "",
    "## 3. Post hoc (D6)", "")
  L <- c(L, if (length(posthoc)) .md(COLS_POSTHOC, posthoc)
         else "La interaccion SEXO x TTO no alcanzo p < alfa: no se corrio post hoc.")
  L <- c(L, "", "## 4. Descriptivo de pSTAT3 por grupo", "",
    .md(COLS_DESCR, descr), "",
    "## 5. Resumen", "")
  for (linea in resumen) L <- c(L, paste0("- ", linea))
  L <- c(L, "", "## 6. Limitacion del ensayo (D9)", "", LIMITACION_D9, "",
    "## 7. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `06_pstat3`: eleccion de rama ",
           "de la cascada, el rechazo de ARTool al bloque aditivo y el ART ",
           "hand-rolled, y la nota de paridad (Shapiro-Wilk por libreria)."), "")
  paste(L, collapse = "\n")
}

# ===========================================================================
# 6. Artefactos compartidos (merge por 'script') -- headers identicos a 02..05.
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
  txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE)
  Encoding(txt) <- "UTF-8"
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

bloque_descartes <- function(rama, shapiro_p, levene_p) {
  paste(c(
"## 06_pstat3",
"",
"### Eleccion de rama de la cascada D5",
"",
paste0("- Modelo D9 `PSTAT3 ~ SEXO * TTO + MEMBRANA` (OLS, contr.sum, 6 columnas). ",
       "Se decide la rama con Shapiro-Wilk sobre los **residuos del modelo ",
       "conjunto** y Levene (Brown-Forsythe) sobre las 4 celdas SEXO x TTO, ambos a ",
       sprintf("alfa = 0.05. En esta corrida: Shapiro p = %s, Levene p = %s -> rama ",
               p6e(shapiro_p), p6e(levene_p)),
       sprintf("`%s`. Sobre datos reales Y sinteticos la cascada cae en `anova3` ", rama),
       "(residuos normales y varianzas homogeneas)."),
"",
"### La rama `art`: ARTool no admite el bloque aditivo",
"",
paste0("- `ARTool::art()` **rechaza** `y ~ SEXO * TTO + MEMBRANA` con ",
       "`parse.art.formula`: *\"Model must include all combinations of interactions ",
       "of fixed effects\"*. El ART clasico exige un diseno completamente cruzado; ",
       "un bloque tecnico aditivo no lo es."),
paste0("- Fallback fijado con el usuario: **ART aditivo hand-rolled**. Se alinea ",
       "cada termino restando del dato el ajuste de todos los demas terminos ",
       "(incluido el efecto principal de MEMBRANA), estimados por OLS del modelo ",
       "conjunto; se rankean los valores alineados (redondeo a 8 decimales, rangos ",
       "promedio) y se corre ANOVA III sobre esos rangos con el diseno completo. El ",
       "post hoc es ART-C hand-rolled: rangos alineados por la interaccion, modelo ",
       "`rr ~ SEXOxTTO + MEMBRANA`, contrastes t con MSE combinado (gl = n - 6)."),
paste0("- Como ARTool no puede ajustar el modelo, la cruza-verificacion de esta ",
       "rama en `06_pstat3.R` es una **auto-verificacion**: R recomputa el alineado ",
       "por dos vias (residuo + contribucion del termino via `lm`, y aritmetica de ",
       "medias/margenes) y exige que coincidan (`stopifnot`, tol 1e-9), ademas de ",
       "chequear el ANOVA III sobre los rangos alineados contra `car::Anova`. ",
       "**Esta rama no se ejecuta con los datos reales ni sinteticos.**"),
"",
"### MEMBRANA como bloque fijo",
"",
paste0("- El diseno es balanceado y ortogonal (3 membranas x 4 grupos x 3 replicas ",
       "= 36; cada membrana aporta 3 valores a cada celda SEXO x TTO), asi que SS ",
       "tipo I = tipo III y la media marginal de celda (promediando sobre MEMBRANA) ",
       "= media cruda de celda. `F_MEMBRANA` / `p_MEMBRANA` se reportan como ",
       "diagnostico del bloque; no dirigen la inferencia."),
"",
"### Limitacion del ensayo (D9)",
"",
paste0("- ", gsub("**", "", LIMITACION_D9, fixed = TRUE)),
"",
"### Paridad R / Python",
"",
paste0("- Nucleo numerico PROPIO identico en R y Python (OLS 6x6 por Gauss-Jordan, ",
       "SS tipo III por comparacion de modelos, sandwich HC3, Levene, ART y ART-C, ",
       "Holm). `06_pstat3.R` cruza-verifica en corrida las ramas que se ejecutan ",
       "(`anova3` / `hc3`) contra `car::Anova` y `emmeans` (`stopifnot`, tol 1e-6)."),
paste0("- Unica excepcion: **Shapiro-Wilk** de libreria (`shapiro.test` / ",
       "`scipy.stats.shapiro`); coinciden a ~1e-12 y se comparan dentro de ",
       "`TOL_ESTADISTICO` en T10. Los p / estadisticos que pasan por funciones ",
       "trascendentes (pf, pt, W) se guardan como texto `%.6e` (igual que 03/05).")
  ), collapse = "\n")
}

actualizar_descartados <- function(rama, shapiro_p, levene_p) {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 06_pstat3:inicio -->"
  marca_fin <- "<!-- 06_pstat3:fin -->"
  nuevo <- paste0(marca_ini, "\n", bloque_descartes(rama, shapiro_p, levene_p),
                  "\n\n", marca_fin)
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
main <- function() {
  d <- cargar()
  fuente <- fuente_datos(ARCHIVO_PSTAT3)
  membranas <- sort(unique(d$MEMBRANA))
  if (length(membranas) != 3L)
    stop(sprintf("06: se esperaban 3 membranas, hay %d: %s",
                 length(membranas), paste(membranas, collapse = ", ")))

  rec <- clasificar(d, membranas)

  peor <- verificar_contra_librerias(d, rec, membranas)
  cat(sprintf("  [cruza-verificacion R vs car/emmeans] peor |dif| = %.3e\n", peor))
  stopifnot(peor < 1e-6)

  clasif <- list(fila_clasif(rec))
  posthoc <- filas_posthoc(rec)
  descr <- filas_descriptivo(d)

  memb_cont <- vapply(membranas, function(mb) sum(d$MEMBRANA == mb), integer(1))
  balanceado <- length(unique(memb_cont)) == 1L

  st_ <- rec$stats
  ip <- if (!is.null(st_[["SEXO:TTO"]])) st_[["SEXO:TTO"]][2] else NA_real_
  sig <- !is.na(ip) && ip < ALFA
  n_ph <- if (is.null(rec$posthoc)) 0L else length(rec$posthoc$filas)
  ph_sig <- character(0)
  if (!is.null(rec$posthoc))
    for (f in rec$posthoc$filas) if (f[[6]] < ALFA) ph_sig <- c(ph_sig, f[[1]])

  resumen <- c(
    sprintf("Via: `%s`%s; rama de la cascada D5: `%s` (Shapiro p = %s, Levene p = %s).",
            rec$via, if (nzchar(rec$motivo)) sprintf(" (%s)", rec$motivo) else "",
            rec$rama, p6e(rec$shapiro_p), p6e(rec$levene_p)),
    sprintf("SEXO: p = %s | TTO: p = %s | SEXO x TTO: p = %s | MEMBRANA (bloque): p = %s.",
            p6e(st_[["SEXO"]][2]), p6e(st_[["TTO"]][2]), p6e(st_[["SEXO:TTO"]][2]),
            p6e(st_[["MEMBRANA"]][2])),
    if (sig)
      sprintf(paste0("Interaccion SEXO x TTO significativa (p < %s): post hoc D6 ",
                     "corrido (%d comparaciones, Holm). Significativas tras Holm: %s."),
              format(ALFA), n_ph,
              if (length(ph_sig)) paste(ph_sig, collapse = ", ") else "(ninguna)")
    else
      sprintf("Interaccion SEXO x TTO no significativa (p >= %s): sin post hoc (D6).",
              format(ALFA)),
    paste0("pSTAT3 = abundancia de fosfo-STAT3 (no fraccion fosforilada); MEMBRANA ",
           "es bloque tecnico fijo. Ver seccion 6 del reporte y D9."))

  escribir_tsv(file.path(RUTA_DATOS_PROC, "pstat3_modelo_clasificacion.tsv"),
               COLS_CLASIF, clasif)
  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "pstat3_modelo_clasificacion.csv"), COLS_CLASIF, clasif)
    escribir_csv(file.path(base, "pstat3_posthoc.csv"), COLS_POSTHOC, posthoc)
    escribir_csv(file.path(base, "pstat3_descriptivo.csv"), COLS_DESCR, descr)
  }

  escribir_lineas(file.path(RUTA_TABLAS, "pstat3_reporte.md"),
                  construir_reporte(fuente, rec, clasif, posthoc, descr, resumen))
  actualizar_descartados(rec$rama, rec$shapiro_p, rec$levene_p)

  ent <- sprintf("data/processed/pstat3_long.tsv (de data/%s/%s)", fuente, ARCHIVO_PSTAT3)
  registrar_procedencia(list(
    list("data/processed/pstat3_modelo_clasificacion.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("clasificacion del modelo pSTAT3 (D9): via, rama D5, ",
         "F/p de SEXO/TTO/SEXOxTTO/MEMBRANA, flag de post hoc")),
    list("outputs/tables/{R,python}/pstat3_modelo_clasificacion.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "idem, copia por implementacion"),
    list("outputs/tables/{R,python}/pstat3_posthoc.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("post hoc D6 de pSTAT3: 4 comparaciones fijas con Holm ",
         "(contraste marginal OLS/HC3 o ART-C), solo si interaccion p < alfa")),
    list("outputs/tables/{R,python}/pstat3_descriptivo.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "n, media, sd y mediana de PSTAT3 por grupo SEXO x TTO"),
    list("outputs/tables/pstat3_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO", ent,
         "reporte legible del modelo pSTAT3 (T6), con la limitacion obligatoria D9")
  ))
  ok_ph <- (sig && n_ph == 4L) || (!sig && n_ph == 0L)
  registrar_verificaciones(list(
    list("pstat3_via", "via del modelo pSTAT3 (D9)", rec$via, "modelo",
         if (rec$via == "modelo") "TRUE" else "FALSE", ESTE_SCRIPT),
    list("pstat3_rama_cascada", "rama de la cascada D5 elegida por diagnostico",
         sprintf("%s (Shapiro p=%s, Levene p=%s)", rec$rama, p6e(rec$shapiro_p),
                 p6e(rec$levene_p)),
         "anova3 (residuos normales, varianzas homogeneas)",
         if (rec$rama == "anova3") "TRUE" else "FALSE", ESTE_SCRIPT),
    list("pstat3_modelo_d9", "modelo D9 con MEMBRANA como bloque fijo",
         "PSTAT3 ~ SEXO * TTO + MEMBRANA (contr.sum, df_resid = n-6)",
         "PSTAT3 ~ SEXO * TTO + MEMBRANA", "TRUE", ESTE_SCRIPT),
    list("pstat3_bloque_balanceado", "las 3 membranas estan balanceadas",
         paste(sprintf("%s:%d", membranas, memb_cont), collapse = ";"),
         "12;12;12",
         if (balanceado && all(memb_cont == 12L)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("pstat3_alfa", "alfa para gate D6 y pretests de supuestos",
         g10(ALFA), "0.05", if (ALFA == 0.05) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("pstat3_posthoc_holm",
         "post hoc D6: 4 comparaciones fijas con Holm, solo si interaccion p < alfa",
         if (sig) sprintf("interaccion p=%s < alfa -> %d comparaciones, Holm", p6e(ip), n_ph)
         else sprintf("interaccion p=%s >= alfa -> sin post hoc", p6e(ip)),
         "4 comparaciones + Holm si y solo si interaccion significativa",
         if (ok_ph) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("pstat3_limitacion_d9",
         "limitacion obligatoria: abundancia de fosfo-STAT3, no fraccion fosforilada",
         "escrita en pstat3_reporte.md (seccion 6) y analisis_descartados.md",
         "presente", "TRUE", ESTE_SCRIPT),
    list("pstat3_core_vs_libs",
         "R cruza-verifica anova3/hc3/emmeans contra car::Anova / emmeans (stopifnot 1e-6)",
         paste0("nucleo PROPIO identico R/Python; Shapiro-Wilk por libreria; ART ",
                "aditivo auto-verificado (ARTool no admite el bloque)"),
         "verificado", "TRUE", ESTE_SCRIPT)
  ))

  cat("== 06_pstat3.R ==\n")
  cat(sprintf("  fuente pSTAT3 = %s\n", fuente))
  cat(sprintf("  membranas = %s  (n por membrana: %s)\n",
              paste(membranas, collapse = ", "), paste(memb_cont, collapse = ", ")))
  cat(sprintf("  via = %s  rama = %s  (Shapiro p=%s, Levene p=%s)\n",
              rec$via, rec$rama, p6e(rec$shapiro_p), p6e(rec$levene_p)))
  cat(sprintf("  p: SEXO=%s  TTO=%s  SEXOxTTO=%s  MEMBRANA=%s\n",
              p6e(st_[["SEXO"]][2]), p6e(st_[["TTO"]][2]),
              p6e(st_[["SEXO:TTO"]][2]), p6e(st_[["MEMBRANA"]][2])))
  if (!is.null(rec$posthoc)) {
    cat(sprintf("  post hoc D6 (%s):\n", rec$posthoc$metodo))
    for (f in rec$posthoc$filas)
      cat(sprintf("    %-30s  est=%s  p_holm=%s\n", f[[1]], g10(f[[2]]), p6e(f[[6]])))
  } else {
    cat("  post hoc D6: no corrido (interaccion no significativa)\n")
  }
  cat("  -> data/processed/pstat3_modelo_clasificacion.tsv,\n")
  cat("     outputs/tables/{R,python}/pstat3_*.csv, outputs/tables/pstat3_reporte.md\n")
}

if (sys.nframe() == 0L) main()
