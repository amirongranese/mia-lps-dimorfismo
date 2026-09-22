# 05_qpcr_modelos.R -- Modelado de la qPCR: cascada D5, post hoc D6, Fisher D7,
#                      columna suplementaria BH D12.
#
# Por que existe este archivo: sobre `neg_ddCt` de 04_qpcr_cuantificacion, ajusta
# UN modelo por gen x tejido y clasifica cada uno segun el diagnostico de
# supuestos, sin elegir el metodo por conveniencia.
#
#   * D5 / seccion 4.1 -- cascada de supuestos (alfa = 0.05 en todo):
#       modelo base  neg_ddCt ~ SEXO * TTO   (OLS, contrastes suma-cero).
#       - Shapiro-Wilk sobre los RESIDUOS del modelo conjunto del gen x tejido
#         (no por celda) y Levene (Brown-Forsythe, centro = mediana) sobre las 4
#         celdas SEXO x TTO.
#       - Shapiro p >= .05  y  Levene p >= .05  -> ANOVA III  (rama "anova3")
#       - Shapiro p >= .05  y  Levene p <  .05  -> Wald III con HC3 (rama "hc3")
#       - Shapiro p <  .05  (con o sin Levene)  -> ART (rama "art")
#   * Piso de celda: si alguna de las 4 celdas SEXO x TTO tiene < 5 valores
#     DETECTADOS, no se ajusta el modelo -> via "descriptivo_n_bajo".
#   * D7 -- gen x tejido con calibrador HEMBRA_CONTROL 0/detectados
#     (`cuantificable == FALSE` en 04): fuera del modelo -> via "D7_deteccion"
#     (Fisher exacto 2x2 Control vs LPS dentro de cada sexo + tabla 2x4).
#   * D6 -- post hoc SOLO si la interaccion SEXO x TTO es significativa (p < .05).
#     4 comparaciones fijas + correccion Holm. Metodo por rama:
#       anova3 -> contraste de medias marginales (vcov OLS)
#       hc3    -> idem con vcov HC3
#       art    -> ART-C (Elkin et al. 2021)
#   * D12 -- columna suplementaria: p ajustado por Benjamini-Hochberg dentro de
#     cada tejido entre los genes MODELADOS, por termino. No dirige la inferencia.
#
# PARIDAD R/Python: el nucleo numerico (OLS 4x4 por eliminacion de Gauss, SS tipo
# III por comparacion de modelos, sandwich HC3, Levene, ART y ART-C, Holm, BH) es
# codigo PROPIO identico en ambos lenguajes. Aca (solo R) se cruza-verifica en
# corrida contra car::Anova / emmeans / ARTool con stopifnot (< 1e-6). Unico
# componente de libreria en el resultado: Shapiro-Wilk (shapiro.test / scipy).
# Los estadisticos y p dependientes de trascendentes se guardan como texto
# "%.6e" (p6e), como en 03_elisa.
#
# Este script NO grafica.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({
  library(car); library(emmeans); library(sandwich); library(ARTool)
})

ESTE_SCRIPT <- "05_qpcr_modelos"

ALFA        <- 0.05
PISO_CELDA  <- 5L

GRUPOS_4 <- c("HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")
CELDAS_4 <- list(c("HEMBRA", "CONTROL"), c("HEMBRA", "LPS"),
                 c("MACHO", "CONTROL"), c("MACHO", "LPS"))
SET_TRANSP <- GENES_TRANSPORTADORES

# D6: (etiqueta, celda_a, celda_b) -> estima a - b
PARES_D6 <- list(
  list("HEMBRA_CONTROL-HEMBRA_LPS",   c("HEMBRA", "CONTROL"), c("HEMBRA", "LPS")),
  list("MACHO_CONTROL-MACHO_LPS",     c("MACHO", "CONTROL"),  c("MACHO", "LPS")),
  list("HEMBRA_LPS-MACHO_LPS",        c("HEMBRA", "LPS"),     c("MACHO", "LPS")),
  list("HEMBRA_CONTROL-MACHO_CONTROL", c("HEMBRA", "CONTROL"), c("MACHO", "CONTROL"))
)

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02/03/04.
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

f_sf   <- function(x, d1, d2) pf(x, d1, d2, lower.tail = FALSE)
t_sf2  <- function(x, df) 2 * pt(-abs(x), df)

# ---------------------------------------------------------------------------
# Diseno suma-cero [1, s, t, s*t]; s=+1 HEMBRA/-1 MACHO, t=+1 CONTROL/-1 LPS.
# ---------------------------------------------------------------------------
fila_diseno <- function(sexo, tto) {
  s <- if (sexo == "HEMBRA") 1 else -1
  t <- if (tto == "CONTROL") 1 else -1
  c(1, s, t, s * t)
}
diseno <- function(m) {
  X <- lapply(seq_len(nrow(m)), function(i) fila_diseno(m$SEXO[i], m$TTO[i]))
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

anova3_terminos <- function(X, y) {
  n <- length(y)
  ff <- ajustar(X, y, 1:4); df <- n - 4L; mse <- ff$sse / df
  out <- list()
  for (nq in list(c("SEXO", 2), c("TTO", 3), c("SEXO:TTO", 4))) {
    cols <- setdiff(1:4, as.integer(nq[[2]]))
    fr <- ajustar(X, y, cols)
    Fv <- (fr$sse - ff$sse) / mse
    out[[nq[[1]]]] <- c(Fv, f_sf(Fv, 1, df))
  }
  list(stats = out, df = df)
}

.vcov_hc3 <- function(X, resid) {
  n <- length(X)
  Xm <- matrix(unlist(X), n, 4L, byrow = TRUE)
  XtX <- t(Xm) %*% Xm
  XtXi <- invertir(XtX)
  h <- vapply(seq_len(n), function(i)
    as.numeric(Xm[i, , drop = FALSE] %*% XtXi %*% t(Xm[i, , drop = FALSE])), numeric(1))
  meat <- matrix(0, 4L, 4L)
  for (i in seq_len(n)) {
    w <- (resid[i]^2) / ((1 - h[i])^2)
    meat <- meat + w * (t(Xm[i, , drop = FALSE]) %*% Xm[i, , drop = FALSE])
  }
  XtXi %*% meat %*% XtXi
}

hc3_terminos <- function(X, y) {
  n <- length(y)
  ff <- ajustar(X, y, 1:4)
  V <- .vcov_hc3(X, ff$resid)
  df <- n - 4L
  out <- list()
  for (ni in list(c("SEXO", 2), c("TTO", 3), c("SEXO:TTO", 4))) {
    idx <- as.integer(ni[[2]])
    Fv <- ff$b[idx]^2 / V[idx, idx]
    out[[ni[[1]]]] <- c(Fv, f_sf(Fv, 1, df))
  }
  list(stats = out, df = df)
}

levene_bf <- function(m) {
  z <- c(); g <- c()
  celdas <- list()
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
  ff <- ajustar(X, y, 1:4)
  sw <- shapiro.test(ff$resid)
  c(as.numeric(sw$statistic), as.numeric(sw$p.value))
}

art_terminos <- function(m) {
  y <- m$y; n <- length(y); gran <- promedio(y)
  sm <- vapply(NIVELES_SEXO, function(l) promedio(m$y[m$SEXO == l]), numeric(1))
  tm <- vapply(NIVELES_TTO,  function(l) promedio(m$y[m$TTO == l]),  numeric(1))
  cm <- list()
  for (k in CELDAS_4) cm[[paste(k, collapse = "\r")]] <-
    promedio(m$y[m$SEXO == k[1] & m$TTO == k[2]])
  out <- list()
  for (term in c("SEXO", "TTO", "SEXO:TTO")) {
    alin <- numeric(n)
    for (i in seq_len(n)) {
      key <- paste(c(m$SEXO[i], m$TTO[i]), collapse = "\r")
      resid <- m$y[i] - cm[[key]]
      efecto <- if (term == "SEXO") sm[[m$SEXO[i]]] - gran
                else if (term == "TTO") tm[[m$TTO[i]]] - gran
                else gran - sm[[m$SEXO[i]]] - tm[[m$TTO[i]]] + cm[[key]]
      alin[i] <- round(resid + efecto, 8)
    }
    rr <- rangos_promedio(alin)
    mr <- data.frame(SEXO = m$SEXO, TTO = m$TTO, y = rr, stringsAsFactors = FALSE)
    d <- diseno(mr)
    st_ <- anova3_terminos(d$X, d$y)$stats
    out[[term]] <- st_[[term]]
  }
  list(stats = out, df = n - 4L)
}

artc_pares <- function(m) {
  y <- m$y; n <- length(y); gran <- promedio(y)
  rr <- rangos_promedio(round(y - gran, 8))
  idx <- lapply(CELDAS_4, function(k) which(m$SEXO == k[1] & m$TTO == k[2]))
  names(idx) <- vapply(CELDAS_4, paste, character(1), collapse = "\r")
  media <- vapply(idx, function(ii) promedio(rr[ii]), numeric(1))
  fitted <- vapply(seq_len(n), function(i)
    media[[paste(c(m$SEXO[i], m$TTO[i]), collapse = "\r")]], numeric(1))
  sse <- .suma((rr - fitted)^2); df <- n - 4L; mse <- sse / df
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

emmeans_pares <- function(m, robusto) {
  d <- diseno(m); X <- d$X; y <- d$y; n <- length(y)
  ff <- ajustar(X, y, 1:4); df <- n - 4L; mse <- ff$sse / df
  if (robusto) {
    V <- .vcov_hc3(X, ff$resid)
  } else {
    Xm <- matrix(unlist(X), n, 4L, byrow = TRUE)
    V <- mse * invertir(t(Xm) %*% Xm)
  }
  media <- list()
  for (k in CELDAS_4) media[[paste(k, collapse = "\r")]] <-
    promedio(m$y[m$SEXO == k[1] & m$TTO == k[2]])
  out <- list()
  for (par in PARES_D6) {
    xa <- fila_diseno(par[[2]][1], par[[2]][2])
    xb <- fila_diseno(par[[3]][1], par[[3]][2])
    cv <- xa - xb
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
bh <- function(pvals) {
  m <- length(pvals); if (m == 0L) return(numeric(0))
  ord <- order(pvals); aj <- numeric(m); prev <- 1
  for (rank in seq(m, 1L)) {
    i <- ord[rank]
    prev <- min(prev, pvals[i] * m / rank)
    aj[i] <- min(1, prev)
  }
  aj
}

fisher_2x2 <- function(a, b, cc, dd) {
  r1 <- a + b; r2 <- cc + dd; c1 <- a + cc; n <- a + b + cc + dd
  or_h <- ((a + 0.5) * (dd + 0.5)) / ((b + 0.5) * (cc + 0.5))
  if (r1 == 0 || r2 == 0 || c1 == 0 || (n - c1) == 0) return(list(or_h = or_h, p = 1))
  lp <- function(k)
    lgamma(r1 + 1) + lgamma(r2 + 1) + lgamma(c1 + 1) + lgamma(n - c1 + 1) -
    lgamma(n + 1) - lgamma(k + 1) - lgamma(r1 - k + 1) - lgamma(c1 - k + 1) -
    lgamma(r2 - c1 + k + 1)
  lo <- max(0, c1 - r2); hi <- min(r1, c1)
  p_obs <- exp(lp(a)); tot <- 0
  for (k in lo:hi) { pk <- exp(lp(k)); if (pk <= p_obs * (1 + 1e-7)) tot <- tot + pk }
  list(or_h = or_h, p = min(1, tot))
}

# ===========================================================================
# 1. Carga
# ===========================================================================
cargar <- function() {
  ruta <- file.path(RUTA_DATOS_PROC, "qpcr_cuantificacion_long.tsv")
  if (!file.exists(ruta))
    stop(sprintf("05: falta %s (correr 04_qpcr_cuantificacion primero)", ruta))
  d <- read.delim(ruta, sep = "\t", quote = "", stringsAsFactors = FALSE,
                  colClasses = "character", check.names = FALSE, encoding = "UTF-8",
                  na.strings = character(0))
  d$y <- ifelse(d$neg_ddCt == "", NA_real_, suppressWarnings(as.numeric(d$neg_ddCt)))
  d$detectado <- d$no_detectado != "TRUE"
  d
}

celdas_detectadas <- function(d, tej, gen) {
  vapply(CELDAS_4, function(k)
    sum(d$TEJIDO == tej & d$GEN == gen & d$SEXO == k[1] & d$TTO == k[2] &
        !is.na(d$y)), integer(1))
}

# ===========================================================================
# 2. Clasificacion y modelado por gen x tejido
# ===========================================================================
clasificar <- function(d) {
  registros <- list()
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    sub <- d[d$TEJIDO == tej & d$GEN == gen, , drop = FALSE]
    cuant <- nrow(sub) > 0 && sub$cuantificable[1] == "TRUE"
    ncel <- celdas_detectadas(d, tej, gen)
    det <- sub[!is.na(sub$y), , drop = FALSE]
    rec <- list(TEJIDO = tej, GEN = gen, es_transportador = gen %in% SET_TRANSP,
                n_celda = ncel, n_detectado_total = nrow(det),
                via = NA_character_, motivo = "",
                shapiro_p = NA_real_, levene_p = NA_real_, rama = "",
                stats = list(), df = NA_integer_, posthoc = NULL, bh = list())
    if (!cuant) {
      rec$via <- "D7_deteccion"
      rec$motivo <- "calibrador HEMBRA_CONTROL 0/detectados (D7)"
    } else if (min(ncel) < PISO_CELDA) {
      rec$via <- "descriptivo_n_bajo"
      rec$motivo <- sprintf("celda SEXO x TTO con < %d detectados", PISO_CELDA)
    } else {
      rec$via <- "modelo"
      rec <- ajustar_cascada(rec, det)
    }
    registros[[length(registros) + 1L]] <- rec
  }
  registros <- bh_por_tejido(registros)
  registros
}

ajustar_cascada <- function(rec, det) {
  m <- data.frame(SEXO = det$SEXO, TTO = det$TTO, y = det$y, stringsAsFactors = FALSE)
  d <- diseno(m)
  sw <- shapiro_residuos(d$X, d$y)
  lv <- levene_bf(m)
  rec$shapiro_p <- sw[2]; rec$levene_p <- lv[2]
  if (sw[2] >= ALFA && lv[2] >= ALFA) {
    rec$rama <- "anova3"; st_ <- anova3_terminos(d$X, d$y)
  } else if (sw[2] >= ALFA && lv[2] < ALFA) {
    rec$rama <- "hc3"; st_ <- hc3_terminos(d$X, d$y)
  } else {
    rec$rama <- "art"; st_ <- art_terminos(m)
  }
  rec$stats <- st_$stats; rec$df <- st_$df
  ip <- st_$stats[["SEXO:TTO"]][2]
  if (ip < ALFA) {
    if (rec$rama == "art") {
      pr <- artc_pares(m); metodo <- "ART-C"
    } else {
      pr <- emmeans_pares(m, robusto = (rec$rama == "hc3"))
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

bh_por_tejido <- function(registros) {
  for (tej in TEJIDOS_E15) {
    idx <- which(vapply(registros, function(r) r$TEJIDO == tej && r$via == "modelo",
                        logical(1)))
    if (!length(idx)) next
    for (term in c("SEXO", "TTO", "SEXO:TTO")) {
      praw <- vapply(idx, function(j) registros[[j]]$stats[[term]][2], numeric(1))
      adj <- bh(praw)
      for (k in seq_along(idx)) registros[[idx[k]]]$bh[[term]] <- adj[k]
    }
  }
  registros
}

# ===========================================================================
# 3. D7 -- deteccion (il6 @ BRAIN_E15)
# ===========================================================================
deteccion_d7 <- function(d, rec) {
  tej <- rec$TEJIDO; gen <- rec$GEN
  sub <- d[d$TEJIDO == tej & d$GEN == gen, , drop = FALSE]
  fisher_filas <- list(); tabla_filas <- list()
  for (sexo in NIVELES_SEXO) {
    gc <- sub[sub$SEXO == sexo & sub$TTO == "CONTROL", , drop = FALSE]
    gl <- sub[sub$SEXO == sexo & sub$TTO == "LPS", , drop = FALSE]
    a <- sum(gc$detectado); b <- nrow(gc) - a
    cc <- sum(gl$detectado); dd <- nrow(gl) - cc
    fr <- fisher_2x2(a, b, cc, dd)
    fisher_filas[[length(fisher_filas) + 1L]] <-
      list(tej, gen, sexo, nrow(gc), a, nrow(gl), cc, g10(fr$or_h), p6e(fr$p))
  }
  for (j in seq_along(CELDAS_4)) {
    k <- CELDAS_4[[j]]; grp <- GRUPOS_4[j]
    cel <- sub[sub$SEXO == k[1] & sub$TTO == k[2], , drop = FALSE]
    nd <- sum(cel$detectado)
    tabla_filas[[length(tabla_filas) + 1L]] <-
      list(tej, gen, grp, nrow(cel), nd, nrow(cel) - nd)
  }
  list(fisher = fisher_filas, tabla = tabla_filas)
}

# ===========================================================================
# 4. Cruza-verificacion en corrida contra car / emmeans / ARTool (solo R)
# ===========================================================================
verificar_contra_librerias <- function(d, registros) {
  op <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(op))
  peor <- 0
  for (rec in registros) {
    if (rec$via != "modelo") next
    sub <- d[d$TEJIDO == rec$TEJIDO & d$GEN == rec$GEN & !is.na(d$y), , drop = FALSE]
    s <- data.frame(y = sub$y,
                    SEXO = factor(sub$SEXO, levels = c("HEMBRA", "MACHO")),
                    TTO = factor(sub$TTO, levels = c("CONTROL", "LPS")))
    m <- lm(y ~ SEXO * TTO, data = s)
    # Levene
    lv <- car::leveneTest(y ~ SEXO * TTO, data = s)[1, "Pr(>F)"]
    peor <- max(peor, abs(lv - rec$levene_p))
    # omnibus por rama
    if (rec$rama == "anova3") {
      a <- car::Anova(m, type = 3)
      for (t in c("SEXO", "TTO", "SEXO:TTO"))
        peor <- max(peor, abs(a[t, "Pr(>F)"] - rec$stats[[t]][2]))
    } else if (rec$rama == "hc3") {
      a <- car::Anova(m, type = 3, white.adjust = "hc3")
      for (t in c("SEXO", "TTO", "SEXO:TTO"))
        peor <- max(peor, abs(a[t, "Pr(>F)"] - rec$stats[[t]][2]))
    } else {
      am <- ARTool::art(y ~ SEXO * TTO, data = s)
      aa <- anova(am); rownames(aa) <- aa$Term
      for (t in c("SEXO", "TTO", "SEXO:TTO"))
        peor <- max(peor, abs(aa[t, "Pr(>F)"] - rec$stats[[t]][2]))
    }
    # post hoc
    if (!is.null(rec$posthoc)) {
      mapc <- c("HEMBRA_CONTROL-HEMBRA_LPS"    = "HEMBRA,CONTROL - HEMBRA,LPS",
                "MACHO_CONTROL-MACHO_LPS"      = "MACHO,CONTROL - MACHO,LPS",
                "HEMBRA_LPS-MACHO_LPS"         = "HEMBRA,LPS - MACHO,LPS",
                "HEMBRA_CONTROL-MACHO_CONTROL" = "HEMBRA,CONTROL - MACHO,CONTROL")
      mape <- c("HEMBRA_CONTROL-HEMBRA_LPS"    = "HEMBRA CONTROL - HEMBRA LPS",
                "MACHO_CONTROL-MACHO_LPS"      = "MACHO CONTROL - MACHO LPS",
                "HEMBRA_LPS-MACHO_LPS"         = "HEMBRA LPS - MACHO LPS",
                "HEMBRA_CONTROL-MACHO_CONTROL" = "HEMBRA CONTROL - MACHO CONTROL")
      if (rec$rama == "art") {
        am <- ARTool::art(y ~ SEXO * TTO, data = s)
        cc <- summary(ARTool::art.con(am, "SEXO:TTO", adjust = "none"))
        pv <- setNames(cc$p.value, gsub("\\s+", " ", as.character(cc$contrast)))
        for (fila in rec$posthoc$filas) {
          e <- fila[[1]]; ref <- pv[[mapc[[e]]]]
          peor <- max(peor, abs(ref - fila[[5]]))
        }
      } else {
        vv <- if (rec$rama == "hc3") sandwich::vcovHC(m, type = "HC3") else vcov(m)
        emm <- emmeans::emmeans(m, ~ SEXO * TTO, vcov. = vv)
        pr <- summary(emmeans::contrast(emm, method = "pairwise", adjust = "none"))
        pv <- setNames(pr$p.value, gsub("\\s+", " ", as.character(pr$contrast)))
        for (fila in rec$posthoc$filas) {
          e <- fila[[1]]; ref <- pv[[mape[[e]]]]
          peor <- max(peor, abs(ref - fila[[5]]))
        }
      }
    }
  }
  peor
}

# ===========================================================================
# 5. Armado de tablas
# ===========================================================================
COLS_CLASIF <- c(
  "TEJIDO", "GEN", "es_transportador",
  "n_HEMBRA_CONTROL", "n_HEMBRA_LPS", "n_MACHO_CONTROL", "n_MACHO_LPS",
  "n_detectado_total", "via", "motivo",
  "shapiro_p_residuos", "levene_p", "rama_cascada",
  "F_SEXO", "p_SEXO", "F_TTO", "p_TTO", "F_SEXOxTTO", "p_SEXOxTTO",
  "p_SEXO_BH", "p_TTO_BH", "p_SEXOxTTO_BH",
  "interaccion_significativa", "post_hoc_corrido")

fila_clasif <- function(rec) {
  gv <- function(term, k) if (!is.null(rec$stats[[term]])) rec$stats[[term]][k] else NA_real_
  bv <- function(term) if (!is.null(rec$bh[[term]])) rec$bh[[term]] else NA_real_
  ip <- if (!is.null(rec$stats[["SEXO:TTO"]])) rec$stats[["SEXO:TTO"]][2] else NA_real_
  sig <- !is.na(ip) && ip < ALFA
  list(
    rec$TEJIDO, rec$GEN, rec$es_transportador,
    rec$n_celda[1], rec$n_celda[2], rec$n_celda[3], rec$n_celda[4],
    rec$n_detectado_total, rec$via, rec$motivo,
    p6e(rec$shapiro_p), p6e(rec$levene_p), rec$rama,
    p6e(gv("SEXO", 1)), p6e(gv("SEXO", 2)),
    p6e(gv("TTO", 1)), p6e(gv("TTO", 2)),
    p6e(gv("SEXO:TTO", 1)), p6e(gv("SEXO:TTO", 2)),
    p6e(bv("SEXO")), p6e(bv("TTO")), p6e(bv("SEXO:TTO")),
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

COLS_NOMODELO <- c("TEJIDO", "GEN", "via", "GRUPO", "n_detectado",
                   "mean_neg_ddCt", "sd_neg_ddCt", "median_neg_ddCt")

filas_nomodelo <- function(d, rec) {
  lapply(seq_along(CELDAS_4), function(j) {
    k <- CELDAS_4[[j]]; grp <- GRUPOS_4[j]
    vals <- d$y[d$TEJIDO == rec$TEJIDO & d$GEN == rec$GEN &
                d$SEXO == k[1] & d$TTO == k[2] & !is.na(d$y)]
    n <- length(vals)
    mu <- if (n) promedio(vals) else NA_real_
    sd <- if (n >= 2) desvio_muestral(vals, mu) else NA_real_
    md <- if (n) mediana(vals) else NA_real_
    list(rec$TEJIDO, rec$GEN, rec$via, grp, n,
         if (is.na(mu)) "" else g10(mu),
         if (is.na(sd)) "" else g10(sd),
         if (is.na(md)) "" else g10(md))
  })
}

# ===========================================================================
# 6. Reporte legible
# ===========================================================================
.md <- function(header, filas) {
  l1 <- paste0("| ", paste(header, collapse = " | "), " |")
  l2 <- paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |")
  cuerpo <- vapply(filas, function(f)
    paste0("| ", paste(vapply(f, .fmt, character(1)), collapse = " | "), " |"),
    character(1))
  paste(c(l1, l2, cuerpo), collapse = "\n")
}

construir_reporte <- function(fuente, clasif, posthoc, fisher_filas, tabla_filas,
                              resumen, d7_nota) {
  L <- c(
    "# Reporte de modelado qPCR (T5)", "",
    "Generado por `05_qpcr_modelos` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. Metodo", "",
    BULLETS_METODO,
    "", "## 2. Clasificacion por gen x tejido", "",
    .md(COLS_CLASIF, clasif), "",
    "## 3. Post hoc (D6) -- comparaciones con interaccion significativa", "")
  L <- c(L, if (length(posthoc)) .md(COLS_POSTHOC, posthoc)
         else "Ninguna interaccion SEXO x TTO alcanzo p < alfa: no se corrio post hoc.")
  L <- c(L, "", "## 4. D7 -- il6 @ BRAIN_E15: proporcion de deteccion", "",
    "Fisher exacto 2x2 (deteccion Control vs LPS) dentro de cada sexo:", "",
    .md(c("TEJIDO", "GEN", "SEXO", "n_CONTROL", "det_CONTROL", "n_LPS",
          "det_LPS", "OR_Haldane", "p_fisher"), fisher_filas), "",
    "Tabla 2x4 (grupo x deteccion):", "",
    .md(c("TEJIDO", "GEN", "GRUPO", "n_total", "n_detectado", "n_no_detectado"),
        tabla_filas), "",
    d7_nota, "",
    "## 5. Resumen", "")
  for (linea in resumen) L <- c(L, paste0("- ", linea))
  L <- c(L, "", "## 6. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `05_qpcr_modelos`: eleccion de ",
           "rama por gen x tejido, il6R @ BRAIN_E15 fuera del modelo por piso de ",
           "celda, il6 @ BRAIN_E15 degenerado (0 detectados en los 4 grupos sobre ",
           "datos reales), y la nota de paridad (Shapiro-Wilk por libreria)."), "")
  paste(L, collapse = "\n")
}

# ===========================================================================
# 7. Artefactos compartidos (merge por 'script') -- headers identicos a 02/03/04.
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
  header <- c("id", "tipo", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script")
  merge_por_script(file.path(RUTA_TABLAS, "verificaciones.csv"), header, filas_nuevas,
                   "script", function(fs) order(
                     vapply(fs, `[[`, character(1), 7), vapply(fs, `[[`, character(1), 1),
                     method = "radix"))
}

# Bloques de prosa -- literales IDENTICOS a los de python/05_qpcr_modelos.py
# (paridad byte a byte de los .md). Una linea por bullet.
BULLETS_METODO <- c(
paste0("- Modelo por gen x tejido: `neg_ddCt ~ SEXO * TTO` (OLS, contrastes ",
       "suma-cero). alfa = 0.05 para la interaccion (gate D6) y para los pretests ",
       "de supuestos."),
paste0("- Cascada D5 / 4.1: Shapiro-Wilk sobre los **residuos del modelo ",
       "conjunto**; Levene (Brown-Forsythe, centro = mediana) sobre las 4 celdas. ",
       "Shapiro>=.05 & Levene>=.05 -> `anova3` (ANOVA III); Shapiro>=.05 & ",
       "Levene<.05 -> `hc3` (Wald III con sandwich HC3); Shapiro<.05 -> `art` ",
       "(Aligned Rank Transform)."),
paste0("- Piso de celda: si alguna celda SEXO x TTO tiene < 5 detectados, no se ",
       "modela (`descriptivo_n_bajo`)."),
paste0("- D7: gen x tejido con calibrador HEMBRA_CONTROL 0/detectados ",
       "(`cuantificable = FALSE`) fuera del modelo -> Fisher exacto 2x2 de ",
       "deteccion Control vs LPS dentro de cada sexo + tabla 2x4."),
paste0("- Post hoc D6 (solo si interaccion p < alfa): 4 comparaciones fijas, ",
       "correccion Holm. `anova3` -> contraste de medias marginales (vcov OLS); ",
       "`hc3` -> idem con vcov HC3; `art` -> ART-C (Elkin et al. 2021)."),
paste0("- D12: columna suplementaria `*_BH` = p ajustado por Benjamini-Hochberg ",
       "dentro de cada tejido entre los genes modelados, por termino. No dirige ",
       "la inferencia.")
)

bloque_descartes <- function(resumen_rama, resumen_bh) {
  paste(c(
"## 05_qpcr_modelos",
"",
"### D13 -- MADRE (camada) no se incluye en el modelo",
"",
paste0("- Se evaluo incluir `MADRE` como efecto aleatorio (18 madres, 2 fetos ",
       "por camada -- 1 hembra + 1 macho). **No se incorporo**: ",
       "<<< COMPLETAR: evidencia de analisis previos >>>. Se asume ",
       "independencia entre fetos para el analisis (D13, AGENTS.md); la ",
       "limitacion se declara en el informe (Seccion 7)."),
"",
"### Eleccion de rama de la cascada D5 por gen x tejido",
"",
paste0("- Se ajusta `neg_ddCt ~ SEXO * TTO` (OLS, contr.sum) y se decide la rama ",
       "con Shapiro-Wilk sobre los **residuos del modelo conjunto** del gen x ",
       "tejido y Levene (Brown-Forsythe) sobre las 4 celdas, ambos a alfa = 0.05. ",
       sprintf("Reparto sobre datos reales: %s. La rama usada queda en ", resumen_rama),
       "`rama_cascada` de `qpcr_modelos_clasificacion.csv` y en `procedencia.csv`."),
paste0("- La rama `hc3` reporta un Wald tipo III (1 gl) con vcov sandwich HC3; su ",
       "post hoc usa el mismo vcov HC3. La rama `art` reporta ART (Wobbrock et al. ",
       "2011) y su post hoc es ART-C (Elkin et al. 2021), NUNCA emmeans directo ",
       "sobre el modelo ART."),
"",
"### il6R @ BRAIN_E15: fuera del modelo por piso de celda",
"",
paste0("- Tras descartar no detectados quedan 3/7/3/3 valores en las celdas ",
       "HEMBRA_CONTROL / HEMBRA_LPS / MACHO_CONTROL / MACHO_LPS. Con celdas de 3 el ",
       "modelo factorial 2x2 y sus pretests de supuestos no son defendibles: via ",
       "`descriptivo_n_bajo` (solo n y descriptivo de neg_ddCt en ",
       "`qpcr_modelos_nomodelo_descriptivo.csv`). El calibrador de il6R @ BRAIN_E15 ",
       "SI tiene detectados (3/9), asi que no es un caso D7."),
"",
"### il6 @ BRAIN_E15: D7 -- solo proporcion de deteccion",
"",
paste0("- Calibrador HEMBRA_CONTROL 0/9 detectados -> `cuantificable = FALSE` en ",
       "04_qpcr_cuantificacion -> via `D7_deteccion` (D7): sin dCt de calibrador no ",
       "hay ddCt, asi que il6 en cerebro no entra al modelo y se analiza solo como ",
       "proporcion de deteccion. Se emite el Fisher exacto 2x2 (deteccion Control ",
       "vs LPS) dentro de cada sexo y la tabla 2x4 grupo x deteccion ",
       "(`qpcr_il6_brain_fisher.csv`, `qpcr_il6_brain_tabla2x4.csv`); el reporte ",
       "dice si algun Fisher alcanza p < alfa. Es el unico gen x tejido via D7; ",
       "se detecta programaticamente (calibrador HEMBRA_CONTROL 0/detectados)."),
"",
"### Correccion entre genes (D12) -- columna suplementaria",
"",
paste0("- Se agrega `p_SEXO_BH`, `p_TTO_BH`, `p_SEXOxTTO_BH` = p ajustado por ",
       "Benjamini-Hochberg dentro de cada tejido a traves de los genes modelados, ",
       "por termino. Es suplementaria: no cambia ninguna conclusion. ", resumen_bh),
"",
"### Paridad R / Python",
"",
paste0("- Todo el nucleo numerico (OLS 4x4 por eliminacion de Gauss, SS tipo III ",
       "por comparacion de modelos, sandwich HC3, Levene, ART y ART-C, Holm, BH) es ",
       "codigo PROPIO identico en R y Python; R lo cruza-verifica en corrida contra ",
       "`car::Anova`, `emmeans` y `ARTool` (`stopifnot`, tol 1e-6)."),
paste0("- Unica excepcion: **Shapiro-Wilk** se toma de libreria (`shapiro.test` en ",
       "R, `scipy.stats.shapiro` en Python). Coinciden a ~1e-12; se comparan dentro ",
       "de `TOL_ESTADISTICO` en T10, no por identidad de bytes. Los p / ",
       "estadisticos que pasan por funciones trascendentes (pf, pt, lgamma, W) se ",
       "guardan como texto `%.6e` (igual que 03_elisa).")
  ), collapse = "\n")
}

actualizar_descartados <- function(resumen_rama, resumen_bh) {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 05_qpcr_modelos:inicio -->"
  marca_fin <- "<!-- 05_qpcr_modelos:fin -->"
  nuevo <- paste0(marca_ini, "\n", bloque_descartes(resumen_rama, resumen_bh),
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
  fuente <- fuente_datos(ARCHIVO_QPCR)
  registros <- clasificar(d)

  ord <- order(vapply(registros, function(r) match(r$TEJIDO, TEJIDOS_E15), integer(1)),
               vapply(registros, function(r) match(r$GEN, GENES), integer(1)))
  registros <- registros[ord]

  peor <- verificar_contra_librerias(d, registros)
  cat(sprintf("  [cruza-verificacion R vs car/emmeans/ARTool] peor |dif| = %.3e\n", peor))
  stopifnot(peor < 1e-6)

  clasif <- lapply(registros, fila_clasif)
  posthoc <- do.call(c, lapply(registros, filas_posthoc))
  if (is.null(posthoc)) posthoc <- list()
  nomodelo <- do.call(c, lapply(Filter(function(r) r$via == "descriptivo_n_bajo",
                                       registros), function(r) filas_nomodelo(d, r)))
  if (is.null(nomodelo)) nomodelo <- list()

  fisher_filas <- list(); tabla_filas <- list()
  for (rec in registros) if (rec$via == "D7_deteccion") {
    dd <- deteccion_d7(d, rec)
    fisher_filas <- c(fisher_filas, dd$fisher)
    tabla_filas <- c(tabla_filas, dd$tabla)
  }

  via_cont <- c(modelo = 0L, D7_deteccion = 0L, descriptivo_n_bajo = 0L)
  for (r in registros) via_cont[r$via] <- via_cont[r$via] + 1L
  rama_cont <- c(anova3 = 0L, hc3 = 0L, art = 0L)
  for (r in registros) if (nzchar(r$rama)) rama_cont[r$rama] <- rama_cont[r$rama] + 1L
  n_posthoc <- sum(vapply(registros, function(r) !is.null(r$posthoc), logical(1)))
  resumen_rama <- paste(sprintf("%s=%d", names(rama_cont), rama_cont), collapse = "; ")

  bh_txt <- character(0)
  for (tej in TEJIDOS_E15) {
    mods <- Filter(function(r) r$TEJIDO == tej && r$via == "modelo", registros)
    crudas <- Filter(function(r) r$stats[["SEXO:TTO"]][2] < ALFA, mods)
    sobre <- Filter(function(r) r$bh[["SEXO:TTO"]] < ALFA, crudas)
    bh_txt <- c(bh_txt, sprintf(paste0("%s: %d/%d interacciones significativas ",
                                       "crudas sobreviven BH"),
                                tej, length(sobre), length(crudas)))
  }
  resumen_bh <- paste0(paste(bh_txt, collapse = "; "), ".")

  resumen <- c(
    sprintf(paste0("gen x tejido: %d modelados, %d por D7 (deteccion), %d ",
                   "descriptivo por piso de celda."),
            via_cont[["modelo"]], via_cont[["D7_deteccion"]],
            via_cont[["descriptivo_n_bajo"]]),
    sprintf("Reparto de ramas D5: %s.", resumen_rama),
    sprintf("Post hoc D6 corrido en %d gen x tejido (interaccion p < %s).",
            n_posthoc, format(ALFA)),
    resumen_bh)

  if (length(fisher_filas)) {
    n_fis <- length(fisher_filas)
    k_fis <- sum(vapply(fisher_filas, function(f)
      nzchar(f[[9]]) && as.numeric(f[[9]]) < ALFA, logical(1)))
    d7_nota <- sprintf(paste0("De los %d Fisher exactos por sexo, %d alcanzan ",
                              "p < %s (Control vs LPS, deteccion de il6 en cerebro)."),
                       n_fis, k_fis, format(ALFA))
  } else {
    d7_nota <- "No hay gen x tejido via D7 en esta corrida."
  }

  escribir_tsv(file.path(RUTA_DATOS_PROC, "qpcr_modelos_clasificacion.tsv"),
               COLS_CLASIF, clasif)
  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "qpcr_modelos_clasificacion.csv"), COLS_CLASIF, clasif)
    escribir_csv(file.path(base, "qpcr_modelos_posthoc.csv"), COLS_POSTHOC, posthoc)
    escribir_csv(file.path(base, "qpcr_modelos_nomodelo_descriptivo.csv"),
                 COLS_NOMODELO, nomodelo)
    escribir_csv(file.path(base, "qpcr_il6_brain_fisher.csv"),
                 c("TEJIDO", "GEN", "SEXO", "n_CONTROL", "det_CONTROL",
                   "n_LPS", "det_LPS", "OR_Haldane", "p_fisher"), fisher_filas)
    escribir_csv(file.path(base, "qpcr_il6_brain_tabla2x4.csv"),
                 c("TEJIDO", "GEN", "GRUPO", "n_total", "n_detectado",
                   "n_no_detectado"), tabla_filas)
  }

  escribir_lineas(file.path(RUTA_TABLAS, "qpcr_modelos_reporte.md"),
                  construir_reporte(fuente, clasif, posthoc, fisher_filas,
                                    tabla_filas, resumen, d7_nota))
  actualizar_descartados(resumen_rama, resumen_bh)

  ent <- sprintf("data/processed/qpcr_cuantificacion_long.tsv (de data/%s/%s)",
                 fuente, ARCHIVO_QPCR)
  registrar_procedencia(list(
    list("data/processed/qpcr_modelos_clasificacion.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("clasificacion por gen x tejido: via, rama D5, F/p ",
         "de SEXO/TTO/SEXOxTTO, columna BH (D12), flag de post hoc")),
    list("outputs/tables/{R,python}/qpcr_modelos_clasificacion.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "idem, copia por implementacion"),
    list("outputs/tables/{R,python}/qpcr_modelos_posthoc.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("post hoc D6: 4 comparaciones fijas con Holm, por ",
         "rama (contraste marginal OLS/HC3 o ART-C)")),
    list("outputs/tables/{R,python}/qpcr_modelos_nomodelo_descriptivo.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "n y descriptivo de neg_ddCt por grupo para gen x tejido no modelados"),
    list("outputs/tables/{R,python}/qpcr_il6_brain_fisher.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "D7: Fisher exacto 2x2 de deteccion Control vs LPS por sexo"),
    list("outputs/tables/{R,python}/qpcr_il6_brain_tabla2x4.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "D7: tabla 2x4 grupo x deteccion para il6 @ BRAIN_E15"),
    list("outputs/tables/qpcr_modelos_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del modelado qPCR (T5)")
  ))
  ok_via <- via_cont[["modelo"]] == 18L && via_cont[["D7_deteccion"]] == 1L &&
            via_cont[["descriptivo_n_bajo"]] == 1L
  registrar_verificaciones(list(
    list("modelos_via_conteo", "recalculo",
         "20 gen x tejido repartidos en via modelo / D7 / descriptivo",
         sprintf("modelo=%d;D7_deteccion=%d;descriptivo_n_bajo=%d",
                 via_cont[["modelo"]], via_cont[["D7_deteccion"]],
                 via_cont[["descriptivo_n_bajo"]]),
         "modelo=18;D7_deteccion=1;descriptivo_n_bajo=1",
         if (ok_via) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("modelos_rama_conteo", "recalculo",
         "cada gen x tejido modelado usa la rama que corresponde a sus propios Shapiro/Levene",
         resumen_rama,
         "rama = art si Shapiro<alfa; hc3 si Shapiro>=alfa & Levene<alfa; anova3 si ambos >=alfa",
         if (!length(Filter(function(r) r$via == "modelo" &&
             r$rama != (if (r$shapiro_p < ALFA) "art"
                        else if (r$levene_p < ALFA) "hc3" else "anova3"),
             registros))) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("modelos_alfa", "declaracion", "alfa para gate D6 y pretests de supuestos",
         g10(ALFA), "0.05", if (ALFA == 0.05) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("modelos_piso_celda", "declaracion",
         "minimo de detectados por celda SEXO x TTO para modelar",
         as.character(PISO_CELDA), "5", if (PISO_CELDA == 5L) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("modelos_posthoc_holm", "declaracion",
         "post hoc D6: 4 comparaciones fijas con correccion Holm",
         sprintf("%d gen x tejido; 4 comparaciones c/u; Holm", n_posthoc),
         "solo si interaccion p < alfa", "TRUE", ESTE_SCRIPT),
    list("modelos_art_core_vs_artool", "declaracion",
         paste0("ART y ART-C PROPIOS; R cruza-verifica en corrida contra ARTool ",
                "(stopifnot 1e-6)"),
         "nucleo PROPIO identico R/Python; Shapiro-Wilk por libreria",
         "verificado", "TRUE", ESTE_SCRIPT),
    list("modelos_bh_d12", "declaracion", "columna suplementaria BH entre genes por tejido (D12)",
         resumen_bh, "no cambia conclusiones", "TRUE", ESTE_SCRIPT),
    list("modelos_d7_il6_brain", "recalculo",
         "il6 @ BRAIN_E15 fuera del modelo (D7), solo deteccion",
         sprintf("%d gen x tejido via D7", via_cont[["D7_deteccion"]]),
         "il6@BRAIN_E15", if (via_cont[["D7_deteccion"]] == 1L) "TRUE" else "FALSE",
         ESTE_SCRIPT)
  ))

  cat("== 05_qpcr_modelos.R ==\n")
  cat(sprintf("  fuente qPCR = %s\n", fuente))
  cat(sprintf("  via: modelo=%d  D7=%d  descriptivo_n_bajo=%d\n",
              via_cont[["modelo"]], via_cont[["D7_deteccion"]],
              via_cont[["descriptivo_n_bajo"]]))
  cat(sprintf("  ramas D5: %s\n", resumen_rama))
  cat(sprintf("  post hoc D6 en %d gen x tejido\n", n_posthoc))
  for (r in registros) {
    if (r$via != "modelo") {
      cat(sprintf("    [%-18s] %-12s %-9s %s\n", r$via, r$TEJIDO, r$GEN, r$motivo))
    } else {
      ip <- r$stats[["SEXO:TTO"]][2]
      marca <- if (ip < ALFA) " *interaccion*" else ""
      cat(sprintf("    [%-6s] %-12s %-9s pS=%s pT=%s pSxT=%s%s\n",
                  r$rama, r$TEJIDO, r$GEN, p6e(r$stats[["SEXO"]][2]),
                  p6e(r$stats[["TTO"]][2]), p6e(ip), marca))
    }
  }
  cat("  -> data/processed/qpcr_modelos_clasificacion.tsv,\n")
  cat("     outputs/tables/{R,python}/qpcr_modelos_*.csv, qpcr_il6_brain_*.csv,\n")
  cat("     outputs/tables/qpcr_modelos_reporte.md\n")
}

if (sys.nframe() == 0L) main()
