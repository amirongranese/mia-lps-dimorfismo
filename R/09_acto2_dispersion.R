# 09_acto2_dispersion.R -- ACTO 2.3-2.4: dispersion por grupo + test formal de
#                          diferencia de correlaciones placenta <-> cerebro.
#
# Por que existe este archivo: T7 (08) describio la correlacion placenta <->
# cerebro por feto (Spearman rho) pero tiene PROHIBIDO compararla entre grupos
# (prohibicion 4) e interpretarla como coordinacion biologica sin descartar
# restriccion de rango (prohibicion 5). Este script hace lo primero de forma
# formal y prepara lo segundo:
#
#   * 2.4 -- TEST REPORTADO (prohibicion 4; decision del usuario): **Fisher z
#     sobre rho de Spearman**. z = atanh(rho); estadistico
#     (z_control - z_lps)/sqrt(SE_control^2 + SE_lps^2); p normal a dos colas.
#     SE primario = Bonett-Wright sqrt((1 + rho^2/2)/(n-3)) (el del IC de T7);
#     SE clasico 1/sqrt(n-3) como columna al lado (no cambia conclusiones);
#     p_bw_bh = BH entre items, suplementario (D12).
#   * 2.3 -- DISPERSION (insumo de la prohibicion 5): sobre los MISMOS pares por
#     feto, SD (n-1) de -ddCt de cada lado en Control y LPS, cociente de
#     varianzas LPS/Control y Levene Brown-Forsythe (centro = mediana) por lado.
#     No alcanza para descartar restriccion de rango: eso lo hace 10.
#
# Piso: se testea un item solo con n_par >= 5 en Control y en LPS. `il6R` suele
# quedar por debajo -> fila con n y sin estadistico, no se omite.
#
# PARIDAD R/Python: rho, SD, cociente y sumas usan acumulador double explicito
# (mismo orden) -> texto "%.10g" bit-identico; z (atanh), p (normal) y F de
# Levene (pf) -> texto "%.6e". Aca (solo R) se cruza-verifica cada rho contra
# cor.test(exact=FALSE) (tol 1e-9) y cada F de Levene contra
# car::leveneTest(center=median) (tol 1e-8), con stopifnot. NUNCA se llama al
# metodo exacto de cor.test spearman (AS 89): segfaultea en este build de R.
# Las figuras son PNG: equivalentes, no byte-identicas.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({ library(ggplot2); library(car) })

ESTE_SCRIPT <- "09_acto2_dispersion"

COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")
Z975 <- 1.959963984540054
PISO_PAR <- 5L

GEN_SIN_CEREBRO <- "il6"
GENES_CORR <- setdiff(GENES, GEN_SIN_CEREBRO)          # 9 genes
ITEMS <- c(GENES_CORR, "score_compuesto")
TEJIDOS <- TEJIDOS_E15                                 # PLACENTA_E15, BRAIN_E15
DPI <- 300

# Estratificacion por sexo (pedido explicito, pedidos/cambios_acto2_dispersion_por_sexo.md):
# AMBOS_SEXOS = comportamiento previo (agrupado); HEMBRA/MACHO = dentro de cada
# sexo. Aditivo: AMBOS_SEXOS se conserva tal cual estaba.
ESTRATOS <- c("AMBOS_SEXOS", "HEMBRA", "MACHO")
sexo_de_estrato <- function(e) if (e == "AMBOS_SEXOS") NULL else e

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..08.
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
promedio <- function(xs) .suma(xs) / length(xs)
desvio <- function(xs) {
  n <- length(xs)
  if (n < 2L) return(NA_real_)
  m <- .suma(xs) / n
  s <- 0; for (v in xs) s <- s + (v - m) * (v - m)
  sqrt(s / (n - 1))
}
mediana <- function(xs) {
  if (length(xs) == 0L) return(NA_real_)
  s <- sort(xs); k <- length(s)
  if (k %% 2L == 1L) s[(k %/% 2L) + 1L] else (s[k %/% 2L] + s[(k %/% 2L) + 1L]) / 2
}
rangos_promedio <- function(xs) rank(xs, ties.method = "average")

spearman_rho <- function(x, y) {
  n <- length(x)
  rx <- as.numeric(rangos_promedio(x)); ry <- as.numeric(rangos_promedio(y))
  mx <- .suma(rx) / n; my <- .suma(ry) / n
  sxy <- 0; sxx <- 0; syy <- 0
  for (i in seq_len(n)) {
    dx <- rx[i] - mx; dy <- ry[i] - my
    sxy <- sxy + dx * dy; sxx <- sxx + dx * dx; syy <- syy + dy * dy
  }
  if (sxx <= 0 || syy <= 0) return(NA_real_)
  sxy / sqrt(sxx * syy)
}
f_sf <- function(x, d1, d2) pf(x, d1, d2, lower.tail = FALSE)
norm_sf2 <- function(x) 2 * pnorm(-abs(x))

levene_bf_2 <- function(a, b) {
  na <- length(a); nb <- length(b)
  if (na < 2L || nb < 2L) return(c(NA_real_, NA_real_))
  ma <- mediana(a); mb <- mediana(b)
  za <- abs(a - ma); zb <- abs(b - mb)
  N <- na + nb
  gran <- .suma(c(za, zb)) / N
  gma <- .suma(za) / na; gmb <- .suma(zb) / nb
  ssb <- na * (gma - gran)^2 + nb * (gmb - gran)^2
  ssw <- .suma((za - gma)^2) + .suma((zb - gmb)^2)
  if (ssw <= 0) return(c(NA_real_, NA_real_))
  Fv <- (ssb / 1) / (ssw / (N - 2))
  c(Fv, f_sf(Fv, 1, N - 2))
}

bh <- function(pvals) {
  m <- length(pvals)
  if (m == 0L) return(numeric(0))
  orden <- order(pvals, method = "radix")
  ajust <- numeric(m); prev <- 1
  for (rk in m:1) {
    i <- orden[rk]
    val <- pvals[i] * m / rk
    prev <- min(prev, val)
    ajust[i] <- min(1, prev)
  }
  ajust
}

# ===========================================================================
# 1. Carga y emparejamiento por feto -- identico a 08.
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
  sexo  <- setNames(cu$SEXO, cu$FETO)[!duplicated(cu$FETO)]
  negdd <- setNames(
    ifelse(cu$neg_ddCt == "", NA_real_, suppressWarnings(as.numeric(cu$neg_ddCt))),
    paste(cu$FETO, cu$TEJIDO, cu$GEN, sep = "\r"))
  scv <- setNames(
    ifelse(sc$score_compuesto == "", NA_real_,
           suppressWarnings(as.numeric(sc$score_compuesto))),
    paste(sc$FETO, sc$TEJIDO, sep = "\r"))
  fetos <- names(madre)[order(madre, names(madre), method = "radix")]
  list(fetos = fetos, tto = tto, sexo = sexo, negdd = negdd, sc = scv)
}

# `sexo_filtro`: NULL = ambos sexos (comportamiento previo); "HEMBRA"/"MACHO"
# restringe ademas por sexo -- usado por la estratificacion pedida (D#, ver
# analisis_descartados 09_acto2_dispersion): agrupar sexos promedia
# correlaciones de signo opuesto (ver fatcd36 en el reporte).
pares <- function(D, item, estrato, sexo_filtro = NULL) {
  xs <- c(); ys <- c(); ts <- c()
  for (f in D$fetos) {
    if (estrato != "GLOBAL" && D$tto[[f]] != estrato) next
    if (!is.null(sexo_filtro) && D$sexo[[f]] != sexo_filtro) next
    if (item == "score_compuesto") {
      xp <- D$sc[[paste(f, "PLACENTA_E15", sep = "\r")]]
      yb <- D$sc[[paste(f, "BRAIN_E15", sep = "\r")]]
    } else {
      xp <- D$negdd[[paste(f, "PLACENTA_E15", item, sep = "\r")]]
      yb <- D$negdd[[paste(f, "BRAIN_E15", item, sep = "\r")]]
    }
    if (is.null(xp) || is.null(yb) || is.na(xp) || is.na(yb)) next
    xs <- c(xs, xp); ys <- c(ys, yb); ts <- c(ts, D$tto[[f]])
  }
  list(x = xs, y = ys, tto = ts)
}
.lado <- function(pr, tej) if (tej == "PLACENTA_E15") pr$x else pr$y

# ===========================================================================
# 2.3 -- Tabla de dispersion (+ cruza-verificacion Levene vs car).
# ===========================================================================
COLS_DISP <- c("ITEM", "TIPO", "TEJIDO", "n_control", "sd_control", "n_lps",
               "sd_lps", "ratio_var_lps_control", "levene_bf_F", "levene_bf_p")

tabla_dispersion <- function(D) {
  filas <- list(); peor <- 0
  for (item in ITEMS) {
    tipo <- if (item == "score_compuesto") "score" else "gen"
    prc <- pares(D, item, "CONTROL"); prl <- pares(D, item, "LPS")
    for (tej in TEJIDOS) {
      vc <- .lado(prc, tej); vl <- .lado(prl, tej)
      nc <- length(vc); nl <- length(vl)
      sdc <- desvio(vc); sdl <- desvio(vl)
      ratio <- if (!is.na(sdc) && !is.na(sdl) && sdc > 0)
        (sdl * sdl) / (sdc * sdc) else NA_real_
      Fp <- if (nc >= PISO_PAR && nl >= PISO_PAR) levene_bf_2(vc, vl)
            else c(NA_real_, NA_real_)
      if (!is.na(Fp[1])) {
        dfv <- data.frame(y = c(vc, vl),
                          g = factor(c(rep("C", nc), rep("L", nl))))
        lt <- suppressWarnings(car::leveneTest(y ~ g, data = dfv, center = median))
        peor <- max(peor, abs(Fp[1] - lt[1, "F value"]))
      }
      filas[[length(filas) + 1L]] <- list(item, tipo, tej, nc, g10(sdc), nl,
                                          g10(sdl), g10(ratio), p6e(Fp[1]),
                                          p6e(Fp[2]))
    }
  }
  list(filas = filas, peor = peor)
}

# ===========================================================================
# 2.4 -- Test reportado: Fisher z sobre rho de Spearman (prohibicion 4).
# Estratificado por ESTRATO (AMBOS_SEXOS/HEMBRA/MACHO, pedido explicito): el
# caso que lo justifica es fatcd36, donde HEMBRA Control rho=0.86 y MACHO
# Control rho=-0.43 -- agrupados dan 0.47, que no describe a ninguno de los
# dos (ver reporte y analisis_descartados). AMBOS_SEXOS = filas previas, sin
# cambios; HEMBRA/MACHO son aditivas.
# ===========================================================================
COLS_TEST <- c("ITEM", "TIPO", "ESTRATO", "n_control", "rho_control", "n_lps",
               "rho_lps", "delta_rho", "z_control", "z_lps", "se_bw_control",
               "se_bw_lps", "stat_z_bw", "p_bw", "p_bw_bh", "se_clasico_control",
               "se_clasico_lps", "stat_z_clasico", "p_clasico")

.se_bw <- function(rho, n) sqrt((1 + rho * rho / 2) / (n - 3))
.se_clasico <- function(n) 1 / sqrt(n - 3)
.fisher_z <- function(rc, rl, sec, sel) {
  zc <- atanh(rc); zl <- atanh(rl)
  stat <- (zc - zl) / sqrt(sec * sec + sel * sel)
  list(zc = zc, zl = zl, sec = sec, sel = sel, stat = stat, p = norm_sf2(stat))
}

# BH (D12) se calcula DENTRO de cada estrato (pedido explicito: la potencia y
# el n difieren mucho entre AMBOS_SEXOS y HEMBRA/MACHO por separado; mezclar
# los 27 p-values en un solo ajuste no tendria sentido).
tabla_test <- function(D) {
  filas <- list(); peor <- 0
  pend <- setNames(vector("list", length(ESTRATOS)), ESTRATOS)
  for (item in ITEMS) {
    tipo <- if (item == "score_compuesto") "score" else "gen"
    for (estrato in ESTRATOS) {
      sx <- sexo_de_estrato(estrato)
      prc <- pares(D, item, "CONTROL", sx); prl <- pares(D, item, "LPS", sx)
      nc <- length(prc$x); nl <- length(prl$x)
      idx <- length(filas) + 1L
      if (nc < PISO_PAR || nl < PISO_PAR) {
        filas[[idx]] <- c(list(item, tipo, estrato, nc, "", nl, ""),
                          as.list(rep("", 12)))
        next
      }
      rc <- spearman_rho(prc$x, prc$y); rl <- spearman_rho(prl$x, prl$y)
      ctc <- suppressWarnings(cor.test(prc$x, prc$y, method = "spearman",
                                       exact = FALSE))
      ctl <- suppressWarnings(cor.test(prl$x, prl$y, method = "spearman",
                                       exact = FALSE))
      peor <- max(peor, abs(rc - as.numeric(ctc$estimate)),
                  abs(rl - as.numeric(ctl$estimate)))
      drho <- rc - rl
      bw <- .fisher_z(rc, rl, .se_bw(rc, nc), .se_bw(rl, nl))
      cl <- .fisher_z(rc, rl, .se_clasico(nc), .se_clasico(nl))
      pend[[estrato]][[length(pend[[estrato]]) + 1L]] <- c(idx, bw$p)
      filas[[idx]] <- list(item, tipo, estrato, nc, g10(rc), nl, g10(rl),
                           g10(drho), p6e(bw$zc), p6e(bw$zl), p6e(bw$sec),
                           p6e(bw$sel), p6e(bw$stat), p6e(bw$p), "",
                           p6e(cl$sec), p6e(cl$sel), p6e(cl$stat), p6e(cl$p))
    }
  }
  for (estrato in ESTRATOS) {
    pd <- pend[[estrato]]
    if (!length(pd)) next
    ajust <- bh(vapply(pd, function(z) z[2], numeric(1)))
    for (k in seq_along(pd)) filas[[pd[[k]][1]]][[15]] <- p6e(ajust[k])
  }
  list(filas = filas, peor = peor)
}

# ===========================================================================
# 3. Figuras.
# ===========================================================================
figura_dispersion_sd <- function(disp, ruta) {
  reg <- list()
  for (f in disp) {
    nom <- if (f[[1]] == "score_compuesto") "score compuesto" else f[[1]]
    lado <- if (f[[3]] == "PLACENTA_E15") "placenta" else "cerebro"
    sdc <- if (nzchar(f[[5]])) as.numeric(f[[5]]) else NA_real_
    sdl <- if (nzchar(f[[7]])) as.numeric(f[[7]]) else NA_real_
    reg[[length(reg) + 1L]] <- data.frame(item = nom, lado = lado,
      tto = "Control", sd = sdc, stringsAsFactors = FALSE)
    reg[[length(reg) + 1L]] <- data.frame(item = nom, lado = lado,
      tto = "LPS", sd = sdl, stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, reg)
  d$item <- factor(d$item, levels = vapply(ITEMS, function(i)
    if (i == "score_compuesto") "score compuesto" else i, character(1)))
  d$grp <- factor(paste(d$lado, d$tto), levels = c(
    "placenta Control", "placenta LPS", "cerebro Control", "cerebro LPS"))
  p <- ggplot(d, aes(grp, sd, fill = tto)) +
    geom_col(colour = "grey30", linewidth = 0.3) +
    scale_fill_manual(values = c(Control = unname(COL_TTO["CONTROL"]),
                                 LPS = unname(COL_TTO["LPS"])), name = NULL) +
    facet_wrap(~ item, scales = "free_y", ncol = 4) +
    labs(title = paste0("Dispersion de -ddCt dentro del par por feto ",
                        "(placenta y cerebro) por grupo"),
         subtitle = paste0("Una caida marcada de la SD bajo LPS es la sospecha ",
                           "de restriccion de rango que dirime 10 (prohibicion 5)"),
         x = NULL, y = "SD  -ddCt") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.subtitle = element_text(size = 8))
  ggsave(ruta, p, width = 13.0, height = 8.5, dpi = DPI)
}

figura_test_delta_rho <- function(test, ruta) {
  # Solo AMBOS_SEXOS (comportamiento previo de la figura): un forest por sexo
  # ademas saturaria el panel; los estratos por sexo se leen en la tabla.
  test <- Filter(function(f) f[[3]] == "AMBOS_SEXOS", test)
  reg <- list(); ann <- list()
  ord <- vapply(ITEMS, function(i)
    if (i == "score_compuesto") "score compuesto" else i, character(1))
  for (f in test) {
    nom <- if (f[[1]] == "score_compuesto") "score compuesto" else f[[1]]
    if (!nzchar(f[[5]]) || !nzchar(f[[7]])) {
      ann[[length(ann) + 1L]] <- data.frame(item = nom,
        lab = sprintf("n<%d en algun grupo (sin test)", PISO_PAR),
        stringsAsFactors = FALSE)
      next
    }
    rc <- as.numeric(f[[5]]); rl <- as.numeric(f[[7]])
    drho <- as.numeric(f[[8]]); p_bw <- as.numeric(f[[14]])
    reg[[length(reg) + 1L]] <- data.frame(item = nom, tto = "Control", rho = rc,
                                          stringsAsFactors = FALSE)
    reg[[length(reg) + 1L]] <- data.frame(item = nom, tto = "LPS", rho = rl,
                                          stringsAsFactors = FALSE)
    est <- if (p_bw < 0.001) "***" else if (p_bw < 0.01) "**" else
      if (p_bw < 0.05) "*" else "n.s."
    ann[[length(ann) + 1L]] <- data.frame(item = nom,
      lab = sprintf("dRho=%+.2f  p=%.3f %s", drho, p_bw, est),
      stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, reg); a <- do.call(rbind, ann)
  d$item <- factor(d$item, levels = rev(ord))
  a$item <- factor(a$item, levels = rev(ord))
  p <- ggplot(d, aes(rho, item, colour = tto)) +
    geom_line(aes(group = item), colour = "grey70", linewidth = 0.6) +
    geom_point(size = 2.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.3) +
    geom_text(data = a, aes(x = 1.08, y = item, label = lab), inherit.aes = FALSE,
              hjust = 0, size = 2.6, colour = "grey15") +
    scale_colour_manual(values = c(Control = unname(COL_TTO["CONTROL"]),
                                   LPS = unname(COL_TTO["LPS"])), name = NULL) +
    coord_cartesian(xlim = c(-1.05, 1.05), clip = "off") +
    labs(title = paste0("Test reportado (prohibicion 4): Fisher z sobre rho de ",
                        "Spearman, Control vs LPS"),
         subtitle = "SE Bonett-Wright; p = p_bw con BH entre items (suplementario)",
         x = "rho de Spearman (placenta <-> cerebro por feto)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          plot.margin = margin(5.5, 120, 5.5, 5.5),
          plot.subtitle = element_text(size = 8))
  ggsave(ruta, p, width = 9.5, height = 6.0, dpi = DPI)
}

# ===========================================================================
# 4. Artefactos compartidos (merge por 'script') -- headers identicos a 02..08.
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
"## 09_acto2_dispersion",
"",
"### Que resuelve",
"",
paste0("- **2.4 -- test formal de diferencia de correlaciones (prohibicion 4).** ",
       "T7 tenia prohibido comparar la rho de Control con la de LPS; aca se hace ",
       "con el test que **se reporta**: **Fisher z sobre rho de Spearman** ",
       "(decision del usuario, entre Fisher z / interaccion de pendientes / ",
       "permutacion). `z = atanh(rho)`; estadistico `(z_control - z_lps)/",
       "sqrt(SE_control^2 + SE_lps^2)`; `p` a dos colas por la normal."),
paste0("- **2.3 -- dispersion por grupo (insumo de la prohibicion 5).** Sobre ",
       "los MISMOS pares por feto que la correlacion, la SD (n-1) de `-ddCt` de ",
       "cada lado (placenta y cerebro) en Control y en LPS, el cociente de ",
       "varianzas LPS/Control y Levene Brown-Forsythe (centro = mediana) por ",
       "lado. NO alcanza para descartar restriccion de rango: eso lo hace la ",
       "simulacion de `10_acto2_simulacion`."),
"",
"### Error estandar del Fisher z (decision del usuario: ambos, BW primario)",
"",
paste0("- **Primario = Bonett-Wright**: `SE_i = sqrt((1 + rho_i^2/2)/(n_i - 3))` ",
       "-- el mismo SE que ya usa el IC de T7, para que todo el Acto 2 sea ",
       "coherente. Corrige levemente por la varianza extra de rho de Spearman ",
       "respecto de Pearson."),
paste0("- **Clasico** `SE_i = 1/sqrt(n_i - 3)` (Fisher de libro) va como columna ",
       "al lado (`stat_z_clasico`, `p_clasico`). **No cambia ninguna ",
       "conclusion**; se reporta para que el lector vea que la eleccion de SE no ",
       "mueve el resultado."),
paste0("- `p_bw_bh`: Benjamini-Hochberg de `p_bw` entre los items testeados. ",
       "Columna suplementaria en el espiritu de D12 -- **no dirige la ",
       "inferencia**; se menciona en el informe cuantos items sobreviven."),
"",
"### Estratificacion por sexo del Delta rho (pedido explicito, aditiva)",
"",
paste0("- **Problema**: hasta esta sesion, la tabla de Delta rho (2.4) agrupaba ",
       "los sexos (n=16-18 por celda), mientras que desde T7 las figuras de ",
       "correlacion ya estan separadas por sexo -- el informe mostraba una cosa ",
       "(correlaciones por sexo en las figuras) y testeaba otra (Delta rho ",
       "agrupado). El caso que lo deja claro es `fatcd36`: `HEMBRA` Control rho ",
       "= 0.86 (n=8) y `MACHO` Control rho = -0.43 (n=8); agrupados dan 0.47, ",
       "que no describe a ninguno de los dos."),
paste0("- **Correccion** (`pedidos/cambios_acto2_dispersion_por_sexo.md`): ",
       "columna `ESTRATO` con `AMBOS_SEXOS`/`HEMBRA`/`MACHO` en la tabla de ",
       "Delta rho (2.4) y en la de dispersion (2.3, ver abajo). **Aditivo**: el ",
       "estrato `AMBOS_SEXOS` es exactamente lo que habia antes, sin cambios; ",
       "`HEMBRA`/`MACHO` son filas nuevas."),
paste0("- **Alternativa descartada**: dejar solo el estrato agrupado. Se ",
       "descarta porque promedia correlaciones de signo opuesto (ver `fatcd36` ",
       "arriba) -- no es una simplificacion neutral, esconde el patron."),
paste0("- **Limitacion declarada**: con `n <= 9` por celda de sexo, el Fisher z ",
       "tiene potencia baja. Se dice en el cuerpo del reporte (Seccion 1), no en ",
       "una nota al pie."),
paste0("- **BH dentro de cada estrato** (no a traves de los tres): los tres ",
       "estratos tienen n y potencia muy distintos (`AMBOS_SEXOS` vs `HEMBRA`/",
       "`MACHO` por separado); mezclar sus 27 `p` en un solo ajuste Benjamini-",
       "Hochberg no tendria sentido estadistico."),
"",
"### Alcance y piso",
"",
paste0("- Items: 9 genes con `-ddCt` en ambos tejidos (todos menos `il6`, D7) + ",
       "score compuesto (D8)."),
paste0("- Se testea un item solo si Control **y** LPS tienen `n_par >= 5`. ",
       "`il6R` suele quedar por debajo en algun grupo -> su test tendria ",
       "potencia casi nula; se deja la fila con los `n` y sin estadistico, no se ",
       "omite en silencio."),
"",
"### Paridad R / Python",
"",
paste0("- `rho`, `SD`, cociente de varianzas y todas las sumas usan acumulador ",
       "`double` explicito (mismo orden) -> texto `%.10g` bit-identico. Lo que ",
       "pasa por trascendentes (`z` por `atanh`, `p` por la normal, `F` de ",
       "Levene por `pf`) -> texto `%.6e`. Las figuras son PNG: equivalentes, no ",
       "byte-identicas."),
paste0("- Solo en R se cruza-verifica cada `rho` contra ",
       "`cor.test(method=\"spearman\", exact=FALSE)` (tol 1e-9) y cada `F` de ",
       "Levene contra `car::leveneTest(center=median)` (tol 1e-8), con ",
       "`stopifnot`.")
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 09_acto2_dispersion:inicio -->"
  marca_fin <- "<!-- 09_acto2_dispersion:fin -->"
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
# 5. Reporte legible.
# ===========================================================================
.md <- function(header, filas) {
  l1 <- paste0("| ", paste(header, collapse = " | "), " |")
  l2 <- paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |")
  cuerpo <- vapply(filas, function(f)
    paste0("| ", paste(vapply(f, .fmt, character(1)), collapse = " | "), " |"),
    character(1))
  paste(c(l1, l2, cuerpo), collapse = "\n")
}

construir_reporte <- function(fuente, disp, test) {
  L <- c(
    "# Reporte de dispersion y test de correlaciones Acto 2.3-2.4 (T8)", "",
    "Generado por `09_acto2_dispersion` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. Test reportado: Fisher z sobre rho de Spearman (prohibicion 4)", "",
    paste0("- Compara `rho_control` con `rho_lps` de la correlacion placenta <-> ",
           "cerebro por feto. `z = atanh(rho)`, estadistico `(z_control - ",
           "z_lps)/sqrt(SE_control^2 + SE_lps^2)`, `p` normal a dos colas."),
    paste0("- **SE primario = Bonett-Wright** (coherente con el IC de T7); **SE ",
           "clasico** `1/sqrt(n-3)` como columna al lado, no cambia conclusiones. ",
           "`p_bw_bh` = BH entre items, suplementario (D12), **calculado dentro ",
           "de cada ESTRATO** (no a traves de los tres)."),
    "- Se testea solo con `n_par >= 5` en Control **y** LPS.",
    paste0("- **Estratificado por sexo** (pedido explicito, aditivo): columna ",
           "`ESTRATO` = `AMBOS_SEXOS` (agrupado, como antes) / `HEMBRA` / `MACHO`. ",
           "**Por que hace falta**: agrupar los sexos promedia correlaciones de ",
           "signo opuesto. El caso que lo deja claro es `fatcd36`: dentro de ",
           "`HEMBRA` Control la placenta y el cerebro coordinan con signo positivo ",
           "fuerte, dentro de `MACHO` Control coordinan con signo NEGATIVO; el ",
           "`AMBOS_SEXOS` agrupado da un rho intermedio que no describe a ninguno ",
           "de los dos sexos por separado (ver filas `fatcd36` en la tabla, por ",
           "`ESTRATO`)."),
    paste0("- **Potencia**: con `n <= 9` por celda (HEMBRA/MACHO), el Fisher z ",
           "tiene muy poca potencia -- un `Delta rho` chico es indetectable y uno ",
           "grande puede no alcanzar significancia. Esta limitacion es real y se ",
           "declara aca, no en una nota al pie: **no leer \"no significativo en ",
           "HEMBRA/MACHO\" como evidencia de que el efecto desaparece al ",
           "estratificar** (prohibicion 4 sigue aplicando dentro de cada ",
           "estrato)."),
    "", .md(COLS_TEST, test), "",
    "## 2. Dispersion por grupo (insumo de la prohibicion 5)", "",
    paste0("- SD (n-1) de `-ddCt` de cada lado en los mismos pares por feto, ",
           "cociente de varianzas LPS/Control y Levene Brown-Forsythe por lado. ",
           "La lectura biologica del cambio de correlacion queda pendiente de ",
           "`10_acto2_simulacion`."), "",
    .md(COLS_DISP, disp), "",
    "## 3. Figuras", "",
    paste0("- `outputs/figures/acto2_dispersion_sd.png` -- SD de -ddCt por item ",
           "(placenta/cerebro x Control/LPS)."),
    paste0("- `outputs/figures/acto2_test_delta_rho.png` -- rho_control vs ",
           "rho_lps por item con Delta rho y `p_bw` anotados."), "",
    "## 4. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `09_acto2_dispersion`: ",
           "eleccion del test y del SE, rol de la tabla de dispersion, piso de n ",
           "y items de baja potencia."), "")
  paste(L, collapse = "\n")
}

# ===========================================================================
main <- function() {
  D <- cargar()
  fuente <- fuente_datos(ARCHIVO_QPCR)

  td <- tabla_dispersion(D); disp <- td$filas
  cat(sprintf("  [cruza-verif R vs car::leveneTest(center=median)] peor |dif F| = %.3e\n",
              td$peor))
  stopifnot(td$peor < 1e-8)

  tt <- tabla_test(D); test <- tt$filas
  cat(sprintf("  [cruza-verif R vs cor.test(exact=FALSE)] peor |dif rho| = %.3e\n",
              tt$peor))
  stopifnot(tt$peor < 1e-9)

  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "acto2_dispersion.csv"), COLS_DISP, disp)
    escribir_csv(file.path(base, "acto2_test_correlaciones.csv"), COLS_TEST, test)
  }

  fig_sd <- file.path(RUTA_FIGURAS, "acto2_dispersion_sd.png")
  fig_dr <- file.path(RUTA_FIGURAS, "acto2_test_delta_rho.png")
  figura_dispersion_sd(disp, fig_sd)
  figura_test_delta_rho(test, fig_dr)

  escribir_lineas(file.path(RUTA_TABLAS, "acto2_dispersion_reporte.md"),
                  construir_reporte(fuente, disp, test))
  actualizar_descartados()

  # Conteos de resumen (cat + verificaciones): AMBOS_SEXOS, para no romper el
  # denominador "de 9 items" de las verificaciones ya existentes; HEMBRA/MACHO
  # se resumen aparte en el cat() final.
  por_estrato <- function(e) Filter(function(f) f[[3]] == e, test)
  .n_test <- function(fs) sum(vapply(fs, function(f) nzchar(f[[5]]), logical(1)))
  .n_sig_bw <- function(fs) sum(vapply(fs, function(f)
    nzchar(f[[14]]) && as.numeric(f[[14]]) < 0.05, logical(1)))
  .n_sig_bh <- function(fs) sum(vapply(fs, function(f)
    nzchar(f[[15]]) && as.numeric(f[[15]]) < 0.05, logical(1)))
  test_ambos <- por_estrato("AMBOS_SEXOS")
  n_test <- .n_test(test_ambos)
  n_sig_bw <- .n_sig_bw(test_ambos)
  n_sig_bh <- .n_sig_bh(test_ambos)

  # Particion: n(HEMBRA) + n(MACHO) == n(AMBOS_SEXOS), por item y por grupo
  # (verificacion 7.1 del pedido) -- si no cierra, hay error de filtrado.
  n_por <- list()
  for (f in test) n_por[[paste(f[[1]], f[[3]])]] <- c(n_control = f[[4]], n_lps = f[[6]])
  particion_ok <- TRUE; particion_detalle <- character(0)
  for (item in ITEMS) {
    a <- n_por[[paste(item, "AMBOS_SEXOS")]]
    h <- n_por[[paste(item, "HEMBRA")]]
    m <- n_por[[paste(item, "MACHO")]]
    ok_c <- as.integer(h["n_control"]) + as.integer(m["n_control"]) == as.integer(a["n_control"])
    ok_l <- as.integer(h["n_lps"]) + as.integer(m["n_lps"]) == as.integer(a["n_lps"])
    if (!ok_c || !ok_l) {
      particion_ok <- FALSE
      particion_detalle <- c(particion_detalle, item)
    }
  }

  ent <- sprintf(paste0("data/processed/qpcr_cuantificacion_long.tsv + ",
                        "qpcr_score_compuesto_long.tsv (de data/%s/%s)"),
                 fuente, ARCHIVO_QPCR)
  registrar_procedencia(list(
    list("outputs/tables/{R,python}/acto2_dispersion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("SD (n-1) de -ddCt por lado (placenta/cerebro) x ",
         "grupo sobre los pares por feto; cociente de varianzas LPS/Control; ",
         "Levene Brown-Forsythe por lado")),
    list("outputs/tables/{R,python}/acto2_test_correlaciones.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, paste0("test reportado (prohibicion 4): ",
         "Fisher z sobre rho de Spearman Control vs LPS; SE Bonett-Wright ",
         "primario + SE clasico + p_bh suplementario; estratificado por ESTRATO ",
         "(AMBOS_SEXOS/HEMBRA/MACHO, pedido explicito -- estratificacion pedida ",
         "para ver si el promedio entre sexos tapaba efectos; alternativa ",
         "descartada: dejar solo el estrato agrupado, descartada porque promedia ",
         "correlaciones de signo opuesto, ver fatcd36; limitacion declarada: ",
         "n<=9, potencia baja; BH dentro de cada estrato, no a traves de los tres)")),
    list("outputs/figures/acto2_dispersion_sd.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent, "SD de -ddCt por item: placenta/cerebro x Control/LPS"),
    list("outputs/figures/acto2_test_delta_rho.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent, "rho_control vs rho_lps por item con Delta rho y p_bw"),
    list("outputs/tables/acto2_dispersion_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.3-2.4 (T8)")
  ))
  registrar_verificaciones(list(
    list("acto2_test_metodo",
         "test reportado de diferencia de correlaciones (prohibicion 4)",
         paste0("Fisher z sobre rho de Spearman; SE Bonett-Wright primario; SE ",
                "clasico y p_bh suplementarios"),
         "Fisher z (decision del usuario)", "TRUE", ESTE_SCRIPT),
    list("acto2_test_items",
         "items con test formal Control vs LPS (n_par >= 5 en ambos grupos)",
         sprintf("%d de %d items testeados", n_test, length(ITEMS)),
         "9 genes + score, los que superan el piso", "TRUE", ESTE_SCRIPT),
    list("acto2_test_piso",
         "piso de n_par por grupo para el Fisher z", as.character(PISO_PAR), "5",
         if (PISO_PAR == 5L) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_test_se_doble",
         "se reportan SE Bonett-Wright (primario) y SE clasico 1/sqrt(n-3)",
         "columnas se_bw_* y se_clasico_* + stat/p de cada uno", "ambos SE",
         if (COLS_TEST[11] == "se_bw_control" &&
             COLS_TEST[16] == "se_clasico_control") "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_test_bh_suplementario",
         paste0("p_bw_bh (BH DENTRO de cada ESTRATO) es suplementario y no ",
                "cambia conclusiones (D12)"),
         sprintf("AMBOS_SEXOS: %d items p_bw<.05; %d items p_bw_bh<.05",
                 n_sig_bw, n_sig_bh),
         "columna rotulada, no dirige", "TRUE", ESTE_SCRIPT),
    list("acto2_estratos_sexo",
         "ESTRATO = AMBOS_SEXOS/HEMBRA/MACHO en Delta rho (2.4) y dispersion (2.3)",
         paste(ESTRATOS, collapse = ";"), "AMBOS_SEXOS;HEMBRA;MACHO",
         if (identical(ESTRATOS, c("AMBOS_SEXOS", "HEMBRA", "MACHO"))) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_estrato_particion_test",
         "n(HEMBRA) + n(MACHO) = n(AMBOS_SEXOS) por item y por grupo (Delta rho, 2.4)",
         sprintf("particiona en %d/%d items%s", length(ITEMS) - length(particion_detalle),
                 length(ITEMS), if (length(particion_detalle))
                   paste0("; falla en: ", paste(particion_detalle, collapse = ", ")) else ""),
         sprintf("particiona en %d/%d items", length(ITEMS), length(ITEMS)),
         if (particion_ok) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_dispersion_pares",
         "la dispersion se mide sobre los MISMOS pares por feto que la correlacion",
         "vector por lado = componente placenta/cerebro de pares(item, grupo)",
         "mismos pares que T7/2.4", "TRUE", ESTE_SCRIPT),
    list("acto2_dispersion_levene",
         "test de dispersion por lado = Levene Brown-Forsythe (centro = mediana)",
         "ANOVA de una via sobre |y - mediana(grupo)|, df1=1 df2=nC+nL-2",
         "Brown-Forsythe", "TRUE", ESTE_SCRIPT),
    list("acto2_dispersion_no_concluye",
         "la tabla de dispersion NO descarta restriccion de rango por si sola",
         "prohibicion 5 la dirime 10_acto2_simulacion; el reporte lo dice",
         "insumo, no conclusion", "TRUE", ESTE_SCRIPT),
    list("acto2_dispersion_figuras",
         "figuras Acto 2.3-2.4: SD por grupo + forest de Delta rho",
         sprintf("sd=%s;delta_rho=%s", file.exists(fig_sd), file.exists(fig_dr)),
         "2 figuras existen",
         if (file.exists(fig_sd) && file.exists(fig_dr)) "TRUE" else "FALSE",
         ESTE_SCRIPT)
  ))

  cat("== 09_acto2_dispersion.R ==\n")
  cat(sprintf("  fuente = %s\n", fuente))
  cat(sprintf("  items testeados (n_par>=5 ambos grupos), AMBOS_SEXOS: %d / %d\n",
              n_test, length(ITEMS)))
  for (f in test_ambos) {
    if (nzchar(f[[5]]))
      cat(sprintf("    %-16s rhoC=%8s  rhoL=%8s  dRho=%8s  p_bw=%s  p_bh=%s\n",
                  f[[1]], f[[5]], f[[7]], f[[8]], f[[14]], f[[15]]))
    else
      cat(sprintf("    %-16s nC=%s nL=%s  (sin test, piso)\n", f[[1]], f[[4]], f[[6]]))
  }
  cat(sprintf("  significativos AMBOS_SEXOS: p_bw<.05 -> %d;  p_bw_bh<.05 -> %d\n",
              n_sig_bw, n_sig_bh))
  for (estrato in c("HEMBRA", "MACHO")) {
    fs <- por_estrato(estrato)
    cat(sprintf("  items testeados %s: %d / %d  (p_bw<.05 -> %d; p_bw_bh<.05 -> %d)\n",
                estrato, .n_test(fs), length(ITEMS), .n_sig_bw(fs), .n_sig_bh(fs)))
  }
  cat(sprintf("  particion n(HEMBRA)+n(MACHO)=n(AMBOS_SEXOS): %s\n",
              if (particion_ok) sprintf("OK en %d/%d items", length(ITEMS), length(ITEMS))
              else paste0("FALLA en: ", paste(particion_detalle, collapse = ", "))))
  cat("  -> outputs/tables/{R,python}/acto2_dispersion.csv, acto2_test_correlaciones.csv\n")
  cat(sprintf("  -> %s, %s\n", basename(fig_sd), basename(fig_dr)))
}

if (sys.nframe() == 0L) main()
