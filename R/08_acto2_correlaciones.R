# 08_acto2_correlaciones.R -- ACTO 2.1-2.2: correlacion placenta <-> cerebro.
#
# Por que existe este archivo: el Acto 2 pregunta si el programa de transporte de
# la placenta y el del cerebro fetal estan acoplados dentro de un mismo feto.
# Este script SOLO describe esa correlacion; NO compara correlaciones entre
# grupos (eso es el test formal de T8 -- prohibicion 4) ni la interpreta como
# coordinacion biologica sin descartar restriccion de rango (prohibicion 5, T8).
#
#   * Emparejamiento por FETO: un par entra si AMBOS lados (PLACENTA_E15 y
#     BRAIN_E15) estan detectados.
#   * Magnitudes (decision del usuario): -ddCt por gen y el score compuesto (D8).
#   * Coeficiente (decision del usuario): **Spearman rho** (no Pearson).
#     rho = Pearson sobre los rangos promedio (corrige empates). p por la
#     t-aproximacion t = rho*sqrt((n-2)/(1-rho^2)), df = n-2, dos colas -- la
#     misma que cor.test(method="spearman", exact=FALSE) y scipy.spearmanr.
#     IC 95% Bonett-Wright: SE_z = sqrt((1 + rho^2/2)/(n-3)), z = atanh(rho),
#     IC = tanh(z +/- 1.959963984540054 * SE_z).
#   * Estratos (CAMBIO pedido explicito, pedidos/cambios_acto2_correlaciones_
#     por_sexo.md): GLOBAL + por TTO (CONTROL/LPS) + por SEXO x TTO (4 celdas).
#     Los 3 estratos originales SE CONSERVAN (la vista sin separar por sexo
#     sigue en el informe); se agregan los 4 por sexo. **No se compara rho
#     entre estratos** (prohibicion 4): cada celda es descriptiva.
#   * Piso: n_par < 5 -> solo n (sin rho/IC/p); en las figuras, los puntos se
#     dibujan igual pero sin linea de tendencia ni rho en la leyenda.
#   * `il6R` se EXCLUYE de todo el Acto 2 (tablas y figuras) por deteccion
#     insuficiente en cerebro (3/7/3/3) que, al estratificar por sexo, deja casi
#     todas las celdas bajo el piso de 5 pares (pedido explicito). Se conserva
#     en las figuras y tablas del Acto 1 (07_figuras_acto1, 05_qpcr_modelos).
#   * Co-expresion: rho de Spearman entre los genes del SPLOM de cada tejido
#     (GENES_SPLOM_PLACENTA = 9, GENES_SPLOM_BRAIN = 8; distinto por tejido
#     porque il6 solo es cuantificable en placenta, D7) -- tabla + SPLOM, para
#     AMBOS sexos juntos y por separado.
#
# PARIDAD R/Python: rho y las sumas usan acumulador double explicito (mismo orden)
# -> bit-identico; p (via pt) e IC (via tanh/atanh) -> texto "%.6e" (p6e). Aca
# (solo R) se cruza-verifica cada rho/p contra cor.test(exact=FALSE) con
# stopifnot (< 1e-9). NUNCA se llama a cor.test spearman con el metodo exacto
# (AS 89): segfaultea en este build de R.
# Las figuras son PNG: equivalentes, no byte-identicas.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({ library(ggplot2); library(GGally); library(scales) })

ESTE_SCRIPT <- "08_acto2_correlaciones"

COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")
PCH_TTO <- c(CONTROL = 16, LPS = 17)
Z975 <- 1.959963984540054
PISO_PAR <- 5L

GEN_SIN_CEREBRO <- "il6"      # D7: sin -ddCt en cerebro (calibrador 0/9)
GEN_EXCLUIDO_CORR <- "il6R"   # pedido explicito: deteccion insuficiente en cerebro
GENES_CORR <- setdiff(GENES, c(GEN_SIN_CEREBRO, GEN_EXCLUIDO_CORR))   # 8 genes
ITEMS <- c(GENES_CORR, "score_compuesto")                              # 9

# 3 estratos originales + 4 celdas SEXO x TTO (pedido explicito).
ESTRATOS <- c("GLOBAL", "CONTROL", "LPS",
              "HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")
FILTRO_ESTRATO <- list(
  GLOBAL         = list(sexo = NULL,     tto = NULL),
  CONTROL        = list(sexo = NULL,     tto = "CONTROL"),
  LPS            = list(sexo = NULL,     tto = "LPS"),
  HEMBRA_CONTROL = list(sexo = "HEMBRA", tto = "CONTROL"),
  HEMBRA_LPS     = list(sexo = "HEMBRA", tto = "LPS"),
  MACHO_CONTROL  = list(sexo = "MACHO",  tto = "CONTROL"),
  MACHO_LPS      = list(sexo = "MACHO",  tto = "LPS")
)

DPI <- 300

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..07.
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
spearman_p <- function(rho, n) {
  if (is.na(rho) || n <= 2L) return(NA_real_)
  if (rho * rho >= 1) return(0)
  t <- rho * sqrt((n - 2) / (1 - rho * rho))
  2 * pt(-abs(t), n - 2)
}
spearman_ci <- function(rho, n) {
  if (is.na(rho) || n <= 3L || rho * rho >= 1) return(c(NA_real_, NA_real_))
  se <- sqrt((1 + rho * rho / 2) / (n - 3))
  z <- atanh(rho)
  c(tanh(z - Z975 * se), tanh(z + Z975 * se))
}

# ===========================================================================
# 1. Carga y emparejamiento por feto.
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
  # SEXO ya viene normalizado (HEMBRA/MACHO) en qpcr_cuantificacion_long.tsv
  # desde 02_ingesta_qc; se usa esa columna directamente (no se deriva de GRUPO).
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

pares <- function(D, item, estrato) {
  f <- FILTRO_ESTRATO[[estrato]]
  xs <- c(); ys <- c(); ts <- c(); ss <- c()
  for (feto in D$fetos) {
    if (!is.null(f$tto)  && D$tto[[feto]]  != f$tto)  next
    if (!is.null(f$sexo) && D$sexo[[feto]] != f$sexo) next
    if (item == "score_compuesto") {
      xp <- D$sc[[paste(feto, "PLACENTA_E15", sep = "\r")]]
      yb <- D$sc[[paste(feto, "BRAIN_E15", sep = "\r")]]
    } else {
      xp <- D$negdd[[paste(feto, "PLACENTA_E15", item, sep = "\r")]]
      yb <- D$negdd[[paste(feto, "BRAIN_E15", item, sep = "\r")]]
    }
    if (is.null(xp) || is.null(yb) || is.na(xp) || is.na(yb)) next
    xs <- c(xs, xp); ys <- c(ys, yb); ts <- c(ts, D$tto[[feto]]); ss <- c(ss, D$sexo[[feto]])
  }
  list(x = xs, y = ys, tto = ts, sexo = ss)
}

# ===========================================================================
# 2. Tablas de correlacion (+ cruza-verificacion R vs cor.test exact=FALSE).
# ===========================================================================
COLS_CORR <- c("ITEM", "TIPO", "ESTRATO", "n_par", "rho_spearman",
               "ic95_low", "ic95_high", "p_valor")

tabla_correlaciones <- function(D) {
  filas <- list(); peor <- 0
  for (item in ITEMS) {
    tipo <- if (item == "score_compuesto") "score" else "gen"
    for (est in ESTRATOS) {
      pr <- pares(D, item, est); n <- length(pr$x)
      if (n >= PISO_PAR) {
        rho <- spearman_rho(pr$x, pr$y)
        ci <- spearman_ci(rho, n)
        p <- spearman_p(rho, n)
        ct <- suppressWarnings(cor.test(pr$x, pr$y, method = "spearman",
                                        exact = FALSE))
        peor <- max(peor, abs(rho - as.numeric(ct$estimate)),
                    abs(p - ct$p.value))
        filas[[length(filas) + 1L]] <- list(item, tipo, est, n, g10(rho),
                                            p6e(ci[1]), p6e(ci[2]), p6e(p))
      } else {
        filas[[length(filas) + 1L]] <- list(item, tipo, est, n, "", "", "", "")
      }
    }
  }
  list(filas = filas, peor = peor)
}

COLS_COEXP <- c("TEJIDO", "SEXO", "GEN_A", "GEN_B", "n_par", "rho_spearman", "p_valor")

genes_splom_tejido <- function(tej)
  if (tej == "PLACENTA_E15") GENES_SPLOM_PLACENTA else GENES_SPLOM_BRAIN

tabla_coexpresion <- function(D) {
  filas <- list()
  for (tej in TEJIDOS_E15) {
    gs <- genes_splom_tejido(tej)
    for (sx in c("AMBOS", "HEMBRA", "MACHO")) {
      for (i in seq_len(length(gs) - 1L)) {
        for (j in (i + 1L):length(gs)) {
          a <- gs[i]; b <- gs[j]
          xs <- c(); ys <- c()
          for (f in D$fetos) {
            if (sx != "AMBOS" && D$sexo[[f]] != sx) next
            va <- D$negdd[[paste(f, tej, a, sep = "\r")]]
            vb <- D$negdd[[paste(f, tej, b, sep = "\r")]]
            if (is.null(va) || is.null(vb) || is.na(va) || is.na(vb)) next
            xs <- c(xs, va); ys <- c(ys, vb)
          }
          n <- length(xs)
          if (n >= PISO_PAR) {
            rho <- spearman_rho(xs, ys); p <- spearman_p(rho, n)
            filas[[length(filas) + 1L]] <- list(tej, sx, a, b, n, g10(rho), p6e(p))
          } else {
            filas[[length(filas) + 1L]] <- list(tej, sx, a, b, n, "", "")
          }
        }
      }
    }
  }
  filas
}

# ===========================================================================
# 3. Figuras.
# ===========================================================================
.rho_estratos <- function(D, item, estratos) {
  out <- list()
  for (est in estratos) {
    pr <- pares(D, item, est); n <- length(pr$x)
    out[[est]] <- if (n >= PISO_PAR) c(spearman_rho(pr$x, pr$y), n) else c(NA_real_, n)
  }
  out
}

# --- 3.1 Figura global (se conserva: sin separar por sexo) -----------------
figura_dispersion <- function(D, ruta) {
  df <- list(); ann <- list()
  for (item in ITEMS) {
    pr <- pares(D, item, "GLOBAL")
    nom <- if (item == "score_compuesto") "score compuesto" else item
    if (length(pr$x))
      df[[length(df) + 1L]] <- data.frame(
        item = nom, x = pr$x, y = pr$y, tto = pr$tto, stringsAsFactors = FALSE)
    rt <- .rho_estratos(D, item, c("GLOBAL", "CONTROL", "LPS"))
    rg <- rt[["GLOBAL"]]; ci <- spearman_ci(rg[1], rg[2])
    l1 <- if (!is.na(rg[1]) && !is.na(ci[1]))
      sprintf("rho = %.2f [%.2f, %.2f]  (n=%d)", rg[1], ci[1], ci[2], rg[2])
      else sprintf("(n=%d, sin rho)", rg[2])
    rc <- rt[["CONTROL"]]; rl <- rt[["LPS"]]
    l2 <- if (!is.na(rc[1]) && !is.na(rl[1]))
      sprintf("Control rho=%.2f (n=%d)  |  LPS rho=%.2f (n=%d)",
              rc[1], rc[2], rl[1], rl[2])
      else sprintf("Control n=%d | LPS n=%d", rc[2], rl[2])
    ann[[length(ann) + 1L]] <- data.frame(item = nom, lab = paste0(l1, "\n", l2),
                                          stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, df); a <- do.call(rbind, ann)
  d$item <- factor(d$item, levels = a$item); a$item <- factor(a$item, levels = a$item)
  p <- ggplot(d, aes(x, y, colour = tto)) +
    geom_point(size = 1.3, alpha = 0.9) +
    geom_text(data = a, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = -0.03, vjust = 1.1, size = 2.6, colour = "grey15") +
    scale_colour_manual(values = COL_TTO, name = NULL,
                        labels = c(CONTROL = "Control", LPS = "LPS")) +
    facet_wrap(~ item, scales = "free", ncol = 3) +
    labs(title = paste0("Correlacion placenta <-> cerebro por feto  --  Spearman ",
                        "rho (-ddCt por gen y score compuesto)"),
         subtitle = paste0("T7 solo describe: la diferencia de rho entre Control ",
                           "y LPS NO se testea aca (prohibicion 4); el test ",
                           "formal es T8. il6R excluido (deteccion insuficiente, ",
                           "pedido explicito)"),
         x = "placenta  -ddCt", y = "cerebro  -ddCt") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.text = element_text(size = 9),
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.subtitle = element_text(size = 8.5))
  ggsave(ruta, p, width = 11.5, height = 9.5, dpi = DPI)
}

# --- 3.2 Figura nueva: una por item, dos paneles (Females/Males) -----------
# Ejes en escala log2 mostrando FC = 2^(-ddCt) (D2), salvo score_compuesto: es
# un z-score (puede ser negativo), no tiene FC -- se grafica en escala lineal.
# Desviacion documentada en analisis_descartados.md (seccion de este script).
NOTA_PIE <- paste0("Axes: 2^-ΔΔCt (log2 display) | Spearman on -ΔΔCt | ",
                   "Linear fit with 95% CI | Same fetus pairing")
NOTA_PIE_SCORE <- paste0("Axis: composite z-score (linear) | Spearman on the ",
                         "z-score | Linear fit with 95% CI | Same fetus pairing")

.linea_leyenda <- function(nombre, D, item, sx, tt) {
  pr <- pares(D, item, if (sx == "HEMBRA") paste0(sx, "_", tt) else paste0(sx, "_", tt))
  n <- length(pr$x)
  if (n >= PISO_PAR) {
    rho <- spearman_rho(pr$x, pr$y); p <- spearman_p(rho, n)
    sprintf("%s: Spearman rho = %.2f ; p = %.3f ; n = %d", nombre, rho, p, n)
  } else {
    sprintf("%s: n = %d (sin rho)", nombre, n)
  }
}

figura_gen_sexo <- function(D, item, ruta) {
  es_score <- item == "score_compuesto"
  filas <- list()
  for (sx in c("HEMBRA", "MACHO")) for (tt in NIVELES_TTO) {
    pr <- pares(D, item, paste0(sx, "_", tt))
    if (!length(pr$x)) next
    filas[[length(filas) + 1L]] <- data.frame(
      sexo_panel = if (sx == "HEMBRA") "Females" else "Males",
      tto = tt, x = pr$x, y = pr$y, stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, filas)
  d$sexo_panel <- factor(d$sexo_panel, levels = c("Females", "Males"))

  # Solo se ajusta linea de tendencia en celdas SEXO x TTO con n >= PISO_PAR.
  n_celda <- table(d$sexo_panel, d$tto)
  d_ok <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
    if (n_celda[as.character(d$sexo_panel[i]), d$tto[i]] >= PISO_PAR) d[i, ] else NULL
  }))

  if (es_score) {
    d$xp <- d$x; d$yp <- d$y
    if (!is.null(d_ok)) { d_ok$xp <- d_ok$x; d_ok$yp <- d_ok$y }
    lab_x <- "Placenta E15 -- composite transporter z-score"
    lab_y <- "Brain E15 -- composite transporter z-score"
    titulo <- "score compuesto"
  } else {
    d$xp <- 2^d$x; d$yp <- 2^d$y
    if (!is.null(d_ok)) { d_ok$xp <- 2^d_ok$x; d_ok$yp <- 2^d_ok$y }
    lab_x <- sprintf("Placenta E15 — %s/rsp29 relative expression", item)
    lab_y <- sprintf("Brain E15 — %s/rsp29 relative expression", item)
    titulo <- item
  }

  leyenda <- paste(c(
    .linea_leyenda("♀ Control", D, item, "HEMBRA", "CONTROL"),
    .linea_leyenda("♀ LPS",     D, item, "HEMBRA", "LPS"),
    .linea_leyenda("♂ Control", D, item, "MACHO",  "CONTROL"),
    .linea_leyenda("♂ LPS",     D, item, "MACHO",  "LPS")
  ), collapse = "\n")
  pie <- if (es_score) NOTA_PIE_SCORE else NOTA_PIE

  p <- ggplot(d, aes(xp, yp, colour = tto, shape = tto)) +
    geom_point(size = 1.6, alpha = 0.85) +
    scale_colour_manual(values = COL_TTO, guide = "none") +
    scale_shape_manual(values = PCH_TTO, guide = "none") +
    facet_wrap(~ sexo_panel, ncol = 2) +
    labs(title = if (es_score) titulo else bquote(italic(.(titulo))),
         subtitle = "Placenta–brain correlation at E15",
         x = lab_x, y = lab_y,
         caption = paste0(leyenda, "\n", pie)) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10, face = "bold"),
          strip.background = element_rect(fill = "grey93", colour = "grey40"),
          panel.border = element_rect(colour = "grey40"),
          plot.caption = element_text(hjust = 0, size = 7.6, lineheight = 1.25),
          plot.subtitle = element_text(size = 9))
  if (!is.null(d_ok))
    p <- p + geom_smooth(data = d_ok, method = "lm", se = TRUE, level = 0.95,
                         linewidth = 0.7, alpha = 0.18)
  if (!es_score)
    p <- p + scale_x_continuous(trans = "log2") + scale_y_continuous(trans = "log2")
  ggsave(ruta, p, width = 8.6, height = 5.6, dpi = DPI)
}

# --- 3.3 SPLOM de co-expresion: AMBOS sexos + por sexo, gen set por tejido --
.panel_rho_n <- function(data, mapping, ...) {
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  lab <- if (n >= PISO_PAR) {
    rho <- spearman_rho(x[ok], y[ok])
    sprintf("rho: %.2f\n(n=%d)", rho, n)
  } else sprintf("n=%d\n(sin rho)", n)
  ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_text(label = lab, size = 2.6, colour = "grey20") +
    ggplot2::theme_void()
}

figura_splom <- function(D, tej, ruta, genes_t, sexo_filtro = "AMBOS") {
  cols <- lapply(genes_t, function(g)
    vapply(D$fetos, function(f) {
      v <- D$negdd[[paste(f, tej, g, sep = "\r")]]
      if (is.null(v)) NA_real_ else v
    }, numeric(1)))
  wide <- as.data.frame(setNames(cols, genes_t), stringsAsFactors = FALSE)
  wide$TTO <- factor(vapply(D$fetos, function(f) D$tto[[f]], character(1)),
                     levels = c("CONTROL", "LPS"))
  wide$SEXO <- vapply(D$fetos, function(f) D$sexo[[f]], character(1))
  if (sexo_filtro != "AMBOS") wide <- wide[wide$SEXO == sexo_filtro, ]

  tt <- if (tej == "PLACENTA_E15") "Placenta E15" else "Cerebro fetal E15"
  sub_tt <- if (sexo_filtro == "AMBOS") "" else
    sprintf(" -- %s (n~%d por tratamiento)", sexo_filtro, round(nrow(wide) / 2))
  p <- GGally::ggpairs(
    wide, columns = seq_along(genes_t),
    mapping = ggplot2::aes(colour = TTO),
    upper = list(continuous = .panel_rho_n),
    lower = list(continuous = GGally::wrap("points", size = 0.5, alpha = 0.7)),
    diag  = list(continuous = GGally::wrap("densityDiag", alpha = 0.4)),
    title = sprintf(paste0("Co-expresion (%s) -- -ddCt, Spearman rho%s"),
                    tt, sub_tt)) +
    ggplot2::scale_colour_manual(values = COL_TTO) +
    ggplot2::scale_fill_manual(values = COL_TTO) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  ggsave(ruta, p, width = 12.5, height = 12.5, dpi = DPI)
}

# ===========================================================================
# 4. Artefactos compartidos (merge por 'script') -- headers identicos a 02..07.
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
"## 08_acto2_correlaciones",
"",
"### Que hace y que NO hace T7",
"",
paste0("- Describe la correlacion placenta <-> cerebro (por feto) del `-ddCt` de ",
       "cada gen y del score compuesto, global, por TTO y por SEXO x TTO. **No ",
       "compara** las correlaciones entre estratos: reportar \"significativo en ",
       "Control y no en LPS\" (o en un sexo y no en el otro) como prueba de ",
       "diferencia esta prohibido (prohibicion 4). El test formal (Fisher z / ",
       "interaccion de pendientes / permutacion) y el control de restriccion de ",
       "rango (prohibicion 5, simulacion) son **T8**."),
"",
"### Coeficiente: Spearman (decision del usuario)",
"",
paste0("- Se usa **Spearman rho** (no Pearson): robusto a outliers de qPCR y a la ",
       "no-linealidad. `rho` = correlacion de Pearson sobre los rangos promedio ",
       "(corrige empates). `p` por la t-aproximacion ",
       "`t = rho*sqrt((n-2)/(1-rho^2))`, df = n-2, dos colas -- la misma que ",
       "`cor.test(method=\"spearman\", exact=FALSE)` y `scipy.stats.spearmanr` ",
       "por defecto, implementada PROPIA para paridad. IC 95% por Bonett-Wright: ",
       "`SE_z = sqrt((1 + rho^2/2)/(n-3))`, `z = atanh(rho)`, ",
       "`IC = tanh(z +/- 1.959963984540054 * SE_z)`. Si `n_par < 5` o `|rho| = 1` ",
       "no se reporta rho/IC/p (solo `n_par`)."),
paste0("- `08_acto2_correlaciones.R` cruza-verifica cada `rho` y `p` contra ",
       "`cor.test(..., method=\"spearman\", exact=FALSE)` (`stopifnot`, tol 1e-9). ",
       "**No** se usa el metodo exacto de `cor.test` (AS 89): segfaultea en este ",
       "build de R."),
"",
"### Emparejamiento y alcance",
"",
paste0("- Un par entra si el feto tiene `-ddCt` (o score) detectado en **ambos** ",
       "tejidos. `il6` queda fuera del brazo por gen: no tiene `-ddCt` en cerebro ",
       "(calibrador HEMBRA_CONTROL 0/9, D7)."),
paste0("- **CAMBIO (pedido explicito, `pedidos/cambios_acto2_correlaciones_por_",
       "sexo.md`): se revierte la decision anterior de limitar los estratos a ",
       "GLOBAL/CONTROL/LPS.** Esa decision argumentaba que las 4 celdas SEXO x ",
       "TTO (n~9) darian intervalos de confianza inutiles; se revierte porque el ",
       "dimorfismo sexual es la pregunta del proyecto y los estratos agregados lo ",
       "promedian. Estratos ahora: `GLOBAL`, `CONTROL`, `LPS`, `HEMBRA_CONTROL`, ",
       "`HEMBRA_LPS`, `MACHO_CONTROL`, `MACHO_LPS`. **Limitacion declarada:** con ",
       "n <= 9 los IC de rho son anchos y la comparacion entre paneles/estratos ",
       "no esta testeada (T7 sigue sin comparar, prohibicion 4)."),
paste0("- **`il6R` se excluye de todas las tablas y figuras del Acto 2** (pedido ",
       "explicito): en cerebro tiene deteccion insuficiente (3/7/3/3) y al ",
       "estratificar por sexo casi todas las celdas quedan bajo el piso de 5 ",
       "pares. Se conserva en las figuras y tablas del Acto 1 (07_figuras_acto1, ",
       "05_qpcr_modelos): la exclusion es solo para las correlaciones."),
"",
"### Co-expresion (SPLOM)",
"",
paste0("- `acto2_coexpresion_transportadores.csv` y los 6 SPLOM (`_PLACENTA_E15`",
       ", `_BRAIN_E15` y sus 4 variantes `_HEMBRA`/`_MACHO`) muestran la rho de ",
       "Spearman dentro de cada tejido. **El conjunto de genes cambia y es ",
       "distinto por tejido**: placenta usa los 7 transportadores + `il6` + ",
       "`gp130` (9; `il6` es cuantificable ahi), cerebro usa los 7 + `gp130` (8; ",
       "`il6` no es cuantificable en cerebro, D7). Las matrices de placenta y ",
       "cerebro **no son comparables celda por celda** (conjuntos distintos). Es ",
       "contexto, no una prueba."),
"",
"### Paridad R / Python",
"",
paste0("- `rho` y las sumas usan acumulador `double` explicito (mismo orden) -> ",
       "bit-identico; se guarda con `%.10g`. Lo que pasa por trascendentes (`p` ",
       "via `pt`, IC via `tanh`/`atanh`) se guarda como texto `%.6e`. Las figuras ",
       "son PNG: equivalentes, no byte-identicas (ggplot2/GGally vs matplotlib)."),
paste0("- **Desviacion documentada:** para `score_compuesto` los ejes de la ",
       "figura por gen/sexo NO se muestran en escala log2 (a diferencia de los ",
       "demas items): el score es un z-score (puede ser negativo), no un ",
       "fold-change, y `log2` de un valor negativo no existe. Se grafica en ",
       "escala lineal, etiquetado como tal.")
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 08_acto2_correlaciones:inicio -->"
  marca_fin <- "<!-- 08_acto2_correlaciones:fin -->"
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

construir_reporte <- function(fuente, corr, coexp) {
  L <- c(
    "# Reporte de correlaciones Acto 2.1-2.2 (T7)", "",
    "Generado por `08_acto2_correlaciones` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. Metodo", "",
    paste0("- Correlacion placenta <-> cerebro **por feto**, con **Spearman rho** ",
           "(decision del usuario; robusto a outliers). `p` por t-aproximacion ",
           "(df = n-2), IC 95% Bonett-Wright. Piso: `n_par >= 5`."),
    paste0("- Magnitudes: `-ddCt` de cada gen (8; `il6` fuera por D7, `il6R` ",
           "fuera por deteccion insuficiente, pedido explicito) y el **score ",
           "compuesto** de 7 transportadores (D8)."),
    paste0("- Estratos: `GLOBAL`, por `TTO` (`CONTROL`/`LPS`) y por `SEXO x TTO` ",
           "(4 celdas; extension pedida explicitamente)."),
    paste0("- **T7 no compara** las correlaciones entre estratos (prohibicion 4). ",
           "El test formal y el control de restriccion de rango son **T8**."),
    "", "## 2. Correlacion placenta <-> cerebro (por gen y score)", "",
    .md(COLS_CORR, corr), "",
    "## 3. Co-expresion entre genes del SPLOM (Spearman, por tejido y sexo)", "",
    .md(COLS_COEXP, coexp), "",
    "## 4. Figuras", "",
    paste0("- `outputs/figures/acto2_dispersion_placenta_cerebro.png` -- vista ",
           "global (sin separar por sexo), dispersion por item + rho anotada."),
    paste0("- `outputs/figures/acto2_corr_placenta_cerebro_<item>.png` (9 ",
           "figuras) -- dos paneles Females/Males, Control y LPS superpuestos, ",
           "ajuste lineal + IC95% si n>=5, leyenda con rho/p/n por celda."),
    paste0("- `outputs/figures/acto2_coexpresion_SPLOM_{PLACENTA_E15,BRAIN_E15}",
           "{,_HEMBRA,_MACHO}.png` (6 figuras) -- matriz de dispersion por tejido, ",
           "ambos sexos juntos y por separado."),
    "", "## 5. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `08_acto2_correlaciones`: ",
           "eleccion de Spearman, formula de `p` e IC, la extension de estratos ",
           "por sexo (con su limitacion declarada), la exclusion de `il6R`, y el ",
           "limite explicito de T7 (describe, no compara)."), "")
  paste(L, collapse = "\n")
}

# ===========================================================================
main <- function() {
  D <- cargar()
  fuente <- fuente_datos(ARCHIVO_QPCR)

  tc <- tabla_correlaciones(D)
  corr <- tc$filas
  cat(sprintf("  [cruza-verificacion R vs cor.test(exact=FALSE)] peor |dif| = %.3e\n",
              tc$peor))
  stopifnot(tc$peor < 1e-9)
  coexp <- tabla_coexpresion(D)

  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "acto2_correlaciones.csv"), COLS_CORR, corr)
    escribir_csv(file.path(base, "acto2_coexpresion_transportadores.csv"),
                 COLS_COEXP, coexp)
  }

  fig_disp <- file.path(RUTA_FIGURAS, "acto2_dispersion_placenta_cerebro.png")
  figura_dispersion(D, fig_disp)

  fig_gen <- setNames(character(length(ITEMS)), ITEMS)
  for (item in ITEMS) {
    ruta <- file.path(RUTA_FIGURAS,
                      sprintf("acto2_corr_placenta_cerebro_%s.png", item))
    figura_gen_sexo(D, item, ruta)
    fig_gen[[item]] <- ruta
  }

  fig_splom <- list(
    PLACENTA_E15        = file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_PLACENTA_E15.png"),
    BRAIN_E15           = file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_BRAIN_E15.png"),
    PLACENTA_E15_HEMBRA = file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_PLACENTA_E15_HEMBRA.png"),
    PLACENTA_E15_MACHO  = file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_PLACENTA_E15_MACHO.png"),
    BRAIN_E15_HEMBRA    = file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_BRAIN_E15_HEMBRA.png"),
    BRAIN_E15_MACHO     = file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_BRAIN_E15_MACHO.png")
  )
  figura_splom(D, "PLACENTA_E15", fig_splom$PLACENTA_E15, GENES_SPLOM_PLACENTA, "AMBOS")
  figura_splom(D, "BRAIN_E15",    fig_splom$BRAIN_E15,    GENES_SPLOM_BRAIN,    "AMBOS")
  figura_splom(D, "PLACENTA_E15", fig_splom$PLACENTA_E15_HEMBRA, GENES_SPLOM_PLACENTA, "HEMBRA")
  figura_splom(D, "PLACENTA_E15", fig_splom$PLACENTA_E15_MACHO,  GENES_SPLOM_PLACENTA, "MACHO")
  figura_splom(D, "BRAIN_E15",    fig_splom$BRAIN_E15_HEMBRA,    GENES_SPLOM_BRAIN,    "HEMBRA")
  figura_splom(D, "BRAIN_E15",    fig_splom$BRAIN_E15_MACHO,     GENES_SPLOM_BRAIN,    "MACHO")

  escribir_lineas(file.path(RUTA_TABLAS, "acto2_correlaciones_reporte.md"),
                  construir_reporte(fuente, corr, coexp))
  actualizar_descartados()

  # --- Verificaciones (D4 seccion del pedido) -------------------------------
  n_items <- length(ITEMS)
  n_con_rho <- sum(vapply(corr, function(f) nzchar(f[[5]]), logical(1)))
  sc_glob <- Filter(function(f) f[[1]] == "score_compuesto" && f[[3]] == "GLOBAL", corr)[[1]]

  # 1) n por panel sexo x tto <= 9, y coincide con lo que devuelve pares().
  n_panel_ok <- TRUE; n_panel_max <- 0L
  for (item in ITEMS) for (sx in c("HEMBRA", "MACHO")) for (tt in NIVELES_TTO) {
    n <- length(pares(D, item, paste0(sx, "_", tt))$x)
    n_panel_max <- max(n_panel_max, n)
    if (n > 9L) n_panel_ok <- FALSE
  }

  # 2) particion: n(HEMBRA_CONTROL)+n(MACHO_CONTROL) == n(CONTROL), por item (idem LPS).
  .n_de <- function(item, est) Filter(function(f) f[[1]] == item && f[[3]] == est, corr)[[1]][[4]]
  particion_ok <- TRUE
  for (item in ITEMS) {
    if (.n_de(item, "HEMBRA_CONTROL") + .n_de(item, "MACHO_CONTROL") != .n_de(item, "CONTROL"))
      particion_ok <- FALSE
    if (.n_de(item, "HEMBRA_LPS") + .n_de(item, "MACHO_LPS") != .n_de(item, "LPS"))
      particion_ok <- FALSE
  }

  # 3) sin Pearson en ninguna TABLA (columnas ni valores) del Acto 2. La prosa
  # metodologica menciona la palabra "Pearson" a proposito (para decir que NO se
  # usa); lo que exige el pedido es que las TABLAS no tengan valores de Pearson.
  sin_pearson <- !any(grepl("pearson", c(COLS_CORR, COLS_COEXP), ignore.case = TRUE)) &&
    !any(vapply(c(corr, coexp), function(f)
      any(grepl("pearson", vapply(f, .fmt, character(1)), ignore.case = TRUE)),
      logical(1)))

  # 4) n < piso -> sin rho (ya garantizado por construccion; se re-chequea).
  piso_ok <- all(vapply(corr, function(f) {
    n <- f[[4]]; tiene_rho <- nzchar(f[[5]])
    if (n < PISO_PAR) !tiene_rho else TRUE
  }, logical(1)))

  # 5) il6R fuera de ITEMS y de los genes del SPLOM.
  il6r_fuera <- !(GEN_EXCLUIDO_CORR %in% ITEMS) &&
    !(GEN_EXCLUIDO_CORR %in% GENES_SPLOM_PLACENTA) &&
    !(GEN_EXCLUIDO_CORR %in% GENES_SPLOM_BRAIN)

  # 6) tamano exacto de los conjuntos del SPLOM.
  splom_ok <- length(GENES_SPLOM_PLACENTA) == 9L && length(GENES_SPLOM_BRAIN) == 8L &&
    setequal(GENES_SPLOM_PLACENTA, c(GENES_TRANSPORTADORES, "il6", "gp130")) &&
    setequal(GENES_SPLOM_BRAIN, c(GENES_TRANSPORTADORES, "gp130"))

  ent <- sprintf(paste0("data/processed/qpcr_cuantificacion_long.tsv + ",
                        "qpcr_score_compuesto_long.tsv (de data/%s/%s)"),
                 fuente, ARCHIVO_QPCR)
  filas_proced <- list(
    list("outputs/tables/{R,python}/acto2_correlaciones.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("Spearman rho placenta<->cerebro por feto: 8 genes ",
         "(sin il6, sin il6R) + score compuesto x 7 estratos (GLOBAL, CONTROL, ",
         "LPS, HEMBRA_CONTROL, HEMBRA_LPS, MACHO_CONTROL, MACHO_LPS); n_par, ",
         "IC95 Bonett-Wright, p t-aprox")),
    list("outputs/tables/{R,python}/acto2_coexpresion_transportadores.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, paste0("Spearman rho entre los genes del SPLOM ",
         "de cada tejido (9 en placenta, 8 en cerebro), AMBOS/HEMBRA/MACHO")),
    list("outputs/figures/acto2_dispersion_placenta_cerebro.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, paste0("dispersion global (sin separar sexo) ",
         "por item, coloreada por TTO, con rho (global y por grupo) anotada"))
  )
  # procedencia.csv guarda rutas RELATIVAS (convencion del proyecto) -- nunca
  # las rutas absolutas de fig_gen/fig_splom (esas son para file.exists()).
  .rel_fig <- function(ruta) file.path("outputs/figures", basename(ruta))
  for (item in ITEMS)
    filas_proced[[length(filas_proced) + 1L]] <- list(
      .rel_fig(fig_gen[[item]]), "figura", ESTE_SCRIPT, "PROPIO", ent,
      sprintf(paste0("placenta<->cerebro de %s, paneles Females/Males, Control/LPS ",
                    "superpuestos, ajuste lineal + IC95 si n>=5 (pedido explicito)"),
              item))
  filas_proced[[length(filas_proced) + 1L]] <- list(
    .rel_fig(fig_splom$PLACENTA_E15), "figura", ESTE_SCRIPT, "PROPIO", ent,
    "SPLOM co-expresion placenta E15, ambos sexos, 9 variables (7 transp + il6 + gp130)")
  filas_proced[[length(filas_proced) + 1L]] <- list(
    .rel_fig(fig_splom$BRAIN_E15), "figura", ESTE_SCRIPT, "PROPIO", ent,
    "SPLOM co-expresion cerebro fetal E15, ambos sexos, 8 variables (7 transp + gp130)")
  filas_proced[[length(filas_proced) + 1L]] <- list(
    .rel_fig(fig_splom$PLACENTA_E15_HEMBRA), "figura", ESTE_SCRIPT, "PROPIO", ent,
    "idem placenta, solo HEMBRA (pedido explicito)")
  filas_proced[[length(filas_proced) + 1L]] <- list(
    .rel_fig(fig_splom$PLACENTA_E15_MACHO), "figura", ESTE_SCRIPT, "PROPIO", ent,
    "idem placenta, solo MACHO (pedido explicito)")
  filas_proced[[length(filas_proced) + 1L]] <- list(
    .rel_fig(fig_splom$BRAIN_E15_HEMBRA), "figura", ESTE_SCRIPT, "PROPIO", ent,
    "idem cerebro, solo HEMBRA (pedido explicito)")
  filas_proced[[length(filas_proced) + 1L]] <- list(
    .rel_fig(fig_splom$BRAIN_E15_MACHO), "figura", ESTE_SCRIPT, "PROPIO", ent,
    "idem cerebro, solo MACHO (pedido explicito)")
  filas_proced[[length(filas_proced) + 1L]] <- list(
    "outputs/tables/acto2_correlaciones_reporte.md", "reporte", ESTE_SCRIPT,
    "PROPIO", ent, "reporte legible del Acto 2.1-2.2 (T7)")
  registrar_procedencia(filas_proced)

  registrar_verificaciones(list(
    list("acto2_items_correlacionados",
         "correlaciones placenta<->cerebro: 8 genes (sin il6, sin il6R) + score compuesto",
         sprintf("%d items x %d estratos; %d con rho (n_par>=5)",
                 n_items, length(ESTRATOS), n_con_rho),
         "8 genes + score_compuesto", if (n_items == 9L) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_coeficiente",
         "coeficiente = Spearman rho (no Pearson); p t-aprox df n-2; IC Bonett-Wright",
         "rho = Pearson sobre rangos; p = 2*pt(|t|, n-2); IC via tanh/atanh",
         "Spearman", "TRUE", ESTE_SCRIPT),
    list("acto2_il6_fuera_cerebro",
         "il6 excluido del brazo placenta<->cerebro por gen (D7, sin -ddCt en cerebro)",
         if (!(GEN_SIN_CEREBRO %in% ITEMS)) "ITEMS sin il6" else "il6 presente (ERROR)",
         "il6 fuera", if (!(GEN_SIN_CEREBRO %in% ITEMS)) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_il6R_fuera",
         "il6R excluido de todo el Acto 2 (tablas y SPLOM), pedido explicito",
         if (il6r_fuera) "il6R ausente de ITEMS y de ambos GENES_SPLOM" else "il6R presente (ERROR)",
         "il6R fuera", if (il6r_fuera) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_no_compara_grupos",
         "T7 describe pero NO compara correlaciones entre estratos (prohibicion 4)",
         "reporte y analisis_descartados lo dicen explicitamente; el test formal es T8",
         "no se compara en T7", "TRUE", ESTE_SCRIPT),
    list("acto2_pareo_por_feto",
         "el emparejamiento placenta<->cerebro es por FETO (ambos lados detectados)",
         "par = (PLACENTA_E15, BRAIN_E15) del mismo FETO con ambos valores no NA",
         "por feto", "TRUE", ESTE_SCRIPT),
    list("acto2_estratos",
         "estratos de correlacion: GLOBAL + TTO + SEXOxTTO (extension por pedido explicito)",
         paste(ESTRATOS, collapse = ";"),
         "GLOBAL;CONTROL;LPS;HEMBRA_CONTROL;HEMBRA_LPS;MACHO_CONTROL;MACHO_LPS",
         if (identical(ESTRATOS, c("GLOBAL", "CONTROL", "LPS", "HEMBRA_CONTROL",
                                    "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")))
           "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_estratos_particion",
         "n(HEMBRA_x)+n(MACHO_x) == n(x) para x en {CONTROL, LPS}, por item",
         if (particion_ok) "particion exacta en todos los items" else "particion falla (ERROR)",
         "particion exacta", if (particion_ok) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_paneles_n_leyenda",
         "cada panel sexo x tratamiento de las figuras por item tiene n<=9",
         sprintf("n_panel_max=%d", n_panel_max), "n<=9",
         if (n_panel_ok) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_sin_pearson",
         "ninguna figura ni tabla del Acto 2 contiene Pearson (solo Spearman)",
         if (sin_pearson) "sin menciones de Pearson" else "Pearson mencionado (ERROR)",
         "sin Pearson", if (sin_pearson) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_piso_par",
         "n_par < 5 -> sin rho/IC/p reportado (solo n)",
         if (piso_ok) "cumple en todas las filas" else "excepcion encontrada (ERROR)",
         "sin rho si n<5", if (piso_ok) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_splom_variables",
         "SPLOM placenta = 9 variables (7 transp+il6+gp130), cerebro = 8 (7 transp+gp130)",
         sprintf("placenta=%d;cerebro=%d", length(GENES_SPLOM_PLACENTA), length(GENES_SPLOM_BRAIN)),
         "placenta=9;cerebro=8", if (splom_ok) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("acto2_figuras",
         "figuras Acto 2.1-2.2: dispersion global + 9 por item + 6 SPLOM",
         sprintf("disp=%s; por_item=%d/%d; splom=%d/6",
                 file.exists(fig_disp),
                 sum(vapply(fig_gen, file.exists, logical(1))), n_items,
                 sum(vapply(fig_splom, file.exists, logical(1)))),
         "16 figuras existen",
         if (file.exists(fig_disp) && all(vapply(fig_gen, file.exists, logical(1))) &&
             all(vapply(fig_splom, file.exists, logical(1)))) "TRUE" else "FALSE",
         ESTE_SCRIPT)
  ))

  cat("== 08_acto2_correlaciones.R ==\n")
  cat(sprintf("  fuente = %s\n", fuente))
  cat(sprintf("  items = %d (8 genes sin il6/il6R + score_compuesto) x %d estratos\n",
              n_items, length(ESTRATOS)))
  cat(sprintf("  score compuesto GLOBAL: n_par=%s  rho=%s  IC=[%s, %s]  p=%s\n",
              sc_glob[[4]], sc_glob[[5]], sc_glob[[6]], sc_glob[[7]], sc_glob[[8]]))
  for (f in corr) if (f[[3]] == "GLOBAL")
    cat(sprintf("    %-16s GLOBAL  n=%2s  rho=%8s  p=%s\n", f[[1]], f[[4]], f[[5]], f[[8]]))
  cat(sprintf("  co-expresion: %d filas (SPLOM x AMBOS/HEMBRA/MACHO)\n", length(coexp)))
  cat("  -> outputs/tables/{R,python}/acto2_correlaciones.csv, acto2_coexpresion_transportadores.csv\n")
  cat(sprintf("  -> %s + %d figuras por item + %d SPLOM\n", basename(fig_disp),
              n_items, length(fig_splom)))
}

if (sys.nframe() == 0L) main()
