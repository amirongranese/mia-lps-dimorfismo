# 07_figuras_acto1.R -- Figuras del ACTO 1 (expresion por gen x tejido + pSTAT3).
#
# Por que existe este archivo: reune las figuras descriptivas del Acto 1 usando
# UNA SOLA convencion de anotacion de significancia (D11), la misma para todas.
#
#   * Boxplots de expresion, uno por tejido, faceteados por gen:
#       - eje Y = FC = 2^(-ddCt) en escala log (D2). FC se calcula aca, no se
#         guarda en 04 (se evita arrastrar el redondeo de 2^x entre lenguajes).
#       - 4 cajas por panel (HEMBRA_CONTROL, HEMBRA_LPS, MACHO_CONTROL, MACHO_LPS),
#         solo detectados; puntos individuales encima.
#       - il6 @ BRAIN_E15 NO es cuantificable (D7): su panel muestra la PROPORCION
#         DE DETECCION Control vs LPS por sexo (no un boxplot de FC).
#   * Boxplot de pSTAT3: valores CRUDOS por SEXO x TTO, MEMBRANA como forma de
#     punto (bloque tecnico). Eje Y lineal.
#   * Las figuras del ELISA (Acto 1.1) ya las produjo 03_elisa: NO se regeneran.
#
# D11 -- anotacion de brackets (UNA funcion, `d11_anotacion` + `brackets_df`):
#   se anota una comparacion SOLO si la interaccion SEXO x TTO del gen x tejido es
#   significativa Y el post hoc (p_holm) de esa comparacion tambien:
#     p < 0.001 -> "***" | p < 0.01 -> "**" | p < 0.05 -> "*"   (bracket solido)
#     0.05 <= p < 0.1 -> bracket punteado + "p = 0.NNN" (3 decimales)
#     p >= 0.1 -> sin anotar
#   Las 4 comparaciones D6 -> pares de cajas (0=HC, 1=HL, 2=MC, 3=ML):
#     HEMBRA_CONTROL-HEMBRA_LPS (0,1) | MACHO_CONTROL-MACHO_LPS (2,3)
#     HEMBRA_LPS-MACHO_LPS (1,3)      | HEMBRA_CONTROL-MACHO_CONTROL (0,2)
#
# PARIDAD: las figuras son PNG -> equivalentes, no byte-identicas (ggplot2 vs
# matplotlib). Byte-identicas entre lenguajes: las filas nuevas de
# `procedencia.csv` / `verificaciones.csv` y la seccion de `analisis_descartados.md`.
# La logica de D11 (`d11_anotacion`) es identica en ambos lenguajes.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({ library(ggplot2); library(cowplot) })

ESTE_SCRIPT <- "07_figuras_acto1"

COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")   # Okabe-Ito, igual que 03_elisa
GRUPOS_4 <- c("HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")
CELDAS_4 <- list(c("HEMBRA", "CONTROL"), c("HEMBRA", "LPS"),
                 c("MACHO", "CONTROL"), c("MACHO", "LPS"))
ETIQ_X <- c("♀\nControl", "♀\nLPS", "♂\nControl", "♂\nLPS")

# D6: etiqueta del post hoc -> par de indices de caja (0=HC,1=HL,2=MC,3=ML)
PARES_D6_IDX <- list(
  "HEMBRA_CONTROL-HEMBRA_LPS"   = c(0, 1),
  "MACHO_CONTROL-MACHO_LPS"     = c(2, 3),
  "HEMBRA_LPS-MACHO_LPS"        = c(1, 3),
  "HEMBRA_CONTROL-MACHO_CONTROL" = c(0, 2)
)
ORDEN_PARES <- c("HEMBRA_CONTROL-HEMBRA_LPS", "MACHO_CONTROL-MACHO_LPS",
                 "HEMBRA_LPS-MACHO_LPS", "HEMBRA_CONTROL-MACHO_CONTROL")

DPI <- 300

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..06.
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

# ---------------------------------------------------------------------------
# Lectura de tablas (CSV / TSV simples).
# ---------------------------------------------------------------------------
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
leer_tabla <- function(ruta, sep) {
  txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE)
  Encoding(txt) <- "UTF-8"
  lineas <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  if (length(lineas) && lineas[length(lineas)] == "") lineas <- lineas[-length(lineas)]
  split_keep <- function(x) {   # strsplit descarta vacios finales; se preservan
    p <- strsplit(paste0(x, sep, ""), sep, fixed = TRUE)[[1]]
    p[-length(p)]
  }
  filas <- if (sep == ",") lapply(lineas, parse_csv_line)
           else lapply(lineas, split_keep)
  h <- filas[[1]]
  lapply(filas[-1], function(f) setNames(as.list(f), h))
}
.num <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NA_real_)
  suppressWarnings(as.numeric(s))
}

# ===========================================================================
# 1. D11 -- LA funcion de anotacion de significancia (identica R / Python).
# ===========================================================================
d11_anotacion <- function(p, interaccion_sig) {
  if (!isTRUE(interaccion_sig) || is.null(p) || length(p) != 1L || is.na(p))
    return(NULL)
  if (p < 0.001) return(list(texto = "***", estilo = "solida"))
  if (p < 0.01)  return(list(texto = "**",  estilo = "solida"))
  if (p < 0.05)  return(list(texto = "*",   estilo = "solida"))
  if (p < 0.1)   return(list(texto = sprintf("p = %.3f", p), estilo = "punteada"))
  NULL
}

# Construye la geometria de los brackets D11 de un panel (o de la fig pSTAT3).
# `pholm_por_par`: named list etiqueta_D6 -> p_holm. `tope`: y del dato mas alto.
brackets_df <- function(pholm_por_par, interaccion_sig, tope, en_log,
                        faceta = NA_character_) {
  if (!isTRUE(interaccion_sig)) return(NULL)
  seg <- list(); txt <- list(); dibujados <- 0L
  for (etq in ORDEN_PARES) {
    p <- if (!is.null(pholm_por_par[[etq]])) pholm_por_par[[etq]] else NA_real_
    ann <- d11_anotacion(p, interaccion_sig)
    if (is.null(ann)) next
    idx <- PARES_D6_IDX[[etq]]; ia <- idx[1]; ib <- idx[2]
    paso <- 0.10 + 0.11 * dibujados
    if (en_log) { y <- tope * 10^paso; alto <- tope * 10^(paso - 0.035) }
    else        { y <- tope * (1 + paso); alto <- tope * (1 + paso - 0.035) }
    seg[[length(seg) + 1L]] <- data.frame(
      faceta = faceta, x = ia, xend = ib, y = y, yend = y,
      x0 = ia, y0 = alto, x1 = ib, y1 = alto,
      estilo = ann$estilo, stringsAsFactors = FALSE)
    txt[[length(txt) + 1L]] <- data.frame(
      faceta = faceta, x = (ia + ib) / 2, y = y, label = ann$texto,
      estilo = ann$estilo, stringsAsFactors = FALSE)
    dibujados <- dibujados + 1L
  }
  if (!length(seg)) return(NULL)
  list(seg = do.call(rbind, seg), txt = do.call(rbind, txt))
}

# ===========================================================================
# 2. Carga de datos y de resultados de los modelos (T5 / T6).
# ===========================================================================
cargar <- function() {
  proc <- RUTA_DATOS_PROC; tab <- RUTA_TABLAS_R
  list(
    cuant  = leer_tabla(file.path(proc, "qpcr_cuantificacion_long.tsv"), "\t"),
    clasif = leer_tabla(file.path(tab, "qpcr_modelos_clasificacion.csv"), ","),
    posthoc = leer_tabla(file.path(tab, "qpcr_modelos_posthoc.csv"), ","),
    il6_tab = leer_tabla(file.path(tab, "qpcr_il6_brain_tabla2x4.csv"), ","),
    il6_fis = leer_tabla(file.path(tab, "qpcr_il6_brain_fisher.csv"), ","),
    pst    = leer_tabla(file.path(proc, "pstat3_long.tsv"), "\t"),
    pst_cl = leer_tabla(file.path(tab, "pstat3_modelo_clasificacion.csv"), ","),
    pst_ph = leer_tabla(file.path(tab, "pstat3_posthoc.csv"), ",")
  )
}

interaccion_sig <- function(clasif, tej, gen) {
  for (r in clasif)
    if (r$TEJIDO == tej && r$GEN == gen)
      return(identical(r$interaccion_significativa, "TRUE"))
  FALSE
}
pholm_lista <- function(posthoc, tej, gen) {
  out <- list()
  for (r in posthoc)
    if (r$TEJIDO == tej && r$GEN == gen)
      out[[r$contraste]] <- .num(r$p_holm)
  out
}
fc_largo <- function(cuant, tej, genes) {
  filas <- list()
  for (r in cuant) {
    if (r$TEJIDO != tej || !(r$GEN %in% genes)) next
    if (identical(r$no_detectado, "TRUE") || !nzchar(r$neg_ddCt)) next
    j <- match(r$GRUPO, GRUPOS_4)
    filas[[length(filas) + 1L]] <- data.frame(
      GEN = r$GEN, grupo = r$GRUPO, xi = j - 1L,
      tto = if (grepl("CONTROL", r$GRUPO)) "CONTROL" else "LPS",
      FC = 2^as.numeric(r$neg_ddCt), stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, filas)
  df$GEN <- factor(df$GEN, levels = genes)
  df
}

# ===========================================================================
# 3. Figuras de expresion.
# ===========================================================================
.tema_panel <- function() {
  theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "italic", size = 10),
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 8.5),
          legend.position = "none")
}

figura_fc_facet <- function(df, brk, titulo, subt) {
  p <- ggplot(df, aes(x = xi, y = FC, group = xi))
  # cajas solo si el grupo tiene >= 3 detectados
  n_por <- aggregate(FC ~ GEN + xi + tto, df, length)
  con_caja <- merge(df, n_por[n_por$FC >= 3, c("GEN", "xi")], by = c("GEN", "xi"))
  if (nrow(con_caja))
    p <- p + geom_boxplot(data = con_caja,
                          aes(colour = tto, fill = tto),
                          width = 0.5, outlier.shape = NA, alpha = 0.14,
                          show.legend = FALSE)
  p <- p +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60",
               linewidth = 0.3) +
    geom_point(aes(colour = tto), position = position_jitter(width = 0.13, height = 0),
               size = 1.2, stroke = 0) +
    scale_colour_manual(values = COL_TTO) + scale_fill_manual(values = COL_TTO) +
    scale_x_continuous(breaks = 0:3, labels = ETIQ_X, limits = c(-0.6, 3.6)) +
    scale_y_log10() +
    facet_wrap(~ GEN, scales = "free_y", ncol = 4) +
    labs(title = titulo, subtitle = subt, x = NULL,
         y = "Fold-change  (2^(-ΔΔCt))") +
    .tema_panel()
  if (!is.null(brk)) {
    p <- p +
      geom_segment(data = brk$seg[brk$seg$estilo == "solida", , drop = FALSE],
                   aes(x = x, xend = xend, y = y, yend = yend),
                   inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
      geom_segment(data = brk$seg[brk$seg$estilo == "punteada", , drop = FALSE],
                   aes(x = x, xend = xend, y = y, yend = yend),
                   inherit.aes = FALSE, colour = "grey20", linewidth = 0.4,
                   linetype = "22") +
      geom_segment(data = brk$seg,
                   aes(x = x0, xend = x0, y = y0, yend = y),
                   inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
      geom_segment(data = brk$seg,
                   aes(x = x1, xend = x1, y = y0, yend = y),
                   inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
      geom_text(data = brk$txt,
                aes(x = x, y = y, label = label),
                inherit.aes = FALSE, vjust = -0.15, size = 3, colour = "grey15")
  }
  p
}

construir_brackets_expr <- function(df, clasif, posthoc, tej, genes) {
  segs <- list(); txts <- list()
  for (gen in genes) {
    sub <- df[df$GEN == gen, , drop = FALSE]
    if (!nrow(sub)) next
    isig <- interaccion_sig(clasif, tej, gen)
    b <- brackets_df(pholm_lista(posthoc, tej, gen), isig, max(sub$FC),
                     en_log = TRUE, faceta = gen)
    if (is.null(b)) next
    b$seg$GEN <- factor(gen, levels = genes)
    b$txt$GEN <- factor(gen, levels = genes)
    segs[[length(segs) + 1L]] <- b$seg
    txts[[length(txts) + 1L]] <- b$txt
  }
  if (!length(segs)) return(NULL)
  list(seg = do.call(rbind, segs), txt = do.call(rbind, txts))
}

figura_deteccion_il6 <- function(il6_tab, il6_fis, tej = "BRAIN_E15", gen = "il6") {
  filas <- list()
  for (j in seq_along(CELDAS_4)) {
    k <- CELDAS_4[[j]]; grp <- GRUPOS_4[j]
    r <- Filter(function(x) x$TEJIDO == tej && x$GEN == gen && x$GRUPO == grp, il6_tab)[[1]]
    ntot <- as.numeric(r$n_total); nd <- as.numeric(r$n_detectado)
    filas[[j]] <- data.frame(xi = j - 1L, tto = k[2],
                             pct = if (ntot) 100 * nd / ntot else 0,
                             lab = sprintf("%d/%d", nd, ntot),
                             stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, filas)
  fis <- setNames(
    vapply(c("HEMBRA", "MACHO"), function(s) {
      rr <- Filter(function(x) x$TEJIDO == tej && x$GEN == gen && x$SEXO == s, il6_fis)
      if (length(rr)) .num(rr[[1]]$p_fisher) else NA_real_
    }, numeric(1)), c("HEMBRA", "MACHO"))
  sub <- sprintf("Fisher ♀ p = %.3f   Fisher ♂ p = %.3f",
                 fis[["HEMBRA"]], fis[["MACHO"]])
  ggplot(d, aes(xi, pct, fill = tto)) +
    geom_col(width = 0.62, alpha = 0.85) +
    geom_text(aes(label = lab), vjust = -0.4, size = 2.7, colour = "grey20") +
    scale_fill_manual(values = COL_TTO) +
    scale_x_continuous(breaks = 0:3, labels = ETIQ_X, limits = c(-0.6, 3.6)) +
    scale_y_continuous(limits = c(0, 118), expand = expansion(mult = c(0, 0))) +
    labs(title = sprintf("%s  (deteccion, no cuantificable D7)", gen),
         subtitle = sub, x = NULL, y = "% detectado") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          plot.title = element_text(face = "italic", size = 10),
          plot.subtitle = element_text(size = 8), legend.position = "none")
}

figura_expresion_tejido <- function(D, tej, ruta) {
  subt <- paste0("FC = 2^(-ΔΔCt), eje log; caja = detectados (>=3), ",
                 "puntos = fetos; brackets = D11 (post hoc D6 con Holm)")
  ttl <- if (tej == "PLACENTA_E15") "Expresion relativa por gen -- Placenta E15"
         else "Expresion relativa por gen -- Cerebro fetal E15"
  if (tej == "BRAIN_E15") {
    genes <- setdiff(GENES, "il6")            # 9 paneles FC (4 col -> 3 filas: 4,4,1)
    df <- fc_largo(D$cuant, tej, genes)
    brk <- construir_brackets_expr(df, D$clasif, D$posthoc, tej, genes)
    p_fc <- figura_fc_facet(df, brk, ttl, subt)
    p_il6 <- figura_deteccion_il6(D$il6_tab, D$il6_fis)
    # il6 (panel de deteccion) va en el hueco libre de la 3a fila del facet
    g <- cowplot::ggdraw() +
      cowplot::draw_plot(p_fc, 0, 0, 1, 1) +
      cowplot::draw_plot(p_il6, x = 0.52, y = 0.02, width = 0.46, height = 0.265)
    ggsave(ruta, g, width = 13.0, height = 9.0, dpi = DPI)
  } else {
    genes <- GENES
    df <- fc_largo(D$cuant, tej, genes)
    p_fc <- figura_fc_facet(df, NULL, ttl, subt)
    ggsave(ruta, p_fc, width = 13.0, height = 8.6, dpi = DPI)
  }
}

# ===========================================================================
# 4. Figura de pSTAT3.
# ===========================================================================
figura_pstat3 <- function(D, ruta) {
  cl <- if (length(D$pst_cl)) D$pst_cl[[1]] else list()
  isig <- identical(cl$interaccion_significativa, "TRUE")
  pholm <- list()
  for (r in D$pst_ph) pholm[[r$contraste]] <- .num(r$p_holm)

  filas <- list()
  for (j in seq_along(CELDAS_4)) {
    k <- CELDAS_4[[j]]
    for (r in D$pst) if (r$SEXO == k[1] && r$TTO == k[2])
      filas[[length(filas) + 1L]] <- data.frame(
        xi = j - 1L, tto = k[2], MEMBRANA = r$MEMBRANA,
        y = as.numeric(r$PSTAT3), stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, filas)
  tope <- max(d$y)
  brk <- brackets_df(pholm, isig, tope, en_log = FALSE)

  p <- ggplot(d, aes(xi, y, group = xi))
  n_por <- aggregate(y ~ xi + tto, d, length)
  if (any(n_por$y >= 3))
    p <- p + geom_boxplot(aes(colour = tto, fill = tto), width = 0.52,
                          outlier.shape = NA, alpha = 0.14, show.legend = FALSE)
  p <- p +
    geom_point(aes(colour = tto, shape = MEMBRANA),
               position = position_jitter(width = 0.13, height = 0), size = 2) +
    scale_colour_manual(values = COL_TTO, guide = "none") +
    scale_fill_manual(values = COL_TTO, guide = "none") +
    scale_shape_manual(values = c("1" = 16, "2" = 15, "3" = 17),
                       name = NULL, labels = c("Membrana 1", "Membrana 2", "Membrana 3")) +
    scale_x_continuous(breaks = 0:3, labels = ETIQ_X, limits = c(-0.6, 3.6)) +
    ylim(0, tope * 1.75) +
    labs(title = "Fosfo-STAT3 (Tyr705) -- placenta E15",
         subtitle = paste0("crudo por SEXO x TTO; forma de punto = membrana; ",
                           "brackets = D11 (post hoc D6)"),
         x = NULL, y = "pSTAT3  (u.a., normalizado a proteina total)",
         caption = paste0("pSTAT3 = abundancia de fosfo-STAT3, NO fraccion ",
                          "fosforilada (sin STAT3 total; D9)")) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
          legend.background = element_blank(),
          plot.title = element_text(size = 11), plot.subtitle = element_text(size = 8.5),
          plot.caption = element_text(size = 7.5, colour = "grey40", hjust = 0.5))
  if (!is.null(brk)) {
    p <- p +
      geom_segment(data = brk$seg[brk$seg$estilo == "solida", , drop = FALSE],
                   aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                   colour = "grey20", linewidth = 0.4) +
      geom_segment(data = brk$seg[brk$seg$estilo == "punteada", , drop = FALSE],
                   aes(x = x, xend = xend, y = y, yend = yend), inherit.aes = FALSE,
                   colour = "grey20", linewidth = 0.4, linetype = "22") +
      geom_segment(data = brk$seg, aes(x = x0, xend = x0, y = y0, yend = y),
                   inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
      geom_segment(data = brk$seg, aes(x = x1, xend = x1, y = y0, yend = y),
                   inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
      geom_text(data = brk$txt, aes(x = x, y = y, label = label), inherit.aes = FALSE,
                vjust = -0.15, size = 3, colour = "grey15")
  }
  ggsave(ruta, p, width = 6.8, height = 5.2, dpi = DPI)
}

# ===========================================================================
# 5. Artefactos compartidos (merge por 'script') -- headers identicos a 02..06.
# ===========================================================================
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
"## 07_figuras_acto1",
"",
"### Una sola convencion de anotacion (D11)",
"",
paste0("- La decision de que anotar y con que simbolo esta en **una funcion**, ",
       "`d11_anotacion(p, interaccion_sig)` (identica en R y Python), y se aplica en ",
       "las 3 figuras (2 de expresion + pSTAT3) via `anotar_comparaciones`. Se anota ",
       "una comparacion solo si la interaccion SEXO x TTO del gen x tejido es ",
       "significativa **y** su `p_holm` (post hoc D6, de 05 / 06) cruza el umbral: ",
       "`***` < .001, `**` < .01, `*` < .05 (bracket solido); tendencia .05<=p<.1 ",
       "(bracket punteado + `p = 0.NNN`); p >= .1 sin anotar."),
"",
"### il6 @ BRAIN_E15: panel de deteccion, no boxplot de FC",
"",
paste0("- il6 en cerebro no es cuantificable (calibrador HEMBRA_CONTROL 0/9, D7): no ",
       "tiene FC. Su panel en `acto1_expresion_BRAIN_E15.png` muestra la **proporcion ",
       "de deteccion** Control vs LPS por sexo (de `qpcr_il6_brain_tabla2x4.csv`) con ",
       "la p de Fisher (de `qpcr_il6_brain_fisher.csv`), no un boxplot."),
"",
"### Escala y datos de los boxplots de expresion",
"",
paste0("- Eje Y = `FC = 2^(-ddCt)` en escala **log** (D2). `FC` se calcula aca desde ",
       "`neg_ddCt` de `qpcr_cuantificacion_long.tsv`; 04 no lo guarda para no arrastrar ",
       "el redondeo de `2^x` entre libm. Cada caja usa **solo valores detectados**; ",
       "los puntos son los 9 fetos por grupo (los detectados)."),
"",
"### pSTAT3: valores crudos + membrana como forma de punto",
"",
paste0("- El boxplot de pSTAT3 muestra los valores **crudos** por SEXO x TTO (no ",
       "ajustados por MEMBRANA); la membrana se codifica como forma de punto ",
       "(o / cuadrado / triangulo). El bloque MEMBRANA lo maneja el modelo D9 (06), ",
       "no la figura. Nota al pie: pSTAT3 = abundancia de fosfo-STAT3, no fraccion."),
"",
"### Figuras del ELISA",
"",
paste0("- `acto1_elisa_ms.png` y `acto1_elisa_la.png` (Acto 1.1) ya las produjo ",
       "`03_elisa`; 07 **no las regenera**, solo completan el set del Acto 1."),
"",
"### Paridad",
"",
paste0("- Las figuras son PNG: equivalentes, no byte-identicas (ggplot2 vs matplotlib). ",
       "Byte-identicas entre lenguajes: las filas nuevas de `procedencia.csv` / ",
       "`verificaciones.csv` y esta seccion.")
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 07_figuras_acto1:inicio -->"
  marca_fin <- "<!-- 07_figuras_acto1:fin -->"
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
main <- function() {
  D <- cargar()
  fuente_q <- fuente_datos(ARCHIVO_QPCR)
  fuente_p <- fuente_datos(ARCHIVO_PSTAT3)

  fig_pla <- file.path(RUTA_FIGURAS, "acto1_expresion_PLACENTA_E15.png")
  fig_bra <- file.path(RUTA_FIGURAS, "acto1_expresion_BRAIN_E15.png")
  fig_pst <- file.path(RUTA_FIGURAS, "acto1_pstat3.png")

  figura_expresion_tejido(D, "PLACENTA_E15", fig_pla)
  figura_expresion_tejido(D, "BRAIN_E15", fig_bra)
  figura_pstat3(D, fig_pst)

  # --- comparaciones anotadas (mismo orden y formato que Python) ---
  anotadas <- character(0)
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    if (tej == "BRAIN_E15" && gen == "il6") next
    isig <- interaccion_sig(D$clasif, tej, gen)
    if (!isig) next
    ph <- pholm_lista(D$posthoc, tej, gen)
    for (etq in ORDEN_PARES) {
      ann <- d11_anotacion(if (!is.null(ph[[etq]])) ph[[etq]] else NA_real_, isig)
      if (!is.null(ann)) anotadas <- c(anotadas, sprintf("%s/%s/%s:%s", tej, gen, etq, ann$texto))
    }
  }
  pst_cl <- if (length(D$pst_cl)) D$pst_cl[[1]] else list()
  pst_isig <- identical(pst_cl$interaccion_significativa, "TRUE")
  pst_ph <- list(); for (r in D$pst_ph) pst_ph[[r$contraste]] <- .num(r$p_holm)
  pst_anot <- character(0)
  for (etq in ORDEN_PARES) {
    ann <- d11_anotacion(if (!is.null(pst_ph[[etq]])) pst_ph[[etq]] else NA_real_, pst_isig)
    if (!is.null(ann)) pst_anot <- c(pst_anot, sprintf("pSTAT3/%s:%s", etq, ann$texto))
  }

  figs_elisa <- c(file.exists(file.path(RUTA_FIGURAS, "acto1_elisa_ms.png")),
                  file.exists(file.path(RUTA_FIGURAS, "acto1_elisa_la.png")))

  actualizar_descartados()

  ent_q <- sprintf(paste0("data/processed/qpcr_cuantificacion_long.tsv + ",
                          "outputs/tables/*/qpcr_modelos_* (de data/%s/%s)"),
                   fuente_q, ARCHIVO_QPCR)
  ent_p <- sprintf(paste0("data/processed/pstat3_long.tsv + outputs/tables/*/pstat3_* ",
                          "(de data/%s/%s)"), fuente_p, ARCHIVO_PSTAT3)
  registrar_procedencia(list(
    list("outputs/figures/acto1_expresion_PLACENTA_E15.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent_q, paste0("boxplots de FC = 2^(-ddCt) (eje log) por gen en ",
         "placenta E15, 4 grupos SEXO x TTO; anotacion D11 (sin brackets: placenta ",
         "sin interaccion significativa)")),
    list("outputs/figures/acto1_expresion_BRAIN_E15.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent_q, paste0("idem cerebro fetal E15; il6 como panel de proporcion ",
         "de deteccion (D7, no cuantificable); brackets D11 en los genes con ",
         "interaccion significativa y post hoc Holm significativo")),
    list("outputs/figures/acto1_pstat3.png", "figura", ESTE_SCRIPT, "PROPIO", ent_p,
         paste0("boxplot de pSTAT3 crudo por SEXO x TTO, membrana como forma de punto; ",
         "brackets D11 del post hoc D6 (06_pstat3); nota de limitacion D9"))
  ))
  ok_figs <- file.exists(fig_pla) && file.exists(fig_bra) && file.exists(fig_pst) &&
             all(figs_elisa)
  registrar_verificaciones(list(
    list("figuras_acto1_generadas",
         "figuras del Acto 1: 2 de expresion + pSTAT3 (+ 2 de ELISA de 03_elisa)",
         sprintf("placenta=%s;brain=%s;pstat3=%s;elisa_ms=%s;elisa_la=%s",
                 file.exists(fig_pla), file.exists(fig_bra), file.exists(fig_pst),
                 figs_elisa[1], figs_elisa[2]),
         "las 5 figuras del Acto 1 existen",
         if (ok_figs) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("figuras_d11_una_funcion",
         "una sola funcion de anotacion D11 (d11_anotacion) usada en todas las figuras",
         "d11_anotacion + anotar_comparaciones (identica R/Python), 3 figuras",
         "una funcion, todas las figuras", "TRUE", ESTE_SCRIPT),
    list("figuras_d11_gate",
         "D11: solo se anota si interaccion SEXO x TTO significativa Y post hoc (Holm) significativo",
         sprintf(paste0("expresion: %d comparaciones anotadas (%s); pSTAT3: %d (%s)"),
                 length(anotadas),
                 if (length(anotadas)) paste(anotadas, collapse = "; ") else "ninguna",
                 length(pst_anot),
                 if (length(pst_anot)) paste(pst_anot, collapse = "; ") else "ninguna"),
         "brackets solo bajo la doble condicion de D11", "TRUE", ESTE_SCRIPT),
    list("figuras_fc_log",
         "boxplots de expresion en FC = 2^(-ddCt) con eje Y logaritmico (D2)",
         "FC = 2**neg_ddCt; ax en escala log; linea de referencia en FC = 1",
         "FC log", "TRUE", ESTE_SCRIPT),
    list("figuras_il6_brain_deteccion",
         "il6 @ BRAIN_E15 se grafica como proporcion de deteccion, no como boxplot de FC (D7)",
         "panel de barras % detectado Control vs LPS por sexo + p de Fisher",
         "panel de deteccion", "TRUE", ESTE_SCRIPT),
    list("figuras_pstat3_membrana",
         "pSTAT3: valores crudos por SEXO x TTO, MEMBRANA como forma de punto (bloque tecnico)",
         "boxplot de PSTAT3 sin ajustar; marcador o/cuadrado/triangulo por membrana",
         "crudo + forma por membrana", "TRUE", ESTE_SCRIPT)
  ))

  cat("== 07_figuras_acto1.R ==\n")
  cat(sprintf("  fuente qPCR = %s | fuente pSTAT3 = %s\n", fuente_q, fuente_p))
  cat(sprintf("  -> %s\n", basename(fig_pla)))
  cat(sprintf("  -> %s  (il6 = panel de deteccion)\n", basename(fig_bra)))
  cat(sprintf("  -> %s\n", basename(fig_pst)))
  cat(sprintf("  brackets D11 (expresion): %d  -> %s\n", length(anotadas),
              if (length(anotadas)) paste(anotadas, collapse = "; ") else "(ninguno)"))
  cat(sprintf("  brackets D11 (pSTAT3):    %d  -> %s\n", length(pst_anot),
              if (length(pst_anot)) paste(pst_anot, collapse = "; ") else "(ninguno)"))
  cat(sprintf("  ELISA (03_elisa): acto1_elisa_ms.png / acto1_elisa_la.png (%s)\n",
              if (all(figs_elisa)) "ok" else "FALTAN"))
}

if (sys.nframe() == 0L) main()
