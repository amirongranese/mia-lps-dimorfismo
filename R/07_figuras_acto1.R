# 07_figuras_acto1.R -- Figuras del ACTO 1 (expresion por gen x tejido + pSTAT3).
#
# Por que existe este archivo: reune las figuras descriptivas del Acto 1 usando
# UNA SOLA convencion de anotacion de significancia (D11), la misma para todas.
#
#   * Boxplots de expresion, uno por tejido: en **R base**, NO ggplot2 (cambio
#     pedido explicitamente, pedidos/boxplots_acto1_base_R.R). El estilo pedido
#     necesita controlar caja, bigote y tope del bigote por separado (bigote
#     punteado en tono claro, tope solido, borde de caja fino oscuro) y
#     `geom_boxplot` no expone esos elementos; en R base son argumentos de
#     `boxplot()`/`bxp()`. **El Acto 2 sigue en ggplot2/GGally**, esto es solo
#     para estos paneles. pSTAT3 tambien sigue en ggplot2 (el pedido de estilo
#     de trazo aplicaba solo a "la parte de figuras de expresion").
#       - eje Y = FC = 2^(-ddCt) en escala log (D2). FC se calcula aca, no se
#         guarda en 04 (se evita arrastrar el redondeo de 2^x entre lenguajes).
#       - 4 cajas por panel (HEMBRA_CONTROL, HEMBRA_LPS, MACHO_CONTROL, MACHO_LPS),
#         solo detectados (caja solo si >=3 detectados); puntos individuales.
#       - Filas del panel agrupadas por via metabolica (lipidos/glucosa/
#         aminoacidos/IL-6), no por el orden crudo de GENES -- estilo pedido.
#         Color: Control por sexo (celeste), LPS por sexo x via (una paleta
#         por via metabolica) -- para que el ojo asocie color con la ruta.
#       - il6 @ BRAIN_E15 NO es cuantificable (D7): su panel muestra la PROPORCION
#         DE DETECCION Control vs LPS por sexo (no un boxplot de FC).
#   * Boxplot de pSTAT3 (ggplot2, sin cambios de estilo): valores CRUDOS por
#     SEXO x TTO, MEMBRANA como forma de punto (bloque tecnico). Eje Y lineal.
#   * Las figuras del ELISA (Acto 1.1) ya las produjo 03_elisa: NO se regeneran.
#
# D11 -- anotacion de brackets, AMPLIADA (cambio pedido explicitamente,
# ver AGENTS.md 4.2; revierte la restriccion previa "solo si la interaccion es
# significativa"). Cascada de 3 ramas, en este orden -- solo se entra a UNA:
#   (a) interaccion SEXO x TTO significativa -> brackets por PAR del post hoc
#       D6 (Holm/ART-C), como antes.
#   (b) interaccion NO significativa y efecto principal de TTO significativo o
#       en tendencia -> UN bracket que abarca los 4 grupos, etiqueta
#       "Control vs LPS" + estrellas/p de `p_TTO`. El modelo no sostiene que el
#       efecto difiera entre sexos, asi que marcar pares sugeriria un
#       dimorfismo no sostenido.
#   (c) interaccion NO significativa y efecto principal de SEXO significativo o
#       en tendencia -> UN bracket entre los centros de cada sexo, etiqueta
#       "♀ vs ♂" + estrellas/p de `p_SEXO`. (b) y (c) no son excluyentes entre si.
# Simbolos (las 3 ramas): p<0.001 -> "***" | p<0.01 -> "**" | p<0.05 -> "*"
# (bracket solido); 0.05<=p<0.1 -> bracket punteado + "p = 0.NNN"; p>=0.1 -> nada.
# UNA sola funcion decide que anotar (`d11_brackets_especificacion`, identica en
# R y Python); la geometria de dibujo es propia de cada motor (R base para
# expresion, ggplot2 para pSTAT3, matplotlib en Python para ambas).
#
# PARIDAD: las figuras son PNG -> equivalentes, no byte-identicas (R base +
# ggplot2 vs matplotlib). Byte-identicas entre lenguajes: las filas nuevas de
# `procedencia.csv` / `verificaciones.csv` y la seccion de `analisis_descartados.md`.
# La logica de D11 (`d11_texto` + `d11_brackets_especificacion`) es identica.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({ library(ggplot2) })   # pSTAT3 sigue en ggplot2; expresion en R base

ESTE_SCRIPT <- "07_figuras_acto1"

# --- Paleta y estilo de los boxplots de expresion (estilo pedido) ----------
# COL_CTRL / COL_LPS viven en 00_config.R (fuente unica, no se reescriben
# aca); pSTAT3 (mas abajo) reusa COL_LPS$IL6 -- pSTAT3 es senalizacion de IL-6.

GPATH <- c(fatcd36 = "LIPIDOS", fatp1 = "LIPIDOS", fatp4 = "LIPIDOS",
           glut1 = "GLUCOSA", glut3 = "GLUCOSA",
           slc38a1 = "AMINOACIDOS", slc38a2 = "AMINOACIDOS",
           il6 = "IL6", il6R = "IL6", gp130 = "IL6")
GDISP <- c(fatcd36 = "CD36", fatp1 = "FATP1", fatp4 = "FATP4",
           glut1 = "GLUT1", glut3 = "GLUT3",
           slc38a1 = "SLC38A1", slc38a2 = "SLC38A2",
           il6 = "il6", il6R = "il6R", gp130 = "gp130")
# Orden de filas del panel (estilo pedido): lipidos, glucosa, aminoacidos, IL6
# -- no el orden crudo de GENES. il6R sale del panel de cerebro (deteccion
# insuficiente en ese tejido, D7 solo aplica a il6 pero el pedido tambien
# retira il6R de ESTA figura por prolijidad visual; se conserva en la tabla de
# modelos y en las correlaciones... salvo il6R que ya esta fuera del Acto 2).
FILAS_VIA <- list(
  LIPIDOS     = c("fatcd36", "fatp1", "fatp4"),
  GLUCOSA     = c("glut1", "glut3"),
  AMINOACIDOS = c("slc38a1", "slc38a2"),
  IL6         = c("il6R", "gp130", "il6")
)
gpath <- function(g) GPATH[[g]]
gdisp <- function(g) GDISP[[g]]
color_for <- function(gen, sexo, tto)
  if (tto == "CONTROL") COL_CTRL[[sexo]] else COL_LPS[[gpath(gen)]][[sexo]]
# pSTAT3 (misma logica, sin "gen": usa directo la paleta IL6 -- pSTAT3 es
# senalizacion rio abajo de IL-6, coherente con el color de esa via).
color_pstat3 <- function(sexo, tto)
  if (tto == "CONTROL") COL_CTRL[[sexo]] else COL_LPS$IL6[[sexo]]

# Estilo de trazo (un solo lugar, para que las 2 figuras de expresion coincidan)
EST <- list(
  borde_caja   = "grey20",   # borde fino oscuro
  lwd_caja     = 0.9,
  ancho_caja   = 0.62,       # boxwex: cajas anchas
  col_bigote   = "grey60",   # mas claro que el borde
  lty_bigote   = 2,          # punteado
  lwd_bigote   = 0.9,
  col_tope     = "grey45",   # staple: linea fina solida
  lwd_tope     = 1.1,
  col_mediana  = "grey10",
  lwd_mediana  = 2.2,
  borde_punto  = "grey20",
  lwd_punto    = 0.6,
  cex_punto    = 1.05,
  alfa_relleno = 0.55        # relleno de caja mas tenue que el punto
)

GRUPOS  <- list(c("HEMBRA", "CONTROL"), c("HEMBRA", "LPS"),
                c("MACHO", "CONTROL"),  c("MACHO", "LPS"))
GLAB    <- c("♀ C", "♀ LPS", "♂ C", "♂ LPS")
GRUPOS_4 <- c("HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")
CELDAS_4 <- GRUPOS   # alias (pSTAT3 usa el mismo orden)
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
    # agregar sep + centinela (no vacio) antes de dividir, para que strsplit
    # no descarte un ultimo campo vacio real (le pasaria si solo agregaramos
    # sep): el segmento final termina en el centinela, nunca en "", asi que
    # sobrevive al auto-descarte de strsplit y despues se remueve a mano.
    p <- strsplit(paste0(x, sep, "\001"), sep, fixed = TRUE)[[1]]
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
fila_de <- function(tabla, tej, gen) {
  for (r in tabla) if (r$TEJIDO == tej && r$GEN == gen) return(r)
  NULL
}

# ===========================================================================
# 1. D11 -- LA funcion de anotacion (AMPLIADA: cascada de 3 ramas, AGENTS 4.2).
#    Identica en R y Python. Decide QUE anotar; la geometria de apilado es de
#    cada motor de dibujo (base R / ggplot2 / matplotlib).
# ===========================================================================
d11_texto <- function(p) {
  if (is.null(p) || length(p) != 1L || is.na(p)) return(NULL)
  if (p < 0.001) return(list(texto = "***", estilo = "solida"))
  if (p < 0.01)  return(list(texto = "**",  estilo = "solida"))
  if (p < 0.05)  return(list(texto = "*",   estilo = "solida"))
  if (p < 0.1)   return(list(texto = sprintf("p = %.3f", p), estilo = "punteada"))
  NULL
}

# `fila_clasif`: fila de qpcr_modelos_clasificacion.csv / pstat3_modelo_
# clasificacion.csv (columnas reales: interaccion_significativa, p_TTO, p_SEXO
# -- OJO, NO "p_interaccion"/"metodo": esas columnas no existen, ver
# analisis_descartados.md). `pholm_por_par`: named list etiqueta_D6 -> p_holm.
d11_brackets_especificacion <- function(fila_clasif, pholm_por_par) {
  out <- list()
  isig <- isTRUE(identical(fila_clasif$interaccion_significativa, "TRUE"))
  if (isig) {
    # (a) post hoc D6 -- un bracket por par cuyo p_holm cruce el umbral.
    for (etq in ORDEN_PARES) {
      p <- pholm_por_par[[etq]]
      ann <- d11_texto(if (is.null(p)) NA_real_ else p)
      if (is.null(ann)) next
      idx <- PARES_D6_IDX[[etq]]
      out[[length(out) + 1L]] <- list(x1 = idx[1], x2 = idx[2],
                                      texto = ann$texto, estilo = ann$estilo)
    }
    return(out)
  }
  # (b) efecto principal de TTO -- bracket unico sobre los 4 grupos (0..3).
  p_tto <- .num(fila_clasif$p_TTO)
  ann_t <- d11_texto(p_tto)
  if (!is.null(ann_t)) {
    texto <- if (ann_t$estilo == "solida") sprintf("Control vs LPS %s", ann_t$texto)
             else sprintf("Control vs LPS  %s", ann_t$texto)
    out[[length(out) + 1L]] <- list(x1 = 0, x2 = 3, texto = texto, estilo = ann_t$estilo)
  }
  # (c) efecto principal de SEXO -- bracket entre los centros de cada sexo.
  p_sexo <- .num(fila_clasif$p_SEXO)
  ann_s <- d11_texto(p_sexo)
  if (!is.null(ann_s)) {
    texto <- if (ann_s$estilo == "solida") sprintf("♀ vs ♂ %s", ann_s$texto)
             else sprintf("♀ vs ♂  %s", ann_s$texto)
    out[[length(out) + 1L]] <- list(x1 = 0.5, x2 = 2.5, texto = texto, estilo = ann_s$estilo)
  }
  out
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
pholm_lista <- function(posthoc, tej, gen) {
  out <- list()
  for (r in posthoc)
    if (r$TEJIDO == tej && r$GEN == gen)
      out[[r$contraste]] <- .num(r$p_holm)
  out
}
fc_por_grupo <- function(cuant, tej, gen) {
  out <- setNames(vector("list", 4L), GRUPOS_4)
  for (g in GRUPOS_4) out[[g]] <- numeric(0)
  for (r in cuant) {
    if (r$TEJIDO != tej || r$GEN != gen) next
    if (identical(r$no_detectado, "TRUE") || !nzchar(r$neg_ddCt)) next
    out[[r$GRUPO]] <- c(out[[r$GRUPO]], 2^as.numeric(r$neg_ddCt))
  }
  out
}

# ===========================================================================
# 3. Primitivas de dibujo en R base (boxplots de expresion).
# ===========================================================================
bracket_base <- function(x1, x2, y, etiqueta, solido) {
  # Linea horizontal simple, SIN las perpendiculares en los extremos (pedido
  # explicito: es mas limpio y mas correcto -- marca un efecto que abarca los
  # 4 grupos o los dos sexos, no una comparacion puntual entre dos extremos).
  col <- if (solido) "black" else "grey40"
  lty <- if (solido) 1 else 2
  lines(c(x1, x2), c(y, y), col = col, lty = lty, lwd = 0.9)
  text((x1 + x2) / 2, y, etiqueta, pos = 3, col = col, cex = 0.85)
}

# El eje Y lo fijan los DATOS (D2: escala log), no los brackets -- antes era al
# reves (el techo se estiraba con `ymax * 1.30^(n_niveles+1)`, siempre, incluso
# sin ningun bracket). Si hay brackets, se reserva una FRACCION FIJA del alto
# total del panel en log10 (`frac_reservada`, constante, no crece con la
# cantidad de niveles) repartida en franjas iguales; si no hay brackets, no se
# reserva nada -- el protagonista del panel es el boxplot.
rango_eje_paneles <- function(ymin_datos, ymax_datos, n_niveles,
                              frac_reservada = 0.25, pad_datos = 1.15) {
  piso <- ymin_datos * 0.6
  if (n_niveles <= 0L)
    return(list(piso = piso, techo = ymax_datos * pad_datos, niveles = numeric(0)))
  log_piso <- log10(piso)
  log_datos_top <- log10(ymax_datos * pad_datos)
  rango_datos <- log_datos_top - log_piso
  log_techo <- log_piso + rango_datos / (1 - frac_reservada)
  franja <- (log_techo - log_datos_top) / n_niveles
  niveles <- vapply(seq_len(n_niveles), function(k)
    10 ^ (log_datos_top + (k - 1 + 0.25) * franja), numeric(1))
  list(piso = piso, techo = 10 ^ log_techo, niveles = niveles)
}

# Panel de un gen x tejido. Eje Y SIEMPRE logaritmico (D2): no se decide segun
# el rango de los datos.
panel_gen <- function(fc, gen, tejido, fila_clasif, pholm_por_par, con_titulo = TRUE) {
  dat <- fc[GRUPOS_4]
  todos <- unlist(dat); todos <- todos[is.finite(todos) & todos > 0]
  if (!length(todos)) { plot.new(); return(invisible(NULL)) }

  especs <- d11_brackets_especificacion(fila_clasif, pholm_por_par)
  geo <- rango_eje_paneles(min(todos), max(todos), length(especs))
  ylim <- c(geo$piso, geo$techo)

  rellenos <- vapply(seq_along(GRUPOS), function(j) {
    g <- GRUPOS[[j]]
    grDevices::adjustcolor(color_for(gen, g[1], g[2]), EST$alfa_relleno)
  }, character(1))

  # Caja solo si el grupo tiene >= 3 detectados (los puntos se dibujan igual).
  dat_caja <- lapply(dat, function(v) if (length(v) >= 3L) v else numeric(0))

  par(mar = c(3.6, 4.4, if (con_titulo) 3.8 else 1.2, 0.8))
  boxplot(dat_caja, xaxt = "n", outline = FALSE, log = "y", ylim = ylim,
          ylab = expression(FC == 2^{-Delta*Delta*Ct}),
          boxwex   = EST$ancho_caja,
          col      = rellenos,
          border   = EST$borde_caja,
          boxlwd   = EST$lwd_caja,
          whisklty = EST$lty_bigote,
          whiskcol = EST$col_bigote,
          whisklwd = EST$lwd_bigote,
          staplelty = 1,
          staplecol = EST$col_tope,
          staplelwd = EST$lwd_tope,
          medcol   = EST$col_mediana,
          medlwd   = EST$lwd_mediana)

  axis(1, at = 1:4, labels = GLAB, tick = FALSE, line = -0.4)
  abline(h = 1, lty = 2, col = "grey70", lwd = 0.8)

  for (j in seq_along(GRUPOS)) {
    g <- GRUPOS[[j]]; v <- dat[[GRUPOS_4[j]]]
    if (length(v))
      points(j + (stats::runif(length(v)) - 0.5) * 0.26, v,
             pch = 21, bg = color_for(gen, g[1], g[2]),
             col = EST$borde_punto, lwd = EST$lwd_punto, cex = EST$cex_punto)
    mtext(sprintf("n=%d", length(v)), side = 1, line = 1.5, at = j, cex = 0.62)
  }

  for (k in seq_along(especs)) {
    b <- especs[[k]]
    bracket_base(b$x1 + 1, b$x2 + 1, geo$niveles[k], b$texto, b$estilo == "solida")
  }

  if (con_titulo) {
    p_int <- .num(fila_clasif$p_SEXOxTTO)
    metodo <- if (is.null(fila_clasif$rama_cascada)) "" else fila_clasif$rama_cascada
    title(main = sprintf("%s · %s", gdisp(gen), tejido),
          cex.main = 1.0, font.main = 2, line = 2.4)
    title(main = sprintf("%s · p SEXOxTTO = %s · n = %d",
                         metodo, if (is.na(p_int)) "NA" else sprintf("%.3f", p_int),
                         length(todos)),
          cex.main = 0.72, font.main = 1, line = 1.1)
  }
  invisible(NULL)
}

# Panel de deteccion (il6 @ BRAIN_E15, D7). Barras + simbolos individuales:
# relleno = detectado, vacio = no detectado (misma convencion que ELISA).
panel_deteccion <- function(il6_tab, il6_fis, tej = "BRAIN_E15", gen = "il6") {
  nd <- numeric(4); nt <- numeric(4)
  for (j in seq_along(GRUPOS)) {
    g <- GRUPOS[[j]]
    r <- Filter(function(x) x$TEJIDO == tej && x$GEN == gen &&
                x$GRUPO == GRUPOS_4[j], il6_tab)[[1]]
    nd[j] <- as.numeric(r$n_detectado); nt[j] <- as.numeric(r$n_total)
  }
  pct <- ifelse(nt > 0, 100 * nd / nt, 0)

  par(mar = c(3.6, 4.4, 3.8, 0.8))
  bp <- barplot(pct, ylim = c(0, 125), names.arg = GLAB,
                ylab = "% detectado", border = EST$borde_caja,
                col = vapply(seq_along(GRUPOS), function(j) {
                  g <- GRUPOS[[j]]
                  grDevices::adjustcolor(color_for(gen, g[1], g[2]), EST$alfa_relleno)
                }, character(1)),
                width = 0.62, space = 0.6, las = 1)

  for (j in seq_along(GRUPOS)) {
    g <- GRUPOS[[j]]
    if (nt[j] == 0) next
    xs <- bp[j] + (seq_len(nt[j]) - (nt[j] + 1) / 2) * 0.075
    ys <- rep(pct[j] + 9, nt[j])
    det <- seq_len(nt[j]) <= nd[j]
    points(xs, ys, pch = 21,
           bg = ifelse(det, color_for(gen, g[1], g[2]), NA),
           col = EST$borde_punto, lwd = EST$lwd_punto, cex = 0.85)
    text(bp[j], pct[j] + 18, sprintf("%d/%d\n(%.0f%%)", nd[j], nt[j], pct[j]),
         cex = 0.65, col = "grey20")
  }

  pf <- vapply(c("HEMBRA", "MACHO"), function(s) {
    rr <- Filter(function(x) x$TEJIDO == tej && x$GEN == gen && x$SEXO == s, il6_fis)
    if (length(rr)) .num(rr[[1]]$p_fisher) else NA_real_
  }, numeric(1))

  title(main = sprintf("%s · %s", gdisp(gen), tej), cex.main = 1.0,
        font.main = 2, line = 2.4)
  title(main = sprintf("no cuantificable (D7) · Fisher ♀ p = %.3f · ♂ p = %.3f",
                       pf[["HEMBRA"]], pf[["MACHO"]]),
        cex.main = 0.72, font.main = 1, line = 1.1)
  invisible(NULL)
}

# Figura global de un tejido: filas por via metabolica (estilo pedido).
# Fila 1 lipidos (3) | fila 2 glucosa (2) | fila 3 aminoacidos (2) | fila 4
# via IL-6 (hasta 3, il6@BRAIN_E15 es el panel de deteccion).
#
# Centrado (pedido explicito, las filas de 2 quedaban corridas): grilla de 6
# columnas donde cada panel ocupa 2 -- la fila de 3 usa las columnas 1-2, 3-4
# y 5-6; la de 2 usa 2-3 y 4-5 (centrada dentro del ancho de la fila de 3).
# `layout()` fusiona celdas contiguas con el mismo indice en una sola region.
INICIOS_FILA <- list(`3` = c(1L, 3L, 5L), `2` = c(2L, 4L))
figura_tejido <- function(D, tejido, ruta) {
  filas <- FILAS_VIA
  if (tejido == "BRAIN_E15") filas$IL6 <- setdiff(filas$IL6, "il6R")

  m <- matrix(0L, nrow = length(filas), ncol = 6L)
  k <- 0L
  for (i in seq_along(filas)) {
    gs <- filas[[i]]
    inicios <- INICIOS_FILA[[as.character(length(gs))]]
    if (is.null(inicios)) stop(sprintf("fila con %d paneles no soportada", length(gs)))
    for (j in seq_along(gs)) {
      k <- k + 1L
      m[i, inicios[j] + 0:1] <- k
    }
  }

  png(ruta, width = 3 * 1150, height = length(filas) * 1000, res = DPI)
  on.exit(dev.off(), add = TRUE)
  set.seed(SEMILLA)   # jitter reproducible
  layout(m)
  par(oma = c(0, 0, 3, 0))

  for (i in seq_along(filas)) for (gen in filas[[i]]) {
    if (tejido == "BRAIN_E15" && gen == "il6") {
      panel_deteccion(D$il6_tab, D$il6_fis); next
    }
    fc <- fc_por_grupo(D$cuant, tejido, gen)
    fila_clasif <- fila_de(D$clasif, tejido, gen)
    ph <- pholm_lista(D$posthoc, tejido, gen)
    panel_gen(fc, gen, tejido, fila_clasif, ph)
  }

  mtext(sprintf("Expresion relativa por gen -- %s  (FC = 2^(-ΔΔCt), eje log; brackets = D11)",
                if (tejido == "PLACENTA_E15") "Placenta E15" else "Cerebro fetal E15"),
        outer = TRUE, cex = 1.1, font = 2, line = 0.8)
  invisible(ruta)
}

# ===========================================================================
# 4. Figura de pSTAT3 (ggplot2, sin cambios de estilo -- solo la cascada D11).
# ===========================================================================
apilar_brackets_ggplot <- function(especs, tope, en_log, faceta = NA_character_) {
  if (!length(especs)) return(NULL)
  seg <- list(); txt <- list()
  for (k in seq_along(especs)) {
    b <- especs[[k]]
    paso <- 0.10 + 0.11 * (k - 1)
    if (en_log) { y <- tope * 10^paso; alto <- tope * 10^(paso - 0.035) }
    else        { y <- tope * (1 + paso); alto <- tope * (1 + paso - 0.035) }
    seg[[k]] <- data.frame(faceta = faceta, x = b$x1, xend = b$x2, y = y, yend = y,
                           x0 = b$x1, y0 = alto, x1 = b$x2, y1 = alto,
                           estilo = b$estilo, stringsAsFactors = FALSE)
    txt[[k]] <- data.frame(faceta = faceta, x = (b$x1 + b$x2) / 2, y = y,
                           label = b$texto, estilo = b$estilo, stringsAsFactors = FALSE)
  }
  list(seg = do.call(rbind, seg), txt = do.call(rbind, txt))
}

figura_pstat3 <- function(D, ruta) {
  cl <- if (length(D$pst_cl)) D$pst_cl[[1]] else list()
  pholm <- list()
  for (r in D$pst_ph) pholm[[r$contraste]] <- .num(r$p_holm)

  filas <- list()
  for (j in seq_along(CELDAS_4)) {
    k <- CELDAS_4[[j]]
    for (r in D$pst) if (r$SEXO == k[1] && r$TTO == k[2])
      filas[[length(filas) + 1L]] <- data.frame(
        xi = j - 1L, sexo = k[1], tto = k[2], MEMBRANA = r$MEMBRANA,
        y = as.numeric(r$PSTAT3), stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, filas)
  d$grupo <- paste(d$sexo, d$tto)
  tope <- max(d$y)
  especs <- d11_brackets_especificacion(cl, pholm)
  brk <- apilar_brackets_ggplot(especs, tope, en_log = FALSE)

  # Paleta por grupo SEXO x TTO (misma logica que los boxplots de expresion):
  # celeste (COL_CTRL) para Control, violeta de la via IL-6 (COL_LPS$IL6) para
  # LPS -- pSTAT3 es senalizacion rio abajo de IL-6.
  pal4 <- setNames(
    vapply(CELDAS_4, function(g) color_pstat3(g[1], g[2]), character(1)),
    vapply(CELDAS_4, function(g) paste(g[1], g[2]), character(1)))

  p <- ggplot(d, aes(xi, y, group = xi))
  n_por <- aggregate(y ~ xi + tto, d, length)
  if (any(n_por$y >= 3))
    p <- p + geom_boxplot(aes(colour = grupo, fill = grupo), width = 0.52,
                          outlier.shape = NA, alpha = 0.14, show.legend = FALSE)
  p <- p +
    geom_point(aes(colour = grupo, shape = MEMBRANA),
               position = position_jitter(width = 0.13, height = 0), size = 2) +
    scale_colour_manual(values = pal4, guide = "none") +
    scale_fill_manual(values = pal4, guide = "none") +
    scale_shape_manual(values = c("1" = 16, "2" = 15, "3" = 17),
                       name = NULL, labels = c("Membrana 1", "Membrana 2", "Membrana 3")) +
    scale_x_continuous(breaks = 0:3, labels = ETIQ_X, limits = c(-0.6, 3.6)) +
    ylim(0, tope * 1.75) +
    labs(title = "Fosfo-STAT3 (Tyr705) -- placenta E15",
         subtitle = paste0("crudo por SEXO x TTO; forma de punto = membrana; ",
                           "brackets = D11 (cascada ampliada)"),
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
"### D11 ampliada: cascada de 3 ramas (cambio pedido explicitamente)",
"",
paste0("- Hasta esta sesion, D11 solo anotaba si la interaccion SEXO x TTO era ",
       "significativa; sin interaccion, el panel quedaba sin ninguna marca aunque ",
       "hubiera un efecto principal fuerte (p. ej. `il6@PLACENTA_E15` con ",
       "`p_TTO` = 6.0e-06 no llevaba bracket). **Se revierte esa restriccion por ",
       "pedido explicito** (`pedidos/boxplots_acto1_base_R.R`, tratado como ",
       "especificacion de estilo y logica, no como codigo a copiar). Cascada ",
       "nueva en `d11_brackets_especificacion` (identica R/Python), documentada ",
       "en AGENTS.md 4.2: (a) interaccion significativa -> brackets por par del ",
       "post hoc D6, igual que antes; (b) sin interaccion y `TTO` significativo/",
       "tendencia -> un bracket sobre los 4 grupos, `Control vs LPS`; (c) sin ",
       "interaccion y `SEXO` significativo/tendencia -> un bracket entre los ",
       "centros de cada sexo, `♀ vs ♂`. (b) y (c) no son excluyentes."),
paste0("- **Nombres de columna corregidos** contra los CSV reales al adaptar el ",
       "pedido: la referencia asumia `p_interaccion` y `metodo` en ",
       "`qpcr_modelos_clasificacion.csv`; las columnas reales son ",
       "`p_SEXOxTTO` (usada solo para mostrarla en el subtitulo del panel; el ",
       "gate sigue siendo `interaccion_significativa`) y `rama_cascada`."),
"",
"### Boxplots de expresion: R base, no ggplot2 (cambio pedido explicitamente)",
"",
paste0("- El estilo pedido necesita bigote (punteado, mas claro) y tope del ",
       "bigote (solido) con trazo distinto del borde de la caja; `geom_boxplot` ",
       "no expone esos tres trazos por separado, `boxplot()`/`bxp()` de R base si ",
       "(`whisklty/whiskcol`, `staplelty/staplecol`, `border`). **Se abandona ",
       "ggplot2 solo para estos paneles**; el Acto 2 (08) sigue en ggplot2/GGally, ",
       "y pSTAT3 (mismo script) tambien sigue en ggplot2. matplotlib (Python) ya ",
       "permite estilar los tres trazos por separado (`whiskerprops`/`capprops`/",
       "`boxprops`), asi que este lenguaje no tuvo que cambiar de libreria: solo ",
       "replica el mismo estilo."),
paste0("- **Filas agrupadas por via metabolica** (lipidos / glucosa / aminoacidos ",
       "/ IL-6), no por el orden crudo de `GENES` -- estilo pedido, ayuda a leer ",
       "el panel por sistema biologico. Color: Control por sexo (celeste, 2 ",
       "tonos), LPS por sexo x via metabolica (una paleta por via) -- il6/il6R/",
       "gp130 comparten paleta (via IL-6)."),
"",
"### il6 @ BRAIN_E15: panel de deteccion, no boxplot de FC",
"",
paste0("- il6 en cerebro no es cuantificable (calibrador HEMBRA_CONTROL 0/9, D7): no ",
       "tiene FC. Su panel en `acto1_expresion_BRAIN_E15.png` muestra la **proporcion ",
       "de deteccion** Control vs LPS por sexo (de `qpcr_il6_brain_tabla2x4.csv`) con ",
       "la p de Fisher (de `qpcr_il6_brain_fisher.csv`), no un boxplot. `il6R` sale ",
       "del panel de cerebro por prolijidad visual (deteccion insuficiente en ese ",
       "tejido); se conserva en la tabla de modelos (05) y en el panel de placenta."),
"",
"### Escala y datos de los boxplots de expresion",
"",
paste0("- Eje Y = `FC = 2^(-ddCt)` en escala **log** (D2). `FC` se calcula aca desde ",
       "`neg_ddCt` de `qpcr_cuantificacion_long.tsv`; 04 no lo guarda para no arrastrar ",
       "el redondeo de `2^x` entre libm. Caja solo si el grupo tiene **>=3 detectados** ",
       "(los puntos se dibujan igual, sin caja, si son menos)."),
"",
"### Brackets: linea simple + eje fijado por los datos (pedido post-cierre)",
"",
paste0("- **Bug reportado**: el techo del eje crecia en proporcion a la cantidad de ",
       "brackets, siempre, incluso en paneles sin un solo bracket (`n_niveles` tenia un ",
       "piso de 1) -- `acto1_expresion_PLACENTA_E15.png` llegaba a 10^4 con datos que no ",
       "pasan de 10, `il6R` (sin brackets) a 10^5. **Se separan las dos responsabilidades**: ",
       "la funcion que arma el rango del eje calcula el techo/piso **solo a partir de los ",
       "datos** (mismo padding de siempre, `pad_datos = 1.15`), y **solo si hay brackets** ",
       "reserva una fraccion fija del alto total del panel en log10 (`frac_reservada = 0.25`, ",
       "constante, no crece con la cantidad de niveles), repartida en franjas iguales; ",
       "sin brackets no se reserva nada."),
paste0("- El dibujo del bracket deja de trazar el rectangulo con perpendiculares en los ",
       "extremos: ahora es una **linea horizontal simple** con el texto encima. Mas ",
       "limpio y mas correcto: el bracket marca un efecto que abarca los 4 grupos o los ",
       "dos sexos (ramas (b)/(c) de D11), no una comparacion puntual entre dos extremos."),
"",
"### pSTAT3: valores crudos + membrana como forma de punto",
"",
paste0("- El boxplot de pSTAT3 muestra los valores **crudos** por SEXO x TTO (no ",
       "ajustados por MEMBRANA); la membrana se codifica como forma de punto ",
       "(o / cuadrado / triangulo). El bloque MEMBRANA lo maneja el modelo D9 (06), ",
       "no la figura. Nota al pie: pSTAT3 = abundancia de fosfo-STAT3, no fraccion. ",
       "Sigue en ggplot2; solo cambia la cascada D11 que decide los brackets."),
"",
"### Paleta de pSTAT3 alineada a los boxplots de expresion (pedido post-cierre)",
"",
paste0("- pSTAT3 dejo de usar la paleta Okabe-Ito (`COL_TTO`, azul/naranja generica) ",
       "y pasa a seguir la **misma logica que los boxplots de expresion**: celeste ",
       "(`COL_CTRL`, mas oscuro en macho) para Control, y el violeta de la via IL-6 ",
       "(`COL_LPS$IL6`/`COL_LPS[\"IL6\"]`, `#C5A3E0` hembra / `#7B4EA8` macho) para LPS ",
       "-- coherente porque pSTAT3 es senalizacion rio abajo de IL-6. Helper nuevo ",
       "`color_pstat3(sexo, tto)`, misma logica que `color_for()` pero sin el argumento ",
       "`gen` (siempre usa la via IL6)."),
paste0("- **`COL_CTRL` y `COL_LPS` se mueven a `00_config.{R,py}`** (antes vivian ",
       "hardcodeados en `07_figuras_acto1`): unica fuente de estos colores, sin ",
       "reescribirlos a mano en el script de figuras. `07_figuras_acto1` los referencia ",
       "(`cfg.COL_CTRL`/`cfg.COL_LPS` en Python; `source()` los deja en el mismo entorno ",
       "en R)."),
"",
"### Figuras del ELISA",
"",
paste0("- `acto1_elisa_ms.png` y `acto1_elisa_la.png` (Acto 1.1) ya las produjo ",
       "`03_elisa`; 07 **no las regenera**, solo completan el set del Acto 1."),
"",
"### Paridad",
"",
paste0("- Las figuras son PNG: equivalentes, no byte-identicas (R base + ggplot2 vs ",
       "matplotlib). Byte-identicas entre lenguajes: las filas nuevas de ",
       "`procedencia.csv` / `verificaciones.csv` y esta seccion. La cascada D11 ",
       "(`d11_texto` + `d11_brackets_especificacion`) es identica en ambos.")
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

  figura_tejido(D, "PLACENTA_E15", fig_pla)
  figura_tejido(D, "BRAIN_E15", fig_bra)
  figura_pstat3(D, fig_pst)

  # --- comparaciones/brackets anotados (mismo orden y formato que Python) ---
  .resumen_brackets <- function(tej, gen, fila_clasif, ph) {
    especs <- d11_brackets_especificacion(fila_clasif, ph)
    vapply(especs, function(b) sprintf("%s/%s/[%s-%s]:%s", tej, gen, b$x1, b$x2, b$texto),
           character(1))
  }
  anotadas <- character(0)
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    if (tej == "BRAIN_E15" && gen == "il6") next
    fc <- fila_de(D$clasif, tej, gen)
    if (is.null(fc)) next
    ph <- pholm_lista(D$posthoc, tej, gen)
    anotadas <- c(anotadas, .resumen_brackets(tej, gen, fc, ph))
  }
  pst_cl <- if (length(D$pst_cl)) D$pst_cl[[1]] else list()
  pst_ph <- list(); for (r in D$pst_ph) pst_ph[[r$contraste]] <- .num(r$p_holm)
  pst_especs <- d11_brackets_especificacion(pst_cl, pst_ph)
  pst_anot <- vapply(pst_especs, function(b)
    sprintf("pSTAT3/[%s-%s]:%s", b$x1, b$x2, b$texto), character(1))

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
         "PROPIO", ent_q, paste0("boxplots (R base) de FC = 2^(-ddCt) (eje log) por ",
         "gen en placenta E15, filas por via metabolica, 4 grupos SEXO x TTO; ",
         "cascada D11 ampliada (brackets por efecto principal si no hay interaccion)")),
    list("outputs/figures/acto1_expresion_BRAIN_E15.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent_q, paste0("idem cerebro fetal E15; il6 como panel de proporcion ",
         "de deteccion (D7, no cuantificable), il6R fuera de este panel; cascada D11 ",
         "ampliada (por par si hay interaccion, por efecto principal si no)")),
    list("outputs/figures/acto1_pstat3.png", "figura", ESTE_SCRIPT, "PROPIO", ent_p,
         paste0("boxplot (ggplot2) de pSTAT3 crudo por SEXO x TTO, membrana como forma ",
         "de punto; cascada D11 ampliada del post hoc D6 (06_pstat3); nota de ",
         "limitacion D9"))
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
         "una sola funcion de anotacion D11 (d11_brackets_especificacion) usada en todas las figuras",
         "d11_texto + d11_brackets_especificacion (identica R/Python), 3 figuras",
         "una funcion, todas las figuras", "TRUE", ESTE_SCRIPT),
    list("figuras_d11_cascada",
         "D11 ampliada: (a) interaccion->por par; (b) TTO principal->bracket 0-3; (c) SEXO principal->bracket 0.5-2.5",
         sprintf(paste0("expresion: %d brackets (%s); pSTAT3: %d (%s)"),
                 length(anotadas),
                 if (length(anotadas)) paste(anotadas, collapse = "; ") else "ninguno",
                 length(pst_anot),
                 if (length(pst_anot)) paste(pst_anot, collapse = "; ") else "ninguno"),
         "cascada de 3 ramas, (b)/(c) no excluyentes", "TRUE", ESTE_SCRIPT),
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
         "crudo + forma por membrana", "TRUE", ESTE_SCRIPT),
    list("figuras_expresion_r_base",
         "boxplots de expresion en R base (bxp/boxplot), no ggplot2 -- pedido explicito",
         "boxplot() con whisklty/staplelty/border distintos; Acto 2 sigue en ggplot2",
         "R base para expresion", "TRUE", ESTE_SCRIPT)
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
