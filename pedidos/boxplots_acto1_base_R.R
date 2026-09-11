# =============================================================================
# BOXPLOTS DEL ACTO 1 -- implementacion de referencia en R base.
#
# Reemplaza la parte de figuras de expresion de 07_figuras_acto1.R.
# Se abandona ggplot2 para estos paneles por una razon concreta: el estilo
# pedido necesita controlar caja, bigote y tope del bigote por separado
# (bigote punteado en tono claro, tope solido, borde de caja fino oscuro), y
# geom_boxplot no expone esos elementos. En R base son argumentos de bxp.
# Las figuras de correlacion del Acto 2 siguen en ggplot2.
#
# NO PROBADO: escrito sin acceso a R ni a los datos. Verificar los nombres de
# columna marcados abajo contra los CSV reales antes de correr.
# =============================================================================

# --- A CONFIRMAR contra outputs/tables/R/qpcr_modelos_clasificacion.csv ------
COL_P_INT    <- "p_interaccion"
COL_P_TTO    <- "p_TTO"
COL_P_SEXO   <- "p_SEXO"
COL_METODO   <- "metodo"
COL_INT_SIG  <- "interaccion_significativa"

ALPHA        <- 0.05
ALPHA_TREND  <- 0.10

# =============================================================================
# 1. Paleta y clasificacion de genes.
#    Estas constantes van en 00_config.R, no aca: las usan 07 y 08.
# =============================================================================
COL_CTRL <- c(HEMBRA = "#AEDCF0", MACHO = "#6BAED6")
COL_LPS  <- list(
  GLUCOSA     = c(HEMBRA = "#F09EC8", MACHO = "#D6317F"),
  AMINOACIDOS = c(HEMBRA = "#8FD9B6", MACHO = "#2E9E6B"),
  LIPIDOS     = c(HEMBRA = "#FBC98A", MACHO = "#E08214"),
  IL6         = c(HEMBRA = "#C5A3E0", MACHO = "#7B4EA8")
)

GPATH <- c(fatcd36 = "LIPIDOS", fatp1 = "LIPIDOS", fatp4 = "LIPIDOS",
           glut1 = "GLUCOSA", glut3 = "GLUCOSA",
           slc38a1 = "AMINOACIDOS", slc38a2 = "AMINOACIDOS",
           il6 = "IL6", il6R = "IL6", gp130 = "IL6")

GDISP <- c(fatcd36 = "CD36", fatp1 = "FATP1", fatp4 = "FATP4",
           glut1 = "GLUT1", glut3 = "GLUT3",
           slc38a1 = "SLC38A1", slc38a2 = "SLC38A2",
           il6 = "il6", il6R = "il6R", gp130 = "gp130")

# Orden de filas del panel global (pedido): lipidos, glucosa, aminoacidos, IL6.
FILAS_VIA <- list(
  LIPIDOS     = c("fatcd36", "fatp1", "fatp4"),
  GLUCOSA     = c("glut1", "glut3"),
  AMINOACIDOS = c("slc38a1", "slc38a2"),
  IL6         = c("il6R", "gp130", "il6")
)

gpath <- function(g) if (tolower(g) %in% tolower(names(GPATH)))
  GPATH[[names(GPATH)[tolower(names(GPATH)) == tolower(g)]]] else "OTRO"
gdisp <- function(g) if (tolower(g) %in% tolower(names(GDISP)))
  GDISP[[names(GDISP)[tolower(names(GDISP)) == tolower(g)]]] else g

color_for <- function(gen, sexo, tto) {
  s <- toupper(sexo)
  if (toupper(tto) == "CONTROL") COL_CTRL[[s]] else COL_LPS[[gpath(gen)]][[s]]
}

# --- Estilo de trazo (un solo lugar, para que las 4 figuras coincidan) -------
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
GLAB    <- c("\u2640 C", "\u2640 LPS", "\u2642 C", "\u2642 LPS")

# =============================================================================
# 2. Anotacion de significancia (D11 ampliada).
#
# Cascada, en este orden -- solo se entra a UNA rama:
#   (a) interaccion SEXO x TTO significativa -> brackets por par del post hoc
#       D6 (ART-C con Holm), como estaba.
#   (b) interaccion NO significativa y efecto principal de TTO significativo ->
#       UN bracket que abarca los cuatro grupos (1 a 4, texto en 2.5). El
#       modelo dice que el efecto no difiere entre sexos, asi que marcar los
#       pares por separado sugeriria un dimorfismo que no esta sostenido.
#   (c) interaccion NO significativa y efecto principal de SEXO significativo ->
#       bracket de 1.5 a 3.5 (centro de cada sexo), texto en 2.5.
# Tendencia (0.05 <= p < 0.10): bracket punteado y el valor de p, sin estrella.
# =============================================================================
estrellas <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Dibuja un bracket. x1/x2 en coordenadas de caja; y es la base; solido si sig.
bracket <- function(x1, x2, y, alto, etiqueta, solido) {
  col <- if (solido) "black" else "grey40"
  lty <- if (solido) 1 else 2
  lines(c(x1, x1, x2, x2), c(y, y + alto, y + alto, y), col = col, lty = lty,
        lwd = 0.9, xpd = NA)
  text((x1 + x2) / 2, y + alto, etiqueta, pos = 3, col = col, cex = 0.85,
       xpd = NA)
}

# Devuelve la lista de brackets a dibujar para un gen x tejido.
# `fila`: la fila de clasificacion. `ph`: data.frame de post hoc del gen.
brackets_para <- function(fila, ph) {
  int_sig <- identical(as.character(fila[[COL_INT_SIG]]), "TRUE")
  out <- list()

  if (int_sig && nrow(ph)) {
    # (a) post hoc D6 -- mismo mapeo contraste -> par de cajas que en 07.
    PARES <- list("HEMBRA_CONTROL-HEMBRA_LPS"    = c(1, 2),
                  "MACHO_CONTROL-MACHO_LPS"      = c(3, 4),
                  "HEMBRA_LPS-MACHO_LPS"         = c(2, 4),
                  "HEMBRA_CONTROL-MACHO_CONTROL" = c(1, 3))
    for (i in order(ph$p_holm)) {
      p <- ph$p_holm[i]
      if (is.na(p) || p >= ALPHA_TREND) next
      pr <- PARES[[as.character(ph$contraste[i])]]
      if (is.null(pr)) next
      sig <- p < ALPHA
      out[[length(out) + 1L]] <- list(
        x1 = pr[1], x2 = pr[2], solido = sig,
        etiqueta = if (sig) estrellas(p) else sprintf("p = %.3f", p))
    }
    return(out)
  }

  # (b) efecto principal de tratamiento -- bracket centrado sobre los 4 grupos.
  p_tto <- suppressWarnings(as.numeric(fila[[COL_P_TTO]]))
  if (!is.na(p_tto) && p_tto < ALPHA_TREND) {
    sig <- p_tto < ALPHA
    out[[length(out) + 1L]] <- list(
      x1 = 1, x2 = 4, solido = sig,
      etiqueta = if (sig) sprintf("Control vs LPS %s", estrellas(p_tto))
                 else sprintf("Control vs LPS  p = %.3f", p_tto))
  }

  # (c) efecto principal de sexo -- entre los centros de cada sexo.
  p_sex <- suppressWarnings(as.numeric(fila[[COL_P_SEXO]]))
  if (!is.na(p_sex) && p_sex < ALPHA_TREND) {
    sig <- p_sex < ALPHA
    out[[length(out) + 1L]] <- list(
      x1 = 1.5, x2 = 3.5, solido = sig,
      etiqueta = if (sig) sprintf("\u2640 vs \u2642 %s", estrellas(p_sex))
                 else sprintf("\u2640 vs \u2642  p = %.3f", p_sex))
  }
  out
}

# =============================================================================
# 3. Panel de un gen x tejido.
#    Eje Y SIEMPRE logaritmico (D2): no se decide segun el rango de los datos.
# =============================================================================
panel_gen <- function(datos_gen, gen, tejido, fila_clasif, ph, con_titulo = TRUE) {
  dat <- lapply(GRUPOS, function(g)
    datos_gen$FC[datos_gen$SEXO == g[1] & datos_gen$TTO == g[2]])
  todos <- unlist(dat); todos <- todos[is.finite(todos) & todos > 0]
  if (!length(todos)) { plot.new(); return(invisible(NULL)) }

  brk <- brackets_para(fila_clasif, ph)
  n_niveles <- max(length(brk), 1L)

  # Espacio arriba para los brackets apilados (en escala log: factor por nivel).
  ymin <- min(todos) * 0.6
  ymax <- max(todos)
  ylim <- c(ymin, ymax * (1.30 ^ (n_niveles + 1)))

  rellenos <- vapply(GRUPOS, function(g)
    grDevices::adjustcolor(color_for(gen, g[1], g[2]), EST$alfa_relleno),
    character(1))

  par(mar = c(3.6, 4.4, if (con_titulo) 3.6 else 1.2, 0.8))
  boxplot(dat, xaxt = "n", outline = FALSE, log = "y", ylim = ylim,
          ylab = expression(FC == 2^{-Delta*Delta*Ct}),
          boxwex   = EST$ancho_caja,
          col      = rellenos,
          border   = EST$borde_caja,   # default para todos los trazos
          boxlwd   = EST$lwd_caja,
          whisklty = EST$lty_bigote,   # bigote punteado...
          whiskcol = EST$col_bigote,   # ...y mas claro que el borde
          whisklwd = EST$lwd_bigote,
          staplelty = 1,               # tope del bigote: solido
          staplecol = EST$col_tope,
          staplelwd = EST$lwd_tope,
          medcol   = EST$col_mediana,
          medlwd   = EST$lwd_mediana)

  axis(1, at = 1:4, labels = GLAB, tick = FALSE, line = -0.4)
  abline(h = 1, lty = 2, col = "grey70", lwd = 0.8)

  # Puntos individuales: circulo relleno del color del grupo, borde fino.
  for (j in seq_along(GRUPOS)) {
    g <- GRUPOS[[j]]; v <- dat[[j]]
    if (length(v))
      points(j + (runif(length(v)) - 0.5) * 0.26, v,
             pch = 21, bg = color_for(gen, g[1], g[2]),
             col = EST$borde_punto, lwd = EST$lwd_punto, cex = EST$cex_punto)
    mtext(sprintf("n=%d", length(v)), side = 1, line = 1.5, at = j, cex = 0.62)
  }

  # Brackets apilados de abajo hacia arriba.
  for (k in seq_along(brk)) {
    b <- brk[[k]]
    y <- ymax * (1.30 ^ k)
    bracket(b$x1, b$x2, y, y * 0.04, b$etiqueta, b$solido)
  }

  if (con_titulo) {
    p_int <- suppressWarnings(as.numeric(fila_clasif[[COL_P_INT]]))
    metodo <- as.character(fila_clasif[[COL_METODO]])
    title(main = sprintf("%s \u00b7 %s", gdisp(gen), tejido),
          cex.main = 1.0, font.main = 2, line = 2.2)
    title(main = sprintf("%s \u00b7 escala -\u0394\u0394Ct = log2(FC) \u00b7 p int = %s \u00b7 n = %d",
                         metodo,
                         if (is.na(p_int)) "NA" else sprintf("%.3f", p_int),
                         length(todos)),
          cex.main = 0.72, font.main = 1, line = 1.0)
  }
  invisible(NULL)
}

# =============================================================================
# 4. Panel de deteccion (il6 @ BRAIN_E15, D7).
#    Barras + simbolos individuales: relleno = detectado, vacio = no detectado,
#    con la misma convencion que las figuras de ELISA.
# =============================================================================
panel_deteccion <- function(tab2x4, fisher, gen = "il6", tejido = "BRAIN_E15") {
  nd <- numeric(4); nt <- numeric(4)
  for (j in seq_along(GRUPOS)) {
    g <- GRUPOS[[j]]
    r <- tab2x4[tab2x4$TEJIDO == tejido & tab2x4$GEN == gen &
                tab2x4$SEXO == g[1] & tab2x4$TTO == g[2], ][1, ]
    nd[j] <- as.numeric(r$n_detectado); nt[j] <- as.numeric(r$n_total)
  }
  pct <- ifelse(nt > 0, 100 * nd / nt, 0)

  par(mar = c(3.6, 4.4, 3.6, 0.8))
  bp <- barplot(pct, ylim = c(0, 125), names.arg = GLAB,
                ylab = "% detectado", border = EST$borde_caja,
                col = vapply(GRUPOS, function(g)
                  grDevices::adjustcolor(color_for(gen, g[1], g[2]),
                                         EST$alfa_relleno), character(1)),
                width = 0.62, space = 0.6, las = 1)

  # Un simbolo por feto: relleno si detectado, vacio si no.
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
    rr <- fisher[fisher$TEJIDO == tejido & fisher$GEN == gen & fisher$SEXO == s, ]
    if (nrow(rr)) as.numeric(rr$p_fisher[1]) else NA_real_
  }, numeric(1))

  title(main = sprintf("%s \u00b7 %s", gdisp(gen), tejido),
        cex.main = 1.0, font.main = 2, line = 2.2)
  title(main = sprintf("no cuantificable (D7) \u00b7 Fisher \u2640 p = %.3f \u00b7 \u2642 p = %.3f",
                       pf[["HEMBRA"]], pf[["MACHO"]]),
        cex.main = 0.72, font.main = 1, line = 1.0)
  invisible(NULL)
}

# =============================================================================
# 5. Figura global de un tejido: filas por naturaleza del gen.
#    Fila 1 lipidos (3) | fila 2 glucosa (2) | fila 3 aminoacidos (2) |
#    fila 4 via IL-6 (hasta 3). Las filas de 2 se centran dejando en blanco la
#    tercera columna (layout admite 0 = celda vacia).
# =============================================================================
figura_tejido <- function(datos, clasif, posthoc, tab2x4, fisher, tejido, ruta) {
  # il6R sale del panel de cerebro (deteccion insuficiente); en placenta queda.
  filas <- FILAS_VIA
  if (tejido == "BRAIN_E15")
    filas$IL6 <- setdiff(filas$IL6, "il6R")

  m <- matrix(0L, nrow = length(filas), ncol = 3)
  k <- 0L
  for (i in seq_along(filas)) {
    gs <- filas[[i]]
    despl <- if (length(gs) == 2L) 1L else 0L   # centra las filas de 2
    for (j in seq_along(gs)) { k <- k + 1L; m[i, j + despl] <- k }
  }

  png(ruta, width = 3 * 1150, height = length(filas) * 1000, res = 300)
  on.exit(dev.off(), add = TRUE)
  set.seed(SEMILLA)   # jitter reproducible
  layout(m)
  par(oma = c(0, 0, 3, 0))

  for (i in seq_along(filas)) for (gen in filas[[i]]) {
    if (tejido == "BRAIN_E15" && gen == "il6") {
      panel_deteccion(tab2x4, fisher); next
    }
    d <- datos[datos$TEJIDO == tejido & datos$GEN == gen & !is.na(datos$FC), ]
    fc <- clasif[clasif$TEJIDO == tejido & clasif$GEN == gen, ][1, ]
    ph <- posthoc[posthoc$TEJIDO == tejido & posthoc$GEN == gen, ]
    panel_gen(d, gen, tejido, fc, ph)
  }

  mtext(sprintf("Expresion relativa por gen -- %s",
                if (tejido == "PLACENTA_E15") "Placenta E15" else "Cerebro fetal E15"),
        outer = TRUE, cex = 1.15, font = 2, line = 0.6)
  invisible(ruta)
}
