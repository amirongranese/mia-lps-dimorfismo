# 03_elisa.R -- Analisis del ELISA de IL-6 con censura a izquierda (D10). Acto 1.1.
#
# Por que existe este archivo: valida el modelo MIA mostrando el efecto de LPS
# sobre IL-6 en suero materno (MS) y liquido amniotico (LA). El ELISA tiene
# censura a izquierda (Conc < 0 == "< LOD"); D10 fija el tratamiento:
#
#   1. Reportar % de deteccion / censura por grupo ANTES de cualquier estadistico.
#   2. Contraste primario = test de rango con los censurados empatados en el rango
#      mas bajo. Se usa Peto-Peto (Fleming-Harrington G-rho = 1) via
#      survival::survdiff(rho=1); ademas se implementa la MISMA formula a mano y
#      se verifica en corrida que coincide con survdiff (stopifnot, tol 1e-8), de
#      modo que R y Python (que implementa la formula a mano) queden alineados.
#   3. Si un grupo tiene censura tan alta que ningun estimador de ubicacion es
#      defendible -> solo proporcion de deteccion. Es el caso de MS: Control tiene
#      4/5 censurados (80%). Decision confirmada con el usuario: MS se analiza
#      SOLO como deteccion + Fisher exacto (Control vs LPS), sin test de ubicacion.
#   4. LA se analiza ESTRATIFICADO POR SEXO: dentro de cada sexo, Peto-Peto
#      Control vs LPS sobre Conc, mas Fisher de deteccion. Sin modelo factorial
#      (censura 0-60% + n 5-9 por celda lo hacen inestable) -- ver
#      analisis_descartados.md.
#
# Descriptivo con censura: mediana del estimador producto-limite (Kaplan-Meier)
# sobre el dato reflejado (t' = M - t pasa la censura izquierda a derecha). NO se
# usa ROS/NADA. LOD primario = 0 (el blanco), como fijo 02_ingesta_qc.
#
# Este script NO cuantifica qPCR y NO toca ninguna decision fuera de D10.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages(library(survival))

ESTE_SCRIPT <- "03_elisa"

# ---------------------------------------------------------------------------
# Grupos y ordenes canonicos.
# ---------------------------------------------------------------------------
GRUPOS_MS <- list(c("MS", "CONTROL", ""), c("MS", "LPS", ""))
GRUPOS_LA <- list(c("LA", "CONTROL", "HEMBRA"), c("LA", "CONTROL", "MACHO"),
                  c("LA", "LPS", "HEMBRA"), c("LA", "LPS", "MACHO"))
GRUPOS_TODOS <- c(GRUPOS_MS, GRUPOS_LA)
SEXOS_LA <- c("HEMBRA", "MACHO")
LOD_ELISA <- 0

COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")  # Okabe-Ito

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02 (paridad R/Python).
# Los p-valores / estadisticos dependientes de trascendentes se guardan como
# texto "%.6e" ya formateado (p6e): asi R y Python escriben el mismo string
# aunque su libm difiera en el ultimo bit.
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
s3e <- function(x) if (is.null(x) || length(x) != 1L || is.na(x)) "NA" else
  sprintf("%.3e", as.numeric(x))
g10 <- function(x) sprintf("%.10g", as.numeric(x))

# ---------------------------------------------------------------------------
# Lectura del intermedio de 02_ingesta_qc (mismo parseo que Python).
# ---------------------------------------------------------------------------
cargar_elisa <- function() {
  ruta <- file.path(RUTA_DATOS_PROC, "elisa_long.tsv")
  if (!file.exists(ruta))
    stop(sprintf("03_elisa: falta %s (correr 02_ingesta_qc primero)", ruta))
  raw <- read.delim(ruta, sep = "\t", quote = "", stringsAsFactors = FALSE,
                    colClasses = "character", check.names = FALSE,
                    encoding = "UTF-8", na.strings = character(0))
  data.frame(
    bloque = raw$bloque, TTO = raw$TTO,
    SEXO = ifelse(is.na(raw$SEXO), "", raw$SEXO),
    ID = raw$ID, Conc = as.numeric(raw$Conc),
    censurado = raw$censurado == "TRUE",
    LOD = ifelse(is.na(raw$LOD) | raw$LOD == "", 0, as.numeric(raw$LOD)),
    stringsAsFactors = FALSE
  )
}

sub_g <- function(df, bloque, tto, sexo) {
  df[df$bloque == bloque & df$TTO == tto &
       (sexo == "" | df$SEXO == sexo), , drop = FALSE]
}

mediana <- function(xs) {
  if (length(xs) == 0L) return(NA_real_)
  s <- sort(xs); k <- length(s)
  if (k %% 2L == 1L) s[(k + 1L) %/% 2L] else (s[k %/% 2L] + s[k %/% 2L + 1L]) / 2
}

# ===========================================================================
# 1. Deteccion / censura por grupo (D10, paso 1)
# ===========================================================================
tabla_deteccion <- function(df) {
  header <- c("bloque", "TTO", "SEXO", "n", "n_detectado", "n_censurado",
              "pct_detectado", "pct_censurado")
  filas <- lapply(GRUPOS_TODOS, function(k) {
    g <- sub_g(df, k[1], k[2], k[3]); n <- nrow(g)
    nc <- sum(g$censurado); nd <- n - nc
    list(k[1], k[2], k[3], n, nd, nc,
         if (n) 100 * nd / n else NA_real_,
         if (n) 100 * nc / n else NA_real_)
  })
  list(header = header, filas = filas)
}

# ===========================================================================
# 2. Descriptivo por grupo: detectados + mediana KM (dato reflejado)
# ===========================================================================
constante_reflexion <- function(df) ceiling(max(df$Conc[!df$censurado])) + 1

km_mediana_reflejada <- function(g, M) {
  t <- ifelse(g$censurado, M - g$LOD, M - g$Conc)
  e <- ifelse(g$censurado, 0L, 1L)
  tev <- sort(unique(t[e == 1L]))
  S <- 1; med_refl <- NA_real_
  for (tt in tev) {
    n_risk <- sum(t >= tt - 1e-12)
    d <- sum(abs(t - tt) <= 1e-12 & e == 1L)
    if (n_risk > 0) S <- S * (n_risk - d) / n_risk
    if (is.na(med_refl) && S <= 0.5 + 1e-12) med_refl <- tt
  }
  if (is.na(med_refl)) NA_real_ else M - med_refl
}

tabla_descriptivo <- function(df, M) {
  header <- c("bloque", "TTO", "SEXO", "n", "n_detectado", "pct_censurado",
              "min_detectado", "mediana_detectada", "max_detectado",
              "km_mediana", "km_nota")
  filas <- lapply(GRUPOS_TODOS, function(k) {
    g <- sub_g(df, k[1], k[2], k[3]); n <- nrow(g)
    dets <- sort(g$Conc[!g$censurado]); nd <- length(dets)
    kmm <- km_mediana_reflejada(g, M)
    nota <- if (is.na(kmm)) "S(t) no baja de 0.5 (censura alta)" else ""
    list(k[1], k[2], k[3], n, nd, if (n) 100 * (n - nd) / n else NA_real_,
         if (nd) dets[1] else NA_real_,
         if (nd) mediana(dets) else NA_real_,
         if (nd) dets[nd] else NA_real_,
         kmm, nota)
  })
  list(header = header, filas = filas)
}

# ===========================================================================
# 3. Peto-Peto (Fleming-Harrington G-rho = 1) para dos grupos.
#
# PROPIO -- misma formula que survival::survdiff (verificado abajo con stopifnot).
# ===========================================================================
peto_peto_2grupos <- function(tiempos, evento, grupo) {
  o <- order(tiempos)
  tiempos <- tiempos[o]; evento <- evento[o]; grupo <- grupo[o]
  gs <- sort(unique(grupo))
  if (length(gs) != 2L) stop("peto_peto_2grupos: se esperan 2 grupos")
  g1 <- gs[1]
  tev <- sort(unique(tiempos[evento == 1L]))
  O1 <- 0; E1 <- 0; V <- 0; S_prev <- 1
  for (t in tev) {
    en_riesgo <- tiempos >= t - 1e-12
    n <- sum(en_riesgo)
    n1 <- sum(en_riesgo & grupo == g1)
    ev_t <- abs(tiempos - t) <= 1e-12 & evento == 1L
    d <- sum(ev_t)
    d1 <- sum(ev_t & grupo == g1)
    w <- S_prev
    O1 <- O1 + w * d1
    E1 <- E1 + w * d * n1 / n
    if (n > 1) V <- V + w * w * d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1)
    S_prev <- S_prev * (n - d) / n
  }
  chisq <- if (V > 0) (O1 - E1)^2 / V else 0
  list(chisq = chisq, df = 1L, p = if (chisq <= 0) 1 else
    pchisq(chisq, 1, lower.tail = FALSE), O1 = O1, E1 = E1, V = V)
}

datos_supervivencia <- function(g, M) {
  list(t = ifelse(g$censurado, M - g$LOD, M - g$Conc),
       e = ifelse(g$censurado, 0L, 1L))
}

tabla_petopeto_la <- function(df, M) {
  header <- c("estrato", "n_control", "n_lps", "cens_control", "cens_lps",
              "chisq", "df", "p_valor", "metodo")
  filas <- lapply(SEXOS_LA, function(sx) {
    gc <- sub_g(df, "LA", "CONTROL", sx); gl <- sub_g(df, "LA", "LPS", sx)
    dc <- datos_supervivencia(gc, M); dl <- datos_supervivencia(gl, M)
    t <- c(dc$t, dl$t); e <- c(dc$e, dl$e)
    grp <- c(rep("CONTROL", nrow(gc)), rep("LPS", nrow(gl)))
    r <- peto_peto_2grupos(t, e, grp)
    # verificacion en corrida contra survival::survdiff(rho=1)
    sd <- survdiff(Surv(t, e) ~ grp, rho = 1)
    p_sd <- pchisq(sd$chisq, length(sd$n) - 1L, lower.tail = FALSE)
    stopifnot(abs(r$chisq - sd$chisq) < 1e-8, abs(r$p - p_sd) < 1e-8)
    list(sx, nrow(gc), nrow(gl), sum(dc$e == 0L), sum(dl$e == 0L),
         g10(r$chisq), 1L, p6e(r$p),
         "Peto-Peto G-rho=1 (survdiff) sobre Conc reflejada")
  })
  list(header = header, filas = filas)
}

# ===========================================================================
# 4. Fisher exacto a dos colas para 2x2 -- misma regla que stats::fisher.test.
#    PROPIO (verificado en corrida contra fisher.test).
# ===========================================================================
fisher_2x2 <- function(a, b, c, d) {
  r1 <- a + b; r2 <- c + d; c1 <- a + c; n <- a + b + c + d
  lp <- function(k)
    lgamma(r1 + 1) + lgamma(r2 + 1) + lgamma(c1 + 1) + lgamma(n - c1 + 1) -
    lgamma(n + 1) - lgamma(k + 1) - lgamma(r1 - k + 1) - lgamma(c1 - k + 1) -
    lgamma(r2 - c1 + k + 1)
  lo <- max(0, c1 - r2); hi <- min(r1, c1)
  p_obs <- exp(lp(a)); tot <- 0
  for (k in lo:hi) { pk <- exp(lp(k)); if (pk <= p_obs * (1 + 1e-7)) tot <- tot + pk }
  p <- min(1, tot)
  or_h <- ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
  # verificacion en corrida
  p_ref <- fisher.test(matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE))$p.value
  stopifnot(abs(p - p_ref) < 1e-9)
  list(or_h = or_h, p = p)
}

fila_fisher_deteccion <- function(df, bloque, sexo) {
  gc <- sub_g(df, bloque, "CONTROL", sexo); gl <- sub_g(df, bloque, "LPS", sexo)
  a <- sum(!gc$censurado); b <- nrow(gc) - a
  cc <- sum(!gl$censurado); dd <- nrow(gl) - cc
  r <- fisher_2x2(a, b, cc, dd)
  list(bloque, if (nzchar(sexo)) sexo else "(sin estrato)",
       a, nrow(gc), cc, nrow(gl), g10(r$or_h), p6e(r$p))
}

tabla_fisher <- function(df) {
  header <- c("bloque", "estrato", "control_detectado", "control_n",
              "lps_detectado", "lps_n", "or_haldane", "p_valor")
  filas <- c(list(fila_fisher_deteccion(df, "MS", "")),
             lapply(SEXOS_LA, function(sx) fila_fisher_deteccion(df, "LA", sx)))
  list(header = header, filas = filas)
}

# ===========================================================================
# 5. Figuras del Acto 1.1 (ggplot2; equivalentes a las de matplotlib en Python)
# ===========================================================================
suppressMessages(library(ggplot2))

.offsets <- function(k) {
  if (k <= 0) return(numeric(0))
  if (k == 1) return(0)
  seq(-0.15, 0.15, length.out = k)
}

.df_puntos <- function(df, bloque, sexos) {
  filas <- list()
  for (sx in sexos) for (i in seq_along(c("CONTROL", "LPS"))) {
    tto <- c("CONTROL", "LPS")[i]
    g <- sub_g(df, bloque, tto, if (bloque == "MS") "" else sx)
    dets <- sort(g$Conc[!g$censurado]); ncen <- sum(g$censurado)
    if (length(dets))
      filas[[length(filas) + 1L]] <- data.frame(
        sexo = sx, tto = tto, x = (i - 1) + .offsets(length(dets)),
        y = dets, tipo = "detectado", stringsAsFactors = FALSE)
    if (ncen)
      filas[[length(filas) + 1L]] <- data.frame(
        sexo = sx, tto = tto, x = (i - 1) + .offsets(ncen),
        y = 0, tipo = "censurado", stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, filas)
  out$tipo <- factor(out$tipo, levels = c("detectado", "censurado"))
  out
}

.df_boxplot <- function(df, bloque, sexos) {
  filas <- list()
  for (sx in sexos) for (i in seq_along(c("CONTROL", "LPS"))) {
    tto <- c("CONTROL", "LPS")[i]
    g <- sub_g(df, bloque, tto, if (bloque == "MS") "" else sx)
    dets <- g$Conc[!g$censurado]
    if (length(dets) >= 3)
      filas[[length(filas) + 1L]] <- data.frame(
        sexo = sx, tto = tto, xc = i - 1, y = dets, stringsAsFactors = FALSE)
  }
  if (length(filas)) do.call(rbind, filas) else NULL
}

.df_anota <- function(df, bloque, sexos, ymax) {
  filas <- list()
  for (sx in sexos) for (i in seq_along(c("CONTROL", "LPS"))) {
    tto <- c("CONTROL", "LPS")[i]
    g <- sub_g(df, bloque, tto, if (bloque == "MS") "" else sx)
    n <- nrow(g); nd <- sum(!g$censurado)
    filas[[length(filas) + 1L]] <- data.frame(
      sexo = sx, tto = tto, x = i - 1, y = ymax * 1.05,
      lab = sprintf("%d/%d det.\n%.0f%%", nd, n, 100 * nd / n),
      stringsAsFactors = FALSE)
  }
  do.call(rbind, filas)
}

.tema_elisa <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "bottom", plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9.5))
}

figura_ms <- function(df, fisher_ms, ruta) {
  ymax <- max(df$Conc[df$bloque == "MS" & !df$censurado])
  pts <- .df_puntos(df, "MS", "MS"); box <- .df_boxplot(df, "MS", "MS")
  ano <- .df_anota(df, "MS", "MS", ymax)
  p <- ggplot()
  if (!is.null(box))
    p <- p + geom_boxplot(data = box, aes(x = xc, y = y, group = xc,
                                          colour = tto, fill = tto),
                          width = 0.46, outlier.shape = NA, alpha = 0.14,
                          show.legend = FALSE)
  p <- p +
    geom_hline(yintercept = LOD_ELISA, linetype = "dashed", colour = "grey60",
               linewidth = 0.3) +
    geom_point(data = pts, aes(x, y, colour = tto, shape = tipo),
               size = 2.5, stroke = 1) +
    geom_text(data = ano, aes(x, y, label = lab), size = 3, colour = "grey25",
              lineheight = 0.9, vjust = 0) +
    scale_colour_manual(values = COL_TTO, guide = "none") +
    scale_fill_manual(values = COL_TTO, guide = "none") +
    scale_shape_manual(values = c(detectado = 19, censurado = 1),
                       name = NULL,
                       labels = c(detectado = "detectado",
                                  censurado = "censurado (< LOD)")) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Control", "LPS"),
                       limits = c(-0.6, 1.6)) +
    ylim(-ymax * 0.06, ymax * 1.22) +
    labs(title = "IL-6 en suero materno (ELISA) -- validacion del modelo MIA",
         subtitle = sprintf("Fisher exacto (deteccion) Control vs LPS:  p = %s",
                            s3e(as.numeric(fisher_ms[[8]]))),
         x = NULL, y = "IL-6  (pg/mL)") +
    .tema_elisa()
  ggsave(ruta, p, width = 4.8, height = 5.0, dpi = 300)
}

figura_la <- function(df, peto_filas, ruta) {
  ymax <- max(df$Conc[df$bloque == "LA" & !df$censurado])
  pts <- .df_puntos(df, "LA", SEXOS_LA); box <- .df_boxplot(df, "LA", SEXOS_LA)
  ano <- .df_anota(df, "LA", SEXOS_LA, ymax)
  pmap <- setNames(vapply(peto_filas, function(r) r[[8]], character(1)), SEXOS_LA)
  etiquetas <- setNames(
    sprintf("%s   (Peto-Peto p = %s)",
            c(HEMBRA = "Hembra", MACHO = "Macho")[SEXOS_LA],
            vapply(SEXOS_LA, function(sx) s3e(as.numeric(pmap[[sx]])), character(1))),
    SEXOS_LA)
  p <- ggplot()
  if (!is.null(box))
    p <- p + geom_boxplot(data = box, aes(x = xc, y = y, group = xc,
                                          colour = tto, fill = tto),
                          width = 0.46, outlier.shape = NA, alpha = 0.14,
                          show.legend = FALSE)
  p <- p +
    geom_hline(yintercept = LOD_ELISA, linetype = "dashed", colour = "grey60",
               linewidth = 0.3) +
    geom_point(data = pts, aes(x, y, colour = tto, shape = tipo, fill = tto),
               size = 2.4, stroke = 0.9) +
    geom_text(data = ano, aes(x, y, label = lab), size = 3, colour = "grey25",
              lineheight = 0.9, vjust = 0) +
    facet_wrap(~ sexo, labeller = as_labeller(etiquetas)) +
    scale_colour_manual(values = COL_TTO, guide = "none") +
    scale_fill_manual(values = COL_TTO, guide = "none") +
    scale_shape_manual(values = c(detectado = 16, censurado = 21),
                       name = NULL,
                       labels = c(detectado = "detectado",
                                  censurado = "censurado (< LOD)")) +
    scale_x_continuous(breaks = c(0, 1), labels = c("Control", "LPS"),
                       limits = c(-0.6, 1.6)) +
    ylim(-ymax * 0.06, ymax * 1.22) +
    labs(title = "IL-6 en liquido amniotico (ELISA) por sexo fetal",
         x = NULL, y = "IL-6  (pg/mL)") +
    .tema_elisa()
  ggsave(ruta, p, width = 7.6, height = 5.0, dpi = 300)
}

# ===========================================================================
# 6. Reporte legible
# ===========================================================================
.md_tabla <- function(header, filas) {
  l1 <- paste0("| ", paste(header, collapse = " | "), " |")
  l2 <- paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |")
  cuerpo <- vapply(filas, function(f)
    paste0("| ", paste(vapply(f, .fmt, character(1)), collapse = " | "), " |"),
    character(1))
  paste(c(l1, l2, cuerpo), collapse = "\n")
}

construir_reporte <- function(fuente, det, desc, fish, peto, M) {
  ms <- fish$filas[[1]]
  L <- c(
    "# Reporte del ELISA de IL-6 (Acto 1.1)", "",
    "Generado por `03_elisa` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. Deteccion / censura por grupo (D10, antes de cualquier estadistico)",
    "", .md_tabla(det$header, det$filas), "",
    "## 2. Descriptivo por grupo", "",
    sprintf(paste0("`km_mediana` = mediana del estimador producto-limite ",
                   "(Kaplan-Meier) sobre el dato reflejado `t' = M - t` (M = %s = ",
                   "ceil(max Conc detectada) + 1). Vacio = S(t) no baja de 0.5 por ",
                   "censura alta -> solo proporcion de deteccion."), .fmt(M)),
    "", .md_tabla(desc$header, desc$filas), "",
    "## 3. Suero materno (MS): solo deteccion + Fisher exacto", "",
    paste0("Control tiene 80% de censura (1 valor detectado): ningun estimador de ",
           "ubicacion es defendible (D10, clausula final). Se reporta solo la ",
           "proporcion de deteccion y el Fisher exacto Control vs LPS."), "",
    sprintf("- Control: %s/%s detectados; LPS: %s/%s detectados.",
            .fmt(ms[[3]]), .fmt(ms[[4]]), .fmt(ms[[5]]), .fmt(ms[[6]])),
    sprintf("- Fisher exacto (2 colas): OR (Haldane) = %s, p = %s.",
            ms[[7]], s3e(as.numeric(ms[[8]]))), "",
    "## 4. Liquido amniotico (LA): Peto-Peto por sexo + Fisher de deteccion", "",
    paste0("Contraste primario: Peto-Peto (Fleming-Harrington G-rho = 1) Control ",
           "vs LPS sobre `Conc`, con los censurados empatados en el rango mas bajo. ",
           "Estratificado por sexo (sin modelo factorial: ver ",
           "`analisis_descartados.md`)."), "",
    .md_tabla(peto$header, peto$filas), "",
    "Fisher exacto sobre deteccion (Control vs LPS) dentro de cada sexo:", "",
    .md_tabla(fish$header, fish$filas), "",
    "## 5. Figuras", "",
    paste0("- `outputs/figures/acto1_elisa_ms.png` -- IL-6 en suero materno por ",
           "tratamiento; censurados dibujados en el LOD (simbolo abierto)."),
    paste0("- `outputs/figures/acto1_elisa_la.png` -- IL-6 en liquido amniotico ",
           "por sexo y tratamiento."), "",
    "## 6. Notas metodologicas", "",
    paste0("Ver `analisis_descartados.md`, seccion `03_elisa`: por que no se uso ",
           "ROS, por que MS no lleva test de ubicacion, por que LA no lleva modelo ",
           "factorial, la reflexion para aplicar survdiff a censura a izquierda, y ",
           "el LOD alternativo (menor estandar de `CURVA IL6`) no aplicado."), ""
  )
  paste(L, collapse = "\n")
}

# ===========================================================================
# 7. Artefactos compartidos (merge por 'script') -- headers identicos a 02.
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

BLOQUE_DESCARTES <- paste(c(
"## 03_elisa",
"",
"### ROS / KM-NADA no usado como estimador de ubicacion",
"",
"- **Que se probo:** reportar media/mediana por grupo con ROS (regression on order",
"  statistics; `NADA` en R) como estimacion puntual con censura.",
"- **Por que no:** (a) agrega una dependencia (`NADA`) fuera del toolchain minimo;",
"  (b) con n 5-9 por grupo y censura 30-80% el ajuste log-normal implicito de ROS",
"  es fragil y sus IC no son defendibles; (c) D10 admite explicitamente el test de",
"  rango con censurados al fondo como alternativa.",
"- **Que se hizo:** contraste primario = Peto-Peto (G-rho = 1). Descriptivo =",
"  mediana del estimador producto-limite (Kaplan-Meier) sobre el dato reflejado;",
"  donde S(t) no baja de 0.5 (censura alta) se reporta vacio + proporcion de",
"  deteccion.",
"",
"### Suero materno (MS): sin test de ubicacion",
"",
"- **Situacion:** Control tiene 4/5 censurados (80%): un solo valor detectado.",
"- **Por que:** ningun estimador de ubicacion (KM, ROS, media) es defendible con",
"  1 dato. D10 (clausula final) -> reportar solo proporcion de deteccion.",
"- **Que se hizo:** MS se analiza solo como deteccion + Fisher exacto Control vs",
"  LPS. La magnitud del efecto en MS se describe (mediana de detectados) pero no",
"  se testea. Decision confirmada con el usuario.",
"",
"### Liquido amniotico (LA): sin modelo factorial SEXO x TTO",
"",
"- **Que se probo:** un modelo factorial con censura (survreg / Tobit) sobre",
"  SEXO * TTO.",
"- **Por que no:** censura 0-60% por celda + n 5-9 lo dejan mal condicionado.",
"- **Que se hizo:** analisis estratificado por sexo -- Peto-Peto Control vs LPS",
"  dentro de hembra y dentro de macho, mas Fisher de deteccion por sexo. La",
"  comparacion formal hembra vs macho del efecto LPS queda fuera de T3.",
"",
"### Peto-Peto por reflexion (censura a izquierda)",
"",
"- `survival::survdiff` maneja censura a DERECHA; el ELISA tiene censura a",
"  IZQUIERDA. Se refleja `t' = M - t` con `M = ceil(max Conc detectada) + 1`: los",
"  \"< LOD\" pasan a censura a derecha en `M - LOD` y quedan empatados en el rango",
"  mas bajo de la escala original, que es lo que pide D10.",
"- La implementacion Python calcula la formula G-rho a mano; la implementacion R",
"  hace lo mismo y ademas se verifica en corrida contra `survival::survdiff(rho=1)`",
"  (`stopifnot`, tolerancia 1e-8).",
"",
"### LOD alternativo no aplicado",
"",
"- D10 menciona como alternativa el menor estandar de la hoja `CURVA IL6`. Aca se",
"  mantiene `LOD = 0` (el blanco de la placa, como en 02_ingesta_qc). El LOD",
"  alternativo queda como analisis de sensibilidad para T9/T10, no aplicado.",
"",
"### p-valores como texto en las tablas",
"",
"- Los p-valores y el `chisq` que dependen de funciones trascendentes (`lgamma`,",
"  `exp`, `erfc` / `pchisq`) se guardan formateados `\"%.6e\"`: asi R y Python",
"  escriben el mismo texto aunque sus librerias matematicas difieran en el ultimo",
"  bit. La comparacion fina R<->Python se hace en T10 con la tolerancia declarada."
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 03_elisa:inicio -->"
  marca_fin <- "<!-- 03_elisa:fin -->"
  nuevo <- paste0(marca_ini, "\n", BLOQUE_DESCARTES, "\n\n", marca_fin)
  if (file.exists(ruta)) {
    txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE)
    Encoding(txt) <- "UTF-8"
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
  df <- cargar_elisa()
  fuente <- fuente_datos(ARCHIVO_ELISA)
  M <- constante_reflexion(df)

  det  <- tabla_deteccion(df)
  desc <- tabla_descriptivo(df, M)
  fish <- tabla_fisher(df)
  peto <- tabla_petopeto_la(df, M)

  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "elisa_deteccion.csv"), det$header, det$filas)
    escribir_csv(file.path(base, "elisa_descriptivo.csv"), desc$header, desc$filas)
    escribir_csv(file.path(base, "elisa_fisher_deteccion.csv"), fish$header, fish$filas)
    escribir_csv(file.path(base, "elisa_petopeto_la.csv"), peto$header, peto$filas)
  }

  figura_ms(df, fish$filas[[1]], file.path(RUTA_FIGURAS, "acto1_elisa_ms.png"))
  figura_la(df, peto$filas, file.path(RUTA_FIGURAS, "acto1_elisa_la.png"))

  escribir_lineas(file.path(RUTA_TABLAS, "elisa_reporte.md"),
                  construir_reporte(fuente, det, desc, fish, peto, M))
  actualizar_descartados()

  ent <- sprintf("data/processed/elisa_long.tsv (de data/%s/%s)", fuente, ARCHIVO_ELISA)
  registrar_procedencia(list(
    list("outputs/tables/{R,python}/elisa_deteccion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "deteccion/censura por grupo (D10 paso 1)"),
    list("outputs/tables/{R,python}/elisa_descriptivo.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "descriptivo por grupo: detectados + mediana KM (dato reflejado)"),
    list("outputs/tables/{R,python}/elisa_fisher_deteccion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "Fisher exacto sobre deteccion: MS y LA por sexo (Control vs LPS)"),
    list("outputs/tables/{R,python}/elisa_petopeto_la.csv", "tabla", ESTE_SCRIPT,
         "PROPIO (R verifica contra survival::survdiff rho=1)", ent,
         "LA por sexo: Peto-Peto Control vs LPS sobre Conc reflejada"),
    list("outputs/tables/elisa_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO", ent,
         "reporte legible del ELISA (Acto 1.1)"),
    list("outputs/figures/acto1_elisa_ms.png", "figura", ESTE_SCRIPT, "PROPIO", ent,
         "IL-6 suero materno: Conc por tratamiento; censurados en el LOD"),
    list("outputs/figures/acto1_elisa_la.png", "figura", ESTE_SCRIPT, "PROPIO", ent,
         "IL-6 liquido amniotico por sexo: Conc por tratamiento")
  ))

  ms_ok <- identical(fish$filas[[1]][[1]], "MS")
  fig_ms <- file.exists(file.path(RUTA_FIGURAS, "acto1_elisa_ms.png"))
  fig_la <- file.exists(file.path(RUTA_FIGURAS, "acto1_elisa_la.png"))
  registrar_verificaciones(list(
    list("elisa_ms_metodo", "MS: solo deteccion + Fisher exacto (D10, censura Control alta)",
         "fisher_deteccion", "fisher_deteccion", if (ms_ok) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("elisa_la_metodo", "LA: Peto-Peto G-rho=1 Control vs LPS, estratificado por sexo",
         "petopeto_por_sexo", "petopeto_por_sexo", "TRUE", ESTE_SCRIPT),
    list("elisa_la_estratos", "estratos de LA", paste(SEXOS_LA, collapse = ";"),
         "HEMBRA;MACHO", if (identical(SEXOS_LA, c("HEMBRA", "MACHO"))) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("elisa_descriptivo_grupos", "grupos en el descriptivo (MSx2 + LAx4)",
         as.character(length(desc$filas)), "6",
         if (length(desc$filas) == 6L) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("elisa_lod", "LOD primario del ELISA (el blanco de la placa)",
         .fmt(LOD_ELISA), "0", if (LOD_ELISA == 0) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("elisa_reflexion_M", "constante de reflexion M = ceil(max Conc detectada)+1",
         .fmt(M), ">0", if (M > 0) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("elisa_figura_ms", "figura Acto 1.1 suero materno existe",
         "acto1_elisa_ms.png", "existe", if (fig_ms) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("elisa_figura_la", "figura Acto 1.1 liquido amniotico existe",
         "acto1_elisa_la.png", "existe", if (fig_la) "TRUE" else "FALSE", ESTE_SCRIPT)
  ))

  cat("== 03_elisa.R ==\n")
  cat(sprintf("  fuente ELISA = %s   M(reflexion) = %s\n", fuente, .fmt(M)))
  for (f in det$filas)
    cat(sprintf("    %-2s %-7s %-6s  n=%2d  det=%2d/%-2d  cens=%d\n",
                f[[1]], f[[2]], f[[3]], f[[4]], f[[5]], f[[4]], f[[6]]))
  cat(sprintf("  MS Fisher deteccion Control vs LPS: OR(Haldane)=%s  p=%s\n",
              fish$filas[[1]][[7]], fish$filas[[1]][[8]]))
  for (f in peto$filas)
    cat(sprintf("  LA %-6s Peto-Peto: chisq=%s  p=%s\n", f[[1]], f[[6]], f[[8]]))
  cat("  -> outputs/tables/{R,python}/elisa_*.csv, outputs/figures/acto1_elisa_*.png,\n")
  cat("     outputs/tables/elisa_reporte.md\n")
}

if (sys.nframe() == 0L) main()
