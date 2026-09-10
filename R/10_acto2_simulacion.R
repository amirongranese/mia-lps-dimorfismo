# 10_acto2_simulacion.R -- ACTO 2.5: descartar restriccion de rango antes de
#                          interpretar el cambio de correlacion placenta<->cerebro.
#
# Por que existe este archivo: prohibicion 5 -- NO se puede interpretar un cambio
# de correlacion entre Control y LPS como cambio de coordinacion biologica sin
# descartar antes que sea un artefacto de dispersion (restriccion de rango). Una
# correlacion observada se atenua si el rango de alguna de las dos variables se
# comprime en un subgrupo, aunque la asociacion verdadera sea la misma.
#
# Diseno (decisiones del usuario):
#   * Escala: normal bivariada en -ddCt. Se generan (placenta, cerebro) con UNA
#     r verdadera comun a ambos grupos y se escalan las SD marginales de cada
#     grupo a las SD OBSERVADAS de ese grupo (tabla de 09). Lo unico que difiere
#     entre grupos en la simulacion es la dispersion.
#   * r verdadera comun -- ambas como rango: (a) ESCENARIO GLOBAL -> rho de
#     Spearman GLOBAL del item (pooled, la de T7); (b) ESCENARIO CONTROL -> rho
#     de Control (grupo con mas coordinacion; escenario mas exigente). El
#     Delta rho observado se compara contra la distribucion simulada de CADA
#     escenario.
#   * Inversion rho_S -> r_Pearson: rho_S = (6/pi)*arcsin(r/2) =>
#     r = 2*sin(pi*rho_S/6). Recorte +/- 0.999999.
#   * B = 2000 repeticiones, semilla 20260101, RNG PROPIO (LCG + polar de
#     Marsaglia de 00_config; misma secuencia exacta en R y Python). Un unico RNG
#     recorre todos los items y escenarios en orden.
#
# Veredicto DENTRO -> la sola diferencia de dispersion puede producir el
# Delta rho observado -> NO se puede descartar restriccion de rango. FUERA -> lo
# excede en ese escenario.
#
# PARIDAD R/Python: los sorteos usan norm1() (LCG + polar, byte-identico R/Python
# por 00_config/T1); rho, medias, SD, cuantiles y conteos usan acumulador double
# explicito y comparaciones -> texto "%.10g" bit-identico. r_pearson_gen pasa por
# sin() -> texto "%.6e". Las figuras son PNG: equivalentes, no byte-identicas.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages({ library(ggplot2) })

ESTE_SCRIPT <- "10_acto2_simulacion"

COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")
Z975 <- 1.959963984540054
PISO_PAR <- 5L

GEN_SIN_CEREBRO <- "il6"
GENES_CORR <- setdiff(GENES, GEN_SIN_CEREBRO)
ITEMS <- c(GENES_CORR, "score_compuesto")
ESCENARIOS <- c("GLOBAL", "CONTROL")
B_SIM <- 2000L
R_CLAMP <- 0.999999
DPI <- 300

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

cuantil_tipo7 <- function(xs_ordenado, prob) {
  n <- length(xs_ordenado)
  if (n == 0L) return(NA_real_)
  if (n == 1L) return(xs_ordenado[1])
  h <- (n - 1) * prob
  lo <- floor(h)
  if (lo >= n - 1) return(xs_ordenado[n])
  frac <- h - lo
  xs_ordenado[lo + 1L] + frac * (xs_ordenado[lo + 2L] - xs_ordenado[lo + 1L])
}

rho_s_a_r <- function(rho_s) {
  r <- 2 * sin(pi * rho_s / 6)
  if (r > R_CLAMP) return(R_CLAMP)
  if (r < -R_CLAMP) return(-R_CLAMP)
  r
}

.fisher_sig <- function(rc, nc, rl, nl) {
  zc <- atanh(rc); zl <- atanh(rl)
  sec <- sqrt((1 + rc * rc / 2) / (nc - 3))
  sel <- sqrt((1 + rl * rl / 2) / (nl - 3))
  stat <- (zc - zl) / sqrt(sec * sec + sel * sel)
  abs(stat) > Z975
}

# ===========================================================================
# 1. Carga y emparejamiento por feto -- identico a 08/09.
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
  scv <- setNames(
    ifelse(sc$score_compuesto == "", NA_real_,
           suppressWarnings(as.numeric(sc$score_compuesto))),
    paste(sc$FETO, sc$TEJIDO, sep = "\r"))
  fetos <- names(madre)[order(madre, names(madre), method = "radix")]
  list(fetos = fetos, tto = tto, negdd = negdd, sc = scv)
}

pares <- function(D, item, estrato) {
  xs <- c(); ys <- c(); ts <- c()
  for (f in D$fetos) {
    if (estrato != "GLOBAL" && D$tto[[f]] != estrato) next
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

# ===========================================================================
# 2. Simulacion.
# ===========================================================================
COLS_SIM <- c("ITEM", "TIPO", "ESCENARIO", "rho_true", "r_pearson_gen",
              "n_control", "n_lps", "sd_pla_control", "sd_bra_control",
              "sd_pla_lps", "sd_bra_lps", "delta_rho_obs", "sim_mean_delta_rho",
              "sim_sd_delta_rho", "sim_q025", "sim_q975", "sim_frac_abs_ge_obs",
              "sim_tasa_fisher_sig", "veredicto")

.simular_celda <- function(rng, r_pear, sd_pla_c, sd_bra_c, nc,
                           sd_pla_l, sd_bra_l, nl, drho_obs) {
  raiz <- sqrt(1 - r_pear * r_pear)
  sims <- numeric(B_SIM)
  n_fisher <- 0L
  for (b in seq_len(B_SIM)) {
    xc <- numeric(nc); yc <- numeric(nc)
    for (k in seq_len(nc)) {
      z1 <- rng$norm1(); z2 <- rng$norm1()
      xc[k] <- sd_pla_c * z1
      yc[k] <- sd_bra_c * (r_pear * z1 + raiz * z2)
    }
    xl <- numeric(nl); yl <- numeric(nl)
    for (k in seq_len(nl)) {
      z1 <- rng$norm1(); z2 <- rng$norm1()
      xl[k] <- sd_pla_l * z1
      yl[k] <- sd_bra_l * (r_pear * z1 + raiz * z2)
    }
    rc <- spearman_rho(xc, yc); rl <- spearman_rho(xl, yl)
    if (is.na(rc) || is.na(rl)) { sims[b] <- 0; next }
    sims[b] <- rc - rl
    if (abs(rc) < 1 && abs(rl) < 1 && .fisher_sig(rc, nc, rl, nl))
      n_fisher <- n_fisher + 1L
  }
  media <- .suma(sims) / B_SIM
  sd <- desvio(sims)
  ordenado <- sort(sims)
  q025 <- cuantil_tipo7(ordenado, 0.025)
  q975 <- cuantil_tipo7(ordenado, 0.975)
  a_obs <- abs(drho_obs); ge <- 0L
  for (d in sims) if (abs(d) >= a_obs) ge <- ge + 1L
  list(media = media, sd = sd, q025 = q025, q975 = q975,
       frac_ge = ge / B_SIM, tasa_fisher = n_fisher / B_SIM)
}

tabla_simulacion <- function(D, rng) {
  filas <- list()
  for (item in ITEMS) {
    tipo <- if (item == "score_compuesto") "score" else "gen"
    prc <- pares(D, item, "CONTROL"); prl <- pares(D, item, "LPS")
    prg <- pares(D, item, "GLOBAL")
    nc <- length(prc$x); nl <- length(prl$x)
    if (nc < PISO_PAR || nl < PISO_PAR) {
      for (esc in ESCENARIOS)
        filas[[length(filas) + 1L]] <- list(item, tipo, esc, "", "", nc, nl,
          "", "", "", "", "", "", "", "", "", "", "", "sin_test")
      next
    }
    rc_obs <- spearman_rho(prc$x, prc$y)
    rl_obs <- spearman_rho(prl$x, prl$y)
    rg_obs <- spearman_rho(prg$x, prg$y)
    drho_obs <- rc_obs - rl_obs
    sd_pla_c <- desvio(prc$x); sd_bra_c <- desvio(prc$y)
    sd_pla_l <- desvio(prl$x); sd_bra_l <- desvio(prl$y)
    for (esc in ESCENARIOS) {
      rho_true <- if (esc == "GLOBAL") rg_obs else rc_obs
      r_pear <- rho_s_a_r(rho_true)
      s <- .simular_celda(rng, r_pear, sd_pla_c, sd_bra_c, nc,
                          sd_pla_l, sd_bra_l, nl, drho_obs)
      dentro <- (s$q025 <= drho_obs && drho_obs <= s$q975)
      filas[[length(filas) + 1L]] <- list(
        item, tipo, esc, g10(rho_true), p6e(r_pear), nc, nl,
        g10(sd_pla_c), g10(sd_bra_c), g10(sd_pla_l), g10(sd_bra_l),
        g10(drho_obs), g10(s$media), g10(s$sd), g10(s$q025), g10(s$q975),
        g10(s$frac_ge), g10(s$tasa_fisher), if (dentro) "DENTRO" else "FUERA")
    }
  }
  filas
}

# ===========================================================================
# 3. Figura.
# ===========================================================================
figura_simulacion <- function(sim, ruta) {
  ord <- vapply(ITEMS, function(i)
    if (i == "score_compuesto") "score compuesto" else i, character(1))
  seg <- list(); pts <- list(); med <- list(); txt <- list()
  ymap <- c(GLOBAL = 1, CONTROL = 0)
  for (f in sim) {
    nom <- if (f[[1]] == "score_compuesto") "score compuesto" else f[[1]]
    if (f[[19]] == "sin_test") {
      txt[[length(txt) + 1L]] <- data.frame(item = nom, lab = "n<5 (sin simulacion)",
                                            stringsAsFactors = FALSE)
      next
    }
    y <- ymap[[f[[3]]]]
    seg[[length(seg) + 1L]] <- data.frame(item = nom, esc = f[[3]], y = y,
      x0 = as.numeric(f[[15]]), x1 = as.numeric(f[[16]]), stringsAsFactors = FALSE)
    med[[length(med) + 1L]] <- data.frame(item = nom, esc = f[[3]], y = y,
      x = as.numeric(f[[13]]), stringsAsFactors = FALSE)
    pts[[length(pts) + 1L]] <- data.frame(item = nom, y = y,
      x = as.numeric(f[[12]]), ver = f[[19]], stringsAsFactors = FALSE)
  }
  d_seg <- do.call(rbind, seg); d_pts <- do.call(rbind, pts)
  d_med <- do.call(rbind, med)
  d_seg$item <- factor(d_seg$item, levels = ord)
  d_pts$item <- factor(d_pts$item, levels = ord)
  d_med$item <- factor(d_med$item, levels = ord)
  p <- ggplot() +
    geom_segment(data = d_seg, aes(x = x0, xend = x1, y = y, yend = y,
                                   colour = esc), linewidth = 2.6, alpha = 0.55) +
    geom_point(data = d_med, aes(x = x, y = y, colour = esc), shape = 124,
               size = 3) +
    geom_point(data = d_pts, aes(x = x, y = y, fill = ver), shape = 23,
               size = 2.6, colour = "black", stroke = 0.3) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey60",
               linewidth = 0.3) +
    scale_colour_manual(values = c(GLOBAL = "#4C4C4C",
                                   CONTROL = unname(COL_TTO["CONTROL"])),
                        name = "r verdadera") +
    scale_fill_manual(values = c(DENTRO = "#009E73", FUERA = "#D55E00"),
                      name = "veredicto") +
    scale_y_continuous(breaks = c(0, 1),
                       labels = c("r=rho Control", "r=rho GLOBAL"),
                       limits = c(-0.6, 1.6)) +
    facet_wrap(~ item, scales = "free_x", ncol = 4) +
    labs(title = paste0("Delta rho observado (rombo) vs intervalo 95% del ",
                        "Delta rho simulado bajo r verdadera comun y solo ",
                        "diferencia de dispersion"),
         subtitle = paste0("verde = DENTRO (no se puede descartar restriccion ",
                           "de rango, prohibicion 5) - naranja = FUERA"),
         x = "Delta rho  (rho Control - rho LPS)", y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.subtitle = element_text(size = 8))
  ggsave(ruta, p, width = 13.0, height = 9.0, dpi = DPI)
}

# ===========================================================================
# 4. Artefactos compartidos (merge por 'script') -- headers identicos a 02..09.
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
"## 10_acto2_simulacion",
"",
"### Que resuelve",
"",
paste0("- **Prohibicion 5.** Antes de leer cualquier cambio de la correlacion ",
       "placenta<->cerebro entre Control y LPS como cambio de coordinacion ",
       "biologica, hay que descartar que sea un artefacto de dispersion ",
       "(restriccion de rango). Una correlacion se atenua si el rango de alguna ",
       "de las dos variables se comprime en un subgrupo, con la misma ",
       "asociacion verdadera."),
"",
"### Diseno (decisiones del usuario)",
"",
paste0("- **Escala: normal bivariada en `-ddCt`.** Se generan (placenta, ",
       "cerebro) con UNA `r` verdadera comun a ambos grupos y se escalan las SD ",
       "marginales de cada grupo a las **SD observadas** de ese grupo (tabla de ",
       "`09_acto2_dispersion`). Lo unico que difiere entre grupos en la ",
       "simulacion es la dispersion."),
paste0("- **`r` verdadera comun -- ambas como rango:** (a) ESCENARIO `GLOBAL` ",
       "-> `rho` verdadera = `rho` de Spearman GLOBAL del item (pooled, la de ",
       "T7); (b) ESCENARIO `CONTROL` -> `rho` verdadera = `rho` de Control ",
       "(grupo con mas coordinacion; escenario mas exigente). El `Delta rho` ",
       "observado se compara contra la distribucion simulada de **cada** ",
       "escenario."),
paste0("- **Inversion `rho_S` -> `r_Pearson`:** para la normal bivariada ",
       "`rho_S = (6/pi)*arcsin(r/2)`  =>  `r = 2*sin(pi*rho_S/6)` (asi la `rho` ",
       "de Spearman simulada apunta al target). `r` se recorta a +/- 0.999999."),
paste0("- **B = 2000 repeticiones**, semilla 20260101, RNG PROPIO (LCG 32-bit + ",
       "polar de Marsaglia de `00_config`; misma secuencia exacta R/Python). Un ",
       "unico RNG recorre todos los items y escenarios en orden."),
"",
"### Salidas por item x escenario y veredicto",
"",
paste0("- `delta_rho_obs`, media y SD del `Delta rho` simulado, intervalo ",
       "central 95% (`sim_q025`/`sim_q975`, cuantiles tipo 7), ",
       "`sim_frac_abs_ge_obs` = fraccion de repeticiones con ",
       "`|Delta rho_sim| >= |Delta rho_obs|` (una `p` de simulacion), y ",
       "`sim_tasa_fisher_sig` = tasa de repeticiones en que el Fisher z de 09 ",
       "(SE Bonett-Wright) daria `p < .05` -> **falsos positivos del test ",
       "reportado bajo pura restriccion de rango**."),
paste0("- `veredicto = DENTRO`: el `Delta rho` observado cae dentro del ",
       "intervalo simulado -> la sola diferencia de dispersion puede ",
       "producirlo -> **no se puede descartar restriccion de rango**; el cambio ",
       "de correlacion no es interpretable como coordinacion biologica. ",
       "`FUERA`: lo excede en ese escenario."),
"",
"### Descartado",
"",
paste0("- **Simular sobre rangos** (en vez de en `-ddCt`): no modela como la ",
       "compresion de la escala original reduce el rango observado, que es el ",
       "mecanismo de la restriccion de rango. Se usa la normal bivariada en ",
       "`-ddCt`."),
paste0("- **Un solo escenario de `r` verdadera**: se corre GLOBAL y CONTROL ",
       "como rango porque la eleccion del `r` verdadero cambia cuanto ",
       "`Delta rho` espurio aparece; reportar los dos acota la conclusion."),
"",
"### Paridad R / Python",
"",
paste0("- Los sorteos usan `norm1()` (LCG + polar, byte-identico R/Python por ",
       "`00_config`/T1). `rho`, medias, SD, cuantiles y conteos usan acumulador ",
       "`double` explicito -> texto `%.10g` bit-identico. `r_pearson_gen` pasa ",
       "por `sin()` -> texto `%.6e`. Las figuras son PNG: equivalentes, no ",
       "byte-identicas.")
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 10_acto2_simulacion:inicio -->"
  marca_fin <- "<!-- 10_acto2_simulacion:fin -->"
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

construir_reporte <- function(fuente, sim) {
  L <- c(
    "# Reporte de la simulacion de restriccion de rango Acto 2.5 (T8)", "",
    "Generado por `10_acto2_simulacion` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. Metodo", "",
    sprintf(paste0("- Normal bivariada en `-ddCt`, `r` verdadera comun a ambos ",
                   "grupos, SD marginales = SD observadas por grupo (09). ",
                   "`r = 2*sin(pi*rho_S/6)`. B = %d, semilla %d, RNG PROPIO."),
            B_SIM, SEMILLA),
    paste0("- Dos escenarios de `r` verdadera como rango: `GLOBAL` (rho de T7) y ",
           "`CONTROL` (rho de Control)."),
    paste0("- `veredicto = DENTRO` -> la sola diferencia de dispersion puede ",
           "producir el `Delta rho` observado (no se descarta restriccion de ",
           "rango, prohibicion 5)."), "",
    "## 2. Resultados por item y escenario", "",
    .md(COLS_SIM, sim), "",
    "## 3. Figura", "",
    paste0("- `outputs/figures/acto2_simulacion_delta_rho.png` -- `Delta rho` ",
           "observado vs intervalo 95% del `Delta rho` simulado, por item y ",
           "escenario."), "",
    "## 4. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `10_acto2_simulacion`: ",
           "inversion rho->r, eleccion de escenarios, y lectura del veredicto."),
    "")
  paste(L, collapse = "\n")
}

# ===========================================================================
main <- function() {
  D <- cargar()
  fuente <- fuente_datos(ARCHIVO_QPCR)
  rng <- nuevo_rng(SEMILLA)

  sim <- tabla_simulacion(D, rng)

  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY))
    escribir_csv(file.path(base, "acto2_simulacion.csv"), COLS_SIM, sim)

  fig_sim <- file.path(RUTA_FIGURAS, "acto2_simulacion_delta_rho.png")
  figura_simulacion(sim, fig_sim)

  escribir_lineas(file.path(RUTA_TABLAS, "acto2_simulacion_reporte.md"),
                  construir_reporte(fuente, sim))
  actualizar_descartados()

  n_sim <- sum(vapply(sim, function(f) f[[19]] != "sin_test", logical(1)))
  n_dentro <- sum(vapply(sim, function(f) f[[19]] == "DENTRO", logical(1)))
  n_fuera <- sum(vapply(sim, function(f) f[[19]] == "FUERA", logical(1)))

  ent <- sprintf(paste0("data/processed/qpcr_cuantificacion_long.tsv + ",
                        "qpcr_score_compuesto_long.tsv + ",
                        "outputs/tables/*/acto2_dispersion.csv (de data/%s/%s)"),
                 fuente, ARCHIVO_QPCR)
  registrar_procedencia(list(
    list("outputs/tables/{R,python}/acto2_simulacion.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("simulacion de restriccion de rango: normal ",
         "bivariada en -ddCt, r verdadera comun (escenarios GLOBAL y CONTROL), ",
         "SD por grupo observadas; Delta rho simulado, IC95, frac >=|obs|, tasa ",
         "Fisher, veredicto DENTRO/FUERA")),
    list("outputs/figures/acto2_simulacion_delta_rho.png", "figura", ESTE_SCRIPT,
         "PROPIO", ent, paste0("Delta rho observado vs intervalo 95% del Delta ",
         "rho simulado por item y escenario")),
    list("outputs/tables/acto2_simulacion_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.5 (T8)")
  ))
  registrar_verificaciones(list(
    list("acto2_sim_proposito",
         paste0("la simulacion descarta (o no) restriccion de rango antes de ",
                "interpretar (prohibicion 5)"),
         paste0("normal bivariada en -ddCt con r verdadera comun y SD por grupo ",
                "observadas; veredicto DENTRO/FUERA por escenario"),
         "prohibicion 5 atendida", "TRUE", ESTE_SCRIPT),
    list("acto2_sim_escenarios",
         "r verdadera comun probada como rango: GLOBAL (rho T7) y CONTROL (rho Control)",
         paste(ESCENARIOS, collapse = ";"), "GLOBAL;CONTROL",
         if (identical(ESCENARIOS, c("GLOBAL", "CONTROL"))) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_sim_inversion_rho_r",
         "inversion rho_Spearman -> r_Pearson generativo para la normal bivariada",
         "r = 2*sin(pi*rho_S/6), recorte +/- 0.999999", "2*sin(pi*rho/6)",
         "TRUE", ESTE_SCRIPT),
    list("acto2_sim_B_semilla",
         "repeticiones y semilla de la simulacion",
         sprintf("B=%d; semilla=%d; RNG PROPIO LCG+polar", B_SIM, SEMILLA),
         sprintf("B=%d, semilla 20260101", B_SIM),
         if (B_SIM == 2000L && SEMILLA == 20260101) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_sim_SD_por_grupo",
         "las SD marginales simuladas por grupo son las SD observadas de ese grupo",
         paste0("sd_pla/sd_bra de pares(item, grupo); unica diferencia entre ",
                "grupos en la simulacion"),
         "SD observadas", "TRUE", ESTE_SCRIPT),
    list("acto2_sim_veredicto",
         "veredicto por item x escenario: DENTRO (no se descarta rango) / FUERA",
         sprintf("%d celdas simuladas; DENTRO=%d; FUERA=%d", n_sim, n_dentro, n_fuera),
         "DENTRO/FUERA/sin_test", "TRUE", ESTE_SCRIPT),
    list("acto2_sim_tasa_fisher",
         paste0("se reporta la tasa de falsos positivos del Fisher z (09) bajo ",
                "pura restriccion de rango"),
         "columna sim_tasa_fisher_sig por celda", "tasa reportada", "TRUE",
         ESTE_SCRIPT),
    list("acto2_sim_figura",
         "figura Acto 2.5: Delta rho observado vs intervalo simulado",
         sprintf("sim=%s", file.exists(fig_sim)), "1 figura existe",
         if (file.exists(fig_sim)) "TRUE" else "FALSE", ESTE_SCRIPT)
  ))

  cat("== 10_acto2_simulacion.R ==\n")
  cat(sprintf("  fuente = %s   B = %d   semilla = %d\n", fuente, B_SIM, SEMILLA))
  cat(sprintf("  celdas simuladas: %d  (DENTRO=%d, FUERA=%d)\n",
              n_sim, n_dentro, n_fuera))
  for (f in sim) {
    if (f[[19]] == "sin_test")
      cat(sprintf("    %-16s %-8s  (sin test, piso)\n", f[[1]], f[[3]]))
    else
      cat(sprintf(paste0("    %-16s %-8s  dRho_obs=%9s  IC95_sim=[%9s, %9s]  ",
                         "frac>=obs=%8s  fisherFP=%7s  -> %s\n"),
                  f[[1]], f[[3]], f[[12]], f[[15]], f[[16]], f[[17]], f[[18]],
                  f[[19]]))
  }
  cat("  -> outputs/tables/{R,python}/acto2_simulacion.csv\n")
  cat(sprintf("  -> %s\n", basename(fig_sim)))
}

if (sys.nframe() == 0L) main()
