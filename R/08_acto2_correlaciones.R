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
#   * Magnitudes (decision del usuario): -ddCt por gen (9 genes: todos menos il6,
#     que no tiene -ddCt en cerebro por D7) y el score compuesto (D8).
#   * Coeficiente (decision del usuario): **Spearman rho** (no Pearson).
#     rho = Pearson sobre los rangos promedio (corrige empates). p por la
#     t-aproximacion t = rho*sqrt((n-2)/(1-rho^2)), df = n-2, dos colas -- la
#     misma que cor.test(method="spearman", exact=FALSE) y scipy.spearmanr.
#     IC 95% Bonett-Wright: SE_z = sqrt((1 + rho^2/2)/(n-3)), z = atanh(rho),
#     IC = tanh(z +/- 1.959963984540054 * SE_z).
#   * Estratos (decision del usuario): GLOBAL + por TTO (CONTROL / LPS).
#   * Piso: n_par < 5 -> solo n (sin rho/IC/p).
#   * Co-expresion: rho de Spearman entre los 7 transportadores dentro de cada
#     tejido (tabla + SPLOM con GGally).
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

suppressMessages({ library(ggplot2); library(GGally) })

ESTE_SCRIPT <- "08_acto2_correlaciones"

COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")
Z975 <- 1.959963984540054
PISO_PAR <- 5L

GEN_SIN_CEREBRO <- "il6"
GENES_CORR <- setdiff(GENES, GEN_SIN_CEREBRO)          # 9 genes
ITEMS <- c(GENES_CORR, "score_compuesto")
ESTRATOS <- c("GLOBAL", "CONTROL", "LPS")
TRANSP <- GENES_TRANSPORTADORES                        # 7, para el SPLOM
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

COLS_COEXP <- c("TEJIDO", "GEN_A", "GEN_B", "n_par", "rho_spearman", "p_valor")

tabla_coexpresion <- function(D) {
  filas <- list()
  for (tej in TEJIDOS_E15) {
    for (i in seq_len(length(TRANSP) - 1L)) {
      for (j in (i + 1L):length(TRANSP)) {
        a <- TRANSP[i]; b <- TRANSP[j]
        xs <- c(); ys <- c()
        for (f in D$fetos) {
          va <- D$negdd[[paste(f, tej, a, sep = "\r")]]
          vb <- D$negdd[[paste(f, tej, b, sep = "\r")]]
          if (is.null(va) || is.null(vb) || is.na(va) || is.na(vb)) next
          xs <- c(xs, va); ys <- c(ys, vb)
        }
        n <- length(xs)
        if (n >= PISO_PAR) {
          rho <- spearman_rho(xs, ys); p <- spearman_p(rho, n)
          filas[[length(filas) + 1L]] <- list(tej, a, b, n, g10(rho), p6e(p))
        } else {
          filas[[length(filas) + 1L]] <- list(tej, a, b, n, "", "")
        }
      }
    }
  }
  filas
}

# ===========================================================================
# 3. Figuras.
# ===========================================================================
.rho_estratos <- function(D, item) {
  out <- list()
  for (est in ESTRATOS) {
    pr <- pares(D, item, est); n <- length(pr$x)
    out[[est]] <- if (n >= PISO_PAR) c(spearman_rho(pr$x, pr$y), n) else c(NA_real_, n)
  }
  out
}

figura_dispersion <- function(D, ruta) {
  df <- list(); ann <- list()
  for (item in ITEMS) {
    pr <- pares(D, item, "GLOBAL")
    nom <- if (item == "score_compuesto") "score compuesto" else item
    if (length(pr$x))
      df[[length(df) + 1L]] <- data.frame(
        item = nom, x = pr$x, y = pr$y, tto = pr$tto, stringsAsFactors = FALSE)
    rt <- .rho_estratos(D, item)
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
    facet_wrap(~ item, scales = "free", ncol = 4) +
    labs(title = paste0("Correlacion placenta <-> cerebro por feto  --  Spearman ",
                        "rho (-ddCt por gen y score compuesto)"),
         subtitle = paste0("T7 solo describe: la diferencia de rho entre Control ",
                           "y LPS NO se testea aca (prohibicion 4); el test ",
                           "formal es T8"),
         x = "placenta  -ddCt", y = "cerebro  -ddCt") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.text = element_text(size = 9),
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.subtitle = element_text(size = 8.5))
  ggsave(ruta, p, width = 13.0, height = 9.0, dpi = DPI)
}

figura_splom <- function(D, tej, ruta) {
  cols <- lapply(TRANSP, function(g)
    vapply(D$fetos, function(f) {
      v <- D$negdd[[paste(f, tej, g, sep = "\r")]]
      if (is.null(v)) NA_real_ else v
    }, numeric(1)))
  wide <- as.data.frame(setNames(cols, TRANSP), stringsAsFactors = FALSE)
  wide$TTO <- factor(vapply(D$fetos, function(f) D$tto[[f]], character(1)),
                     levels = c("CONTROL", "LPS"))
  tt <- if (tej == "PLACENTA_E15") "Placenta E15" else "Cerebro fetal E15"
  p <- GGally::ggpairs(
    wide, columns = seq_along(TRANSP),
    mapping = ggplot2::aes(colour = TTO),
    upper = list(continuous = GGally::wrap("cor", method = "spearman", size = 2.8)),
    lower = list(continuous = GGally::wrap("points", size = 0.5, alpha = 0.7)),
    diag  = list(continuous = GGally::wrap("densityDiag", alpha = 0.4)),
    title = sprintf(paste0("Co-expresion de los 7 transportadores (%s) -- -ddCt, ",
                           "Spearman rho"), tt)) +
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
       "cada gen y del score compuesto, global y por TTO. **No compara** las ",
       "correlaciones entre grupos: reportar \"significativo en Control y no en ",
       "LPS\" como prueba de diferencia esta prohibido (prohibicion 4). El test ",
       "formal (Fisher z / interaccion de pendientes / permutacion) y el control ",
       "de restriccion de rango (prohibicion 5, simulacion) son **T8**."),
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
       "(calibrador HEMBRA_CONTROL 0/9, D7). `il6R` en cerebro tiene n bajo -> sus ",
       "`n_par` por estrato pueden quedar por debajo del piso de 5."),
paste0("- Estratos: `GLOBAL` (n<=36) y por `TTO` (`CONTROL` / `LPS`, n<=18). Las ",
       "4 celdas SEXO x TTO (n~9) darian IC inutiles y no se usan."),
"",
"### Co-expresion (SPLOM)",
"",
paste0("- `acto2_coexpresion_transportadores.csv` y `acto2_coexpresion_SPLOM_",
       "<tejido>.png` muestran la rho de Spearman entre los 7 transportadores ",
       "dentro de cada tejido (21 pares x 2 tejidos). Es contexto, no una prueba."),
"",
"### Paridad R / Python",
"",
paste0("- `rho` y las sumas usan acumulador `double` explicito (mismo orden) -> ",
       "bit-identico; se guarda con `%.10g`. Lo que pasa por trascendentes (`p` ",
       "via `pt`, IC via `tanh`/`atanh`) se guarda como texto `%.6e`. Las figuras ",
       "son PNG: equivalentes, no byte-identicas (ggplot2/GGally vs matplotlib).")
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
    paste0("- Magnitudes: `-ddCt` de cada gen (9; `il6` fuera por D7) y el **score ",
           "compuesto** de 7 transportadores (D8)."),
    "- Estratos: `GLOBAL` y por `TTO` (`CONTROL` / `LPS`).",
    paste0("- **T7 no compara** las correlaciones entre grupos (prohibicion 4). El ",
           "test formal y el control de restriccion de rango son **T8**."),
    "", "## 2. Correlacion placenta <-> cerebro (por gen y score)", "",
    .md(COLS_CORR, corr), "",
    "## 3. Co-expresion entre transportadores (Spearman, por tejido)", "",
    .md(COLS_COEXP, coexp), "",
    "## 4. Figuras", "",
    paste0("- `outputs/figures/acto2_dispersion_placenta_cerebro.png` -- dispersion ",
           "placenta vs cerebro por gen + score, coloreada por TTO, con rho (global ",
           "y por grupo) anotada."),
    paste0("- `outputs/figures/acto2_coexpresion_SPLOM_PLACENTA_E15.png` y ",
           "`..._BRAIN_E15.png` -- matriz de dispersion de los 7 transportadores ",
           "por tejido."),
    "", "## 5. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `08_acto2_correlaciones`: ",
           "eleccion de Spearman, formula de `p` e IC, y el limite explicito de T7 ",
           "(describe, no compara)."), "")
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
  fig_spl_p <- file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_PLACENTA_E15.png")
  fig_spl_b <- file.path(RUTA_FIGURAS, "acto2_coexpresion_SPLOM_BRAIN_E15.png")
  figura_dispersion(D, fig_disp)
  figura_splom(D, "PLACENTA_E15", fig_spl_p)
  figura_splom(D, "BRAIN_E15", fig_spl_b)

  escribir_lineas(file.path(RUTA_TABLAS, "acto2_correlaciones_reporte.md"),
                  construir_reporte(fuente, corr, coexp))
  actualizar_descartados()

  n_items <- length(ITEMS)
  n_con_rho <- sum(vapply(corr, function(f) nzchar(f[[5]]), logical(1)))
  sc_glob <- Filter(function(f) f[[1]] == "score_compuesto" && f[[3]] == "GLOBAL", corr)[[1]]

  ent <- sprintf(paste0("data/processed/qpcr_cuantificacion_long.tsv + ",
                        "qpcr_score_compuesto_long.tsv (de data/%s/%s)"),
                 fuente, ARCHIVO_QPCR)
  registrar_procedencia(list(
    list("outputs/tables/{R,python}/acto2_correlaciones.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, paste0("Spearman rho placenta<->cerebro por feto: 9 genes + ",
         "score compuesto x {GLOBAL, CONTROL, LPS}; n_par, IC95 Bonett-Wright, p ",
         "t-aprox")),
    list("outputs/tables/{R,python}/acto2_coexpresion_transportadores.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, paste0("Spearman rho entre los 7 ",
         "transportadores dentro de cada tejido (21 pares x 2 tejidos)")),
    list("outputs/figures/acto2_dispersion_placenta_cerebro.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, paste0("dispersion placenta vs cerebro por gen ",
         "+ score, coloreada por TTO, con rho (global y por grupo) anotada")),
    list("outputs/figures/acto2_coexpresion_SPLOM_PLACENTA_E15.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, paste0("matriz de dispersion (SPLOM) de los 7 ",
         "transportadores en placenta E15")),
    list("outputs/figures/acto2_coexpresion_SPLOM_BRAIN_E15.png", "figura",
         ESTE_SCRIPT, "PROPIO", ent, "idem en cerebro fetal E15"),
    list("outputs/tables/acto2_correlaciones_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible del Acto 2.1-2.2 (T7)")
  ))
  registrar_verificaciones(list(
    list("acto2_items_correlacionados",
         "correlaciones placenta<->cerebro: 9 genes (sin il6) + score compuesto",
         sprintf("%d items x %d estratos; %d con rho (n_par>=5)",
                 n_items, length(ESTRATOS), n_con_rho),
         "9 genes + score_compuesto", if (n_items == 10L) "TRUE" else "FALSE",
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
    list("acto2_no_compara_grupos",
         "T7 describe pero NO compara correlaciones entre grupos (prohibicion 4)",
         "reporte y analisis_descartados lo dicen explicitamente; el test formal es T8",
         "no se compara en T7", "TRUE", ESTE_SCRIPT),
    list("acto2_pareo_por_feto",
         "el emparejamiento placenta<->cerebro es por FETO (ambos lados detectados)",
         "par = (PLACENTA_E15, BRAIN_E15) del mismo FETO con ambos valores no NA",
         "por feto", "TRUE", ESTE_SCRIPT),
    list("acto2_estratos",
         "estratos de correlacion: GLOBAL + por TTO (Control/LPS); sin celdas SEXOxTTO",
         paste(ESTRATOS, collapse = ";"), "GLOBAL;CONTROL;LPS",
         if (identical(ESTRATOS, c("GLOBAL", "CONTROL", "LPS"))) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_piso_par", "n_par minimo para calcular rho/IC/p",
         as.character(PISO_PAR), "5", if (PISO_PAR == 5L) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("acto2_figuras", "figuras Acto 2.1-2.2: dispersion por gen + 2 SPLOM",
         sprintf("disp=%s;splom_pla=%s;splom_bra=%s", file.exists(fig_disp),
                 file.exists(fig_spl_p), file.exists(fig_spl_b)),
         "3 figuras existen",
         if (file.exists(fig_disp) && file.exists(fig_spl_p) && file.exists(fig_spl_b))
           "TRUE" else "FALSE", ESTE_SCRIPT)
  ))

  cat("== 08_acto2_correlaciones.R ==\n")
  cat(sprintf("  fuente = %s\n", fuente))
  cat(sprintf("  items = %d (9 genes sin il6 + score_compuesto) x %d estratos\n",
              n_items, length(ESTRATOS)))
  cat(sprintf("  score compuesto GLOBAL: n_par=%s  rho=%s  IC=[%s, %s]  p=%s\n",
              sc_glob[[4]], sc_glob[[5]], sc_glob[[6]], sc_glob[[7]], sc_glob[[8]]))
  for (f in corr) if (f[[3]] == "GLOBAL")
    cat(sprintf("    %-16s GLOBAL  n=%2s  rho=%8s  p=%s\n", f[[1]], f[[4]], f[[5]], f[[8]]))
  cat(sprintf("  co-expresion: %d pares (21 x 2 tejidos)\n", length(coexp)))
  cat("  -> outputs/tables/{R,python}/acto2_correlaciones.csv, acto2_coexpresion_transportadores.csv\n")
  cat(sprintf("  -> %s, %s, %s\n", basename(fig_disp), basename(fig_spl_p),
              basename(fig_spl_b)))
}

if (sys.nframe() == 0L) main()
