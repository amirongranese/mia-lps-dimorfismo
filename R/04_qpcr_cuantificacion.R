# 04_qpcr_cuantificacion.R -- Cuantificacion relativa de la qPCR (D1, D2, D8).
#
# Por que existe este archivo: convierte los CT crudos de 02_ingesta_qc en las
# medidas sobre las que modela T5. Aplica -- sin imputar NADA -- las decisiones:
#
#   * D1: dCt = CT_gen - CT_rsp29 por muestra. Calibrador = HEMBRA_CONTROL, por
#     gen x tejido, promediando SOLO valores detectados. ddCt = dCt - dCt_calib.
#   * D2: la medida de analisis es -ddCt. FC = 2^(-ddCt) es solo para graficar
#     (T6/T7): NO se guarda aca para no arrastrar el redondeo de 2^x entre lenguajes.
#   * D7: si el calibrador HEMBRA_CONTROL de un gen x tejido tiene 0 detectados,
#     ese gen x tejido NO es cuantificable (ddCt/-ddCt/z NA, cuantificable=FALSE).
#     Se detecta programaticamente. Hoy: solo il6 @ BRAIN_E15.
#   * D8: z-score de cada gen por separado, dentro de cada tejido, sobre los 36
#     fetos, usando SOLO los detectados; desvio estandar MUESTRAL (n-1). Score
#     compuesto = promedio de los 7 z de transportadores por feto x tejido.
#     PROHIBIDO llamarlo "TONE".
#
# D3: los no detectados quedan NA y no se imputan. Se emite
# qpcr_no_detectados_descriptivo.csv (n y % de NA por gen x tejido x grupo +
# marca de calibrador sin detectados) como sustento numerico de D3 y D7.
#
# Este script NO modela y NO grafica.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

ESTE_SCRIPT <- "04_qpcr_cuantificacion"

CAL_GRUPO   <- "HEMBRA_CONTROL"
GRUPOS_4    <- c("HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")
SET_TRANSP  <- GENES_TRANSPORTADORES

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 01/02/03.
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

# ---------------------------------------------------------------------------
# Aritmetica -- acumulador double explicito, mismo orden de suma que Python
# (R sum() usa long double en Windows/MinGW -> se evita a proposito).
# ---------------------------------------------------------------------------
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

# ===========================================================================
# 1. Lectura + indexado
# ===========================================================================
cargar <- function() {
  ruta <- file.path(RUTA_DATOS_PROC, "qpcr_e15_long.tsv")
  if (!file.exists(ruta))
    stop(sprintf("04: falta %s (correr 02_ingesta_qc primero)", ruta))
  read.delim(ruta, sep = "\t", quote = "", stringsAsFactors = FALSE,
             colClasses = "character", check.names = FALSE, encoding = "UTF-8",
             na.strings = character(0))
}

indexar <- function(raw) {
  # metadata por feto (primer match), ordenada por (MADRE_ID, FETO) codepoint
  key <- !duplicated(raw$FETO)
  meta <- data.frame(MADRE_ID = raw$MADRE_ID[key], FETO = raw$FETO[key],
                     SEXO = raw$SEXO[key], TTO = raw$TTO[key],
                     GRUPO = raw$GRUPO[key], stringsAsFactors = FALSE)
  meta <- meta[order(meta$MADRE_ID, meta$FETO, method = "radix"), , drop = FALSE]

  k <- paste(raw$FETO, raw$TEJIDO, raw$GEN, sep = "\r")
  ct  <- setNames(suppressWarnings(as.numeric(raw$CT)), k)
  nod <- setNames(raw$no_detectado == "TRUE", k)
  kr  <- paste(raw$FETO, raw$TEJIDO, sep = "\r")
  rsp <- suppressWarnings(as.numeric(raw$rsp29))
  rsp <- rsp[!duplicated(kr)]; names(rsp) <- kr[!duplicated(kr)]
  list(meta = meta, ct = ct, nod = nod, rsp = rsp)
}

kg <- function(feto, tej, gen) paste(feto, tej, gen, sep = "\r")
kt <- function(feto, tej) paste(feto, tej, sep = "\r")

# ===========================================================================
# 2. Calibrador (D1) + regla D7
# ===========================================================================
calcular_calibradores <- function(ix) {
  cal_fetos <- ix$meta$FETO[ix$meta$GRUPO == CAL_GRUPO]
  calib <- list()
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    vals <- vapply(cal_fetos, function(f) {
      c <- ix$ct[[kg(f, tej, gen)]]
      if (is.na(c)) NA_real_ else c - ix$rsp[[kt(f, tej)]]
    }, numeric(1))
    det <- vals[!is.na(vals)]
    calib[[paste(tej, gen, sep = "\r")]] <- list(
      n_total = length(vals), n_detectado = length(det),
      dct_calibrador = if (length(det)) promedio(det) else NA_real_,
      cuantificable = length(det) > 0
    )
  }
  calib
}

# ===========================================================================
# 3. dCt / ddCt / -ddCt / z (D2, D8) + score compuesto
# ===========================================================================
cuantificar <- function(ix, calib) {
  meta <- ix$meta
  long <- list(); negdd <- list()
  for (i in seq_len(nrow(meta))) {
    m <- meta[i, ]
    for (tej in TEJIDOS_E15) for (gen in GENES) {
      k <- kg(m$FETO, tej, gen)
      c <- ix$ct[[k]]
      d <- if (is.na(c)) NA_real_ else c - ix$rsp[[kt(m$FETO, tej)]]
      cal <- calib[[paste(tej, gen, sep = "\r")]]
      dcal <- cal$dct_calibrador
      dd <- if (is.na(d) || is.na(dcal)) NA_real_ else d - dcal
      nd <- if (is.na(dd)) NA_real_ else -dd
      negdd[[k]] <- nd
      long[[length(long) + 1L]] <- list(
        MADRE_ID = m$MADRE_ID, FETO = m$FETO, SEXO = m$SEXO, TTO = m$TTO,
        GRUPO = m$GRUPO, TEJIDO = tej, GEN = gen,
        es_transportador = gen %in% SET_TRANSP,
        no_detectado = ix$nod[[k]], cuantificable = cal$cuantificable,
        rsp29 = ix$rsp[[kt(m$FETO, tej)]], CT = c,
        dCt = d, dCt_calibrador = dcal, ddCt = dd, neg_ddCt = nd
      )
    }
  }

  # z-score por gen x tejido sobre los 36 fetos (solo detectados, sd n-1)
  zt <- list(); zstats <- list()
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    vals <- vapply(meta$FETO, function(f) {
      v <- negdd[[kg(f, tej, gen)]]; if (is.null(v)) NA_real_ else v
    }, numeric(1))
    det <- vals[!is.na(vals)]
    if (length(det) >= 2L) {
      mu <- promedio(det); sd <- desvio_muestral(det, mu)
      for (j in seq_along(meta$FETO)) {
        v <- vals[j]
        zt[[kg(meta$FETO[j], tej, gen)]] <-
          if (is.na(v) || is.na(sd) || sd == 0) NA_real_ else (v - mu) / sd
      }
      zstats[[paste(tej, gen, sep = "\r")]] <- c(length(det), mu, sd)
    } else {
      for (f in meta$FETO) zt[[kg(f, tej, gen)]] <- NA_real_
      zstats[[paste(tej, gen, sep = "\r")]] <- c(length(det), NA, NA)
    }
  }
  for (r in seq_along(long))
    long[[r]]$z <- zt[[kg(long[[r]]$FETO, long[[r]]$TEJIDO, long[[r]]$GEN)]]

  # score compuesto por feto x tejido
  score <- list()
  for (i in seq_len(nrow(meta))) {
    m <- meta[i, ]
    for (tej in TEJIDOS_E15) {
      zs <- vapply(SET_TRANSP, function(g) {
        v <- zt[[kg(m$FETO, tej, g)]]; if (is.null(v)) NA_real_ else v
      }, numeric(1))
      zs <- zs[!is.na(zs)]
      score[[length(score) + 1L]] <- list(
        MADRE_ID = m$MADRE_ID, FETO = m$FETO, SEXO = m$SEXO, TTO = m$TTO,
        GRUPO = m$GRUPO, TEJIDO = tej, n_z_disponibles = length(zs),
        score_compuesto = if (length(zs)) promedio(zs) else NA_real_
      )
    }
  }
  list(long = long, score = score, zstats = zstats)
}

# ===========================================================================
# 4. Tabla descriptiva de no detectados (sustento D3 / D7)
# ===========================================================================
tabla_no_detectados <- function(ix, calib) {
  header <- c("TEJIDO", "GEN", "GRUPO", "n_total", "n_no_detectado",
              "pct_no_detectado", "calibrador_cero_detectados")
  meta <- ix$meta
  filas <- list()
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    cal_cero <- !calib[[paste(tej, gen, sep = "\r")]]$cuantificable
    for (grp in c(GRUPOS_4, "TODOS")) {
      fs <- if (grp == "TODOS") meta$FETO else meta$FETO[meta$GRUPO == grp]
      n <- length(fs)
      nnd <- sum(vapply(fs, function(f) isTRUE(ix$nod[[kg(f, tej, gen)]]), logical(1)))
      filas[[length(filas) + 1L]] <- list(tej, gen, grp, n, nnd,
        if (n) 100 * nnd / n else NA_real_, cal_cero)
    }
  }
  list(header = header, filas = filas)
}

# ===========================================================================
# 5. Resumenes
# ===========================================================================
tabla_calibradores <- function(calib) {
  header <- c("TEJIDO", "GEN", "es_transportador", "n_calibrador",
              "n_calibrador_detectado", "dCt_calibrador", "cuantificable")
  filas <- list()
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    c <- calib[[paste(tej, gen, sep = "\r")]]
    filas[[length(filas) + 1L]] <- list(tej, gen, gen %in% SET_TRANSP,
      c$n_total, c$n_detectado, c$dct_calibrador, c$cuantificable)
  }
  list(header = header, filas = filas)
}

tabla_resumen_cuant <- function(long, calib) {
  header <- c("TEJIDO", "GEN", "es_transportador", "cuantificable",
              "n_detectado_de_36", "mean_neg_ddCt", "sd_neg_ddCt",
              "median_neg_ddCt")
  filas <- list()
  for (tej in TEJIDOS_E15) for (gen in GENES) {
    vals <- unlist(lapply(long, function(r)
      if (r$TEJIDO == tej && r$GEN == gen && !is.na(r$neg_ddCt)) r$neg_ddCt else NULL))
    n <- length(vals)
    mu <- if (n) promedio(vals) else NA_real_
    sd <- if (n >= 2L) desvio_muestral(vals, mu) else NA_real_
    md <- if (n) mediana(vals) else NA_real_
    filas[[length(filas) + 1L]] <- list(tej, gen, gen %in% SET_TRANSP,
      calib[[paste(tej, gen, sep = "\r")]]$cuantificable, n, mu, sd, md)
  }
  list(header = header, filas = filas)
}

tabla_resumen_score <- function(score) {
  header <- c("TEJIDO", "GRUPO", "n_fetos", "mean_score_compuesto",
              "sd_score_compuesto")
  filas <- list()
  for (tej in TEJIDOS_E15) for (grp in GRUPOS_4) {
    vals <- unlist(lapply(score, function(r)
      if (r$TEJIDO == tej && r$GRUPO == grp && !is.na(r$score_compuesto))
        r$score_compuesto else NULL))
    n <- length(vals)
    mu <- if (n) promedio(vals) else NA_real_
    sd <- if (n >= 2L) desvio_muestral(vals, mu) else NA_real_
    filas[[length(filas) + 1L]] <- list(tej, grp, n, mu, sd)
  }
  list(header = header, filas = filas)
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

construir_reporte <- function(fuente, cal, nd, res, sc, no_cuant) {
  L <- c(
    "# Reporte de cuantificacion qPCR (T4)", "",
    "Generado por `04_qpcr_cuantificacion` (R y Python producen este archivo identico).",
    sprintf("Fuente de datos en uso: `%s`.", fuente), "",
    "## 1. Metodo (D1, D2, D8) -- sin imputacion (D3)", "",
    paste0("- `dCt = CT_gen - CT_rsp29` por muestra. Calibrador = promedio de ",
           "`dCt` en **HEMBRA_CONTROL**, por gen x tejido, **solo detectados** (D1)."),
    paste0("- `ddCt = dCt - dCt_calibrador`; la medida de analisis es `-ddCt` (D2). ",
           "`FC = 2^(-ddCt)` se calcula al graficar, no se guarda."),
    paste0("- z-score por gen dentro de tejido sobre los 36 fetos, solo detectados, ",
           "**desvio muestral (n-1)** (D8)."),
    paste0("- Score compuesto = promedio de los z de los 7 transportadores por ",
           "feto x tejido (los disponibles si hay <7) (D8)."), "",
    "## 2. Calibradores por gen x tejido (D1) y regla D7", ""
  )
  L <- c(L, if (length(no_cuant))
    sprintf(paste0("**No cuantificable(s) por D7** (calibrador HEMBRA_CONTROL con ",
                   "0 detectados): %s."), paste(no_cuant, collapse = ", "))
    else "Todos los gen x tejido tienen >=1 detectado en el calibrador.")
  L <- c(L, "", .md_tabla(cal$header, cal$filas), "",
    "## 3. No detectados por gen x tejido x grupo (sustento numerico de D3 y D7)", "",
    paste0("Columna `calibrador_cero_detectados` = TRUE marca el gen x tejido donde ",
           "el calibrador HEMBRA_CONTROL no tiene ningun detectado: ahi la imputacion ",
           "MNAR no tendria nada sobre lo que anclar (se argumenta en T10). Tabla ",
           "completa en `qpcr_no_detectados_descriptivo.csv`; aca solo la fila `TODOS` ",
           "por gen x tejido."), "",
    .md_tabla(nd$header, Filter(function(f) f[[3]] == "TODOS", nd$filas)), "",
    "## 4. Resumen de -ddCt por gen x tejido (36 fetos)", "",
    .md_tabla(res$header, res$filas), "",
    "## 5. Score compuesto por grupo x tejido", "",
    .md_tabla(sc$header, sc$filas), "",
    "## 6. Notas", "",
    paste0("Ver `analisis_descartados.md`, seccion `04_qpcr_cuantificacion`: ",
           "il6 @ BRAIN_E15 fuera de la cuantificacion (D7), y el recordatorio de que ",
           "no se imputa (D3) con el puntero a la tabla de no detectados."), "")
  paste(L, collapse = "\n")
}

# ===========================================================================
# 7. Artefactos compartidos (merge por 'script') -- headers identicos a 02/03.
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
"## 04_qpcr_cuantificacion",
"",
"### il6 @ BRAIN_E15: fuera de la cuantificacion (D7)",
"",
"- El calibrador HEMBRA_CONTROL de il6 @ BRAIN_E15 tiene 0/9 detectados (ver",
"  `qpcr_no_detectados_descriptivo.csv`: TEJIDO=BRAIN_E15, GEN=il6,",
"  GRUPO=HEMBRA_CONTROL, n_no_detectado = n_total). Sin dCt de calibrador no hay",
"  ddCt: se marca `cuantificable = FALSE`, sus `ddCt` / `-ddCt` / `z` quedan NA,",
"  y el gen entra al analisis solo como proporcion de deteccion en T5 (D7).",
"- Es el unico gen x tejido donde ocurre. Se detecta programaticamente:",
"  gen x tejido con HEMBRA_CONTROL `n_detectado == 0 & n_total > 0`.",
"",
"### No se imputan los no detectados (D3) -- sustento numerico",
"",
"- Los CT no detectados (`CT_CRUDO == 40` o celda vacia, unificados en",
"  02_ingesta_qc) quedan NA y NO se imputan por ningun metodo. El conteo y % de",
"  NA por gen x tejido x grupo esta en `qpcr_no_detectados_descriptivo.csv`.",
"- El calibrador y todos los promedios/desvios (dCt de calibrador, media y sd",
"  del z-score) se calculan **solo sobre valores detectados** (D1/D8). Ningun",
"  `dCt` se rellena: la cantidad de `dCt` NA es exactamente la de no detectados.",
"- La justificacion completa del descarte de la imputacion MNAR (`nondetects`) se",
"  redacta en T10 apoyandose en esta tabla: la imputacion MNAR estima lo no",
"  observado a partir de un modelo de la censura, y en las celdas donde el",
"  calibrador no tiene detectados (D7) no hay nada sobre lo que anclar la",
"  estimacion, asi que el fold-change quedaria definido contra un valor inventado."
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 04_qpcr_cuantificacion:inicio -->"
  marca_fin <- "<!-- 04_qpcr_cuantificacion:fin -->"
  nuevo <- paste0(marca_ini, "\n", BLOQUE_DESCARTES, "\n\n", marca_fin)
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
filas_de <- function(lst, cols) lapply(lst, function(r) lapply(cols, function(c) r[[c]]))

main <- function() {
  raw <- cargar()
  fuente <- fuente_datos(ARCHIVO_QPCR)
  ix <- indexar(raw)
  calib <- calcular_calibradores(ix)
  q <- cuantificar(ix, calib)
  long <- q$long; score <- q$score

  no_cuant <- character(0)
  for (tej in TEJIDOS_E15) for (gen in GENES)
    if (!calib[[paste(tej, gen, sep = "\r")]]$cuantificable)
      no_cuant <- c(no_cuant, paste0(gen, "@", tej))

  cols_long <- c("MADRE_ID", "FETO", "SEXO", "TTO", "GRUPO", "TEJIDO", "GEN",
                 "es_transportador", "no_detectado", "cuantificable", "rsp29",
                 "CT", "dCt", "dCt_calibrador", "ddCt", "neg_ddCt", "z")
  escribir_tsv(file.path(RUTA_DATOS_PROC, "qpcr_cuantificacion_long.tsv"),
               cols_long, filas_de(long, cols_long))
  cols_sc <- c("MADRE_ID", "FETO", "SEXO", "TTO", "GRUPO", "TEJIDO",
               "n_z_disponibles", "score_compuesto")
  escribir_tsv(file.path(RUTA_DATOS_PROC, "qpcr_score_compuesto_long.tsv"),
               cols_sc, filas_de(score, cols_sc))

  cal <- tabla_calibradores(calib)
  nd  <- tabla_no_detectados(ix, calib)
  res <- tabla_resumen_cuant(long, calib)
  sc  <- tabla_resumen_score(score)
  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "qpcr_calibradores.csv"), cal$header, cal$filas)
    escribir_csv(file.path(base, "qpcr_no_detectados_descriptivo.csv"), nd$header, nd$filas)
    escribir_csv(file.path(base, "qpcr_cuantificacion_resumen.csv"), res$header, res$filas)
    escribir_csv(file.path(base, "qpcr_score_compuesto_resumen.csv"), sc$header, sc$filas)
  }

  escribir_lineas(file.path(RUTA_TABLAS, "qpcr_cuantificacion_reporte.md"),
                  construir_reporte(fuente, cal, nd, res, sc, no_cuant))
  actualizar_descartados()

  n_dct_na <- sum(vapply(long, function(r) is.na(r$dCt), logical(1)))
  n_nodet  <- sum(vapply(long, function(r) isTRUE(r$no_detectado), logical(1)))

  ent <- sprintf("data/processed/qpcr_e15_long.tsv (de data/%s/%s)", fuente, ARCHIVO_QPCR)
  registrar_procedencia(list(
    list("data/processed/qpcr_cuantificacion_long.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "dCt, calibrador, ddCt, -ddCt, z por feto x tejido x gen (D1/D2/D8)"),
    list("data/processed/qpcr_score_compuesto_long.tsv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "score compuesto (promedio z de 7 transportadores) por feto x tejido (D8)"),
    list("outputs/tables/{R,python}/qpcr_calibradores.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent, "dCt de calibrador HEMBRA_CONTROL por gen x tejido + flag cuantificable (D1/D7)"),
    list("outputs/tables/{R,python}/qpcr_no_detectados_descriptivo.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "n y % de no detectados por gen x tejido x grupo + calibrador_cero_detectados (sustento D3/D7)"),
    list("outputs/tables/{R,python}/qpcr_cuantificacion_resumen.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent,
         "n detectados, media/sd/mediana de -ddCt por gen x tejido"),
    list("outputs/tables/{R,python}/qpcr_score_compuesto_resumen.csv", "tabla",
         ESTE_SCRIPT, "PROPIO", ent, "media/sd del score compuesto por grupo x tejido"),
    list("outputs/tables/qpcr_cuantificacion_reporte.md", "reporte", ESTE_SCRIPT,
         "PROPIO", ent, "reporte legible de la cuantificacion qPCR")
  ))
  registrar_verificaciones(list(
    list("cuant_filas_long", "qpcr_cuantificacion_long tiene 36x2x10 = 720 filas",
         as.character(length(long)), "720", if (length(long) == 720L) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("cuant_score_filas", "score compuesto: 36x2 = 72 filas",
         as.character(length(score)), "72", if (length(score) == 72L) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("cuant_no_cuantificables", "gen x tejido no cuantificables por D7",
         if (length(no_cuant)) paste(no_cuant, collapse = ";") else "(ninguno)",
         "il6@BRAIN_E15", if (identical(no_cuant, "il6@BRAIN_E15")) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("cuant_sin_imputacion", "cantidad de dCt NA == cantidad de no detectados (no se imputa)",
         sprintf("%d==%d", n_dct_na, n_nodet), "iguales",
         if (n_dct_na == n_nodet) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("cuant_sd_metodo", "z-score con desvio estandar muestral (n-1)",
         "muestral_n-1", "muestral_n-1", "TRUE", ESTE_SCRIPT),
    list("cuant_transportadores_score", "score compuesto sobre 7 transportadores",
         as.character(length(GENES_TRANSPORTADORES)), "7",
         if (length(GENES_TRANSPORTADORES) == 7L) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("cuant_calibrador_grupo", "calibrador de la cuantificacion relativa (D1)",
         CAL_GRUPO, "HEMBRA_CONTROL", if (CAL_GRUPO == "HEMBRA_CONTROL") "TRUE" else "FALSE",
         ESTE_SCRIPT)
  ))

  cat("== 04_qpcr_cuantificacion.R ==\n")
  cat(sprintf("  fuente qPCR = %s\n", fuente))
  cat(sprintf("  filas long = %d   score compuesto filas = %d\n", length(long), length(score)))
  cat(sprintf("  no cuantificable (D7): %s\n",
              if (length(no_cuant)) paste(no_cuant, collapse = ", ") else "(ninguno)"))
  cat(sprintf("  dCt NA = %d  (== no detectados = %d) -> sin imputar\n", n_dct_na, n_nodet))
  for (f in sc$filas)
    cat(sprintf("    score %-12s %-14s n=%2d  media=%s  sd=%s\n",
                f[[1]], f[[2]], f[[3]], .fmt(f[[4]]), .fmt(f[[5]])))
  cat("  -> data/processed/qpcr_{cuantificacion,score_compuesto}_long.tsv,\n")
  cat("     outputs/tables/{R,python}/qpcr_*.csv, outputs/tables/qpcr_cuantificacion_reporte.md\n")
}

if (sys.nframe() == 0L) main()
