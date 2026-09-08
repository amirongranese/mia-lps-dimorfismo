# 02_ingesta_qc.R -- Ingesta y control de calidad de los tres archivos de entrada.
#
# Por que existe este archivo: es la unica puerta de entrada de los datos al
# pipeline. Lee los tres Excel (crudo real si esta en data/raw/, si no el
# sintetico versionado), los deja en formato largo y limpio, y aplica -- sin
# tocar ninguna otra decision -- las tres reglas de saneamiento de la Seccion 2
# del brief:
#
#   * qPCR: CT_CRUDO == 40 y celda vacia son la MISMA cosa -> "no detectado" ->
#     NA (D3). Nunca se imputa. El tejido BRAIN_P1 se excluye del alcance E15
#     (se analiza en otro informe) y queda registrado en analisis_descartados.md.
#   * ELISA: la hoja mezcla suero materno y liquido amniotico; se parte en dos
#     bloques por la columna TEJIDO. Conc < 0 es censura a izquierda (D10): se
#     marca censurado=TRUE, el valor se guarda como NA y el LOD (= 0, el blanco)
#     se registra aparte. La columna "IL-6" del crudo NO se usa (ver
#     analisis_descartados.md).
#   * pSTAT3: se normaliza SEXO/TTO y se verifica el balanceo de las 3 membranas.
#
# Entregable: el "reporte de QC" -- tablas de n real por archivo x grupo x sexo,
# de no-detectados por gen x tejido x grupo, de censura del ELISA por bloque x
# grupo, y de que madres/fetos faltan en cada bloque respecto del diseno qPCR.
# Este script NO cuantifica, NO modela y NO grafica.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

suppressMessages(library(readxl))

# Identificador para los registros compartidos (procedencia, verificaciones,
# analisis_descartados): SIN extension, para que R y Python -- que producen los
# mismos artefactos -- reemplacen la misma fila y el archivo quede byte-identico.
ESTE_SCRIPT <- "02_ingesta_qc"

# ---------------------------------------------------------------------------
# Ordenes canonicos: todas las salidas se ordenan por estas claves para que la
# comparacion byte a byte con la implementacion Python sea posible.
# ---------------------------------------------------------------------------
GRUPOS_4     <- c("HEMBRA_CONTROL", "HEMBRA_LPS", "MACHO_CONTROL", "MACHO_LPS")
ORDEN_GRUPO  <- setNames(seq_along(GRUPOS_4) - 1L, GRUPOS_4)
ORDEN_TEJIDO <- setNames(seq_along(TEJIDOS_E15) - 1L, TEJIDOS_E15)
ORDEN_GEN    <- setNames(seq_along(GENES) - 1L, GENES)
ORDEN_TTO    <- c(CONTROL = 0L, LPS = 1L)
ORDEN_SEXO   <- setNames(c(0L, 1L, 2L), c("", "HEMBRA", "MACHO"))

grupo_norm <- function(sexo, tto) paste0(sexo, "_", tto)

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a los de 01_generar_sinteticos (paridad
# R/Python): NA -> "", float -> %.10g, entero -> sin decimales, logico -> TRUE/FALSE.
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
  con <- file(ruta, open = "wb")
  writeBin(charToRaw(enc2utf8(texto)), con)
  close(con)
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
s_txt <- function(v) {
  v <- as.character(v)
  v[is.na(v)] <- ""
  trimws(v)
}

to_num <- function(v) {
  x <- as.character(v)
  n <- suppressWarnings(as.numeric(x))
  falta <- is.na(n) & !is.na(x) & trimws(x) != ""
  if (any(falta)) {
    n2 <- suppressWarnings(as.numeric(gsub(",", ".", x[falta], fixed = TRUE)))
    n[falta] <- n2
  }
  n
}

madre_id_de <- function(feto) ifelse(feto == "", "", sub("\\.[0-9]+$", "", feto))

.leer_grid_xlsx <- function(ruta, hoja) {
  g <- suppressMessages(readxl::read_xlsx(ruta, sheet = hoja, col_names = FALSE,
                                          col_types = "text"))
  as.data.frame(g, stringsAsFactors = FALSE)
}

.leer_grid_tsv <- function(ruta) {
  txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE)
  Encoding(txt) <- "UTF-8"
  lineas <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  if (length(lineas) && lineas[length(lineas)] == "") lineas <- lineas[-length(lineas)]
  partes <- strsplit(lineas, "\t", fixed = TRUE)  # strsplit descarta vacios finales
  nc <- max(vapply(partes, length, integer(1)))
  m <- t(vapply(partes, function(p) { length(p) <- nc; p }, character(nc)))
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  df[!is.na(df) & df == ""] <- NA_character_
  df
}

# Grilla cruda (data.frame de texto, sin encabezado asumido).
# - Datos reales: .xlsx de data/raw/.
# - Datos sinteticos: mirror .tsv canonico de data/synthetic/ (forma versionada,
#   byte-identica R/Python; el .xlsx sintetico no se versiona). Fallback al .xlsx
#   sintetico si estuviera y no el .tsv.
leer_grid <- function(nombre_archivo, hoja) {
  if (fuente_datos(nombre_archivo) == "real")
    return(.leer_grid_xlsx(ruta_datos(nombre_archivo), hoja))
  tallo <- sub("\\.[^.]+$", "", nombre_archivo)
  for (tsv in c(file.path(RUTA_DATOS_SINT, paste0(tallo, ".tsv")),
                file.path(RUTA_DATOS_SINT, paste0(tallo, "__", hoja, ".tsv"))))
    if (file.exists(tsv)) return(.leer_grid_tsv(tsv))
  ruta <- ruta_datos(nombre_archivo)
  if (file.exists(ruta)) return(.leer_grid_xlsx(ruta, hoja))
  stop(sprintf("ingesta: no hay fuente para '%s' (hoja '%s')", nombre_archivo, hoja))
}

ordkey <- function(...) do.call(order, list(...))

# ===========================================================================
# 1. qPCR: lectura, 40/vacio -> NA, exclusion de BRAIN_P1
# ===========================================================================
COLS_QPCR <- c("MADRE", "NOMINACION", "FETO", "GRUPO", "SEXO", "TTO", "TEJIDO",
               "rsp29", "GEN", "CT_CRUDO")

ingesta_qpcr <- function() {
  g <- leer_grid(ARCHIVO_QPCR, HOJA_QPCR)
  encab <- s_txt(unlist(g[1, ], use.names = FALSE))
  if (!identical(encab[1:10], COLS_QPCR))
    stop(sprintf("qPCR: encabezado inesperado: %s", paste(encab[1:10], collapse = ", ")))

  raw <- g[-1, , drop = FALSE]
  cel <- lapply(1:10, function(j) s_txt(raw[[j]]))
  allblank <- Reduce(`&`, lapply(cel, function(v) v == ""))
  keep <- !allblank

  MADRE      <- cel[[1]][keep]
  NOMINACION <- cel[[2]][keep]
  FETO       <- cel[[3]][keep]
  GRUPO_RAW  <- cel[[4]][keep]
  SEXO_RAW   <- cel[[5]][keep]
  TTO_RAW    <- cel[[6]][keep]
  TEJIDO     <- cel[[7]][keep]
  rsp29      <- to_num(raw[[8]])[keep]
  GEN        <- cel[[9]][keep]
  cv         <- cel[[10]][keep]
  x          <- to_num(raw[[10]])[keep]

  es_vacio <- cv == ""
  es_40    <- !is.na(x) & x == 40
  n_ct_no_num <- sum(cv != "" & is.na(x))
  no_det <- es_vacio | es_40

  SEXO <- vapply(SEXO_RAW, normalizar_sexo, character(1), USE.NAMES = FALSE)
  TTO  <- vapply(TTO_RAW,  normalizar_tto,  character(1), USE.NAMES = FALSE)

  df <- data.frame(
    MADRE = MADRE, MADRE_ID = madre_id_de(FETO), NOMINACION = NOMINACION,
    FETO = FETO, GRUPO_RAW = GRUPO_RAW, SEXO_RAW = SEXO_RAW, TTO_RAW = TTO_RAW,
    SEXO = SEXO, TTO = TTO, GRUPO = grupo_norm(SEXO, TTO), TEJIDO = TEJIDO,
    rsp29 = rsp29, GEN = GEN, CT_CRUDO = x,
    no_detectado = ifelse(no_det, "TRUE", "FALSE"),
    CT = ifelse(no_det, NA_real_, x),
    stringsAsFactors = FALSE
  )

  es_bp1 <- df$TEJIDO == TEJIDO_EXCLUIDO
  e15 <- df[!es_bp1, , drop = FALSE]
  bp1 <- df[es_bp1, , drop = FALSE]

  stopifnot(all(!is.na(df$rsp29) & df$rsp29 != 40))
  stopifnot(identical(sort(unique(e15$TEJIDO)), sort(TEJIDOS_E15)))
  stopifnot(length(unique(e15$FETO)) == 36L)
  stopifnot(all(nzchar(e15$FETO)))
  stopifnot(nrow(e15) == 36L * length(TEJIDOS_E15) * length(GENES))

  list(e15 = e15, bp1 = bp1, diag = list(
    n_ct40 = sum(es_40), n_ctvacio = sum(es_vacio), n_ct_no_num = n_ct_no_num,
    n_bp1_filas = nrow(bp1), madres_bp1 = sort(unique(bp1$MADRE))
  ))
}

# ===========================================================================
# 2. ELISA: parte suero materno / liquido amniotico, censura a izquierda (D10)
# ===========================================================================
LOD_ELISA <- 0
TEJIDO_MS <- "Suero materno"
TEJIDO_LA <- "Líquido amniótico"

ingesta_elisa <- function() {
  g <- leer_grid(ARCHIVO_ELISA, HOJA_ELISA)
  col6 <- s_txt(g[[6]])
  fila_encab <- which(col6 == "TEJIDO")
  if (length(fila_encab) == 0L)
    stop("ELISA: no se encontro la fila de encabezado (celda 'TEJIDO')")
  fila_encab <- fila_encab[1]
  encab <- s_txt(unlist(g[fila_encab, ], use.names = FALSE))
  esperado <- c("TRATAMIENTO", "Abs 450", "Conc", "IL-6", "IDMADRE", "TEJIDO")
  for (k in seq_along(esperado)) {
    if (!startsWith(encab[k], esperado[k]))
      stop(sprintf("ELISA: columna %d = '%s', se esperaba ~'%s'",
                   k - 1L, encab[k], esperado[k]))
  }

  raw <- g[(fila_encab + 1L):nrow(g), , drop = FALSE]
  trat_raw <- s_txt(raw[[1]])
  abs450   <- to_num(raw[[2]])
  conc     <- to_num(raw[[3]])
  il6col   <- to_num(raw[[4]])
  idm      <- s_txt(raw[[5]])
  tejido   <- s_txt(raw[[6]])

  keep <- tejido %in% c(TEJIDO_MS, TEJIDO_LA)
  trat_raw <- trat_raw[keep]; abs450 <- abs450[keep]; conc <- conc[keep]
  il6col <- il6col[keep]; idm <- idm[keep]; tejido <- tejido[keep]

  bloque <- ifelse(tejido == TEJIDO_MS, "MS", "LA")
  tto <- vapply(strsplit(trat_raw, " ", fixed = TRUE),
                function(p) normalizar_tto(p[1]), character(1))
  sexo <- rep("", length(trat_raw))
  es_la <- bloque == "LA"
  sexo[es_la] <- ifelse(grepl(EMOJI_HEMBRA, trat_raw[es_la], fixed = TRUE), "HEMBRA",
                 ifelse(grepl(EMOJI_MACHO,  trat_raw[es_la], fixed = TRUE), "MACHO", ""))
  if (any(es_la & sexo == ""))
    stop("ELISA-LA: no se pudo leer el sexo de alguna fila")
  if (any(is.na(conc)))
    stop("ELISA: Conc no numerica en alguna fila")

  censurado <- conc < 0
  df <- data.frame(
    bloque = bloque, TRATAMIENTO_RAW = trat_raw, TTO = tto, SEXO = sexo,
    ID = idm, MADRE_ID = ifelse(bloque == "MS", idm, madre_id_de(idm)),
    TEJIDO = tejido, Abs450 = abs450, Conc = conc, IL6_col_crudo = il6col,
    censurado = ifelse(censurado, "TRUE", "FALSE"),
    IL6_pgml = ifelse(censurado, NA_real_, conc), LOD = LOD_ELISA,
    stringsAsFactors = FALSE
  )
  n <- list(MS = sum(df$bloque == "MS"), LA = sum(df$bloque == "LA"),
            MS_cens = sum(df$bloque == "MS" & censurado),
            LA_cens = sum(df$bloque == "LA" & censurado))
  stopifnot(n$MS > 0, n$LA > 0)
  list(filas = df, n = n)
}

# ===========================================================================
# 3. pSTAT3: normaliza SEXO/TTO, verifica balanceo de membranas (bloque D9)
# ===========================================================================
COLS_PSTAT3 <- c("MEMBRANA", "MADRE", "NOMINACION", "FETO", "GRUPO", "SEXO", "TTO",
                 "TEJIDO", "PSTAT3")

ingesta_pstat3 <- function() {
  g <- leer_grid(ARCHIVO_PSTAT3, HOJA_PSTAT3)
  encab <- s_txt(unlist(g[1, ], use.names = FALSE))
  if (!identical(encab[1:9], COLS_PSTAT3))
    stop(sprintf("pSTAT3: encabezado inesperado: %s", paste(encab[1:9], collapse = ", ")))

  raw <- g[-1, , drop = FALSE]
  cel <- lapply(1:9, function(j) s_txt(raw[[j]]))
  keep <- !Reduce(`&`, lapply(cel, function(v) v == ""))

  memb  <- to_num(raw[[1]])[keep]
  MADRE <- cel[[2]][keep]
  NOMIN <- cel[[3]][keep]
  FETO  <- cel[[4]][keep]
  GR_RAW <- cel[[5]][keep]
  SX_RAW <- cel[[6]][keep]
  TT_RAW <- cel[[7]][keep]
  TEJIDO <- cel[[8]][keep]
  pstat  <- to_num(raw[[9]])[keep]

  stopifnot(all(!is.na(memb) & memb == floor(memb)))
  stopifnot(all(!is.na(pstat)))
  SEXO <- vapply(SX_RAW, normalizar_sexo, character(1), USE.NAMES = FALSE)
  TTO  <- vapply(TT_RAW, normalizar_tto,  character(1), USE.NAMES = FALSE)

  df <- data.frame(
    MEMBRANA = as.integer(memb), MADRE = MADRE, MADRE_ID = madre_id_de(FETO),
    NOMINACION = NOMIN, FETO = FETO, GRUPO_RAW = GR_RAW, SEXO_RAW = SX_RAW,
    TTO_RAW = TT_RAW, SEXO = SEXO, TTO = TTO, GRUPO = grupo_norm(SEXO, TTO),
    TEJIDO = TEJIDO, PSTAT3 = pstat, stringsAsFactors = FALSE
  )
  membs <- sort(unique(df$MEMBRANA))
  conteo <- vapply(membs, function(m) sum(df$MEMBRANA == m), integer(1))
  balinstr <- paste(sprintf("%d:%d", membs, conteo), collapse = "|")
  balanceado <- length(unique(conteo)) == 1L && nrow(df) %% 4L == 0L
  list(filas = df, diag = list(membs = membs, conteo = conteo, balinstr = balinstr,
                               balanceado = balanceado, n = nrow(df)))
}

# ===========================================================================
# 4. Tabla de n real por archivo x grupo x sexo
# ===========================================================================
tabla_n <- function(qpcr_e15, elisa, pstat3) {
  filas <- list()
  add <- function(...) filas[[length(filas) + 1L]] <<- list(...)

  for (tej in TEJIDOS_E15) {
    sub <- qpcr_e15[qpcr_e15$TEJIDO == tej, , drop = FALSE]
    for (grp in GRUPOS_4) {
      sg <- sub[sub$GRUPO == grp, , drop = FALSE]
      sx <- sub("_.*$", "", grp); tt <- sub("^.*_", "", grp)
      add(ARCHIVO_QPCR, tej, "feto", sx, tt, grp,
          length(unique(sg$FETO)), nrow(sg))
    }
  }
  for (grp in GRUPOS_4) {
    sg <- pstat3[pstat3$GRUPO == grp, , drop = FALSE]
    sx <- sub("_.*$", "", grp); tt <- sub("^.*_", "", grp)
    add(ARCHIVO_PSTAT3, "", "feto", sx, tt, grp, length(unique(sg$FETO)), nrow(sg))
  }
  ms <- elisa[elisa$bloque == "MS", , drop = FALSE]
  for (tt in c("CONTROL", "LPS")) {
    sg <- ms[ms$TTO == tt, , drop = FALSE]
    add(ARCHIVO_ELISA, "MS (Suero materno)", "madre", "", tt, "",
        length(unique(sg$ID)), nrow(sg))
  }
  la <- elisa[elisa$bloque == "LA", , drop = FALSE]
  for (tt in c("CONTROL", "LPS")) for (sx in c("HEMBRA", "MACHO")) {
    sg <- la[la$TTO == tt & la$SEXO == sx, , drop = FALSE]
    add(ARCHIVO_ELISA, "LA (Liquido amniotico)", "saco", sx, tt,
        grupo_norm(sx, tt), length(unique(sg$ID)), nrow(sg))
  }

  orden_arch <- c(0L, 1L, 2L); names(orden_arch) <- c(ARCHIVO_QPCR, ARCHIVO_PSTAT3, ARCHIVO_ELISA)
  orden_bloque <- setNames(c(0L, 1L, 2L, 3L, 4L),
                           c("PLACENTA_E15", "BRAIN_E15", "",
                             "MS (Suero materno)", "LA (Liquido amniotico)"))
  k <- order(
    orden_arch[vapply(filas, `[[`, character(1), 1)],
    orden_bloque[vapply(filas, `[[`, character(1), 2)],
    ORDEN_TTO[vapply(filas, `[[`, character(1), 5)],
    (match(vapply(filas, `[[`, character(1), 4), c("","HEMBRA","MACHO")) - 1L)
  )
  list(header = c("archivo", "bloque", "unidad", "SEXO", "TTO", "GRUPO",
                  "n_unidades", "n_registros"),
       filas = filas[k])
}

# ===========================================================================
# 5. No-detectados de qPCR por gen x tejido x grupo
# ===========================================================================
tabla_no_detectados <- function(qpcr_e15) {
  filas <- list(); cero_calib <- character(0)
  add <- function(...) filas[[length(filas) + 1L]] <<- list(...)
  for (tej in TEJIDOS_E15) for (gen in GENES) for (grp in GRUPOS_4) {
    sub <- qpcr_e15[qpcr_e15$TEJIDO == tej & qpcr_e15$GEN == gen &
                      qpcr_e15$GRUPO == grp, , drop = FALSE]
    nt <- nrow(sub)
    nd <- sum(sub$no_detectado == "FALSE")
    add(tej, gen, grp, nt, nd, nt - nd, if (nt > 0) nd / nt else NA_real_)
    if (grp == "HEMBRA_CONTROL" && nt > 0 && nd == 0)
      cero_calib <- c(cero_calib, paste0(gen, "@", tej))
  }
  k <- order(
    ORDEN_TEJIDO[vapply(filas, `[[`, character(1), 1)],
    ORDEN_GEN[vapply(filas, `[[`, character(1), 2)],
    ORDEN_GRUPO[vapply(filas, `[[`, character(1), 3)]
  )
  list(header = c("TEJIDO", "GEN", "GRUPO", "n_total", "n_detectado",
                  "n_no_detectado", "prop_detectado"),
       filas = filas[k], cero_calib = cero_calib)
}

# ===========================================================================
# 6. Censura del ELISA por bloque x grupo (D10)
# ===========================================================================
tabla_censura <- function(elisa) {
  filas <- list()
  add <- function(...) filas[[length(filas) + 1L]] <<- list(...)
  ms <- elisa[elisa$bloque == "MS", , drop = FALSE]
  for (tt in c("CONTROL", "LPS")) {
    sg <- ms[ms$TTO == tt, , drop = FALSE]
    nc <- sum(sg$censurado == "TRUE")
    add("MS", tt, "", nrow(sg), nc, nrow(sg) - nc,
        if (nrow(sg) > 0) 100 * nc / nrow(sg) else NA_real_)
  }
  la <- elisa[elisa$bloque == "LA", , drop = FALSE]
  for (tt in c("CONTROL", "LPS")) for (sx in c("HEMBRA", "MACHO")) {
    sg <- la[la$TTO == tt & la$SEXO == sx, , drop = FALSE]
    nc <- sum(sg$censurado == "TRUE")
    add("LA", tt, sx, nrow(sg), nc, nrow(sg) - nc,
        if (nrow(sg) > 0) 100 * nc / nrow(sg) else NA_real_)
  }
  k <- order(
    ifelse(vapply(filas, `[[`, character(1), 1) == "MS", 0L, 1L),
    ORDEN_TTO[vapply(filas, `[[`, character(1), 2)],
    (match(vapply(filas, `[[`, character(1), 3), c("","HEMBRA","MACHO")) - 1L)
  )
  list(header = c("bloque", "TTO", "SEXO", "n", "n_censurado", "n_detectado",
                  "pct_censurado"),
       filas = filas[k])
}

# ===========================================================================
# 7. Que madres/fetos faltan en cada bloque respecto del diseno qPCR (Seccion 2.2)
# ===========================================================================
tabla_faltantes <- function(qpcr_e15, elisa, pstat3) {
  ref_madres <- sort(unique(qpcr_e15$MADRE_ID))
  ref_fetos  <- sort(unique(qpcr_e15$FETO))
  ms_ids     <- sort(unique(elisa$ID[elisa$bloque == "MS"]))
  la_ids     <- sort(unique(elisa$ID[elisa$bloque == "LA"]))
  la_madres  <- sort(unique(elisa$MADRE_ID[elisa$bloque == "LA"]))
  ps_fetos   <- sort(unique(pstat3$FETO))

  fila <- function(bloque, referencia, ref, dato) {
    list(bloque, referencia, length(ref), length(dato),
         length(intersect(ref, dato)),
         paste(sort(setdiff(ref, dato)), collapse = ";"),
         paste(sort(setdiff(dato, ref)), collapse = ";"))
  }
  filas <- list(
    fila("ELISA_MS", "qpcr_E15_madre_id", ref_madres, ms_ids),
    fila("ELISA_LA", "qpcr_E15_feto", ref_fetos, la_ids),
    fila("ELISA_LA", "qpcr_E15_madre_id", ref_madres, la_madres),
    fila("pSTAT3", "qpcr_E15_feto", ref_fetos, ps_fetos)
  )
  list(header = c("bloque", "referencia", "n_ref", "n_dato", "n_match",
                  "ids_ref_sin_dato", "ids_dato_sin_ref"),
       filas = filas)
}

# ===========================================================================
# 8. Resumen maquina-legible + reporte legible
# ===========================================================================
construir_resumen <- function(qd, ed, pd, cero_calib) {
  filas <- list()
  row <- function(chk, val, esp, ok)
    filas[[length(filas) + 1L]] <<- list(chk, .fmt(val), .fmt(esp),
                                         if (ok) "TRUE" else "FALSE")
  row("fuente_qpcr", fuente_datos(ARCHIVO_QPCR), "", TRUE)
  row("fuente_elisa", fuente_datos(ARCHIVO_ELISA), "", TRUE)
  row("fuente_pstat3", fuente_datos(ARCHIVO_PSTAT3), "", TRUE)
  row("qpcr_fetos_e15", qd$n_fetos_e15, 36, qd$n_fetos_e15 == 36)
  row("qpcr_filas_e15", qd$n_filas_e15, 720, qd$n_filas_e15 == 720)
  row("qpcr_tejidos_e15", paste(TEJIDOS_E15, collapse = ";"),
      paste(TEJIDOS_E15, collapse = ";"), TRUE)
  row("qpcr_brain_p1_filas_excluidas", qd$n_bp1_filas, ">0", qd$n_bp1_filas > 0)
  row("qpcr_ct_igual_40", qd$n_ct40, ">=0", TRUE)
  row("qpcr_ct_vacio", qd$n_ctvacio, ">=0", TRUE)
  row("qpcr_ct_no_detectado", qd$n_ct40 + qd$n_ctvacio, ">=0", TRUE)
  row("qpcr_ct_no_numerico", qd$n_ct_no_num, 0, qd$n_ct_no_num == 0)
  row("qpcr_genes_cero_det_en_calibrador",
      if (length(cero_calib)) paste(cero_calib, collapse = ";") else "(ninguno)", "", TRUE)
  row("elisa_ms_n", ed$MS, 14, ed$MS == 14)
  row("elisa_ms_censurado", ed$MS_cens, ">=0", TRUE)
  row("elisa_la_n", ed$LA, 28, ed$LA == 28)
  row("elisa_la_censurado", ed$LA_cens, ">=0", TRUE)
  row("elisa_columna_il6_usada", "no (ver analisis_descartados.md)", "no", TRUE)
  row("pstat3_filas", pd$n, 36, pd$n == 36)
  row("pstat3_membranas", pd$balinstr, "1:12|2:12|3:12", pd$balanceado)
  row("pstat3_balanceado", pd$balanceado, TRUE, pd$balanceado)
  list(header = c("check", "valor", "esperado", "ok"), filas = filas)
}

.trunc <- function(s, n = 60L) if (nchar(s) <= n) s else paste0(substr(s, 1, n - 1L), "…")

.md_tabla <- function(header, filas) {
  l1 <- paste0("| ", paste(header, collapse = " | "), " |")
  l2 <- paste0("| ", paste(rep("---", length(header)), collapse = " | "), " |")
  cuerpo <- vapply(filas, function(f)
    paste0("| ", paste(vapply(f, .fmt, character(1)), collapse = " | "), " |"),
    character(1))
  paste(c(l1, l2, cuerpo), collapse = "\n")
}

construir_reporte_md <- function(ntab, ndet, censt, faltt, rest, qd, ed, pd) {
  L <- c(
    "# Reporte de QC -- ingesta MIA-LPS", "",
    "Generado por `02_ingesta_qc` (R y Python producen este archivo identico).",
    sprintf("Fuentes en uso: qPCR=`%s`, ELISA=`%s`, pSTAT3=`%s`.",
            fuente_datos(ARCHIVO_QPCR), fuente_datos(ARCHIVO_ELISA),
            fuente_datos(ARCHIVO_PSTAT3)),
    "", "## 1. Saneamiento aplicado", "",
    sprintf("- qPCR `CT_CRUDO`: %d celdas `== 40` + %d vacias = %d no-detectados -> `NA` (D3, sin imputar).",
            qd$n_ct40, qd$n_ctvacio, qd$n_ct40 + qd$n_ctvacio),
    sprintf("- qPCR `BRAIN_P1`: %d filas excluidas del alcance E15 (madres %s); ver `analisis_descartados.md`.",
            qd$n_bp1_filas, paste(qd$madres_bp1, collapse = ", ")),
    sprintf("- ELISA: hoja partida en MS (n=%d) y LA (n=%d); `Conc < 0` -> censura a izquierda (D10), valor `NA`, `LOD = %g`.",
            ed$MS, ed$LA, LOD_ELISA),
    "- ELISA columna `IL-6` del crudo: NO usada (ver `analisis_descartados.md`).",
    sprintf("- pSTAT3: %d filas, membranas `%s` (%s).",
            pd$n, pd$balinstr, if (pd$balanceado) "balanceado" else "DESBALANCEADO"),
    "", "## 2. n real por archivo x grupo x sexo", "",
    .md_tabla(ntab$header, ntab$filas),
    "", "## 3. Censura del ELISA por bloque x grupo (D10)", "",
    .md_tabla(censt$header, censt$filas),
    "", "## 4. No-detectados de qPCR: casos con 0 detectados en el calibrador", ""
  )
  cero <- Filter(function(f) f[[3]] == "HEMBRA_CONTROL" && f[[4]] > 0 && f[[5]] == 0,
                 ndet$filas)
  if (length(cero)) {
    L <- c(L, "Genes x tejido no cuantificables por D7 (se excluyen del modelo en 05):",
           "", .md_tabla(ndet$header, cero))
  } else {
    L <- c(L, "Ninguno: todos los gen x tejido tienen >=1 deteccion en HEMBRA_CONTROL.")
  }
  L <- c(L, "", "(Tabla completa gen x tejido x grupo en `qc_no_detectados_qpcr.csv`.)",
         "", "## 5. Madres / fetos faltantes por bloque (respecto del diseno qPCR E15)", "")
  faltt_trunc <- lapply(faltt$filas, function(f)
    c(f[1:5], list(.trunc(f[[6]]), .trunc(f[[7]]))))
  L <- c(L, .md_tabla(faltt$header, faltt_trunc),
         "", "(Listas completas de IDs en `qc_faltantes.csv`.)",
         "", "## 6. Chequeos", "", .md_tabla(rest$header, rest$filas), "")
  paste(L, collapse = "\n")
}

# ===========================================================================
# 9. Artefactos compartidos: analisis_descartados / procedencia / verificaciones
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
  list(header = parse_csv_line(lineas[1]),
       filas = lapply(lineas[-1], parse_csv_line))
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
  todas <- todas[clave_orden(todas)]  # clave_orden usa method="radix" (== codepoint en Python)
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
"## 02_ingesta_qc",
"",
"### BRAIN_P1 excluido de la ingesta",
"",
"- **Que se probo / que trae el crudo:** el archivo `Raw data CTs.xlsx` incluye",
"  filas con `TEJIDO == BRAIN_P1` (cerebro, dia postnatal 1) ademas de",
"  `PLACENTA_E15` y `BRAIN_E15`.",
"- **Por que no se usa aca:** queda fuera del alcance de la Seccion 3 del brief",
"  (placenta y cerebro fetal E15). El P1 se analiza en un informe aparte.",
"- **Que se hizo en su lugar:** `02_ingesta_qc` filtra esas filas antes de",
"  cualquier calculo, deja constancia de cuantas eran y de que madres provienen",
"  (ver `qc_resumen.csv`, fila `qpcr_brain_p1_filas_excluidas`), y el resto del",
"  pipeline opera solo sobre `PLACENTA_E15` + `BRAIN_E15`.",
"",
"### Columna `IL-6` del ELISA: no se usa",
"",
"- **Que se probo / que trae el crudo:** la hoja `Sueros y LA` tiene una columna",
"  `IL-6` (indice 3) ademas de `Conc` (indice 2). En casi todas las filas es",
"  `max(Conc, 0)`; en el bloque de suero materno / LPS es `Conc x 4` (un factor",
"  de dilucion no documentado en la hoja).",
"- **Por que no sirve:** mezcla dos transformaciones (piso en 0 y factor de",
"  dilucion), pisa la censura a izquierda y no es una medida homogenea entre",
"  filas.",
"- **Que se hizo en su lugar:** se ignora la columna `IL-6`. El analisis del",
"  ELISA (03) usa `Conc` (indice 2) con el tratamiento de censura de D10:",
"  `Conc < 0` -> `censurado = TRUE`, valor `NA`, `LOD` registrado aparte.",
"",
"### Imputacion de no-detectados (recordatorio)",
"",
"- D3 ya fija que los no-detectados de qPCR (`CT_CRUDO == 40` o celda vacia) van a",
"  `NA` y **no se imputan por ningun metodo**. El detalle de por que se descarto la",
"  imputacion MNAR (`nondetects`) se documenta con numeros en T4/T10."
), collapse = "\n")

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  marca_ini <- "<!-- 02_ingesta_qc:inicio -->"
  marca_fin <- "<!-- 02_ingesta_qc:fin -->"
  nuevo <- paste0(marca_ini, "\n", BLOQUE_DESCARTES, "\n\n", marca_fin)
  if (file.exists(ruta)) {
    txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE)
    Encoding(txt) <- "UTF-8"
  } else {
    txt <- paste0("# Analisis descartados\n\n",
      "Que se probo, por que no funciono o no se uso, y que se hizo en su",
      " lugar. Una seccion por script.\n")
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
filas_de <- function(df, cols) {
  lapply(seq_len(nrow(df)), function(i) {
    r <- as.list(df[i, cols, drop = FALSE])
    names(r) <- NULL
    lapply(r, function(v) if (length(v) == 1L && is.na(v)) NA else v)
  })
}

main <- function() {
  q <- ingesta_qpcr(); qpcr_e15 <- q$e15; qpcr_bp1 <- q$bp1; qd <- q$diag
  e <- ingesta_elisa(); elisa <- e$filas; ed <- e$n
  p <- ingesta_pstat3(); pstat3 <- p$filas; pd <- p$diag

  qd$n_fetos_e15 <- length(unique(qpcr_e15$FETO))
  qd$n_filas_e15 <- nrow(qpcr_e15)

  # --- intermedios largos (data/processed/, regenerables) --------------
  q_ord <- order(qpcr_e15$MADRE_ID, qpcr_e15$FETO,
                 ORDEN_TEJIDO[qpcr_e15$TEJIDO], ORDEN_GEN[qpcr_e15$GEN])
  cols_e15 <- c("MADRE", "MADRE_ID", "NOMINACION", "FETO", "GRUPO_RAW", "SEXO_RAW",
                "TTO_RAW", "SEXO", "TTO", "GRUPO", "TEJIDO", "rsp29", "GEN",
                "CT_CRUDO", "no_detectado", "CT")
  escribir_tsv(file.path(RUTA_DATOS_PROC, "qpcr_e15_long.tsv"), cols_e15,
               filas_de(qpcr_e15[q_ord, , drop = FALSE], cols_e15))

  b_ord <- order(qpcr_bp1$MADRE, qpcr_bp1$FETO, qpcr_bp1$GEN)
  cols_bp1 <- c("MADRE", "NOMINACION", "FETO", "GRUPO_RAW", "SEXO", "TTO", "TEJIDO",
                "rsp29", "GEN", "CT_CRUDO", "no_detectado")
  escribir_tsv(file.path(RUTA_DATOS_PROC, "qpcr_brain_p1_excluido.tsv"), cols_bp1,
               filas_de(qpcr_bp1[b_ord, , drop = FALSE], cols_bp1))

  el_ord <- order(ifelse(elisa$bloque == "MS", 0L, 1L), ORDEN_TTO[elisa$TTO],
                  (match(elisa$SEXO, c("","HEMBRA","MACHO")) - 1L), elisa$ID)
  cols_el <- c("bloque", "TRATAMIENTO_RAW", "TTO", "SEXO", "ID", "MADRE_ID", "TEJIDO",
               "Abs450", "Conc", "IL6_col_crudo", "censurado", "IL6_pgml", "LOD")
  escribir_tsv(file.path(RUTA_DATOS_PROC, "elisa_long.tsv"), cols_el,
               filas_de(elisa[el_ord, , drop = FALSE], cols_el))

  p_ord <- order(pstat3$MEMBRANA, ORDEN_GRUPO[pstat3$GRUPO], pstat3$FETO)
  cols_ps <- c("MEMBRANA", "MADRE", "MADRE_ID", "NOMINACION", "FETO", "GRUPO_RAW",
               "SEXO_RAW", "TTO_RAW", "SEXO", "TTO", "GRUPO", "TEJIDO", "PSTAT3")
  escribir_tsv(file.path(RUTA_DATOS_PROC, "pstat3_long.tsv"), cols_ps,
               filas_de(pstat3[p_ord, , drop = FALSE], cols_ps))

  # --- tablas de QC (una copia por implementacion) --------------------
  ntab  <- tabla_n(qpcr_e15, elisa, pstat3)
  ndet  <- tabla_no_detectados(qpcr_e15)
  censt <- tabla_censura(elisa)
  faltt <- tabla_faltantes(qpcr_e15, elisa, pstat3)
  rest  <- construir_resumen(qd, ed, pd, ndet$cero_calib)

  for (base in c(RUTA_TABLAS_R, RUTA_TABLAS_PY)) {
    escribir_csv(file.path(base, "qc_n_por_grupo.csv"), ntab$header, ntab$filas)
    escribir_csv(file.path(base, "qc_no_detectados_qpcr.csv"), ndet$header, ndet$filas)
    escribir_csv(file.path(base, "qc_censura_elisa.csv"), censt$header, censt$filas)
    escribir_csv(file.path(base, "qc_faltantes.csv"), faltt$header, faltt$filas)
    escribir_csv(file.path(base, "qc_resumen.csv"), rest$header, rest$filas)
  }

  # --- reporte legible + artefactos compartidos ----------------------
  md <- construir_reporte_md(ntab, ndet, censt, faltt, rest, qd, ed, pd)
  escribir_lineas(file.path(RUTA_TABLAS, "qc_reporte.md"), md)
  actualizar_descartados()

  ent_q <- sprintf("data/%s/%s", fuente_datos(ARCHIVO_QPCR), ARCHIVO_QPCR)
  ent_e <- sprintf("data/%s/%s", fuente_datos(ARCHIVO_ELISA), ARCHIVO_ELISA)
  ent_p <- sprintf("data/%s/%s", fuente_datos(ARCHIVO_PSTAT3), ARCHIVO_PSTAT3)
  q3 <- paste(ent_q, ent_e, ent_p, sep = "; ")
  registrar_procedencia(list(
    list("data/processed/qpcr_e15_long.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_q, "qPCR largo E15, 40/vacio->NA (D3), BRAIN_P1 excluido"),
    list("data/processed/qpcr_brain_p1_excluido.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_q, "filas BRAIN_P1 apartadas del alcance E15"),
    list("data/processed/elisa_long.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_e, "ELISA largo MS+LA, censura a izquierda D10"),
    list("data/processed/pstat3_long.tsv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent_p, "pSTAT3 largo, SEXO/TTO normalizados"),
    list("outputs/tables/{R,python}/qc_n_por_grupo.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         q3, "n real por archivo x grupo x sexo"),
    list("outputs/tables/{R,python}/qc_no_detectados_qpcr.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent_q, "no-detectados por gen x tejido x grupo; marca D7"),
    list("outputs/tables/{R,python}/qc_censura_elisa.csv", "tabla", ESTE_SCRIPT,
         "PROPIO", ent_e, "% censura ELISA por bloque x grupo (D10)"),
    list("outputs/tables/{R,python}/qc_faltantes.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         q3, "IDs de referencia sin dato y viceversa por bloque"),
    list("outputs/tables/{R,python}/qc_resumen.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         q3, "chequeos maquina-legibles de la ingesta"),
    list("outputs/tables/qc_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO",
         q3, "reporte de QC legible")
  ))
  registrar_verificaciones(list(
    list("ingesta_qpcr_fetos_e15", "qPCR E15 tiene 36 fetos",
         .fmt(qd$n_fetos_e15), "36", if (qd$n_fetos_e15 == 36) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("ingesta_qpcr_filas_e15", "qPCR E15 tiene 36x2x10 = 720 filas",
         .fmt(qd$n_filas_e15), "720", if (qd$n_filas_e15 == 720) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("ingesta_qpcr_tejidos", "tejidos E15 = PLACENTA_E15;BRAIN_E15 (sin BRAIN_P1)",
         paste(sort(unique(qpcr_e15$TEJIDO)), collapse = ";"),
         paste(sort(TEJIDOS_E15), collapse = ";"),
         if (identical(sort(unique(qpcr_e15$TEJIDO)), sort(TEJIDOS_E15))) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("ingesta_ct_no_numerico", "no hay CT_CRUDO no numerico fuera de vacio",
         .fmt(qd$n_ct_no_num), "0", if (qd$n_ct_no_num == 0) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("ingesta_pstat3_balanceado", "pSTAT3: 3 membranas con igual n",
         pd$balinstr, "1:12|2:12|3:12", if (pd$balanceado) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("ingesta_elisa_bloques", "ELISA parte en MS y LA con n>0",
         sprintf("MS=%d;LA=%d", ed$MS, ed$LA), "MS>0;LA>0",
         if (ed$MS > 0 && ed$LA > 0) "TRUE" else "FALSE", ESTE_SCRIPT)
  ))

  # --- resumen por consola ------------------------------------------
  cat("== 02_ingesta_qc.R ==\n")
  cat(sprintf("  qPCR   fuente=%-9s E15: %d fetos / %d filas   BRAIN_P1 excluido: %d filas\n",
              fuente_datos(ARCHIVO_QPCR), qd$n_fetos_e15, qd$n_filas_e15, qd$n_bp1_filas))
  cat(sprintf("         CT no detectado: %d (==40) + %d (vacio) = %d\n",
              qd$n_ct40, qd$n_ctvacio, qd$n_ct40 + qd$n_ctvacio))
  cat(sprintf("         no cuantificable por D7: %s\n",
              if (length(ndet$cero_calib)) paste(ndet$cero_calib, collapse = ", ") else "(ninguno)"))
  cat(sprintf("  ELISA  fuente=%-9s MS n=%d (cens %d)   LA n=%d (cens %d)\n",
              fuente_datos(ARCHIVO_ELISA), ed$MS, ed$MS_cens, ed$LA, ed$LA_cens))
  cat(sprintf("  pSTAT3 fuente=%-9s %d filas   membranas %s   %s\n",
              fuente_datos(ARCHIVO_PSTAT3), pd$n, pd$balinstr,
              if (pd$balanceado) "balanceado" else "DESBALANCEADO"))
  cat("  -> data/processed/*.tsv, outputs/tables/{R,python}/qc_*.csv, outputs/tables/qc_reporte.md\n")
  cat("  invariantes OK\n")
}

if (sys.nframe() == 0L) main()
