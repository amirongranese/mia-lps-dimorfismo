# 01_generar_sinteticos.R -- Generador de datos sinteticos.
#
# Por que existe: el objetivo del proyecto es que un tercero, sin acceso a los
# datos crudos (ineditos, no versionados), pueda correr TODO el pipeline. Para
# eso data/synthetic/ contiene tres archivos que imitan la estructura EXACTA de
# los tres crudos -- mismos nombres de hoja y de columna, mismos tipos, y las
# mismas PATOLOGIAS que fuerzan las decisiones del brief:
#
#   * 18 camadas / 36 fetos E15 (una hembra y un macho por madre) + BRAIN_P1
#     como tejido extra que 02_ingesta_qc va a excluir.
#   * CT_CRUDO con no-detectados de dos formas (el valor literal 40 y la celda
#     vacia), igual que en el crudo. Ambas = "no detectado" (D3).
#   * il6 en BRAIN_E15 con 0 detectados en el calibrador (hembra-control):
#     dispara la regla D7 (gen no cuantificable -> solo proporcion de deteccion).
#   * ELISA con la hoja partida en suero materno / liquido amniotico y con Conc
#     negativos = censura a izquierda (D10).
#   * pSTAT3 con 3 membranas balanceadas (bloque fijo de D9).
#
# Los EFECTOS simulados son arbitrarios. La secuencia de numeros aleatorios se
# consume en el MISMO orden que python/01_generar_sinteticos.py (mismo RNG
# propio, misma semilla) -> los .tsv canonicos salen byte a byte identicos entre
# lenguajes (se verifica en data/synthetic/MANIFEST.tsv).

suppressPackageStartupMessages({
  library(openxlsx)   # escritura de .xlsx con control de tipo por celda
  library(digest)     # sha256 de los .tsv canonicos
})

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- file.path(getwd(), "R")
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

# ===========================================================================
# 1. Tabla canonica de fetos (deterministica, sin RNG)
#
# Los identificadores son estructura, no medicion: fijarlos sin tocar el stream
# aleatorio deja toda la secuencia para los valores medidos y hace trivial la
# paridad con Python.
# ===========================================================================
.FECHAS <- c("011025", "021025", "031025", "091025", "101025",
             "171025", "241025", "251025", "311025")

.tabla_fetos <- function() {
  filas <- list()
  for (par in list(c("C", "CONTROL"), c("L", "LPS"))) {
    prefijo <- par[1]; tto <- par[2]
    for (i in 1:9) {
      madre <- paste0(prefijo, i)
      madre_id <- paste0(prefijo, "_", .FECHAS[i], "_1")
      for (k in seq_along(NIVELES_SEXO)) {   # 1=HEMBRA, 2=MACHO
        sexo <- NIVELES_SEXO[k]
        nn <- 3 + 2 * (i - 1) + (k - 1)      # distinto dentro de la madre
        emo <- if (sexo == "HEMBRA") EMOJI_HEMBRA else EMOJI_MACHO
        filas[[length(filas) + 1L]] <- list(
          MADRE = madre, MADRE_ID = madre_id,
          FETO = sprintf("%s.%02d", madre_id, nn),
          NOMINACION = paste0(madre, emo, VS16),
          SEXO = sexo, TTO = tto,
          GRUPO = etiqueta_grupo(sexo, tto)
        )
      }
    }
  }
  filas
}

FETOS <- .tabla_fetos()   # 36 fetos: C1-H, C1-M, C2-H, ... L9-M
.clave_fs <- function(madre, sexo) paste(madre, sexo, sep = "|")
FETO_POR_MADRE_SEXO <- setNames(FETOS,
  vapply(FETOS, function(f) .clave_fs(f$MADRE, f$SEXO), character(1)))

# Madres con tejido BRAIN_P1 (subconjunto): C1..C7 y L1..L7. Tres (C6, C7, L7)
# con la celda FETO vacia -> 22 FETO unicos + 6 vacios = 28 fetos * 10 = 280
# filas, igual que el crudo.
.MADRES_P1 <- c(paste0("C", 1:7), paste0("L", 1:7))
.MADRES_P1_FETO_VACIO <- c("C6", "C7", "L7")

.GENES_CON_NA <- c("fatcd36", "glut1", "glut3", "slc38a1", "slc38a2")
.N_NA_POR_GEN <- 6L

# ===========================================================================
# 2. qPCR: tabla larga de CT
#
# Ct plausible por gen x tejido: la via IL-6 (il6, il6R) alta -> muchos
# no-detectados, il6 en cerebro practicamente ausente; transportadores en rango
# detectable con pocos no-detectados. Un efecto SEXO*TTO minusculo en un par de
# genes evita el ruido puro; su tamano es arbitrario.
# ===========================================================================
.CT_BASE <- c(fatcd36 = 27.0, fatp1 = 32.0, fatp4 = 28.0, glut1 = 26.0,
              glut3 = 31.5, slc38a1 = 33.0, slc38a2 = 29.5, gp130 = 25.5,
              il6R = 36.5, il6 = 37.5)
.CT_SHIFT_TEJIDO <- c(PLACENTA_E15 = 0.0, BRAIN_E15 = 2.0, BRAIN_P1 = 1.5)

.media_ct <- function(gen, tejido, sexo, tto) {
  mu <- .CT_BASE[[gen]] + .CT_SHIFT_TEJIDO[[tejido]]
  if (gen == "il6") {
    mu <- c(PLACENTA_E15 = 35.5, BRAIN_E15 = 40.7, BRAIN_P1 = 41.0)[[tejido]]
  }
  if (gen == "il6R" && tejido != "PLACENTA_E15") mu <- mu + 2.5
  if (gen %in% c("glut1", "slc38a2") && sexo == "MACHO" && tto == "LPS") mu <- mu - 0.6
  mu
}

.sigma_ct <- function(gen, tejido) {
  s <- if (gen %in% c("il6", "il6R")) 1.3 else 0.8
  if (tejido != "PLACENTA_E15") s <- s + 0.2
  s
}

.generar_qpcr <- function(rng) {
  n_max <- 1000L
  MADRE <- character(n_max); NOMINACION <- character(n_max); FETO <- character(n_max)
  GRUPO <- character(n_max); SEXO <- character(n_max); TTO <- character(n_max)
  TEJIDO <- character(n_max); rsp29 <- numeric(n_max); GEN <- character(n_max)
  lat <- numeric(n_max)
  p <- 0L

  bloque <- function(feto, tejido, feto_field) {
    rsp <- min(max(rng$norm1(25.0, 2.3), 19.5), 31.9)
    for (gen in GENES) {
      l <- .media_ct(gen, tejido, feto$SEXO, feto$TTO) +
        rng$norm1(0.0, .sigma_ct(gen, tejido))
      p <<- p + 1L
      MADRE[p] <<- feto$MADRE; NOMINACION[p] <<- feto$NOMINACION
      FETO[p] <<- feto_field; GRUPO[p] <<- feto$GRUPO; SEXO[p] <<- feto$SEXO
      TTO[p] <<- feto$TTO; TEJIDO[p] <<- tejido; rsp29[p] <<- rsp
      GEN[p] <<- gen; lat[p] <<- l
    }
  }

  for (feto in FETOS) for (tejido in TEJIDOS_E15) bloque(feto, tejido, feto$FETO)
  for (feto in FETOS) {
    if (!(feto$MADRE %in% .MADRES_P1)) next
    ff <- if (feto$MADRE %in% .MADRES_P1_FETO_VACIO) "" else feto$FETO
    bloque(feto, TEJIDO_EXCLUIDO, ff)
  }

  df <- data.frame(
    MADRE = MADRE[1:p], NOMINACION = NOMINACION[1:p], FETO = FETO[1:p],
    GRUPO = GRUPO[1:p], SEXO = SEXO[1:p], TTO = TTO[1:p], TEJIDO = TEJIDO[1:p],
    rsp29 = rsp29[1:p], GEN = GEN[1:p], stringsAsFactors = FALSE
  )
  lat <- lat[1:p]

  # --- Deteccion: latente >= 40 -> no detectado -----------------------
  ct <- ifelse(lat >= 40.0, 40.0, round(pmin(pmax(lat, 21.0), 39.9), 2))

  # --- D7: il6 / BRAIN_E15 sin detectados en el calibrador ------------
  # En el crudo il6 en cerebro tiene unos pocos detectados fuera del calibrador
  # pero 0/9 en hembra-control. Se deja lo estocastico y se fuerza el calibrador
  # a no-detectado (condicion exacta que dispara D7).
  m_cal <- df$GEN == "il6" & df$TEJIDO == "BRAIN_E15" &
    df$SEXO == "HEMBRA" & df$TTO == "CONTROL"
  ct[m_cal] <- 40.0

  df$CT_CRUDO <- ct

  # --- Celdas vacias: 6 por gen en 5 genes = 30 NA -------------------
  for (gen in .GENES_CON_NA) {
    cand <- which(df$GEN == gen)                     # posiciones ascendentes
    elegidas <- rng$elegir_sin_reemplazo(cand, .N_NA_POR_GEN)
    df$CT_CRUDO[elegidas] <- NA_real_
  }

  df[, c("MADRE", "NOMINACION", "FETO", "GRUPO", "SEXO", "TTO", "TEJIDO",
         "rsp29", "GEN", "CT_CRUDO")]
}

# ===========================================================================
# 3. ELISA IL-6: hoja "Sueros y LA" (44 x 14) + "CURVA IL6"
#
# El crudo es una planilla con un rotulo suelto en la fila 0, encabezados en la
# fila 1, dos bloques apilados (suero materno / liquido amniotico) y celdas de
# scratch a la derecha. 02 la parte por la columna 5. Relacion Abs<->Conc real:
# los coeficientes 0.0016 y 0.1057 son las celdas M/OO del crudo, y
#   Conc = (Abs - 0.1057) / 0.0016 ,  IL-6 = max(Conc, 0)  [columna a NO usar].
# Conc negativa = censura a izquierda (D10).
# ===========================================================================
.ELISA_SLOPE <- 0.0016
.ELISA_BLANK_ABS <- 0.1057

.abs_conc_il6 <- function(conc_objetivo) {
  abs450 <- round(conc_objetivo * .ELISA_SLOPE + .ELISA_BLANK_ABS, 3)
  conc <- (abs450 - .ELISA_BLANK_ABS) / .ELISA_SLOPE
  c(abs450 = abs450, conc = conc, il6 = max(conc, 0.0))
}

.generar_elisa <- function(rng) {
  ms <- list()
  c15 <- paste0("C", 1:5)
  for (j in seq_along(c15)) {
    madre_id <- FETO_POR_MADRE_SEXO[[.clave_fs(c15[j], "HEMBRA")]]$MADRE_ID
    conc_obj <- if (j <= 4) -(5.0 + abs(rng$norm1(0.0, 15.0))) else
      30.0 + abs(rng$norm1(0.0, 20.0))
    v <- .abs_conc_il6(conc_obj)
    ms[[length(ms) + 1L]] <- list(col0 = "Control", abs = v[["abs450"]],
      conc = v[["conc"]], il6 = v[["il6"]], id = madre_id, tej = "Suero materno")
  }
  for (madre in paste0("L", 1:9)) {
    madre_id <- FETO_POR_MADRE_SEXO[[.clave_fs(madre, "HEMBRA")]]$MADRE_ID
    conc_obj <- 150.0 + abs(rng$norm1(650.0, 380.0))
    v <- .abs_conc_il6(conc_obj)
    ms[[length(ms) + 1L]] <- list(col0 = "LPS", abs = v[["abs450"]],
      conc = v[["conc"]], il6 = v[["il6"]], id = madre_id, tej = "Suero materno")
  }

  la <- list()
  bloque_la <- function(madres, sexo, tto, media_baja) {
    emo <- if (sexo == "MACHO") EMOJI_MACHO else EMOJI_HEMBRA
    etq <- paste0(if (tto == "CONTROL") "Control" else "LPS", " ", emo)
    for (madre in madres) {
      feto <- FETO_POR_MADRE_SEXO[[.clave_fs(madre, sexo)]]
      if (tto == "LPS") {
        pp <- rng$unif1()
        conc_obj <- if (pp < 0.25) 150.0 + abs(rng$norm1(300.0, 250.0)) else
          rng$norm1(media_baja, 22.0)
      } else {
        conc_obj <- rng$norm1(media_baja, 10.0)
      }
      v <- .abs_conc_il6(conc_obj)
      la[[length(la) + 1L]] <<- list(col0 = etq, abs = v[["abs450"]],
        conc = v[["conc"]], il6 = v[["il6"]], id = feto$FETO,
        tej = "Líquido amniótico")
    }
  }
  bloque_la(c15, "MACHO", "CONTROL", 7.0)
  bloque_la(c15, "HEMBRA", "CONTROL", -9.0)
  bloque_la(paste0("L", 1:9), "MACHO", "LPS", 1.0)
  bloque_la(paste0("L", 1:9), "HEMBRA", "LPS", 0.0)
  list(ms = ms, la = la)
}

.CURVA_PGML <- c(0.0, 15.6, 31.25, 62.5, 125.0, 250.0, 500.0, 1000.0)

.generar_curva <- function(rng) {
  filas <- vector("list", length(.CURVA_PGML))
  for (i in seq_along(.CURVA_PGML)) {
    pg <- .CURVA_PGML[i]
    a1 <- round(pg * .ELISA_SLOPE + 0.0695 + rng$norm1(0.0, 0.004), 4)
    a2 <- round(pg * .ELISA_SLOPE + 0.0695 + rng$norm1(0.0, 0.004), 4)
    # media sin redondear: evita el caso "mitad exacta" en el 4o decimal, que
    # round() resuelve distinto en R y en Python. Es una celda de scratch.
    filas[[i]] <- c(a1, a2, (a1 + a2) / 2, pg)
  }
  filas
}

# ===========================================================================
# 4. pSTAT3: 36 filas, 3 membranas balanceadas (bloque fijo de D9)
# ===========================================================================
.MEMBRANAS <- list(
  c("C1", "C2", "C3", "L1", "L2", "L3"),
  c("C4", "C5", "C6", "L4", "L5", "L6"),
  c("C7", "C8", "C9", "L7", "L8", "L9")
)

.generar_pstat3 <- function(rng) {
  filas <- list()
  for (m in seq_along(.MEMBRANAS)) {
    madres_m <- .MEMBRANAS[[m]]
    controles <- madres_m[startsWith(madres_m, "C")]
    lpss <- madres_m[startsWith(madres_m, "L")]
    for (sexo in NIVELES_SEXO) {
      for (par in list(list("CONTROL", controles), list("LPS", lpss))) {
        tto <- par[[1]]; grp_madres <- par[[2]]
        mu <- if (tto == "CONTROL") 2.2 else 6.0
        for (madre in grp_madres) {
          feto <- FETO_POR_MADRE_SEXO[[.clave_fs(madre, sexo)]]
          val <- min(max(rng$norm1(mu, 0.65), 0.85), 7.2)
          filas[[length(filas) + 1L]] <- list(
            MEMBRANA = m, MADRE = madre, NOMINACION = feto$NOMINACION,
            FETO = feto$FETO, GRUPO = etiqueta_grupo(sexo, tto),
            SEXO = if (sexo == "HEMBRA") "Hembra" else "Macho",
            TTO = if (tto == "CONTROL") "Control" else "LPS",
            TEJIDO = "PLACENTA_E15", PSTAT3 = val
          )
        }
      }
    }
  }
  do.call(rbind.data.frame, c(lapply(filas, function(r)
    data.frame(r, stringsAsFactors = FALSE)), make.row.names = FALSE))
}

# ===========================================================================
# 5. Escritura: .xlsx (para el pipeline) + .tsv canonico (paridad R/Python)
# ===========================================================================
.fmt <- function(x) {
  if (is.null(x) || (length(x) == 1L && is.na(x))) return("")
  if (is.character(x)) return(x)
  if (is.logical(x)) return(if (x) "TRUE" else "FALSE")
  # numerico: entero exacto -> sin decimales; si no, %.10g (igual que Python)
  if (x == floor(x) && abs(x) < 1e15) return(sprintf("%d", as.integer(round(x))))
  sprintf("%.10g", x)
}

.escribir_tsv <- function(ruta, filas, encabezado) {
  lineas <- character(0)
  if (!is.null(encabezado)) lineas <- paste(encabezado, collapse = "\t")
  cuerpo <- vapply(filas, function(fila)
    paste(vapply(fila, .fmt, character(1)), collapse = "\t"), character(1))
  texto <- paste0(paste(c(lineas, cuerpo), collapse = "\n"), "\n")
  con <- file(ruta, open = "wb")
  writeBin(charToRaw(enc2utf8(texto)), con)
  close(con)
  digest(charToRaw(enc2utf8(texto)), algo = "sha256", serialize = FALSE)
}

.df_a_filas <- function(df) {
  lapply(seq_len(nrow(df)), function(i) as.list(df[i, , drop = FALSE]))
}

main <- function() {
  rng <- nuevo_rng(SEMILLA)

  # --- orden de consumo del RNG: qpcr -> elisa -> curva -> pstat3 -------
  qpcr <- .generar_qpcr(rng)
  el <- .generar_elisa(rng); ms <- el$ms; la <- el$la
  curva <- .generar_curva(rng)
  pstat3 <- .generar_pstat3(rng)

  dst <- RUTA_DATOS_SINT
  if (!dir.exists(dst)) dir.create(dst, recursive = TRUE)

  # --- qPCR .xlsx + .tsv -------------------------------------------
  wb <- createWorkbook()
  addWorksheet(wb, HOJA_QPCR)
  writeData(wb, HOJA_QPCR, qpcr, keepNA = FALSE)
  saveWorkbook(wb, file.path(dst, ARCHIVO_QPCR), overwrite = TRUE)
  h_qpcr <- .escribir_tsv(file.path(dst, "Raw data CTs.tsv"),
                          .df_a_filas(qpcr), colnames(qpcr))

  # --- pSTAT3 .xlsx + .tsv ----------------------------------------
  wb <- createWorkbook()
  addWorksheet(wb, HOJA_PSTAT3)
  writeData(wb, HOJA_PSTAT3, pstat3, keepNA = FALSE)
  saveWorkbook(wb, file.path(dst, ARCHIVO_PSTAT3), overwrite = TRUE)
  h_pstat3 <- .escribir_tsv(file.path(dst, "pstat3 placenta.tsv"),
                            .df_a_filas(pstat3), colnames(pstat3))

  # --- ELISA .xlsx (2 hojas tipo grilla, sin encabezado real) --------
  ancho <- 14L
  filas_datos <- c(ms, la)
  g_elisa <- vector("list", 2L + length(filas_datos))
  g_elisa[[1]] <- c(list("IL6"), as.list(rep("", ancho - 1L)))
  g_elisa[[2]] <- list("TRATAMIENTO y sexo", "Abs 450 MS", "Conc ", "IL-6 ",
    "IDMADRE/FETO", "TEJIDO", "", "", "", "", "", "Abs", "m", "oo")
  for (k in seq_along(filas_datos)) {
    fd <- filas_datos[[k]]
    fila <- c(list(fd$col0, fd$abs, fd$conc, fd$il6, fd$id, fd$tej),
              as.list(rep("", ancho - 6L)))
    if (k == 1L) { fila[[13]] <- .ELISA_SLOPE; fila[[14]] <- .ELISA_BLANK_ABS }
    g_elisa[[2L + k]] <- fila
  }

  alto_curva <- 22L
  g_curva <- lapply(seq_len(alto_curva), function(i) as.list(rep("", 10L)))
  g_curva[[14]] <- list("", "", "Abs 450", "pg/ml", "", "", "", "", "", "")
  for (r in seq_along(curva)) {
    fila <- curva[[r]]
    g_curva[[14L + r]] <- c(as.list(fila), as.list(rep("", 10L - length(fila))))
  }

  wb <- createWorkbook()
  # orden de hojas como en el crudo: CURVA IL6 primero
  addWorksheet(wb, HOJA_CURVA)
  writeData(wb, HOJA_CURVA,
    x = as.data.frame(list("", "", "Abs 450", "pg/ml"), stringsAsFactors = FALSE),
    startCol = 1, startRow = 14, colNames = FALSE)
  df_curva <- as.data.frame(do.call(rbind, curva))
  names(df_curva) <- c("a1", "a2", "media", "pg")
  writeData(wb, HOJA_CURVA, df_curva, startCol = 1, startRow = 15, colNames = FALSE)

  addWorksheet(wb, HOJA_ELISA)
  # bloques tipados: banner, encabezado, scratch, y datos numericos
  writeData(wb, HOJA_ELISA, "IL6", startCol = 1, startRow = 1, colNames = FALSE)
  writeData(wb, HOJA_ELISA,
    x = as.data.frame(as.list(unlist(g_elisa[[2]])), stringsAsFactors = FALSE),
    startCol = 1, startRow = 2, colNames = FALSE)
  df_elisa <- data.frame(
    col0 = vapply(filas_datos, function(x) x$col0, character(1)),
    abs  = vapply(filas_datos, function(x) x$abs, numeric(1)),
    conc = vapply(filas_datos, function(x) x$conc, numeric(1)),
    il6  = vapply(filas_datos, function(x) x$il6, numeric(1)),
    id   = vapply(filas_datos, function(x) x$id, character(1)),
    tej  = vapply(filas_datos, function(x) x$tej, character(1)),
    stringsAsFactors = FALSE
  )
  writeData(wb, HOJA_ELISA, df_elisa, startCol = 1, startRow = 3, colNames = FALSE)
  writeData(wb, HOJA_ELISA, .ELISA_SLOPE, startCol = 13, startRow = 3, colNames = FALSE)
  writeData(wb, HOJA_ELISA, .ELISA_BLANK_ABS, startCol = 14, startRow = 3, colNames = FALSE)

  saveWorkbook(wb, file.path(dst, ARCHIVO_ELISA), overwrite = TRUE)

  h_elisa <- .escribir_tsv(
    file.path(dst, "ELISA IL6 2026 Dosis 100__Sueros y LA.tsv"), g_elisa, NULL)
  h_curva <- .escribir_tsv(
    file.path(dst, "ELISA IL6 2026 Dosis 100__CURVA IL6.tsv"), g_curva, NULL)

  # --- MANIFEST.tsv: derivado solo del contenido -> byte-identico R/Py --
  manifest <- list(
    list(ARCHIVO_QPCR, HOJA_QPCR, nrow(qpcr), ncol(qpcr), h_qpcr),
    list(ARCHIVO_PSTAT3, HOJA_PSTAT3, nrow(pstat3), ncol(pstat3), h_pstat3),
    list(ARCHIVO_ELISA, HOJA_ELISA, length(g_elisa), ancho, h_elisa),
    list(ARCHIVO_ELISA, HOJA_CURVA, alto_curva, 10L, h_curva)
  )
  .escribir_tsv(file.path(dst, "MANIFEST.tsv"), manifest,
                c("archivo", "hoja", "n_filas", "n_columnas", "sha256"))

  # --- PATOLOGIAS.tsv: contadores que 02+ dan por sentados ------------
  es_e15 <- qpcr$TEJIDO %in% TEJIDOS_E15
  ct40 <- sum(qpcr$CT_CRUDO == 40.0, na.rm = TRUE)
  ctna <- sum(is.na(qpcr$CT_CRUDO))
  m_il6 <- qpcr$GEN == "il6" & qpcr$TEJIDO == "BRAIN_E15"
  il6_det <- sum(m_il6 & !is.na(qpcr$CT_CRUDO) & qpcr$CT_CRUDO < 40)
  m_il6_cal <- m_il6 & qpcr$GRUPO == etiqueta_grupo("HEMBRA", "CONTROL")
  il6_cal_det <- sum(m_il6_cal & !is.na(qpcr$CT_CRUDO) & qpcr$CT_CRUDO < 40)
  ms_neg <- sum(vapply(ms, function(x) x$conc < 0, logical(1)))
  la_neg <- sum(vapply(la, function(x) x$conc < 0, logical(1)))
  memb <- paste(vapply(1:3, function(m)
    sprintf("%d:%d", m, sum(pstat3$MEMBRANA == m)), character(1)), collapse = "|")
  pat <- list(
    list("qpcr_n_filas", nrow(qpcr)),
    list("qpcr_n_madres", length(unique(qpcr$MADRE))),
    list("qpcr_n_fetos_e15", length(unique(qpcr$FETO[es_e15]))),
    list("qpcr_n_filas_e15", sum(es_e15)),
    list("qpcr_ct_igual_40", ct40),
    list("qpcr_ct_vacio", ctna),
    list("qpcr_il6_brain_e15_detectados", il6_det),
    list("qpcr_il6_brain_e15_hembracontrol_detectados", il6_cal_det),
    list("elisa_ms_n", length(ms)),
    list("elisa_ms_conc_neg", ms_neg),
    list("elisa_la_n", length(la)),
    list("elisa_la_conc_neg", la_neg),
    list("pstat3_n_filas", nrow(pstat3)),
    list("pstat3_membranas", memb)
  )
  .escribir_tsv(file.path(dst, "PATOLOGIAS.tsv"), pat, c("patologia", "valor"))

  # --- invariantes duros ------------------------------------------
  stopifnot(
    nrow(qpcr) == 1000L,
    length(unique(qpcr$MADRE)) == 18L,
    length(unique(qpcr$FETO[es_e15])) == 36L,
    sum(es_e15) == 720L,
    ctna == 30L,
    ct40 > 0L,
    il6_cal_det == 0L,
    il6_det >= 1L,
    ms_neg > 0L, la_neg > 0L,
    nrow(pstat3) == 36L,
    all(vapply(1:3, function(m) sum(pstat3$MEMBRANA == m) == 12L, logical(1)))
  )

  cat("== 01_generar_sinteticos.R ==\n")
  cat(sprintf("destino: %s\n", dst))
  cat(sprintf("  %-32s %4d filas  sha256(tsv)=%s\n", ARCHIVO_QPCR, nrow(qpcr),
              substr(h_qpcr, 1, 12)))
  cat(sprintf("  %-32s %4d filas  sha256(tsv)=%s\n", ARCHIVO_PSTAT3, nrow(pstat3),
              substr(h_pstat3, 1, 12)))
  cat(sprintf("  %-32s %4d+%d filas  sha256(tsv)=%s / %s\n", ARCHIVO_ELISA,
              length(g_elisa), alto_curva, substr(h_elisa, 1, 12),
              substr(h_curva, 1, 12)))
  cat(sprintf("  CT==40: %d   CT vacio: %d   il6/BRAIN_E15 detectados: %d\n",
              ct40, ctna, il6_det))
  cat(sprintf("  ELISA Conc<0: MS=%d  LA=%d\n", ms_neg, la_neg))
  cat("  invariantes OK\n")
}

if (sys.nframe() == 0L) {
  main()
}
