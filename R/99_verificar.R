# 99_verificar.R -- T11: verificacion final del reanalisis reproducible MIA-LPS.
#
# Por que existe este archivo: es la compuerta del proyecto. Partiendo de un
# clon limpio, `.\run_all.ps1` corre todo y este script comprueba, sin acceso a
# los datos crudos, que:
#   1. existe y no esta vacio cada item del checklist de AGENTS.md 1;
#   2. las salidas numericas R y Python concuerdan celda a celda dentro de la
#      tolerancia declarada (re-corre la comparacion de 98_comparacion en
#      proceso, sobre outputs/tables/{R,python}/*.csv);
#   3. los `.md` de copia unica de outputs/tables/ y la proyeccion sin figuras
#      de docs/informe.html son BYTE-IDENTICOS entre la corrida R y la Python
#      (via los snapshots que deja 12_informe en outputs/intermediate/render/);
#   4. ninguna fila de verificaciones.csv quedo distinta de TRUE.
# Escribe logs/corrida_<fecha>.txt (fecha, entorno, versiones, semilla, fuente,
# concordancia, cuantas verificaciones pasaron) e imprime
# `TODAS LAS VERIFICACIONES PASARON` si y solo si todo lo duro pasa Y nada quedo
# sin ejecutar.
#
# UNA SOLA IMPLEMENTACION (punto 4.3 de pedidos/cambios_revision_codex.md): quien
# solo tenga R (o solo Python) tiene que poder verificar su mitad. Con
# `run_all.ps1 -Only R|python` no corre 98_comparacion, asi que los chequeos 2 y 3
# (cruce R<->Python) no se pueden ejecutar: quedan en un tercer estado,
# NO_EJECUTADA -- ni pase ni fallo, visible en la salida y en el log. En ese caso
# la salida final dice VERIFICACION PARCIAL, nunca "TODAS LAS VERIFICACIONES
# PASARON". El modo lo informa la variable de entorno MIA_LPS_UNICA_IMPL, que
# escribe run_all.ps1; corriendo 99 a mano, la ausencia de la contraparte en
# outputs/tables/ tiene el mismo efecto.
#
# NO hace analisis. Solo lee salidas. `docs/informe.pdf` es blando: si 12_informe
# no encontro motor de PDF, su ausencia no frena la verificacion (AGENTS 8).
#
# PARIDAD: los helpers de lectura/comparacion son copia de 98_comparacion (la
# convencion del repo es que cada script importe solo 00_config).

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

LANG <- "R"
TOL_EST <- TOL_ESTADISTICO
TOL_P <- TOL_P_ITERATIVO

# Corrida de una sola implementacion: "R", "python" o "" (corrida completa).
UNICA_IMPL <- Sys.getenv("MIA_LPS_UNICA_IMPL", "")

# --- Checklist de AGENTS.md 1 ------------------------------------------
SCRIPTS_NUM <- c(
  "00_config", "01_generar_sinteticos", "02_ingesta_qc", "03_elisa",
  "04_qpcr_cuantificacion", "05_qpcr_modelos", "06_pstat3", "07_figuras_acto1",
  "08_acto2_correlaciones", "09_acto2_dispersion", "10_acto2_simulacion",
  "11_sensibilidad", "98_comparacion", "12_informe", "99_verificar"
)
FIG_ACTO1 <- c(
  "acto1_expresion_PLACENTA_E15.png", "acto1_expresion_BRAIN_E15.png",
  "acto1_pstat3.png", "acto1_elisa_ms.png", "acto1_elisa_la.png"
)
FIG_ACTO2 <- c(
  "acto2_dispersion_placenta_cerebro.png",
  "acto2_coexpresion_SPLOM_PLACENTA_E15.png",
  "acto2_coexpresion_SPLOM_BRAIN_E15.png",
  "acto2_dispersion_sd.png", "acto2_test_delta_rho.png",
  "acto2_simulacion_delta_rho.png",
  "acto2_sensibilidad_eigengene.png", "acto2_sensibilidad_excl_extremo.png"
)
SINTETICOS <- c(
  "Raw data CTs.tsv", "pstat3 placenta.tsv",
  "ELISA IL6 2026 Dosis 100__Sueros y LA.tsv",
  "ELISA IL6 2026 Dosis 100__CURVA IL6.tsv", "MANIFEST.tsv"
)
REPORTES_MD <- c(
  "qc_reporte.md", "elisa_reporte.md", "qpcr_cuantificacion_reporte.md",
  "qpcr_modelos_reporte.md", "pstat3_reporte.md",
  "acto2_correlaciones_reporte.md", "acto2_dispersion_reporte.md",
  "acto2_simulacion_reporte.md", "acto2_sensibilidad_reporte.md",
  "comparacion_reporte.md", "analisis_descartados.md"
)
R_PAQUETES <- c("car", "ARTool", "emmeans", "survival", "sandwich", "lmtest",
                "nlme", "ggplot2", "readxl", "writexl", "NADA", "GGally",
                "patchwork", "digest")

# =========================================================================
# Lectura de CSV -- copia de 98_comparacion.
# =========================================================================
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
  if (!file.exists(ruta)) return(list(header = NULL, filas = list()))
  txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(txt) <- "UTF-8"
  lineas <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  if (length(lineas) && lineas[length(lineas)] == "") lineas <- lineas[-length(lineas)]
  if (!length(lineas)) return(list(header = NULL, filas = list()))
  list(header = parse_csv_line(lineas[1]), filas = lapply(lineas[-1], parse_csv_line))
}
leer_bytes <- function(ruta) {
  if (!file.exists(ruta)) return(NULL)
  readBin(ruta, "raw", n = file.info(ruta)$size)
}

# =========================================================================
# Comparacion numerica celda a celda -- copia de 98_comparacion.
# =========================================================================
.es_num <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NA_real_)
  v <- suppressWarnings(as.numeric(s))
  if (is.na(v) || !is.finite(v)) return(NA_real_)
  v
}
.dentro <- function(a, b, tol) {
  d <- abs(a - b)
  if (d <= tol) return(TRUE)
  m <- max(abs(a), abs(b))
  m > 0 && d / m <= tol
}
comparar_archivo <- function(nombre, ruta_r, ruta_py) {
  r <- leer_csv(ruta_r); p <- leer_csv(ruta_py)
  hr <- r$header; hp <- p$header; fr <- r$filas; fp <- p$filas
  br <- leer_bytes(ruta_r); bp <- leer_bytes(ruta_py)
  byte_id <- !is.null(br) && !is.null(bp) && identical(br, bp)
  header_igual <- identical(hr, hp)
  filas_ok <- length(fr) == length(fp)
  n_num <- 0L; n_dif_txt <- 0L; n_fuera <- 0L; max_abs <- 0; peor <- ""
  if (header_igual && filas_ok && !is.null(hr)) {
    for (i in seq_along(fr)) {
      rr <- fr[[i]]; rp <- fp[[i]]
      for (j in seq_along(hr)) {
        cr <- if (j <= length(rr)) rr[[j]] else ""
        cp <- if (j <= length(rp)) rp[[j]] else ""
        va <- .es_num(cr); vb <- .es_num(cp)
        if (!is.na(va) && !is.na(vb)) {
          n_num <- n_num + 1L
          d <- abs(va - vb)
          if (d > max_abs) { max_abs <- d; peor <- sprintf("fila %d / %s", i + 1L, hr[[j]]) }
          if (!.dentro(va, vb, TOL_EST) && !.dentro(va, vb, TOL_P))
            n_fuera <- n_fuera + 1L
        } else if (!identical(cr, cp)) {
          n_dif_txt <- n_dif_txt + 1L
          if (!nzchar(peor)) peor <- sprintf("fila %d / %s (texto)", i + 1L, hr[[j]])
        }
      }
    }
  }
  ok <- header_igual && filas_ok && n_dif_txt == 0L && n_fuera == 0L && !is.null(hr)
  list(archivo = nombre, ok = ok, n_num = n_num, n_dif_texto = n_dif_txt,
       n_fuera_tol = n_fuera, max_dif_abs = max_abs, peor = peor,
       byte_identico = byte_id)
}

# =========================================================================
# Utilidades.
# =========================================================================
no_vacio <- function(ruta) file.exists(ruta) && file.info(ruta)$size > 0
.ver_paquetes <- function() {
  out <- character(0)
  for (p in R_PAQUETES) {
    v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "?")
    out <- c(out, sprintf("%s %s", p, v))
  }
  out
}

# =========================================================================
main <- function() {
  dur <- character(0); warn <- character(0); sinej <- character(0)
  chk_ok <- 0L; chk_tot <- 0L
  chequear <- function(desc, cond, blando = FALSE) {
    chk_tot <<- chk_tot + 1L
    if (isTRUE(cond)) chk_ok <<- chk_ok + 1L
    else if (blando) warn <<- c(warn, desc)
    else dur <<- c(dur, desc)
    isTRUE(cond)
  }
  # Tercer estado: el chequeo no se pudo ejecutar en esta corrida. No cuenta
  # como pase (no suma a chk_ok) ni como fallo (no suma a dur).
  sin_ejecutar <- function(desc, motivo) {
    chk_tot <<- chk_tot + 1L
    sinej <<- c(sinej, sprintf("%s [%s]", desc, motivo))
    FALSE
  }

  raiz <- RAIZ_REPO
  fuente <- fuente_datos(ARCHIVO_QPCR)

  # --- Se puede cruzar R contra Python en esta corrida? ------------------
  csv_r <- sort(list.files(RUTA_TABLAS_R, pattern = "\\.csv$"), method = "radix")
  csv_py <- sort(list.files(RUTA_TABLAS_PY, pattern = "\\.csv$"), method = "radix")
  solo_una <- nzchar(UNICA_IMPL)
  otro_lang <- if (LANG == "R") "python" else "R"
  # Con -Only, lo que haya en disco de la contraparte es de OTRA corrida: no
  # sirve de comparacion, aunque exista.
  cruce_posible <- !solo_una && length(csv_r) > 0L && length(csv_py) > 0L
  motivo_cruce <- if (solo_una)
    sprintf("corrida de una sola implementacion (-Only %s): 98_comparacion no corrio",
            UNICA_IMPL)
  else
    "no hay salidas de la contraparte en outputs/tables/ para comparar"

  # --- 1. Checklist de AGENTS 1 -------------------------------------
  for (f in c("AGENTS.md", "CLAUDE.md", "README.md", "ESTADO.md",
              "run_all.ps1", "requirements.txt", "renv.lock"))
    chequear(sprintf("existe/no vacio: %s", f), no_vacio(file.path(raiz, f)))

  for (s in SCRIPTS_NUM) {
    chequear(sprintf("script R/%s.R", s),
             no_vacio(file.path(raiz, "R", sprintf("%s.R", s))))
    chequear(sprintf("script python/%s.py", s),
             no_vacio(file.path(raiz, "python", sprintf("%s.py", s))))
  }
  for (t in SINTETICOS)
    chequear(sprintf("dato sintetico: %s", t),
             no_vacio(file.path(RUTA_DATOS_SINT, t)))
  for (fg in c(FIG_ACTO1, FIG_ACTO2))
    chequear(sprintf("figura: %s", fg), no_vacio(file.path(RUTA_FIGURAS, fg)))
  # comparacion_reporte.md y comparacion_R_python.csv los escribe 98_comparacion:
  # si no corrio, lo que haya en disco es de otra corrida -> NO_EJECUTADA.
  for (m in REPORTES_MD) {
    desc <- sprintf("reporte: outputs/tables/%s", m)
    if (m == "comparacion_reporte.md" && solo_una) sin_ejecutar(desc, motivo_cruce)
    else chequear(desc, no_vacio(file.path(RUTA_TABLAS, m)))
  }
  for (t in c("procedencia.csv", "verificaciones.csv", "comparacion_R_python.csv")) {
    desc <- sprintf("auditoria: outputs/tables/%s", t)
    if (t == "comparacion_R_python.csv" && solo_una) sin_ejecutar(desc, motivo_cruce)
    else chequear(desc, no_vacio(file.path(RUTA_TABLAS, t)))
  }

  if (solo_una) {
    chequear(sprintf("outputs/tables/%s/ tiene CSV", LANG),
             length(if (LANG == "R") csv_r else csv_py) > 0L)
    sin_ejecutar(sprintf("outputs/tables/%s/ tiene CSV", otro_lang), motivo_cruce)
  } else {
    chequear("outputs/tables/R/ tiene CSV", length(csv_r) > 0L)
    chequear("outputs/tables/python/ tiene CSV", length(csv_py) > 0L)
  }

  chequear("informe: docs/informe.html", no_vacio(file.path(RUTA_DOCS, "informe.html")))

  pdf_status <- "?"
  for (lg in c("python", "R")) {
    p <- file.path(RUTA_INTERMEDIOS, "render", lg, "pdf_status.txt")
    if (file.exists(p)) {
      pdf_status <- trimws(readChar(p, file.info(p)$size, useBytes = TRUE))
    }
  }
  pdf_ok <- no_vacio(file.path(RUTA_DOCS, "informe.pdf"))
  chequear("informe: docs/informe.pdf (blando si no hay motor de PDF)",
           pdf_ok || pdf_status != "ok", blando = !pdf_ok)

  # --- 2. Concordancia numerica R <-> Python ---------------------
  comps <- list(); n_byte <- 0L; n_fuera <- 0L; n_dif_txt <- 0L; peor_abs <- 0
  if (cruce_posible) {
    comunes <- sort(intersect(csv_r, csv_py), method = "radix")
    solo_r <- sort(setdiff(csv_r, csv_py), method = "radix")
    solo_py <- sort(setdiff(csv_py, csv_r), method = "radix")
    chequear("sin CSV huerfanos entre R/ y python/",
             !length(solo_r) && !length(solo_py))
    comps <- lapply(comunes, function(n)
      comparar_archivo(n, file.path(RUTA_TABLAS_R, n), file.path(RUTA_TABLAS_PY, n)))
    n_byte <- sum(vapply(comps, function(d) d$byte_identico, logical(1)))
    n_fuera <- sum(vapply(comps, function(d) d$n_fuera_tol, integer(1)))
    n_dif_txt <- sum(vapply(comps, function(d) d$n_dif_texto, integer(1)))
    peor_abs <- if (length(comps))
      max(vapply(comps, function(d) d$max_dif_abs, numeric(1))) else 0
    malos <- Filter(nzchar, vapply(comps, function(d)
      if (d$ok) "" else d$archivo, character(1)))
    chequear(sprintf("concordancia numerica R<->Python (%d CSV, tol 1e-6/1e-4)",
                     length(comps)),
             !length(malos) && n_fuera == 0L && n_dif_txt == 0L)
  } else {
    sin_ejecutar("sin CSV huerfanos entre R/ y python/", motivo_cruce)
    sin_ejecutar("concordancia numerica R<->Python", motivo_cruce)
  }

  # --- 3. Paridad de render R <-> Python -----------------------
  ren_r <- file.path(RUTA_INTERMEDIOS, "render", "R")
  ren_py <- file.path(RUTA_INTERMEDIOS, "render", "python")
  paridad_render <- "n/d"
  if (cruce_posible && dir.exists(ren_r) && dir.exists(ren_py)) {
    difs <- character(0)
    a <- leer_bytes(file.path(ren_r, "informe.textonly.html"))
    b <- leer_bytes(file.path(ren_py, "informe.textonly.html"))
    if (is.null(a) || is.null(b) || !identical(a, b))
      difs <- c(difs, "informe.textonly.html")
    for (md in sort(list.files(file.path(ren_r, "tables"), pattern = "\\.md$"),
                    method = "radix")) {
      x <- leer_bytes(file.path(ren_r, "tables", md))
      y <- leer_bytes(file.path(ren_py, "tables", md))
      if (is.null(y) || !identical(x, y)) difs <- c(difs, paste0("tables/", md))
    }
    for (md in sort(list.files(file.path(ren_py, "tables"), pattern = "\\.md$"),
                    method = "radix")) {
      live <- leer_bytes(file.path(RUTA_TABLAS, md))
      snap <- leer_bytes(file.path(ren_py, "tables", md))
      if (is.null(live) || !identical(live, snap))
        difs <- c(difs, paste0("disco != snapshot: ", md))
    }
    paridad_render <- if (!length(difs)) "OK" else
      paste0("DIFIEREN: ", paste(difs, collapse = ", "))
    chequear("paridad R<->Python de .md e informe (sin figuras)", !length(difs))
  } else {
    paridad_render <- "NO_EJECUTADA"
    sin_ejecutar("paridad R<->Python de .md e informe (sin figuras)", motivo_cruce)
  }

  # --- 4. verificaciones.csv: ninguna en FALSE (NO_EJECUTADA no es fallo) ----
  # NO_EJECUTADA marca una verificacion que no corrio en esta corrida (ej.
  # comparacion R<->Python en -Only R, o PDF sin motor disponible) -- no es
  # ni un pase ni un fallo, y no bloquea "TODAS LAS VERIFICACIONES PASARON".
  vf <- leer_csv(file.path(RUTA_TABLAS, "verificaciones.csv"))
  v_no_true <- character(0); v_no_ejec <- character(0)
  if (!is.null(vf$header) && "ok" %in% vf$header) {
    j <- match("ok", vf$header); jd <- match("id", vf$header)
    if (is.na(jd)) jd <- 1L
    ok_col <- vapply(vf$filas, function(f) if (j <= length(f)) f[[j]] else "", character(1))
    id_col <- vapply(vf$filas, function(f) if (jd <= length(f)) f[[jd]] else "", character(1))
    v_no_true <- id_col[ok_col == "FALSE"]
    v_no_ejec <- id_col[ok_col == "NO_EJECUTADA"]
  }
  chequear(sprintf("verificaciones.csv: %d filas, ninguna en FALSE (%d NO_EJECUTADA)",
                   length(vf$filas), length(v_no_ejec)),
           !is.null(vf$header) && !length(v_no_true))

  # --- 5. log de corrida -------------------------------------
  fecha <- format(Sys.Date(), "%Y-%m-%d")
  si <- Sys.info()
  seccion <- c(
    sprintf("## Corrida %s -- %s", LANG, fecha), "",
    sprintf("- fecha            : %s", fecha),
    sprintf("- plataforma       : %s %s (%s)", si[["sysname"]], si[["release"]],
            si[["machine"]]),
    sprintf("- lenguaje         : %s", R.version.string),
    sprintf("- paquetes         : %s", paste(.ver_paquetes(), collapse = "; ")),
    sprintf("- semilla          : %d", SEMILLA),
    sprintf("- fuente de datos  : %s", fuente),
    sprintf(paste0("- CSV R<->Python   : %d comparados; byte-identicos %d; ",
            "fuera de tol %d; dif texto %d; peor |dif| %.2e"),
            length(comps), n_byte, n_fuera, n_dif_txt, peor_abs),
    sprintf("- paridad render   : %s", paridad_render),
    sprintf("- informe.pdf      : %s (estado 12_informe: %s)",
            if (pdf_ok) "presente" else "ausente", pdf_status),
    sprintf("- verificaciones   : %d/%d chequeos duros OK; %d NO_EJECUTADA; avisos: %d",
            chk_ok, chk_tot, length(sinej), length(warn)),
    sprintf("- modo             : %s",
            if (solo_una) sprintf("una sola implementacion (%s)", UNICA_IMPL)
            else "R + Python"),
    sprintf("- filas NO_EJECUTADA de verificaciones.csv: %d%s",
            length(v_no_ejec),
            if (length(v_no_ejec)) paste0(" (", paste(v_no_ejec, collapse = ", "), ")")
            else ""),
    sprintf("- resultado        : %s", .resultado(dur, sinej, v_no_ejec))
  )
  if (length(dur)) seccion <- c(seccion, "", "### Fallos", paste0("  - ", dur))
  if (length(sinej))
    seccion <- c(seccion, "", "### No ejecutadas", paste0("  - ", sinej))
  if (length(warn)) seccion <- c(seccion, "", "### Avisos", paste0("  - ", warn))
  seccion <- c(seccion, "")
  .escribir_log(fecha, paste(seccion, collapse = "\n"))

  # --- salida legible --------------------------------------
  cat("== 99_verificar.R ==\n")
  cat(sprintf("  fuente = %s\n", fuente))
  cat(sprintf("  modo   = %s\n",
              if (solo_una) sprintf("una sola implementacion (%s)", UNICA_IMPL)
              else "R + Python"))
  cat(sprintf("  checklist: %d/%d chequeos duros OK; %d NO_EJECUTADA\n",
              chk_ok, chk_tot, length(sinej)))
  cat(sprintf(paste0("  CSV R<->Python: %d; byte-identicos %d; fuera de tol %d; ",
              "peor |dif| %.2e\n"), length(comps), n_byte, n_fuera, peor_abs))
  cat(sprintf("  paridad render: %s\n", paridad_render))
  cat(sprintf("  informe.pdf: %s (%s)\n", if (pdf_ok) "ok" else "ausente", pdf_status))
  cat(sprintf("  log -> logs/corrida_%s.txt\n", fecha))
  for (w in warn) cat(sprintf("  aviso: %s\n", w))
  for (s in sinej) cat(sprintf("  NO_EJECUTADA: %s\n", s))
  if (length(v_no_ejec))
    cat(sprintf("  NO_EJECUTADA (filas de verificaciones.csv): %s\n",
                paste(v_no_ejec, collapse = ", ")))
  if (length(dur)) {
    cat("\n")
    for (d in dur) cat(sprintf("  FALLA: %s\n", d))
    cat(sprintf("\n  *** VERIFICACION FALLIDA: %d chequeos duros no pasaron ***\n",
                length(dur)))
    quit(status = 1L)
  }
  cat(sprintf("\n  %s\n", .resultado(dur, sinej, v_no_ejec)))
  if (length(sinej) || length(v_no_ejec))
    cat(paste0("  (sin fallos, pero algo quedo sin ejecutar: para la validacion ",
               "completa correr .\\run_all.ps1 sin -Only)\n"))
}

# Una sola implementacion o una verificacion sin ejecutar => VERIFICACION
# PARCIAL. "TODAS LAS VERIFICACIONES PASARON" exige cero fallos y cero
# NO_EJECUTADA, en los chequeos duros y en verificaciones.csv (punto 4.3).
.resultado <- function(dur, sinej, v_no_ejec) {
  if (length(dur)) return(sprintf("VERIFICACION FALLIDA (%d)", length(dur)))
  if (length(sinej) || length(v_no_ejec))
    return(sprintf(paste0("VERIFICACION PARCIAL: 0 fallos, %d chequeo(s) y %d ",
                          "fila(s) NO_EJECUTADA"),
                   length(sinej), length(v_no_ejec)))
  "TODAS LAS VERIFICACIONES PASARON"
}

.escribir_log <- function(fecha, seccion_lang) {
  ruta <- file.path(RUTA_LOGS, sprintf("corrida_%s.txt", fecha))
  cab <- "# Corrida de verificacion -- reanalisis MIA-LPS\n"
  marca_ini <- sprintf("<!-- %s:inicio -->", LANG)
  marca_fin <- sprintf("<!-- %s:fin -->", LANG)
  bloque <- sprintf("%s\n%s\n%s", marca_ini, seccion_lang, marca_fin)
  txt <- if (file.exists(ruta)) {
    t <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(t) <- "UTF-8"; t
  } else paste0(cab, "\n")
  if (grepl(marca_ini, txt, fixed = TRUE) && grepl(marca_fin, txt, fixed = TRUE)) {
    pre <- strsplit(txt, marca_ini, fixed = TRUE)[[1]][1]
    post <- strsplit(txt, marca_fin, fixed = TRUE)[[1]][2]
    txt <- paste0(pre, bloque, post)
  } else {
    if (!endsWith(txt, "\n")) txt <- paste0(txt, "\n")
    txt <- paste0(txt, "\n", bloque, "\n")
  }
  con <- file(ruta, open = "wb"); writeBin(charToRaw(enc2utf8(txt)), con); close(con)
}

if (sys.nframe() == 0L) main()
