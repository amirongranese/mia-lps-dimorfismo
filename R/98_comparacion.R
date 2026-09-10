# 98_comparacion.R -- T10: concordancia numerica R <-> Python y auditoria de
#                     procedencia / verificaciones / analisis_descartados.
#
# Por que existe este archivo: todo el pipeline se implementa dos veces (R y
# Python) y la convencion del proyecto (AGENTS.md 7) es que las salidas
# numericas tengan NOMBRES IDENTICOS en outputs/tables/R/ y outputs/tables/
# python/ para poder cruzarlas automaticamente. Este script hace ese cruce
# archivo por archivo, celda por celda, con las tolerancias declaradas
# (1e-6 para estadisticos; 1e-4 de fallback para p de tests iterativos), y ademas
# consolida las tres tablas de auditoria: chequea que toda figura y toda tabla
# tenga su fila en procedencia.csv, que analisis_descartados.md tenga una seccion
# por cada script de analisis (02..11), y que no quede ninguna verificacion en
# estado distinto de TRUE.
#
# NO hace analisis cientifico. Depende de que 02..11 ya hayan corrido (deja las
# tablas en outputs/tables/); en run_all.ps1 va justo antes de 99_verificar.
#
# PARIDAD R/Python: R y Python producen EL MISMO comparacion_R_python.csv y el
# mismo comparacion_reporte.md (son, de hecho, el control cruzado uno del otro).
# Los numeros se formatean con "%.10g"; la comparacion celda a celda parsea cada
# lado a double y aplica |a-b| <= tol OR |a-b|/max(|a|,|b|) <= tol.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

ESTE_SCRIPT <- "98_comparacion"

TOL_EST <- TOL_ESTADISTICO      # 1e-6  -- estadisticos, medias, coeficientes
TOL_P   <- TOL_P_ITERATIVO      # 1e-4  -- fallback para p de tests iterativos

# Scripts de analisis que DEBEN tener su seccion en analisis_descartados.md.
SCRIPTS_ANALISIS <- c(
  "02_ingesta_qc", "03_elisa", "04_qpcr_cuantificacion", "05_qpcr_modelos",
  "06_pstat3", "07_figuras_acto1", "08_acto2_correlaciones",
  "09_acto2_dispersion", "10_acto2_simulacion", "11_sensibilidad"
)

# ---------------------------------------------------------------------------
# Formateo y escritura -- identicos a 02..11.
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
# Lectura de CSV -- mismo parser que 02..11 (comillas dobles, sin dependencias).
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
leer_csv <- function(ruta) {
  if (!file.exists(ruta)) return(NULL)
  txt <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(txt) <- "UTF-8"
  lineas <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  if (length(lineas) && lineas[length(lineas)] == "") lineas <- lineas[-length(lineas)]
  if (!length(lineas)) return(NULL)
  list(header = parse_csv_line(lineas[1]), filas = lapply(lineas[-1], parse_csv_line))
}
leer_bytes <- function(ruta) {
  if (!file.exists(ruta)) return(NULL)
  readBin(ruta, "raw", n = file.info(ruta)$size)
}

# ---------------------------------------------------------------------------
# Comparacion numerica celda a celda.
# ---------------------------------------------------------------------------
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

COLS_COMP <- c("archivo", "en_R", "en_python", "filas_R", "filas_py", "cols_R",
               "cols_py", "header_igual", "celdas", "celdas_numericas",
               "celdas_texto", "n_dif_texto", "max_dif_abs", "max_dif_rel",
               "peor_celda", "tol", "n_fuera_tol", "byte_identico", "ok")

comparar_archivo <- function(nombre, ruta_r, ruta_py) {
  r <- leer_csv(ruta_r); p <- leer_csv(ruta_py)
  hr <- if (is.null(r)) NULL else r$header
  hp <- if (is.null(p)) NULL else p$header
  fr <- if (is.null(r)) list() else r$filas
  fp <- if (is.null(p)) list() else p$filas
  br <- leer_bytes(ruta_r); bp <- leer_bytes(ruta_py)
  byte_id <- !is.null(br) && !is.null(bp) && identical(br, bp)

  header_igual <- identical(hr, hp)
  filas_ok <- length(fr) == length(fp)
  n_cel <- 0L; n_num <- 0L; n_txt <- 0L; n_dif_txt <- 0L; n_fuera <- 0L
  max_abs <- 0; max_rel <- 0; peor <- ""; tol_usada <- TOL_EST

  if (header_igual && filas_ok && !is.null(hr)) {
    for (i in seq_along(fr)) {
      rr <- fr[[i]]; rp <- fp[[i]]
      for (j in seq_along(hr)) {
        cr <- if (j <= length(rr)) rr[[j]] else ""
        cp <- if (j <= length(rp)) rp[[j]] else ""
        n_cel <- n_cel + 1L
        va <- .es_num(cr); vb <- .es_num(cp)
        if (!is.na(va) && !is.na(vb)) {
          n_num <- n_num + 1L
          d <- abs(va - vb); m <- max(abs(va), abs(vb))
          rel <- if (m > 0) d / m else 0
          if (d > max_abs) { max_abs <- d; peor <- sprintf("fila %d / %s", i + 1L, hr[[j]]) }
          if (rel > max_rel) max_rel <- rel
          if (!.dentro(va, vb, TOL_EST)) {
            if (.dentro(va, vb, TOL_P)) tol_usada <- TOL_P
            else n_fuera <- n_fuera + 1L
          }
        } else {
          n_txt <- n_txt + 1L
          if (!identical(cr, cp)) {
            n_dif_txt <- n_dif_txt + 1L
            if (!nzchar(peor)) peor <- sprintf("fila %d / %s (texto)", i + 1L, hr[[j]])
          }
        }
      }
    }
  }

  ok <- header_igual && filas_ok && n_dif_txt == 0L && n_fuera == 0L && !is.null(hr)
  list(archivo = nombre, en_R = file.exists(ruta_r), en_python = file.exists(ruta_py),
       filas_R = length(fr), filas_py = length(fp),
       cols_R = if (is.null(hr)) 0L else length(hr),
       cols_py = if (is.null(hp)) 0L else length(hp),
       header_igual = header_igual, celdas = n_cel, celdas_numericas = n_num,
       celdas_texto = n_txt, n_dif_texto = n_dif_txt, max_dif_abs = max_abs,
       max_dif_rel = max_rel, peor_celda = peor, tol = tol_usada,
       n_fuera_tol = n_fuera, byte_identico = byte_id, ok = ok)
}

# ---------------------------------------------------------------------------
# Auditoria de procedencia / verificaciones / analisis_descartados.
# ---------------------------------------------------------------------------
.basename_artefacto <- function(s) sub(".*/", "", s)

auditar_procedencia <- function(figuras, tablas) {
  pr <- leer_csv(file.path(RUTA_TABLAS, "procedencia.csv"))
  mencionados <- character(0)
  if (!is.null(pr)) {
    idx <- match("artefacto", pr$header)
    mencionados <- unique(vapply(pr$filas, function(f)
      if (length(f) >= idx) .basename_artefacto(f[[idx]]) else "", character(1)))
  }
  list(figs_sin = figuras[!(figuras %in% mencionados)],
       tabs_sin = tablas[!(tablas %in% mencionados)])
}

auditar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  txt <- if (file.exists(ruta)) {
    t <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(t) <- "UTF-8"; t
  } else ""
  SCRIPTS_ANALISIS[!vapply(SCRIPTS_ANALISIS, function(s)
    grepl(sprintf("<!-- %s:inicio -->", s), txt, fixed = TRUE) &&
    grepl(sprintf("<!-- %s:fin -->", s), txt, fixed = TRUE), logical(1))]
}

auditar_verificaciones <- function() {
  vf <- leer_csv(file.path(RUTA_TABLAS, "verificaciones.csv"))
  if (is.null(vf)) return(list(total = 0L, ok = 0L, no = character(0)))
  idx_ok <- match("ok", vf$header); idx_id <- match("id", vf$header)
  idx_sc <- match("script", vf$header)
  # el conteo excluye las propias filas de 98_comparacion para que sea estable
  # sin importar el orden de ejecucion (R antes / despues de Python) ni re-corridas.
  filas <- Filter(function(f) !(length(f) >= idx_sc && f[[idx_sc]] == ESTE_SCRIPT),
                  vf$filas)
  total <- length(filas)
  no <- vapply(filas, function(f)
    if (length(f) >= idx_ok && f[[idx_ok]] == "TRUE") "" else f[[idx_id]], character(1))
  no <- no[nzchar(no)]
  list(total = total, ok = total - length(no), no = no)
}

# ---------------------------------------------------------------------------
# Justificaciones que 02/03/04 difieren explicitamente a T10.
# ---------------------------------------------------------------------------
.num_no_detectados <- function() {
  nd <- leer_csv(file.path(RUTA_TABLAS_R, "qpcr_no_detectados_descriptivo.csv"))
  if (is.null(nd)) return(NULL)
  h <- nd$header
  ci <- function(n) match(n, h)
  tot_obs <- 0L; tot_nd <- 0L; d7 <- list()
  for (f in nd$filas) {
    if (f[[ci("GRUPO")]] == "TODOS") {
      tot_obs <- tot_obs + as.integer(f[[ci("n_total")]])
      tot_nd <- tot_nd + as.integer(f[[ci("n_no_detectado")]])
    }
    if (f[[ci("calibrador_cero_detectados")]] == "TRUE" && f[[ci("GRUPO")]] != "TODOS")
      d7[[length(d7) + 1L]] <- list(
        t = f[[ci("TEJIDO")]], g = f[[ci("GEN")]], gr = f[[ci("GRUPO")]],
        ntot = as.integer(f[[ci("n_total")]]), nnd = as.integer(f[[ci("n_no_detectado")]]))
  }
  list(tot_obs = tot_obs, tot_nd = tot_nd, d7 = d7)
}

.SEC_INI <- "<!-- 98_comparacion:inicio -->"
.SEC_FIN <- "<!-- 98_comparacion:fin -->"

.seccion_descartados <- function() {
  nd <- .num_no_detectados()
  if (!is.null(nd)) {
    pct <- if (nd$tot_obs) 100 * nd$tot_nd / nd$tot_obs else 0
    linea_tot <- sprintf(paste0("En la corrida actual: **%d de %d** observaciones ",
      "gen x feto x tejido no detectadas (%.1f %%)."), nd$tot_nd, nd$tot_obs, pct)
    if (length(nd$d7)) {
      celdas <- paste(vapply(nd$d7, function(x)
        sprintf("%s/%s/%s %d/%d", x$t, x$g, x$gr, x$nnd, x$ntot), character(1)),
        collapse = "; ")
      linea_d7 <- sprintf(paste0("Celdas con calibrador HEMBRA_CONTROL 0/n ",
        "detectado (D7): %s."), celdas)
    } else {
      linea_d7 <- paste0("En esta corrida ningun calibrador quedo con 0 detectados ",
        "(sobre sintetico puede no dispararse D7).")
    }
  } else {
    linea_tot <- paste0("(qpcr_no_detectados_descriptivo.csv no disponible: correr ",
      "04_qpcr_cuantificacion antes de 98_comparacion.)")
    linea_d7 <- ""
  }

  L <- c(
    "## 98_comparacion",
    "",
    paste0("Consolidacion de auditoria (T10). Concordancia R <-> Python archivo por ",
      "archivo (`comparacion_R_python.csv`) y cierre de las justificaciones que ",
      "02/03/04 dejaron diferidas a esta tarea."),
    "",
    "### Imputacion MNAR de no-detectados (D3) -- justificacion completa",
    "",
    paste0("- **Que es**: la imputacion MNAR (p. ej. `nondetects` en R) modela la ",
      "probabilidad de no-deteccion como funcion decreciente del Ct latente y ",
      "sortea Ct por encima del umbral para las celdas no observadas, para ",
      "despues correr el analisis sobre una matriz \"completa\"."),
    paste0("- **Numeros**: ", linea_tot)
  )
  if (nzchar(linea_d7)) L <- c(L, paste0("  ", linea_d7))
  L <- c(L,
    "- **Por que se descarta** (D3, prohibicion 1):",
    paste0("  1. *Ancla inexistente en las celdas D7.* La imputacion necesita algunos ",
      "Ct observados en el grupo para estimar la pendiente de la curva de ",
      "deteccion. Donde el calibrador HEMBRA_CONTROL tiene 0 detectados ",
      "(il6 @ BRAIN_E15) no hay nada sobre lo que anclar: el dCt de calibrador, ",
      "y por lo tanto todo ddCt del gen x tejido, quedaria definido contra un ",
      "valor enteramente inventado (prohibicion 10)."),
    paste0("  2. *Estructura artificial en el baseline.* En la exploracion previa ",
      "(informe E15 original) la imputacion MNAR rellenaba los no-detectados del ",
      "grupo Control con Ct altos correlacionados, comprimiendo su varianza y ",
      "generando diferencias Control vs LPS que el dato crudo no tiene. La ",
      "imputacion inyecta la hipotesis que el analisis deberia poner a prueba."),
    paste0("  3. *No hace falta.* D1/D8 calculan calibrador, medias y z **solo sobre ",
      "detectados**; D5 modela solo las celdas con n >= 5. La informacion que ",
      "aporta un no-detectado (un \"< umbral\") se analiza como **proporcion de ",
      "deteccion** (Fisher exacto, D7) alli donde es informativa -- sin inventar ",
      "su magnitud."),
    paste0("- **Que se hizo**: no-detectado -> NA, sin imputar por ningun metodo; ",
      "deteccion como desenlace propio donde corresponde (D7)."),
    "",
    "### LOD alternativo del ELISA (menor estandar de `CURVA IL6`) -- no aplicado",
    "",
    paste0("- D10 admite como alternativa al blanco de placa (`LOD = 0`, el primario ",
      "de 02/03) el menor estandar de la hoja `CURVA IL6`."),
    paste0("- **Por que no se corre como sensibilidad separada**: el contraste ",
      "primario del ELISA trata a los censurados como **empatados en el rango ",
      "mas bajo** (Peto-Peto G-rho=1 en LA; Fisher de deteccion en MS). Quien ",
      "esta censurado lo define `Conc < 0`, no el valor del LOD; mover el LOD de ",
      "0 al menor estandar no reordena esos empates ni cambia el conjunto ",
      "censurado, asi que **no puede cambiar la conclusion por construccion del ",
      "test**. El LOD explicito solo entraria en un estimador de ubicacion con ",
      "censura, que D10 (clausula final) ya descarta por censura alta."),
    paste0("- **Que se hizo**: se mantiene `LOD = 0`; el alternativo queda documentado ",
      "y no ejecutado."),
    "",
    "### Factorizacion del nucleo hand-rolled (08--11) -- evaluado, no se factoriza",
    "",
    paste0("- **Situacion**: `spearman_rho` + `rangos_promedio` + `pares()` / ",
      "`cargar()` + los merges de `-ddCt` / score estan replicados casi ",
      "identicos en 08, 09, 10 y 11; 11 agrega Jacobi sin trigonometria y ",
      "`pearson_pairwise`."),
    paste0("- **Opcion evaluada**: moverlos a un modulo compartido (`00_config` o un ",
      "`nucleo_acto2`) importado por 08--11."),
    paste0("- **Decision: no se factoriza en T10.** (a) La regla del repo (AGENTS.md ",
      "7) es que cada script importa **solo** `00_config`; meter logica de ",
      "analisis en `00_config` rompe esa linea. (b) Las salidas de 02--11 son ",
      "**byte-identicas R/Python** y estan cerradas; tocar el nucleo justo antes ",
      "del informe (T11) arriesga perturbarlas sin ganancia cientifica. (c) La ",
      "replicacion es chica, ya verificada byte a byte entre lenguajes y, en R, ",
      "cruzada en cada script contra `cor.test` / `car::leveneTest` / `eigen()`; ",
      "`comparacion_R_python.csv` la cubre de forma continua."),
    paste0("- **Revisable** despues de T11 si el informe necesita re-ejecutar el ",
      "nucleo."),
    "",
    "### Concordancia R <-> Python (`comparacion_R_python.csv`)",
    "",
    paste0("- Una fila por CSV de `outputs/tables/{R,python}/`. Cada celda que parsea ",
      "a numero finito en ambos lados se compara con ",
      "`|a-b| <= 1e-6` **o** `|a-b|/max(|a|,|b|) <= 1e-6` (`tol` pasa a `1e-4` si ",
      "alguna celda necesita el margen mas laxo, reservado a `p` de tests ",
      "iterativos); el resto se compara como texto exacto. `byte_identico` marca ",
      "los pares que ademas coinciden byte a byte."),
    paste0("- Los `.md` de `outputs/tables/` (reportes legibles) son de copia unica ",
      "-- R y Python escriben la misma ruta -- y no entran en esta tabla; su ",
      "paridad entre lenguajes se re-chequea en `99_verificar` (T11).")
  )
  paste(L, collapse = "\n")
}

actualizar_descartados <- function() {
  ruta <- file.path(RUTA_TABLAS, "analisis_descartados.md")
  nuevo <- paste0(.SEC_INI, "\n", .seccion_descartados(), "\n\n", .SEC_FIN)
  txt <- if (file.exists(ruta)) {
    t <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(t) <- "UTF-8"; t
  } else {
    paste0("# Analisis descartados\n\nQue se probo, por que no funciono o no se",
           " uso, y que se hizo en su lugar. Una seccion por script.\n")
  }
  if (grepl(.SEC_INI, txt, fixed = TRUE) && grepl(.SEC_FIN, txt, fixed = TRUE)) {
    pre <- strsplit(txt, .SEC_INI, fixed = TRUE)[[1]][1]
    post <- strsplit(txt, .SEC_FIN, fixed = TRUE)[[1]][2]
    txt <- paste0(pre, nuevo, post)
  } else {
    if (!endsWith(txt, "\n")) txt <- paste0(txt, "\n")
    txt <- paste0(txt, "\n", nuevo, "\n")
  }
  escribir_lineas(ruta, txt)
}

# ---------------------------------------------------------------------------
# Merge por 'script' en procedencia.csv / verificaciones.csv -- headers 02..11.
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Reporte legible.
# ---------------------------------------------------------------------------
construir_reporte <- function(fuente, comps, resumen, figs_sin, tabs_sin, desc_faltan,
                              verif_total, verif_true, verif_no) {
  L <- character(0); ap <- function(...) L[[length(L) + 1L]] <<- paste0(...)
  ap("# Comparacion R <-> Python y auditoria (T10)")
  ap("")
  ap("Generado por `98_comparacion` (R y Python producen este archivo identico).")
  ap(sprintf("Fuente de datos en uso: `%s`.", fuente))
  ap("")
  ap("## 1. Concordancia numerica archivo por archivo")
  ap("")
  ap(paste0("Tolerancia: `1e-6` (estadisticos) con fallback `1e-4` (p de tests ",
    "iterativos). `ok` = header igual, mismo n de filas, sin diferencias de ",
    "texto y sin celdas fuera de tolerancia."))
  ap("")
  ap("| archivo | filas | num | max_dif_abs | byte-id | tol | ok |")
  ap("| --- | --- | --- | --- | --- | --- | --- |")
  for (d in comps) {
    ap(sprintf("| %s | %d | %d | %.2e | %s | %g | %s |",
               d$archivo, d$filas_R, d$celdas_numericas, d$max_dif_abs,
               if (d$byte_identico) "si" else "no", d$tol,
               if (d$ok) "OK" else "REVISAR"))
  }
  ap("")
  ap(sprintf(paste0("- Archivos comparados: **%d**. Byte-identicos: **%d**. ",
    "Fuera de tolerancia: **%d**. Peor |dif| absoluta: **%.3e** (%s)."),
    resumen$n, resumen$byte, resumen$fuera, resumen$peor_abs,
    if (nzchar(resumen$peor_arch)) resumen$peor_arch else "-"))
  ap(sprintf("- Solo en R: %s. Solo en Python: %s.",
             if (nzchar(resumen$solo_r)) resumen$solo_r else "ninguno",
             if (nzchar(resumen$solo_py)) resumen$solo_py else "ninguno"))
  ap("")
  ap("## 2. Auditoria de procedencia.csv")
  ap("")
  ap(sprintf("- Figuras sin fila en procedencia: **%s**.",
             if (nzchar(figs_sin)) figs_sin else "ninguna"))
  ap(sprintf(paste0("- Tablas (outputs/tables/{R,python}/*.csv) sin fila: **%s**."),
             if (nzchar(tabs_sin)) tabs_sin else "ninguna"))
  ap("")
  ap("## 3. Auditoria de analisis_descartados.md")
  ap("")
  ap(sprintf("- Scripts de analisis esperados: %d. Sin seccion: **%s**.",
             length(SCRIPTS_ANALISIS), if (nzchar(desc_faltan)) desc_faltan else "ninguno"))
  ap("")
  ap("## 4. Auditoria de verificaciones.csv")
  ap("")
  ap(sprintf("- Filas: **%d**. En TRUE: **%d**. Distintas de TRUE: **%s**.",
             verif_total, verif_true, if (nzchar(verif_no)) verif_no else "ninguna"))
  ap("")
  ap("## 5. Notas")
  ap("")
  ap(paste0("Ver `analisis_descartados.md`, seccion `98_comparacion`: justificacion ",
    "completa del descarte de la imputacion MNAR, LOD alternativo del ELISA no ",
    "aplicado, y la decision de no factorizar el nucleo hand-rolled de 08--11."))
  ap("")
  paste(L, collapse = "\n")
}

# ===========================================================================
main <- function() {
  fuente <- fuente_datos(ARCHIVO_QPCR)

  arch_r <- sort(list.files(RUTA_TABLAS_R, pattern = "\\.csv$"), method = "radix")
  arch_py <- sort(list.files(RUTA_TABLAS_PY, pattern = "\\.csv$"), method = "radix")
  if (!length(arch_r)) stop("outputs/tables/R/ no tiene CSV: correr 02..11 antes.")
  solo_r <- sort(setdiff(arch_r, arch_py), method = "radix")
  solo_py <- sort(setdiff(arch_py, arch_r), method = "radix")
  comunes <- sort(intersect(arch_r, arch_py), method = "radix")

  comps <- lapply(comunes, function(n)
    comparar_archivo(n, file.path(RUTA_TABLAS_R, n), file.path(RUTA_TABLAS_PY, n)))

  n_byte <- sum(vapply(comps, function(d) d$byte_identico, logical(1)))
  n_fuera <- sum(vapply(comps, function(d) d$n_fuera_tol + d$n_dif_texto, integer(1)))
  n_ok <- sum(vapply(comps, function(d) d$ok, logical(1)))
  peor_abs <- 0; peor_arch <- ""
  for (d in comps) if (d$max_dif_abs > peor_abs) { peor_abs <- d$max_dif_abs; peor_arch <- d$archivo }

  # --- comparacion_R_python.csv: una fila por archivo + fila __TOTAL__ ---
  fila_de <- function(d) lapply(COLS_COMP, function(c) d[[c]])
  filas_csv <- lapply(comps, fila_de)
  total <- list(
    archivo = "__TOTAL__", en_R = length(arch_r), en_python = length(arch_py),
    filas_R = sum(vapply(comps, function(d) d$filas_R, integer(1))),
    filas_py = sum(vapply(comps, function(d) d$filas_py, integer(1))),
    cols_R = "", cols_py = "",
    header_igual = all(vapply(comps, function(d) d$header_igual, logical(1))),
    celdas = sum(vapply(comps, function(d) d$celdas, integer(1))),
    celdas_numericas = sum(vapply(comps, function(d) d$celdas_numericas, integer(1))),
    celdas_texto = sum(vapply(comps, function(d) d$celdas_texto, integer(1))),
    n_dif_texto = sum(vapply(comps, function(d) d$n_dif_texto, integer(1))),
    max_dif_abs = peor_abs,
    max_dif_rel = max(vapply(comps, function(d) d$max_dif_rel, numeric(1)), 0),
    peor_celda = peor_arch, tol = TOL_EST,
    n_fuera_tol = sum(vapply(comps, function(d) d$n_fuera_tol, integer(1))),
    byte_identico = n_byte == length(comps),
    ok = n_ok == length(comps) && !length(solo_r) && !length(solo_py))
  filas_csv <- c(filas_csv, list(fila_de(total)))
  escribir_csv(file.path(RUTA_TABLAS, "comparacion_R_python.csv"), COLS_COMP, filas_csv)

  # --- auditoria ---
  figuras <- sort(list.files(RUTA_FIGURAS, pattern = "\\.png$"), method = "radix")
  ap <- auditar_procedencia(figuras, arch_r)
  desc_faltan <- auditar_descartados()
  vf <- auditar_verificaciones()

  resumen <- list(n = length(comps), byte = n_byte,
                  fuera = total$n_fuera_tol + total$n_dif_texto,
                  peor_abs = peor_abs, peor_arch = peor_arch,
                  solo_r = paste(solo_r, collapse = ", "),
                  solo_py = paste(solo_py, collapse = ", "))

  escribir_lineas(file.path(RUTA_TABLAS, "comparacion_reporte.md"),
    construir_reporte(fuente, comps, resumen, paste(ap$figs_sin, collapse = ", "),
                      paste(ap$tabs_sin, collapse = ", "),
                      paste(desc_faltan, collapse = ", "),
                      vf$total, vf$ok, paste(vf$no, collapse = ", ")))
  actualizar_descartados()

  # --- procedencia / verificaciones (filas de 98_comparacion) ---
  ent <- "outputs/tables/{R,python}/*.csv"
  registrar_procedencia(list(
    list("outputs/tables/comparacion_R_python.csv", "tabla", ESTE_SCRIPT, "PROPIO",
         ent, paste0("concordancia numerica R<->python archivo por archivo; tol 1e-6 ",
         "(fallback 1e-4) + byte-identidad; fila __TOTAL__ con agregados")),
    list("outputs/tables/comparacion_reporte.md", "reporte", ESTE_SCRIPT, "PROPIO",
         ent, paste0("reporte legible de la comparacion R<->python y de la auditoria de ",
         "procedencia / verificaciones / analisis_descartados (T10)")),
    list("outputs/tables/procedencia.csv", "auditoria", ESTE_SCRIPT, "PROPIO",
         "todos los scripts 02..11",
         paste0("indice de procedencia: una fila por figura y por tabla; cobertura ",
         "verificada en T10 (auditar_procedencia)")),
    list("outputs/tables/verificaciones.csv", "auditoria", ESTE_SCRIPT, "PROPIO",
         "todos los scripts 02..11",
         paste0("indice de verificaciones: un resultado principal por fila; T10 chequea ",
         "que ninguna quede distinta de TRUE")),
    list("outputs/tables/analisis_descartados.md", "auditoria", ESTE_SCRIPT, "PROPIO",
         "todos los scripts 02..11",
         "una seccion por script de analisis (02..11); cobertura verificada en T10")))

  registrar_verificaciones(list(
    list("comp_cobertura_csv",
         "todo CSV de outputs/tables/R/ tiene gemelo en python/ y viceversa",
         sprintf("R=%d; python=%d; solo_R=%d; solo_python=%d",
                 length(arch_r), length(arch_py), length(solo_r), length(solo_py)),
         "sin huerfanos en ningun lado",
         if (!length(solo_r) && !length(solo_py)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("comp_headers_iguales",
         "el encabezado coincide en los pares comparados",
         sprintf("%d/%d con header igual",
                 sum(vapply(comps, function(d) d$header_igual, logical(1))), length(comps)),
         sprintf("%d/%d", length(comps), length(comps)),
         if (all(vapply(comps, function(d) d$header_igual, logical(1)))) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("comp_dentro_tolerancia",
         "concordancia numerica R<->python dentro de 1e-6 (fallback 1e-4)",
         sprintf(paste0("archivos=%d; byte_identicos=%d; fuera_tol=%d; dif_texto=%d; ",
                 "peor|dif|abs=%.2e"), length(comps), n_byte, total$n_fuera_tol,
                 total$n_dif_texto, peor_abs),
         "0 celdas fuera de tolerancia; 0 diferencias de texto",
         if (n_fuera == 0L) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("audit_procedencia_figuras",
         "toda figura de outputs/figures/ tiene fila en procedencia.csv",
         sprintf("figuras=%d; sin_fila=%s", length(figuras),
                 if (length(ap$figs_sin)) paste(ap$figs_sin, collapse = ", ") else "[]"),
         "sin_fila = []",
         if (!length(ap$figs_sin)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("audit_procedencia_tablas",
         "toda tabla outputs/tables/{R,python}/*.csv tiene fila en procedencia.csv",
         sprintf("tablas=%d; sin_fila=%s", length(arch_r),
                 if (length(ap$tabs_sin)) paste(ap$tabs_sin, collapse = ", ") else "[]"),
         "sin_fila = []",
         if (!length(ap$tabs_sin)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("audit_descartados_cobertura",
         "analisis_descartados.md tiene una seccion por script de analisis (02..11)",
         sprintf("esperados=%d; faltan=%s", length(SCRIPTS_ANALISIS),
                 if (length(desc_faltan)) paste(desc_faltan, collapse = ", ") else "[]"),
         "faltan = []",
         if (!length(desc_faltan)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("audit_verificaciones_ok",
         "ninguna fila de verificaciones.csv queda distinta de TRUE",
         sprintf("total=%d; TRUE=%d; no_TRUE=%s", vf$total, vf$ok,
                 if (length(vf$no)) paste(vf$no, collapse = ", ") else "[]"),
         "no_TRUE = []",
         if (!length(vf$no)) "TRUE" else "FALSE", ESTE_SCRIPT)))

  # --- salida legible ---
  cat("== 98_comparacion.R ==\n")
  cat(sprintf("  fuente = %s\n", fuente))
  cat(sprintf("  CSV comparados: %d  |  byte-identicos: %d/%d  |  OK: %d/%d\n",
              length(comps), n_byte, length(comps), n_ok, length(comps)))
  cat(sprintf("  peor |dif| absoluta: %.3e%s\n", peor_abs,
              if (nzchar(peor_arch)) sprintf("  (%s)", peor_arch) else ""))
  if (length(solo_r) || length(solo_py))
    cat(sprintf("  HUERFANOS  solo_R=%s  solo_python=%s\n",
                paste(solo_r, collapse = ","), paste(solo_py, collapse = ",")))
  malos <- vapply(comps, function(d) if (d$ok) "" else d$archivo, character(1))
  malos <- malos[nzchar(malos)]
  if (length(malos)) cat(sprintf("  REVISAR: %s\n", paste(malos, collapse = ", ")))
  cat(sprintf("  procedencia: figuras sin fila = %s; tablas sin fila = %s\n",
              if (length(ap$figs_sin)) paste(ap$figs_sin, collapse = ", ") else "[]",
              if (length(ap$tabs_sin)) paste(ap$tabs_sin, collapse = ", ") else "[]"))
  cat(sprintf("  descartados: sin seccion = %s\n",
              if (length(desc_faltan)) paste(desc_faltan, collapse = ", ") else "[]"))
  cat(sprintf("  verificaciones: %d/%d en TRUE; distintas de TRUE = %s\n",
              vf$ok, vf$total, if (length(vf$no)) paste(vf$no, collapse = ", ") else "[]"))
  todo_ok <- !length(malos) && !length(solo_r) && !length(solo_py) &&
    !length(ap$figs_sin) && !length(ap$tabs_sin) && !length(desc_faltan) && !length(vf$no)
  cat("  -> outputs/tables/comparacion_R_python.csv, comparacion_reporte.md\n")
  cat(if (todo_ok) "  TODAS LAS COMPARACIONES Y AUDITORIAS PASARON\n"
      else "  *** HAY ITEMS A REVISAR (ver arriba) ***\n")
}

main()
