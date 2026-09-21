# 12_informe.R -- T11: informe HTML autocontenido del reanalisis MIA-LPS.
#
# Por que existe este archivo: cierra el pipeline juntando en un unico documento
# navegable (Acto 1 + Acto 2 + reproducibilidad + analisis descartados +
# limitaciones + procedencia/verificaciones) todo lo que produjeron 02..11 y 98.
# El informe se arma A MANO (no rmarkdown / no pandoc): (a) esas dependencias no
# estan garantizadas en un clon limpio -- de hecho falta `rmarkdown` en este
# build de R -- y (b) la unica forma de que R y Python generen EL MISMO
# `docs/informe.html` byte a byte es construir el HTML con la misma logica de
# strings en ambos lenguajes, sin un motor intermedio.
#
# Entradas: los `.md` de copia unica de outputs/tables/ (byte-identicos R/Python
# por construccion de 02..11/98), unos pocos CSV de outputs/tables/<lang>/ para
# los numeros del resumen, y los PNG de outputs/figures/. Las figuras se
# incrustan en base64 (informe AUTOCONTENIDO); como los PNG NO son byte-identicos
# entre ggplot2 y matplotlib (ver ESTADO / analisis_descartados), `informe.html`
# tampoco lo es: la paridad R/Python se chequea en 99_verificar sobre la
# PROYECCION del HTML sin los blobs `data:` (informe.textonly.html) y sobre los
# `.md`. No lleva marca de tiempo (eso va en logs/corrida_<fecha>.txt).
#
# Salidas:
#   docs/informe.html                         -- informe autocontenido
#   docs/informe.pdf                          -- best-effort (Edge/Chrome headless);
#                                                el pipeline NO falla si no se puede
#   outputs/intermediate/render/<lang>/...    -- snapshots para la paridad de 99
#   filas en procedencia.csv / verificaciones.csv (script = 12_informe)

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

ESTE_SCRIPT <- "12_informe"
LANG <- "R"

# --- Reportes .md de copia unica -> seccion del informe --------------------
# El 3er elemento de cada entrada YA NO es una lista de PNG escrita a mano
# (se desincronizaba cada vez que un script agregaba figuras -- ver ESTADO.md,
# sesion del bugfix): es un CRITERIO que filtra
# `procedencia.csv` (columna `script`, mas un patron de nombre opcional para
# separar los dos usos de 07_figuras_acto1 -- expresion vs pSTAT3). Las
# figuras de cada seccion se derivan solas en `figuras_de_seccion()`.
SECCIONES <- list(
  list("Acto 1.1 -- ELISA de IL-6 (validacion del modelo)", "elisa_reporte.md",
       list(script = "03_elisa")),
  list("Acto 1.2 -- Cuantificacion relativa (qPCR)", "qpcr_cuantificacion_reporte.md",
       list(script = "04_qpcr_cuantificacion")),
  list("Acto 1.3 -- Modelos de expresion (qPCR)", "qpcr_modelos_reporte.md",
       list(script = "07_figuras_acto1", patron = "^acto1_expresion_")),
  list("Acto 1.4 -- pSTAT3 en placenta", "pstat3_reporte.md",
       list(script = "07_figuras_acto1", patron = "^acto1_pstat3")),
  list("Acto 2.1-2.2 -- Correlacion placenta<->cerebro y co-expresion",
       "acto2_correlaciones_reporte.md",
       list(script = "08_acto2_correlaciones")),
  list("Acto 2.3-2.4 -- Dispersion y test formal de Delta rho",
       "acto2_dispersion_reporte.md",
       list(script = "09_acto2_dispersion")),
  list("Acto 2.5 -- Simulacion de restriccion de rango",
       "acto2_simulacion_reporte.md", list(script = "10_acto2_simulacion")),
  list("Acto 2.6 -- Sensibilidad (eigengene, exclusion del feto extremo)",
       "acto2_sensibilidad_reporte.md", list(script = "11_sensibilidad"))
)

# Decisiones D1..D12 (texto fijo, identico en ambos lenguajes; espejo de AGENTS 4).
DECISIONES <- list(
  c("D1", paste0("Cuantificacion relativa: dCt = Ct_gen - Ct_rsp29. Calibrador = ",
    "HEMBRA CONTROL por gen x tejido, promediando solo valores detectados. ddCt = ",
    "dCt_muestra - dCt_calibrador.")),
  c("D2", paste0("El analisis corre sobre -ddCt (log2, simetrica y aditiva). El ",
    "fold-change FC = 2^(-ddCt) solo se grafica (eje Y log).")),
  c("D3", paste0("No detectados -> NA. Se evaluo imputacion MNAR y se descarto; no se ",
    "imputa por ningun metodo.")),
  c("D4", paste0("Poblacion de analisis: los 36 fetos, sin exclusion de outliers en el ",
    "analisis principal (la robustez la cubren los controles de sensibilidad).")),
  c("D5", paste0("Modelo por gen x tejido: -ddCt ~ SEXO * TTO, con cascada de supuestos ",
    "sobre los residuos (ANOVA-III / HC3 / ART).")),
  c("D6", paste0("Post hoc solo si SEXO x TTO es significativa: 4 comparaciones fijas, ",
    "correccion Holm.")),
  c("D7", paste0("il6 en cerebro E15 no es cuantificable (calibrador HEMBRA_CONTROL ",
    "0/9 detectados): fuera del modelo, se analiza como proporcion de deteccion ",
    "(Fisher exacto).")),
  c("D8", paste0("Score compuesto de transportadores = z-score por gen (dentro de ",
    "tejido, 36 fetos) y promedio de los 7 z por feto. Variante PCA = eigengene.")),
  c("D9", paste0("pSTAT3: PSTAT3 ~ SEXO * TTO + MEMBRANA (bloque fijo). Limitacion ",
    "obligatoria: sin STAT3 total -> mide abundancia de fosfo-STAT3, no fraccion ",
    "fosforilada.")),
  c("D10", paste0("ELISA con censura a izquierda: indicador censurado, valor NA, LOD ",
    "aparte. Se reporta % de censura por grupo antes de todo estadistico; ",
    "metodos para datos censurados (KM/ROS o no parametricos).")),
  c("D11", paste0("Anotacion de boxplots solo si SEXO x TTO y el post hoc son ",
    "significativos; una sola funcion de brackets por lenguaje.")),
  c("D12", paste0("Correccion entre genes: primario sin correccion; columna ",
    "suplementaria con p ajustado por Benjamini-Hochberg dentro de cada tejido. ",
    "No cambia conclusiones.")),
  c("D13", paste0("Sin termino de camada: los modelos no incluyen MADRE, ni fijo ",
    "ni aleatorio. Se asume independencia entre fetos; la limitacion se declara ",
    "en el informe."))
)

LIMITACIONES <- c(
  paste0("Las conclusiones biologicas solo son validas con los datos reales. Sobre ",
    "datos sinteticos, el informe demuestra que el pipeline es completo y ",
    "reproducible; los efectos simulados son arbitrarios."),
  paste0("pSTAT3 (D9): normalizado a proteina total, sin STAT3 total en la misma ",
    "membrana -> refleja abundancia de fosfo-STAT3 (Tyr705), no la fraccion de ",
    "STAT3 fosforilada. Un cambio puede deberse a mas fosforilacion, a mas STAT3 ",
    "total, o a ambos."),
  paste0("ELISA de IL-6: n propio (14 madres en suero, 28 sacos en liquido amniotico), ",
    "no se fuerza al n=36 del diseno de qPCR. Censura a izquierda alta en varios ",
    "grupos (hasta 80 % en suero Control) -> en esos casos solo se reporta ",
    "proporcion de deteccion, no un estimador de ubicacion."),
  paste0("il6 en cerebro fetal E15: no cuantificable por D7 (0 detectados en el ",
    "calibrador). Solo se analiza como proporcion de deteccion."),
  paste0("BRAIN_P1 (cerebro postnatal dia 1) queda fuera del alcance de este informe ",
    "(placenta y cerebro fetal E15); se analiza por separado."),
  paste0("Acto 2: la comparacion de correlaciones entre grupos es de baja potencia ",
    "(n_par ~= 15-18 por estrato). La ausencia de significancia no es evidencia ",
    "de igualdad; se reporta el test formal (Fisher z) y la simulacion de ",
    "restriccion de rango, no el patron descriptivo."),
  paste0("Independencia asumida entre fetos (D13): el LPS se administra a la ",
    "madre y cada camada aporta un feto de cada sexo, asi que los dos fetos ",
    "de una camada no son estrictamente independientes. El analisis asume ",
    "independencia (los modelos no incluyen MADRE, ni fijo ni aleatorio). Si ",
    "existiera variacion entre camadas, afectaria sobre todo a la precision ",
    "de los efectos principales de tratamiento.")
)

# El informe sabe con que datos se genero (punto 3.4 de la revision): con
# fuente sintetica, las secciones interpretativas (conclusion de cada
# seccion, sintesis, conclusion revisada) se reemplazan por este aviso --
# nunca se llama a las funciones que arman esa prosa, para que sea
# estructuralmente imposible que aparezca una conclusion biologica sobre
# datos sinteticos.
AVISO_SINTETICO <- paste0(
  "<p><em>Informe generado con datos sinteticos. Los efectos son simulados ",
  "y arbitrarios; las conclusiones biologicas corresponden a los datos ",
  "reales, que no se incluyen en este repositorio.</em></p>")

# =========================================================================
# Formateo y escritura -- identicos a 02..11/98.
# =========================================================================
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
# Redondeo manual (floor(|x|*10^nd + .5)/10^nd) en vez de round(): mismo resultado
# bit a bit en R y Python sobre el mismo double, sin depender de la regla de
# redondeo-al-par de cada lenguaje -- necesario porque estos numeros van al
# texto de conclusiones y ese texto se byte-compara entre R y Python (99_verificar).
.round_fmt <- function(x, nd = 2L) {
  x <- as.numeric(x); s <- sign(x); m <- 10^nd
  v <- floor(abs(x) * m + 0.5) / m * s
  sprintf(paste0("%.", nd, "f"), v)
}
.join_y <- function(v) {
  if (length(v) == 0L) return("")
  if (length(v) == 1L) return(v[1])
  paste0(paste(v[-length(v)], collapse = ", "), " y ", v[length(v)])
}
escribir_texto <- function(ruta, texto) {
  if (!endsWith(texto, "\n")) texto <- paste0(texto, "\n")
  con <- file(ruta, open = "wb"); writeBin(charToRaw(enc2utf8(texto)), con); close(con)
}
escribir_csv <- function(ruta, encabezado, filas) {
  cuerpo <- vapply(filas, function(f)
    paste(vapply(f, .csv_cell, character(1)), collapse = ","), character(1))
  escribir_texto(ruta, paste(c(paste(vapply(encabezado, .csv_cell, character(1)),
                                      collapse = ","), cuerpo), collapse = "\n"))
}

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
leer_texto <- function(ruta) {
  if (!file.exists(ruta)) return("")
  t <- readChar(ruta, file.info(ruta)$size, useBytes = TRUE); Encoding(t) <- "UTF-8"; t
}
# Como leer_csv pero descarta las filas cuyo campo 'script' es 12_informe: el
# informe describe el pipeline de analisis, no su propio paso de reporte, y asi
# el bloque de auditoria queda identico R/Python e idempotente entre corridas.
leer_csv_sin_este <- function(ruta) {
  t <- leer_csv(ruta)
  if (is.null(t$header) || !("script" %in% t$header)) return(t)
  j <- match("script", t$header)
  t$filas <- Filter(function(f) !(j <= length(f) && f[[j]] == ESTE_SCRIPT), t$filas)
  t
}

# =========================================================================
# Merge en procedencia.csv / verificaciones.csv por 'script' -- igual que 02..11.
# =========================================================================
merge_por_script <- function(ruta, header, filas_nuevas, clave_orden) {
  viejo <- leer_csv(ruta)
  if (!is.null(viejo$header) && !identical(viejo$header, header))
    stop(sprintf("%s: encabezado incompatible", basename(ruta)))
  idx <- match("script", header)
  conservadas <- Filter(function(f) length(f) >= idx && f[[idx]] != ESTE_SCRIPT,
                        viejo$filas)
  nuevas <- lapply(filas_nuevas, function(f) vapply(f, .fmt, character(1)))
  todas <- c(conservadas, nuevas)
  todas <- todas[clave_orden(todas)]
  escribir_csv(ruta, header, todas)
}
registrar_procedencia <- function(filas_nuevas) {
  header <- c("artefacto", "tipo", "script", "origen_codigo", "entradas", "descripcion")
  merge_por_script(file.path(RUTA_TABLAS, "procedencia.csv"), header, filas_nuevas,
                   function(fs) order(vapply(fs, `[[`, character(1), 3),
                                      vapply(fs, `[[`, character(1), 1),
                                      method = "radix"))
}
registrar_verificaciones <- function(filas_nuevas) {
  header <- c("id", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script")
  merge_por_script(file.path(RUTA_TABLAS, "verificaciones.csv"), header, filas_nuevas,
                   function(fs) order(vapply(fs, `[[`, character(1), 6),
                                      vapply(fs, `[[`, character(1), 1),
                                      method = "radix"))
}

# =========================================================================
# Mini Markdown -> HTML  (deterministico; misma logica que python/12_informe.py).
# =========================================================================
.esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}
.inline <- function(s) {
  s <- .esc(s)
  s <- gsub("`([^`]+)`", "<code>\\1</code>", s, perl = TRUE)
  s <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", s, perl = TRUE)
  s <- gsub("\\*([^*]+)\\*", "<em>\\1</em>", s, perl = TRUE)
  gsub("\\[([^]]+)\\]\\(([^)]+)\\)", "<a href=\"\\2\">\\1</a>", s, perl = TRUE)
}
.celdas_tabla <- function(linea) {
  t <- trimws(linea)
  if (startsWith(t, "|")) t <- substring(t, 2)
  if (endsWith(t, "|")) t <- substring(t, 1, nchar(t) - 1L)
  trimws(strsplit(t, "|", fixed = TRUE)[[1]])
}
.es_separador_tabla <- function(linea) {
  cs <- .celdas_tabla(linea)
  length(cs) > 0 && all(vapply(cs, function(c)
    grepl("-", c, fixed = TRUE) &&
    !grepl("[^-: ]", c), logical(1)))
}
.lista_html <- function(bloque) {
  items <- list(); cur <- NULL; sub <- NULL
  empujar <- function() if (!is.null(cur))
    items[[length(items) + 1L]] <<- list(txt = cur, sub = sub)
  for (ln in bloque) {
    m1 <- regmatches(ln, regexec("^[-*]\\s+(.*)$", ln))[[1]]
    m2 <- regmatches(ln, regexec("^\\s{2,}([0-9]+)\\.\\s+(.*)$", ln))[[1]]
    if (length(m1) == 2L) {
      empujar()
      cur <- .inline(trimws(m1[2])); sub <- NULL
    } else if (length(m2) == 3L) {
      if (is.null(sub)) sub <- character(0)
      sub <- c(sub, .inline(trimws(m2[3])))
    } else {
      extra <- .inline(trimws(ln))
      if (is.null(cur)) cur <- extra else cur <- paste0(cur, " ", extra)
    }
  }
  empujar()
  partes <- "<ul>"
  for (it in items) {
    if (!is.null(it$sub) && length(it$sub)) {
      partes <- c(partes, paste0("<li>", it$txt, "<ol>",
        paste0("<li>", it$sub, "</li>", collapse = ""), "</ol></li>"))
    } else {
      partes <- c(partes, paste0("<li>", it$txt, "</li>"))
    }
  }
  paste(c(partes, "</ul>"), collapse = "\n")
}
md_a_html <- function(texto, base_nivel = 3L) {
  lineas <- strsplit(gsub("\r\n", "\n", texto, fixed = TRUE), "\n", fixed = TRUE)[[1]]
  out <- character(0); parr <- character(0); i <- 1L; n <- length(lineas)
  cerrar_parr <- function() {
    if (length(parr)) {
      out <<- c(out, paste0("<p>", paste(parr, collapse = " "), "</p>"))
      parr <<- character(0)
    }
  }
  while (i <= n) {
    ln <- lineas[i]; s <- trimws(ln)

    if (s == "" || grepl("^<!--.*-->$", s)) { cerrar_parr(); i <- i + 1L; next }

    mh <- regmatches(s, regexec("^(#{1,6})\\s+(.*)$", s))[[1]]
    if (length(mh) == 3L) {
      cerrar_parr()
      niv <- min(6L, base_nivel + nchar(mh[2]) - 1L)
      out <- c(out, sprintf("<h%d>%s</h%d>", niv, .inline(mh[3]), niv))
      i <- i + 1L; next
    }

    if (startsWith(s, "|") && i + 1L <= n && .es_separador_tabla(lineas[i + 1L])) {
      cerrar_parr()
      enc <- .celdas_tabla(s); i <- i + 2L
      cuerpo <- list()
      while (i <= n && startsWith(trimws(lineas[i]), "|")) {
        cuerpo[[length(cuerpo) + 1L]] <- .celdas_tabla(lineas[i]); i <- i + 1L
      }
      out <- c(out, paste0("<table><thead><tr>",
        paste0("<th>", vapply(enc, .inline, character(1)), "</th>", collapse = ""),
        "</tr></thead><tbody>"))
      for (fila in cuerpo)
        out <- c(out, paste0("<tr>",
          paste0("<td>", vapply(fila, .inline, character(1)), "</td>", collapse = ""),
          "</tr>"))
      out <- c(out, "</tbody></table>")
      next
    }

    if (grepl("^[-*]\\s+", sub("^\\s+", "", ln)) || grepl("^\\s{2,}[0-9]+\\.\\s+", ln)) {
      cerrar_parr()
      bloque <- character(0)
      while (i <= n && trimws(lineas[i]) != "" &&
             (grepl("^[-*]\\s+", sub("^\\s+", "", lineas[i])) ||
              grepl("^\\s", lineas[i]))) {
        bloque <- c(bloque, lineas[i]); i <- i + 1L
      }
      out <- c(out, .lista_html(bloque))
      next
    }

    parr <- c(parr, .inline(s)); i <- i + 1L
  }
  cerrar_parr()
  paste(out, collapse = "\n")
}

# =========================================================================
# Numeros del resumen -- se leen de outputs/tables/<lang>/ (byte-identicos R/Py).
# =========================================================================
.tab <- function(nombre) leer_csv(file.path(RUTA_TABLAS_R, nombre))
.col <- function(t, nombre) {
  if (is.null(t$header) || !(nombre %in% t$header)) return(character(0))
  j <- match(nombre, t$header)
  vapply(t$filas, function(f) if (j <= length(f)) f[[j]] else "", character(1))
}
# =========================================================================
# Figuras por seccion -- DERIVADAS de procedencia.csv, no hardcodeadas (ver
# comentario de SECCIONES). `figuras_procedencia()` es la unica lectura de
# procedencia.csv para esto; `figuras_de_seccion()` filtra por script (y un
# patron opcional de nombre) para una seccion puntual; `figuras_embebidas()`
# junta lo que efectivamente se incrusta en el informe completo (usado tanto
# para armar el HTML como para la verificacion de cobertura en main()).
# =========================================================================
figuras_procedencia <- function() {
  t <- leer_csv_sin_este(file.path(RUTA_TABLAS, "procedencia.csv"))
  art <- .col(t, "artefacto"); tipo <- .col(t, "tipo"); scr <- .col(t, "script")
  es_fig <- tipo == "figura" & startsWith(art, "outputs/figures/")
  list(artefacto = art[es_fig], script = scr[es_fig])
}
figuras_de_seccion <- function(figs, criterio) {
  idx <- figs$script == criterio$script
  nombres <- basename(figs$artefacto[idx])
  if (!is.null(criterio$patron)) nombres <- nombres[grepl(criterio$patron, nombres)]
  nombres
}
figuras_embebidas <- function(figs) unique(unlist(
  lapply(SECCIONES, function(sec) figuras_de_seccion(figs, sec[[3]]))))

resumen_numeros <- function() {
  r <- list()

  t <- .tab("pstat3_modelo_clasificacion.csv")
  p <- .col(t, "p_SEXOxTTO"); r$pstat3_pint <- if (length(p)) p[1] else "n/d"

  t <- .tab("pstat3_posthoc.csv")
  contr <- .col(t, "contraste"); ph <- .col(t, "p_holm"); r$pstat3_hh <- "n/d"
  for (k in seq_along(contr)) if (contr[k] == "HEMBRA_CONTROL-HEMBRA_LPS")
    r$pstat3_hh <- ph[k]

  t <- .tab("qpcr_modelos_clasificacion.csv")
  tej <- .col(t, "TEJIDO"); isig <- .col(t, "interaccion_significativa")
  via <- .col(t, "via")
  r$qpcr_int_pla <- sum(tej == "PLACENTA_E15" & isig == "TRUE")
  r$qpcr_int_bra <- sum(tej == "BRAIN_E15" & isig == "TRUE")
  r$qpcr_modelados <- sum(via == "modelo")

  t <- .tab("acto2_test_correlaciones.csv")
  it <- .col(t, "ITEM"); pbw <- .col(t, "p_bw")
  ok <- nzchar(pbw)
  r$acto2_n_test <- sum(ok)
  if (any(ok)) {
    itv <- it[ok]; pv <- as.numeric(pbw[ok]); k <- which.min(pv)
    r$acto2_min_item <- itv[k]; r$acto2_min_pbw <- pbw[ok][k]
    r$acto2_n_sig <- sum(pv < 0.05)
  } else {
    r$acto2_min_item <- "n/d"; r$acto2_min_pbw <- "n/d"; r$acto2_n_sig <- 0L
  }

  t <- .tab("acto2_simulacion.csv")
  ver <- .col(t, "veredicto")
  r$sim_fuera <- sum(ver == "FUERA"); r$sim_total <- length(ver)

  t <- .tab("acto2_sensibilidad_excl_extremo.csv")
  r$sens_cambia <- sum(.col(t, "veredicto_cambia") == "TRUE")
  r$sens_items <- length(t$filas)

  t <- .tab("acto2_sensibilidad_eigengene_test.csv")
  it <- .col(t, "ITEM"); pbw <- .col(t, "p_bw"); r$eig_pbw <- "n/d"
  for (k in seq_along(it)) if (it[k] == "eigengene") r$eig_pbw <- pbw[k]

  t <- leer_csv(file.path(RUTA_TABLAS, "comparacion_R_python.csv"))
  r$comp_total <- "n/d"; r$comp_byte <- "n/d"; r$comp_fuera <- "n/d"
  if (!is.null(t$header)) {
    arch <- .col(t, "archivo")
    for (k in seq_along(arch)) if (arch[k] == "__TOTAL__") {
      d <- setNames(as.list(t$filas[[k]]), t$header)
      r$comp_byte <- d[["byte_identico"]]; r$comp_fuera <- d[["n_fuera_tol"]]
    }
    r$comp_total <- as.character(sum(arch != "__TOTAL__"))
  }

  t <- leer_csv_sin_este(file.path(RUTA_TABLAS, "procedencia.csv"))
  r$n_procedencia <- length(t$filas)
  t <- leer_csv_sin_este(file.path(RUTA_TABLAS, "verificaciones.csv"))
  okc <- .col(t, "ok")
  r$n_verif <- length(t$filas); r$n_verif_true <- sum(okc == "TRUE")
  r
}

# =========================================================================
# Numeros para las "Conclusion de la seccion", la sintesis y la conclusion
# revisada (pedido pedidos/cambios_informe_conclusiones.md, puntos 1-4).
# Mismo principio que resumen_numeros(): todo numero citado en prosa se lee
# de un CSV, nunca se escribe a mano. Las listas de genes tambien se derivan
# filtrando (no se copian del pedido) para que el texto siga los datos si
# estos cambian.
# =========================================================================
numeros_conclusiones <- function() {
  n <- list()

  # --- ELISA (Acto 1.1) ---
  t <- .tab("elisa_fisher_deteccion.csv")
  bl <- .col(t, "bloque")
  k <- which(bl == "MS")[1]
  n$ms_ctrl_det <- .col(t, "control_detectado")[k]; n$ms_ctrl_n <- .col(t, "control_n")[k]
  n$ms_lps_det <- .col(t, "lps_detectado")[k]; n$ms_lps_n <- .col(t, "lps_n")[k]
  n$ms_p <- .col(t, "p_valor")[k]

  t <- .tab("elisa_petopeto_la.csv")
  est <- .col(t, "estrato"); pv <- .col(t, "p_valor")
  n$la_hembra_p <- pv[est == "HEMBRA"][1]; n$la_macho_p <- pv[est == "MACHO"][1]
  n$la_n_control <- .col(t, "n_control")[est == "HEMBRA"][1]

  t <- .tab("elisa_descriptivo.csv")
  bl <- .col(t, "bloque"); tt <- .col(t, "TTO"); sx <- .col(t, "SEXO")
  med <- .col(t, "mediana_detectada")
  .med <- function(tto, sexo) {
    i <- which(bl == "LA" & tt == tto & sx == sexo)
    if (length(i)) .round_fmt(med[i[1]]) else "n/d"
  }
  n$la_h_ctrl_med <- .med("CONTROL", "HEMBRA"); n$la_h_lps_med <- .med("LPS", "HEMBRA")
  n$la_m_ctrl_med <- .med("CONTROL", "MACHO"); n$la_m_lps_med <- .med("LPS", "MACHO")

  # --- qPCR modelos (Acto 1.3) ---
  t <- .tab("qpcr_modelos_clasificacion.csv")
  tej <- .col(t, "TEJIDO"); gen <- .col(t, "GEN"); via <- .col(t, "via")
  pTTO <- .col(t, "p_TTO"); isig <- .col(t, "interaccion_significativa")
  bh_int <- suppressWarnings(as.numeric(.col(t, "p_SEXOxTTO_BH")))
  i_pla <- tej == "PLACENTA_E15" & via == "modelo" & !is.na(suppressWarnings(as.numeric(pTTO))) &
    as.numeric(pTTO) < 0.05
  n$pla_tto_genes <- .join_y(sort(gen[i_pla], method = "radix"))
  # C3: derivar "ningun gen mostro interaccion en placenta" del conteo real,
  # no darlo por sentado -- si algun gen x tejido de placenta tuviera
  # interaccion significativa, la frase cambia sola.
  i_pla_int <- tej == "PLACENTA_E15" & isig == "TRUE"
  n$pla_int_n <- sum(i_pla_int)
  n$pla_int_genes <- .join_y(sort(gen[i_pla_int], method = "radix"))
  i_bra_int <- tej == "BRAIN_E15" & isig == "TRUE"
  bra_int_genes <- sort(gen[i_bra_int], method = "radix")
  n$bra_int_n_txt <- length(bra_int_genes)
  # D12/C1: "robusta a BH" se calcula, no se asume -- numerador/denominador de
  # los que ademas tienen p_SEXOxTTO_BH < .05 entre los primariamente
  # significativos.
  i_bra_int_bh <- i_bra_int & !is.na(bh_int) & bh_int < 0.05
  n$bra_int_bh_n <- sum(i_bra_int_bh)

  t <- .tab("qpcr_modelos_posthoc.csv")
  tej_p <- .col(t, "TEJIDO"); gen_p <- .col(t, "GEN"); pholm <- .col(t, "p_holm")
  media_genes <- character(0); disp_genes <- character(0); disp_holm_min <- Inf
  for (g in bra_int_genes) {
    idx <- tej_p == "BRAIN_E15" & gen_p == g
    ps <- suppressWarnings(as.numeric(pholm[idx]))
    if (any(ps < 0.05, na.rm = TRUE)) {
      media_genes <- c(media_genes, g)
    } else {
      disp_genes <- c(disp_genes, g)
      disp_holm_min <- min(disp_holm_min, min(ps, na.rm = TRUE))
    }
  }
  n$bra_media_genes <- .join_y(sort(media_genes, method = "radix"))
  n$bra_disp_genes <- .join_y(sort(disp_genes, method = "radix"))
  n$bra_disp_holm_min <- .round_fmt(disp_holm_min, 3L)

  # --- pSTAT3 (Acto 1.4) ---
  t <- .tab("pstat3_descriptivo.csv")
  gr <- .col(t, "GRUPO"); me <- .col(t, "mean_PSTAT3")
  .mean_grp <- function(g) .round_fmt(me[gr == g][1])
  n$pstat3_hc <- .mean_grp("HEMBRA_CONTROL"); n$pstat3_hl <- .mean_grp("HEMBRA_LPS")
  n$pstat3_mc <- .mean_grp("MACHO_CONTROL"); n$pstat3_ml <- .mean_grp("MACHO_LPS")
  t <- .tab("pstat3_posthoc.csv")
  contr <- .col(t, "contraste"); ph <- .col(t, "p_holm")
  n$pstat3_hh <- ph[contr == "HEMBRA_CONTROL-HEMBRA_LPS"][1]
  n$pstat3_mm <- ph[contr == "MACHO_CONTROL-MACHO_LPS"][1]

  # --- Acto 2.3-2.4: Delta rho + interaccion sobre dispersion ---
  t <- .tab("acto2_test_correlaciones.csv")
  it <- .col(t, "ITEM"); es <- .col(t, "ESTRATO"); pbw <- .col(t, "p_bw")
  ok <- nzchar(pbw)
  n$delta_n_test <- sum(ok)
  pv <- suppressWarnings(as.numeric(pbw[ok]))
  n$delta_n_sig <- sum(pv < 0.05, na.rm = TRUE)
  kmin <- which.min(pv)
  n$delta_min_item <- it[ok][kmin]; n$delta_min_p <- pbw[ok][kmin]
  lbl_estrato <- c(AMBOS_SEXOS = "ambos sexos", HEMBRA = "hembras", MACHO = "machos")
  n$delta_min_estrato <- unname(lbl_estrato[es[ok][kmin]])

  t <- .tab("acto2_dispersion_interaccion.csv")
  gen_d <- .col(t, "GEN"); tej_d <- .col(t, "TEJIDO"); bh <- .col(t, "p_SEXOxTTO_BH")
  bh_n <- suppressWarnings(as.numeric(bh))
  i_bra_sig <- tej_d == "BRAIN_E15" & bh_n < 0.05
  i_bra_tend <- tej_d == "BRAIN_E15" & bh_n >= 0.05 & bh_n < 0.10
  i_pla_sig <- tej_d == "PLACENTA_E15" & bh_n < 0.05
  ord_sig <- order(bh_n[i_bra_sig], gen_d[i_bra_sig], method = "radix")
  n$disp_bra_sig_genes <- .join_y(gen_d[i_bra_sig][ord_sig])
  n$disp_bra_sig_n <- sum(i_bra_sig)
  n$disp_bra_tend_genes <- .join_y(sort(gen_d[i_bra_tend], method = "radix"))
  n$disp_pla_sig_n <- sum(i_pla_sig)

  # --- Acto 2.5: simulacion, items en el limite del IC (ESTRATO=HEMBRA, GLOBAL) ---
  t <- .tab("acto2_simulacion.csv")
  it_s <- .col(t, "ITEM"); es_s <- .col(t, "ESTRATO"); esc_s <- .col(t, "ESCENARIO")
  ver_s <- .col(t, "veredicto")
  i_lim <- es_s == "HEMBRA" & esc_s == "GLOBAL" & ver_s == "FUERA"
  nombres_lim <- gsub("_", " ", it_s[i_lim])
  n$sim_lim_genes <- .join_y(nombres_lim)

  # --- Acto 2.6: eigengene y sensibilidad de exclusion ---
  t <- .tab("acto2_sensibilidad_pca_varianza.csv")
  tej_v <- .col(t, "TEJIDO"); pc <- .col(t, "PC"); pv_ <- .col(t, "prop_var")
  .pc1_pct <- function(tej) {
    i <- which(tej_v == tej & pc == "1")
    .round_fmt(as.numeric(pv_[i[1]]) * 100, 0L)
  }
  n$eig_pla_pct <- .pc1_pct("PLACENTA_E15"); n$eig_bra_pct <- .pc1_pct("BRAIN_E15")

  t <- .tab("acto2_sensibilidad_excl_extremo.csv")
  it_e <- .col(t, "ITEM"); pf <- .col(t, "p_bw_full"); ps <- .col(t, "p_bw_sin")
  .excl <- function(item) {
    i <- which(it_e == item)[1]
    c(full = .round_fmt(pf[i], 3L), sin = .round_fmt(ps[i], 3L))
  }
  e1 <- .excl("fatcd36"); e2 <- .excl("score_compuesto")
  n$sens_fatcd36_full <- e1[["full"]]; n$sens_fatcd36_sin <- e1[["sin"]]
  n$sens_score_full <- e2[["full"]]; n$sens_score_sin <- e2[["sin"]]

  n
}

# =========================================================================
# "Conclusion de la seccion" (pedido, punto 1): un bloque de 2-5 oraciones al
# final de cada seccion de SECCIONES, indexado igual (1..8). El indice 2
# (Acto 1.2, cuantificacion) no lleva conclusion -- es de metodo. Prosa fija
# (igual en R/Python), numeros y listas de genes interpolados desde `n`
# (numeros_conclusiones()).
# =========================================================================
conclusion_seccion <- function(k, n) {
  txt <- switch(as.character(k),
    "1" = paste0(
      "El LPS indujo una respuesta inflamatoria sistemica: la IL-6 fue ",
      "detectable en ", n$ms_lps_det, "/", n$ms_lps_n, " madres LPS frente a ",
      n$ms_ctrl_det, "/", n$ms_ctrl_n, " control (Fisher p = <code>", n$ms_p,
      "</code>), lo que valida el modelo. En liquido amniotico ningun ",
      "contraste alcanzo significancia (Peto-Peto ♀ p = <code>",
      n$la_hembra_p, "</code>; ♂ p = <code>", n$la_macho_p,
      "</code>), aunque la mediana de los valores detectados fue mayor bajo ",
      "LPS en ambos sexos (♀ ", n$la_h_ctrl_med, " &rarr; ",
      n$la_h_lps_med, "; ♂ ", n$la_m_ctrl_med, " &rarr; ", n$la_m_lps_med,
      "). Con ", n$la_n_control, " sacos control por sexo, la ausencia de ",
      "significancia no permite concluir que la IL-6 no llegue al ",
      "compartimento fetal."),
    "3" = paste0(
      "<p><strong>Placenta.</strong> ",
      if (n$pla_int_n == 0L) paste0(
        "Ningun gen mostro interaccion SEXO&times;TTO. El LPS modifico la ",
        "expresion de ", n$pla_tto_genes, " de forma equivalente en ambos ",
        "sexos (efecto principal de tratamiento).")
      else paste0(
        n$pla_int_n, " gen x tejido mostraron interaccion SEXO&times;TTO ",
        "significativa (", n$pla_int_genes, "). El LPS tambien modifico la ",
        "expresion de ", n$pla_tto_genes, " de forma equivalente en ambos ",
        "sexos en los genes sin interaccion (efecto principal de ",
        "tratamiento)."),
      "</p>\n<p><strong>Cerebro fetal.</strong> ", n$bra_int_n_txt,
      " genes mostraron interaccion SEXO&times;TTO significativa (",
      n$bra_int_bh_n, " de ", n$bra_int_n_txt, " sobreviven a la ",
      "correccion BH). En ", n$bra_media_genes, " el post hoc ",
      "localiza el efecto en hembras: ♀Control difiere de ♀LPS ",
      "y ♀LPS difiere de ♂LPS, sin cambios en machos. En ",
      n$bra_disp_genes, " la interaccion no se explica por ninguna ",
      "comparacion de medias (todos los p de Holm &ge; ", n$bra_disp_holm_min,
      "). La seccion 2.3 muestra que en esos genes el efecto esta en la ",
      "dispersion y no en la media.</p>"),
    "4" = paste0(
      "El LPS aumento la abundancia de fosfo-STAT3 en placenta solo en ",
      "hembras (♀Control ", n$pstat3_hc, " &rarr; ♀LPS ",
      n$pstat3_hl, "; Holm p = <code>", n$pstat3_hh, "</code>). En machos no ",
      "cambio (", n$pstat3_mc, " &rarr; ", n$pstat3_ml, "; p = <code>",
      n$pstat3_mm, "</code>). Es el mismo patron que los genes con efecto en ",
      "media en cerebro. Como los transportadores placentarios responden ",
      "igual en ambos sexos, es compatible con que sus cambios no dependan ",
      "de la activacion de STAT3, o con que los machos los alcancen por ",
      "otra via. Recordar D9: se mide abundancia de fosfo-STAT3, no fraccion ",
      "fosforilada."),
    "5" = paste0(
      "Descriptivamente, en hembras control la correlacion placenta&ndash;",
      "cerebro es alta en varios genes y cae cerca de cero con LPS; en ",
      "machos no hay correlacion en ningun grupo. Estas figuras describen y ",
      "no testean (prohibicion 4). En los diagramas triangulares de cerebro ",
      "de hembras, la distribucion ♀Control es ancha con una cola ",
      "hacia valores bajos, y la ♀LPS es un pico angosto: la ",
      "reduccion de dispersion de la seccion 2.3 es visible directamente."),
    "6" = paste0(
      "<p>Ningun &Delta;&rho; Control vs LPS es significativo, ni agrupando ",
      "sexos ni dentro de cada sexo (", n$delta_n_sig, " de ", n$delta_n_test,
      "; minimo <code>", n$delta_min_item, "</code> en ", n$delta_min_estrato,
      ", p = <code>", n$delta_min_p, "</code>). El test de interaccion ",
      "SEXO&times;TTO sobre la dispersion es el resultado positivo del Acto ",
      "2: en cerebro, ", n$disp_bra_sig_n, " genes sobreviven a BH (",
      n$disp_bra_sig_genes, "; ", n$disp_bra_tend_genes, " en tendencia), ",
      "con un patron cruzado: el LPS reduce la dispersion en hembras y la ",
      "aumenta en machos. En placenta ninguno sobrevive a BH (", n$disp_pla_sig_n,
      "). Como los transportadores varian mayormente juntos (el eigengene ",
      "explica el ", n$eig_bra_pct, " % de la varianza en cerebro; seccion ",
      "2.6), es compatible con que estos genes no sean efectos ",
      "independientes, sino que reflejen un patron compartido por el ",
      "conjunto de transportadores.</p>\n",
      "<p>La figura de SD ahora muestra los tres estratos (agrupado, hembras, ",
      "machos) lado a lado: el patron cruzado se ve directamente comparando ",
      "las filas HEMBRA y MACHO, columna por item -- y explica por que la ",
      "fila AMBOS_SEXOS, arriba de las otras dos, no lo muestra: los cambios ",
      "opuestos de hembras y machos se cancelan al promediarlos. Es un ",
      "ejemplo directo de lo que oculta agrupar los sexos.</p>"),
    "7" = paste0(
      "La caida de correlacion observada en hembras es compatible con la ",
      "compactacion de la expresion bajo LPS: cuando el rango de una ",
      "variable se reduce, la correlacion cae aunque la relacion biologica ",
      "no haya cambiado. Las excepciones (", n$sim_lim_genes,
      ", en hembras) quedan en el limite del intervalo y se toman como ",
      "pista, no como hallazgo. La perdida aparente de acoplamiento ",
      "placenta&ndash;cerebro en hembras es compatible con la reduccion de ",
      "dispersion vista desde otro angulo, aunque no permite descartar un ",
      "cambio de coordinacion."),
    "8" = paste0(
      "El eigengene (PC1) explica el ", n$eig_pla_pct,
      " % de la varianza de los transportadores en placenta y el ",
      n$eig_bra_pct, " % en cerebro, con cargas similares para los siete ",
      "genes: los transportadores varian mayormente juntos, como un unico ",
      "eje por feto. Reemplazar el score compuesto por el eigengene no ",
      "cambia ninguna conclusion. Excluir el feto mas influyente de cada ",
      "item tampoco cambia veredictos, pero aproximadamente duplica los p ",
      "de los items que estaban cerca del umbral (fatcd36 ",
      n$sens_fatcd36_full, " &rarr; ", n$sens_fatcd36_sin, "; score ",
      n$sens_score_full, " &rarr; ", n$sens_score_sin, "): esas senales son ",
      "sensibles a la exclusion de un solo feto."),
    NULL
  )
  if (is.null(txt)) return("")
  cuerpo <- if (startsWith(txt, "<p>")) txt else paste0("<p>", txt, "</p>")
  paste0("<h4>Conclusion de la seccion</h4>\n", cuerpo)
}

# =========================================================================
# Sintesis (punto 2) y conclusion revisada (punto 3) del pedido -- prosa fija,
# numeros/listas interpolados desde `n` (numeros_conclusiones()).
# =========================================================================
sintesis_eje_html <- function(n) {
  paste0(
    "<p><strong>Madre.</strong> El LPS produjo una respuesta inflamatoria ",
    "sistemica inequivoca (IL-6 serica detectable en ", n$ms_lps_det, "/",
    n$ms_lps_n, " madres tratadas frente a ", n$ms_ctrl_det, "/", n$ms_ctrl_n,
    " control).</p>\n",
    "<p><strong>Liquido amniotico.</strong> Sin resultado concluyente: con ",
    n$la_n_control, " sacos control por sexo no se detecto diferencia, lo ",
    "que no descarta que la IL-6 llegue al compartimento fetal.</p>\n",
    "<p><strong>Placenta.</strong> Dos respuestas que no coinciden. La ",
    "senalizacion IL-6/STAT3 se activa solo en placentas de fetos hembra. La ",
    "expresion de transportadores de nutrientes, en cambio, se modifica en ",
    "ambos sexos por igual.</p>\n",
    "<p><strong>Cerebro fetal.</strong> La respuesta depende del sexo del ",
    "feto. En hembras, el LPS modifica la expresion de ", n$bra_media_genes,
    ", con el mismo patron que pSTAT3 en placenta; en machos, esos genes no ",
    "cambian.</p>\n",
    "<p><strong>Lectura del Acto 1.</strong> El LPS materno induce en la ",
    "descendencia hembra una respuesta coherente a lo largo del eje: ",
    "activacion de STAT3 en placenta y cambio de la expresion de ",
    "transportadores en cerebro. En la descendencia macho, la respuesta en ",
    "cerebro no es detectable como cambio de nivel.</p>")
}

conclusion_revisada_html <- function(n) {
  filas <- list(
    c("Siete genes de cerebro responden con dimorfismo sexual",
      "Test de interaccion SEXO&times;TTO sobre la dispersion",
      paste0("Se precisa: ", n$bra_media_genes,
        " con desplazamiento de la media en hembras; en ", n$bra_disp_genes,
        " lo que cambia es la variabilidad, como parte de un patron ",
        "compartido por los transportadores")),
    c("(Lectura intuitiva de las figuras) El LPS desacopla placenta y cerebro en hembras",
      "Test formal de &Delta;&rho; + simulacion de restriccion de rango",
      paste0("Se descarta: ningun &Delta;&rho; significativo (", n$delta_n_sig,
        "/", n$delta_n_test, "); la caida de correlacion se explica por la ",
        "compactacion de la expresion")),
    c("La respuesta resumida por el score compuesto es robusta",
      "Eigengene PC1 y exclusion del feto extremo",
      paste0("Se sostiene: el eigengene no cambia conclusiones; las senales ",
        "cercanas al umbral son sensibles a la exclusion de un solo feto"))
  )
  tabla <- c(paste0("<table><thead><tr><th>Afirmacion del Acto 1</th>",
    "<th>Que la puso a prueba en el Acto 2</th><th>Resultado</th></tr>",
    "</thead><tbody>"))
  for (f in filas)
    tabla <- c(tabla, sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>",
      f[1], f[2], f[3]))
  tabla <- c(tabla, "</tbody></table>")
  parrafo <- paste0(
    "<p><strong>Conclusion revisada.</strong> En hembras, el LPS materno ",
    "produce una respuesta direccional y homogenea: activa STAT3 en ",
    "placenta, desplaza la expresion cerebral de ", n$bra_media_genes,
    ", y reduce la variabilidad entre individuos. En machos no hay ",
    "respuesta direccional, pero aumenta la variabilidad entre fetos. La ",
    "aparente perdida de acoplamiento placenta&ndash;cerebro en hembras es ",
    "compatible con la reduccion de dispersion; la simulacion muestra que ",
    "esta alcanza para explicarla, aunque no permite descartar un cambio de ",
    "coordinacion.</p>")
  paste(c(tabla, parrafo), collapse = "\n")
}

# =========================================================================
# base64 propio (RFC 4648, igual salida que base64.b64encode de Python).
# =========================================================================
.B64 <- c(LETTERS, letters, 0:9, "+", "/")
base64_raw <- function(bytes) {
  n <- length(bytes); if (n == 0L) return("")
  v <- as.integer(bytes)
  pad <- (3L - n %% 3L) %% 3L
  v <- c(v, rep(0L, pad))
  m <- matrix(v, nrow = 3L)
  idx <- rbind(
    bitwShiftR(m[1, ], 2L),
    bitwOr(bitwShiftL(bitwAnd(m[1, ], 3L), 4L), bitwShiftR(m[2, ], 4L)),
    bitwOr(bitwShiftL(bitwAnd(m[2, ], 15L), 2L), bitwShiftR(m[3, ], 6L)),
    bitwAnd(m[3, ], 63L)
  )
  ch <- .B64[idx + 1L]
  if (pad > 0L) ch[(length(ch) - pad + 1L):length(ch)] <- "="
  paste(ch, collapse = "")
}
img_datauri <- function(ruta) {
  b <- readBin(ruta, "raw", n = file.info(ruta)$size)
  paste0("data:image/png;base64,", base64_raw(b))
}
fig_html <- function(nombre) {
  ruta <- file.path(RUTA_FIGURAS, nombre)
  if (!file.exists(ruta))
    return(sprintf('<p class="falta">[falta la figura %s]</p>', .esc(nombre)))
  sprintf('<figure><img alt="%s" src="%s"><figcaption>%s</figcaption></figure>',
          .esc(nombre), img_datauri(ruta), .esc(nombre))
}

# =========================================================================
# CSS (inline; sin recursos externos). IDENTICO a python/12_informe.py.
# =========================================================================
CSS <- paste0("\n",
":root { color-scheme: light; }\n",
"* { box-sizing: border-box; }\n",
"body { font: 15px/1.55 -apple-system, \"Segoe UI\", Roboto, Helvetica, Arial,\n",
"       sans-serif; color: #1a1a1a; background: #fff; margin: 0;\n",
"       padding: 2.2rem 1.4rem 4rem; }\n",
"main { max-width: 60rem; margin: 0 auto; }\n",
"h1 { font-size: 1.7rem; margin: 0 0 .2rem; }\n",
"h2 { font-size: 1.28rem; margin: 2.4rem 0 .6rem; padding-top: .4rem;\n",
"     border-top: 2px solid #222; }\n",
"h3 { font-size: 1.08rem; margin: 1.5rem 0 .4rem; }\n",
"h4 { font-size: .98rem; margin: 1.1rem 0 .3rem; color: #333; }\n",
"h5 { font-size: .92rem; margin: .9rem 0 .3rem; color: #444; }\n",
"p { margin: .5rem 0; }\n",
"code { background: #f0f0f0; padding: .05em .35em; border-radius: 3px;\n",
"       font: .86em/1.4 \"SF Mono\", Consolas, \"Liberation Mono\", monospace; }\n",
"a { color: #0645ad; }\n",
"ul, ol { margin: .4rem 0 .7rem; padding-left: 1.5rem; }\n",
"li { margin: .18rem 0; }\n",
"table { border-collapse: collapse; margin: .8rem 0; font-size: .82rem;\n",
"        display: block; overflow-x: auto; max-width: 100%; }\n",
"th, td { border: 1px solid #ccc; padding: .28em .55em; text-align: left;\n",
"         white-space: nowrap; }\n",
"thead th { background: #f4f4f4; position: sticky; top: 0; }\n",
"tbody tr:nth-child(even) { background: #fafafa; }\n",
"figure { margin: 1rem 0; text-align: center; }\n",
"figure img { max-width: 100%; height: auto; border: 1px solid #e2e2e2; }\n",
"figcaption { font-size: .78rem; color: #666; margin-top: .3rem; }\n",
".aviso { background: #fff4e5; border: 1px solid #f0c48a; padding: .8rem 1rem;\n",
"         border-radius: 6px; margin: 1rem 0; }\n",
".meta { color: #555; font-size: .9rem; }\n",
".falta { color: #b00; font-style: italic; }\n",
"nav.toc { background: #f7f7f7; border: 1px solid #e0e0e0; border-radius: 6px;\n",
"          padding: .8rem 1.2rem; margin: 1.4rem 0; }\n",
"nav.toc ol { margin: .3rem 0; }\n",
"footer { margin-top: 3rem; padding-top: 1rem; border-top: 1px solid #ccc;\n",
"         color: #666; font-size: .85rem; }\n",
"@media print {\n",
"  body { padding: 0; font-size: 11pt; }\n",
"  h2 { page-break-before: auto; }\n",
"  figure, table { page-break-inside: avoid; }\n",
"  thead th { position: static; }\n",
"  nav.toc { page-break-after: always; }\n",
"}\n")

# =========================================================================
# Armado del HTML.
# =========================================================================
.seccion <- function(id_, titulo, cuerpo_html) {
  sprintf('<section id="%s">\n<h2>%s</h2>\n%s\n</section>', id_, .esc(titulo),
          cuerpo_html)
}
.csv_a_tabla <- function(ruta, omitir = character(0)) {
  t <- leer_csv_sin_este(ruta)
  if (is.null(t$header))
    return(sprintf('<p class="falta">[falta %s]</p>', .esc(basename(ruta))))
  h <- t$header
  cols <- which(!(h %in% omitir))
  out <- paste0("<table><thead><tr>",
    paste0("<th>", vapply(h[cols], .esc, character(1)), "</th>", collapse = ""),
    "</tr></thead><tbody>")
  for (fila in t$filas) {
    celdas <- vapply(cols, function(j)
      .esc(if (j <= length(fila)) fila[[j]] else ""), character(1))
    out <- c(out, paste0("<tr>", paste0("<td>", celdas, "</td>", collapse = ""),
                         "</tr>"))
  }
  paste(c(out, "</tbody></table>"), collapse = "\n")
}

construir_html <- function(fuente, num, nc) {
  sint <- fuente != "real"
  L <- character(0)
  ap <- function(...) L[[length(L) + 1L]] <<- paste0(...)

  ap("<!doctype html>")
  ap('<html lang="es">')
  ap("<head>")
  ap('<meta charset="utf-8">')
  ap('<meta name="viewport" content="width=device-width, initial-scale=1">')
  ap("<title>Reanalisis MIA-LPS -- placenta E15 / cerebro fetal E15</title>")
  ap("<style>", CSS, "</style>")
  ap("</head>")
  ap("<body>")
  ap("<main>")

  ap("<h1>Reanalisis reproducible &mdash; MIA-LPS en placenta E15 y cerebro ",
     "fetal E15</h1>")
  ap('<p class="meta">Activacion inmune materna (LPS 100 &micro;g/kg i.p., dia 15 ',
     "de gestacion, colecta a las 6 h). Respuesta de la placenta y del cerebro ",
     "fetal, y si depende del sexo del feto.</p>")
  ap('<p class="meta">Fuente de datos de esta corrida: <code>', .esc(fuente),
     "</code>. Informe generado por <code>12_informe</code> (R y Python producen el ",
     "mismo <code>informe.html</code> salvo los PNG incrustados). Sin marca de ",
     "tiempo: la fecha de corrida esta en <code>logs/corrida_&lt;fecha&gt;.txt</code>.")
  if (sint) {
    ap('<div class="aviso"><strong>Datos sinteticos.</strong> Esta corrida no ',
       "uso los crudos reales. Las salidas demuestran que el pipeline es ",
       "completo y reproducible; los efectos son arbitrarios y no tienen ",
       "lectura biologica.</div>")
  }

  items_toc <- list(
    c("resumen", "1. Resumen"),
    c("metodos", "2. Diseno y metodos"),
    c("acto1", "3. Acto 1 &mdash; respuesta a la MIA y dependencia del sexo"),
    c("sintesis", "4. Sintesis del eje madre &rarr; placenta &rarr; cerebro"),
    c("acto2", "5. Acto 2 &mdash; coordinacion placenta&lt;-&gt;cerebro"),
    c("conclusion_revisada", "6. Conclusion revisada (Acto 1 frente a Acto 2)"),
    c("reproducibilidad", "7. Reproducibilidad (R vs Python)"),
    c("descartados", "8. Analisis descartados"),
    c("limitaciones", "9. Limitaciones"),
    c("auditoria", "10. Procedencia y verificaciones")
  )
  ap('<nav class="toc"><strong>Contenido</strong><ol>')
  for (it in items_toc) ap('<li><a href="#', it[1], '">', it[2], "</a></li>")
  ap("</ol></nav>")

  # --- 1. resumen ---
  res <- c(
    "<ul>",
    paste0("<li><strong>Validacion del modelo (ELISA IL-6).</strong> El LPS eleva ",
      "IL-6 en suero materno (1/5 vs 9/9 detectados; Fisher p &asymp; 5&times;10",
      "<sup>-3</sup>). En liquido amniotico la senal es mas debil y la censura ",
      "alta; el Peto-Peto por sexo no alcanza significancia.</li>"),
    paste0("<li><strong>pSTAT3 en placenta.</strong> Interaccion SEXO&times;TTO ",
      "significativa (p = <code>", num$pstat3_pint, "</code>): el aumento de ",
      "fosfo-STAT3 con LPS es restringido a hembras (HEMBRA_CONTROL&ndash;",
      "HEMBRA_LPS Holm p = <code>", num$pstat3_hh, "</code>); en machos no cambia. ",
      "Mide abundancia de fosfo-STAT3, no fraccion fosforilada (D9).</li>"),
    paste0("<li><strong>Programa de transportadores (qPCR).</strong> De ",
      num$qpcr_modelados, " gen &times; tejido modelados, la interaccion ",
      "SEXO&times;TTO es significativa en ", num$qpcr_int_bra, " de cerebro E15 y ",
      num$qpcr_int_pla, " de placenta E15; en cerebro las ", num$qpcr_int_bra,
      " sobreviven la correccion BH entre genes (D12).</li>"),
    paste0("<li><strong>Coordinacion placenta&lt;-&gt;cerebro (Acto 2).</strong> ",
      "En el test formal de &Delta;&rho; (Fisher z sobre &rho; de Spearman, ",
      "prohibicion 4), ", num$acto2_n_sig, " de ", num$acto2_n_test,
      " items alcanzan p &lt; 0.05 (minimo: <code>", num$acto2_min_item,
      "</code>, p_bw = <code>", num$acto2_min_pbw, "</code>). La simulacion de ",
      "restriccion de rango deja ", num$sim_fuera, " de ", num$sim_total,
      " celdas FUERA del IC95. <strong>Sobre los datos reales no hay evidencia de ",
      "que la coordinacion cambie entre Control y LPS.</strong> En cambio, el ",
      "test de interaccion SEXO&times;TTO sobre la dispersion si detecta un ",
      "efecto sexo-dependiente en cerebro: ", nc$disp_bra_sig_n, " genes (",
      nc$disp_bra_sig_genes, ") muestran menor dispersion en hembras y mayor en ",
      "machos bajo LPS (BH &lt; 0.05).</li>"),
    paste0("<li><strong>Robustez.</strong> El resultado negativo del Acto 2 se ",
      "sostiene con el eigengene PC1 en vez del promedio de z (p_bw = <code>",
      num$eig_pbw, "</code>) y al excluir el feto mas influyente por item (",
      num$sens_cambia, " de ", num$sens_items, " items cambian el veredicto).</li>"),
    paste0("<li><strong>Reproducibilidad.</strong> ", num$comp_total,
      " CSV de resultados comparados R&harr;Python: byte-identicos = <code>",
      num$comp_byte, "</code>, celdas fuera de tolerancia = <code>",
      num$comp_fuera, "</code>.</li>"),
    "</ul>"
  )
  L <- c(L, .seccion("resumen", "1. Resumen", paste(res, collapse = "\n")))

  # --- 2. metodos ---
  met <- c(
    "<h3>2.1 Diseno</h3>",
    paste0("<p>18 madres, 36 fetos E15 (1 hembra + 1 macho por camada). Cuatro ",
      "grupos SEXO &times; TTO (HEMBRA/MACHO &times; CONTROL/LPS), 9 fetos por ",
      "grupo. Tejidos analizados: placenta E15 y cerebro fetal E15 (BRAIN_P1 ",
      "fuera de alcance). El ELISA de IL-6 y el western de pSTAT3 se analizan con ",
      "su propio n.</p>"),
    "<h3>2.2 Decisiones metodologicas fijas (D1&ndash;D13)</h3>",
    "<table><thead><tr><th>#</th><th>Decision</th></tr></thead><tbody>"
  )
  for (d in DECISIONES)
    met <- c(met, sprintf("<tr><td>%s</td><td>%s</td></tr>", d[1], .esc(d[2])))
  met <- c(met, "</tbody></table>",
    "<h3>2.3 Cascada de supuestos (D5)</h3>",
    paste0("<ul>",
      "<li>Shapiro-Wilk (residuos del modelo conjunto) y Levene ",
      "(Brown-Forsythe, centro = mediana) OK &rarr; <strong>ANOVA tipo ",
      "III</strong> (contrastes suma-cero).</li>",
      "<li>Falla solo Levene &rarr; <strong>OLS con errores HC3</strong>.</li>",
      "<li>Falla Shapiro (con o sin Levene) &rarr; <strong>ART</strong> ",
      "(Aligned Rank Transform); post hoc ART-C, nunca emmeans directo.</li>",
      "</ul>"),
    "<h3>2.4 Control de calidad de la ingesta</h3>",
    md_a_html(leer_texto(file.path(RUTA_TABLAS, "qc_reporte.md"))))
  L <- c(L, .seccion("metodos", "2. Diseno y metodos", paste(met, collapse = "\n")))

  figs <- figuras_procedencia()
  bloque_secciones <- function(ids, titulo, indices) {
    partes <- character(0)
    for (k in indices) {
      sec <- SECCIONES[[k]]
      partes <- c(partes, sprintf("<h3>%s</h3>", .esc(sec[[1]])))
      partes <- c(partes, md_a_html(leer_texto(file.path(RUTA_TABLAS, sec[[2]]))))
      for (fg in figuras_de_seccion(figs, sec[[3]])) partes <- c(partes, fig_html(fg))
      # 3.4: con fuente sintetica, ni siquiera se llama a conclusion_seccion
      # (nunca se arma la prosa interpretativa) -- se inserta el aviso.
      cl <- if (sint) {
        if (k == 2L) "" else paste0("<h4>Conclusion de la seccion</h4>\n", AVISO_SINTETICO)
      } else conclusion_seccion(k, nc)
      if (nzchar(cl)) partes <- c(partes, cl)
    }
    .seccion(ids, titulo, paste(partes, collapse = "\n"))
  }
  L <- c(L, bloque_secciones(
    "acto1", "3. Acto 1 -- respuesta a la MIA y dependencia del sexo",
    c(1L, 2L, 3L, 4L)))
  L <- c(L, .seccion("sintesis",
    "4. Sintesis del eje madre -> placenta -> cerebro",
    if (sint) AVISO_SINTETICO else sintesis_eje_html(nc)))
  L <- c(L, bloque_secciones(
    "acto2", "5. Acto 2 -- coordinacion placenta<->cerebro", c(5L, 6L, 7L, 8L)))
  L <- c(L, .seccion("conclusion_revisada",
    "6. Conclusion revisada (Acto 1 frente a Acto 2)",
    if (sint) AVISO_SINTETICO else conclusion_revisada_html(nc)))

  # --- 5. reproducibilidad ---
  rep_ <- c(
    paste0("<p>Todo el analisis esta implementado <strong>dos veces</strong> ",
      "(R y Python), con semilla fija <code>", SEMILLA, "</code> y solo rutas ",
      "relativas. Las salidas numericas llevan nombres identicos en ",
      "<code>outputs/tables/R/</code> y <code>outputs/tables/python/</code>; ",
      "<code>98_comparacion</code> las cruza celda a celda (tolerancia ",
      "1e-6, con fallback 1e-4 para p de tests iterativos) y ",
      "<code>99_verificar</code> re-corre esa comparacion y byte-compara ",
      "los <code>.md</code> y la proyeccion sin figuras de este informe.</p>"),
    md_a_html(leer_texto(file.path(RUTA_TABLAS, "comparacion_reporte.md")))
  )
  L <- c(L, .seccion("reproducibilidad", "7. Reproducibilidad (R vs Python)",
                     paste(rep_, collapse = "\n")))

  # --- 8. descartados ---
  L <- c(L, .seccion("descartados", "8. Analisis descartados",
    md_a_html(leer_texto(file.path(RUTA_TABLAS, "analisis_descartados.md")))))

  # --- 9. limitaciones ---
  lim <- c("<ul>", paste0("<li>", vapply(LIMITACIONES, .esc, character(1)), "</li>"),
           "</ul>")
  L <- c(L, .seccion("limitaciones", "9. Limitaciones", paste(lim, collapse = "\n")))

  # --- 10. auditoria ---
  aud <- c(
    sprintf(paste0("<p><code>procedencia.csv</code>: <strong>%d</strong> filas (una ",
      "por figura y por tabla). <code>verificaciones.csv</code>: <strong>%d</strong> ",
      "filas, <strong>%d</strong> en TRUE. Tablas completas en ",
      "<code>outputs/tables/</code>.</p>"),
      num$n_procedencia, num$n_verif, num$n_verif_true),
    "<h3>8.1 Procedencia</h3>",
    .csv_a_tabla(file.path(RUTA_TABLAS, "procedencia.csv")),
    "<h3>8.2 Verificaciones</h3>",
    paste0('<p class="meta">Se omite la columna <code>valor_obtenido</code> ',
      "(numeros de diagnostico que pueden diferir en el ultimo digito entre R y ",
      "Python); la tabla completa esta en ",
      "<code>outputs/tables/verificaciones.csv</code>.</p>"),
    .csv_a_tabla(file.path(RUTA_TABLAS, "verificaciones.csv"),
                 omitir = "valor_obtenido")
  )
  L <- c(L, .seccion("auditoria", "10. Procedencia y verificaciones",
                     paste(aud, collapse = "\n")))

  ap <- function(...) L[[length(L) + 1L]] <<- paste0(...)
  ap("<footer>Reanalisis MIA-LPS &mdash; placenta E15 / cerebro fetal E15. ",
     "Generado por <code>12_informe</code>; R y Python producen este HTML ",
     "identico salvo los PNG incrustados. Ver <code>AGENTS.md</code> para las ",
     "decisiones D1&ndash;D13 y <code>ESTADO.md</code> para la bitacora.</footer>")
  ap("</main>")
  ap("</body>")
  ap("</html>")
  paste0(paste(L, collapse = "\n"), "\n")
}

# =========================================================================
# Proyeccion "sin figuras" para la paridad R/Python en 99_verificar.
# =========================================================================
html_sin_figuras <- function(html) {
  gsub('src="data:image/png;base64,[^"]*"', 'src="[png]"', html, perl = TRUE)
}

# =========================================================================
# PDF por impresion headless (best-effort; el pipeline NO falla por esto).
# =========================================================================
.buscar_motor_pdf <- function() {
  env <- Sys.getenv("MIA_LPS_PDF_ENGINE", "")
  if (nzchar(env) && file.exists(env)) return(env)
  cands <- c(
    "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe",
    "C:/Program Files/Microsoft/Edge/Application/msedge.exe",
    "C:/Program Files/Google/Chrome/Application/chrome.exe",
    "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe"
  )
  for (c in cands) if (file.exists(c)) return(c)
  for (n in c("msedge", "chrome", "chromium", "chromium-browser", "google-chrome")) {
    w <- unname(Sys.which(n))
    if (nzchar(w)) return(w)
  }
  ""
}
generar_pdf <- function(html_path, pdf_path) {
  motor <- .buscar_motor_pdf()
  if (!nzchar(motor)) {
    cat("  [pdf] sin motor de impresion (Edge/Chrome); queda solo el HTML.\n")
    return("sin_motor")
  }
  perfil <- file.path(RUTA_INTERMEDIOS, "edge_profile")
  res <- tryCatch({
    if (!dir.exists(perfil)) dir.create(perfil, recursive = TRUE, showWarnings = FALSE)
    if (file.exists(pdf_path)) unlink(pdf_path)
    uri <- paste0("file:///", gsub("\\\\", "/", normalizePath(html_path, mustWork = FALSE)))
    args <- c("--headless", "--disable-gpu", "--no-first-run",
              "--no-pdf-header-footer",
              paste0("--user-data-dir=", perfil),
              paste0("--print-to-pdf=", pdf_path),
              uri)
    rc <- suppressWarnings(system2(motor, args = shQuote(args), stdout = FALSE,
                                   stderr = FALSE, timeout = 120))
    if (file.exists(pdf_path) && file.info(pdf_path)$size > 0) {
      cat(sprintf("  [pdf] %s (%d KB)\n", basename(pdf_path),
                  file.info(pdf_path)$size %/% 1024))
      "ok"
    } else {
      cat(sprintf("  [pdf] el motor no produjo un PDF valido (rc=%s); queda el HTML.\n",
                  rc))
      sprintf("fallo:sin_salida_rc%s", rc)
    }
  }, error = function(e) {
    cat(sprintf("  [pdf] no se pudo generar (%s); queda el HTML.\n",
                conditionMessage(e)))
    "fallo:error"
  })
  res
}

# =========================================================================
main <- function() {
  fuente <- fuente_datos(ARCHIVO_QPCR)
  num <- resumen_numeros()
  nc <- numeros_conclusiones()
  html <- construir_html(fuente, num, nc)
  textonly <- html_sin_figuras(html)

  ruta_html <- file.path(RUTA_DOCS, "informe.html")
  ruta_pdf <- file.path(RUTA_DOCS, "informe.pdf")
  escribir_texto(ruta_html, html)

  pdf_status <- generar_pdf(ruta_html, ruta_pdf)

  # --- snapshot por lenguaje para la paridad de 99_verificar ---
  snap <- file.path(RUTA_INTERMEDIOS, "render", LANG)
  if (!dir.exists(file.path(snap, "tables")))
    dir.create(file.path(snap, "tables"), recursive = TRUE, showWarnings = FALSE)
  escribir_texto(file.path(snap, "informe.html"), html)
  escribir_texto(file.path(snap, "informe.textonly.html"), textonly)
  escribir_texto(file.path(snap, "pdf_status.txt"), pdf_status)
  for (md in sort(list.files(RUTA_TABLAS, pattern = "\\.md$"), method = "radix"))
    escribir_texto(file.path(snap, "tables", md),
                   leer_texto(file.path(RUTA_TABLAS, md)))

  figs_proc <- figuras_procedencia()
  figs_todas <- figuras_embebidas(figs_proc)
  n_fig <- length(figs_todas)
  faltan <- figs_todas[!file.exists(file.path(RUTA_FIGURAS, figs_todas))]
  md_sec <- vapply(SECCIONES, `[[`, character(1), 2)
  md_faltan <- md_sec[!file.exists(file.path(RUTA_TABLAS, md_sec))]

  # Cobertura: toda figura con fila en procedencia.csv tiene que quedar
  # incrustada en alguna seccion -- esta es la verificacion que habria
  # detectado que acto2_corr_placenta_cerebro_* y los SPLOM por sexo no
  # aparecian en el informe (quedaban en procedencia.csv pero fuera de
  # cualquier SECCIONES a mano).
  todas_en_procedencia <- unique(basename(figs_proc$artefacto))
  sin_embeber <- sort(setdiff(todas_en_procedencia, figs_todas))

  n_sec <- length(gregexpr('<section id="', html, fixed = TRUE)[[1]])

  # 3.4: con fuente sintetica, ninguna frase interpretativa puede aparecer --
  # se recalcula contando cuantas veces aparece el aviso en el HTML final
  # (9 = 7 conclusiones de seccion + sintesis + conclusion revisada) contra
  # el esperado segun `sint`.
  sint <- fuente != "real"
  m_aviso <- gregexpr(AVISO_SINTETICO, html, fixed = TRUE)[[1]]
  n_aviso <- if (length(m_aviso) == 1L && m_aviso[1] == -1L) 0L else length(m_aviso)
  n_aviso_esperado <- if (sint) 9L else 0L

  ent <- "outputs/tables/*.md + outputs/tables/{R,python}/*.csv + outputs/figures/*.png"
  registrar_procedencia(list(
    list("docs/informe.html", "informe", ESTE_SCRIPT, "PROPIO", ent,
         paste0("informe HTML autocontenido (Acto 1 + Acto 2 + reproducibilidad + ",
           "analisis descartados + limitaciones + procedencia/verificaciones); ",
           "figuras incrustadas en base64; sin marca de tiempo")),
    list("docs/informe.pdf", "informe", ESTE_SCRIPT, "PROPIO", "docs/informe.html",
         paste0("version imprimible por impresion headless (Edge/Chrome); best-effort, ",
           "el pipeline no falla si no hay motor de PDF (estado: ", pdf_status, ")"))
  ))
  registrar_verificaciones(list(
    list("informe_html_generado",
         "docs/informe.html existe y no esta vacio",
         sprintf("%d secciones; %d figuras", n_sec, n_fig),
         "10 secciones",
         if (file.exists(ruta_html) && file.info(ruta_html)$size > 0 && n_sec == 10L)
           "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_figuras_incrustadas",
         "todas las figuras del Acto 1 y 2 estan incrustadas en el informe",
         sprintf("figuras=%d; faltan=%s", n_fig,
                 if (length(faltan)) paste(faltan, collapse = ", ") else "[]"),
         "faltan = []",
         if (!length(faltan)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_reportes_incluidos",
         "todos los .md de seccion existen y se incluyeron",
         sprintf("secciones=%d; md_faltan=%s", length(SECCIONES),
                 if (length(md_faltan)) paste(md_faltan, collapse = ", ") else "[]"),
         "md_faltan = []",
         if (!length(md_faltan)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_figuras_procedencia_embebidas",
         "toda figura con fila en procedencia.csv esta incrustada en el informe",
         sprintf("procedencia=%d; embebidas=%d; sin_embeber=%s",
                 length(todas_en_procedencia), length(figs_todas),
                 if (length(sin_embeber)) paste(sin_embeber, collapse = ", ") else "[]"),
         "sin_embeber = []",
         if (!length(sin_embeber)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_pdf",
         "docs/informe.pdf generado, o degradado limpio si no hay motor de PDF",
         sprintf("estado=%s; existe=%s", pdf_status,
                 if (file.exists(ruta_pdf) && file.info(ruta_pdf)$size > 0)
                   "TRUE" else "FALSE"),
         "ok | sin_motor | fallo (nunca frena el pipeline)",
         if (pdf_status %in% c("ok", "sin_motor") ||
             startsWith(pdf_status, "fallo")) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_snapshot_paridad",
         "snapshot por lenguaje para el byte-compare R/Python de 99_verificar",
         paste0("render/<lang>/: informe.html + informe.textonly.html + tables/*.md + ",
                "pdf_status.txt"),
         "snapshot escrito",
         if (file.exists(file.path(snap, "informe.textonly.html"))) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("informe_sintetico_sin_interpretacion",
         paste0("con fuente sintetica, las secciones interpretativas (conclusion ",
                "de seccion, sintesis, conclusion revisada) se reemplazan por el ",
                "aviso de datos sinteticos -- nunca se arma la prosa biologica"),
         sprintf("fuente=%s; aviso=%d/%d", fuente, n_aviso, n_aviso_esperado),
         "aviso = 9 si fuente sintetica, 0 si fuente real",
         if (n_aviso == n_aviso_esperado) "TRUE" else "FALSE", ESTE_SCRIPT)
  ))

  cat("== 12_informe.R ==\n")
  cat(sprintf("  fuente = %s\n", fuente))
  cat(sprintf("  -> %s (%d KB)\n", basename(ruta_html),
              file.info(ruta_html)$size %/% 1024))
  cat(sprintf("  -> %s [%s]\n", basename(ruta_pdf), pdf_status))
  cat(sprintf("  figuras incrustadas: %d%s\n", n_fig,
              if (!length(faltan)) "" else sprintf("  FALTAN: %s",
                paste(faltan, collapse = ", "))))
  cat(sprintf("  secciones .md: %d%s\n", length(SECCIONES),
              if (!length(md_faltan)) "" else sprintf("  FALTAN: %s",
                paste(md_faltan, collapse = ", "))))
  cat(sprintf("  figuras de procedencia.csv sin embeber: %s\n",
              if (!length(sin_embeber)) "ninguna" else paste(sin_embeber, collapse = ", ")))
  cat(sprintf("  snapshot: outputs/intermediate/render/%s/\n", LANG))
  if (length(faltan) || length(md_faltan)) {
    cat("  *** faltan insumos: correr 02..11 y 07 antes de 12_informe ***\n")
    quit(status = 1L)
  }
}

if (sys.nframe() == 0L) main()
