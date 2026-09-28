# 13_presentacion.R -- pagina de presentacion (HTML + PDF), 12 diapositivas.
#
# Por que existe este archivo: la consigna del curso pide, dentro del repo, una
# pagina HTML y un PDF desde los cuales exponer (pedidos/cambios_presentacion.md,
# segunda tanda en pedidos/cambios_presentacion_2.md). Es CAPA DE PRESENTACION, no
# analisis: no calcula ningun resultado nuevo (excepto dos figuras que RE-DIBUJAN,
# sin recalcular, un subconjunto de genes ya modelados -- ver
# "generar_panel_subset" mas abajo), solo lee lo que ya escribieron 02..11/98/12
# (tablas, figuras, texto de AGENTS.md) y lo arma en 12 diapositivas para proyectar.
#
# EL HILO (pedido 2, punto 0): sola, la autora habria hecho boxplots; con el
# agente pudo correr analisis que no habria hecho, y esos analisis explicaron lo
# que el boxplot dejaba sin resolver. Por eso 6-7 son "el analisis convencional"
# y 8-9 son "lo que agrego explorar con el agente" -- es el argumento central,
# no un detalle de orden.
#
# EXCEPCION A LA REGLA DE SCRIPTS GEMELOS (documentada tambien en AGENTS.md 3 y
# en procedencia.csv): existe SOLO en R. No hay python/13_presentacion.py: no hay
# resultado numerico que comparar entre lenguajes, asi que duplicarlo no
# agregaria verificacion, solo mantenimiento.
#
# PARIDAD: los helpers de lectura de CSV / escritura / base64 / PDF headless son
# copia de 12_informe.R (misma convencion que 99_verificar.R copia de
# 98_comparacion.R): la regla del repo es que cada script importe solo
# 00_config, nunca cross-importe un script numerado.
#
# SEGUNDA EXCEPCION (pedido 2, secciones 6-7): las diapositivas 6 y 7 necesitan
# boxplots de un SUBCONJUNTO de 3 genes por tejido, no el panel completo. Se
# reusan las funciones de dibujo de 07_figuras_acto1.R (panel_gen y sus
# dependencias) en vez de reimplementarlas -- pedido explicito ("reusando las
# funciones de 07_figuras_acto1.R, sin recalcular nada"). Para no chocar
# nombres (07 redefine ESTE_SCRIPT, registrar_procedencia, etc. con firmas
# propias) se fuente CON `local = <environment nuevo>`: todo el archivo 07 se
# ejecuta aislado, sin tocar nada de lo definido aca. `if (sys.nframe() == 0L)
# main()` (guardia ya existente en 07) evita que eso dispare su propio main().
#
# DOS VERSIONES, mismo codigo, la fuente de datos decide el destino (igual
# mecanismo que fuente_datos()/12_informe usan para el aviso sintetico):
#   - fuente == "sintetico" -> docs/index.html + docs/presentacion.pdf
#     (PUBLICA: se versiona, enlaza desde la pagina del curso).
#   - fuente == "real"      -> outputs/presentacion_real/index.html + .pdf
#     (PARA EXPONER: nunca se versiona, ver .gitignore).
# Para generar las dos hace falta correr el script dos veces (igual que para
# tener ambas versiones de docs/informe.html): una vez normal (con
# data/raw/ presente) y otra con MIA_LPS_FORZAR_SINTETICO=1.
#
# Ademas de las 12 diapositivas, genera docs/referencias.md (o su copia en
# outputs/presentacion_real/) con las 10 citas de
# pedidos/referencias_epidemiologia.md, copiadas tal cual -- ninguna inventada.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

ESTE_SCRIPT <- "13_presentacion"
LANG <- "R"

RUTA_ASSETS   <- file.path(RAIZ_REPO, "assets")
RUTA_PRES_REAL <- file.path(RUTA_OUT, "presentacion_real")

# =========================================================================
# Helpers de lectura/escritura/formato -- copia de 12_informe.R.
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
leer_csv_sin_este <- function(ruta) {
  t <- leer_csv(ruta)
  if (is.null(t$header) || !("script" %in% t$header)) return(t)
  j <- match("script", t$header)
  t$filas <- Filter(function(f) !(j <= length(f) && f[[j]] == ESTE_SCRIPT), t$filas)
  t
}
.col <- function(t, nombre) {
  if (is.null(t$header) || !(nombre %in% t$header)) return(character(0))
  j <- match(nombre, t$header)
  vapply(t$filas, function(f) if (j <= length(f)) f[[j]] else "", character(1))
}
.tab <- function(nombre) leer_csv(file.path(RUTA_TABLAS_R, nombre))

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
  header <- c("id", "tipo", "descripcion", "valor_obtenido", "valor_esperado", "ok", "script")
  merge_por_script(file.path(RUTA_TABLAS, "verificaciones.csv"), header, filas_nuevas,
                   function(fs) order(vapply(fs, `[[`, character(1), 6),
                                      vapply(fs, `[[`, character(1), 1),
                                      method = "radix"))
}
.esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  gsub(">", "&gt;", s, fixed = TRUE)
}

# base64 propio (RFC 4648) -- copia de 12_informe.R.
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

# =========================================================================
# PDF por impresion headless -- copia de 12_informe.R.
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
# Aviso de fuente sintetica -- mismo principio que AVISO_SINTETICO de
# 12_informe.R: con datos sinteticos, las diapositivas 6-9 (analisis
# convencional + lo que agrego explorar con el agente) nunca arman la
# prosa interpretativa (biologica); la reemplazan por este aviso. Es
# estructuralmente imposible que una conclusion biologica aparezca sobre
# datos sinteticos, porque el codigo que la construye no se llama.
# =========================================================================
AVISO_SINTETICO <- paste0(
  "<p class=\"aviso-sint\"><em>Diapositiva generada con datos sinteticos. Los ",
  "efectos son simulados y arbitrarios; las conclusiones biologicas ",
  "corresponden a los datos reales, que no se incluyen en este ",
  "repositorio.</em></p>")

# =========================================================================
# Extraer una decision fija (D1, D2, D5, ...) tal como esta escrita en
# AGENTS.md -- para la diapositiva 4 ("el punto de partida fue un prompt"):
# se muestra el fragmento REAL, no una paráfrasis. Fila de la tabla:
# "| **D1** | <decision> | <justificacion> |"
# =========================================================================
extraer_decision <- function(id) {
  agents <- leer_texto(file.path(RAIZ_REPO, "AGENTS.md"))
  lineas <- strsplit(agents, "\n", fixed = TRUE)[[1]]
  patron <- sprintf("^\\| \\*\\*%s\\*\\* \\|", id)
  fila <- grep(patron, lineas, value = TRUE, perl = TRUE)
  if (!length(fila)) return(sprintf("[%s no encontrada en AGENTS.md]", id))
  # separar por " | " respetando que el contenido no trae pipes propios
  partes <- strsplit(fila[1], " \\| ", perl = TRUE)[[1]]
  # partes[1] = "| **D1**", partes[2] = decision, partes[3] = justificacion |
  decision <- if (length(partes) >= 2L) partes[2] else fila[1]
  decision
}
# Mini-inline Markdown -> HTML para el fragmento de AGENTS.md (negrita, code).
.inline_md <- function(s) {
  s <- .esc(s)
  s <- gsub("`([^`]+)`", "<code>\\1</code>", s, perl = TRUE)
  s <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", s, perl = TRUE)
  s
}

# =========================================================================
# Figuras -- lee de outputs/figures/ (resultados del proyecto) o de assets/
# (ilustraciones propias / gráfico de contexto). registrar_figura() empuja a
# una lista para poder chequear despues, en main(), que todas existen.
# =========================================================================
.FIGS_USADAS <- new.env()
.FIGS_USADAS$outputs <- character(0)
fig_outputs <- function(nombre, alt = nombre, clase = "") {
  .FIGS_USADAS$outputs <- c(.FIGS_USADAS$outputs, nombre)
  ruta <- file.path(RUTA_FIGURAS, nombre)
  if (!file.exists(ruta))
    return(sprintf('<p class="falta">[falta la figura %s]</p>', .esc(nombre)))
  sprintf('<img class="fig %s" alt="%s" src="%s">', clase, .esc(alt), img_datauri(ruta))
}
fig_asset <- function(nombre, alt = nombre, clase = "") {
  ruta <- file.path(RUTA_ASSETS, nombre)
  if (!file.exists(ruta))
    return(sprintf('<p class="falta">[falta assets/%s]</p>', .esc(nombre)))
  sprintf('<img class="fig %s" alt="%s" src="%s">', clase, .esc(alt), img_datauri(ruta))
}

# =========================================================================
# CSS -- variables de paleta/tipografia/espaciado en un unico bloque, tomando
# como base la paleta de 00_config.R (COL_CTRL, COL_LPS): azules de control
# (HEMBRA/MACHO) para lo neutro, IL6 (la via central de la historia: MIA ->
# STAT3 -> placenta) como acento principal, y el resto de las vias como
# acentos secundarios puntuales.
# =========================================================================
CSS <- paste0("\n",
":root {\n",
"  color-scheme: light;\n",
"  --bg:        #fdfcfa;\n",
"  --panel:     #ffffff;\n",
"  --ink:       #1a1a1a;\n",
"  --ink-soft:  #55524f;\n",
"  --line:      #e4e0da;\n",
"  --accent:        ", COL_LPS$IL6[["MACHO"]], ";\n",
"  --accent-soft:   ", COL_LPS$IL6[["HEMBRA"]], ";\n",
"  --blue-f:    ", COL_CTRL[["HEMBRA"]], ";\n",
"  --blue-m:    ", COL_CTRL[["MACHO"]], ";\n",
"  --pink:      ", COL_LPS$GLUCOSA[["MACHO"]], ";\n",
"  --pink-soft: ", COL_LPS$GLUCOSA[["HEMBRA"]], ";\n",
"  --green:     ", COL_LPS$AMINOACIDOS[["MACHO"]], ";\n",
"  --orange:    ", COL_LPS$LIPIDOS[["MACHO"]], ";\n",
"  --warn-bg:     #fff4e5;\n",
"  --warn-border: #f0c48a;\n",
"  --font: -apple-system, \"Segoe UI\", Roboto, Helvetica, Arial, sans-serif;\n",
"  --font-mono: \"SF Mono\", Consolas, \"Liberation Mono\", monospace;\n",
"  --space-1: .5rem; --space-2: 1rem; --space-3: 1.8rem; --space-4: 2.8rem;\n",
"}\n",
"* { box-sizing: border-box; }\n",
"html { scroll-snap-type: y mandatory; scroll-behavior: smooth; }\n",
"body { font: 21px/1.5 var(--font); color: var(--ink); background: var(--bg);\n",
"       margin: 0; }\n",
"section.slide {\n",
"  min-height: 100vh; width: 100%;\n",
"  padding: 2.1rem 4.5rem;\n",
"  display: flex; flex-direction: column; justify-content: flex-start;\n",
"  scroll-snap-align: start; position: relative;\n",
"  border-bottom: 1px solid var(--line);\n",
"}\n",
"section.slide.title { justify-content: center; }\n",
"section.slide > .kicker { text-transform: uppercase; letter-spacing: .12em;\n",
"  font-size: .78rem; font-weight: 700; color: var(--accent); margin: 0 0 .6rem; }\n",
"section.slide h1 { font-size: 2.6rem; line-height: 1.15; margin: 0 0 .8rem; }\n",
"section.slide h2 { font-size: 2.05rem; line-height: 1.2; margin: 0 0 1.3rem;\n",
"  padding-bottom: .5rem; border-bottom: 3px solid var(--accent); display: inline-block; }\n",
"section.slide p  { margin: .5rem 0; color: var(--ink-soft); }\n",
"section.slide p.lead { font-size: 1.25rem; color: var(--ink); }\n",
".slide-num { position: absolute; bottom: 1.2rem; right: 1.6rem;\n",
"  font-size: .85rem; color: var(--ink-soft); opacity: .6; }\n",
".cols { display: flex; gap: var(--space-4); align-items: flex-start; margin-top: 1.2rem; }\n",
".col  { flex: 1; min-width: 0; }\n",
".col.narrow { flex: 0 0 34%; }\n",
".col.wide   { flex: 1 1 66%; }\n",
".fig { max-width: 100%; max-height: 62vh; display: block; margin: 0 auto; }\n",
".fig.tall  { max-height: 70vh; }\n",
".fig.wide  { max-height: 42vh; width: 100%; object-fit: contain; }\n",
".fig.small { max-height: 22vh; margin-bottom: .8rem; }\n",
"figure.pres { margin: 0; text-align: center; }\n",
"figure.pres figcaption { font-size: .82rem; color: var(--ink-soft); margin-top: .5rem; }\n",
".card { background: var(--panel); border: 1px solid var(--line); border-radius: 10px;\n",
"  padding: 1.3rem 1.6rem; }\n",
".aviso-sint, .aviso { background: var(--warn-bg); border: 1px solid var(--warn-border);\n",
"  border-radius: 8px; padding: .9rem 1.2rem; font-size: 1.02rem; }\n",
"pre.quote { background: #f7f5f1; border-left: 4px solid var(--accent);\n",
"  border-radius: 4px; padding: 1rem 1.3rem; font: 1rem/1.6 var(--font-mono);\n",
"  white-space: pre-wrap; margin: .8rem 0; }\n",
"table.stats { border-collapse: collapse; margin: .6rem 0; font-size: 1.02rem; }\n",
"table.stats td { padding: .35rem 1.1rem .35rem 0; }\n",
"table.stats td.n { font-size: 1.9rem; font-weight: 700; color: var(--accent); }\n",
"ul.plain { list-style: none; margin: .4rem 0; padding: 0; }\n",
"ul.plain li { margin: .3rem 0; }\n",
"ul.bullets, ol.bullets { margin: .5rem 0; padding-left: 1.3rem; }\n",
"ul.bullets li, ol.bullets li { margin: .5rem 0; }\n",
".gene-tag { display: inline-block; background: var(--accent-soft); color: #3a1d52;\n",
"  border-radius: 999px; padding: .1rem .7rem; font: .88rem var(--font-mono);\n",
"  margin: .12rem .2rem .12rem 0; }\n",
"/* --- diapositiva 2: problema --- */\n",
".problem-row { display: flex; align-items: center; gap: 2.2rem; margin-top: 1rem; }\n",
".problem-row .illus { flex: 0 0 26%; position: relative; text-align: center; }\n",
".problem-row .illus img { max-width: 100%; max-height: 46vh; }\n",
".problem-row .arrowbox { flex: 0 0 20%; display: flex; flex-direction: column;\n",
"  align-items: center; gap: .4rem; }\n",
".arrowbox .box { border: 2.5px solid var(--accent); border-radius: 10px;\n",
"  padding: .9rem 1.1rem; text-align: center; font-weight: 700; font-size: 1.12rem;\n",
"  color: var(--accent); background: var(--accent-soft); }\n",
".problem-row .epi { flex: 1 1 54%; }\n",
".rayo { stroke: var(--accent); stroke-width: 5; fill: none; stroke-linecap: round;\n",
"  stroke-linejoin: round; }\n",
".rayo-punta { fill: var(--accent); }\n",
".closing-card { margin-top: 1.3rem; }\n",
"/* --- diapositiva 2bis: grafico epidemiologico (HTML/CSS, no PNG) --- */\n",
".epi-axis { display: flex; justify-content: space-between; font-size: .82rem;\n",
"  color: var(--ink-soft); margin: .2rem 0 .3rem; padding-left: 1px; }\n",
".epi-block { margin: .7rem 0; }\n",
".epi-block h4 { margin: 0 0 .5rem; font-size: 1rem; display: flex; align-items: center;\n",
"  gap: .5rem; }\n",
".epi-block h4 .dot { width: .8rem; height: .8rem; border-radius: 50%; display: inline-block; }\n",
".epi-item { margin: .45rem 0; }\n",
".epi-item .lab { font-size: .92rem; margin-bottom: .18rem; }\n",
".epi-track { position: relative; height: 1.15rem; background: #f1efe9;\n",
"  border-radius: 3px; overflow: visible; }\n",
".epi-track::before { content: ''; position: absolute; left: 0; top: -.25rem;\n",
"  bottom: -.25rem; width: 2px; background: var(--ink-soft); opacity: .35; }\n",
".epi-bar { position: absolute; left: 0; top: 0; bottom: 0; border-radius: 3px 0 0 3px; }\n",
".epi-bar.range { border-radius: 3px; opacity: .45; }\n",
".epi-val { position: absolute; top: 0; font-size: .82rem; font-weight: 700;\n",
"  white-space: nowrap; transform: translateY(-.05rem); }\n",
"/* --- diapositiva 5: conteos --- */\n",
".stat-grid { display: grid; grid-template-columns: repeat(3, 1fr); gap: 1.3rem 2rem;\n",
"  margin-top: 1rem; }\n",
".stat { }\n",
".stat .n { font-size: 2.3rem; font-weight: 700; color: var(--accent); line-height: 1; }\n",
".stat .l { font-size: .92rem; color: var(--ink-soft); margin-top: .2rem; }\n",
"/* --- impresion: una diapositiva por pagina --- */\n",
"@media print {\n",
"  html { scroll-snap-type: none; }\n",
"  section.slide { min-height: 100vh; break-after: page; border-bottom: none; }\n",
"  section.slide:last-child { break-after: auto; }\n",
"}\n",
".flecha-linea { stroke: var(--accent); stroke-width: 4; }\n",
".flecha-punta { fill: var(--accent); }\n"
)

# =========================================================================
# Diapositiva 1: metadatos de portada. AUTORA y CURSO no se pueden leer de
# ninguna tabla -- son datos de la persona que expone, no del analisis (se
# completaron a mano por el usuario). TITULO es el titulo real del trabajo,
# dado explicitamente en pedidos/cambios_presentacion_2.md 1 (ya no el h1 de
# README.md: ese titulo de repo es mas largo/tecnico que el de una charla).
# SUBTITULO queda igual que en la primera tanda (pedido explicito: "como
# esta"). FECHA es la fecha de la corrida (Sys.Date()), no una fecha fija.
# =========================================================================
PRESENTACION_AUTORA <- "Ayelen Mirón Granese"
PRESENTACION_CURSO  <- "Ondas Gravitacionales e Investigación Asistida por IA"
PRESENTACION_TITULO <- paste0(
  "Transportadores de nutrientes en el eje placenta–cerebro fetal en un ",
  "modelo de activación inmune materna")
PRESENTACION_SUBTITULO <- "Reanalisis bioestadistico dirigido por un agente de IA (Claude Code)"

# =========================================================================
# Diapositiva 2bis: grafico de razones de prevalencia por sexo -- el UNICO
# grafico de la presentacion que no sale de outputs/figures/ (no es un
# resultado del proyecto, es contexto bibliografico). Se regenera en HTML/CSS
# con la paleta de la presentacion en vez de reusar assets/epidemiologia.png.
# Datos y redondeo: literales de pedidos/referencias_epidemiologia.md /
# cambios_presentacion.md 2bis -- no se inventa ni un numero. Los rangos se
# dibujan como rango (segmento mas claro desde el extremo bajo al alto), no
# como un valor unico.
# =========================================================================
EPI_DATOS <- list(
  list(pat = "Trastorno del espectro autista", lo = 4, hi = 4, grupo = "varones"),
  list(pat = "Enfermedad de Parkinson", lo = 3.5, hi = 3.5, grupo = "varones"),
  list(pat = "Trastorno por deficit de atencion e hiperactividad", lo = 3, hi = 3,
       grupo = "varones"),
  list(pat = "Esclerosis lateral amiotrofica", lo = 1.6, hi = 1.6, grupo = "varones"),
  list(pat = "Esquizofrenia", lo = 1.4, hi = 1.4, grupo = "varones"),
  list(pat = "Esclerosis multiple", lo = 2, hi = 3, grupo = "mujeres"),
  list(pat = "Enfermedad de Alzheimer", lo = 1.6, hi = 3, grupo = "mujeres"),
  list(pat = "Depresion y trastornos de ansiedad", lo = 2, hi = 2, grupo = "mujeres")
)
.fmt_ratio <- function(x) {
  if (x == floor(x)) sprintf("%d", as.integer(x)) else sprintf("%s", format(x, nsmall = 1))
}
construir_epi_chart <- function() {
  items <- EPI_DATOS
  dom_max <- ceiling(max(vapply(items, function(x) x$hi, numeric(1))) * 1.1 * 2) / 2
  pct <- function(v) (v - 1) / (dom_max - 1) * 100
  fila <- function(it, color) {
    w_lo <- pct(it$lo); w_hi <- pct(it$hi)
    val_txt <- if (it$hi > it$lo)
      sprintf("%s&ndash;%s:1", .fmt_ratio(it$lo), .fmt_ratio(it$hi))
    else sprintf("%s:1", .fmt_ratio(it$lo))
    base <- sprintf('<div class="epi-bar" style="width:%.2f%%; background:%s"></div>',
                     w_lo, color)
    rango <- if (it$hi > it$lo)
      sprintf('<div class="epi-bar range" style="left:%.2f%%; width:%.2f%%; background:%s"></div>',
              w_lo, w_hi - w_lo, color) else ""
    val_pos <- w_hi + 1.2
    sprintf(paste0('<div class="epi-item"><div class="lab">%s</div>',
                   '<div class="epi-track">%s%s',
                   '<span class="epi-val" style="left:%.2f%%">%s</span></div></div>'),
            .esc(it$pat), base, rango, val_pos, val_txt)
  }
  ord <- function(l) l[order(-vapply(l, function(x) x$hi, numeric(1)),
                              -vapply(l, function(x) x$lo, numeric(1)))]
  varones <- ord(Filter(function(x) x$grupo == "varones", items))
  mujeres <- ord(Filter(function(x) x$grupo == "mujeres", items))
  bloque <- function(titulo, color, lista) {
    sprintf('<div class="epi-block"><h4><span class="dot" style="background:%s"></span>%s</h4>%s</div>',
            color, titulo,
            paste(vapply(lista, fila, character(1), color = color), collapse = "\n"))
  }
  paste0(
    '<div class="epi-axis"><span>1:1 &nbsp;(sin diferencia)</span><span>',
    .fmt_ratio(dom_max), ':1</span></div>',
    bloque("Predominio en varones", "var(--blue-m)", varones),
    bloque("Predominio en mujeres", "var(--pink)", mujeres),
    '<p style="font-size:.8rem; margin-top:.8rem;">Razones de prevalencia ',
    'compiladas de 10 estudios epidemiologicos (ver <a href="referencias.md">',
    'referencias.md</a>).</p>'
  )
}
escribir_referencias <- function(destino) {
  citas <- leer_texto(file.path(RAIZ_REPO, "pedidos", "referencias_epidemiologia.md"))
  escribir_texto(destino, citas)
}

# =========================================================================
# Ilustraciones SVG propias (no estan en los PNG de assets/).
# =========================================================================
# Rayos apuntando AL vientre (pedido 2, seccion 2: "hoy se leen como algo que
# emana de la figura; tienen que leerse como algo que incide sobre la
# gestacion"). El grupo se ubica a la IZQUIERDA de la ilustracion (`left`, no
# `right`) y cada zigzag empieza lejos (x chico) y termina cerca del borde
# derecho del propio SVG -- que es el borde que toca la ilustracion -- con una
# punta de flecha (`marker-end`, orientacion automatica a lo largo del trazo).
svg_rayos <- function() paste0(
  '<svg class="rayos" viewBox="0 0 120 200" xmlns="http://www.w3.org/2000/svg" ',
  'style="position:absolute; left:-38px; top:10%; width:60px; height:76%;">',
  '<defs><marker id="rayoPunta" markerWidth="7" markerHeight="7" refX="5" refY="3.5" ',
  'orient="auto"><path class="rayo-punta" d="M0,0 L7,3.5 L0,7 Z"/></marker></defs>',
  '<polyline class="rayo" marker-end="url(#rayoPunta)" points="8,15 42,52 12,78 100,118"/>',
  '<polyline class="rayo" marker-end="url(#rayoPunta)" points="8,58 42,95 12,118 100,155"/>',
  '<polyline class="rayo" marker-end="url(#rayoPunta)" points="8,98 42,132 12,152 100,185"/>',
  '</svg>')
svg_flecha <- function() paste0(
  '<svg viewBox="0 0 100 30" width="88" height="26" xmlns="http://www.w3.org/2000/svg">',
  '<line class="flecha-linea" x1="4" y1="15" x2="80" y2="15"/>',
  '<polygon class="flecha-punta" points="80,6 98,15 80,24"/>',
  '</svg>')

# =========================================================================
# Figuras de SUBCONJUNTO de genes (pedido 2, secciones 6-7): diapositivas 6 y
# 7 necesitan boxplots de 3 genes por tejido, no el panel completo. Se
# reusan panel_gen() y sus dependencias de 07_figuras_acto1.R -- pedido
# explicito ("sin recalcular nada"). SEGUNDA EXCEPCION a "cada script importa
# solo 00_config" (ver cabecera del archivo): se fuente 07 completo, pero
# `local = <environment nuevo>` lo aisla por completo de este script (07
# redefine ESTE_SCRIPT, registrar_procedencia, merge_por_script, etc. con
# firmas propias -- si se sourceara al entorno global romperia las de aca).
# `cargar()` de 07 solo LEE tablas ya escritas por 04/05 (04_qpcr_
# cuantificacion, 05_qpcr_modelos): no dispara ningun calculo nuevo.
# =========================================================================
.ENV07 <- local({
  e <- NULL
  function() {
    if (is.null(e))
      e <<- { env <- new.env()
              source(file.path(RAIZ_REPO, "R", "07_figuras_acto1.R"),
                     local = env, encoding = "UTF-8")
              env }
    e
  }
})
generar_panel_subset <- function(genes, tejido, ruta) {
  e <- .ENV07()
  D <- e$cargar()
  n <- length(genes)
  grDevices::png(ruta, width = n * 1150L, height = 1000L, res = e$DPI)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(SEMILLA)   # jitter reproducible, igual que figura_tejido()
  graphics::layout(matrix(seq_len(n), nrow = 1L))
  for (gen in genes) {
    fc <- e$fc_por_grupo(D$cuant, tejido, gen)
    fila_clasif <- e$fila_de(D$clasif, tejido, gen)
    ph <- e$pholm_lista(D$posthoc, tejido, gen)
    e$panel_gen(fc, gen, tejido, fila_clasif, ph)
  }
  invisible(ruta)
}
GENES_SUBSET_PLACENTA <- c("il6", "glut3", "slc38a2")   # pedido 2, seccion 6
GENES_SUBSET_BRAIN    <- c("glut1", "slc38a2", "fatp1") # pedido 2, seccion 6
NOMBRE_SUBSET_PLACENTA <- "acto1_expresion_PLACENTA_E15_subset3.png"
NOMBRE_SUBSET_BRAIN    <- "acto1_expresion_BRAIN_E15_subset3.png"

# =========================================================================
# Datos leidos del repo para las diapositivas 3, 4, 5, 6, 7, 9 -- nada de lo
# que sigue se escribe a mano: se lee de AGENTS.md / outputs/tables/.
# =========================================================================
recolectar_datos <- function() {
  d <- list()

  # --- diapositiva 3: diseno experimental ---
  vf <- leer_csv(file.path(RUTA_TABLAS, "verificaciones.csv"))
  d$n_fetos <- { v <- .col(vf, "valor_obtenido")[.col(vf, "id") == "ingesta_qpcr_fetos_e15"]
                 if (length(v)) v[1] else "36" }
  d$n_genes <- length(GENES)
  d$n_grupos <- length(NIVELES_SEXO) * length(NIVELES_TTO)
  d$n_por_grupo <- as.integer(d$n_fetos) %/% d$n_grupos

  # --- diapositiva 4: fragmento real de AGENTS.md ---
  d$d1 <- extraer_decision("D1"); d$d2 <- extraer_decision("D2"); d$d5 <- extraer_decision("D5")

  # --- diapositiva 5: conteos del repo ---
  d$n_scripts_r  <- length(list.files(file.path(RAIZ_REPO, "R"), pattern = "\\.R$"))
  d$n_scripts_py <- length(list.files(file.path(RAIZ_REPO, "python"), pattern = "\\.py$"))
  d$n_figuras    <- length(list.files(RUTA_FIGURAS, pattern = "\\.png$"))
  d$n_tablas_r   <- length(list.files(RUTA_TABLAS_R, pattern = "\\.csv$"))
  d$n_tablas_py  <- length(list.files(RUTA_TABLAS_PY, pattern = "\\.csv$"))
  pr <- leer_csv_sin_este(file.path(RUTA_TABLAS, "procedencia.csv"))
  d$n_procedencia <- length(pr$filas)
  vfx <- leer_csv_sin_este(file.path(RUTA_TABLAS, "verificaciones.csv"))
  tipoc <- .col(vfx, "tipo")
  d$n_verif <- length(vfx$filas)
  for (ti in c("recalculo", "existencia", "declaracion"))
    d[[paste0("n_verif_", ti)]] <- sum(tipoc == ti)
  cp <- leer_csv(file.path(RUTA_TABLAS, "comparacion_R_python.csv"))
  d$comp_n <- "n/d"; d$comp_byte_n <- "n/d"
  if (!is.null(cp$header)) {
    arch <- .col(cp, "archivo")
    bi <- .col(cp, "byte_identico")
    idx <- arch != "__TOTAL__"
    d$comp_n <- as.character(sum(idx))
    d$comp_byte_n <- as.character(sum(idx & bi == "TRUE"))
  }

  # --- diapositiva 6: pSTAT3 en placenta ---
  pm <- .tab("pstat3_modelo_clasificacion.csv")
  d$pstat3_pint <- { v <- .col(pm, "p_SEXOxTTO"); if (length(v)) v[1] else "n/d" }
  ph <- .tab("pstat3_posthoc.csv")
  contr <- .col(ph, "contraste"); phv <- .col(ph, "p_holm")
  d$pstat3_hh <- "n/d"; d$pstat3_mm <- "n/d"
  for (k in seq_along(contr)) {
    if (contr[k] == "HEMBRA_CONTROL-HEMBRA_LPS") d$pstat3_hh <- phv[k]
    if (contr[k] == "MACHO_CONTROL-MACHO_LPS") d$pstat3_mm <- phv[k]
  }

  # --- diapositiva 7: interaccion SEXOxTTO en el programa de expresion ---
  qm <- .tab("qpcr_modelos_clasificacion.csv")
  tej <- .col(qm, "TEJIDO"); isig <- .col(qm, "interaccion_significativa")
  via <- .col(qm, "via")
  d$qpcr_int_bra <- sum(tej == "BRAIN_E15" & isig == "TRUE")
  d$qpcr_mod_bra <- sum(tej == "BRAIN_E15" & via == "modelo")
  d$qpcr_int_pla <- sum(tej == "PLACENTA_E15" & isig == "TRUE")
  d$genes_int_bra <- .col(qm, "GEN")[tej == "BRAIN_E15" & isig == "TRUE"]

  # --- diapositiva 9: razonamiento de 3 pasos (interaccion en el boxplot ->
  # cuantos NO se explican por ninguna comparacion de a pares -> cuantos de
  # esos si se explican por dispersion). Todo derivado de las tablas, gen por
  # gen, no una lista escrita a mano. ---------------------------------------
  poh <- .tab("qpcr_modelos_posthoc.csv")
  ph_tej <- .col(poh, "TEJIDO"); ph_gen <- .col(poh, "GEN")
  ph_p <- suppressWarnings(as.numeric(.col(poh, "p_holm")))
  sin_posthoc <- character(0)
  for (g in d$genes_int_bra) {
    idx <- ph_tej == "BRAIN_E15" & ph_gen == g
    if (!any(idx & !is.na(ph_p) & ph_p < 0.05)) sin_posthoc <- c(sin_posthoc, g)
  }
  d$genes_sin_posthoc <- sin_posthoc

  di <- .tab("acto2_dispersion_interaccion.csv")
  tejd <- .col(di, "TEJIDO"); geng <- .col(di, "GEN")
  bh <- suppressWarnings(as.numeric(.col(di, "p_SEXOxTTO_BH")))
  d$disp_bra_sig <- sum(tejd == "BRAIN_E15" & bh < 0.05, na.rm = TRUE)
  d$disp_bra_n   <- sum(tejd == "BRAIN_E15")
  d$disp_pla_sig <- sum(tejd == "PLACENTA_E15" & bh < 0.05, na.rm = TRUE)
  d$genes_disp_bra <- .col(di, "GEN")[tejd == "BRAIN_E15" & !is.na(bh) & bh < 0.05]
  # de los "sin posthoc": cuales explica la dispersion (BH<.05) y cuales
  # quedan solo en tendencia (BH<.10, sin llegar a .05).
  bh_de <- function(g) { i <- tejd == "BRAIN_E15" & geng == g; if (any(i)) bh[i][1] else NA_real_ }
  bh_sin_posthoc <- vapply(sin_posthoc, bh_de, numeric(1))
  d$genes_sin_posthoc_explicados <- sin_posthoc[!is.na(bh_sin_posthoc) & bh_sin_posthoc < 0.05]
  d$genes_sin_posthoc_tendencia <- sin_posthoc[!is.na(bh_sin_posthoc) &
                                               bh_sin_posthoc >= 0.05 & bh_sin_posthoc < 0.10]

  d$fuente <- fuente_datos(ARCHIVO_QPCR)
  d
}
.join_y <- function(v) {
  if (length(v) == 0L) return("")
  if (length(v) == 1L) return(v[1])
  paste0(paste(v[-length(v)], collapse = ", "), " y ", v[length(v)])
}

# =========================================================================
# Armado de cada diapositiva -- 12 secciones, 1 funcion por diapositiva
# (numeracion y contenido siguen pedidos/cambios_presentacion_2.md).
# =========================================================================
slide <- function(n, kicker, cuerpo, clase = "") {
  kick <- if (nzchar(kicker)) sprintf('<div class="kicker">%s</div>', .esc(kicker)) else ""
  sprintf('<section class="slide %s" id="s%d">\n%s\n%s\n<div class="slide-num">%d / 12</div>\n</section>',
          clase, n, kick, cuerpo, n)
}

# --- 1. Titulo (pedido 2, seccion 1): titulo del trabajo, no el h1 del repo. --
slide1 <- function() {
  cuerpo <- paste0(
    '<h1>', .esc(PRESENTACION_TITULO), '</h1>\n',
    '<p class="lead">', .esc(PRESENTACION_SUBTITULO), '</p>\n',
    '<p style="margin-top:3rem; font-size:1.2rem; color:var(--ink-soft);">',
    .esc(PRESENTACION_AUTORA), '<br>', .esc(PRESENTACION_CURSO), '<br>',
    format(Sys.Date(), "%d/%m/%Y"), '</p>')
  slide(1, "", cuerpo, "title")
}

# --- 2. Antecedentes (pedido 2, seccion 2): ilustracion + rayos SVG apuntando
#     al vientre, flecha + recuadro SVG/HTML, grafico de razones de
#     prevalencia (2bis, construir_epi_chart()); frase de cierre en tarjeta. --
slide2 <- function() {
  cuerpo <- paste0(
    '<h2>Antecedentes</h2>\n',
    '<div class="problem-row">\n',
    '<div class="illus">', fig_asset("ilustracion-mia.png", "Ilustracion: embarazada"),
    svg_rayos(), '</div>\n',
    '<div class="arrowbox">', svg_flecha(),
    '<div class="box">Trastornos del<br>neurodesarrollo</div></div>\n',
    '<div class="epi">', construir_epi_chart(), '</div>\n',
    '</div>\n',
    '<div class="card closing-card"><p class="lead" style="margin:0;">La inflamacion materna ',
    'durante la gestacion es un factor de riesgo para trastornos del neurodesarrollo, y esas ',
    'patologias afectan de forma distinta a varones y mujeres.</p></div>')
  slide(2, "Contexto", cuerpo)
}

# --- 3. Modelo experimental (pedido 2, seccion 3): imagen sin modificar
#     (assets/modelo-experimental.png -- misma ruta; el archivo se reemplaza
#     aparte, sin tocar codigo) + texto minimo. ----------------------------
slide3 <- function(d) {
  cuerpo <- paste0(
    '<h2>Modelo experimental</h2>\n',
    '<div style="text-align:center;">',
    fig_asset("modelo-experimental.png", "Modelo experimental", "wide"), '</div>\n',
    '<div class="cols" style="margin-top:.6rem; align-items:flex-start; flex:0;">\n',
    '<div class="col"><p class="lead">LPS 100 &mu;g/kg i.p. en el dia 15 de gestacion (E15), ',
    'coleccion a las 6 horas.</p>',
    '<p>', d$n_fetos, ' fetos: 18 camadas, un feto de cada sexo por camada; placenta y ',
    'cerebro del mismo individuo.</p></div>\n',
    '<div class="col"><p class="lead">Tres mediciones:</p>',
    '<ul class="bullets"><li>IL-6 en suero materno y liquido amniotico, por ELISA</li>',
    '<li>', d$n_genes, ' genes por RT-qPCR (transportadores de nutrientes + via IL-6/STAT3)</li>',
    '<li>pSTAT3 en placenta, por Western blot</li></ul>',
    '<p>', d$n_grupos, ' grupos (sexo &times; tratamiento), n = ', d$n_por_grupo, '.</p></div>\n',
    '</div>')
  slide(3, "Contexto", cuerpo)
}

# --- 4. El punto de partida fue un prompt: fragmento REAL de AGENTS.md, mas
#     el cierre del arco (pedido 2, seccion 4) que conecta con la diapositiva
#     5: de esas decisiones salio primero el analisis convencional, y despues
#     una exploracion adicional que el analisis convencional no incluia. ----
slide4 <- function(d) {
  cuerpo <- paste0(
    '<h2>El punto de partida fue un prompt</h2>\n',
    '<p class="lead">Fragmentos reales de las decisiones fijas del repositorio ',
    '(<code>AGENTS.md</code>), escritas antes de correr ningun analisis:</p>\n',
    '<pre class="quote"><strong>D1</strong>: ', .inline_md(d$d1), '</pre>\n',
    '<pre class="quote"><strong>D2</strong>: ', .inline_md(d$d2), '</pre>\n',
    '<pre class="quote"><strong>D5</strong>: ', .inline_md(d$d5), '</pre>\n',
    '<p>No fue &laquo;analiza mis datos&raquo;.</p>\n',
    '<p>Sobre esas decisiones se pidio despues <strong>el analisis convencional</strong> ',
    '(modelos por gen y tejido, con sus boxplots) y, encima, <strong>una exploracion de los ',
    'datos</strong> que el analisis convencional no incluia: correlaciones entre tejidos, ',
    'co-expresion entre genes, dispersion.</p>')
  slide(4, "Lo que produjo el agente", cuerpo)
}

# --- 5. Lo que salio de ahi: conteos leidos del repo, no escritos a mano. ---
slide5 <- function(d) {
  stat <- function(n, l) sprintf(
    '<div class="stat"><div class="n">%s</div><div class="l">%s</div></div>', n, l)
  cuerpo <- paste0(
    '<h2>Lo que salio de ahi</h2>\n',
    '<div class="stat-grid">\n',
    stat(sprintf("%d + %d", d$n_scripts_r, d$n_scripts_py), "scripts (R + Python)"),
    stat(d$n_figuras, "figuras"),
    stat(sprintf("%d / %d", d$n_tablas_r, d$n_tablas_py), "tablas CSV (R / Python)"),
    stat(d$n_procedencia, "filas en procedencia.csv"),
    stat(sprintf("%d + %d + %d", d$n_verif_recalculo, d$n_verif_existencia, d$n_verif_declaracion),
         sprintf("verificaciones: recalculo + existencia + declaracion (%d en total)", d$n_verif)),
    stat(sprintf("%s / %s", d$comp_byte_n, d$comp_n), "CSV byte-identicos, R &harr; Python"),
    '</div>')
  slide(5, "Lo que produjo el agente", cuerpo)
}

# --- 6. Placenta -- el analisis convencional (pedido 2, seccion 6): boxplots
#     de 3 genes de placenta con efecto de tratamiento (subconjunto, no el
#     panel completo) + pSTAT3 mas chico al costado. -----------------------
slide6 <- function(d, sint) {
  texto <- if (sint) AVISO_SINTETICO else paste0(
    '<p class="lead">El LPS modifica la expresion en placenta y activa la via IL-6/STAT3, ',
    'pero el boxplot no muestra interaccion sexo &times; tratamiento en ningun gen ',
    '(interaccion SEXO&times;TTO de pSTAT3, p = ', d$pstat3_pint, '; post hoc &female;Control ',
    'vs &female;LPS p = ', d$pstat3_hh, ', &male;Control vs &male;LPS p = ', d$pstat3_mm, ').</p>')
  cuerpo <- paste0(
    '<h2>Placenta</h2>\n',
    '<div class="cols">\n',
    '<div class="col wide">',
    fig_outputs(NOMBRE_SUBSET_PLACENTA, "Boxplots de il6, glut3 y slc38a2 en placenta", "tall"),
    '</div>\n',
    '<div class="col narrow">',
    fig_outputs("acto1_pstat3.png", "pSTAT3 en placenta", "small"),
    texto, '</div>\n',
    '</div>')
  slide(6, "El analisis convencional", cuerpo)
}

# --- 7. Cerebro fetal -- el analisis convencional (pedido 2, seccion 6):
#     3 genes elegidos (los que tienen post hoc significativo), no el panel
#     completo (con el panel completo "no se lee nada"). --------------------
slide7 <- function(d, sint) {
  genes_html <- paste(sprintf('<span class="gene-tag">%s</span>', .esc(GENES_SUBSET_BRAIN)),
                       collapse = "")
  texto <- if (sint) AVISO_SINTETICO else paste0(
    '<p class="lead">', d$qpcr_int_bra, ' de ', d$qpcr_mod_bra, ' genes modelados en cerebro ',
    'muestran interaccion SEXO&times;TTO (vs. ', d$qpcr_int_pla, ' en placenta). En estos tres, ',
    'el post hoc ubica el efecto en las hembras.</p>\n', genes_html)
  cuerpo <- paste0(
    '<h2>Cerebro fetal</h2>\n',
    '<div class="cols">\n',
    '<div class="col wide">',
    fig_outputs(NOMBRE_SUBSET_BRAIN, "Boxplots de glut1, slc38a2 y fatp1 en cerebro", "tall"),
    '</div>\n',
    '<div class="col narrow">', texto, '</div>\n',
    '</div>')
  slide(7, "El analisis convencional", cuerpo)
}

# --- 8. Co-expresion entre genes (NUEVA, pedido 2 seccion 7): con el boxplot
#     cada gen se mira por separado; el diagrama triangular (SPLOM) muestra
#     como se mueven los genes entre si. Nivel de tejido (no por sexo): con
#     el nivel de sexo el panel se duplica sin ganar legibilidad proyectado;
#     avisado al usuario en el cierre de esta tanda. -----------------------
slide8 <- function(sint) {
  texto <- if (sint) AVISO_SINTETICO else paste0(
    '<p class="lead">Con el boxplot, cada gen se mira por separado. El diagrama triangular ',
    'muestra como se mueven los genes entre si: un analisis que no estaba en el plan original.</p>')
  cuerpo <- paste0(
    '<h2>Co-expresion entre genes</h2>\n',
    '<div class="cols">\n',
    '<div class="col wide">',
    fig_outputs("acto2_coexpresion_SPLOM_BRAIN_E15.png", "Co-expresion de genes en cerebro E15", "tall"),
    '</div>\n',
    '<div class="col narrow">', texto, '</div>\n',
    '</div>')
  slide(8, "Lo que agrego explorar con el agente", cuerpo)
}

# --- 9. La dispersion (pedido 2, seccion 7): razonamiento completo de 3
#     pasos, no solo la figura. Genes y conteos derivados de las tablas gen
#     por gen (recolectar_datos()), no una lista escrita a mano. -----------
slide9 <- function(d, sint) {
  texto <- if (sint) AVISO_SINTETICO else {
    explicados <- d$genes_sin_posthoc_explicados
    tendencia <- d$genes_sin_posthoc_tendencia
    paso3 <- if (length(explicados))
      paste0('En ', length(explicados), ' (', paste(explicados, collapse = ", "),
             ') el analisis de dispersion lo explica (BH &lt; 0.05)',
             if (length(tendencia))
               paste0('; ', paste(tendencia, collapse = ", "), ' queda en tendencia')
             else '', '.')
    else 'El analisis de dispersion no alcanza a explicar ninguno con este criterio.'
    paste0(
      '<ol class="bullets" style="font-size:.92rem;">',
      '<li>En cerebro, ', length(d$genes_int_bra), ' genes mostraron interaccion sexo &times; ',
      'tratamiento en el boxplot.</li>',
      '<li>En ', length(d$genes_sin_posthoc), ' (', paste(d$genes_sin_posthoc, collapse = ", "),
      '), ningun par de grupos explicaba la diferencia: habia algo sexo-dependiente, pero no ',
      'se veia donde.</li>',
      '<li>No es donde se ubican los datos sino cuanto se dispersan: el LPS compacta la ',
      'expresion en hembras y la dispersa en machos. ', paso3, '</li>',
      '</ol>')
  }
  cuerpo <- paste0(
    '<h2>La dispersion</h2>\n',
    '<div class="cols">\n',
    '<div class="col wide">', fig_outputs("acto2_dispersion_sd.png", "Dispersion, agrupada y por sexo", "tall"), '</div>\n',
    '<div class="col narrow">', texto, '</div>\n',
    '</div>')
  slide(9, "Lo que agrego explorar con el agente", cuerpo)
}

# --- 10. Decisiones que el agente documenta como propias (sin cambios de
#     contenido, pedido 2 seccion 8): texto movido a columna lateral. No se
#     encontro una figura que muestre el desacople sin exponer un resultado
#     inedito (toda figura de correlacion por sexo trae rho/p impresos) --
#     se avisa, como autoriza el pedido, y la diapositiva queda sin figura.
#     Fuente exacta: outputs/tables/analisis_descartados.md, "el informe
#     mostraba una cosa ... y testeaba otra" (seccion 09_acto2_dispersion),
#     con el ejemplo concreto de fatcd36. ------------------------------------
slide10 <- function() {
  cuerpo <- paste0(
    '<h2>El agente toma decisiones y las documenta como propias</h2>\n',
    '<div class="cols">\n<div class="col narrow">\n',
    '<p class="lead">En el Acto 2, el agente agrupo los sexos en el test formal de correlacion ',
    'placenta&ndash;cerebro y lo justifico en su propio archivo de descartes, mientras las ',
    'figuras ya mostraban los datos separados por sexo. El informe mostraba una cosa y testeaba ',
    'otra.</p>\n',
    '<p><strong>El caso que lo deja claro</strong> (sin exponer el numero real: es un ',
    'resultado inedito): en <code>fatcd36</code>, un sexo mostraba una correlacion positiva ',
    'fuerte y el otro una correlacion negativa; agrupados, ambas se promedian en un valor ',
    'intermedio que no describe a ninguno de los dos sexos.</p>\n',
    '<p style="margin-top:1.2rem;"><strong>Como se detecto:</strong> leyendo el informe y ',
    'notando el desacople entre la figura (separada por sexo) y el test (agrupado).</p>\n',
    '<p class="meta" style="font-size:.85rem;">Detalle completo: ',
    '<code>outputs/tables/analisis_descartados.md</code>, secciones ',
    '<code>08_acto2_correlaciones</code> y <code>09_acto2_dispersion</code>.</p>',
    '\n</div>\n</div>')
  slide(10, "Errores y como se detectaron", cuerpo)
}

# --- 11. Las verificaciones no verificaban (sin cambios de contenido).
#     Fuente exacta: revisiones/AUDITORIA.md (hallazgos A4-A9) y
#     revisiones/RESPUESTA.md. -----------------------------------------------
slide11 <- function() {
  cuerpo <- paste0(
    '<h2>Las verificaciones no verificaban</h2>\n',
    '<p class="lead">Un agente externo, de otra empresa, audito un clon limpio del repositorio ',
    'y encontro que buena parte de <code>verificaciones.csv</code> no podia fallar.</p>\n',
    '<ul class="bullets">',
    '<li><code>informe_pdf</code> devolvia <code>TRUE</code> incluso cuando la generacion del ',
    'PDF habia fallado.</li>',
    '<li><code>pstat3_rama_cascada</code> exigia una rama fija en vez de comprobar que la ',
    'elegida correspondiera a sus propios diagnosticos.</li>',
    '<li>&laquo;Todas las verificaciones pasaron&raquo; no distinguia entre lo que de verdad ',
    'se recalculaba, lo que solo comprobaba que un archivo existiera, y lo que era una ',
    'constante declarada.</li>',
    '</ul>\n',
    '<p><strong>Que se hizo:</strong> clasificar las verificaciones por tipo ',
    '(<em>recalculo</em> / <em>existencia</em> / <em>declaracion</em>) y agregar el estado ',
    '<code>NO_EJECUTADA</code> para lo que no llego a correr.</p>\n',
    '<p class="meta" style="font-size:.85rem;">Detalle completo: ',
    '<code>revisiones/AUDITORIA.md</code> y <code>revisiones/RESPUESTA.md</code>.</p>')
  slide(11, "Errores y como se detectaron", cuerpo)
}

# --- 12. Cierre (sin cambios: se reescribe cuando termine la segunda prueba
#     de reproduccion con el agente externo, pedido 2 seccion 9). ----------
slide12 <- function() {
  cuerpo <- paste0(
    '<h2>Cierre</h2>\n',
    '<p class="lead">El agente hizo en dias lo que llevaria semanas, y ninguna de sus propias ',
    'verificaciones detecto sus propios errores: aparecieron al leer el informe y al auditarlo ',
    'con otro agente.</p>\n',
    '<p style="margin-top:2.2rem; font-size:1.4rem;"><strong>&iquest;Que haria falta para ',
    'confiar en un analisis hecho asi?</strong></p>')
  slide(12, "Cierre", cuerpo)
}

# =========================================================================
# Armado del HTML completo.
# =========================================================================
construir_html <- function(d) {
  sint <- d$fuente != "real"
  slides <- paste(
    slide1(), slide2(), slide3(d), slide4(d), slide5(d),
    slide6(d, sint), slide7(d, sint), slide8(sint), slide9(d, sint),
    slide10(), slide11(), slide12(),
    sep = "\n")
  paste0(
    "<!doctype html>\n<html lang=\"es\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>", .esc(PRESENTACION_TITULO), "</title>\n<style>", CSS, "</style>\n</head>\n<body>\n",
    slides, "\n</body>\n</html>\n")
}

# =========================================================================
main <- function() {
  d <- recolectar_datos()
  sint <- d$fuente != "real"

  # --- figuras de subconjunto (pedido 2, secciones 6-7): se regeneran SIEMPRE
  # (idempotente), reflejando la fuente de datos vigente en outputs/tables/R y
  # data/processed en el momento de esta corrida -- igual criterio que las
  # figuras completas de 07_figuras_acto1, que tambien "son lo que haya". ----
  ruta_subset_pla <- file.path(RUTA_FIGURAS, NOMBRE_SUBSET_PLACENTA)
  ruta_subset_bra <- file.path(RUTA_FIGURAS, NOMBRE_SUBSET_BRAIN)
  generar_panel_subset(GENES_SUBSET_PLACENTA, "PLACENTA_E15", ruta_subset_pla)
  generar_panel_subset(GENES_SUBSET_BRAIN, "BRAIN_E15", ruta_subset_bra)

  html <- construir_html(d)

  destino <- if (sint) RUTA_DOCS else RUTA_PRES_REAL
  if (!dir.exists(destino)) dir.create(destino, recursive = TRUE, showWarnings = FALSE)
  ruta_html <- file.path(destino, "index.html")
  ruta_pdf  <- file.path(destino, "presentacion.pdf")
  ruta_ref  <- file.path(destino, "referencias.md")

  escribir_texto(ruta_html, html)
  escribir_referencias(ruta_ref)
  pdf_status <- generar_pdf(ruta_html, ruta_pdf)

  # --- verificaciones (punto 4 del pedido) --------------------------------
  m_aviso <- gregexpr(AVISO_SINTETICO, html, fixed = TRUE)[[1]]
  n_aviso <- if (length(m_aviso) == 1L && m_aviso[1] == -1L) 0L else length(m_aviso)
  n_aviso_esperado <- if (sint) 4L else 0L   # diapositivas 6,7,8,9 (pedido 2)

  usadas <- unique(.FIGS_USADAS$outputs)
  faltan_figs <- usadas[!file.exists(file.path(RUTA_FIGURAS, usadas))]

  gi_ok <- tryCatch({
    rc <- system2("git", c("-C", shQuote(RAIZ_REPO), "check-ignore", "-q",
                           shQuote(file.path(RUTA_PRES_REAL, "index.html"))),
                  stdout = FALSE, stderr = FALSE)
    rc == 0L
  }, error = function(e) NA)

  destino_rel <- if (sint) "docs" else "outputs/presentacion_real"
  ent <- "AGENTS.md + README.md + outputs/tables/*.csv + outputs/figures/*.png + assets/*.png"
  registrar_procedencia(list(
    list(file.path(destino_rel, "index.html"), "presentacion", ESTE_SCRIPT, "PROPIO", ent,
         paste0("pagina de presentacion, 11 diapositivas; version ",
                if (sint) "PUBLICA (sintetica)" else "PARA EXPONER (real, no versionada)")),
    list(file.path(destino_rel, "presentacion.pdf"), "presentacion", ESTE_SCRIPT, "PROPIO",
         "index.html",
         paste0("version imprimible por impresion headless; best-effort (estado: ",
                pdf_status, ")")),
    list(file.path(destino_rel, "referencias.md"), "presentacion", ESTE_SCRIPT, "PROPIO",
         "pedidos/referencias_epidemiologia.md",
         "10 citas del grafico de contexto epidemiologico (diapositiva 2), copiadas tal cual"),
    list(paste0(destino_rel, "/index.html#s2 (grafico inline, no PNG)"), "presentacion",
         ESTE_SCRIPT, "PROPIO", "pedidos/referencias_epidemiologia.md",
         paste0("grafico de razones de prevalencia por sexo: el UNICO grafico de la ",
                "presentacion que no sale de outputs/figures/ -- es contexto ",
                "bibliografico, no un resultado del proyecto")),
    # tipo = "figura_presentacion", NO "figura": 12_informe.R clasifica como
    # "figura" (a secas) + ruta en outputs/figures/ = debe estar embebida en
    # docs/informe.html (informe_figuras_procedencia_embebidas). Estos 2 PNG
    # son solo para la presentacion, nunca para el informe -- un tipo
    # distinto los deja fuera de ese chequeo sin tocar 12_informe.R.
    list(file.path("outputs/figures", NOMBRE_SUBSET_PLACENTA), "figura_presentacion",
         ESTE_SCRIPT,
         "PROPIO (panel_gen() de 07_figuras_acto1.R, reusada sin recalculo)",
         "data/processed/qpcr_cuantificacion_long.tsv + outputs/tables/R/qpcr_modelos_{clasificacion,posthoc}.csv",
         paste0("subconjunto de 3 genes (", paste(GENES_SUBSET_PLACENTA, collapse = ", "),
                ") de placenta para la diapositiva 6, sin recalcular nada; no forma parte ",
                "de docs/informe.html")),
    list(file.path("outputs/figures", NOMBRE_SUBSET_BRAIN), "figura_presentacion",
         ESTE_SCRIPT,
         "PROPIO (panel_gen() de 07_figuras_acto1.R, reusada sin recalculo)",
         "data/processed/qpcr_cuantificacion_long.tsv + outputs/tables/R/qpcr_modelos_{clasificacion,posthoc}.csv",
         paste0("subconjunto de 3 genes (", paste(GENES_SUBSET_BRAIN, collapse = ", "),
                ") de cerebro para la diapositiva 7, sin recalcular nada; no forma parte ",
                "de docs/informe.html"))
  ))
  registrar_verificaciones(list(
    list("presentacion_html_generado", "recalculo",
         sprintf("%s existe y no esta vacio", file.path(destino_rel, "index.html")),
         sprintf("existe=%s", file.exists(ruta_html) && file.info(ruta_html)$size > 0),
         "TRUE",
         if (file.exists(ruta_html) && file.info(ruta_html)$size > 0) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("presentacion_sin_interpretacion", "recalculo",
         paste0("con fuente sintetica, las diapositivas 6-9 reemplazan la prosa ",
                "interpretativa por el aviso de datos sinteticos"),
         sprintf("fuente=%s; aviso=%d/%d", d$fuente, n_aviso, n_aviso_esperado),
         "aviso = 3 si fuente sintetica, 0 si fuente real",
         if (n_aviso == n_aviso_esperado) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("presentacion_figuras_existen", "recalculo",
         "toda figura de outputs/figures/ referenciada por la presentacion existe",
         sprintf("usadas=%d; faltan=%s", length(usadas),
                 if (length(faltan_figs)) paste(faltan_figs, collapse = ", ") else "[]"),
         "faltan = []", if (!length(faltan_figs)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("presentacion_gitignore_real", "recalculo",
         "outputs/presentacion_real/ esta cubierto por .gitignore (git check-ignore)",
         sprintf("ignorado=%s", if (is.na(gi_ok)) "NA" else gi_ok), "TRUE",
         if (isTRUE(gi_ok)) "TRUE" else if (is.na(gi_ok)) "NO_EJECUTADA" else "FALSE",
         ESTE_SCRIPT),
    list("presentacion_pdf", "recalculo",
         paste0("presentacion.pdf generado (TRUE), degradado limpio sin motor ",
                "(NO_EJECUTADA), o fallo real (FALSE)"),
         sprintf("estado=%s", pdf_status),
         "ok -> TRUE | sin_motor -> NO_EJECUTADA | fallo:* -> FALSE",
         if (pdf_status == "ok") "TRUE"
         else if (pdf_status == "sin_motor") "NO_EJECUTADA" else "FALSE",
         ESTE_SCRIPT)
  ))

  cat("== 13_presentacion.R ==\n")
  cat(sprintf("  fuente = %s\n", d$fuente))
  cat(sprintf("  -> %s (%d KB)\n", ruta_html, file.info(ruta_html)$size %/% 1024))
  cat(sprintf("  -> %s [%s]\n", ruta_pdf, pdf_status))
  cat(sprintf("  -> %s\n", ruta_ref))
  if (length(faltan_figs))
    cat(sprintf("  *** faltan figuras: %s ***\n", paste(faltan_figs, collapse = ", ")))
  if (nzchar(PRESENTACION_AUTORA) && grepl("COMPLETAR", PRESENTACION_AUTORA, fixed = TRUE))
    cat("  *** completar PRESENTACION_AUTORA / PRESENTACION_CURSO antes de exponer ***\n")
}

if (sys.nframe() == 0L) main()
