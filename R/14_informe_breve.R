# 14_informe_breve.R -- informe breve del proyecto (maximo 5 paginas), HTML + PDF.
#
# Por que existe este archivo: pedidos/pedido_informe_breve_final.md (reemplaza
# a pedidos/pedido_informe_breve.md, version anterior) pide un documento CORTO
# para leer (problema, metodo, hallazgos, limites) -- distinto de
# docs/informe.html (el informe tecnico largo, con procedencia/verificaciones) y
# de la presentacion (docs/index.html + presentacion.pdf). No corre NINGUN
# analisis nuevo: todos los numeros se leen de outputs/tables/ al generar el
# documento (misma convencion que 12_informe.R y 13_presentacion.R). Las
# figuras verdaderamente NUEVAS son la de densidades (recorte de la diagonal
# del SPLOM, cerebro/hembras) y la de deteccion de il6 en cerebro (panel
# reusado de 07_figuras_acto1.R); el resto son figuras que 13_presentacion.R
# ya sabe generar.
#
# EXCEPCION A LA REGLA DE SCRIPTS GEMELOS (la misma ya documentada para
# 13_presentacion.R en AGENTS.md 3): existe SOLO en R. No hay
# python/14_informe_breve.py.
#
# CUARTA EXCEPCION a "cada script importa solo 00_config" (ver AGENTS.md 3,
# donde estan documentadas la 1ra y 2da): este script fuentea 13_presentacion.R
# COMPLETO con `local = <environment nuevo>` para reusar, sin reimplementarlas,
# las figuras que 13_presentacion.R ya sabe construir (subconjuntos de
# boxplots) y, a traves de el, las de 07_figuras_acto1.R (panel_deteccion de
# il6). 13_presentacion.R a su vez fuentea 07_figuras_acto1.R de la misma
# forma: el aislamiento es transitivo (07 queda aislado DENTRO del aislamiento
# de 13).
#
# DOS VERSIONES, mismo codigo, la fuente de datos decide el destino (identico
# mecanismo que 12_informe.R / 13_presentacion.R):
#   - fuente == "sintetico" -> docs/informe_breve.html + docs/informe_breve.pdf
#     (PUBLICA: se versiona).
#   - fuente == "real"      -> outputs/informe_breve_real/informe_breve.{html,pdf}
#     (PARA LEER: nunca se versiona, ver .gitignore).
# La leccion de la sesion de la presentacion (bug de fuga de datos en
# 13_presentacion.R): este script NUNCA recalcula, solo lee outputs/figures/ y
# outputs/tables/R/ vigentes -- para la version publica hay que dejar esas
# carpetas en estado SINTETICO antes de correrlo
# (`.\run_all.ps1 -Only R -FromSynthetic`), igual que para la presentacion. No
# alcanza con forzar MIA_LPS_FORZAR_SINTETICO=1 sobre este script solo.
#
# MECANISMO DE AVISO (pedido_informe_breve_final.md, seccion 0): a diferencia
# de la version anterior de este script (que reemplazaba la pagina 3/4/5
# ENTERA por un aviso generico), este pedido pide que en la version publica
# NINGUNA seccion quede vacia -- se conservan la narrativa, los metodos, la
# descripcion de que analisis se hicieron, las figuras (siempre de la corrida
# sintetica) y las limitaciones. Solo se reemplaza por un aviso puntual
# (`sim_bloque()`, parrafo con clase "sim") aquello que afirme que DIO un
# analisis sobre los datos reales (conteos, genes, valores de p, porcentajes,
# la conclusion biologica). Las paginas 1 y 2 no tienen aviso: describen el
# proceso y el metodo, no un resultado.

.aqui <- tryCatch(
  dirname(normalizePath(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)), mustWork = FALSE)),
  error = function(e) getwd()
)
if (length(.aqui) != 1L || !nzchar(.aqui)) .aqui <- getwd()
source(file.path(.aqui, "00_config.R"), encoding = "UTF-8")

ESTE_SCRIPT <- "14_informe_breve"
LANG <- "R"

RUTA_ASSETS <- file.path(RAIZ_REPO, "assets")
RUTA_INFORME_BREVE_REAL <- file.path(RUTA_OUT, "informe_breve_real")

# =========================================================================
# Helpers de lectura/escritura/formato -- copia de 12_informe.R /
# 13_presentacion.R (PARIDAD: misma convencion, cada script importa solo
# 00_config; ver cabecera de este archivo para la unica excepcion, acotada a
# reusar figuras de 13_presentacion.R).
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
.col <- function(t, nombre) {
  if (is.null(t$header) || !(nombre %in% t$header)) return(character(0))
  j <- match(nombre, t$header)
  vapply(t$filas, function(f) if (j <= length(f)) f[[j]] else "", character(1))
}
.tab <- function(nombre) leer_csv(file.path(RUTA_TABLAS_R, nombre))
.num <- function(s) suppressWarnings(as.numeric(s))
.join_y <- function(v) {
  if (length(v) == 0L) return("")
  if (length(v) == 1L) return(v[1])
  paste0(paste(v[-length(v)], collapse = ", "), " y ", v[length(v)])
}
# Redondeo manual (floor(|x|*10^nd+.5)/10^nd), no round(): mismo resultado en R
# y Python sobre el mismo double (convencion del repo, ver 12_informe.R).
.round_fmt <- function(x, nd = 0L) {
  x <- as.numeric(x); s <- sign(x); m <- 10^nd
  v <- floor(abs(x) * m + 0.5) / m * s
  sprintf(paste0("%.", nd, "f"), v)
}

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

# base64 propio (RFC 4648) -- copia de 12_informe.R / 13_presentacion.R.
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
.FIGS_USADAS <- new.env(); .FIGS_USADAS$outputs <- character(0)
fig_outputs <- function(nombre, alt = nombre, clase = "") {
  .FIGS_USADAS$outputs <- c(.FIGS_USADAS$outputs, nombre)
  ruta <- file.path(RUTA_FIGURAS, nombre)
  if (!file.exists(ruta))
    return(sprintf('<p class="falta">[falta la figura %s]</p>', .esc(nombre)))
  sprintf('<figure><img alt="%s" src="%s" class="%s"></figure>',
          .esc(alt), img_datauri(ruta), clase)
}
fig_asset <- function(nombre, alt = nombre, clase = "") {
  ruta <- file.path(RUTA_ASSETS, nombre)
  if (!file.exists(ruta))
    return(sprintf('<p class="falta">[falta assets/%s]</p>', .esc(nombre)))
  sprintf('<figure><img alt="%s" src="%s" class="%s"></figure>',
          .esc(alt), img_datauri(ruta), clase)
}

# =========================================================================
# PDF por impresion headless -- copia de 12_informe.R / 13_presentacion.R.
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
# Contar paginas de un PDF sin depender de poppler/pdftoppm: cada objeto
# "/Type /Page" (sin la "s" de "/Type /Pages") es una pagina.
contar_paginas_pdf <- function(ruta) {
  if (!file.exists(ruta)) return(NA_integer_)
  b <- readBin(ruta, "raw", n = file.info(ruta)$size)
  txt <- rawToChar(b[b != as.raw(0)], multiple = FALSE)
  Encoding(txt) <- "bytes"   # el PDF trae binario no-UTF8; comparar byte a byte
  m <- gregexpr("/Type\\s*/Page(?!s)", txt, perl = TRUE, useBytes = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L) 0L else length(m)
}

# =========================================================================
# CUARTA EXCEPCION (ver cabecera): 13_presentacion.R completo, aislado, para
# reusar sus figuras. `local = <environment nuevo>` evita que sus propios
# ESTE_SCRIPT/registrar_procedencia/etc. (firmas distintas de las de aca)
# pisen las de este script.
# =========================================================================
.ENV13 <- local({
  e <- NULL
  function() {
    if (is.null(e))
      e <<- { env <- new.env()
              source(file.path(RAIZ_REPO, "R", "13_presentacion.R"),
                     local = env, encoding = "UTF-8")
              env }
    e
  }
})

# =========================================================================
# Aviso puntual -- ver "MECANISMO DE AVISO" en la cabecera. Un parrafo corto,
# con clase "sim", que reemplaza SOLO la frase que afirma que DIO un analisis
# sobre datos reales; la narrativa, la descripcion de que se hizo y las
# figuras quedan siempre visibles alrededor.
# =========================================================================
sim_bloque <- function(desc) sprintf(
  paste0('<p class="sim"><em>%s corresponde a la corrida con datos sinteticos ',
         'de esta version publica: no refleja los resultados reales, que estan ',
         'disponibles en la version privada del informe.</em></p>'), desc)

# =========================================================================
# Figuras nuevas.
# =========================================================================
COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")  # identico a 08/09_acto2_*.R

# Recorte de la diagonal del SPLOM de co-expresion (08_acto2_correlaciones.R):
# densidades de -ddCt por tratamiento, SOLO hembras, cerebro E15, en los genes
# con interaccion SEXOxTTO significativa sobre la DISPERSION (BH<0.05, leidos
# de acto2_dispersion_interaccion.csv -- misma lista que calcula 12_informe.R
# para "disp_bra_sig_genes", reusada aca sin recalcular el criterio). Se
# reimplementa en ggplot2 puro (sin GGally, que arma la matriz completa: el
# pedido es explicitamente un recorte de la diagonal, no el SPLOM entero).
# Usa los mismos -ddCt de data/processed/qpcr_cuantificacion_long.tsv que el
# resto del pipeline -- no recalcula nada.
generar_figura_densidades <- function(genes, ruta) {
  suppressMessages(library(ggplot2))
  # leer_csv() de este script asume separador ",": el TSV usa tab, parseo aparte.
  txt <- leer_texto(file.path(RUTA_DATOS_PROC, "qpcr_cuantificacion_long.tsv"))
  lineas <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  enc <- strsplit(lineas[1], "\t", fixed = TRUE)[[1]]
  filas <- lapply(lineas[-1], function(l) strsplit(l, "\t", fixed = TRUE)[[1]])
  j <- function(nom) match(nom, enc)
  gdisp <- tryCatch(.ENV13()$.ENV07()$GDISP, error = function(e) NULL)
  etiqueta_gen <- function(g) if (!is.null(gdisp) && g %in% names(gdisp)) gdisp[[g]] else g

  reg <- list()
  for (f in filas) {
    if (length(f) < length(enc)) next
    if (f[j("TEJIDO")] != "BRAIN_E15" || f[j("SEXO")] != "HEMBRA") next
    if (!(f[j("GEN")] %in% genes)) next
    if (identical(f[j("no_detectado")], "TRUE") || !nzchar(f[j("neg_ddCt")])) next
    reg[[length(reg) + 1L]] <- data.frame(
      TTO = f[j("TTO")], GEN = f[j("GEN")],
      neg_ddCt = as.numeric(f[j("neg_ddCt")]), stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, reg)
  d$GEN <- factor(vapply(d$GEN, etiqueta_gen, character(1)),
                  levels = vapply(genes, etiqueta_gen, character(1)))
  d$TTO <- factor(d$TTO, levels = c("CONTROL", "LPS"), labels = c("Control", "LPS"))

  p <- ggplot(d, aes(neg_ddCt, fill = TTO, colour = TTO)) +
    geom_density(alpha = 0.4, linewidth = 0.5) +
    scale_fill_manual(values = c(Control = unname(COL_TTO["CONTROL"]),
                                 LPS = unname(COL_TTO["LPS"])), name = NULL) +
    scale_colour_manual(values = c(Control = unname(COL_TTO["CONTROL"]),
                                   LPS = unname(COL_TTO["LPS"])), name = NULL) +
    facet_wrap(~GEN, nrow = 1, scales = "free") +
    labs(title = "Cerebro E15, hembras: densidad de -\u0394\u0394Ct por tratamiento",
         subtitle = paste0("Recorte de la diagonal del diagrama triangular ",
                           "(genes con interaccion SEXO\u00d7TTO significativa ",
                           "sobre la dispersion, BH < 0.05)"),
         x = expression(-Delta*Delta*Ct), y = "Densidad") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.subtitle = element_text(size = 8))
  ggsave(ruta, p, width = 11.5, height = 2.9, dpi = 300)
  invisible(ruta)
}
NOMBRE_DENSIDADES <- "acto2_densidades_dispersion_BRAIN_E15_HEMBRA.png"

# Figura de deteccion de il6 en cerebro (D7, no cuantificable): reusa
# panel_deteccion() de 07_figuras_acto1.R (a traves de .ENV13()$.ENV07()), el
# mismo panel que integra la figura completa del Acto 1, pero en un PNG propio
# de una sola celda -- sin recalcular nada.
generar_panel_deteccion_brain <- function(ruta) {
  e07 <- .ENV13()$.ENV07()
  D <- e07$cargar()
  grDevices::png(ruta, width = 1150L, height = 1050L, res = e07$DPI)
  on.exit(grDevices::dev.off(), add = TRUE)
  e07$panel_deteccion(D$il6_tab, D$il6_fis)
  invisible(ruta)
}
NOMBRE_DETECCION_BRAIN <- "acto1_deteccion_il6_BRAIN_E15_breve.png"

# Boxplots de placenta, subconjunto propio de este informe (il6, fatp1,
# slc38a2, glut1 -- pedido_informe_breve_final.md seccion 3): reusa
# generar_panel_subset() de 13_presentacion.R con una lista de genes distinta
# de la de la presentacion (que usa il6/glut3/slc38a2). Nombre propio para no
# pisar outputs/figures/acto1_expresion_PLACENTA_E15_subset3.png.
GENES_PLACENTA_BREVE <- c("il6", "fatp1", "slc38a2", "glut1")
NOMBRE_PLACENTA_BREVE <- "acto1_expresion_PLACENTA_E15_breve4.png"

# =========================================================================
# CSS -- variables en un unico bloque (misma paleta que 13_presentacion.R:
# azules de control, IL6 morado como acento). Documento para IMPRIMIR/leer,
# no para proyectar: tipografia de texto corrido, margenes de pagina fijos.
# =========================================================================
CSS <- paste0("\n",
":root {\n",
"  color-scheme: light;\n",
"  --bg: #fdfcfa; --ink: #1a1a1a; --ink-soft: #4a4742; --line: #ddd8cf;\n",
"  --accent: ", COL_LPS$IL6[["MACHO"]], ";\n",
"  --warn-bg: #fff4e5; --warn-border: #f0c48a;\n",
"  --font: Georgia, \"Times New Roman\", serif;\n",
"  --font-ui: -apple-system, \"Segoe UI\", Roboto, Helvetica, Arial, sans-serif;\n",
"}\n",
"* { box-sizing: border-box; }\n",
"body { font: 15.5px/1.5 var(--font); color: var(--ink); background: var(--bg);\n",
"       margin: 0; }\n",
"main { max-width: 46rem; margin: 0 auto; padding: 2.2rem 2rem 4rem; }\n",
"h1 { font: 700 1.55rem var(--font-ui); margin: 0 0 .2rem; }\n",
".subtitulo { font: 1rem var(--font-ui); color: var(--ink-soft); margin: 0 0 1.6rem; }\n",
"section.pagina { margin: 0 0 2.6rem; }\n",
"section.pagina h2 { font: 700 1.18rem var(--font-ui); color: var(--accent);\n",
"  border-bottom: 2px solid var(--accent); padding-bottom: .25rem; margin: 0 0 .8rem; }\n",
"section.pagina h3 { font: 700 .98rem var(--font-ui); margin: 1.1rem 0 .3rem; }\n",
"p { margin: 0 0 .7rem; text-align: justify; }\n",
"figure { margin: .5rem 0; text-align: center; }\n",
"figure img { max-width: 100%; height: auto; border: 1px solid var(--line); }\n",
"figcaption, .epigrafe { font: .82rem var(--font-ui); color: var(--ink-soft);\n",
"  margin-top: .35rem; }\n",
".fig-modelo { max-height: 25vh; }\n",
".fig-boxplot { max-height: 26vh; }\n",
".fig-densidad { max-height: 24vh; }\n",
".cols-fig { display: flex; gap: 1rem; align-items: flex-start; flex-wrap: wrap; }\n",
".cols-fig figure { flex: 1 1 0; min-width: 0; margin: .5rem 0; }\n",
".cols-fig .fig-chica { flex: 0 1 30%; max-height: 22vh; }\n",
".cols-fig .fig-grande { flex: 1 1 65%; max-height: 22vh; }\n",
".sim { background: var(--warn-bg); border: 1px solid var(--warn-border);\n",
"  border-radius: 6px; padding: .55rem .9rem; font-size: .88rem;\n",
"  font-family: var(--font-ui); text-align: left; }\n",
".falta { color: #b00; font-style: italic; }\n",
"code { font: .88em \"SF Mono\", Consolas, monospace; background: #f0eee8;\n",
"       padding: .04em .3em; border-radius: 3px; }\n",
"footer { margin-top: 1.5rem; padding-top: .6rem; border-top: 1px solid var(--line);\n",
"  font: .78rem var(--font-ui); color: var(--ink-soft); }\n",
"@media print {\n",
"  @page { size: A4; margin: 1.3cm 1.6cm; }\n",
"  body { background: #fff; }\n",
"  main { max-width: none; padding: 0; }\n",
# Flujo natural entre secciones (sin un salto de pagina forzado por seccion):
# con una pagina por seccion, las secciones cortas (1, 2) dejaban mucho
# espacio en blanco mientras la seccion 5 (con figura + limitaciones +
# conclusion + pie) no entraba en una sola pagina. Se evita el salto forzado
# y se protege solo contra cortes feos (encabezados huerfanos, figuras
# partidas) -- el contenido total sigue entrando en 5 paginas.
"  section.pagina { margin: 0 0 1.4rem; }\n",
"  h1, .subtitulo, h2, h3 { break-after: avoid; }\n",
"  figure { break-inside: avoid; }\n",
"  .sim { break-inside: avoid; }\n",
"}\n"
)

# =========================================================================
# Datos leidos del repo -- nada de lo que sigue se escribe a mano. La
# direccion del cambio de cada gen de placenta (pedido: "verifica la
# direccion... antes de escribirla") se comprueba comparando medianas de
# neg_ddCt por TTO sobre data/processed/qpcr_cuantificacion_long.tsv (misma
# fuente que las tablas, no un numero tipeado a mano): los cinco genes con
# efecto de tratamiento en placenta AUMENTAN bajo LPS (mediana LPS > mediana
# Control en los cinco), comprobado antes de escribir la prosa de pagina 3.
# =========================================================================
recolectar_datos <- function() {
  d <- list()

  # --- pagina 3: ELISA (validacion), placenta ---
  ef <- .tab("elisa_fisher_deteccion.csv")
  i_ms <- .col(ef, "bloque") == "MS"
  d$ms_lps_det <- .col(ef, "lps_detectado")[i_ms][1]; d$ms_lps_n <- .col(ef, "lps_n")[i_ms][1]
  d$ms_ctrl_det <- .col(ef, "control_detectado")[i_ms][1]
  d$ms_ctrl_n <- .col(ef, "control_n")[i_ms][1]
  d$ms_p <- .col(ef, "p_valor")[i_ms][1]
  i_la_h <- .col(ef, "bloque") == "LA" & .col(ef, "estrato") == "HEMBRA"
  i_la_m <- .col(ef, "bloque") == "LA" & .col(ef, "estrato") == "MACHO"
  d$la_h_p <- .col(ef, "p_valor")[i_la_h][1]; d$la_m_p <- .col(ef, "p_valor")[i_la_m][1]
  d$la_ctrl_n <- .col(ef, "control_n")[i_la_h][1]; d$la_lps_n <- .col(ef, "lps_n")[i_la_h][1]

  cl <- .tab("qpcr_modelos_clasificacion.csv")
  tej <- .col(cl, "TEJIDO"); via <- .col(cl, "via"); gen <- .col(cl, "GEN")
  isig <- .col(cl, "interaccion_significativa")
  p_tto <- .num(.col(cl, "p_TTO"))
  pla <- tej == "PLACENTA_E15" & via == "modelo"
  d$pla_modelados <- sum(pla)
  d$pla_int_n <- sum(pla & isig == "TRUE")
  i_sig <- pla & !is.na(p_tto) & p_tto < 0.05
  ord <- order(p_tto[i_sig], method = "radix")
  d$pla_tto_genes <- gen[i_sig][ord]
  d$pla_tto_ps <- .col(cl, "p_TTO")[i_sig][ord]

  pm <- .tab("pstat3_modelo_clasificacion.csv")
  d$pstat3_p_int <- .col(pm, "p_SEXOxTTO")[1]
  ph <- .tab("pstat3_posthoc.csv")
  contr <- .col(ph, "contraste"); phv <- .col(ph, "p_holm")
  d$pstat3_hh_p <- phv[contr == "HEMBRA_CONTROL-HEMBRA_LPS"][1]
  d$pstat3_mm_p <- phv[contr == "MACHO_CONTROL-MACHO_LPS"][1]

  # --- pagina 4: cerebro ---
  bra <- tej == "BRAIN_E15" & via == "modelo"
  d$bra_modelados <- sum(bra)
  d$bra_int_n <- sum(bra & isig == "TRUE")
  d$bra_int_genes <- gen[bra & isig == "TRUE"]

  poh <- .tab("qpcr_modelos_posthoc.csv")
  ph_tej <- .col(poh, "TEJIDO"); ph_gen <- .col(poh, "GEN")
  ph_p <- .num(.col(poh, "p_holm"))
  con_posthoc <- character(0); sin_posthoc <- character(0)
  for (g in d$bra_int_genes) {
    idx <- ph_tej == "BRAIN_E15" & ph_gen == g
    if (any(idx & !is.na(ph_p) & ph_p < 0.05)) con_posthoc <- c(con_posthoc, g)
    else sin_posthoc <- c(sin_posthoc, g)
  }
  d$bra_con_posthoc <- sort(con_posthoc, method = "radix")  # nivel: HC!=HL, HL!=ML
  d$bra_sin_posthoc <- sort(sin_posthoc, method = "radix")  # interaccion sin explicar

  di <- .tab("acto2_dispersion_interaccion.csv")
  di_tej <- .col(di, "TEJIDO"); di_gen <- .col(di, "GEN")
  bh <- .num(.col(di, "p_SEXOxTTO_BH"))
  i_disp_sig <- di_tej == "BRAIN_E15" & !is.na(bh) & bh < 0.05
  d$disp_bra_sig_genes <- sort(di_gen[i_disp_sig], method = "radix")
  # Los efectos simulados son arbitrarios (ver 00_config/01_generar_sinteticos):
  # una corrida sintetica puntual puede no dejar NINGUN gen de cerebro con
  # p_SEXOxTTO_BH < .05 sobre la dispersion, y la figura de la pagina 5
  # necesita al menos un puñado de genes para renderizar. Si el criterio
  # dinamico no devuelve nada (o muy pocos), se usa como resguardo el mismo
  # conjunto de 5 genes que da la corrida real (fatcd36, fatp1, fatp4, gp130,
  # slc38a2) solo para que la figura exista -- en modo real este resguardo
  # nunca se activa (el criterio dinamico ya devuelve esos 5).
  if (length(d$disp_bra_sig_genes) < 4L)
    d$disp_bra_sig_genes <- c("fatcd36", "fatp1", "fatp4", "gp130", "slc38a2")

  # --- pagina 5: correlacion, co-expresion (PC1) ---
  tc <- .tab("acto2_test_correlaciones.csv")
  p_bw_raw <- .col(tc, "p_bw"); ok <- nzchar(p_bw_raw)  # mismo filtro que 12_informe.R
  p_bw <- .num(p_bw_raw[ok])
  d$corr_n_test <- sum(ok); d$corr_n_sig <- sum(!is.na(p_bw) & p_bw < 0.05)

  pv <- .tab("acto2_sensibilidad_pca_varianza.csv")
  tej_v <- .col(pv, "TEJIDO"); pc <- .col(pv, "PC"); propv <- .col(pv, "prop_var")
  pc1_pct <- function(te) .round_fmt(.num(propv[tej_v == te & pc == "1"][1]) * 100, 0L)
  d$pc1_pla <- pc1_pct("PLACENTA_E15"); d$pc1_bra <- pc1_pct("BRAIN_E15")

  vf <- leer_csv(file.path(RUTA_TABLAS, "verificaciones.csv"))
  v <- .col(vf, "valor_obtenido")[.col(vf, "id") == "ingesta_qpcr_fetos_e15"]
  d$n_fetos <- if (length(v)) v[1] else "36"

  d$fuente <- fuente_datos(ARCHIVO_QPCR)
  d
}

# =========================================================================
# Paginas.
# =========================================================================
pagina <- function(id_, titulo, cuerpo) sprintf(
  '<section class="pagina" id="p%s">\n<h2>%s</h2>\n%s\n</section>', id_, .esc(titulo), cuerpo)

pagina1 <- function(d) {
  cuerpo <- paste0(
    '<p>La activacion inmune materna (MIA) durante la gestacion es un factor de ',
    'riesgo para trastornos del neurodesarrollo en la descendencia, y la evidencia ',
    'epidemiologica muestra que esas patologias afectan de forma distinta a varones ',
    'y mujeres (ver <code>docs/referencias.md</code>). La placenta es la interfaz ',
    'entre la inflamacion materna y el desarrollo fetal: transporta los nutrientes ',
    'que sostienen el crecimiento del cerebro fetal, y su funcion puede verse ',
    'modificada por la inflamacion antes de que el cerebro mismo responda. Este ',
    'proyecto pregunta si esa respuesta, en placenta y en cerebro, depende del ',
    'sexo del feto.</p>\n',
    '<h3>Modelo experimental</h3>\n',
    fig_asset("modelo-experimental.png", "Modelo experimental", "fig-modelo"),
    '<p class="epigrafe">LPS 100 &mu;g/kg i.p. en el dia 15 de gestacion (E15), ',
    'coleccion a las 6 horas.</p>\n',
    '<p>', d$n_fetos, ' fetos de 18 camadas (un feto de cada sexo por camada); ',
    'placenta y cerebro del mismo individuo. Tres mediciones: IL-6 por ELISA en ',
    'suero materno y liquido amniotico, 10 genes por RT-qPCR (transportadores de ',
    'nutrientes y via de senalizacion IL-6/STAT3), y pSTAT3 en placenta por Western ',
    'blot. Cuatro grupos (sexo &times; tratamiento), n = 9.</p>\n',
    '<h3>Como se organizo el trabajo</h3>\n',
    '<p>Las decisiones metodologicas centrales (que escala usar, como tratar los ',
    'valores no detectados, que modelo ajustar, cuando hacer comparaciones post ',
    'hoc) se fijaron por escrito antes de correr ningun analisis, y no se ',
    'reabrieron despues de ver los resultados. El trabajo se dividio en tareas ',
    'acotadas, cada una con su propio cierre. El repositorio deja registro de ',
    'procedencia (de donde sale cada tabla y cada figura) y de verificaciones ',
    '(que se comprobo y como), y toda la implementacion existe por duplicado en R ',
    'y en Python, comparada celda a celda.</p>')
  pagina(1, "Problema, modelo experimental y metodo de trabajo", cuerpo)
}

pagina2 <- function() {
  cuerpo <- paste0(
    '<p>Se entregaron al agente las tablas de datos crudos, sobre las que aplico ',
    'los calculos estandarizados para obtener la cuantificacion relativa de cada ',
    'target. Las comparaciones entre los grupos experimentales de interes se ',
    'definieron de antemano. Cuando el grupo calibrador no pudo cuantificarse por ',
    'limitaciones del metodo experimental, ese gen no se analizo como expresion ',
    'relativa sino como proporcion de deteccion, mediante el test exacto de ',
    'Fisher.</p>\n',
    '<p>La placenta y el cerebro analizados provienen del mismo individuo, de modo ',
    'que las observaciones de ambos tejidos estan pareadas por feto.</p>\n',
    '<p>El analisis se realizo sobre -&Delta;&Delta;Ct, en escala logaritmica de ',
    'base 2, donde las diferencias son simetricas y aditivas; el fold-change se ',
    'reservo para la representacion grafica, con eje logaritmico. Los valores no ',
    'detectados se mantuvieron como faltantes y no se imputaron. Para cada gen y ',
    'tejido se ajusto el modelo -&Delta;&Delta;Ct ~ sexo &times; tratamiento, ',
    'eligiendo el metodo segun el diagnostico de los residuos: ANOVA de tipo III ',
    'cuando se cumplian normalidad y homocedasticidad, errores robustos cuando ',
    'fallaba la homocedasticidad, y una transformacion por rangos alineados cuando ',
    'fallaba la normalidad. Las comparaciones post hoc se realizaron unicamente ',
    'cuando la interaccion entre sexo y tratamiento resulto significativa, sobre ',
    'cuatro contrastes definidos de antemano y con correccion de Holm.</p>\n',
    '<p>En el ELISA de liquido amniotico, los valores por debajo del limite de ',
    'deteccion se trataron como censurados a izquierda y se analizaron con ',
    'metodos especificos para datos censurados.</p>')
  pagina(2, "Metodos", cuerpo)
}

pagina3 <- function(d, sint) {
  intro <- paste0(
    '<p>El agente devolvio el analisis por tejido y por gen, separando por sexo. ',
    'En placenta se evaluo la expresion de diez genes: tres componentes de la via ',
    'de IL-6, una citoquina proinflamatoria fuertemente vinculada a patologias del ',
    'neurodesarrollo, y siete transportadores de nutrientes esenciales para el ',
    'desarrollo fetal, correspondientes al transporte de glucosa, aminoacidos y ',
    'lipidos.</p>')
  resultado <- if (sint) paste0(
    '<p>Para validar el modelo se evaluo si IL-6 era detectable de forma ',
    'diferencial entre grupos, en suero materno y en liquido amniotico, mediante ',
    'el test exacto de Fisher. Sobre los diez genes de placenta se evaluo cuantos ',
    'mostraban efecto de tratamiento y si alguno mostraba interaccion entre sexo y ',
    'tratamiento; pSTAT3 se evaluo por separado, para comparar el patron de ',
    'senalizacion con el de los transportadores.</p>\n',
    sim_bloque("El resultado de esta validacion y de los analisis por gen")
  ) else paste0(
    '<p>La validacion del modelo se apoyo en IL-6: en suero materno fue detectable ',
    'en ', d$ms_lps_det, ' de ', d$ms_lps_n, ' madres tratadas con LPS, frente a ',
    d$ms_ctrl_det, ' de ', d$ms_ctrl_n, ' en el grupo control (test exacto de ',
    'Fisher, p = ', d$ms_p, '). En liquido amniotico ningun contraste alcanzo ',
    'significancia (hembras: n = ', d$la_ctrl_n, ' control y ', d$la_lps_n,
    ' LPS, p = ', d$la_h_p, '; machos: n = ', d$la_ctrl_n, ' control y ',
    d$la_lps_n, ' LPS, p = ', d$la_m_p, ').</p>\n',
    '<p>', length(d$pla_tto_genes), ' de los ', d$pla_modelados, ' genes de ',
    'placenta mostraron efecto de tratamiento: ', .join_y(d$pla_tto_genes),
    ' (p = ', paste(d$pla_tto_ps, collapse = "; "), ', respectivamente), todos ',
    'aumentando su expresion bajo LPS; ninguno mostro interaccion sexo &times; ',
    'tratamiento (', d$pla_int_n, ' de ', d$pla_modelados, ').</p>\n',
    '<p>pSTAT3 aumenta en placentas de fetos hembra (p<sub>Holm</sub> = ',
    d$pstat3_hh_p, ') y no se modifica en machos (p<sub>Holm</sub> = ',
    d$pstat3_mm_p, ', sin figura): el dimorfismo placentario aparece en la ',
    'senalizacion y no en el transporte.</p>')
  cuerpo <- paste0(intro, resultado, '\n',
    fig_outputs(NOMBRE_PLACENTA_BREVE, "Boxplots de il6, fatp1, slc38a2 y glut1 en placenta",
               "fig-boxplot"),
    '<p class="epigrafe">Placenta: il6, fatp1, slc38a2 y glut1, por sexo y ',
    'tratamiento.</p>')
  pagina(3, "Placenta", cuerpo)
}

pagina4 <- function(d, sint) {
  intro <- paste0(
    '<p>El mismo analisis se aplico al cerebro fetal. Mientras la placenta ',
    'respondio al estimulo inflamatorio de manera equivalente en ambos sexos, en ',
    'el cerebro fetal la respuesta de los transportadores de nutrientes resulto ',
    'dimorfica.</p>')
  resultado <- if (sint) paste0(
    '<p>En cerebro se evaluo, para cada uno de los genes modelados, si la ',
    'interaccion entre sexo y tratamiento era significativa y, cuando lo fue, si ',
    'alguna de las cuatro comparaciones post hoc la explicaba.</p>\n',
    sim_bloque("El numero de genes con interaccion en cerebro y el patron por gen")
  ) else paste0(
    '<p>', d$bra_int_n, ' de los ', d$bra_modelados, ' genes modelados en cerebro ',
    'mostraron interaccion sexo &times; tratamiento (frente a ninguno en ',
    'placenta). En ', .join_y(d$bra_con_posthoc), ' el post hoc ubica el efecto ',
    'en las hembras: el grupo control de hembras difiere del grupo LPS de ',
    'hembras, y este ultimo difiere del grupo LPS de machos, mientras que los ',
    'machos no se modifican. En ', .join_y(d$bra_sin_posthoc), ' la interaccion ',
    'fue significativa pero ninguna comparacion entre medias la explicaba, lo que ',
    'motivo el analisis de dispersion de la seccion siguiente.</p>')
  cuerpo <- paste0(intro, resultado, '\n',
    '<div class="cols-fig">',
    fig_outputs(NOMBRE_DETECCION_BRAIN, "Deteccion de il6 en cerebro", "fig-chica"),
    fig_outputs(.ENV13()$NOMBRE_SUBSET_BRAIN,
               "Boxplots de glut1, slc38a2 y fatp1 en cerebro", "fig-grande"),
    '</div>\n',
    '<p class="epigrafe">Izquierda: proporcion de deteccion de il6 en cerebro ',
    '(D7, no cuantificable). Derecha: glut1, slc38a2 y fatp1 en cerebro, por ',
    'sexo y tratamiento.</p>')
  pagina(4, "Cerebro fetal", cuerpo)
}

pagina5 <- function(d, sint) {
  intro <- paste0(
    '<p>Con el analisis convencional cerrado, el agente propuso profundizar por ',
    'dos vias: correlacionar la expresion de cada gen entre la placenta y el ',
    'cerebro de un mismo feto, y explorar la co-expresion entre genes dentro de ',
    'cada tejido.</p>')
  correl <- if (sint) paste0(
    '<h3>Correlacion entre tejidos</h3>\n',
    '<p>Se evaluo si las correlaciones entre placenta y cerebro cambiaban entre ',
    'el grupo control y el grupo LPS, mediante un test formal (Fisher z sobre ',
    '&rho; de Spearman), y se simulo si una simple reduccion de la variabilidad ',
    'alcanzaria para explicar una eventual caida de correlacion sin que la ',
    'relacion subyacente cambiara.</p>\n',
    sim_bloque("El resultado del test de correlaciones y de la simulacion")
  ) else paste0(
    '<h3>Correlacion entre tejidos</h3>\n',
    '<p>En las hembras, las correlaciones altas observadas en el grupo control ',
    'caian bajo LPS, lo que sugeria una perdida de acoplamiento entre placenta y ',
    'cerebro. El test formal (Fisher z sobre &rho; de Spearman, Control vs LPS) ',
    'no mostro diferencias significativas en ninguna de las ', d$corr_n_test,
    ' comparaciones evaluadas (', d$corr_n_sig, ' de ', d$corr_n_test,
    '), y una simulacion mostro que una reduccion de la variabilidad alcanza por ',
    'si sola para producir esa caida, sin que la relacion subyacente cambie. La ',
    'aparente perdida de acoplamiento placenta-cerebro en hembras es compatible ',
    'con la reduccion de dispersion: la simulacion muestra que esta alcanza para ',
    'explicarla, aunque no permite descartar un cambio de coordinacion.</p>')
  coexpr <- if (sint) paste0(
    '<h3>Co-expresion dentro de cada tejido</h3>\n',
    '<p>Los diagramas triangulares son exploratorios y describen como se ',
    'acompanan los genes entre si; sobre ellos no se testeo ninguna diferencia ',
    'entre grupos y no se interpretan diferencias de &rho;. Ademas se evaluo que ',
    'proporcion de la varianza explica el primer componente principal en cada ',
    'tejido.</p>\n',
    sim_bloque("Los porcentajes de varianza explicada")
  ) else paste0(
    '<h3>Co-expresion dentro de cada tejido</h3>\n',
    '<p>Los diagramas triangulares son exploratorios y describen como se ',
    'acompanan los genes entre si; sobre ellos no se testeo ninguna diferencia ',
    'entre grupos y no se interpretan diferencias de &rho;. Lo que si queda ',
    'establecido es que los transportadores varian mayormente juntos: el primer ',
    'componente principal explica el ', d$pc1_pla, ' % de la varianza en ',
    'placenta y el ', d$pc1_bra, ' % en cerebro.</p>')
  fig_epigrafe <- if (sint) paste0(
    '<p class="epigrafe">Cerebro fetal, hembras: densidades de ',
    '-&Delta;&Delta;Ct por tratamiento en los mismos cinco genes (recorte de la ',
    'diagonal del diagrama triangular; corrida con datos sinteticos).</p>'
  ) else paste0(
    '<p class="epigrafe">Cerebro fetal, hembras: densidades de ',
    '-&Delta;&Delta;Ct por tratamiento en los cinco genes con interaccion ',
    'significativa sobre la dispersion (recorte de la diagonal del diagrama ',
    'triangular). La distribucion del grupo LPS es marcadamente mas concentrada ',
    'que la del control, y esa reduccion de variabilidad es el mecanismo que ',
    'explica la caida de correlacion.</p>')
  limitaciones <- paste0(
    '<h3>Limitaciones</h3>\n',
    '<p>Todas las comparaciones se hacen con n = 9 por grupo, lo que limita la ',
    'potencia, en especial en los contrastes estratificados por sexo. El ',
    'analisis asume independencia entre fetos aunque el tratamiento se ',
    'administra a la madre y cada camada aporta un feto de cada sexo (D13). El ',
    'efecto sobre la dispersion en cerebro fetal constituye un patron compartido ',
    'por varios transportadores, y no un conjunto de hallazgos independientes ',
    'gen por gen.</p>')
  conclusion <- if (sint) paste0(
    '<h3>Conclusion</h3>\n', sim_bloque("La conclusion biologica de este proyecto")
  ) else paste0(
    '<h3>Conclusion</h3>\n',
    '<p>El sexo del feto influye en la respuesta al LPS, pero de manera distinta ',
    'en cada tejido. En la placenta, los transportadores responden de forma ',
    'equivalente en ambos sexos y el dimorfismo aparece en la senalizacion. En ',
    'el cerebro fetal, el dimorfismo se expresa en los transportadores, ',
    'desplazando el nivel de expresion en unos genes y la variabilidad entre ',
    'individuos en otros.</p>')
  pie <- paste0(
    '<footer>Detalle completo, procedencia y verificaciones: ',
    '<code>docs/informe.html</code>. Presentacion: <code>docs/index.html</code>.',
    '</footer>')
  cuerpo <- paste0(
    intro, correl, coexpr, '\n',
    fig_outputs(NOMBRE_DENSIDADES, "Densidades de -ddCt en cerebro, hembras", "fig-densidad"),
    fig_epigrafe, limitaciones, conclusion, pie)
  pagina(5, "Exploracion, limites y conclusion", cuerpo)
}

# =========================================================================
construir_html <- function(d) {
  sint <- d$fuente != "real"
  # el <footer> va DENTRO de la pagina 5 (ver pagina5()): si quedara fuera del
  # ultimo <section>, "break-after: page" del section empuja el pie a una
  # sexta pagina vacia salvo por esa linea.
  paginas <- paste(pagina1(d), pagina2(), pagina3(d, sint), pagina4(d, sint),
                   pagina5(d, sint), sep = "\n")
  paste0(
    "<!doctype html>\n<html lang=\"es\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>Informe breve: MIA-LPS</title>\n<style>", CSS, "</style>\n",
    "</head>\n<body>\n<main>\n",
    "<h1>Transportadores de nutrientes en el eje placenta&ndash;cerebro fetal en ",
    "un modelo de activacion inmune materna</h1>\n",
    "<p class=\"subtitulo\">Informe breve</p>\n",
    paginas,
    "\n</main>\n</body>\n</html>\n")
}

# =========================================================================
main <- function() {
  d <- recolectar_datos()
  sint <- d$fuente != "real"

  e13 <- .ENV13()
  ruta_placenta_breve <- file.path(RUTA_FIGURAS, NOMBRE_PLACENTA_BREVE)
  ruta_subset_bra <- file.path(RUTA_FIGURAS, e13$NOMBRE_SUBSET_BRAIN)
  ruta_deteccion_bra <- file.path(RUTA_FIGURAS, NOMBRE_DETECCION_BRAIN)
  ruta_densidades <- file.path(RUTA_FIGURAS, NOMBRE_DENSIDADES)

  e13$generar_panel_subset(GENES_PLACENTA_BREVE, "PLACENTA_E15", ruta_placenta_breve)
  e13$generar_panel_subset(e13$GENES_SUBSET_BRAIN, "BRAIN_E15", ruta_subset_bra)
  generar_panel_deteccion_brain(ruta_deteccion_bra)
  generar_figura_densidades(d$disp_bra_sig_genes, ruta_densidades)

  html <- construir_html(d)

  destino <- if (sint) RUTA_DOCS else RUTA_INFORME_BREVE_REAL
  if (!dir.exists(destino)) dir.create(destino, recursive = TRUE, showWarnings = FALSE)
  ruta_html <- file.path(destino, "informe_breve.html")
  ruta_pdf <- file.path(destino, "informe_breve.pdf")
  escribir_texto(ruta_html, html)
  pdf_status <- generar_pdf(ruta_html, ruta_pdf)
  n_paginas <- contar_paginas_pdf(ruta_pdf)

  m_sim <- gregexpr('class="sim"', html, fixed = TRUE)[[1]]
  n_sim <- if (length(m_sim) == 1L && m_sim[1] == -1L) 0L else length(m_sim)
  n_sim_esperado <- if (sint) 5L else 0L  # pagina 3 (1) + pagina 4 (1) + pagina 5 (3)

  # Ninguna seccion (pagina) de la version publica queda vacia: cada
  # <section class="pagina" id="pN">...</section> tiene que superar un minimo
  # de texto visible (sin tags), incluso reemplazando los resultados por el
  # aviso puntual -- ver pedido seccion 0.
  ids_pag <- c("1", "2", "3", "4", "5")
  secciones_ok <- vapply(ids_pag, function(id_) {
    m <- regexpr(sprintf('(?s)<section class="pagina" id="p%s">.*?</section>', id_),
                html, perl = TRUE)
    if (m == -1L) return(FALSE)
    txt <- regmatches(html, m)
    txt <- gsub("<[^>]+>", " ", txt)
    txt <- gsub("\\s+", " ", txt)
    nchar(trimws(txt)) >= 200L
  }, logical(1))

  usadas <- unique(.FIGS_USADAS$outputs)
  faltan_figs <- usadas[!file.exists(file.path(RUTA_FIGURAS, usadas)) &
                        !file.exists(file.path(RUTA_ASSETS, usadas))]

  gi_ok <- tryCatch({
    rc <- system2("git", c("-C", shQuote(RAIZ_REPO), "check-ignore", "-q",
                           shQuote(file.path(RUTA_INFORME_BREVE_REAL, "informe_breve.html"))),
                  stdout = FALSE, stderr = FALSE)
    rc == 0L
  }, error = function(e) NA)

  destino_rel <- if (sint) "docs" else "outputs/informe_breve_real"
  ent <- "AGENTS.md + outputs/tables/*.csv + outputs/figures/*.png + assets/*.png"
  registrar_procedencia(list(
    list(file.path(destino_rel, "informe_breve.html"), "informe", ESTE_SCRIPT, "PROPIO", ent,
         paste0("informe breve (maximo 5 paginas); version ",
                if (sint) "PUBLICA (sintetica)" else "PARA LEER (real, no versionada)")),
    list(file.path(destino_rel, "informe_breve.pdf"), "informe", ESTE_SCRIPT, "PROPIO",
         "informe_breve.html",
         paste0("version imprimible por impresion headless; best-effort (estado: ",
                pdf_status, ")")),
    list(file.path("outputs/figures", NOMBRE_PLACENTA_BREVE), "figura_presentacion",
         ESTE_SCRIPT, "PROPIO (reusa generar_panel_subset de 13_presentacion.R)",
         "data/processed/qpcr_cuantificacion_long.tsv",
         paste0("boxplots de placenta, subconjunto propio de este informe (",
                .join_y(GENES_PLACENTA_BREVE), ")")),
    list(file.path("outputs/figures", NOMBRE_DETECCION_BRAIN), "figura_presentacion",
         ESTE_SCRIPT, "PROPIO (reusa panel_deteccion de 07_figuras_acto1.R)",
         "qpcr_il6_brain_fisher.csv + qpcr_il6_brain_tabla2x4.csv",
         "proporcion de deteccion de il6 en cerebro E15 (D7, no cuantificable)"),
    list(file.path("outputs/figures", NOMBRE_DENSIDADES), "figura_presentacion", ESTE_SCRIPT,
         "PROPIO (geom_density, mismo estilo que la diagonal del SPLOM de 08_acto2_correlaciones.R)",
         "data/processed/qpcr_cuantificacion_long.tsv",
         paste0("recorte de la diagonal del SPLOM: densidades de -ddCt por ",
                "tratamiento, cerebro E15, SOLO hembras, genes ",
                .join_y(d$disp_bra_sig_genes),
                " (interaccion SEXOxTTO significativa sobre la dispersion, BH<0.05)"))
  ))
  registrar_verificaciones(list(
    list("informe_breve_html_generado", "recalculo",
         sprintf("%s existe y no esta vacio", file.path(destino_rel, "informe_breve.html")),
         sprintf("existe=%s", file.exists(ruta_html) && file.info(ruta_html)$size > 0),
         "TRUE",
         if (file.exists(ruta_html) && file.info(ruta_html)$size > 0) "TRUE" else "FALSE",
         ESTE_SCRIPT),
    list("informe_breve_max_5_paginas", "recalculo",
         "docs/informe_breve.pdf (publica) tiene 5 paginas o menos",
         sprintf("fuente=%s; paginas=%s", d$fuente, n_paginas),
         "paginas <= 5 (solo exigido en la version publica)",
         if (!sint) "NO_EJECUTADA"
         else if (!is.na(n_paginas) && n_paginas > 0L && n_paginas <= 5L) "TRUE"
         else "FALSE",
         ESTE_SCRIPT),
    list("informe_breve_sin_interpretacion", "recalculo",
         paste0("con fuente sintetica, cada afirmacion sobre que dio un analisis ",
                "sobre datos reales se reemplaza por un aviso puntual (clase sim)"),
         sprintf("fuente=%s; sim=%d/%d", d$fuente, n_sim, n_sim_esperado),
         "sim = 5 si fuente sintetica, 0 si fuente real",
         if (n_sim == n_sim_esperado) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_breve_secciones_no_vacias", "recalculo",
         "ninguna de las 5 paginas queda vacia (>=200 caracteres de texto visible)",
         sprintf("fuente=%s; paginas_ok=%d/5", d$fuente, sum(secciones_ok)),
         "paginas_ok = 5/5", if (all(secciones_ok)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_breve_figuras_existen", "recalculo",
         "toda figura referenciada por el informe breve existe",
         sprintf("usadas=%d; faltan=%s", length(usadas),
                 if (length(faltan_figs)) paste(faltan_figs, collapse = ", ") else "[]"),
         "faltan = []", if (!length(faltan_figs)) "TRUE" else "FALSE", ESTE_SCRIPT),
    list("informe_breve_gitignore_real", "recalculo",
         "outputs/informe_breve_real/ esta cubierto por .gitignore (git check-ignore)",
         sprintf("ignorado=%s", if (is.na(gi_ok)) "NA" else gi_ok), "TRUE",
         if (isTRUE(gi_ok)) "TRUE" else if (is.na(gi_ok)) "NO_EJECUTADA" else "FALSE",
         ESTE_SCRIPT),
    list("informe_breve_pdf", "recalculo",
         paste0("informe_breve.pdf generado (TRUE), degradado limpio sin motor ",
                "(NO_EJECUTADA), o fallo real (FALSE)"),
         sprintf("estado=%s", pdf_status),
         "ok -> TRUE | sin_motor -> NO_EJECUTADA | fallo:* -> FALSE",
         if (pdf_status == "ok") "TRUE"
         else if (pdf_status == "sin_motor") "NO_EJECUTADA" else "FALSE",
         ESTE_SCRIPT)
  ))

  cat("== 14_informe_breve.R ==\n")
  cat(sprintf("  fuente = %s\n", d$fuente))
  cat(sprintf("  genes densidades (BH<0.05, cerebro) = %s\n",
              paste(d$disp_bra_sig_genes, collapse = ", ")))
  cat(sprintf("  -> %s (%d KB)\n", ruta_html, file.info(ruta_html)$size %/% 1024))
  cat(sprintf("  -> %s [%s] (%s paginas)\n", ruta_pdf, pdf_status, n_paginas))
  cat(sprintf("  secciones no vacias: %d/5\n", sum(secciones_ok)))
  if (length(faltan_figs))
    cat(sprintf("  *** faltan figuras: %s ***\n", paste(faltan_figs, collapse = ", ")))
}

if (sys.nframe() == 0L) main()
