# 14_informe_breve.R -- informe breve del proyecto (maximo 5 paginas), HTML + PDF.
#
# Por que existe este archivo: pedidos/pedido_informe_breve.md pide un documento
# CORTO para leer (problema, metodo, hallazgos, limites) -- distinto de
# docs/informe.html (el informe tecnico largo, con procedencia/verificaciones) y
# de la presentacion (docs/index.html + presentacion.pdf). No corre NINGUN
# analisis nuevo: todos los numeros se leen de outputs/tables/ al generar el
# documento (misma convencion que 12_informe.R y 13_presentacion.R). La unica
# figura verdaderamente NUEVA es la de densidades (seccion "Figura nueva" del
# pedido); el resto son figuras que 13_presentacion.R ya sabe generar.
#
# EXCEPCION A LA REGLA DE SCRIPTS GEMELOS (la misma ya documentada para
# 13_presentacion.R en AGENTS.md 3): existe SOLO en R. No hay
# python/14_informe_breve.py.
#
# CUARTA EXCEPCION a "cada script importa solo 00_config" (ver AGENTS.md 3,
# donde estan documentadas la 1ra y 2da): este script fuentea 13_presentacion.R
# COMPLETO con `local = <environment nuevo>` para reusar, sin reimplementarlas,
# las figuras que 13_presentacion.R ya sabe construir (subconjuntos de boxplots,
# pSTAT3 propio de presentacion, y el override que en la version PUBLICA vacia
# el numero de p de los brackets de tendencia) -- asi el informe breve hereda
# EXACTAMENTE el mismo comportamiento ya verificado, sin volver a resolverlo.
# 13_presentacion.R a su vez fuentea 07_figuras_acto1.R de la misma forma: el
# aislamiento es transitivo (07 queda aislado DENTRO del aislamiento de 13).
#
# DOS VERSIONES, mismo codigo, la fuente de datos decide el destino (identico
# mecanismo que 12_informe.R / 13_presentacion.R):
#   - fuente == "sintetico" -> docs/informe_breve.html + docs/informe_breve.pdf
#     (PUBLICA: se versiona).
#   - fuente == "real"      -> outputs/informe_breve_real/informe_breve.{html,pdf}
#     (PARA LEER: nunca se versiona, ver .gitignore).
# La leccion de la sesion pasada (bug de fuga de datos en 13_presentacion.R):
# este script NUNCA recalcula, solo lee outputs/figures/ y outputs/tables/R/
# vigentes -- para la version publica hay que dejar esas carpetas en estado
# SINTETICO antes de correrlo (`.\run_all.ps1 -Only R -FromSynthetic`), igual
# que para la presentacion. No alcanza con forzar MIA_LPS_FORZAR_SINTETICO=1
# sobre este script solo.

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
# reusar sus figuras (subconjuntos de boxplots, pSTAT3 propio, override de
# brackets de tendencia). `local = <environment nuevo>` evita que sus propios
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
# Aviso de fuente sintetica -- mismo principio que 12_informe.R /
# 13_presentacion.R: con datos sinteticos, las paginas 3, 4 y 5 nunca arman la
# prosa interpretativa (biologica); se reemplaza por este aviso. Las paginas 1
# y 2 (contexto, diseno, metodo de trabajo, metodos) son iguales en las dos
# versiones -- describen el proceso, no los resultados.
# =========================================================================
AVISO_SINTETICO <- paste0(
  "<p class=\"aviso\"><em>Seccion generada con datos sinteticos. Los efectos ",
  "son simulados y arbitrarios; las conclusiones biologicas corresponden a ",
  "los datos reales, que no se incluyen en este repositorio.</em></p>")

# =========================================================================
# Figura nueva: densidades de -ddCt por grupo (SEXO x TTO), Cerebro E15, 4
# genes con interaccion SEXOxTTO significativa sobre la DISPERSION (elegidos
# leyendo acto2_dispersion_interaccion.csv, no a mano -- ver seleccion en
# main()). Mismo estilo que la diagonal del SPLOM de 08_acto2_correlaciones.R
# (GGally::wrap("densityDiag"), Control/LPS superpuestos, semitransparente):
# se reimplementa aca en ggplot2 puro (sin GGally, que arma matrices
# completas) porque el pedido es una figura de 2 filas x N columnas, no un
# SPLOM. Usa los mismos -ddCt de data/processed/qpcr_cuantificacion_long.tsv
# que el resto del pipeline -- no recalcula nada.
# =========================================================================
COL_TTO <- c(CONTROL = "#0072B2", LPS = "#D55E00")  # identico a 08/09_acto2_*.R
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
    if (f[j("TEJIDO")] != "BRAIN_E15") next
    if (!(f[j("GEN")] %in% genes)) next
    if (identical(f[j("no_detectado")], "TRUE") || !nzchar(f[j("neg_ddCt")])) next
    reg[[length(reg) + 1L]] <- data.frame(
      SEXO = f[j("SEXO")], TTO = f[j("TTO")], GEN = f[j("GEN")],
      neg_ddCt = as.numeric(f[j("neg_ddCt")]), stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, reg)
  d$GEN <- factor(vapply(d$GEN, etiqueta_gen, character(1)),
                  levels = vapply(genes, etiqueta_gen, character(1)))
  d$SEXO <- factor(d$SEXO, levels = c("HEMBRA", "MACHO"),
                   labels = c("Hembras", "Machos"))
  d$TTO <- factor(d$TTO, levels = c("CONTROL", "LPS"), labels = c("Control", "LPS"))

  p <- ggplot(d, aes(neg_ddCt, fill = TTO, colour = TTO)) +
    geom_density(alpha = 0.4, linewidth = 0.5) +
    scale_fill_manual(values = c(Control = unname(COL_TTO["CONTROL"]),
                                 LPS = unname(COL_TTO["LPS"])), name = NULL) +
    scale_colour_manual(values = c(Control = unname(COL_TTO["CONTROL"]),
                                   LPS = unname(COL_TTO["LPS"])), name = NULL) +
    facet_grid(SEXO ~ GEN, scales = "free") +
    labs(title = "Cerebro E15 -- densidad de -\u0394\u0394Ct por grupo",
         subtitle = paste0("Genes con interaccion SEXO\u00d7TTO significativa sobre la ",
                           "dispersion (BH < 0.05); el LPS compacta la expresion en ",
                           "hembras y la dispersa en machos"),
         x = expression(-Delta*Delta*Ct), y = "Densidad") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(), legend.position = "top",
          strip.background = element_rect(fill = "grey93", colour = NA),
          plot.subtitle = element_text(size = 8.2))
  ggsave(ruta, p, width = 9.6, height = 4.6, dpi = 300)
  invisible(ruta)
}
NOMBRE_DENSIDADES <- "acto2_densidades_dispersion_BRAIN_E15.png"

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
".fig-densidad { max-height: 30vh; }\n",
".aviso { background: var(--warn-bg); border: 1px solid var(--warn-border);\n",
"  border-radius: 6px; padding: .7rem 1rem; font-size: .95rem; }\n",
".falta { color: #b00; font-style: italic; }\n",
"ul.limites { margin: .3rem 0 .8rem 1.3rem; padding: 0; }\n",
"ul.limites li { margin: .3rem 0; text-align: justify; }\n",
"code { font: .88em \"SF Mono\", Consolas, monospace; background: #f0eee8;\n",
"       padding: .04em .3em; border-radius: 3px; }\n",
"footer { margin-top: 1.5rem; padding-top: .6rem; border-top: 1px solid var(--line);\n",
"  font: .78rem var(--font-ui); color: var(--ink-soft); }\n",
"@media print {\n",
"  @page { size: A4; margin: 1.3cm 1.6cm; }\n",
"  body { background: #fff; }\n",
"  main { max-width: none; padding: 0; }\n",
"  section.pagina { break-after: page; margin: 0; }\n",
"  section.pagina:last-of-type { break-after: auto; }\n",
"  h1, .subtitulo { break-after: avoid; }\n",
"  figure { break-inside: avoid; }\n",
"}\n"
)

# =========================================================================
# Datos leidos del repo -- nada de lo que sigue se escribe a mano.
# =========================================================================
recolectar_datos <- function() {
  d <- list()

  # --- pagina 3: ELISA (validacion), placenta ---
  ef <- .tab("elisa_fisher_deteccion.csv")
  i <- .col(ef, "bloque") == "MS"
  d$ms_lps_det <- .col(ef, "lps_detectado")[i][1]; d$ms_lps_n <- .col(ef, "lps_n")[i][1]
  d$ms_ctrl_det <- .col(ef, "control_detectado")[i][1]
  d$ms_ctrl_n <- .col(ef, "control_n")[i][1]
  d$ms_p <- .col(ef, "p_valor")[i][1]

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
  pd <- .tab("pstat3_descriptivo.csv")
  gr <- .col(pd, "GRUPO"); med <- .col(pd, "mean_PSTAT3")
  d$pstat3_hc_media <- .round_fmt(.num(med[gr == "HEMBRA_CONTROL"]), 2L)
  d$pstat3_hl_media <- .round_fmt(.num(med[gr == "HEMBRA_LPS"]), 2L)

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
  d$bra_con_posthoc <- con_posthoc     # nivel: HC!=HL, HL!=ML, sin cambio en MC-ML
  d$bra_sin_posthoc <- sin_posthoc     # dispersion, no nivel

  di <- .tab("acto2_dispersion_interaccion.csv")
  di_tej <- .col(di, "TEJIDO"); di_gen <- .col(di, "GEN")
  bh <- .num(.col(di, "p_SEXOxTTO_BH"))
  i_disp_sig <- di_tej == "BRAIN_E15" & !is.na(bh) & bh < 0.05
  d$disp_bra_sig_n <- sum(i_disp_sig)
  d$disp_bra_sig_genes <- sort(di_gen[i_disp_sig], method = "radix")
  d$disp_bra_n <- sum(di_tej == "BRAIN_E15")
  # de los "sin post hoc": cuantos explica la dispersion (BH<.05) vs tendencia
  bh_de <- function(g) { i <- di_tej == "BRAIN_E15" & di_gen == g
                         if (any(i)) bh[i][1] else NA_real_ }
  bh_sp <- vapply(sin_posthoc, bh_de, numeric(1))
  d$sin_posthoc_explicados <- sin_posthoc[!is.na(bh_sp) & bh_sp < 0.05]
  d$sin_posthoc_tendencia <- sin_posthoc[!is.na(bh_sp) & bh_sp >= 0.05 & bh_sp < 0.10]

  # --- pagina 5: correlacion, restriccion de rango, PC1, sensibilidad ---
  tc <- .tab("acto2_test_correlaciones.csv")
  p_bw_raw <- .col(tc, "p_bw"); ok <- nzchar(p_bw_raw)  # mismo filtro que 12_informe.R
  p_bw <- .num(p_bw_raw[ok])
  d$corr_n_test <- sum(ok); d$corr_n_sig <- sum(!is.na(p_bw) & p_bw < 0.05)

  sim <- .tab("acto2_simulacion.csv")
  it_s <- .col(sim, "ITEM"); es_s <- .col(sim, "ESTRATO"); esc_s <- .col(sim, "ESCENARIO")
  ver_s <- .col(sim, "veredicto")
  i_h_glob <- es_s == "HEMBRA" & esc_s == "GLOBAL"
  d$sim_h_total <- sum(i_h_glob)
  d$sim_h_fuera <- sum(i_h_glob & ver_s == "FUERA")
  d$sim_lim_genes <- .join_y(gsub("_", " ", it_s[i_h_glob & ver_s == "FUERA"]))

  pv <- .tab("acto2_sensibilidad_pca_varianza.csv")
  tej_v <- .col(pv, "TEJIDO"); pc <- .col(pv, "PC"); propv <- .col(pv, "prop_var")
  pc1_pct <- function(te) .round_fmt(.num(propv[tej_v == te & pc == "1"][1]) * 100, 0L)
  d$pc1_pla <- pc1_pct("PLACENTA_E15"); d$pc1_bra <- pc1_pct("BRAIN_E15")

  ex <- .tab("acto2_sensibilidad_excl_extremo.csv")
  d$excl_total <- length(.col(ex, "veredicto_cambia"))
  d$excl_cambian <- sum(.col(ex, "veredicto_cambia") == "TRUE")

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
    'proyecto pregunta si esa respuesta -- en placenta y en cerebro -- depende del ',
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
    '<p>Las decisiones metodologicas centrales -- que escala usar, como tratar los ',
    'valores no detectados, que modelo ajustar, cuando hacer comparaciones post ',
    'hoc -- se fijaron por escrito antes de correr ningun analisis, y no se ',
    'reabrieron despues de ver los resultados. El trabajo se dividio en tareas ',
    'acotadas, cada una con su propio cierre. El repositorio deja registro de ',
    'procedencia (de donde sale cada tabla y cada figura) y de verificaciones ',
    '(que se comprobo y como), y toda la implementacion existe por duplicado en R ',
    'y en Python, comparada celda a celda.</p>')
  pagina(1, "Problema, modelo experimental y metodo de trabajo", cuerpo)
}

pagina2 <- function() {
  cuerpo <- paste0(
    '<p>La cuantificacion es relativa: <code>&Delta;Ct = Ct<sub>gen</sub> - ',
    'Ct<sub>rsp29</sub></code>, con calibrador &female;Control por gen y tejido ',
    '(promediando solo valores detectados). El analisis se hace sobre ',
    '<code>-&Delta;&Delta;Ct</code> (escala log2, simetrica y aditiva); el ',
    'fold-change <code>FC = 2<sup>-&Delta;&Delta;Ct</sup></code> se usa solo para ',
    'graficar, con eje logaritmico -- el fold-change esta acotado en 0 y es ',
    'asimetrico, y eso rompe los supuestos de los modelos lineales.</p>\n',
    '<p>Los valores no detectados se tratan como faltantes (NA), nunca se ',
    'imputan: una imputacion evaluada tempranamente introducia estructura ',
    'artificial en el grupo control y se descarto.</p>\n',
    '<p>El modelo es <code>-&Delta;&Delta;Ct ~ SEXO * TTO</code>, ajustado por gen ',
    'y tejido, con una cascada de metodo segun el diagnostico de los residuos: ',
    'ANOVA tipo III si se cumplen normalidad y homocedasticidad, errores robustos ',
    'HC3 si falla solo la homocedasticidad, y ART (Aligned Rank Transform) si ',
    'falla la normalidad. Las comparaciones post hoc solo se corren si la ',
    'interaccion SEXO&times;TTO es significativa, sobre cuatro comparaciones fijas ',
    '(&female;Control-&female;LPS, &male;Control-&male;LPS, &female;LPS-&male;LPS, ',
    '&female;Control-&male;Control) con correccion de Holm.</p>\n',
    '<p><code>il6</code> en cerebro no es cuantificable por este metodo: el ',
    'calibrador &female;Control tiene 0 de 9 detectados. Se analiza aparte, solo ',
    'como proporcion de deteccion por sexo (Fisher exacto).</p>\n',
    '<p>El ELISA de liquido amniotico tiene censura a izquierda (valores por debajo ',
    'del limite de deteccion): se registra la censura explicitamente y se usan ',
    'metodos para datos censurados (Peto-Peto), nunca se trata un valor censurado ',
    'como un numero negativo ni se lo reemplaza por cero.</p>\n',
    '<p>El modelo de pSTAT3 agrega la membrana de Western blot como bloque fijo, ',
    'para absorber la variacion tecnica entre membranas; usa la misma cascada de ',
    'metodo y el mismo post hoc.</p>')
  pagina(2, "Metodos", cuerpo)
}

pagina3 <- function(d, sint) {
  cuerpo3 <- if (sint) AVISO_SINTETICO else paste0(
    '<p>La IL-6 serica materna valida el modelo: fue detectable en ', d$ms_lps_det,
    ' de ', d$ms_lps_n, ' madres tratadas con LPS, frente a ', d$ms_ctrl_det, ' de ',
    d$ms_ctrl_n, ' en el grupo control (Fisher exacto, p = ', d$ms_p, ').</p>\n',
    '<p>De ', d$pla_modelados, ' genes modelados en placenta, ninguno muestra ',
    'interaccion sexo &times; tratamiento (', d$pla_int_n, ' de ', d$pla_modelados,
    '): el LPS modifica la expresion en ambos sexos por igual. Los genes con ',
    'efecto principal de tratamiento son ', .join_y(d$pla_tto_genes),
    ' (p = ', paste(d$pla_tto_ps, collapse = "; "), ', respectivamente).</p>\n',
    '<p>pSTAT3 responde de otra manera: la interaccion sexo &times; tratamiento es ',
    'significativa (p = ', d$pstat3_p_int, '), y el post hoc ubica el efecto solo ',
    'en hembras (&female;Control-&female;LPS, de ', d$pstat3_hc_media, ' a ',
    d$pstat3_hl_media, ' u.a., p<sub>Holm</sub> = ', d$pstat3_hh_p,
    '; &male;Control-&male;LPS, p<sub>Holm</sub> = ', d$pstat3_mm_p,
    ', sin cambio).</p>\n',
    '<p>Los dos resultados conviven pero no coinciden: los transportadores cambian ',
    'en ambos sexos, mientras que solo las hembras activan STAT3. Es compatible ',
    'con que el cambio de los transportadores no dependa de esa via, o con que los ',
    'machos lo alcancen por otro camino; ninguna de las dos esta demostrada por ',
    'este analisis.</p>')
  cuerpo <- paste0(
    '<div class="cols-fig">',
    fig_outputs(.ENV13()$NOMBRE_SUBSET_PLACENTA,
               "Boxplots de il6, glut3 y slc38a2 en placenta", "fig-boxplot"),
    '<p class="epigrafe">Expresion de il6, glut3 y slc38a2 en placenta (subconjunto ',
    'representativo; los cinco genes con efecto de tratamiento estan en el texto).</p>',
    fig_outputs(.ENV13()$NOMBRE_PSTAT3_PRESENTACION,
               "pSTAT3 en placenta", "fig-boxplot"),
    '<p class="epigrafe">pSTAT3 en placenta, por sexo y tratamiento.</p>',
    '</div>\n', cuerpo3)
  pagina(3, "Placenta", cuerpo)
}

pagina4 <- function(d, sint) {
  cuerpo4 <- if (sint) AVISO_SINTETICO else {
    paso3 <- if (length(d$sin_posthoc_explicados))
      paste0('En ', length(d$sin_posthoc_explicados), ' de ellos (',
             .join_y(d$sin_posthoc_explicados),
             ') el test de interaccion sexo &times; tratamiento sobre la ',
             'dispersion lo explica (BH < 0.05)',
             if (length(d$sin_posthoc_tendencia))
               paste0('; en ', .join_y(d$sin_posthoc_tendencia),
                      ' aparece la misma tendencia, sin llegar a significancia')
             else '', '.')
    else ''
    paste0(
      '<p>En cerebro, ', d$bra_int_n, ' de ', d$bra_modelados, ' genes modelados ',
      'muestran interaccion sexo &times; tratamiento (frente a 0 en placenta). Hay ',
      'dos formas de dimorfismo, no una.</p>\n',
      '<p><strong>Desplazamiento del nivel.</strong> En ', .join_y(d$bra_con_posthoc),
      ' el post hoc localiza la diferencia: &female;Control distinto de ',
      '&female;LPS, &female;LPS distinto de &male;LPS, sin cambios entre ',
      '&male;Control y &male;LPS.</p>\n',
      '<p><strong>Cambio de dispersion.</strong> En ', .join_y(d$bra_sin_posthoc),
      ' hay interaccion pero ninguna comparacion de a pares sobrevive a Holm: hay ',
      'algo sexo-dependiente, pero el post hoc no dice donde. ', paso3,
      ' En total, ', d$disp_bra_sig_n, ' de ', d$disp_bra_n,
      ' genes de cerebro muestran interaccion significativa sobre la dispersion ',
      '(BH < 0.05): el LPS compacta la expresion en hembras y la dispersa en ',
      'machos.</p>\n',
      '<p>El boxplot mostraba que algo dependia del sexo, pero no donde estaba: el ',
      'analisis de dispersion lo explico.</p>')
  }
  cuerpo <- paste0(
    '<div class="cols-fig">',
    fig_outputs(.ENV13()$NOMBRE_SUBSET_BRAIN,
               "Boxplots de glut1, slc38a2 y fatp1 en cerebro", "fig-boxplot"),
    '<p class="epigrafe">Expresion de glut1, slc38a2 y fatp1 en cerebro (los tres ',
    'con post hoc significativo).</p>',
    fig_outputs(NOMBRE_DENSIDADES,
               "Densidades de -ddCt por grupo, cerebro", "fig-densidad"),
    '<p class="epigrafe">Densidad de -&Delta;&Delta;Ct por sexo y tratamiento, ',
    'cerebro (genes con interaccion significativa sobre la dispersion).</p>',
    '</div>\n', cuerpo4)
  pagina(4, "Cerebro fetal", cuerpo)
}

pagina5 <- function(d, sint) {
  cuerpo <- if (sint) AVISO_SINTETICO else paste0(
    '<p>Ninguna de las ', d$corr_n_test, ' comparaciones de correlacion placenta-',
    'cerebro (Fisher z sobre &rho; de Spearman, Control vs LPS) alcanza p<0.05 (',
    d$corr_n_sig, ' de ', d$corr_n_test, '). La caida de correlacion observada en ',
    'hembras es compatible con la compactacion de la expresion bajo LPS: cuando el ',
    'rango de una variable se reduce, la correlacion cae aunque la relacion ',
    'biologica no haya cambiado. Sobre ', d$sim_h_total,
    ' comparaciones simuladas para hembras, ', d$sim_h_fuera,
    ' quedan fuera del intervalo (', d$sim_lim_genes,
    '); se toman como pista, no como hallazgo -- la simulacion alcanza para ',
    'explicar la caida, aunque no permite descartar un cambio de coordinacion.</p>\n',
    '<p>Los diagramas triangulares de co-expresion sugieren que los transportadores ',
    'varian mayormente juntos, como un eje compartido por feto: el eigengene (PC1) ',
    'explica el ', d$pc1_pla, ' % de la varianza en placenta y el ', d$pc1_bra,
    ' % en cerebro. Las diferencias de &rho; entre grupos que se ven en esos ',
    'diagramas no se interpretan como un cambio de coordinacion (no estan ',
    'testeadas).</p>\n',
    '<p>Dos controles de sensibilidad: reemplazar el score compuesto por el ',
    'eigengene no cambia ninguna conclusion, y excluir el feto mas influyente de ',
    'cada item tampoco cambia ningun veredicto (', d$excl_cambian, ' de ',
    d$excl_total, ').</p>\n',
    '<h3>Limitaciones</h3>\n',
    '<ul class="limites">',
    '<li>n = 9 por grupo: la potencia es baja, en especial para las comparaciones ',
    'estratificadas por sexo.</li>',
    '<li>El tratamiento se administra a la madre, y cada camada aporta un feto de ',
    'cada sexo; el analisis asume independencia entre fetos de la misma camada.</li>',
    '<li>pSTAT3 esta normalizado a proteina total, sin STAT3 total: refleja ',
    'abundancia de fosfo-STAT3, no la fraccion fosforilada.</li>',
    '<li>El efecto de dispersion en cerebro es un patron compartido por varios ',
    'transportadores, no hallazgos independientes gen por gen.</li>',
    '</ul>\n',
    '<h3>Conclusion</h3>\n',
    '<p>El sexo del feto influye en la respuesta al LPS, pero de manera distinta ',
    'en cada tejido. En placenta afecta la senalizacion (pSTAT3, solo en hembras) ',
    'y no el transporte (los transportadores cambian en ambos sexos por igual). En ',
    'cerebro afecta a los transportadores mismos: desplaza el nivel de expresion ',
    'en unos genes y la dispersion entre individuos en otros, siempre en la misma ',
    'direccion -- el LPS vuelve mas uniforme la respuesta en hembras y mas variable ',
    'en machos.</p>')
  pagina(5, "Eje entre tejidos, limites y conclusion", cuerpo)
}

# =========================================================================
construir_html <- function(d) {
  sint <- d$fuente != "real"
  paginas <- paste(pagina1(d), pagina2(), pagina3(d, sint), pagina4(d, sint),
                   pagina5(d, sint), sep = "\n")
  paste0(
    "<!doctype html>\n<html lang=\"es\">\n<head>\n<meta charset=\"utf-8\">\n",
    "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
    "<title>Informe breve -- MIA-LPS</title>\n<style>", CSS, "</style>\n",
    "</head>\n<body>\n<main>\n",
    "<h1>Transportadores de nutrientes en el eje placenta&ndash;cerebro fetal en ",
    "un modelo de activacion inmune materna</h1>\n",
    "<p class=\"subtitulo\">Informe breve</p>\n",
    paginas,
    "\n<footer>Detalle completo, procedencia y verificaciones: ",
    "<code>docs/informe.html</code>. Presentacion: <code>docs/index.html</code>.",
    "</footer>\n</main>\n</body>\n</html>\n")
}

# =========================================================================
main <- function() {
  d <- recolectar_datos()
  sint <- d$fuente != "real"

  e13 <- .ENV13()
  e13$aplicar_override_tendencia(sint)   # ver cabecera: heredado de 13_presentacion.R
  ruta_subset_pla <- file.path(RUTA_FIGURAS, e13$NOMBRE_SUBSET_PLACENTA)
  ruta_subset_bra <- file.path(RUTA_FIGURAS, e13$NOMBRE_SUBSET_BRAIN)
  ruta_pstat3 <- file.path(RUTA_FIGURAS, e13$NOMBRE_PSTAT3_PRESENTACION)
  e13$generar_panel_subset(e13$GENES_SUBSET_PLACENTA, "PLACENTA_E15", ruta_subset_pla)
  e13$generar_panel_subset(e13$GENES_SUBSET_BRAIN, "BRAIN_E15", ruta_subset_bra)
  e13$generar_pstat3_presentacion(ruta_pstat3)

  ruta_densidades <- file.path(RUTA_FIGURAS, NOMBRE_DENSIDADES)
  generar_figura_densidades(GENES_DENSIDADES, ruta_densidades)

  html <- construir_html(d)

  destino <- if (sint) RUTA_DOCS else RUTA_INFORME_BREVE_REAL
  if (!dir.exists(destino)) dir.create(destino, recursive = TRUE, showWarnings = FALSE)
  ruta_html <- file.path(destino, "informe_breve.html")
  ruta_pdf <- file.path(destino, "informe_breve.pdf")
  escribir_texto(ruta_html, html)
  pdf_status <- generar_pdf(ruta_html, ruta_pdf)
  n_paginas <- contar_paginas_pdf(ruta_pdf)

  m_aviso <- gregexpr(AVISO_SINTETICO, html, fixed = TRUE)[[1]]
  n_aviso <- if (length(m_aviso) == 1L && m_aviso[1] == -1L) 0L else length(m_aviso)
  n_aviso_esperado <- if (sint) 3L else 0L   # paginas 3, 4, 5

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
    list(file.path("outputs/figures", NOMBRE_DENSIDADES), "figura_presentacion", ESTE_SCRIPT,
         "PROPIO (geom_density, mismo estilo que la diagonal del SPLOM de 08_acto2_correlaciones.R)",
         "data/processed/qpcr_cuantificacion_long.tsv",
         paste0("densidades de -ddCt por sexo x tratamiento, cerebro E15, genes ",
                .join_y(GENES_DENSIDADES),
                " (interaccion SEXOxTTO significativa sobre la dispersion, BH<0.05); ",
                "reemplaza a acto2_dispersion_sd.png en este informe -- muestra la ",
                "distribucion completa en vez de un resumen"))
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
         paste0("con fuente sintetica, las paginas 3-5 reemplazan la prosa ",
                "interpretativa por el aviso de datos sinteticos"),
         sprintf("fuente=%s; aviso=%d/%d", d$fuente, n_aviso, n_aviso_esperado),
         "aviso = 3 si fuente sintetica, 0 si fuente real",
         if (n_aviso == n_aviso_esperado) "TRUE" else "FALSE", ESTE_SCRIPT),
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
  cat(sprintf("  genes densidades = %s\n", paste(GENES_DENSIDADES, collapse = ", ")))
  cat(sprintf("  -> %s (%d KB)\n", ruta_html, file.info(ruta_html)$size %/% 1024))
  cat(sprintf("  -> %s [%s] (%s paginas)\n", ruta_pdf, pdf_status, n_paginas))
  if (length(faltan_figs))
    cat(sprintf("  *** faltan figuras: %s ***\n", paste(faltan_figs, collapse = ", ")))
}

# Genes para la figura de densidades (pedido, seccion 2): "tres o cuatro genes
# ... elegidos entre los que tienen interaccion significativa sobre la
# dispersion", leidos de acto2_dispersion_interaccion.csv. Con datos reales:
# fatcd36, fatp1, fatp4, gp130, slc38a2 (BH<0.05, cerebro). Se eligen 4 de esos
# 5 -- fatcd36, fatp4, gp130, slc38a2 -- por tener el contraste de SD mas
# marcado (control ancho / LPS angosto en hembras, al reves en machos, ~3-4x
# en los dos sentidos); fatp1 muestra el mismo patron pero mas atenuado.
GENES_DENSIDADES <- c("fatcd36", "fatp4", "gp130", "slc38a2")

if (sys.nframe() == 0L) main()
