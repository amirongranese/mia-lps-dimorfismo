# 00_config.R -- Constantes, rutas y utilidades deterministas compartidas.
#
# Por que existe este archivo: el pipeline se implementa dos veces (R y Python) y
# todo script arranca haciendo source() de esta configuracion. Aca se fijan (a)
# las rutas, siempre relativas a la raiz del repo resuelta desde la ubicacion del
# script, (b) la semilla global, (c) las listas de genes / grupos / tejidos, y
# (d) un generador pseudoaleatorio propio que produce EXACTAMENTE la misma
# secuencia que python/00_config.py, para que los datos sinteticos salgan
# identicos entre lenguajes.
#
# Este archivo NO hace analisis. Su unico efecto colateral es crear las carpetas
# de salida si no existen (idempotente).

# ---------------------------------------------------------------------------
# Raiz del repo: se deduce del propio archivo (prohibido hardcodear rutas
# absolutas). Funciona con Rscript (--file=) y con source() desde otro script.
# ---------------------------------------------------------------------------
.detectar_raiz_repo <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m) == 1L) {
    return(dirname(dirname(normalizePath(sub("^--file=", "", m), mustWork = FALSE))))
  }
  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (!is.null(of)) {
      return(dirname(dirname(normalizePath(of, mustWork = FALSE))))
    }
  }
  wd <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(wd) %in% c("R", "python")) return(dirname(wd))
  wd
}

RAIZ_REPO <- .detectar_raiz_repo()

# --- Rutas de datos y salidas (todas relativas a RAIZ_REPO) -----------------
RUTA_DATOS_RAW   <- file.path(RAIZ_REPO, "data", "raw")
RUTA_DATOS_SINT  <- file.path(RAIZ_REPO, "data", "synthetic")
RUTA_DATOS_PROC  <- file.path(RAIZ_REPO, "data", "processed")
RUTA_OUT         <- file.path(RAIZ_REPO, "outputs")
RUTA_FIGURAS     <- file.path(RUTA_OUT, "figures")
RUTA_TABLAS      <- file.path(RUTA_OUT, "tables")
RUTA_TABLAS_R    <- file.path(RUTA_TABLAS, "R")
RUTA_TABLAS_PY   <- file.path(RUTA_TABLAS, "python")
RUTA_INTERMEDIOS <- file.path(RUTA_OUT, "intermediate")
RUTA_DOCS        <- file.path(RAIZ_REPO, "docs")
RUTA_LOGS        <- file.path(RAIZ_REPO, "logs")

# Crear lo regenerable si falta (un clon limpio trae solo data/synthetic).
for (.d in c(RUTA_DATOS_SINT, RUTA_DATOS_PROC, RUTA_FIGURAS, RUTA_TABLAS_R,
             RUTA_TABLAS_PY, RUTA_INTERMEDIOS, RUTA_DOCS, RUTA_LOGS)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE, showWarnings = FALSE)
}

# --- Archivos de entrada y hojas ------------------------------------------
ARCHIVO_QPCR   <- "Raw data CTs.xlsx";              HOJA_QPCR   <- "Sheet1"
ARCHIVO_ELISA  <- "ELISA IL6 2026 Dosis 100.xlsx";  HOJA_ELISA  <- "Sueros y LA"
HOJA_CURVA     <- "CURVA IL6"
ARCHIVO_PSTAT3 <- "pstat3 placenta.xlsx";           HOJA_PSTAT3 <- "Sheet1"

# run_all.ps1 -FromSynthetic exporta MIA_LPS_FORZAR_SINTETICO=1 para correr el
# pipeline sobre data/synthetic/ aunque existan los crudos (util para el chequeo
# de reproducibilidad "corre de punta a punta sin datos reales").
.forzar_sintetico <- function() {
  !(trimws(Sys.getenv("MIA_LPS_FORZAR_SINTETICO", "")) %in% c("", "0"))
}

# ruta_datos(): usa el crudo real de data/raw/ si existe; si no, el sintetico
# versionado de data/synthetic/. Nunca escribe en data/raw/. Un clon limpio
# (sin datos crudos) cae automaticamente al sintetico.
ruta_datos <- function(nombre_archivo) {
  real <- file.path(RUTA_DATOS_RAW, nombre_archivo)
  if (file.exists(real) && !.forzar_sintetico()) return(real)
  file.path(RUTA_DATOS_SINT, nombre_archivo)
}

fuente_datos <- function(nombre_archivo) {
  if (.forzar_sintetico()) return("sintetico")
  if (file.exists(file.path(RUTA_DATOS_RAW, nombre_archivo))) "real" else "sintetico"
}

# --- Semilla global -----------------------------------------------------
SEMILLA <- 20260101L

# --- Genes y su clasificacion (orden canonico del archivo crudo) -----------
GENES <- c("fatcd36", "fatp1", "fatp4", "glut1", "glut3",
           "gp130", "il6", "il6R", "slc38a1", "slc38a2")
GENES_TRANSPORTADORES <- c("fatcd36", "fatp1", "fatp4", "glut1", "glut3",
                           "slc38a1", "slc38a2")  # 7 transportadores de nutrientes
GENES_VIA_IL6 <- c("gp130", "il6", "il6R")         # 3 de la via IL-6
GEN_HOUSEKEEPING <- "rsp29"

# Conjuntos de variables para los SPLOM de co-expresion del Acto 2 (T7,
# extension pedida en pedidos/cambios_acto2_correlaciones_por_sexo.md): el
# conjunto es distinto por tejido porque il6 es cuantificable en placenta pero
# no en cerebro (calibrador HEMBRA_CONTROL 0/9, D7); il6R no entra en ningun
# SPLOM (deteccion insuficiente: 3/7/3/3 en cerebro).
GENES_SPLOM_PLACENTA <- c(GENES_TRANSPORTADORES, "il6", "gp130")   # 9 variables
GENES_SPLOM_BRAIN    <- c(GENES_TRANSPORTADORES, "gp130")          # 8 variables

stopifnot(
  setequal(GENES, c(GENES_TRANSPORTADORES, GENES_VIA_IL6)),
  length(GENES_TRANSPORTADORES) == 7L,
  length(GENES_VIA_IL6) == 3L,
  length(unique(GENES)) == 10L,
  length(GENES_SPLOM_PLACENTA) == 9L, length(GENES_SPLOM_BRAIN) == 8L,
  setequal(GENES_SPLOM_PLACENTA, c(GENES_TRANSPORTADORES, "il6", "gp130")),
  setequal(GENES_SPLOM_BRAIN, c(GENES_TRANSPORTADORES, "gp130"))
)

# --- Tejidos ----------------------------------------------------------
TEJIDOS_E15     <- c("PLACENTA_E15", "BRAIN_E15")
TEJIDO_EXCLUIDO <- "BRAIN_P1"   # se excluye en 02_ingesta_qc (se analiza en otro informe)
TEJIDOS_TODOS   <- c(TEJIDOS_E15, TEJIDO_EXCLUIDO)

# --- Diseno de grupos -------------------------------------------------
NIVELES_SEXO <- c("HEMBRA", "MACHO")
NIVELES_TTO  <- c("CONTROL", "LPS")
# D1: el calibrador de la cuantificacion relativa es HEMBRA CONTROL.
CALIBRADOR <- c(SEXO = "HEMBRA", TTO = "CONTROL")

# --- Paleta de los boxplots de expresion / pSTAT3 (Acto 1, 07_figuras_acto1) --
# Unica fuente de estos colores -- 07_figuras_acto1.{R,py} los referencia, no
# los reescribe a mano. Control por sexo (celeste, 2 tonos); LPS por sexo x
# via metabolica (una paleta por via) para que el ojo asocie color con la
# ruta; pSTAT3 (senalizacion de IL-6) reusa la paleta IL6 para su LPS.
COL_CTRL <- c(HEMBRA = "#AEDCF0", MACHO = "#6BAED6")
COL_LPS  <- list(
  GLUCOSA     = c(HEMBRA = "#F09EC8", MACHO = "#D6317F"),
  AMINOACIDOS = c(HEMBRA = "#8FD9B6", MACHO = "#2E9E6B"),
  LIPIDOS     = c(HEMBRA = "#FBC98A", MACHO = "#E08214"),
  IL6         = c(HEMBRA = "#C5A3E0", MACHO = "#7B4EA8")
)

EMOJI_HEMBRA <- "♀"
EMOJI_MACHO  <- "♂"
VS16         <- "️"   # variation selector que traen las NOMINACION de los crudos

normalizar_sexo <- function(x) {
  u <- toupper(trimws(as.character(x)))
  if (startsWith(u, "H") || startsWith(u, "F") || u == EMOJI_HEMBRA) return("HEMBRA")
  if (startsWith(u, "M") || u == EMOJI_MACHO) return("MACHO")
  stop(sprintf("SEXO no reconocido: %s", x))
}

normalizar_tto <- function(x) {
  u <- toupper(trimws(as.character(x)))
  if (startsWith(u, "C")) return("CONTROL")
  if (startsWith(u, "L")) return("LPS")
  stop(sprintf("TTO no reconocido: %s", x))
}

# Etiqueta GRUPO tal como aparece en los crudos, reproduciendo el espaciado
# inconsistente real: 'HEMBRA'+'CONTROL' va pegado, el resto lleva un espacio.
etiqueta_grupo <- function(sexo, tto) {
  s <- normalizar_sexo(sexo)
  t <- normalizar_tto(tto)
  emo <- if (s == "HEMBRA") EMOJI_HEMBRA else EMOJI_MACHO
  tt  <- if (t == "CONTROL") "Control" else "LPS"
  sep <- if (s == "HEMBRA" && t == "CONTROL") "" else " "
  paste0(emo, sep, tt)
}

# --- Tolerancias de comparacion R <-> Python (convencion del proyecto) -----
TOL_ESTADISTICO <- 1e-6   # estadisticos, medias, coeficientes
TOL_P_ITERATIVO <- 1e-4   # p-valores de tests iterativos / permutacion

# ---------------------------------------------------------------------------
# Generador pseudoaleatorio propio -- PARIDAD EXACTA R <-> Python
#
# PROPIO. No se usa para nada cientifico: unicamente para fabricar los datos
# sinteticos de 01_generar_sinteticos de forma identica en ambos lenguajes.
# LCG de 32 bits (constantes de Numerical Recipes) para el stream uniforme, y
# metodo polar de Marsaglia para los normales (solo sqrt y log, correctamente
# redondeados en R y Python -> los bits coinciden; se evita sin/cos por eso).
# Nota: 1664525 * (2^32 - 1) < 2^53, asi que la aritmetica en doubles es exacta.
# ---------------------------------------------------------------------------
.LCG_A <- 1664525
.LCG_C <- 1013904223
.LCG_M <- 4294967296   # 2^32

nuevo_rng <- function(semilla = SEMILLA) {
  estado <- as.numeric(semilla) %% .LCG_M
  cache <- 0
  tiene_cache <- FALSE

  unif1 <- function() {
    estado <<- (.LCG_A * estado + .LCG_C) %% .LCG_M
    estado / .LCG_M
  }

  norm1 <- function(media = 0, sd = 1) {
    if (tiene_cache) {
      tiene_cache <<- FALSE
      return(media + sd * cache)
    }
    repeat {
      v1 <- 2 * unif1() - 1
      v2 <- 2 * unif1() - 1
      w <- v1 * v1 + v2 * v2
      if (w > 0 && w < 1) break
    }
    f <- sqrt(-2 * log(w) / w)
    cache <<- v2 * f
    tiene_cache <<- TRUE
    media + sd * (v1 * f)
  }

  entero <- function(lo, hi) lo + floor(unif1() * (hi - lo + 1))

  # k elementos distintos de 'opciones', por Fisher-Yates parcial. La aritmetica
  # replica el indice 0-based de Python: paso R i (1..k) == paso Python i-1.
  elegir_sin_reemplazo <- function(opciones, k) {
    n <- length(opciones)
    idx <- seq_len(n)
    for (i in seq_len(k)) {
      j0 <- (i - 1) + floor(unif1() * (n - (i - 1)))  # objetivo 0-based en [i-1, n-1]
      j <- j0 + 1L
      tmp <- idx[i]; idx[i] <- idx[j]; idx[j] <- tmp
    }
    opciones[idx[seq_len(k)]]
  }

  list(unif1 = unif1, norm1 = norm1, entero = entero,
       elegir_sin_reemplazo = elegir_sin_reemplazo)
}

# ---------------------------------------------------------------------------
.resumen_config <- function() {
  cat("== 00_config.R ==\n")
  cat(sprintf("RAIZ_REPO        : %s\n", RAIZ_REPO))
  cat(sprintf("SEMILLA          : %d\n", SEMILLA))
  cat(sprintf("GENES (%d)       : %s\n", length(GENES), paste(GENES, collapse = ", ")))
  cat(sprintf("  transportadores: %s\n", paste(GENES_TRANSPORTADORES, collapse = ", ")))
  cat(sprintf("  via IL-6       : %s\n", paste(GENES_VIA_IL6, collapse = ", ")))
  cat(sprintf("  housekeeping   : %s\n", GEN_HOUSEKEEPING))
  cat(sprintf("TEJIDOS_E15      : %s   (excluido: %s)\n",
              paste(TEJIDOS_E15, collapse = ", "), TEJIDO_EXCLUIDO))
  for (s in NIVELES_SEXO) for (t in NIVELES_TTO) {
    cat(sprintf("  grupo %-6s %-7s -> '%s'\n", s, t, etiqueta_grupo(s, t)))
  }
  for (f in c(ARCHIVO_QPCR, ARCHIVO_ELISA, ARCHIVO_PSTAT3)) {
    cat(sprintf("  %-32s fuente actual: %s\n", f, fuente_datos(f)))
  }
  r <- nuevo_rng(SEMILLA)
  u <- round(c(r$unif1(), r$unif1(), r$unif1()), 10)
  cat(sprintf("RNG check unif1 x3: [%s]\n", paste(sprintf("%.10g", u), collapse = ", ")))
}

if (sys.nframe() == 0L) {
  .resumen_config()
}
