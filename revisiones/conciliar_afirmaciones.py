# conciliar_afirmaciones.py -- herramienta de auditoria de un solo uso, NO parte
# del pipeline (por eso no tiene gemelo en R, no escribe nada en outputs/ y no
# figura en el checklist de AGENTS.md 1).
#
# Punto 5 de pedidos/cambios_revision_codex.md: revisiones/AUDITORIA.md rastreo
# cinco afirmaciones del informe hasta su tabla, pero no pudo comprobar los
# numeros porque las tablas no estan en el clon. Este script los comprueba:
# recalcula cada numero desde outputs/tables/R/*.csv y verifica que la frase
# exacta que quedo en el informe (proyeccion sin figuras del render R) contenga
# ese numero, con el mismo redondeo.
#
# Uso, despues de una corrida con datos reales:
#     python revisiones/conciliar_afirmaciones.py
# El resultado esta transcripto en revisiones/conciliacion_auditoria.md.
#
import csv
import html
import math
import pathlib
import re

T = pathlib.Path("outputs/tables/R")
_raw = pathlib.Path("outputs/intermediate/render/R/informe.textonly.html").read_text(
    encoding="utf-8")
# Mismas sustituciones que hace el navegador al mostrar el informe: se quitan las
# etiquetas, se resuelven las entidades y se normaliza el espacio, para poder
# buscar el texto tal como se lee.
TXT = html.unescape(re.sub(r"<[^>]+>", " ", _raw))
TXT = TXT.replace("×", "x").replace("Δ", "Delta").replace("ρ", "rho")
TXT = re.sub(r"\s+", " ", TXT)


def tab(n):
    with open(T / n, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def half_up(x, nd=0):
    m = 10 ** nd
    return math.floor(abs(x) * m + 0.5) / m * (1 if x >= 0 else -1)


def dice(frag):
    """El informe contiene literalmente este fragmento?"""
    return frag in TXT


print("=" * 78)
print("1. Fisher de deteccion de IL-6 en suero materno (MS)")
f = [r for r in tab("elisa_fisher_deteccion.csv") if r["bloque"] == "MS"]
assert len(f) == 1
r = f[0]
print("   fila CSV : bloque=%s estrato=%s control=%s/%s lps=%s/%s p_valor=%s"
      % (r["bloque"], r["estrato"], r["control_detectado"], r["control_n"],
         r["lps_detectado"], r["lps_n"], r["p_valor"]))
frag = ("detectable en %s/%s madres LPS frente a %s/%s control (Fisher p = %s"
        % (r["lps_detectado"], r["lps_n"], r["control_detectado"], r["control_n"],
           r["p_valor"]))
print("   informe  :", frag)
print("   VERDICTO :", "COINCIDE" if dice(frag) else "NO COINCIDE")
print("   redondeo : ninguno; el informe imprime el string del CSV tal cual (%.6e)")

print("=" * 78)
print("2. Liquido amniotico: ningun contraste alcanza significancia")
pp = tab("elisa_petopeto_la.csv")
fi = [x for x in tab("elisa_fisher_deteccion.csv") if x["bloque"] == "LA"]
for x in pp:
    print("   Peto-Peto %-7s p_valor=%s  (chisq=%s df=%s)"
          % (x["estrato"], x["p_valor"], x["chisq"], x["df"]))
for x in fi:
    print("   Fisher    %-7s p_valor=%s" % (x["estrato"], x["p_valor"]))
todos = [float(x["p_valor"]) for x in pp] + [float(x["p_valor"]) for x in fi]
h = [x for x in pp if x["estrato"] == "HEMBRA"][0]["p_valor"]
m = [x for x in pp if x["estrato"] == "MACHO"][0]["p_valor"]
frag2 = "ningun contraste alcanzo significancia (Peto-Peto"
print("   informe  : ... %s ... p = %s ... p = %s" % (frag2, h, m))
ok2 = (dice(frag2) and dice("p = %s" % h) and dice("p = %s" % m)
       and all(p >= 0.05 for p in todos))
print("   contrastes LA con p < 0.05: %d de %d" % (sum(p < 0.05 for p in todos),
                                                   len(todos)))
print("   VERDICTO :", "COINCIDE" if ok2 else "NO COINCIDE")

print("=" * 78)
print("3. Interaccion SEXO x TTO por tejido")
cl = tab("qpcr_modelos_clasificacion.csv")
mod = [r for r in cl if r["via"] == "modelo"]
bra = [r for r in mod if r["TEJIDO"] == "BRAIN_E15"]
pla = [r for r in mod if r["TEJIDO"] == "PLACENTA_E15"]
bra_int = [r for r in bra if r["interaccion_significativa"] == "TRUE"]
pla_int = [r for r in pla if r["interaccion_significativa"] == "TRUE"]
print("   gen x tejido modelados (via=modelo): %d  (placenta %d, cerebro %d)"
      % (len(mod), len(pla), len(bra)))
print("   interaccion_significativa=TRUE: cerebro %d (%s); placenta %d"
      % (len(bra_int), ", ".join(sorted(x["GEN"] for x in bra_int)), len(pla_int)))
lim = [r for r in bra if r["interaccion_significativa"] != "TRUE" and r["p_SEXOxTTO"]]
print("   cerebro modelado SIN interaccion: %s"
      % ", ".join("%s p=%s" % (x["GEN"], x["p_SEXOxTTO"]) for x in lim))
bh_ok = [r for r in bra_int if float(r["p_SEXOxTTO_BH"]) < 0.05]
print("   BH (p_SEXOxTTO_BH < .05) entre los %d: %d  -> valores %s"
      % (len(bra_int), len(bh_ok), sorted({x["p_SEXOxTTO_BH"] for x in bra_int})))
frags3 = [
    "De %d gen x tejido modelados, la interaccion SEXOxTTO es significativa en %d de "
    "cerebro E15 y %d de placenta E15" % (len(mod), len(bra_int), len(pla_int)),
    "%d genes mostraron interaccion SEXOxTTO significativa (%d de %d sobreviven a la "
    "correccion BH)" % (len(bra_int), len(bh_ok), len(bra_int)),
    "Placenta. Ningun gen mostro interaccion SEXOxTTO.",
]
for x in frags3:
    print("   informe  : [%s] %s" % ("OK" if dice(x) else "FALTA", x))
print("   VERDICTO :", "COINCIDE" if all(dice(x) for x in frags3) else "NO COINCIDE")

print("=" * 78)
print("4. Delta rho Control vs LPS: ninguno significativo")
tc = tab("acto2_test_correlaciones.csv")
items = sorted({r["ITEM"] for r in tc})
estr = sorted({r["ESTRATO"] for r in tc})
ps = [(float(r["p_bw"]), r["ITEM"], r["ESTRATO"]) for r in tc]
sig = [x for x in ps if x[0] < 0.05]
mn = min(ps)
print("   filas=%d  = %d items x %d estratos %s"
      % (len(tc), len(items), len(estr), estr))
print("   items: %s" % ", ".join(items))
print("   p_bw < 0.05: %d de %d" % (len(sig), len(ps)))
fila_min = [r for r in tc if r["ITEM"] == mn[1] and r["ESTRATO"] == mn[2]][0]
print("   minimo: ITEM=%s ESTRATO=%s p_bw=%s (delta_rho=%s)"
      % (mn[1], mn[2], fila_min["p_bw"], fila_min["delta_rho"]))
frag4 = ("(%d de %d; minimo %s en hembras, p = %s"
         % (len(sig), len(ps), mn[1], fila_min["p_bw"]))
print("   informe  : [%s] %s" % ("OK" if dice(frag4) else "FALTA", frag4))
print("   VERDICTO :", "COINCIDE" if dice(frag4) and mn[2] == "HEMBRA"
      else "NO COINCIDE")

print("=" * 78)
print("5. PC1: varianza explicada por tejido")
pv = tab("acto2_sensibilidad_pca_varianza.csv")
for tej in ("PLACENTA_E15", "BRAIN_E15"):
    r = [x for x in pv if x["TEJIDO"] == tej and x["PC"] == "1"][0]
    p = float(r["prop_var"])
    print("   %-13s prop_var=%s -> x100 = %.6f -> half-up(0) = %d"
          % (tej, r["prop_var"], p * 100, int(half_up(p * 100))))
pla = float([x for x in pv if x["TEJIDO"] == "PLACENTA_E15" and x["PC"] == "1"][0]
            ["prop_var"]) * 100
bra = float([x for x in pv if x["TEJIDO"] == "BRAIN_E15" and x["PC"] == "1"][0]
            ["prop_var"]) * 100
frag5 = ("El eigengene (PC1) explica el %d %% de la varianza de los transportadores "
         "en placenta y el %d %% en cerebro"
         % (int(half_up(pla)), int(half_up(bra))))
print("   informe  : [%s] %s" % ("OK" if dice(frag5) else "FALTA", frag5))
print("   VERDICTO :", "COINCIDE" if dice(frag5) else "NO COINCIDE")
print("=" * 78)
