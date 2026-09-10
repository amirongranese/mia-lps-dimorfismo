<#
    run_all.ps1 -- Corre el pipeline completo del reanalisis MIA-LPS, R y Python.

    Estructura (T11): para poder byte-comparar de verdad R contra Python -- no solo
    los CSV numericos (eso lo cubre 98_comparacion) sino tambien los `.md` de copia
    unica y la proyeccion sin figuras del informe -- cada lenguaje tiene que dejar
    SUS PROPIOS `.md` en disco cuando corre `12_informe`. Con la pasada intercalada
    (R, luego Python, por paso) eso no pasa: Python pisa siempre lo de R. Por eso el
    modo `both` corre en pasadas separadas:

        0. PRIME  (python 00..11)   -- puebla outputs/tables/python/ para que el 98
                                       de la pasada R tenga con quien compararse.
        1. R      (R 00..11, 98, 12) -- deja los .md, docs/informe.html y el snapshot
                                       outputs/intermediate/render/R/ en version R.
        2. PYTHON (python 00..11, 98, 12) -- idem en version Python; docs/informe.html
                                       final queda en la version Python.
        3. VERIFY (R 99, luego python 99) -- 99_verificar compara render/R vs
                                       render/python, re-corre la concordancia de los
                                       CSV, chequea el checklist de AGENTS 1 y escribe
                                       logs/corrida_<fecha>.txt.

    `-Only R|python` hace UNA sola pasada de ese lenguaje (00..11 + 12), sin 98 ni
    99: sirve para iterar una implementacion cuando la otra ya dejo sus salidas.
    La validacion completa (98 + 99 + "TODAS LAS VERIFICACIONES PASARON") necesita
    `both`.

    Uso:
        .\run_all.ps1                     # todo, R y Python, con verificacion final
        .\run_all.ps1 -Only python        # solo Python (analisis + informe)
        .\run_all.ps1 -Only R             # solo R (analisis + informe)
        .\run_all.ps1 -FromSynthetic      # forzar data/synthetic/ aunque haya crudos
        .\run_all.ps1 -SkipReport         # saltear 12_informe y 99_verificar
        .\run_all.ps1 -SkipVerify         # correr todo menos 99_verificar

    Requisitos: ver README.md 2. No usa make ni utilidades Unix.
#>

[CmdletBinding()]
param(
    [ValidateSet('both', 'R', 'python')]
    [string]$Only = 'both',
    [switch]$FromSynthetic,
    [switch]$SkipReport,
    [switch]$SkipVerify
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

# --- Consola en UTF-8: los scripts imprimen '♀'/'→' etc.; sin esto, con la
#     salida redirigida (CI, `| tee`) Python cae en cp1252 y explota. ---
$env:PYTHONUTF8 = '1'
$env:PYTHONIOENCODING = 'utf-8'
try {
    [Console]::OutputEncoding = [System.Text.Encoding]::UTF8
    $OutputEncoding = [System.Text.Encoding]::UTF8
}
catch { }

# --- Resolucion de rutas (todo relativo a la ubicacion de este script) ---
$RepoRoot = $PSScriptRoot

$rscriptCmd = Get-Command Rscript.exe -ErrorAction SilentlyContinue
$RscriptCandidates = @()
if ($rscriptCmd) { $RscriptCandidates += $rscriptCmd.Source }
$RscriptCandidates += 'C:\Program Files\R\R-4.6.1\bin\Rscript.exe'
$RscriptCandidates += 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe'
$Rscript = $RscriptCandidates | Where-Object { $_ -and (Test-Path $_) } | Select-Object -First 1

$pythonCmd = Get-Command python -ErrorAction SilentlyContinue
$Python = if ($pythonCmd) { $pythonCmd.Source } else { $null }

if (-not $Rscript) { throw "No se encontro Rscript.exe (ver README.md 2)." }
if (-not $Python) { throw "No se encontro python en PATH (ver README.md 2)." }

if ($FromSynthetic) { $env:MIA_LPS_FORZAR_SINTETICO = '1' }
else { Remove-Item Env:\MIA_LPS_FORZAR_SINTETICO -ErrorAction SilentlyContinue }

Write-Host "Repo:    $RepoRoot"
Write-Host "Rscript: $Rscript"
Write-Host "Python:  $Python"
Write-Host ("Modo:    Only={0}  FromSynthetic={1}  SkipReport={2}  SkipVerify={3}" -f `
        $Only, [bool]$FromSynthetic, [bool]$SkipReport, [bool]$SkipVerify)

# --- Orden de los scripts (sin extension) ---
$StepsBase = @(
    '00_config', '01_generar_sinteticos', '02_ingesta_qc', '03_elisa',
    '04_qpcr_cuantificacion', '05_qpcr_modelos', '06_pstat3', '07_figuras_acto1',
    '08_acto2_correlaciones', '09_acto2_dispersion', '10_acto2_simulacion',
    '11_sensibilidad'
)

function Invoke-Lang {
    param(
        [ValidateSet('R', 'python')] [string]$Lang,
        [string[]]$Steps
    )
    foreach ($step in $Steps) {
        if ($Lang -eq 'R') {
            $path = Join-Path $RepoRoot "R\$step.R"
            $exe = $Rscript
            $tag = ">> R      $step"
            $col = 'Cyan'
        }
        else {
            $path = Join-Path $RepoRoot "python\$step.py"
            $exe = $Python
            $tag = ">> python $step"
            $col = 'Green'
        }
        if (-not (Test-Path $path)) {
            Write-Host "   (no existe $path)" -ForegroundColor DarkGray
            continue
        }
        Write-Host $tag -ForegroundColor $col
        & $exe $path
        if ($LASTEXITCODE -ne 0) { throw "$Lang / $step termino con codigo $LASTEXITCODE" }
    }
}

$sw = [System.Diagnostics.Stopwatch]::StartNew()

if ($Only -eq 'both') {
    $stepsFull = $StepsBase + @('98_comparacion')
    if (-not $SkipReport) { $stepsFull += '12_informe' }

    Write-Host ""
    Write-Host "== Pasada 0/3: PRIME (python 00..11) ==" -ForegroundColor Magenta
    Invoke-Lang -Lang python -Steps $StepsBase

    Write-Host ""
    Write-Host "== Pasada 1/3: R ==" -ForegroundColor Magenta
    Invoke-Lang -Lang R -Steps $stepsFull

    Write-Host ""
    Write-Host "== Pasada 2/3: PYTHON ==" -ForegroundColor Magenta
    Invoke-Lang -Lang python -Steps $stepsFull

    if (-not $SkipReport -and -not $SkipVerify) {
        Write-Host ""
        Write-Host "== Pasada 3/3: VERIFY (99_verificar) ==" -ForegroundColor Magenta
        Invoke-Lang -Lang R -Steps @('99_verificar')
        Invoke-Lang -Lang python -Steps @('99_verificar')
    }
    else {
        Write-Host ""
        Write-Host "   (99_verificar salteado)" -ForegroundColor DarkYellow
    }
}
else {
    # -Only R | python : una sola pasada, sin 98 ni 99.
    $steps = $StepsBase
    if (-not $SkipReport) { $steps += '12_informe' }
    Write-Host ""
    Write-Host "== Pasada unica: $Only (sin 98_comparacion ni 99_verificar) ==" `
        -ForegroundColor Magenta
    Invoke-Lang -Lang $Only -Steps $steps
    Write-Host ""
    Write-Host "   Modo -Only: sin comparacion R<->Python ni verificacion final." `
        -ForegroundColor DarkYellow
    Write-Host "   Para la validacion completa: .\run_all.ps1  (sin -Only)." `
        -ForegroundColor DarkYellow
}

$sw.Stop()
Write-Host ""
Write-Host ("Pipeline terminado en {0:n1} min." -f $sw.Elapsed.TotalMinutes) `
    -ForegroundColor White
Write-Host "Revisar: outputs/tables/verificaciones.csv, docs/informe.html y logs/corrida_<fecha>.txt"
