<#
    run_all.ps1 — Corre el pipeline completo del reanálisis MIA-LPS, en orden, R y Python.

    ESQUELETO (T0): define el orden y el andamiaje. El cuerpo real de cada paso se completa
    a medida que se implementan los scripts (T1..T11). Cada bloque numerado ejecuta primero
    la versión R y después la Python del mismo script.

    Uso:
        .\run_all.ps1                     # todo, R y Python, con datos sintéticos si no hay crudos
        .\run_all.ps1 -Only python        # solo la implementación Python
        .\run_all.ps1 -Only R             # solo la implementación R
        .\run_all.ps1 -FromSynthetic      # forzar datos sintéticos aunque existan los crudos
        .\run_all.ps1 -SkipReport         # saltear 12_informe (HTML/PDF)

    Requisitos: ver README.md §2. No usa make ni utilidades Unix.
#>

[CmdletBinding()]
param(
    [ValidateSet('both', 'R', 'python')]
    [string]$Only = 'both',
    [switch]$FromSynthetic,
    [switch]$SkipReport
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

# --- Resolución de rutas (todo relativo a la ubicación de este script) ---
$RepoRoot = $PSScriptRoot
$RscriptCandidates = @(
    (Get-Command Rscript.exe -ErrorAction SilentlyContinue).Source,
    'C:\Program Files\R\R-4.6.1\bin\Rscript.exe'
) | Where-Object { $_ -and (Test-Path $_) }
$Rscript = $RscriptCandidates | Select-Object -First 1
$Python = (Get-Command python -ErrorAction SilentlyContinue).Source

Write-Host "Repo:    $RepoRoot"
Write-Host "Rscript: $Rscript"
Write-Host "Python:  $Python"
Write-Host ("Modo:    Only={0}  FromSynthetic={1}  SkipReport={2}" -f $Only, $FromSynthetic, $SkipReport)

# --- Orden de los scripts (sin extensión) ---
$Steps = @(
    '00_config',
    '01_generar_sinteticos',
    '02_ingesta_qc',
    '03_elisa',
    '04_qpcr_cuantificacion',
    '05_qpcr_modelos',
    '06_pstat3',
    '07_figuras_acto1',
    '08_acto2_correlaciones',
    '09_acto2_dispersion',
    '10_acto2_simulacion',
    '11_sensibilidad',
    '12_informe',
    '99_verificar'
)

function Invoke-Step {
    param([string]$Step)

    if ($SkipReport -and $Step -eq '12_informe') {
        Write-Host "  (salteado: $Step)" -ForegroundColor DarkYellow
        return
    }

    if ($Only -in @('both', 'R')) {
        $rPath = Join-Path $RepoRoot "R\$Step.R"
        if (Test-Path $rPath) {
            Write-Host ">> R      $Step" -ForegroundColor Cyan
            & $Rscript $rPath
            if ($LASTEXITCODE -ne 0) { throw "R/$Step.R terminó con código $LASTEXITCODE" }
        }
        else {
            Write-Host "   (aún no existe R\$Step.R)" -ForegroundColor DarkGray
        }
    }

    if ($Only -in @('both', 'python')) {
        $pyPath = Join-Path $RepoRoot "python\$Step.py"
        if (Test-Path $pyPath) {
            Write-Host ">> python $Step" -ForegroundColor Green
            & $Python $pyPath
            if ($LASTEXITCODE -ne 0) { throw "python/$Step.py terminó con código $LASTEXITCODE" }
        }
        else {
            Write-Host "   (aún no existe python\$Step.py)" -ForegroundColor DarkGray
        }
    }
}

$sw = [System.Diagnostics.Stopwatch]::StartNew()
foreach ($step in $Steps) { Invoke-Step -Step $step }
$sw.Stop()

Write-Host ""
Write-Host ("Pipeline terminado en {0:n1} min." -f $sw.Elapsed.TotalMinutes) -ForegroundColor White
Write-Host "Revisar: outputs/tables/verificaciones.csv y logs/corrida_<fecha>.txt"
