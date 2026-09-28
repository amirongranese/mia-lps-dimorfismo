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

    `-Only R|python` hace UNA sola pasada de ese lenguaje (00..11 + 12 + 99), sin
    98_comparacion: sirve para iterar una implementacion, y para que quien SOLO
    tenga R (o solo Python) pueda verificar su mitad del pipeline. En ese modo
    99_verificar marca los chequeos R<->Python como NO_EJECUTADA (via la variable
    de entorno MIA_LPS_UNICA_IMPL) y NO imprime "TODAS LAS VERIFICACIONES
    PASARON", porque parte quedo sin ejecutar. La validacion completa necesita
    `both`.

    Interpretes: no hay rutas fijas. Rscript y python se buscan en el PATH, en el
    registro de R (R-core), en el venv del repo (.venv) y en las instalaciones
    estandar de Windows; ver README 2. Se exige SOLO el interprete que el modo
    elegido va a usar (`-Only R` no necesita Python).

    Uso:
        .\run_all.ps1                     # todo, R y Python, con verificacion final
        .\run_all.ps1 -Only python        # solo Python (analisis + informe + 99)
        .\run_all.ps1 -Only R             # solo R (analisis + informe + 99)
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

# --- Resolucion de interpretes -------------------------------------------
# Sin rutas fijas (punto 4.4 de pedidos/cambios_revision_codex.md): las rutas del
# README son EJEMPLOS. Se busca, en orden: PATH -> registro/venv -> instalaciones
# estandar de Windows. Cada candidato se valida ejecutandolo, asi los alias de
# app-execution de la Store ("python3" en WindowsApps, que solo abre la Store) no
# se toman por interpretes.

function Test-Interprete {
    param([string]$Exe, [string[]]$Argumentos)
    if (-not $Exe -or -not (Test-Path -LiteralPath $Exe)) { return $false }
    try {
        $null = & $Exe @Argumentos
        return ($LASTEXITCODE -eq 0)
    }
    catch { return $false }
}

function Resolve-Rscript {
    $cands = @()
    foreach ($n in @('Rscript.exe', 'Rscript')) {
        $c = Get-Command $n -ErrorAction SilentlyContinue
        if ($c -and $c.Source) { $cands += $c.Source }
    }
    # Registro de R-core: lo escribe el instalador oficial de R en Windows.
    foreach ($k in @('HKLM:\SOFTWARE\R-core\R', 'HKCU:\SOFTWARE\R-core\R',
            'HKLM:\SOFTWARE\WOW6432Node\R-core\R')) {
        if (-not (Test-Path $k)) { continue }
        $prop = Get-ItemProperty -Path $k -ErrorAction SilentlyContinue
        if ($prop -and $prop.PSObject.Properties['InstallPath'] -and $prop.InstallPath) {
            $cands += (Join-Path $prop.InstallPath 'bin\Rscript.exe')
        }
    }
    # Instalaciones estandar. Orden por nombre descendente => version mas nueva
    # primero mientras los numeros tengan un digito (R-4.6.1 antes que R-4.4.1).
    foreach ($r in @("$env:ProgramFiles\R", "${env:ProgramFiles(x86)}\R",
            "$env:LOCALAPPDATA\Programs\R")) {
        if (-not (Test-Path $r)) { continue }
        $cands += (Get-ChildItem -LiteralPath $r -Directory -ErrorAction SilentlyContinue |
            Sort-Object Name -Descending |
            ForEach-Object { Join-Path $_.FullName 'bin\Rscript.exe' })
    }
    foreach ($c in $cands) {
        if (Test-Interprete -Exe $c -Argumentos @('-e', 'invisible(1)')) { return $c }
    }
    return $null
}

function Resolve-Python {
    $cands = @()
    # 1. venv del repo: es lo que documenta el README 2.2.
    $cands += (Join-Path $RepoRoot '.venv\Scripts\python.exe')
    # 2. PATH.
    foreach ($n in @('python.exe', 'python', 'python3.exe', 'python3')) {
        $c = Get-Command $n -ErrorAction SilentlyContinue
        if ($c -and $c.Source) { $cands += $c.Source }
    }
    # 3. Lanzador `py`: sabe donde esta el interprete real aunque no este en PATH.
    $py = Get-Command 'py.exe' -ErrorAction SilentlyContinue
    if ($py -and $py.Source) {
        try {
            $exe = & $py.Source -3 -c 'import sys; print(sys.executable)'
            if ($LASTEXITCODE -eq 0 -and $exe) { $cands += ([string]$exe).Trim() }
        }
        catch { }
    }
    # 4. Instalaciones estandar (per-user y per-machine).
    foreach ($r in @("$env:LOCALAPPDATA\Programs\Python", "$env:ProgramFiles\Python",
            "${env:ProgramFiles(x86)}\Python", 'C:\Python')) {
        if (-not (Test-Path $r)) { continue }
        $cands += (Get-ChildItem -LiteralPath $r -Directory -ErrorAction SilentlyContinue |
            Sort-Object Name -Descending |
            ForEach-Object { Join-Path $_.FullName 'python.exe' })
    }
    foreach ($c in $cands) {
        if (Test-Interprete -Exe $c -Argumentos @('-c', 'import sys; sys.exit(0)')) { return $c }
    }
    return $null
}

# Solo se exige el interprete que este modo va a usar (punto 4.2).
$NecesitaR = ($Only -eq 'both' -or $Only -eq 'R')
$NecesitaPy = ($Only -eq 'both' -or $Only -eq 'python')

$Rscript = if ($NecesitaR) { Resolve-Rscript } else { $null }
$Python = if ($NecesitaPy) { Resolve-Python } else { $null }

if ($NecesitaR -and -not $Rscript) {
    throw ("No se encontro Rscript. Instalar R >= 4.4 y dejarlo en el PATH o en la " +
        "instalacion estandar de Windows (ver README.md 2.1).")
}
if ($NecesitaPy -and -not $Python) {
    throw ("No se encontro python. Instalar Python >= 3.11 y dejarlo en el PATH, o " +
        "crear el venv del repo en .venv\ (ver README.md 2.2).")
}

if ($FromSynthetic) { $env:MIA_LPS_FORZAR_SINTETICO = '1' }
else { Remove-Item Env:\MIA_LPS_FORZAR_SINTETICO -ErrorAction SilentlyContinue }

$SinUso = '(no se usa en este modo)'
Write-Host "Repo:    $RepoRoot"
Write-Host ("Rscript: {0}" -f $(if ($Rscript) { $Rscript } else { $SinUso }))
Write-Host ("Python:  {0}" -f $(if ($Python) { $Python } else { $SinUso }))
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
    Remove-Item Env:\MIA_LPS_UNICA_IMPL -ErrorAction SilentlyContinue
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
    # -Only R | python : una sola pasada, sin 98_comparacion, pero CON 99_verificar
    # (punto 4.3). MIA_LPS_UNICA_IMPL le avisa a 99 que los chequeos R<->Python no
    # se pueden ejecutar en esta corrida: quedan NO_EJECUTADA, ni pase ni fallo.
    $env:MIA_LPS_UNICA_IMPL = $Only
    $steps = $StepsBase
    if (-not $SkipReport) { $steps += '12_informe' }
    Write-Host ""
    Write-Host "== Pasada unica: $Only (sin 98_comparacion) ==" -ForegroundColor Magenta
    Invoke-Lang -Lang $Only -Steps $steps

    if (-not $SkipReport -and -not $SkipVerify) {
        Write-Host ""
        Write-Host "== Verificacion de una sola implementacion: $Only 99_verificar ==" `
            -ForegroundColor Magenta
        Invoke-Lang -Lang $Only -Steps @('99_verificar')
    }
    else {
        Write-Host ""
        Write-Host "   (99_verificar salteado)" -ForegroundColor DarkYellow
    }

    Write-Host ""
    Write-Host "   Modo -Only: sin comparacion R<->Python; esos chequeos quedan NO_EJECUTADA." `
        -ForegroundColor DarkYellow
    Write-Host "   Para la validacion completa: .\run_all.ps1  (sin -Only)." `
        -ForegroundColor DarkYellow
}

$sw.Stop()
Write-Host ""
Write-Host ("Pipeline terminado en {0:n1} min." -f $sw.Elapsed.TotalMinutes) `
    -ForegroundColor White
Write-Host "Revisar: outputs/tables/verificaciones.csv, docs/informe.html y logs/corrida_<fecha>.txt"
