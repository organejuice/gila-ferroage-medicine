param(
  [string]$ProjectRoot = ".",
  [string]$Rscript = "D:\R\R-4.3.3\bin\Rscript.exe",
  [string]$Python = "D:\conda\python.exe",
  [string]$RLibrary = "C:\Users\administer\.codex\visualizations\2026\09\23\01a0ce11-eaa9-7cb1-aa30-49773cbedb70\g6_R_library"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$Rscript = (Resolve-Path -LiteralPath $Rscript).Path
$Python = (Resolve-Path -LiteralPath $Python).Path
$env:R_LIBS_USER = $RLibrary
$analysisExit = 0
$reportExit = 0
$expectedPreflightStop = $false
Push-Location $ProjectRoot
try {
  & $Rscript --vanilla "src/73_run_g6_primary_analysis.R" "."
  $analysisExit = $LASTEXITCODE
  & $Python "src/76_write_g6_report.py" "."
  $reportExit = $LASTEXITCODE
  $preflightPath = Join-Path $ProjectRoot "results/2026-09-24-disease-side-v3.5-g6-primary/G6_effect_blind_preflight.json"
  $decisionPath = Join-Path $ProjectRoot "results/2026-09-24-disease-side-v3.5-g6-primary/G6_decision.json"
  if ($analysisExit -ne 0 -and (Test-Path -LiteralPath $preflightPath) -and (Test-Path -LiteralPath $decisionPath)) {
    $preflightResult = Get-Content -LiteralPath $preflightPath -Raw | ConvertFrom-Json
    $decisionResult = Get-Content -LiteralPath $decisionPath -Raw | ConvertFrom-Json
    $expectedPreflightStop = ($preflightResult.preflight_pass -eq $false -and
                              $preflightResult.effects_calculated -eq $false -and
                              $decisionResult.decision -eq "G6_NOT_EXECUTABLE" -and
                              $decisionResult.disease_effects_calculated -eq $false)
  }
}
finally {
  Pop-Location
}
if ($reportExit -ne 0) { throw "G6 report generation failed with R exit code $reportExit" }
if ($analysisExit -ne 0 -and -not $expectedPreflightStop) { throw "G6 primary analysis failed with R exit code $analysisExit; the report records the fail-closed state." }
if ($expectedPreflightStop) { Write-Output "Effect-blind gate failed as recorded; G6_NOT_EXECUTABLE report and QC are complete without group effects." }
