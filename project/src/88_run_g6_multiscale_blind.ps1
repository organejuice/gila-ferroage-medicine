param(
  [ValidateSet('Build','Fit','Summarize')][string]$Phase,
  [string]$Rscript = 'D:\R\R-4.3.3\bin\Rscript.exe',
  [string]$RLibrary = 'C:\Users\administer\.codex\visualizations\2026\09\23\01a0ce11-eaa9-7cb1-aa30-49773cbedb70\g6_R_library'
)
$ErrorActionPreference = 'Stop'
$root = (Get-Location).Path
if (-not (Test-Path -LiteralPath (Join-Path $root 'configs\g6_multiscale_pilot_001.yaml'))) {
  throw 'Run from the project root.'
}
if (-not (Test-Path -LiteralPath $Rscript)) { throw 'Rscript is missing.' }
$env:R_LIBS_USER = $RLibrary
$r = (Resolve-Path -LiteralPath $Rscript).Path
if ($Phase -eq 'Build') {
  & $r 'src/87_prepare_g6_multiscale_conservative_reference.R' .
  if ($LASTEXITCODE -ne 0) { throw 'Conservative reference preparation failed.' }
  foreach ($width in @(50, 110)) {
    & $r 'src/83_build_g6_multiscale_blind.R' . $width
    if ($LASTEXITCODE -ne 0) { throw "Raw bin$width reconstruction failed." }
  }
  return
}
if ($Phase -eq 'Fit') {
  $manifestPath = Join-Path $root 'results\2026-09-24-disease-side-v3.5-g6-input-audit\registries\GSE269906_tissue_CNGB_file_manifest.tsv'
  $manifest = Import-Csv -LiteralPath $manifestPath -Delimiter "`t" | Select-Object GSM, donor_file
  if ($manifest.Count -ne 12) { throw 'Expected 12 tissue GEF records.' }
  $stageDir = Join-Path $env:TEMP ('g6_ms_' + [guid]::NewGuid().ToString('N'))
  New-Item -ItemType Directory -Path $stageDir | Out-Null
  try {
    foreach ($width in @(50, 110)) {
      foreach ($row in $manifest) {
        $signature = Join-Path $root ("data\interim\G6_GSE269906_multiscale_001\bin{0}\signature\{1}_bin{0}_tissue_signature_counts.rds" -f $width, $row.GSM)
        if (-not (Test-Path -LiteralPath $signature)) { throw "Missing scale reconstruction: $signature" }
        $raw = Join-Path $root (Join-Path 'data\raw\GSE269906_tissue' $row.donor_file)
        $stage = Join-Path $stageDir $row.donor_file
        Copy-Item -LiteralPath $raw -Destination $stage -ErrorAction Stop
        try {
          & $r 'src/84_extract_g6_multiscale_rctd_inputs.R' . $stage $row.GSM $width
          if ($LASTEXITCODE -ne 0) { throw "Full-count extraction failed: $($row.GSM) bin$width" }
        } finally {
          if (Test-Path -LiteralPath $stage) { Remove-Item -LiteralPath $stage -Force }
        }
        & $r 'src/85_fit_g6_multiscale_rctd_blind.R' . $row.GSM $width
        if ($LASTEXITCODE -ne 0) { throw "RCTD failed: $($row.GSM) bin$width" }
      }
    }
  } finally {
    if (Test-Path -LiteralPath $stageDir) { Remove-Item -LiteralPath $stageDir }
  }
  return
}
& $r 'src/86_summarize_g6_multiscale_blind.R' .
if ($LASTEXITCODE -ne 0) { throw 'Blind summary failed.' }
