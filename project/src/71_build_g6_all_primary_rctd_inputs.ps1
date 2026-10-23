$ErrorActionPreference = 'Stop'
$projectRoot = (Get-Location).Path
$stageDir = 'C:\Users\administer\.codex\visualizations\2026\09\23\01a0ce11-eaa9-7cb1-aa30-49773cbedb70\g6_hdf5_stage'
$library = 'C:\Users\administer\.codex\visualizations\2026\09\23\01a0ce11-eaa9-7cb1-aa30-49773cbedb70\g6_R_library'
$rscript = 'D:\R\R-4.3.3\bin\Rscript.exe'
$manifestPath = 'results/2026-09-24-disease-side-v3.5-g6-input-audit/registries/GSE269906_tissue_CNGB_file_manifest.tsv'
$manifest = Import-Csv -LiteralPath $manifestPath -Delimiter "`t"
$outDir = 'data/interim/G6_GSE269906_bin50_full'
New-Item -ItemType Directory -Path $outDir -Force | Out-Null
$env:R_LIBS_USER = $library

foreach ($row in $manifest) {
    $outPath = Join-Path $outDir ($row.GSM + '_bin50_centroid_supported_counts.rds')
    if (Test-Path -LiteralPath $outPath) {
        Write-Host "Already extracted $($row.GSM); skipping."
        continue
    }
    $sourcePath = Join-Path 'data/raw/GSE269906_tissue' $row.donor_file
    $stagedPath = Join-Path $stageDir $row.donor_file
    if (-not (Test-Path -LiteralPath $sourcePath)) { throw "Missing raw source: $sourcePath" }
    Copy-Item -LiteralPath $sourcePath -Destination $stagedPath -Force
    try {
        & $rscript 'src/70_extract_g6_bin50_rctd_inputs.R' '.' $stagedPath $row.GSM
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed for $($row.GSM), exit=$LASTEXITCODE" }
    }
    finally {
        if (Test-Path -LiteralPath $stagedPath) { Remove-Item -LiteralPath $stagedPath -Force }
    }
}

$rows = Get-ChildItem -LiteralPath $outDir -Filter '*_bin50_centroid_supported_qc.csv' | ForEach-Object { Import-Csv -LiteralPath $_.FullName }
if ($rows.Count -ne 12) { throw "Expected 12 primary donor input matrices; found $($rows.Count)" }
$rows | Export-Csv -LiteralPath (Join-Path $outDir 'GSE269906_full_bin50_input_manifest.csv') -NoTypeInformation
Write-Host "All 12 primary donor matrices passed extraction checks."
