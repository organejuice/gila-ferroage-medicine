param(
  [string]$ProjectRoot = ".",
  [string]$Rscript = "D:\R\R-4.3.3\bin\Rscript.exe",
  [string]$RLibrary = "C:\Users\administer\.codex\visualizations\2026\09\23\01a0ce11-eaa9-7cb1-aa30-49773cbedb70\g6_R_library"
)

$ErrorActionPreference = "Stop"
$ProjectRoot = (Resolve-Path -LiteralPath $ProjectRoot).Path
$Rscript = (Resolve-Path -LiteralPath $Rscript).Path
$env:R_LIBS_USER = $RLibrary
$auditRoot = Join-Path $ProjectRoot "results\2026-09-24-disease-side-v3.5-g6-input-audit"
$manifest = Import-Csv -LiteralPath (Join-Path $auditRoot "registries\GSE269906_donor_section_manifest.tsv") -Delimiter "`t"
$countInputs = Join-Path $ProjectRoot "data\interim\G6_GSE269906_bin50_full"
$weightsRoot = Join-Path $ProjectRoot "data\interim\G6_GSE269906_RCTD"
$logRoot = Join-Path $auditRoot "logs"
New-Item -ItemType Directory -Path $logRoot -Force | Out-Null

function Get-AvailableMemoryMB {
  try {
    $counter = Get-Counter "\Memory\Available MBytes" -ErrorAction Stop
    return [int][math]::Floor($counter.CounterSamples[0].CookedValue)
  }
  catch {
    return 0
  }
}

$work = foreach ($donor in $manifest) {
  $gsm = $donor.GSM
  $countPath = Join-Path $countInputs "${gsm}_bin50_centroid_supported_counts.rds"
  $qcPath = Join-Path $countInputs "${gsm}_bin50_centroid_supported_qc.csv"
  $weightPath = Join-Path $weightsRoot "${gsm}_RCTD_weights_effect_blind.rds"
  $sumPath = Join-Path $weightsRoot "${gsm}_RCTD_weight_sum_qc.csv"
  $markerPath = Join-Path $weightsRoot "${gsm}_RCTD_holdout_marker_qc.csv"
  if (!(Test-Path -LiteralPath $countPath) -or !(Test-Path -LiteralPath $qcPath)) {
    throw "Missing locked full-count input or QC table for $gsm"
  }
  $nBins = [int](Import-Csv -LiteralPath $qcPath).n_bins
  [pscustomobject]@{
    GSM = $gsm
    nBins = $nBins
    complete = (Test-Path -LiteralPath $weightPath) -and (Test-Path -LiteralPath $sumPath) -and (Test-Path -LiteralPath $markerPath)
  }
}

$todo = @($work | Sort-Object nBins, GSM | Where-Object { !$_.complete })
while ($todo.Count -gt 0) {
  $availableMB = Get-AvailableMemoryMB
  $batchSize = if ($availableMB -ge 6000) { 2 } else { 1 }
  $batch = @($todo | Select-Object -First $batchSize)
  $todo = @($todo | Select-Object -Skip $batch.Count)
  Write-Host "Available memory=${availableMB}MB; starting batch of $($batch.Count) donor fit(s)."

  $jobs = @()
  foreach ($item in $batch) {
    if ($jobs.Count -gt 0 -and (Get-AvailableMemoryMB) -lt 4500) {
      $todo = @($item) + $todo
      Write-Host "Memory guard deferred $($item.GSM); remaining fits will run after this batch."
      break
    }
    $stdoutPath = Join-Path $logRoot "$($item.GSM)_RCTD.log"
    $stderrPath = Join-Path $logRoot "$($item.GSM)_RCTD.stderr.log"
    Write-Host "START effect-blind RCTD: $($item.GSM) ($($item.nBins) bins); stdout=$stdoutPath; stderr=$stderrPath"
    $process = Start-Process -FilePath $Rscript `
      -ArgumentList @("--vanilla", "src/72_fit_g6_rctd_effect_blind.R", ".", $item.GSM) `
      -WorkingDirectory $ProjectRoot -WindowStyle Hidden -PassThru `
      -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath
    $jobs += [pscustomobject]@{ Item = $item; Process = $process; Stdout = $stdoutPath; Stderr = $stderrPath }
  }

  foreach ($job in $jobs) { $job.Process.WaitForExit() }
  foreach ($job in $jobs) {
    $job.Process.Refresh()
    if ($job.Process.ExitCode -ne 0) {
      throw "RCTD failed for $($job.Item.GSM) with exit code $($job.Process.ExitCode). See $($job.Stdout) and $($job.Stderr)"
    }
    foreach ($requiredPath in @(
      (Join-Path $weightsRoot "$($job.Item.GSM)_RCTD_weights_effect_blind.rds"),
      (Join-Path $weightsRoot "$($job.Item.GSM)_RCTD_weight_sum_qc.csv"),
      (Join-Path $weightsRoot "$($job.Item.GSM)_RCTD_holdout_marker_qc.csv")
    )) {
      if (!(Test-Path -LiteralPath $requiredPath)) { throw "RCTD exited successfully but did not produce $requiredPath" }
    }
    Write-Host "COMPLETE effect-blind RCTD: $($job.Item.GSM)"
  }
}

Write-Host "All 12 primary donors have complete effect-blind RCTD outputs."
