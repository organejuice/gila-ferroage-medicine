param(
    [string]$ProjectRoot = "."
)

$ErrorActionPreference = "Stop"
$metaDir = Join-Path $ProjectRoot "data/raw/G6_source_metadata_2026-09-23"
$targetDir = Join-Path $ProjectRoot "data/interim/GSE269906_tissue_reacquire"
$outputDir = Join-Path $ProjectRoot "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries"
$fileListPath = Join-Path $metaDir "GSE269906_filelist.txt"
$chunkSize = 33554432L
$maxParallel = 2
New-Item -ItemType Directory -Force -Path $targetDir, $outputDir | Out-Null

$geoRows = Get-Content -LiteralPath $fileListPath | Select-Object -Skip 1 |
    ConvertFrom-Csv -Delimiter "`t" -Header Kind, Name, Time, Size, Type
$items = foreach ($row in $geoRows) {
    if ($row.Kind -ne "File" -or $row.Name -notmatch "\.tissue\.gef$") { continue }
    $gsm = [regex]::Match($row.Name, "^GSM\d+").Value
    [pscustomobject]@{
        Name = $row.Name
        GSM = $gsm
        Expected = [long]($row.Size -replace ",", "")
        Url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8330nnn/$gsm/suppl/$($row.Name)"
    }
}
if (@($items).Count -ne 12) { throw "Expected 12 tissue GEF files in the official GEO file listing; found $(@($items).Count)." }
$totalBytes = [long](($items | Measure-Object -Property Expected -Sum).Sum)
$freeBytes = (Get-PSDrive -Name ([System.IO.Path]::GetPathRoot((Resolve-Path $ProjectRoot).Path).TrimEnd('\').Substring(0,1))).Free
if ($freeBytes -lt ($totalBytes * 1.5)) { throw "Insufficient free space for tissue GEF staging: need at least 1.5x $totalBytes bytes." }
Write-Output "GEO_TISSUE_TOTAL_BYTES $totalBytes"

$partManifest = [System.Collections.Generic.List[object]]::new()
foreach ($item in $items) {
    $assembled = Join-Path $targetDir ($item.Name + ".reassembled")
    if ((Test-Path -LiteralPath $assembled) -and (Get-Item -LiteralPath $assembled).Length -eq $item.Expected) {
        $partManifest.Add([pscustomobject]@{ GSM=$item.GSM; Name=$item.Name; Start=0L; End=$item.Expected-1L; Bytes=$item.Expected; Status="existing_exact_reassembled" })
        continue
    }

    $plan = [System.Collections.Generic.List[object]]::new()
    for ($start = 0L; $start -lt $item.Expected; $start += $chunkSize) {
        $end = [Math]::Min($item.Expected - 1L, $start + $chunkSize - 1L)
        $expectedPart = $end - $start + 1L
        $partPath = Join-Path $targetDir ("$($item.GSM).range.$start-$end.part")
        if ((Test-Path -LiteralPath $partPath) -and (Get-Item -LiteralPath $partPath).Length -eq $expectedPart) {
            $plan.Add([pscustomobject]@{ GSM=$item.GSM; Name=$item.Name; Url=$item.Url; Start=$start; End=$end; Bytes=$expectedPart; Part=$partPath; Status="reused_exact_length" })
        } else {
            if (Test-Path -LiteralPath $partPath) { Remove-Item -LiteralPath $partPath }
            $plan.Add([pscustomobject]@{ GSM=$item.GSM; Name=$item.Name; Url=$item.Url; Start=$start; End=$end; Bytes=$expectedPart; Part=$partPath; Status="pending" })
        }
    }

    $pending = @($plan | Where-Object { $_.Status -eq "pending" })
    for ($i = 0; $i -lt $pending.Count; $i += 8) {
        $batchEnd = [Math]::Min($pending.Count - 1, $i + 7)
        $batch = @($pending[$i..$batchEnd])
        $curlArgs = [System.Collections.Generic.List[string]]::new()
        $curlArgs.Add("--fail"); $curlArgs.Add("--location"); $curlArgs.Add("--silent"); $curlArgs.Add("--show-error")
        $curlArgs.Add("--connect-timeout"); $curlArgs.Add("30"); $curlArgs.Add("--max-time"); $curlArgs.Add("1800")
        $curlArgs.Add("--parallel"); $curlArgs.Add("--parallel-immediate"); $curlArgs.Add("--parallel-max"); $curlArgs.Add([string]$maxParallel)
        for ($j = 0; $j -lt $batch.Count; $j++) {
            if ($j -gt 0) { $curlArgs.Add("--next") }
            $curlArgs.Add("--range"); $curlArgs.Add("$($batch[$j].Start)-$($batch[$j].End)")
            $curlArgs.Add("--output"); $curlArgs.Add($batch[$j].Part)
            $curlArgs.Add($batch[$j].Url)
        }
        Write-Output "GEO_TISSUE_BATCH $($item.GSM) parts=$($i+1)-$($batchEnd+1) of $($pending.Count)"
        & curl.exe @curlArgs
        $curlExit = $LASTEXITCODE
        foreach ($request in $batch) {
            $actual = if (Test-Path -LiteralPath $request.Part) { (Get-Item -LiteralPath $request.Part).Length } else { 0L }
            if ($actual -eq $request.Bytes) {
                $request.Status = "downloaded_exact_length"
            } else {
                if (Test-Path -LiteralPath $request.Part) { Remove-Item -LiteralPath $request.Part }
                $succeeded = $false
                for ($attempt = 1; $attempt -le 3 -and -not $succeeded; $attempt++) {
                    Write-Output "GEO_TISSUE_RETRY $($item.GSM) $($request.Start)-$($request.End) attempt=$attempt curl_exit=$curlExit"
                    $response = & curl.exe --fail --location --silent --show-error --connect-timeout 30 --max-time 1800 --range "$($request.Start)-$($request.End)" --output $request.Part --write-out "%{http_code}|%{size_download}" $request.Url
                    $oneExit = $LASTEXITCODE
                    $oneBytes = if (Test-Path -LiteralPath $request.Part) { (Get-Item -LiteralPath $request.Part).Length } else { 0L }
                    if ($oneExit -eq 0 -and "$response" -eq "206|$($request.Bytes)" -and $oneBytes -eq $request.Bytes) { $succeeded = $true }
                    elseif (Test-Path -LiteralPath $request.Part) { Remove-Item -LiteralPath $request.Part }
                }
                if (-not $succeeded) { throw "Failed to retrieve exact official GEO range $($request.Start)-$($request.End) for $($item.Name)." }
                $request.Status = "retried_206_exact_length"
            }
        }
    }

    $assemblePart = Join-Path $targetDir ($item.Name + ".reassembled.partial")
    $output = [System.IO.File]::Open($assemblePart, [System.IO.FileMode]::Create, [System.IO.FileAccess]::Write, [System.IO.FileShare]::None)
    try {
        foreach ($request in ($plan | Sort-Object Start)) {
            $input = [System.IO.File]::OpenRead($request.Part)
            try { $input.CopyTo($output) } finally { $input.Dispose() }
        }
    } finally { $output.Dispose() }
    $actualTotal = (Get-Item -LiteralPath $assemblePart).Length
    if ($actualTotal -ne $item.Expected) { throw "Reassembled size mismatch for $($item.Name): $actualTotal vs $($item.Expected)." }
    Move-Item -LiteralPath $assemblePart -Destination $assembled -Force
    foreach ($request in $plan) { $partManifest.Add($request) }
    $sha = (Get-FileHash -Algorithm SHA256 -LiteralPath $assembled).Hash
    [pscustomobject]@{ file=$item.Name; expected_bytes=$item.Expected; actual_bytes=$actualTotal; sha256=$sha } | ConvertTo-Json -Compress
}

$partManifest | Select-Object GSM,Name,Start,End,Bytes,Status | Export-Csv -LiteralPath (Join-Path $outputDir "GSE269906_tissue_GEO_range_retrieval_manifest.tsv") -Delimiter "`t" -NoTypeInformation -Encoding utf8
