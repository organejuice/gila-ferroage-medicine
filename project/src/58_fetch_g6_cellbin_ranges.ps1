param(
    [string]$ProjectRoot = "."
)

$ErrorActionPreference = "Stop"
$metaDir = Join-Path $ProjectRoot "data/raw/G6_source_metadata_2026-09-23"
$targetDir = Join-Path $ProjectRoot "data/interim/GSE269906_cellbin_reacquire"
$outputDir = Join-Path $ProjectRoot "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries"
$fileListPath = Join-Path $metaDir "GSE269906_filelist.txt"
$chunkSize = 33554432L
$maxParallel = 2
New-Item -ItemType Directory -Force -Path $targetDir, $outputDir | Out-Null

$geoRows = Get-Content -LiteralPath $fileListPath | Select-Object -Skip 1 |
    ConvertFrom-Csv -Delimiter "`t" -Header Kind, Name, Time, Size, Type
$requested = @(
    "GSM8330067_B01809C2.cellbin.gef",
    "GSM8330068_B01809A4.cellbin.gef",
    "GSM8330070_D02175A6.cellbin.gef"
)
$items = foreach ($name in $requested) {
    $row = $geoRows | Where-Object { $_.Kind -eq "File" -and $_.Name -eq $name } | Select-Object -First 1
    if (-not $row) { throw "File is absent from the official GEO listing: $name" }
    [pscustomobject]@{
        Name = $name
        GSM = [regex]::Match($name, "^GSM\d+").Value
        Expected = [long]($row.Size -replace ",", "")
        Url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8330nnn/$([regex]::Match($name, '^GSM\d+').Value)/suppl/$name"
    }
}

$partManifest = [System.Collections.Generic.List[object]]::new()
foreach ($item in $items) {
    $assembled = Join-Path $targetDir ($item.Name + ".reassembled")
    if ((Test-Path -LiteralPath $assembled) -and (Get-Item -LiteralPath $assembled).Length -eq $item.Expected) {
        $partManifest.Add([pscustomobject]@{ GSM=$item.GSM; Name=$item.Name; Start=0L; End=$item.Expected-1L; Bytes=$item.Expected; Status="existing_exact_reassembled" })
        continue
    }

    $basePath = Join-Path $targetDir $item.Name
    $baseLength = 0L
    if (Test-Path -LiteralPath $basePath) { $baseLength = (Get-Item -LiteralPath $basePath).Length }
    if ($baseLength -ge $item.Expected) { $baseLength = 0L }
    $plan = [System.Collections.Generic.List[object]]::new()
    for ($start = $baseLength; $start -lt $item.Expected; $start += $chunkSize) {
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
    for ($i = 0; $i -lt $pending.Count; $i += $maxParallel) {
        $batchEnd = [Math]::Min($pending.Count - 1, $i + $maxParallel - 1)
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
        Write-Output "GEO_RANGE_BATCH $($item.GSM) $($batch[0].Start)-$($batch[-1].End) requests=$($batch.Count)"
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
                    Write-Output "GEO_RANGE_RETRY $($item.GSM) $($request.Start)-$($request.End) attempt=$attempt curl_exit=$curlExit"
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
        if ($baseLength -gt 0) {
            $input = [System.IO.File]::OpenRead($basePath)
            try { $input.CopyTo($output) } finally { $input.Dispose() }
        }
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

$partManifest | Select-Object GSM,Name,Start,End,Bytes,Status | Export-Csv -LiteralPath (Join-Path $outputDir "GSE269906_GEO_range_retrieval_manifest.tsv") -Delimiter "`t" -NoTypeInformation -Encoding utf8
