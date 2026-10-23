param(
    [string]$ProjectRoot = "."
)

$ErrorActionPreference = "Stop"
$metaDir = Join-Path $ProjectRoot "data/raw/G6_source_metadata_2026-09-23"
$targetDir = Join-Path $ProjectRoot "data/interim/GSE269906_tissue_cngb_reacquire"
$outputDir = Join-Path $ProjectRoot "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries"
$fileListPath = Join-Path $metaDir "GSE269906_filelist.txt"
$chunkSize = 8388608L
$maxParallel = 8
$batchSize = 32
New-Item -ItemType Directory -Force -Path $targetDir, $outputDir | Out-Null

$geoRows = Get-Content -LiteralPath $fileListPath | Select-Object -Skip 1 |
    ConvertFrom-Csv -Delimiter "`t" -Header Kind, Name, Time, Size, Type
$items = @(
    foreach ($row in $geoRows | Where-Object { $_.Kind -eq "File" -and $_.Name -match "\.tissue\.gef$" } | Sort-Object Name) {
        $gsm = [regex]::Match($row.Name, "^GSM\d+").Value
        $gsmNumber = [int]$gsm.Substring(3)
        $sampleOffset = $gsmNumber - 8330060
        if ($sampleOffset -lt 0 -or $sampleOffset -gt 11) { throw "Unexpected GSM in tissue file listing: $gsm" }
        $sampleId = "STSP{0:D7}" -f (10540 + $sampleOffset)
        [pscustomobject]@{
            Name = $row.Name
            GSM = $gsm
            STOmicsSample = $sampleId
            Expected = [long]($row.Size -replace ",", "")
            Url = "ftp://ftp.cngb.org/pub/SciRAID/stomics/STDS0000318/stomics/$sampleId/$gsm/$($row.Name)"
            DirectoryUrl = "ftp://ftp.cngb.org/pub/SciRAID/stomics/STDS0000318/stomics/$sampleId/$gsm/"
        }
    }
)
if ($items.Count -ne 12) { throw "Expected 12 tissue GEF files in official GEO listing; found $($items.Count)." }
$totalBytes = [long](($items | Measure-Object -Property Expected -Sum).Sum)
$driveName = [System.IO.Path]::GetPathRoot((Resolve-Path $ProjectRoot).Path).TrimEnd('\').Substring(0,1)
$freeBytes = (Get-PSDrive -Name $driveName).Free
if ($freeBytes -lt ($totalBytes * 1.5)) { throw "Insufficient free space for tissue GEF staging: need at least 1.5x $totalBytes bytes." }
Write-Output "CNGB_TISSUE_TOTAL_BYTES $totalBytes"

# Validate the public mirror's GSM/sample/file crosswalk before transferring.
foreach ($item in $items) {
    $listing = & curl.exe --fail --silent --show-error --connect-timeout 15 --max-time 45 --list-only $item.DirectoryUrl
    if ($LASTEXITCODE -ne 0 -or $listing -notcontains $item.Name) {
        throw "The STOmics mirror does not map $($item.GSM) to $($item.STOmicsSample)/$($item.Name)."
    }
}
Write-Output "CNGB_SAMPLE_FILE_CROSSWALK_VERIFIED $($items.Count)"

$partManifest = [System.Collections.Generic.List[object]]::new()
$fileManifest = [System.Collections.Generic.List[object]]::new()
foreach ($item in $items) {
    $assembled = Join-Path $targetDir ($item.Name + ".reassembled")
    $plan = [System.Collections.Generic.List[object]]::new()
    if ((Test-Path -LiteralPath $assembled) -and (Get-Item -LiteralPath $assembled).Length -eq $item.Expected) {
        $sha = (Get-FileHash -Algorithm SHA256 -LiteralPath $assembled).Hash.ToLowerInvariant()
        $fileManifest.Add([pscustomobject]@{
            accession="GSE269906"; GSM=$item.GSM; donor_file=$item.Name; sample_id=$item.STOmicsSample
            GEO_expected_bytes=$item.Expected; CNGB_actual_bytes=$item.Expected; sha256=$sha
            source_database="NCBI GEO mirrored by CNGB STOmicsDB STDS0000318"
            GEO_source_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8330nnn/$($item.GSM)/suppl/$($item.Name)"
            mirror_url=$item.Url; status="existing_exact_reassembled"
        })
        continue
    }
    if (Test-Path -LiteralPath $assembled) { Remove-Item -LiteralPath $assembled }

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
    for ($i = 0; $i -lt $pending.Count; $i += $batchSize) {
        $batchEnd = [Math]::Min($pending.Count - 1, $i + $batchSize - 1)
        $batch = @($pending[$i..$batchEnd])
        $curlArgs = [System.Collections.Generic.List[string]]::new()
        $curlArgs.Add("--fail"); $curlArgs.Add("--silent"); $curlArgs.Add("--show-error")
        $curlArgs.Add("--connect-timeout"); $curlArgs.Add("30"); $curlArgs.Add("--max-time"); $curlArgs.Add("1800")
        $curlArgs.Add("--parallel"); $curlArgs.Add("--parallel-immediate"); $curlArgs.Add("--parallel-max"); $curlArgs.Add([string]$maxParallel)
        for ($j = 0; $j -lt $batch.Count; $j++) {
            if ($j -gt 0) { $curlArgs.Add("--next") }
            $curlArgs.Add("--range"); $curlArgs.Add("$($batch[$j].Start)-$($batch[$j].End)")
            $curlArgs.Add("--output"); $curlArgs.Add($batch[$j].Part)
            $curlArgs.Add($batch[$j].Url)
        }
        Write-Output "CNGB_TISSUE_BATCH $($item.GSM) parts=$($i+1)-$($batchEnd+1) of $($pending.Count)"
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
                    Write-Output "CNGB_TISSUE_RETRY $($item.GSM) $($request.Start)-$($request.End) attempt=$attempt curl_exit=$curlExit"
                    $response = & curl.exe --fail --silent --show-error --connect-timeout 30 --max-time 1800 --range "$($request.Start)-$($request.End)" --output $request.Part --write-out "%{size_download}" $request.Url
                    $oneExit = $LASTEXITCODE
                    $oneBytes = if (Test-Path -LiteralPath $request.Part) { (Get-Item -LiteralPath $request.Part).Length } else { 0L }
                    if ($oneExit -eq 0 -and "$response" -eq "$($request.Bytes)" -and $oneBytes -eq $request.Bytes) { $succeeded = $true }
                    elseif (Test-Path -LiteralPath $request.Part) { Remove-Item -LiteralPath $request.Part }
                }
                if (-not $succeeded) { throw "Failed to retrieve exact STOmics/CNGB range $($request.Start)-$($request.End) for $($item.Name)." }
                $request.Status = "retried_exact_length"
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
    foreach ($request in $plan) {
        $partManifest.Add([pscustomobject]@{
            accession="GSE269906"; GSM=$item.GSM; file=$item.Name; STOmics_sample_id=$item.STOmicsSample
            start=$request.Start; end=$request.End; bytes=$request.Bytes; status=$request.Status
            source_url=$request.Url; source_database="CNGB STOmicsDB STDS0000318"
        })
    }
    $sha = (Get-FileHash -Algorithm SHA256 -LiteralPath $assembled).Hash.ToLowerInvariant()
    $fileManifest.Add([pscustomobject]@{
        accession="GSE269906"; GSM=$item.GSM; donor_file=$item.Name; sample_id=$item.STOmicsSample
        GEO_expected_bytes=$item.Expected; CNGB_actual_bytes=$actualTotal; sha256=$sha
        source_database="NCBI GEO mirrored by CNGB STOmicsDB STDS0000318"
        GEO_source_url="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8330nnn/$($item.GSM)/suppl/$($item.Name)"
        mirror_url=$item.Url; status="complete_exact_length"
    })
    [pscustomobject]@{ file=$item.Name; sample_id=$item.STOmicsSample; expected_bytes=$item.Expected; actual_bytes=$actualTotal; sha256=$sha } | ConvertTo-Json -Compress
}

$partManifest | Export-Csv -LiteralPath (Join-Path $outputDir "GSE269906_tissue_CNGB_range_retrieval_manifest.tsv") -Delimiter "`t" -NoTypeInformation -Encoding utf8
$fileManifest | Export-Csv -LiteralPath (Join-Path $outputDir "GSE269906_tissue_CNGB_file_manifest.tsv") -Delimiter "`t" -NoTypeInformation -Encoding utf8
