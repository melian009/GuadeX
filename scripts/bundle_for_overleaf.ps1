#requires -Version 7.0
<#
.SYNOPSIS
    Builds a small, self-contained folder that can be uploaded to Overleaf.

.DESCRIPTION
    Copies a main .tex file, every .tex it \input's or \include's, the .bib and
    .bst files it needs, and only the figures referenced by \includegraphics.
    Paths are preserved relative to the repository root so the existing
    relative paths (e.g. ../results/.../fig.png) still resolve on Overleaf.

.PARAMETER Main
    Main .tex file, relative to the repo root or absolute.
    Default: docs/FinalReportOctober2026.tex

.PARAMETER OutDir
    Output folder (relative to repo root or absolute).
    Default: overleaf_bundle

.PARAMETER Zip
    Also produce a .zip of the bundle, which is what you upload to Overleaf.

.PARAMETER Force
    Overwrite an existing bundle folder that was not created by this script.

.EXAMPLE
    pwsh scripts/bundle_for_overleaf.ps1 -Zip

.EXAMPLE
    pwsh scripts/bundle_for_overleaf.ps1 -Main docs/WorkingPaper/ABC3D.tex -OutDir overleaf_abc3d -Zip
#>
[CmdletBinding()]
param(
    [string]$Main = 'docs/FinalReportOctober2026.tex',
    [string]$OutDir = 'overleaf_bundle',
    [switch]$Zip,
    [switch]$Force
)

$ErrorActionPreference = 'Stop'

$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$Marker   = '.generated-by-bundle_for_overleaf'

function Resolve-RepoPath([string]$Path) {
    if ([System.IO.Path]::IsPathRooted($Path)) { return $Path }
    Join-Path $RepoRoot $Path
}

$MainPath = Resolve-RepoPath $Main
if (-not (Test-Path -LiteralPath $MainPath)) {
    throw "Main .tex file not found: $MainPath"
}
$MainPath = (Resolve-Path -LiteralPath $MainPath).Path

$OutDir = Resolve-RepoPath $OutDir
if (Test-Path -LiteralPath $OutDir) {
    if (Test-Path -LiteralPath (Join-Path $OutDir $Marker)) {
        Remove-Item -LiteralPath $OutDir -Recurse -Force
    }
    elseif ($Force) {
        Remove-Item -LiteralPath $OutDir -Recurse -Force
    }
    else {
        throw "Output folder already exists and was not created by this script: $OutDir (use -Force to overwrite)"
    }
}
New-Item -ItemType Directory -Path $OutDir -Force | Out-Null
Set-Content -LiteralPath (Join-Path $OutDir $Marker) -Value "Generated $(Get-Date -Format o)"

$script:Copied   = [System.Collections.Generic.List[object]]::new()
$script:Seen     = [System.Collections.Generic.HashSet[string]]::new([StringComparer]::OrdinalIgnoreCase)
$script:Warnings = [System.Collections.Generic.List[string]]::new()

function Get-DisplayPath([string]$FullPath) {
    [System.IO.Path]::GetRelativePath($RepoRoot, $FullPath).Replace('\', '/')
}

function Copy-IntoBundle([string]$FullPath, [string]$Reason) {
    if (-not $script:Seen.Add($FullPath)) { return }
    $rel     = Get-DisplayPath $FullPath
    $dest    = Join-Path $OutDir $rel
    $destDir = Split-Path -Parent $dest
    New-Item -ItemType Directory -Path $destDir -Force | Out-Null
    Copy-Item -LiteralPath $FullPath -Destination $dest -Force
    $script:Copied.Add([pscustomobject]@{
        Reason = $Reason
        Path   = $rel
        KB     = [math]::Round((Get-Item -LiteralPath $FullPath).Length / 1KB, 1)
    })
}

function Resolve-WithExtensions([string]$Spec, [string]$BaseDir, [string[]]$Extensions) {
    $candidates = [System.Collections.Generic.List[string]]::new()
    $candidates.Add($Spec)
    if ([string]::IsNullOrEmpty([System.IO.Path]::GetExtension($Spec))) {
        foreach ($ext in $Extensions) { $candidates.Add("$Spec$ext") }
    }
    foreach ($candidate in $candidates) {
        $p = if ([System.IO.Path]::IsPathRooted($candidate)) { $candidate } else { Join-Path $BaseDir $candidate }
        if (Test-Path -LiteralPath $p -PathType Leaf) { return (Resolve-Path -LiteralPath $p).Path }
    }
    return $null
}

$queue = [System.Collections.Generic.Queue[object]]::new()
$queue.Enqueue([pscustomobject]@{ Path = $MainPath; Reason = 'main document' })

while ($queue.Count -gt 0) {
    $item    = $queue.Dequeue()
    $texPath = $item.Path
    $docDir  = Split-Path -Parent $texPath
    $text    = Get-Content -LiteralPath $texPath -Raw

    Copy-IntoBundle $texPath $item.Reason

    foreach ($m in [regex]::Matches($text, '\\includegraphics\*?(?:\[[^\]]*\])?\{([^}]+)\}')) {
        $specs = $m.Groups[1].Value.Split(',') | ForEach-Object { $_.Trim().Trim('"') } | Where-Object { $_ }
        foreach ($spec in $specs) {
            $found = Resolve-WithExtensions $spec $docDir @('.pdf', '.png', '.jpg', '.jpeg', '.eps')
            if ($found) {
                Copy-IntoBundle $found 'figure'
                if ([System.IO.Path]::GetExtension($found) -ieq '.eps') {
                    $script:Warnings.Add("EPS figure included: $(Get-DisplayPath $found). Overleaf's pdfLaTeX cannot use EPS; switch the project to LaTeX or convert it to PDF.")
                }
            }
            else {
                $script:Warnings.Add("Missing figure '$spec' referenced from $(Get-DisplayPath $texPath)")
            }
        }
    }

    foreach ($m in [regex]::Matches($text, '\\(?:input|include)\{([^}]+)\}')) {
        $spec = $m.Groups[1].Value.Trim().Trim('"')
        $found = Resolve-WithExtensions $spec $docDir @('.tex')
        if ($found) {
            $queue.Enqueue([pscustomobject]@{ Path = $found; Reason = 'included .tex' })
        }
        else {
            $script:Warnings.Add("Missing included file '$spec' referenced from $(Get-DisplayPath $texPath)")
        }
    }

    foreach ($m in [regex]::Matches($text, '\\bibliography\{([^}]+)\}')) {
        foreach ($spec in ($m.Groups[1].Value.Split(',') | ForEach-Object { $_.Trim().Trim('"') } | Where-Object { $_ })) {
            $found = Resolve-WithExtensions $spec $docDir @('.bib')
            if ($found) {
                Copy-IntoBundle $found 'bibliography'
            }
            else {
                $script:Warnings.Add("Missing bibliography '$spec' referenced from $(Get-DisplayPath $texPath)")
            }
        }
    }

    foreach ($m in [regex]::Matches($text, '\\bibliographystyle\{([^}]+)\}')) {
        $spec = $m.Groups[1].Value.Trim().Trim('"')
        $found = Resolve-WithExtensions $spec $docDir @('.bst')
        if ($found) {
            Copy-IntoBundle $found 'bibliography style'
        }
        elseif ($spec -notin @('plain', 'unsrt', 'alpha', 'abbrv')) {
            $script:Warnings.Add("Missing .bst '$spec' referenced from $(Get-DisplayPath $texPath) (Overleaf may have it built in)")
        }
    }
}

$totalKB = ($script:Copied | Measure-Object -Property KB -Sum).Sum
$totalMB = [math]::Round($totalKB / 1KB, 2)

Write-Output ''
Write-Output ("Bundle: {0}" -f (Get-DisplayPath $OutDir))
Write-Output ("Main document: {0}" -f (Get-DisplayPath $MainPath))
Write-Output ("Files: {0}   Total: {1} MB" -f $script:Copied.Count, $totalMB)
Write-Output ''
$script:Copied | Sort-Object Reason, Path | Format-Table -AutoSize | Out-String | Write-Output

foreach ($w in $script:Warnings) { Write-Warning $w }

if ($totalMB -gt 50) {
    Write-Warning ("Bundle is {0} MB, above the 50 MB Overleaf limit. Compress the figures (e.g. PNG -> JPEG/WebP) before uploading." -f $totalMB)
}

if ($Zip) {
    $zipPath = "$OutDir.zip"
    if (Test-Path -LiteralPath $zipPath) { Remove-Item -LiteralPath $zipPath -Force }
    Compress-Archive -Path (Join-Path $OutDir '*') -DestinationPath $zipPath -Force
    Write-Output ("Zip: {0} ({1} MB)" -f (Get-DisplayPath $zipPath), [math]::Round((Get-Item -LiteralPath $zipPath).Length / 1MB, 2))
    Write-Output 'Upload this zip to Overleaf, then set the main document to the path shown above.'
}
