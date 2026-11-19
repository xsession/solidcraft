<#
.SYNOPSIS
    Pack every crates.io dependency of SolidCraft and its companion into one zip, so the
    workspace builds with no network (cargo --offline) on another machine.
.DESCRIPTION
    Runs `cargo vendor --locked` for the workspace and companion/cfd-solvers-rs, writes the
    source-replacement settings next to the crates, and zips both. Unzip at the repository root
    and copy vendor-crates\cargo-config.toml to .cargo\config.toml (or pass it with --config)
    to build offline. Nothing in the repository is changed.
.PARAMETER Out
    The zip to write (default: solidcraft-vendor.zip in the repository root).
#>
[CmdletBinding()]
param([string]$Out = "solidcraft-vendor.zip")
$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Push-Location $repo
try {
    $vendor = Join-Path $repo "vendor-crates"
    if (Test-Path $vendor) { Remove-Item -Recurse -Force $vendor }
    $config = cargo vendor --locked --versioned-dirs `
        --sync companion/cfd-solvers-rs/Cargo.toml vendor-crates
    if ($LASTEXITCODE -ne 0) { throw "cargo vendor failed" }
    $config | Out-File -Encoding utf8 (Join-Path $vendor "cargo-config.toml")
    $count = (Get-ChildItem $vendor -Directory).Count
    if (Test-Path $Out) { Remove-Item -Force $Out }
    Compress-Archive -Path $vendor -DestinationPath $Out -CompressionLevel Optimal
    $size = [math]::Round((Get-Item $Out).Length / 1MB, 1)
    Write-Host "Packed $count crates into $Out ($size MB)."
    Write-Host "If the zip is too large to upload in one piece, split it, e.g. with 7-Zip: 7z a -v95m solidcraft-vendor.7z vendor-crates"
}
finally {
    Pop-Location
}
