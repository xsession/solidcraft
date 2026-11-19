<#
.SYNOPSIS
    Build, test and lint everything the Multiphysics workspace changes touch, and keep one log.
.DESCRIPTION
    Runs, in order, and records every step's output and exit code in multiphysics-check.log:
      1. cargo fmt --check for the three crates
      2. cargo check for solidcraft-multiphysics, solidcraft-engine and solidcraft-ui-egui
         (all targets, so tests compile too)
      3. cargo clippy on the same crates (-D warnings, as CI)
      4. cargo test -p solidcraft-multiphysics
      5. cargo test -p solidcraft-engine multiphysics
      6. cargo test -p solidcraft-ui-egui multiphysics
      7. cargo test in companion/cfd-solvers-rs
    Later steps still run when one fails, so a single log shows every problem. Send the log
    back to have the errors fixed.
.PARAMETER Offline
    Pass --offline to cargo (after unpacking a vendor archive from offline-vendor.ps1).
#>
[CmdletBinding()]
param([switch]$Offline)
$repo = Split-Path -Parent $PSScriptRoot
$log = Join-Path $repo "multiphysics-check.log"
$net = @(); if ($Offline) { $net = @("--offline") }
$crates = @("-p", "solidcraft-multiphysics", "-p", "solidcraft-engine", "-p", "solidcraft-ui-egui")
"multiphysics-check $(Get-Date -Format o)" | Out-File -Encoding utf8 $log
"rustc: $(rustc --version)" | Out-File -Append -Encoding utf8 $log
$failed = @()
function Step([string]$name, [string]$dir, [string[]]$arguments) {
    "`n==== $name ====`n> cargo $($arguments -join ' ')" | Out-File -Append -Encoding utf8 $log
    Push-Location $dir
    try {
        & cargo @arguments 2>&1 | ForEach-Object { "$_" } | Out-File -Append -Encoding utf8 $log
        $code = $LASTEXITCODE
    }
    finally { Pop-Location }
    "exit code: $code" | Out-File -Append -Encoding utf8 $log
    Write-Host ("{0,-40} {1}" -f $name, $(if ($code -eq 0) { "ok" } else { "FAILED ($code)" }))
    if ($code -ne 0) { $script:failed += $name }
}
Step "format" $repo (@("fmt", "--check") + $crates)
Step "check (all targets)" $repo (@("check", "--all-targets") + $net + $crates)
Step "clippy" $repo (@("clippy", "--all-targets") + $net + $crates + @("--", "-D", "warnings"))
Step "test multiphysics crate" $repo (@("test") + $net + @("-p", "solidcraft-multiphysics"))
Step "test engine (multiphysics)" $repo (@("test") + $net + @("-p", "solidcraft-engine", "multiphysics"))
Step "test ui (multiphysics)" $repo (@("test") + $net + @("-p", "solidcraft-ui-egui", "multiphysics"))
Step "test companion" (Join-Path $repo "companion/cfd-solvers-rs") (@("test", "--locked") + $net)
if ($failed.Count -eq 0) { Write-Host "All steps passed. Log: $log" }
else { Write-Host "Failed: $($failed -join ', '). Send $log to have them fixed." }
