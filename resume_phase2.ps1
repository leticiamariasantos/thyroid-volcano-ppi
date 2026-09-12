param(
  [int]$From = 25,
  [switch]$Force
)
# Retomada do pipeline Phase 2 + upgrade 2026.
# IMPORTANTE: rodar via PowerShell (não via bash MSYS2), pois Rscript chamado
# diretamente do MSYS2 causa segfault em fwrite/saveRDS (data.table 1.18.4/R 4.6.1).
# Uso: powershell -NoProfile -ExecutionPolicy Bypass -File resume_phase2.ps1 -From 25
$ErrorActionPreference = "Continue"
$Rscript = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$stdout = "logs/phase2_resume_${stamp}_stdout.log"
$stderr = "logs/phase2_resume_${stamp}_stderr.log"
$argList = @("scripts/run_phase2.R", "--from=$From")
if ($Force) { $argList += "--force" }
$p = Start-Process -FilePath $Rscript -ArgumentList $argList `
  -RedirectStandardOutput $stdout -RedirectStandardError $stderr `
  -WindowStyle Hidden -PassThru
Set-Content -Path "logs/last_run_pid.txt" -Value $p.Id
Write-Output "PID=$($p.Id)"
Write-Output "STDOUT=$stdout"
Write-Output "STDERR=$stderr"
Write-Output "Para monitorar: tail -c 2000 $stdout"
