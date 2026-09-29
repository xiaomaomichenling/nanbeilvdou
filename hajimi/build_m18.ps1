$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
if (-not (Get-Command g++ -ErrorAction SilentlyContinue)) {
    $env:PATH = [Environment]::GetEnvironmentVariable('Path', 'Machine') + ';' + [Environment]::GetEnvironmentVariable('Path', 'User')
}
& g++ -O3 -march=native -std=c++20 -static -pthread -o .\m18_complete_hard_search.exe .\m18_complete_hard_search.cpp
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
