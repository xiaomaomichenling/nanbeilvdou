$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
& .\m18_complete_hard_search.exe @args
exit $LASTEXITCODE
