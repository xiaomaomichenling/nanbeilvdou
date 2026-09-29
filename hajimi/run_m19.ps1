$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath $PSScriptRoot
& .\m18_complete_hard_search.exe --M 19 @args
exit $LASTEXITCODE
