@echo off
setlocal EnableExtensions
set "PROJECT=%~1"
if "%PROJECT%"=="" set "PROJECT=%~dp0..\data\raw_public_inputs"
set "REPO=%~dp0.."
call "%~dp0_environment.cmd"
if errorlevel 1 exit /b 1
"%PY%" -c "import numpy,pandas,openpyxl" >nul 2>&1
if errorlevel 1 (
  echo Installing required Python packages...
  "%PY%" -m pip install -r "%REPO%\requirements.txt"
  if errorlevel 1 exit /b 1
)
for /f %%T in ('powershell -NoProfile -Command "Get-Date -Format yyyyMMdd_HHmmss"') do set "STAMP=%%T"
set "OUT=%PROJECT%\analysis_exports_%STAMP%"

echo ============================================================
echo EXPORT SUPPLEMENTARY WORKBOOKS AND MACHINE-READABLE SOURCE DATA
echo ============================================================
echo Project : %PROJECT%
echo Python  : %PY%
echo Output  : %OUT%
echo.

"%PY%" "%REPO%\scripts\04_export_supplementary_outputs.py" ^
  --project-dir "%PROJECT%" ^
  --output-dir "%OUT%"
if errorlevel 1 (
  echo ERROR: Supplementary-output export failed.
  exit /b 1
)

echo.
echo Completed: %OUT%
explorer "%OUT%"
endlocal
