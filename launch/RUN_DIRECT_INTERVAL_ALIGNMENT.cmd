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
set "OUT=%PROJECT%\revision_validation_outputs\ALIGNMENT_FINAL_%STAMP%"
set "FINAL=%OUT%\final"
set "LOG=%OUT%\alignment_console.log"
mkdir "%OUT%" 2>nul

echo ============================================================
echo DIRECT CHROMOSOME 4 INTERVAL ALIGNMENT
echo ============================================================
echo Project : %PROJECT%
echo Python  : %PY%
echo Output  : %OUT%
echo.

"%PY%" "%REPO%\scripts\run_and_tee.py" --log "%LOG%" "%PY%" -X faulthandler "%REPO%\scripts\02_direct_interval_alignment.py" ^
  --project-dir "%PROJECT%" ^
  --output-dir "%OUT%" ^
  --intervals-csv "%REPO%\config\analysis_intervals.csv" ^
  --threads 4
if errorlevel 1 (
  echo ERROR: Direct alignment failed. Review %LOG%
  exit /b 1
)

"%PY%" "%REPO%\scripts\03_finalize_alignment_annotations.py" ^
  --project-dir "%PROJECT%" ^
  --alignment-dir "%OUT%" ^
  --output-dir "%FINAL%"
if errorlevel 1 (
  echo ERROR: Alignment annotation finalization failed.
  exit /b 1
)

echo.
echo Completed: %FINAL%
explorer "%FINAL%"
endlocal
