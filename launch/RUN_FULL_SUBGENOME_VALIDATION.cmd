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
set "OUT=%PROJECT%\revision_validation_outputs\FULL_SENSITIVITY_FILTERED_%STAMP%"
set "LOG=%OUT%\full_sensitivity_console.log"
mkdir "%OUT%" 2>nul

echo ============================================================
echo FULL STRICT SUBGENOME-FILTERED VALIDATION
echo ============================================================
echo Project : %PROJECT%
echo Python  : %PY%
echo Output  : %OUT%
echo.

"%PY%" "%REPO%\scripts\native_stack_preflight.py"
if errorlevel 1 (
  echo ERROR: Native numerical stack preflight failed.
  exit /b 1
)

"%PY%" "%REPO%\scripts\run_and_tee.py" --log "%LOG%" "%PY%" -X faulthandler "%REPO%\scripts\01_subgenome_filtered_sensitivity.py" ^
  --project-dir "%PROJECT%" ^
  --output-dir "%OUT%" ^
  --panel-sizes "6000,12000,24000,48000" ^
  --seeds "20250416,20250417,20250418,20250419,20250420" ^
  --submitted-panel-size 12000 ^
  --submitted-seed 20250416 ^
  --min-site-call-rate 0.80 ^
  --maf-threshold 0.05 ^
  --window-size-bp 1000000 ^
  --top-markers 250 ^
  --n-permutations 500 ^
  --min-window-markers 5 ^
  --n-components 5 ^
  --progress-every 250000 ^
  --batch-size 4096 ^
  --contig-filter-mode matching_pseudomolecules_only
if errorlevel 1 (
  echo ERROR: Full sensitivity analysis failed. Review %LOG%
  exit /b 1
)

echo.
echo Completed: %OUT%
explorer "%OUT%"
endlocal
