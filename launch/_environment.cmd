@echo off
rem Resolve a usable Python interpreter without invoking conda or mamba activation.
set "PY="
set "ENV_ROOT=%USERPROFILE%\anaconda3\envs\gwas_env"
if exist "%ENV_ROOT%\python.exe" (
  set "PY=%ENV_ROOT%\python.exe"
  set "PATH=%ENV_ROOT%;%ENV_ROOT%\Library\mingw-w64\bin;%ENV_ROOT%\Library\usr\bin;%ENV_ROOT%\Library\bin;%ENV_ROOT%\Scripts;%ENV_ROOT%\bin;%ENV_ROOT%\DLLs;%PATH%"
)
if not defined PY (
  for /f "delims=" %%P in ('where python 2^>nul') do if not defined PY set "PY=%%P"
)
if not defined PY (
  echo ERROR: Python was not found. Install the packages in requirements.txt or update launch\_environment.cmd.
  exit /b 1
)
set "PYTHONUTF8=1"
set "PYTHONIOENCODING=utf-8"
set "PYTHONUNBUFFERED=1"
set "PYTHONFAULTHANDLER=1"
set "OMP_NUM_THREADS=1"
set "OPENBLAS_NUM_THREADS=1"
set "MKL_NUM_THREADS=1"
set "NUMEXPR_NUM_THREADS=1"
set "VECLIB_MAXIMUM_THREADS=1"
exit /b 0
