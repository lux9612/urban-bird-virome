@echo off
setlocal
cd /d "%~dp0"
if exist .venv\Scripts\python.exe goto run
where py >nul 2>nul
if errorlevel 1 (python -m venv .venv) else (py -3.12 -m venv .venv)
if errorlevel 1 goto error
:run
.venv\Scripts\python.exe -c "import sys; sys.exit(0 if sys.version_info >= (3, 12) else 1)"
if errorlevel 1 goto error
.venv\Scripts\python.exe -m pip install -r requirements.txt
if errorlevel 1 goto error
.venv\Scripts\python.exe plot_fig7b.py
if errorlevel 1 goto error
echo Finished. See the output folder.
pause
exit /b 0
:error
echo Python 3.12 or newer is required. See the error above.
pause
exit /b 1
