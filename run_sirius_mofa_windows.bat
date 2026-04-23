@echo off
setlocal EnableExtensions
cd /d "%~dp0"

set "OLLAMA_HOST=http://127.0.0.1:11434"

rem --- Cles API optionnelles pour Claude / DeepSeek ---
rem --- Decommentez et remplissez selon le(s) provider(s) utilise(s) ---
rem set "ANTHROPIC_API_KEY=sk-ant-xxxxxxxxxxxxxxxxxxxx"
rem set "DEEPSEEK_API_KEY=sk-xxxxxxxxxxxxxxxxxxxx"

echo AI Sirius MOFA launcher

set "VENV_PY=%~dp0.venv\Scripts\python.exe"
set "VENV_PIP=%~dp0.venv\Scripts\pip.exe"
set "REQ=%~dp0requirements.txt"

if exist "%VENV_PY%" (
    echo Using venv Python: "%VENV_PY%"
    "%VENV_PY%" --version

    rem --- Pour reinstaller les dependances, decommentez la ligne ci-dessous ---
    rem "%VENV_PIP%" install -r "%REQ%" --quiet --timeout 120

    echo Starting AI Sirius MOFA...
    "%VENV_PY%" "%~dp0ai_sirius_mofa.py"
    goto :done
)

rem --- Fallback: system Python ---
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: python not found in PATH and no .venv found.
    pause
    exit /b 1
)

echo Starting AI Sirius MOFA...
python "%~dp0ai_sirius_mofa.py"

:done
pause
