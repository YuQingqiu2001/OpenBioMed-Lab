@echo off
REM Visible Chrome launcher for bioRxiv/medRxiv browser extraction.
REM This deliberately uses a normal visible Chrome window, not headless Chrome.
REM You may need to manually complete Cloudflare/Turnstile once in the opened window.

setlocal
set DEBUG_PORT=9223
if defined LITERATURE_KG_ROOT (
  set PROFILE_DIR=%LITERATURE_KG_ROOT%\chrome_cdp_profile
) else (
  set PROFILE_DIR=%~dp0..\literature-workspace\chrome_cdp_profile
)
set START_URL=https://www.biorxiv.org/

set CHROME=%ProgramFiles%\Google\Chrome\Application\chrome.exe
if not exist "%CHROME%" set CHROME=%ProgramFiles(x86)%\Google\Chrome\Application\chrome.exe
if not exist "%CHROME%" set CHROME=%LocalAppData%\Google\Chrome\Application\chrome.exe

if not exist "%CHROME%" (
  echo [ERROR] Cannot find chrome.exe. Edit this file and set CHROME manually.
  pause
  exit /b 1
)

if not exist "%PROFILE_DIR%" mkdir "%PROFILE_DIR%"

echo [bioRxiv CDP] Starting visible Chrome on http://127.0.0.1:%DEBUG_PORT%
echo [bioRxiv CDP] Profile: %PROFILE_DIR%
echo [bioRxiv CDP] If Cloudflare appears, complete it manually in this Chrome window.

start "bioRxiv Chrome CDP" "%CHROME%" ^
  --remote-debugging-port=%DEBUG_PORT% ^
  --remote-debugging-address=127.0.0.1 ^
  --user-data-dir="%PROFILE_DIR%" ^
  --no-first-run ^
  --no-default-browser-check ^
  --disable-sync ^
  --disable-features=ChromeSignin,SigninIntercept,AccountConsistency,EnableEphemeralGuestProfilesOnDesktop ^
  "%START_URL%"

echo.
echo Verify from Git-Bash/Hermes terminal:
echo   curl http://127.0.0.1:%DEBUG_PORT%/json/version
echo.
echo Then run, for example:
echo   node "%~dp0extract_biorxiv_cdp.mjs" --doi 10.64898/2026.05.31.727600 --port %DEBUG_PORT%
echo.
endlocal
