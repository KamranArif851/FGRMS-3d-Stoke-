@echo off
setlocal
cd /d "%~dp0"
set "CC=E:\install_c\bin\gcc.exe"
if not exist "%CC%" set "CC=gcc"
"%CC%" -O3 -std=c99 -Wall -Wextra -Wpedantic -Werror -Wl,--large-address-aware iteration.c -o iteration.exe -lm
if errorlevel 1 exit /b 1
"%CC%" -O3 -std=c99 -Wall -Wextra -Wpedantic -Werror -Wl,--large-address-aware comparison_plots.c -o comparison_plots.exe -lm
exit /b %errorlevel%
