@echo off
cd /d "%~dp0"
call build.bat
if errorlevel 1 goto :error
if not exist results mkdir results
iteration.exe input.cfg results/main
if errorlevel 1 goto :error
comparison_plots.exe results/main_input.cfg results/main
if errorlevel 1 goto :error
start "" gallery.html
pause
exit /b 0
:error
 echo Failed. Read the message above.
 pause
 exit /b 1
