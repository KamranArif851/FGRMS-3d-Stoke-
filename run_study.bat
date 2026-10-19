@echo off
cd /d "%~dp0"
call build.bat
if errorlevel 1 goto :error
if not exist results mkdir results
for %%C in (q4_n2 q4_n3 q4_443) do (
 iteration.exe cases/%%C.cfg results/%%C
 if errorlevel 1 goto :error
)
comparison_plots.exe
if errorlevel 1 goto :error
copy /y gallery.html study_gallery.html >nul
start "" gallery.html
pause
exit /b 0
:error
 echo Failed. Read the message above.
 pause
 exit /b 1
