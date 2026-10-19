@echo off
cd /d "%~dp0"
if not exist gallery.html (
 echo No gallery found. Run run.bat first.
 pause
 exit /b 1
)
start "" "%~dp0gallery.html"
