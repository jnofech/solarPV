SETLOCAL EnableExtensions
set EXE=%1
FOR /F %%x IN ('tasklist /NH /FI "IMAGENAME eq %EXE%"') DO IF %%x == %EXE% goto FOUND
echo Not running!
goto FIN
:FOUND
echo Running!!
:FIN