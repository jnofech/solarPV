@Set "EX=%1"
@tasklist /fi "status eq not responding" /nh | find "%EX%" >NUL
@if %errorlevel% == 0 (
    echo Application not responding!!
)