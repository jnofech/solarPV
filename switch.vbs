Dim WshShell
Set WshShell = WScript.CreateObject("WScript.Shell")
Dim ARGS
set ARGS=WScript.Arguments
WshShell.AppActivate(ARGS.Item(0))
' WshShell.SendKeys("~")
' WScript.Echo(ARGS.Item(0) + " activated")
WScript.Quit(0)