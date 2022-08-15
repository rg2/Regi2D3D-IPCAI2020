echo on

REM Create a .zip of the compiled binaries

REM Use the drive where this script currently resides
SET "CUR_DRIVE=%~d0"
ECHO CUR_DRIVE = %CUR_DRIVE%

REM Everything will be installed here using the standard
REM bin, include, lib layout
SET "INSTALL_ROOT=%CUR_DRIVE%\usr\local"
ECHO INSTALL_ROOT = %INSTALL_ROOT%

SET "XREG_NAME=regi2d3d-ipcai2020-win64"
MKDIR "%XREG_NAME%" || EXIT /b
MKDIR "%XREG_NAME%\bin" || EXIT /b

XCOPY /E "%INSTALL_ROOT%\bin\OpenCL.dll" "%XREG_NAME%\bin\"
XCOPY /E "%INSTALL_ROOT%\bin\tbb.dll" "%XREG_NAME%\bin\"
XCOPY /E "%INSTALL_ROOT%\bin\setup-xreg-vars.bat" "%XREG_NAME%\bin\"
XCOPY /E "%INSTALL_ROOT%\bin\ffmpeg.exe" "%XREG_NAME%\bin\"
XCOPY /E "%INSTALL_ROOT%\bin\xreg-ipcai-*" "%XREG_NAME%\bin\"

XCOPY /E  regi2d3d-ipcai2020\README.md "%XREG_NAME%\bin\"
XCOPY /E  regi2d3d-ipcai2020\LICENSE "%XREG_NAME%\bin\"

REM GitHub runners have 7zip installed, so we will use that to create the .zip
7z a -tzip "%XREG_NAME%.zip" "%XREG_NAME%"
