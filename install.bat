@echo off

echo Library installer script
echo ------------------------
echo[

set LIBRARY_NAME=PelvicLandmarkIdentification
set VERSION=2.0.0

if "%MATLAB_LIBS%" == "" (
    echo Error: The environment variable 'MATLAB_LIBS' is not set.
    pause
    exit -1
) else (
    echo Installing "%LIBRARY_NAME%-%VERSION%" in "%MATLAB_LIBS%".
    echo Close this window to abort the installation.
    timeout /t 15

    if not exist "%MATLAB_LIBS%" (
        echo Creating folder "%MATLAB_LIBS%"...
        mkdir "%MATLAB_LIBS%"
    )

    if not exist "%MATLAB_LIBS%\%LIBRARY_NAME%-%VERSION%" (
        echo Creating folder "%LIBRARY_NAME%-%VERSION%" in "%MATLAB_LIBS%"...
        mkdir "%MATLAB_LIBS%\%LIBRARY_NAME%-%VERSION%"
    )

    echo Copying the data...
    if exist exclude.txt del exclude.txt
    echo %~dp0.git >> exclude.txt
    echo %~dp0tests\cache >> exclude.txt
    echo exclude.txt >> exclude.txt
    xcopy "%~dp0." "%MATLAB_LIBS%\%LIBRARY_NAME%-%VERSION%\" /i /f /y /e /EXCLUDE:exclude.txt
    del exclude.txt

    pause
)
