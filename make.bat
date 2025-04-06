@echo off

REM SPDX-FileCopyrightText: 2024-2025 Authors (see AUTHORS.txt)
REM
REM SPDX-License-Identifier: Apache-2.0


REM Convenience wrapper for CMake commands

REM Script command (1st parameter)
set COMMAND=%1


if "%LUX_PYTHON%" == "" (
    set LUX_PYTHON=python
)

set LUX_CMAKE="%LUX_PYTHON%" -u cmake\cmake.py

if "%COMMAND%" == "" (
    call :Config
    call :BuildAndInstall luxcore
    call :BuildAndInstall luxcoreui
    call :BuildAndInstall luxcoreconsole
    call :BuildAndInstall pyluxcore
) else if "%COMMAND%" == "luxcore" (
    call :Config
    call :BuildAndInstall luxcore
) else if "%COMMAND%" == "pyluxcore" (
    call :Config
    call :BuildAndInstall luxcore
    call :BuildAndInstall pyluxcore
) else if "%COMMAND%" == "luxcoreui" (
    call :Config
    call :BuildAndInstall luxcore
    call :BuildAndInstall luxcoreui
) else if "%COMMAND%" == "luxcoreconsole" (
    call :Config
    call :BuildAndInstall luxcore
    call :BuildAndInstall luxcoreconsole
) else if "%COMMAND%" == "config" (
    call :Config
) else if "%COMMAND%" == "package" (
    call :BuildAndInstall package
) else if "%COMMAND%" == "install" (
    call :Install
) else if "%COMMAND%" == "clean" (
    call :Clean
) else if "%COMMAND%" == "clear" (
    call :Clear
) else if "%COMMAND%" == "deps" (
    call :Deps
) else if "%COMMAND%" == "list-presets" (
    call :ListPresets
) else (
    echo Command "%COMMAND%" unknown
)
exit /B

:Deps
call %LUX_CMAKE% deps
goto :EOF

:ListPresets
call %LUX_CMAKE% list-presets
goto :EOF

:Config
call %LUX_CMAKE% config
goto :EOF

:BuildAndInstall
call %LUX_CMAKE% build-and-install %1
goto :EOF

:Install
IF "%~1" == "" (
    %LUX_CMAKE% all
) else (
    %LUX_CMAKE% %1
)
goto :EOF

:Clean
call %LUX_CMAKE% clean
goto :EOF

:Clear
REM rmdir /S /Q
call %LUX_CMAKE% clear
goto :EOF

:EOF
