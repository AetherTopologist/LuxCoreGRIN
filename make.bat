@echo off

REM Convenience wrapper for CMake commands

REM Script command (1st parameter)
set COMMAND=%1

for /f "delims=" %%A in ('python cmake\get_preset.py') do set "CONAN_PRESET=%%A"
echo CMake preset: %CONAN_PRESET%

echo Build cmake args: %BUILD_CMAKE_ARGS%

if "%BUILD_DIR%" == "" (
    set BUILD_DIR=.\build
)

if "%INSTALL_DIR%" == "" (
    set INSTALL_DIR=.\build
)

set SOURCE_DIR=%cd%

if "%LUX_PYTHON%" == "" (
    set LUX_PYTHON=python.exe
)

if "%COMMAND%" == "" (
    call :Luxcore
    call :PyLuxcore
    call :LuxcoreUI
    call :LuxcoreConsole
) else if "%COMMAND%" == "luxcore" (
    call :Luxcore
) else if "%COMMAND%" == "pyluxcore" (
    call :PyLuxcore
) else if "%COMMAND%" == "luxcoreui" (
    call :LuxcoreUI
) else if "%COMMAND%" == "luxcoreconsole" (
    call :LuxcoreConsole
) else if "%COMMAND%" == "config" (
    call :Config
) else if "%COMMAND%" == "install" (
    call :InvokeCMakeInstall
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

:InvokeCMake
setlocal
set PRESET=%1
cmake --build --preset %PRESET%
endlocal
goto :EOF

:InvokeCMakeConfig
setlocal
set PRESET=%1
set TARGET=%2
cmake %BUILD_CMAKE_ARGS% -D%TARGET%=ON --preset %PRESET% -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR% -S %SOURCE_DIR%
endlocal
goto :EOF

:InvokeCMakeInstall
cmake --install %BUILD_DIR%/cmake --prefix %BUILD_DIR%
goto :EOF

:Clean
call :InvokeCMake %CONAN_PRESET% clean
goto :EOF

:Config
set TARGET=%1
call :InvokeCMakeConfig %CONAN_PRESET% %TARGET%
goto :EOF

:Luxcore
call :Config LUXCORE_LIBONLY
call :InvokeCMake %CONAN_PRESET%
goto :EOF

:PyLuxcore
call :Config LUXCORE_PYLUXCORE
call :InvokeCMake %CONAN_PRESET%
goto :EOF

:LuxcoreUI
call :Config LUXCORE_UI
call :InvokeCMake %CONAN_PRESET%
goto :EOF

:LuxcoreConsole
call :Config LUXCORE_CON
call :InvokeCMake %CONAN_PRESET%
goto :EOF

:Clear
rmdir /S /Q %BUILD_DIR%
goto :EOF

:Deps
%LUX_PYTHON% -u cmake\make_deps.py
goto :EOF

:ListPresets
cmake --list-presets
goto :EOF

:EOF
