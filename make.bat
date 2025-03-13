@echo off

REM Convenience wrapper for CMake commands

REM Script command (1st parameter)
set COMMAND=%1

echo Build cmake args: %BUILD_CMAKE_ARGS%

if "%OUTPUT_DIR%" == "" (
    set OUTPUT_DIR=out
)

set BUILD_DIR=%OUTPUT_DIR%\build
set INSTALL_DIR=%OUTPUT_DIR%\install

set SOURCE_DIR=%cd%

if "%LUX_PYTHON%" == "" (
    set LUX_PYTHON=python
)

if "%COMMAND%" == "" (
    call :Luxcore
    call :PyLuxcore
    call :LuxcoreUI
    call :LuxcoreConsole
    call :InvokeCMakeInstall
) else if "%COMMAND%" == "luxcore" (
    call :Luxcore
    call :InvokeCMakeInstall luxcore
) else if "%COMMAND%" == "pyluxcore" (
    call :PyLuxcore
    call :InvokeCMakeInstall pyluxcore
) else if "%COMMAND%" == "luxcoreui" (
    call :LuxcoreUI
    call :InvokeCMakeInstall luxcore
    call :InvokeCMakeInstall luxcoreui
) else if "%COMMAND%" == "luxcoreconsole" (
    call :LuxcoreConsole
    call :InvokeCMakeInstall luxcore
    call :InvokeCMakeInstall luxcoreconsole
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

:InvokeCMakeConfig
setlocal
for /f "delims=" %%A in ('python cmake\get_preset.py') do set "PRESET=%%A"
set PRESET=conan-release
echo CMake preset: %PRESET%
cmake %BUILD_CMAKE_ARGS% --preset %PRESET% -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR% -S %SOURCE_DIR% -G Ninja
endlocal
goto :EOF

:InvokeCMakeBuild
setlocal
for /f "delims=" %%A in ('python cmake\get_preset.py') do set "PRESET=%%A"
set PRESET=conan-release
echo CMake preset: %PRESET%
set TARGET=%1
cmake --build --preset %PRESET% --target %TARGET% %BUILD_CMAKE_ARGS%
endlocal
goto :EOF

:InvokeCMakeInstall
IF "%~1" == "" (
    cmake --install %BUILD_DIR%/Release --prefix %BUILD_DIR%
) else (
    cmake --install %BUILD_DIR%/Release --prefix %BUILD_DIR% --component %1
)
goto :EOF

:Clean
call :InvokeCMakeBuild clean
goto :EOF

:Config
call :InvokeCMakeConfig %CONAN_PRESET%
goto :EOF

:Luxcore
call :Config
call :InvokeCMakeBuild luxcore
goto :EOF

:PyLuxcore
call :Config
call :InvokeCMakeBuild pyluxcore
goto :EOF

:LuxcoreUI
call :Config
call :InvokeCMakeBuild luxcoreui
goto :EOF

:LuxcoreConsole
call :Config
call :InvokeCMakeBuild luxcoreconsole
goto :EOF

:Clear
rmdir /S /Q %BUILD_DIR%
goto :EOF

:Deps
call python -u cmake\make_deps.py
goto :EOF

:ListPresets
cmake --list-presets
goto :EOF

:EOF
