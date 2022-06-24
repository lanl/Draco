echo off
rem -----------------------------------------------------------------------------------------------
rem File  : gitlab-nightly-regress.bat
rem Date  : Friday, Aug 14, 2020, 14:20 pm
rem Author: Kelly Thompson <kgt@lanl.gov>
rem Note  : Copyright (C) 2020-2022, Triad National Security, LLC., All rights are reserved.
rem -----------------------------------------------------------------------------------------------

echo This is windows...
rem set

echo ---------------------------------------
echo SITE_ID      = %SITE_ID%
echo DEPLOY       = %DEPLOY%
echo VCVARS       = %VCVARS%
echo MINGW64PATH  = %MINGW64PATH%
echo DRACO_SOURCE_DIR = %DRACO_SOURCE_DIR%
echo DRACO_BINARY_DIR = %DRACO_BINARY_DIR%
echo NUMBER_OF_PROCS  = %NUMBER_OF_PROCESSORS%
echo GENERATOR        = %GENERATOR%
echo CTEST_BUILD_NAME     = %CTEST_BUILD_NAME%
echo CMAKE_TOOLCHAIN_FILE = %CMAKE_TOOLCHAIN_FILE%
echo DEPLOY_DIR           = %DEPLOY_DIR%
echo DRACO_INSTALL_DIR    = %DRACO_INSTALL_DIR%
echo CMAKE_C_FLAGS        = %CMAKE_C_FLAGS%
echo CMAKE_Fortran_FLAGS  = %CMAKE_Fortran_FLAGS%
echo ---------------------------------------

rem Fix LANL proxy issues
rem - git must use a proxy (settings found in %USERPROFILE%/.gitconfig
rem - To allow cdash to post results, the environment variables HTTP[S] must not be set
git.exe config --global http.proxy http://proxyout.lanl.gov:8080
git.exe config --global https.proxy https://proxyout.lanl.gov:8080
set HTTP_PROXY=
set HTTPS_PROXY=

rem gfortran runtimes
set PATH=%PATH%;%MINGW64PATH%;=%
rem set PATH=%PATH%;C:\msys64\mingw64\bin;=%
rem numdiff
set PATH=%PATH%;C:\work\vendors64\bin;=%

echo PATH = %PATH%

rem -----------------------------------------------------------------------------------------------
rem Environment setup and checking
rem -----------------------------------------------------------------------------------------------

:setupmsvccommenv

dir %VCVARS% > %TEMP%\glr.log
if %ERRORLEVEL% == 0 (echo   vcvars64.bat   : found) else ( goto errnovcvars )
@call %VCVARS% %* > %TEMP%\glr.log
dir %CMAKE_TOOLCHAIN_FILE% > %TEMP%\glr.log
if %ERRORLEVEL% == 0 (echo   toolchain.cmake: found) else ( goto errnotoolchain )
dir %VCPKGLOC%  > %TEMP%\glr.log
if %ERRORLEVEL% == 0 (echo   vcpkg.cmake    : found) else ( goto errnovcpkg )
where /q cmake
if %ERRORLEVEL% == 0 (echo   cmake          : found) else ( goto errnocmake )
where /q numdiff
if %ERRORLEVEL% == 0 (echo   numdiff        : found) else ( goto errnonumdiff )

rem -----------------------------------------------------------------------------------------------
rem Main tests
rem -----------------------------------------------------------------------------------------------

:maintests

set modes=Configure,Build,Test,Submit
set script=%DRACO_SOURCE_DIR%\.gitlab\ci\draco-nightly.cmake

echo ctest -VV -C %CMAKE_BUILD_TYPE% -S %script%,%modes%
ctest -VV -C %CMAKE_BUILD_TYPE% -S %script%,%modes%

goto :done

rem -----------------------------------------------------------------------------------------------
rem Error handling
rem -----------------------------------------------------------------------------------------------

:errnovcvars
echo FATAL ERROR :: VCVARS = %VCVARS%
echo                was not found.  This batch file is required to setup the Visual Studio build
echo                environment.
goto:eof

:errnotoolchain
echo FATAL ERROR :: CMAKE_TOOLCHAIN_FILE = %CMAKE_TOOLCHAIN_FILE%
echo                was not found. This file is required.
goto:eof

:errnovcpkg
echo FATAL ERROR :: VCPKGLOC = %VCPKGLOC%
echo                was not found. This file is required.
goto:eof

:errnocmake
echo FATAL ERROR :: cmake was not found in PATH.
goto:eof

:errnonumdiff
echo FATAL ERROR :: numdiff was not found in PATH.
goto:eof

:done
echo
echo ======== end gitlab-nightly-regress.bat ==========
echo

rem -----------------------------------------------------------------------------------------------
rem End .gitlab/ci/gitlab-nightly-regress.bat
rem -----------------------------------------------------------------------------------------------
