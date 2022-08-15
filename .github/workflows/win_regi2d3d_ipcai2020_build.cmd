@echo off

setlocal

set "CUR_DRIVE=%~d0" || EXIT /b

set "PATH=%CUR_DRIVE%\tmp\cmake-3.22.2-windows-x86_64\bin;%CUR_DRIVE%\tmp\ninja-bin;%PATH%" || EXIT /b

set "INSTALL_ROOT=%CUR_DRIVE%\usr\local" || EXIT /b

set "XREG_BUILD_PATH=%CUR_DRIVE%\tmp\xreg_build" || EXIT /b

cd regi2d3d-ipcai2020 || EXIT /b

mkdir build || EXIT /b

cd build || EXIT /b

cmake -G Ninja ^
    -Dxreg_DIR:PATH="%XREG_BUILD_PATH%" ^
    -DCMAKE_PREFIX_PATH:PATH="%INSTALL_ROOT%" ^
    -DCMAKE_INSTALL_PREFIX:PATH="%INSTALL_ROOT%" ^
    -DCMAKE_CXX_STANDARD:STRING="11" ^
    -DCMAKE_BUILD_TYPE:STRING="Release" ^
    -DBUILD_SHARED_LIBS:BOOL="OFF" ^
    ..  || EXIT /b

cmake --build . --config Release || EXIT /b

cmake --install . || EXIT /b
