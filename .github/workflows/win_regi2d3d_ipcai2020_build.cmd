@echo off

setlocal

set CUR_DRIVE="%~d0" || EXIT /b

set "PATH=%CUR_DRIVE%\tmp\cmake-3.22.2-windows-x86_64\bin;%PATH%" || EXIT /b

cd regi2d3d-ipcai2020 || EXIT /b

mkdir build || EXIT /b

cd build || EXIT /b

cmake -G Ninja -DCMAKE_MAKE_PROGRAM:PATH="%CUR_DRIVE%/tmp/ninja-bin/ninja.exe" ^
    -Dxreg_DIR:PATH="%CUR_DRIVE%/tmp/xreg_build:" ^
    -DCMAKE_PREFIX_PATH:PATH="%CUR_DRIVE%/usr/local" ^
    -DCMAKE_CXX_STANDARD:STRING="11" ^
    -DCMAKE_BUILD_TYPE:STRING="Release" ^
    -DBUILD_SHARED_LIBS:BOOL="OFF" ^
    ..  || EXIT /b

cmake --build . --config Release || EXIT /b
