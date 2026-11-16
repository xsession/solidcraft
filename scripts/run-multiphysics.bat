@echo off
setlocal
for %%I in ("%~dp0..") do set "MP_REPO=%%~fI"
set "MP_COMPANION=%MP_REPO%\companion\cfd-solvers-rs"
pushd "%MP_COMPANION%" || exit /b 1
call cargo build --release --locked --target-dir "%MP_COMPANION%\target"
if errorlevel 1 (
    popd
    exit /b 1
)
popd
set "SOLIDCRAFT_CFD_SOLVE=%MP_COMPANION%\target\release\cfd-solve.exe"
pushd "%MP_REPO%" || exit /b 1
call cargo run --release -p solidcraft -- %*
set "MP_EXIT=%ERRORLEVEL%"
popd
exit /b %MP_EXIT%
