@echo off
setlocal
pushd "%~dp0"
where cargo >nul 2>nul
if errorlevel 1 (
  echo Install Rust 1.85.1 and the Windows C++ Build Tools linker, then run this script again.
  popd
  exit /b 2
)
cargo test --locked --all-features
if errorlevel 1 goto failed
cargo build --locked --release --all-features
if errorlevel 1 goto failed
echo Built target\release\cfd-solve.exe
popd
exit /b 0
:failed
popd
exit /b 1
