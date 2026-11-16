#!/usr/bin/env sh
set -eu
cd "$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
if ! command -v cargo >/dev/null 2>&1; then
    printf '%s\n' 'Install Rust 1.85.1 or newer, including your platform linker, then run this script again.' >&2
    exit 2
fi
cargo test --locked --all-features
cargo build --locked --release --all-features
printf '%s\n' 'Built target/release/cfd-solve. Run: target/release/cfd-solve --help'
