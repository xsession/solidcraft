#!/usr/bin/env sh
# Build the separate companion before launching SolidCraft with its explicit path.
set -eu
mp_repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
mp_companion="$mp_repo/companion/cfd-solvers-rs"
(
    cd "$mp_companion"
    cargo build --release --locked --target-dir "$mp_companion/target"
)
SOLIDCRAFT_CFD_SOLVE="$mp_companion/target/release/cfd-solve"
export SOLIDCRAFT_CFD_SOLVE
cd "$mp_repo"
exec cargo run --release -p solidcraft -- "$@"
