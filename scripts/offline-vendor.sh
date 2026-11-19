#!/usr/bin/env sh
# Pack every crates.io dependency of SolidCraft and its companion into one archive, so the
# workspace builds with no network (cargo --offline) on another machine. Unpack at the
# repository root and use vendor-crates/cargo-config.toml as .cargo/config.toml.
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
out=${1:-solidcraft-vendor.tar.gz}
cd "$repo"
rm -rf vendor-crates
cargo vendor --locked --versioned-dirs --sync companion/cfd-solvers-rs/Cargo.toml vendor-crates > vendor-crates.config.tmp
mv vendor-crates.config.tmp vendor-crates/cargo-config.toml
tar -czf "$out" vendor-crates
echo "Packed $(ls vendor-crates | wc -l) crates into $out ($(du -h "$out" | cut -f1))."
echo "Too large for one upload? split -b 95m $out $out.part-"
