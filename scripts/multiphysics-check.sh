#!/usr/bin/env sh
# Build, test and lint everything the Multiphysics workspace changes touch, and keep one log
# (multiphysics-check.log). Every step runs even if an earlier one fails. Pass --offline after
# unpacking a vendor archive from offline-vendor.sh.
set -u
repo=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
log="$repo/multiphysics-check.log"
net=""
[ "${1:-}" = "--offline" ] && net="--offline"
crates="-p solidcraft-multiphysics -p solidcraft-engine -p solidcraft-ui-egui"
echo "multiphysics-check $(date -u +%FT%TZ)" > "$log"
echo "rustc: $(rustc --version)" >> "$log"
failed=""
step() {
    name=$1; dir=$2; shift 2
    printf '\n==== %s ====\n> cargo %s\n' "$name" "$*" >> "$log"
    (cd "$dir" && cargo "$@") >> "$log" 2>&1
    code=$?
    echo "exit code: $code" >> "$log"
    if [ "$code" -eq 0 ]; then printf '%-40s ok\n' "$name"; else printf '%-40s FAILED (%s)\n' "$name" "$code"; failed="$failed; $name"; fi
}
# shellcheck disable=SC2086
step "format" "$repo" fmt --check $crates
step "check (all targets)" "$repo" check --all-targets $net $crates
step "clippy" "$repo" clippy --all-targets $net $crates -- -D warnings
step "test multiphysics crate" "$repo" test $net -p solidcraft-multiphysics
step "test engine (multiphysics)" "$repo" test $net -p solidcraft-engine multiphysics
step "test ui (multiphysics)" "$repo" test $net -p solidcraft-ui-egui multiphysics
step "test companion" "$repo/companion/cfd-solvers-rs" test --locked $net
if [ -z "$failed" ]; then echo "All steps passed. Log: $log"; else echo "Failed:${failed#;}. Send $log to have them fixed."; fi
