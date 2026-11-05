#!/usr/bin/env bash
# Compile each example ODE of test/Conversion/ODEToACR/examples to a reaction
# network, simulate it with the reservoir emulator and print how far the
# network's trajectories are from the ODE solution (MSE per state).
#
#   utils/check-odes.sh            # all examples
#   utils/check-odes.sh FILE...    # selected ode.system files
#
# Needs chemcomp-opt / chemcomp-translate in $CHEMCOMP_BUILD_DIR/bin and the
# emulator package (pixi run check-odes in the default environment).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="${CHEMCOMP_BUILD_DIR:-$ROOT/build/pixi-default}/bin"
CHECK="$ROOT/tools/chemcomp-reservoir-emulator/scripts/check_realization.py"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

files=("$@")
[[ ${#files[@]} -gt 0 ]] || files=("$ROOT"/test/Conversion/ODEToACR/examples/*.mlir)

status=0
for file in "${files[@]}"; do
  # The simulated horizon is the --t-end of the file's own RUN line.
  t_end="$(sed -n 's/.*--t-end \([0-9.]*\).*/\1/p' "$file" | head -1)"
  "$BIN/chemcomp-opt" "$file" --ode-to-crn \
    | "$BIN/chemcomp-translate" -mlir-to-json > "$WORK/crn.json"
  "$BIN/chemcomp-translate" -ode-to-json "$file" > "$WORK/ode.json"
  python -s "$CHECK" --ode "$WORK/ode.json" --crn "$WORK/crn.json" \
    --t-end "${t_end:-10}" | grep -v '^Created output directory' || status=1
done
exit $status
