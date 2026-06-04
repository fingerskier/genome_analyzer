#!/usr/bin/env bash
set -uo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$DIR/.." && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
export REPO TMP
# shellcheck disable=SC1091
source "$DIR/lib.sh"
for t in "$DIR"/test_*.sh; do
  [[ -e "$t" ]] || continue
  printf '\n=== %s ===\n' "$(basename "$t")"
  # Isolate each test file: a test may source a script that runs `set -e`
  # (e.g. analyze.sh), and errexit would otherwise leak into later files and
  # abort the runner on the first nonzero top-level command. Reset to the
  # runner's baseline (nounset + pipefail, no errexit) before each file.
  set +e -uo pipefail
  # shellcheck disable=SC1090
  source "$t"
done
finish
