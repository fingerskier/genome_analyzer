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
  # shellcheck disable=SC1090
  source "$t"
done
finish
