# test/lib.sh — zero-dependency assertions + fixture builder. Sourced by run-tests.sh.
PASS=0; FAIL=0

assert_eq() {  # expected actual message
  if [[ "$1" == "$2" ]]; then PASS=$((PASS+1)); printf '  ok   %s\n' "$3"
  else FAIL=$((FAIL+1)); printf '  FAIL %s\n    expected: [%s]\n    actual:   [%s]\n' "$3" "$1" "$2"; fi
}

assert_contains() {  # haystack needle message
  if printf '%s' "$1" | grep -qF -- "$2"; then PASS=$((PASS+1)); printf '  ok   %s\n' "$3"
  else FAIL=$((FAIL+1)); printf '  FAIL %s\n    missing: [%s]\n' "$3" "$2"; fi
}

# make_fixture <fixtures/name.vcf.txt> -> prints path to a bgzipped+indexed VCF in $TMP
make_fixture() {
  local src="$1" base out
  base="$(basename "${src%.vcf.txt}")"
  out="$TMP/${base}.vcf.gz"
  bgzip -c "$src" > "$out"
  tabix -f -p vcf "$out"
  printf '%s' "$out"
}

finish() { printf '\n%d passed, %d failed\n' "$PASS" "$FAIL"; [[ "$FAIL" -eq 0 ]]; }
