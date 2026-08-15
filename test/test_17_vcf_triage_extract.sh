# test_17_vcf_triage_extract.sh — vcf-triage extract: CLI, sections, thresholds, degradation.
XT="$REPO/.claude/skills/vcf-triage/extract-triage.sh"
XDIR="$TMP/triage"; mkdir -p "$XDIR"

# --- CLI ---
out="$(bash "$XT" 2>&1)"; rc=$?
assert_eq "1" "$rc" "triage: no args exits 1"
assert_contains "$out" "Usage" "triage: prints usage on no args"

# --- Missing everything: INPUTS notes each absent file, exit 0 ---
out="$(bash "$XT" "$XDIR" nothere 2>&1)"; rc=$?
assert_eq "0" "$rc" "triage: missing inputs still exits 0"
assert_contains "$out" "== INPUTS ==" "triage: INPUTS section present"
assert_contains "$out" "clinvar: MISSING" "triage: notes missing clinvar table"
assert_contains "$out" "gwas: MISSING" "triage: notes missing gwas table"
assert_contains "$out" "xref-traits.sh" "triage: missing note names the producing script"
assert_contains "$out" "== PHARMGKB ==" "triage: PHARMGKB section present"
assert_contains "$out" "absent" "triage: pharmgkb absent note"

# --- Single-artifact path form derives dir/base ---
printf 'x\n' > "$XDIR/mybase.traits.clinvar.tsv"
out="$(bash "$XT" "$XDIR/mybase.traits.clinvar.tsv" 2>&1)"
assert_contains "$out" "clinvar: found" "triage: artifact-path form finds sibling table"
rm "$XDIR/mybase.traits.clinvar.tsv"
