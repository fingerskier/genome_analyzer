# test_16_xref_entry.sh — xref-traits.sh CLI + graceful degradation (no DBs).
sample="$(make_fixture "$REPO/test/fixtures/xref_sample.vcf.txt")"

# Missing input -> nonzero exit, usage shown.
out="$("$REPO/xref-traits.sh" 2>&1)"; rc=$?
assert_eq "1" "$rc" "entry: no -i exits 1"
assert_contains "$out" "Usage" "entry: prints usage when -i missing"

# Both DBs absent -> still writes a SUMMARY with skip notes, exit 0.
od="$TMP/entry_out"
"$REPO/xref-traits.sh" -i "$sample" -o "$od" -g "$TMP/nope.tsv" -c "$TMP/nope.vcf.gz" >/dev/null 2>&1
rc=$?
assert_eq "0" "$rc" "entry: missing DBs still exits 0"
sum="$(cat "$od/xref_sample.traits.SUMMARY.md" 2>/dev/null)"
assert_contains "$sum" "fetch-gwas.sh"    "entry: GWAS skip note guides the user to fetch-gwas.sh"
assert_contains "$sum" "fetch-clinvar.sh" "entry: ClinVar skip note guides the user to fetch-clinvar.sh"
