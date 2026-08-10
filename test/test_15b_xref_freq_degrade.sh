# test_15b_xref_freq_degrade.sh — frequency degradation: legacy 7-col GWAS
# cache (pre-raf) and a ClinVar VCF without AF_* tags must not break the run.
( source "$REPO/common.sh"; source "$REPO/lib/trait-xref.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/xref_sample.vcf.txt")"
  CLINVAR_VCF="$(make_fixture "$REPO/test/fixtures/clinvar_noaf.vcf.txt")"
  GWAS_TBL="$REPO/test/fixtures/gwas_norm_legacy.tsv"
  OUTDIR="$TMP/xref_degrade_out"; mkdir -p "$OUTDIR"; BASE="xrefd"
  run_trait_xref )
rc=$?
assert_eq "0" "$rc" "freq degrade: run exits 0"

dsum="$(cat "$TMP/xref_degrade_out/xrefd.traits.SUMMARY.md")"
dgtsv="$(cat "$TMP/xref_degrade_out/xrefd.traits.gwas.tsv")"
dctsv="$(cat "$TMP/xref_degrade_out/xrefd.traits.clinvar.tsv")"

assert_eq "0" "$(printf '%s\n' "$dgtsv" | awk -F'\t' '$11!="NA"' | grep -c .)" "freq degrade: legacy table -> every raf is NA"
assert_contains "$dsum" "fetch-gwas.sh -f" "freq degrade: summary tells the user how to refresh the catalog"
assert_contains "$dctsv" "NoAF condition" "freq degrade: AF-less ClinVar hit still reported"
assert_contains "$dctsv" "NA	NA" "freq degrade: AF-less ClinVar VCF -> pop_af/source NA"
