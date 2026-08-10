# test_15_xref.sh — end-to-end trait xref on fixtures (GWAS + ClinVar).
( source "$REPO/common.sh"; source "$REPO/lib/trait-xref.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/xref_sample.vcf.txt")"
  CLINVAR_VCF="$(make_fixture "$REPO/test/fixtures/clinvar.vcf.txt")"
  GWAS_TBL="$REPO/test/fixtures/gwas_norm.tsv"
  OUTDIR="$TMP/xref_out"; mkdir -p "$OUTDIR"; BASE="xref"
  run_trait_xref )

sum="$(cat "$TMP/xref_out/xref.traits.SUMMARY.md")"
gtsv="$(cat "$TMP/xref_out/xref.traits.gwas.tsv")"
ctsv="$(cat "$TMP/xref_out/xref.traits.clinvar.tsv")"

assert_contains "$gtsv" "Height"            "gwas: rs100 Height association reported"
assert_contains "$gtsv" "Coffee consumption" "gwas: rs200 Coffee association reported"
assert_eq "0" "$(printf '%s\n' "$gtsv" | grep -c 'Weak association')" "gwas: p=1e-6 below 5e-8 threshold excluded"
assert_contains "$sum" "**T/T**"            "gwas summary: homozygous risk genotype highlighted"
assert_contains "$sum" "strand-ambiguous"   "gwas summary: mentions the ambiguous-SNP bucket"
assert_contains "$gtsv" "ambiguous"         "gwas tsv: palindromic rs300 flagged ambiguous"

assert_contains "$ctsv" "Example condition" "clinvar: pathogenic hit reported (underscores->spaces)"
assert_contains "$ctsv" "homozygous"        "clinvar: zygosity reported"
assert_contains "$ctsv" "GENEX"             "clinvar: gene symbol reported"
assert_eq "0" "$(printf '%s\n' "$ctsv" | grep -c 'Nothing notable')" "clinvar: benign variant excluded"
assert_contains "$sum" "not a diagnosis"    "summary: carries research-grade caveat"

# Multi-condition ClinVar entries (CLNDN pipe-delimited) must not break the
# markdown table: pipes collapse to "; " and the list caps at 3 conditions.
assert_contains "$sum" "Cond one; Cond two; Cond three (+1 more)" "clinvar summary: pipe-delimited conditions collapsed + capped"
assert_eq "0" "$(printf '%s\n' "$sum" | grep -c 'Cond four')" "clinvar summary: 4th condition capped out of the rendered cell"
assert_eq "0" "$(printf '%s\n' "$sum" | grep -c 'Cond one|Cond two')" "clinvar summary: no raw pipe characters bleed into the table"

# Frequency context: raf is appended as GTSV column 11 and rendered as a
# percent in the summary GWAS table.
assert_contains "$gtsv" "GENEH	111	ok	0.23" "gwas tsv: raf appended after flag (col 11)"
assert_contains "$sum" "Risk-allele freq" "gwas summary: frequency column header present"
assert_contains "$sum" "| 23% |"          "gwas summary: raf rendered as rounded percent"

# ClinVar population frequency: AF_EXAC harvested with source, absent AF ->
# NA/unknown, buckets rendered in the summary with a plain-language bullet.
assert_contains "$ctsv" "ExAC"            "clinvar tsv: AF source recorded for the ExAC-backed hit"
assert_contains "$ctsv" "NA	NA"           "clinvar tsv: AF-less record -> pop_af NA, source NA"
assert_contains "$sum" "How common?"      "clinvar summary: frequency column header present"
assert_contains "$sum" "common (~30%)"    "clinvar summary: AF_EXAC=0.30 -> common bucket with percent"
assert_contains "$sum" "| unknown |"      "clinvar summary: AF-less hit shows unknown"
assert_contains "$sum" "large fraction of the population" "summary: plain-language common-variant caveat present"
