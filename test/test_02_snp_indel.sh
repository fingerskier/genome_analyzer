( source "$REPO/common.sh"
  source "$REPO/lib/snp-indel.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/snp-indel.vcf.txt")"
  OUTDIR="$TMP/snp_out"; mkdir -p "$OUTDIR"
  BASE="snp-indel"; THREADS=1; ANNOTATE=none; BUILD=GRCh38; REF=""
  run_snp_indel
)
sum="$TMP/snp_out/snp-indel.SUMMARY.md"
[[ -s "$sum" ]] && e=yes || e=no
assert_eq "yes" "$e" "snp-indel emits a non-empty SUMMARY.md"
assert_contains "$(cat "$sum")" "SNPs" "snp-indel summary reports SNP count"
