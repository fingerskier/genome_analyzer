( source "$REPO/common.sh"; source "$REPO/lib/cnv.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/cnv.vcf.txt")"
  OUTDIR="$TMP/cnv_out"; mkdir -p "$OUTDIR"; BASE="cnv"; THREADS=1; GENES_BED=""
  run_cnv )
sum="$(cat "$TMP/cnv_out/cnv.SUMMARY.md")"
assert_contains "$sum" "Losses | 2" "CNV summary: 2 PASS losses"
assert_contains "$sum" "Gains | 1" "CNV summary: 1 PASS gain"
assert_contains "$sum" "homozygous deletion" "CNV summary: flags CN0 homozygous deletion"
assert_contains "$sum" "L10kb" "CNV summary: notes filtered calls"
