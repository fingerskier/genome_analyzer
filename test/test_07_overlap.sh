if ! command -v bedtools >/dev/null 2>&1; then
  printf '  skip bedtools not installed; gene_overlap tested manually\n'
else
  source "$REPO/common.sh"
  # events BED: chrom \t start \t end \t id  (1-based POS-1 already applied by caller)
  printf '1\t4999\t5500\tD1\n2\t6999\t8500\tU1\n' > "$TMP/ev.bed"
  gene_overlap "$TMP/ev.bed" "$REPO/test/fixtures/genes.bed" "$TMP/ov.tsv"
  out="$(cat "$TMP/ov.tsv")"
  assert_contains "$out" "GENEA" "overlap: chr1 event hits GENEA"
  assert_contains "$out" "GENEB" "overlap: chr2 event hits GENEB"
fi
