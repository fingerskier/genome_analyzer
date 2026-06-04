# test_08_degrade.sh — gene-overlap: degradation + positive overlap paths.

# (1) Degradation: bedtools masked via have() override (PATH-masking is unreliable on
# macOS because setting PATH="nobin" inside a subshell blocks mkdir/ln before the nobin
# dir can be populated). The plan's documented fallback: override have() to return 1.
# SV must still emit a summary whose Gene-overlap section mentions "bedtools".
( source "$REPO/common.sh"; source "$REPO/lib/sv.sh"
  have(){ return 1; }
  INPUT="$(make_fixture "$REPO/test/fixtures/sv.vcf.txt")"
  OUTDIR="$TMP/deg_out"; mkdir -p "$OUTDIR"; BASE="sv"; GENES_BED=""
  run_sv ) >/dev/null 2>&1
assert_contains "$(cat "$TMP/deg_out/sv.SUMMARY.md")" "bedtools" "degraded SV summary explains how to enable overlap"

# (2) Positive SV overlap: with bedtools + the fixture gene BED, run_sv writes a
# populated genes.tsv. GENEA proves chr1 overlap; GENEE proves BND endpoint parsing.
( source "$REPO/common.sh"; source "$REPO/lib/sv.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/sv.vcf.txt")"
  OUTDIR="$TMP/sv_ov"; mkdir -p "$OUTDIR"; BASE="sv"; GENES_BED="$REPO/test/fixtures/genes.bed"
  run_sv ) >/dev/null 2>&1
svg="$(cat "$TMP/sv_ov/sv.genes.tsv" 2>/dev/null)"
assert_contains "$svg" "GENEA" "SV overlap: chr1 events hit GENEA"
assert_contains "$svg" "GENEE" "SV overlap: BND chr5 endpoint hits GENEE (both-endpoints)"

# (3) Positive CNV overlap.
( source "$REPO/common.sh"; source "$REPO/lib/cnv.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/cnv.vcf.txt")"
  OUTDIR="$TMP/cnv_ov"; mkdir -p "$OUTDIR"; BASE="cnv"; GENES_BED="$REPO/test/fixtures/genes.bed"
  run_cnv ) >/dev/null 2>&1
cng="$(cat "$TMP/cnv_ov/cnv.genes.tsv" 2>/dev/null)"
assert_contains "$cng" "GENEC" "CNV overlap: chr3 gain hits GENEC"
