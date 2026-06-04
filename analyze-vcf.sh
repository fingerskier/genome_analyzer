#!/usr/bin/env bash
#
# analyze-vcf.sh — interrogate a single-sample germline VCF and emit a tidy summary.
#
# Pipeline:  validate/index -> stats -> assay detection -> PASS filter -> [annotate] -> SUMMARY.md
#
# The core path (stats + filter + summary) has only one hard dependency: bcftools (+ tabix/htslib).
# Annotation is opt-in because VEP/SnpEff carry heavy cache/database setup.
#
# Usage:
#   ./analyze-vcf.sh -i sample.vcf.gz [-o out_dir] [-a none|snpeff|vep] [-b GRCh38] [-r ref.fa]
#
# Examples:
#   ./analyze-vcf.sh -i her.vcf.gz                       # stats + PASS filter + summary
#   ./analyze-vcf.sh -i her.vcf.gz -a snpeff -b GRCh38   # + SnpEff annotation
#
set -euo pipefail

# ---- defaults ---------------------------------------------------------------
INPUT=""
OUTDIR="out"
ANNOTATE="none"          # none | snpeff | vep
BUILD="GRCh38"           # genome build for annotation DB / VEP cache
REF=""                   # optional reference FASTA (VEP can use it for HGVS)
THREADS="$(nproc 2>/dev/null || echo 4)"

# ---- helpers ----------------------------------------------------------------
log()  { printf '\033[1;36m[%(%H:%M:%S)T]\033[0m %s\n' -1 "$*" >&2; }
die()  { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

usage() { sed -n '2,20p' "$0"; exit "${1:-0}"; }

# ---- arg parsing ------------------------------------------------------------
while getopts ":i:o:a:b:r:t:h" opt; do
  case "$opt" in
    i) INPUT="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    a) ANNOTATE="$OPTARG" ;;
    b) BUILD="$OPTARG" ;;
    r) REF="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage 0 ;;
    :) die "option -$OPTARG requires an argument" ;;
    \?) die "unknown option -$OPTARG" ;;
  esac
done

[[ -n "$INPUT" ]] || { usage 1; }
[[ -f "$INPUT" ]] || die "input not found: $INPUT"
have bcftools || die "bcftools not on PATH (install htslib/bcftools first)"
have tabix    || die "tabix not on PATH (part of htslib)"

mkdir -p "$OUTDIR"
BASE="$(basename "${INPUT%.vcf.gz}")"
BASE="${BASE%.vcf}"
STATS="$OUTDIR/${BASE}.stats.txt"
PASS="$OUTDIR/${BASE}.pass.vcf.gz"
ANNO="$OUTDIR/${BASE}.annotated.vcf.gz"
SUMMARY="$OUTDIR/${BASE}.SUMMARY.md"

# ---- stage 1: validate + index --------------------------------------------
log "Stage 1/5  validate + index"
# bgzip sanity: gzip-but-not-bgzip VCFs break tabix. quickcheck catches truncation.
if ! bcftools view "$INPUT" -h >/dev/null 2>&1; then
  die "$INPUT is not a readable VCF (is it bgzipped, not plain gzip?). Try: zcat in.vcf.gz | bgzip > out.vcf.gz"
fi
if [[ ! -f "${INPUT}.tbi" && ! -f "${INPUT}.csi" ]]; then
  log "  no index found, creating .tbi"
  tabix -p vcf "$INPUT"
fi

# ---- stage 2: stats ---------------------------------------------------------
log "Stage 2/5  bcftools stats"
bcftools stats "$INPUT" > "$STATS"

SAMPLES="$(bcftools query -l "$INPUT" | paste -sd, -)"
NSAMP="$(bcftools query -l "$INPUT" | wc -l | tr -d ' ')"
get_sn() { awk -F'\t' -v k="$1" '$1=="SN" && $3==k {print $4; exit}' "$STATS"; }
RECORDS="$(get_sn 'number of records:')"
SNPS="$(get_sn 'number of SNPs:')"
INDELS="$(get_sn 'number of indels:')"
MNPS="$(get_sn 'number of MNPs:')"
TITV="$(awk -F'\t' '$1=="TSTV"{print $5; exit}' "$STATS")"
: "${RECORDS:=0}"; : "${SNPS:=0}"; : "${INDELS:=0}"; : "${TITV:=NA}"

# ---- stage 3: assay detection ----------------------------------------------
# Heuristic on SNP count. Rough but reliable enough to set expectations.
if   (( SNPS > 1000000 )); then ASSAY="Whole genome (WGS)"; TITV_OK="~2.0–2.1"
elif (( SNPS > 20000   )); then ASSAY="Whole exome (WES) / large panel"; TITV_OK="~3.0–3.3"
else                            ASSAY="Targeted panel / genotyping array"; TITV_OK="varies"
fi
log "Stage 3/5  assay looks like: $ASSAY  (SNPs=$SNPS, Ti/Tv=$TITV)"

# ---- stage 4: PASS filter ---------------------------------------------------
log "Stage 4/5  PASS-only filter"
# Keep records where FILTER is PASS or '.'(unfiltered). Drop everything else.
bcftools view -f PASS,. "$INPUT" -Oz -o "$PASS"
tabix -f -p vcf "$PASS"
PASS_N="$(bcftools index -n "$PASS" 2>/dev/null || bcftools stats "$PASS" | awk -F'\t' '$1=="SN" && $3=="number of records:"{print $4}')"

# ---- stage 5: annotation (optional) ----------------------------------------
FINAL="$PASS"
ANNO_NOTE="skipped (run with -a snpeff|vep)"
case "$ANNOTATE" in
  none) log "Stage 5/5  annotation: skipped" ;;
  snpeff)
    have snpeff || die "snpeff not on PATH"
    log "Stage 5/5  SnpEff annotation ($BUILD)"
    snpeff -noStats "$BUILD" "$PASS" | bgzip > "$ANNO"
    tabix -f -p vcf "$ANNO"; FINAL="$ANNO"; ANNO_NOTE="SnpEff / $BUILD" ;;
  vep)
    have vep || die "vep not on PATH"
    log "Stage 5/5  VEP annotation ($BUILD)"
    REF_ARG=(); [[ -n "$REF" ]] && REF_ARG=(--fasta "$REF" --hgvs)
    vep -i "$PASS" --cache --offline --assembly "$BUILD" --everything \
        --vcf --compress_output bgzip --fork "$THREADS" \
        "${REF_ARG[@]}" -o "$ANNO"
    tabix -f -p vcf "$ANNO"; FINAL="$ANNO"; ANNO_NOTE="Ensembl VEP / $BUILD" ;;
  *) die "unknown annotation engine: $ANNOTATE (use none|snpeff|vep)" ;;
esac

# ---- emit summary -----------------------------------------------------------
log "Writing $SUMMARY"
cat > "$SUMMARY" <<EOF
# VCF Summary — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by analyze-vcf.sh_

## Overview
| Field | Value |
|---|---|
| Source file | \`$(basename "$INPUT")\` |
| Sample(s) | $NSAMP — \`${SAMPLES}\` |
| Likely assay | **$ASSAY** |
| Total records | $RECORDS |
| SNPs | $SNPS |
| Indels | $INDELS |
| Ti/Tv ratio | **$TITV** (healthy: $TITV_OK) |
| PASS records | $PASS_N |
| Annotation | $ANNO_NOTE |

## Quality read
- **Ti/Tv** $TITV against an expected $TITV_OK for this assay type. A value well below range hints at false-positive enrichment; treat downstream calls with more skepticism.
- $(( PASS_N < RECORDS )) && echo "Filtering removed $(( RECORDS - PASS_N )) non-PASS records." || echo "All records were PASS/unfiltered."

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$STATS")\` | Full \`bcftools stats\` dump |
| \`$(basename "$PASS")\` | PASS-filtered VCF (+ index) |
$( [[ "$FINAL" == "$ANNO" ]] && echo "| \`$(basename "$ANNO")\` | Annotated VCF (+ index) |" )

## Next step
Hand \`$(basename "$FINAL")\` to the **Genome Analyzer** Claude Code skill for SNPedia / GWAS
cross-referencing, or query specific variants via **BioMCP**. Remember: associations are
research-grade, not clinical findings.
EOF

log "Done. Open: $SUMMARY"
echo "$SUMMARY"
