# lib/snp-indel.sh — SNP/indel analyzer. Sourced by analyze.sh. Expects globals:
# INPUT OUTDIR BASE THREADS ANNOTATE BUILD REF ; helpers from common.sh.

run_snp_indel() {
  local STATS="$OUTDIR/${BASE}.stats.txt"
  local PASS="$OUTDIR/${BASE}.pass.vcf.gz"
  local ANNO="$OUTDIR/${BASE}.annotated.vcf.gz"
  local SUMMARY="$OUTDIR/${BASE}.SUMMARY.md"

  log "stats"
  bcftools stats "$INPUT" > "$STATS"
  local SAMPLES NSAMP RECORDS SNPS INDELS TITV
  SAMPLES="$(bcftools query -l "$INPUT" | paste -sd, -)"
  NSAMP="$(bcftools query -l "$INPUT" | wc -l | tr -d ' ')"
  get_sn() { awk -F'\t' -v k="$1" '$1=="SN" && $3==k {print $4; exit}' "$STATS"; }
  RECORDS="$(get_sn 'number of records:')"
  SNPS="$(get_sn 'number of SNPs:')"
  INDELS="$(get_sn 'number of indels:')"
  TITV="$(awk -F'\t' '$1=="TSTV"{print $5; exit}' "$STATS")"
  : "${RECORDS:=0}"; : "${SNPS:=0}"; : "${INDELS:=0}"; : "${TITV:=NA}"

  local ASSAY TITV_OK
  if   (( SNPS > 1000000 )); then ASSAY="Whole genome (WGS)"; TITV_OK="~2.0-2.1"
  elif (( SNPS > 20000   )); then ASSAY="Whole exome (WES) / large panel"; TITV_OK="~3.0-3.3"
  else                            ASSAY="Targeted panel / genotyping array"; TITV_OK="varies"
  fi
  log "assay looks like: $ASSAY  (SNPs=$SNPS, Ti/Tv=$TITV)"

  log "PASS-only filter"
  bcftools view -f PASS,. "$INPUT" -Oz -o "$PASS"
  tabix -f -p vcf "$PASS"
  local PASS_N
  PASS_N="$(bcftools index -n "$PASS" 2>/dev/null || bcftools stats "$PASS" | awk -F'\t' '$1=="SN" && $3=="number of records:"{print $4}')"

  local FINAL="$PASS" ANNO_NOTE="skipped (run with -a snpeff|vep)"
  case "$ANNOTATE" in
    none) log "annotation: skipped" ;;
    snpeff)
      have snpeff || die "snpeff not on PATH"
      log "SnpEff annotation ($BUILD)"
      snpeff -noStats "$BUILD" "$PASS" | bgzip > "$ANNO"
      tabix -f -p vcf "$ANNO"; FINAL="$ANNO"; ANNO_NOTE="SnpEff / $BUILD" ;;
    vep)
      have vep || die "vep not on PATH"
      log "VEP annotation ($BUILD)"
      local REF_ARG=(); [[ -n "$REF" ]] && REF_ARG=(--fasta "$REF" --hgvs)
      vep -i "$PASS" --cache --offline --assembly "$BUILD" --everything \
          --vcf --compress_output bgzip --fork "$THREADS" \
          "${REF_ARG[@]}" -o "$ANNO"
      tabix -f -p vcf "$ANNO"; FINAL="$ANNO"; ANNO_NOTE="Ensembl VEP / $BUILD" ;;
    *) die "unknown annotation engine: $ANNOTATE (use none|snpeff|vep)" ;;
  esac

  log "writing $SUMMARY"
  cat > "$SUMMARY" <<EOF
# VCF Summary — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by analyze.sh (snp-indel)_

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

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$STATS")\` | Full \`bcftools stats\` dump |
| \`$(basename "$PASS")\` | PASS-filtered VCF (+ index) |
$( [[ "$FINAL" == "$ANNO" ]] && echo "| \`$(basename "$ANNO")\` | Annotated VCF (+ index) |" )

## Next step
Hand \`$(basename "$FINAL")\` to the **Genome Analyzer** skill for SNPedia / GWAS
cross-referencing. Associations are research-grade, not clinical findings.
EOF

  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
