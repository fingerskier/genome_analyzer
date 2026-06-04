# lib/cnv.sh — Canvas copy-number-variant analyzer. Sourced by analyze.sh.
# Expects globals: INPUT OUTDIR BASE GENES_BED ; helpers from common.sh.

: "${SCRIPT_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

run_cnv() {
  local PASS="$OUTDIR/${BASE}.pass.vcf.gz"
  local SUMMARY="$OUTDIR/${BASE}.SUMMARY.md"
  local GENES_TSV="$OUTDIR/${BASE}.genes.tsv"

  # Real events only: Canvas:REF segments carry no SVTYPE and ALT='.'; exclude them.
  # Keep records whose ALT is a symbolic <...>. Then PASS-filter.
  log "drop reference segments + PASS filter"
  bcftools view -f PASS,. "$INPUT" \
    | awk -F'\t' '/^#/ || $5 ~ /^</' \
    | bcftools view -Oz -o "$PASS"
  tabix -f -p vcf "$PASS"

  # Per-event: CN \t CNVLEN
  local Q
  Q="$(bcftools query -f '[%CN]\t%INFO/CNVLEN\n' "$PASS")"

  local n_loss n_gain n_homdel loss_bp gain_bp
  n_loss="$(printf '%s\n' "$Q"  | awk -F'\t' '$1!="." && $1<2' | grep -c . || true)"
  n_gain="$(printf '%s\n' "$Q"  | awk -F'\t' '$1!="." && $1>2' | grep -c . || true)"
  n_homdel="$(printf '%s\n' "$Q"| awk -F'\t' '$1=="0"' | grep -c . || true)"
  loss_bp="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="." && $1<2 && $2!="."{s+=$2} END{print s+0}')"
  gain_bp="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="." && $1>2 && $2!="."{s+=$2} END{print s+0}')"
  local loss_mb gain_mb
  loss_mb="$(awk -v b="$loss_bp" 'BEGIN{printf "%.2f", b/1000000}')"
  gain_mb="$(awk -v b="$gain_bp" 'BEGIN{printf "%.2f", b/1000000}')"

  local reasons
  reasons="$(bcftools view -H "$INPUT" | awk -F'\t' '$5 ~ /^</ && $7!="PASS" && $7!="."{print $7}' | sort | uniq -c | awk '{printf "%s (%s), ", $2, $1}' | sed 's/, $//')"
  : "${reasons:=none}"

  # Per-chromosome tally (PASS events) and the largest CNVs by spanned length.
  local perchrom largest
  perchrom="$(bcftools query -f '%CHROM\n' "$PASS" | sort -V | uniq -c | awk '{printf "| %s | %s |\n", $2, $1}')"
  largest="$(bcftools query -f '%CHROM\t%POS\t[%CN]\t%INFO/CNVLEN\n' "$PASS" \
    | awk -F'\t' '$4!="."{dir=($3<2?"loss":"gain"); print $4"\t"$1":"$2"\t"dir" CN="$3}' \
    | sort -rn | awk -F'\t' 'NR<=10{printf "| %s | %s | %s bp |\n", $3, $2, $1}')"

  local homdel_note=""
  (( n_homdel > 0 )) && homdel_note=" (incl. $n_homdel homozygous deletion(s), CN=0)"

  local gene_note bed="${GENES_BED:-$SCRIPT_DIR/data/genes.GRCh38.bed}"
  if ! have bedtools; then
    gene_note="_skipped: bedtools not found. \`brew install bedtools && ./fetch-genes.sh\`_"
  elif [[ ! -s "$bed" ]]; then
    gene_note="_skipped: no gene BED at \`$bed\`. Run \`./fetch-genes.sh\` (or pass \`-g\`)._"
  else
    local ev="$OUTDIR/${BASE}.events.bed"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' "$PASS" \
      | awk -F'\t' 'BEGIN{OFS="\t"} $3!="."{print $1, $2-1, $3, $4}' > "$ev"
    gene_overlap "$ev" "$bed" "$GENES_TSV"
    gene_note="$(awk -F'\t' 'END{print NR}' "$GENES_TSV") events overlap genes — see \`$(basename "$GENES_TSV")\`."
  fi

  log "writing $SUMMARY"
  cat > "$SUMMARY" <<EOF
# CNV Summary — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by analyze.sh (cnv / Canvas)_

## Events (PASS, reference segments excluded)
| Class | Count |
|---|---|
| Losses | $n_loss$homdel_note |
| Gains | $n_gain |

## Genome burden (PASS)
| Direction | Mb |
|---|---|
| Lost | $loss_mb |
| Gained | $gain_mb |

## Per chromosome (PASS)
| Chrom | Events |
|---|---|
$perchrom

## Largest events
| Direction | Location | Size (bp) |
|---|---|---|
$largest

## Filtering
Non-PASS events dropped by: $reasons

## Gene overlap
$gene_note

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$PASS")\` | PASS-filtered CNV VCF, reference segments removed (+ index) |
| \`$(basename "$GENES_TSV")\` | Event → genes-hit table (when enabled) |
EOF

  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
