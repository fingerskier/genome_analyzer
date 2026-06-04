# lib/sv.sh — Manta structural-variant analyzer. Sourced by analyze.sh.
# Expects globals: INPUT OUTDIR BASE GENES_BED ; helpers from common.sh.

run_sv() {
  local PASS="$OUTDIR/${BASE}.pass.vcf.gz"
  local SUMMARY="$OUTDIR/${BASE}.SUMMARY.md"
  local GENES_TSV="$OUTDIR/${BASE}.genes.tsv"

  log "PASS-only filter"
  bcftools view -f PASS,. "$INPUT" -Oz -o "$PASS"
  tabix -f -p vcf "$PASS"

  # Per-record: SVTYPE \t SVLEN \t EVENT \t ID   (PASS only)
  local Q
  Q="$(bcftools query -f '%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/EVENT\t%ID\n' "$PASS")"

  # Counts by type. BND counted as distinct EVENT (fallback: records/2 rounded up).
  local typecounts
  typecounts="$(printf '%s\n' "$Q" | awk -F'\t' '
    $1=="BND"{ ev=($3=="."?$4:$3); if(!(ev in seen)){seen[ev]=1; bnd++} next }
    { c[$1]++ }
    END{
      printf "DEL %d\n", c["DEL"]+0
      printf "DUP %d\n", c["DUP"]+0
      printf "INS %d\n", c["INS"]+0
      printf "INV %d\n", c["INV"]+0
      printf "BND %d\n", bnd+0
    }')"
  local n_del n_dup n_ins n_inv n_bnd
  n_del="$(awk '$1=="DEL"{print $2}' <<<"$typecounts")"
  n_dup="$(awk '$1=="DUP"{print $2}' <<<"$typecounts")"
  n_ins="$(awk '$1=="INS"{print $2}' <<<"$typecounts")"
  n_inv="$(awk '$1=="INV"{print $2}' <<<"$typecounts")"
  n_bnd="$(awk '$1=="BND"{print $2}' <<<"$typecounts")"

  # Size stats from abs(SVLEN), non-BND, non-empty.
  local maxlen medlen
  maxlen="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="BND" && $2!="."{v=($2<0?-$2:$2); if(v>m)m=v} END{print m+0}')"
  medlen="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="BND" && $2!="."{v=($2<0?-$2:$2); print v}' | sort -n | awk '{a[NR]=$1} END{if(NR==0){print 0} else if(NR%2){print a[(NR+1)/2]} else {print int((a[NR/2]+a[NR/2+1])/2)}}')"

  # Per-chromosome tally (PASS) and the largest non-BND events.
  local perchrom largest
  perchrom="$(bcftools query -f '%CHROM\n' "$PASS" | sort -V | uniq -c | awk '{printf "| %s | %s |\n", $2, $1}')"
  largest="$(bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$PASS" \
    | awk -F'\t' '$3!="BND" && $4!="."{v=($4<0?-$4:$4); print v"\t"$1":"$2"\t"$3}' \
    | sort -rn | head -10 | awk -F'\t' '{printf "| %s | %s | %s |\n", $3, $2, $1}')"

  # Non-PASS filter reason tally.
  local reasons
  reasons="$(bcftools view -H "$INPUT" | awk -F'\t' '$7!="PASS" && $7!="."{print $7}' | sort | uniq -c | awk '{printf "%s (%s), ", $2, $1}' | sed 's/, $//')"
  : "${reasons:=none}"

  # Gene overlap placeholder (wired in Task 9).
  local gene_note="_gene overlap added in a later step_"

  log "writing $SUMMARY"
  cat > "$SUMMARY" <<EOF
# SV Summary — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by analyze.sh (sv / Manta)_

## Events by type (PASS)
| Type | Count |
|---|---|
| DEL | $n_del |
| DUP | $n_dup |
| INS | $n_ins |
| INV | $n_inv |
| BND | $n_bnd |

## Size (non-BND, PASS)
| Metric | bp |
|---|---|
| Median | $medlen |
| Largest | $maxlen |

## Per chromosome (PASS)
| Chrom | Events |
|---|---|
$perchrom

## Largest events (non-BND)
| Type | Location | Size (bp) |
|---|---|---|
$largest

## Filtering
Non-PASS records dropped by: $reasons

## Gene overlap
$gene_note

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$PASS")\` | PASS-filtered SV VCF (+ index) |
| \`$(basename "$GENES_TSV")\` | Event → genes-hit table (when enabled) |
EOF

  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
