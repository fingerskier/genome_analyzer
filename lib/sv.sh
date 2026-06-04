# lib/sv.sh — Manta structural-variant analyzer. Sourced by analyze.sh.
# Expects globals: INPUT OUTDIR BASE GENES_BED ; helpers from common.sh.

: "${SCRIPT_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

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
  # Primary contigs (1-22, X, Y, MT/M) listed individually in karyotype order;
  # alt/random/decoy contigs collapsed into one "other contigs" row to keep the
  # table readable on real GRCh38 data (which has dozens of decoy contigs).
  perchrom="$(bcftools query -f '%CHROM\n' "$PASS" | sort | uniq -c \
    | awk '{
        cnt=$1; chr=$2; sub(/^chr/,"",chr)
        if(chr ~ /^([1-9]|1[0-9]|2[0-2]|X|Y|MT|M)$/){ prim[chr]=cnt }
        else { other_ev+=cnt; other_n++ }
      }
      END{
        n=split("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT M", ord, " ")
        for(i=1;i<=n;i++){ c=ord[i]; if(c in prim) printf "| %s | %s |\n", c, prim[c] }
        if(other_n>0) printf "| other contigs (%d) | %d |\n", other_n, other_ev
      }')"
  largest="$(bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$PASS" \
    | awk -F'\t' '$3!="BND" && $4!="."{v=($4<0?-$4:$4); print v"\t"$1":"$2"\t"$3}' \
    | sort -rn | awk -F'\t' 'NR<=10{printf "| %s | %s | %s |\n", $3, $2, $1}')"

  # Non-PASS filter reason tally.
  local reasons
  reasons="$(bcftools view -H "$INPUT" | awk -F'\t' '$7!="PASS" && $7!="."{print $7}' | sort | uniq -c | awk '{printf "%s (%s), ", $2, $1}' | sed 's/, $//')"
  : "${reasons:=none}"

  # Gene overlap (graceful: stats stand alone if bedtools / gene BED are absent).
  local gene_note bed="${GENES_BED:-$SCRIPT_DIR/data/genes.GRCh38.bed}"
  if ! have bedtools; then
    gene_note="_skipped: bedtools not found. \`brew install bedtools && ./fetch-genes.sh\`_"
  elif [[ ! -s "$bed" ]]; then
    gene_note="_skipped: no gene BED at \`$bed\`. Run \`./fetch-genes.sh\` (or pass \`-g\`)._"
  else
    # Build events BED. Non-BND: POS-1..END. BND: each endpoint as a 1bp interval
    # (this record's POS, plus the mate position parsed from ALT N[chr:pos[ / ]chr:pos]N).
    local ev="$OUTDIR/${BASE}.events.bed"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\t%ALT\n' "$PASS" \
      | awk -F'\t' 'BEGIN{OFS="\t"}
          $4!="BND" && $3!="."{ print $1, $2-1, $3, $5 }
          $4=="BND"{
            print $1, $2-1, $2, $5
            alt=$6
            if(match(alt, /[][][0-9XYMT]+:[0-9]+/)){
              loc=substr(alt, RSTART, RLENGTH); sub(/^[][]/,"",loc)
              split(loc, m, ":"); print m[1], m[2]-1, m[2], $5
            }
          }' > "$ev"
    gene_overlap "$ev" "$bed" "$GENES_TSV"
    gene_note="$(awk -F'\t' 'END{print NR}' "$GENES_TSV") events overlap genes — see \`$(basename "$GENES_TSV")\`."
  fi

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
