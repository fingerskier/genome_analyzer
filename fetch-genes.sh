#!/usr/bin/env bash
# fetch-genes.sh — download GENCODE basic GRCh38 gene spans, collapse to one
# interval per gene symbol, normalize chr-prefix, cache to data/genes.GRCh38.bed.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# GENCODE basic annotation (release 45, GRCh38). Override with GENCODE_URL env var.
GENCODE_URL="${GENCODE_URL:-https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz}"

# gtf_to_bed: GENCODE GTF on stdin -> 'chrom<TAB>start0<TAB>end<TAB>gene_name', chr stripped.
gtf_to_bed() {
  awk -F'\t' 'BEGIN{OFS="\t"}
    $0 ~ /^#/ {next}
    $3=="gene"{
      name=""
      if(match($9, /gene_name "[^"]+"/)){ name=substr($9, RSTART+11, RLENGTH-12) }
      if(name==""){ next }
      chrom=$1; sub(/^chr/, "", chrom)
      print chrom, $4-1, $5, name
    }'
}

# collapse_genes: BED4 'chrom<TAB>start<TAB>end<TAB>name' on stdin -> one interval
# per (gene name, chromosome), widest span, genomic-sorted. Keying on name AND
# chromosome keeps pseudo-autosomal genes (e.g. SHOX on chrX and chrY) separate.
collapse_genes() {
  sort -k4,4 -k1,1 -k2,2n \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        { if($4!=g || $1!=gc){ if(g!=""){print c,s,e,g}; g=$4; gc=$1; c=$1; s=$2; e=$3 }
          else { if($2<s)s=$2; if($3>e)e=$3 } }
        END{ if(g!="") print c,s,e,g }' \
    | sort -k1,1 -k2,2n
}

main() {
  local out="$SCRIPT_DIR/data/genes.GRCh38.bed"
  mkdir -p "$SCRIPT_DIR/data"
  command -v bedtools >/dev/null 2>&1 || echo "note: bedtools not installed; analyzers need it for overlap" >&2
  echo "downloading $GENCODE_URL" >&2
  curl -fsSL "$GENCODE_URL" | gunzip -c | gtf_to_bed | collapse_genes > "$out"
  echo "wrote $out ($(wc -l < "$out" | tr -d ' ') genes)" >&2
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
