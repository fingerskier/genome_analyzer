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

main() {
  local out="$SCRIPT_DIR/data/genes.GRCh38.bed"
  mkdir -p "$SCRIPT_DIR/data"
  command -v bedtools >/dev/null 2>&1 || echo "note: bedtools not installed; analyzers need it for overlap" >&2
  echo "downloading $GENCODE_URL" >&2
  # One interval per gene symbol: keep the widest span seen for each name.
  curl -fsSL "$GENCODE_URL" | gunzip -c | gtf_to_bed \
    | sort -k4,4 -k1,1 -k2,2n \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        { if($4!=g){ if(g!=""){print c,s,e,g}; g=$4;c=$1;s=$2;e=$3 }
          else { if($2<s)s=$2; if($3>e)e=$3 } }
        END{ if(g!="") print c,s,e,g }' \
    | sort -k1,1 -k2,2n > "$out"
  echo "wrote $out ($(wc -l < "$out" | tr -d ' ') genes)" >&2
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
