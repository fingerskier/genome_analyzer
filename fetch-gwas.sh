#!/usr/bin/env bash
# fetch-gwas.sh — download GWAS Catalog associations, normalize to a compact
# join-ready table at data/gwas-catalog.tsv:
#   rsid risk trait pval orbeta gene pmid raf
# Only single-rsID rows with a single-base risk allele are kept (v1 scope).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# "ontology-annotated, alternative, full" associations zip. Override with GWAS_URL.
GWAS_URL="${GWAS_URL:-https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip}"

# normalize_gwas: GWAS Catalog associations TSV on stdin -> normalized table.
# Columns read (1-based): 2 PUBMEDID, 8 DISEASE/TRAIT, 15 MAPPED_GENE,
# 21 STRONGEST SNP-RISK ALLELE (e.g. "rs1234-A"), 27 RISK ALLELE FREQUENCY,
# 28 P-VALUE, 31 OR or BETA. Col 27 is sanity-checked against the header: on
# drift we warn and emit NA rather than harvest a wrong column.
normalize_gwas() {
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1{
      if($27 != "RISK ALLELE FREQUENCY"){
        badhdr=1
        print "WARNING: header col 27 is not RISK ALLELE FREQUENCY; raf column set to NA" > "/dev/stderr"
      }
      next
    }
    {
      sra=$21
      p=0; for(i=length(sra);i>=1;i--){ if(substr(sra,i,1)=="-"){p=i; break} }
      if(p==0) next
      rsid=substr(sra,1,p-1); risk=substr(sra,p+1)
      if(rsid !~ /^rs[0-9]+$/) next
      if(risk !~ /^[ACGT]$/)   next
      trait=$8; gsub(/[\t\r]/," ",trait)
      raf=$27; gsub(/\r/,"",raf)
      if(badhdr || raf !~ /^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ || raf+0>1) raf="NA"
      print rsid, risk, trait, $28, $31, $15, $2, raf
    }'
}

main() {
  local out="$SCRIPT_DIR/data/gwas-catalog.tsv"
  mkdir -p "$SCRIPT_DIR/data"
  if [[ -s "$out" && "${1:-}" != "-f" && "${1:-}" != "--force" ]]; then
    echo "already present: $out ($(wc -l < "$out" | tr -d ' ') rows). Use -f to refresh." >&2
    return 0
  fi
  command -v curl  >/dev/null 2>&1 || { echo "ERROR: curl required"  >&2; exit 1; }
  command -v unzip >/dev/null 2>&1 || { echo "ERROR: unzip required" >&2; exit 1; }
  echo "downloading $GWAS_URL" >&2
  local zip; zip="$(mktemp)"
  curl -fsSL "$GWAS_URL" -o "$zip"
  unzip -p "$zip" | normalize_gwas | LC_ALL=C sort -k1,1 > "$out"
  rm -f "$zip"
  echo "wrote $out ($(wc -l < "$out" | tr -d ' ') associations)" >&2
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
