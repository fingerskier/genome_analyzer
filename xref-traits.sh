#!/usr/bin/env bash
# xref-traits.sh — cross-reference a SNP/indel PASS VCF against the GWAS Catalog
# and ClinVar, reporting carried trait/risk and pathogenic alleles with zygosity.
# Usage: ./xref-traits.sh -i sample.snp-indel.genome.pass.vcf.gz [-o out]
#                         [-g data/gwas-catalog.tsv] [-c data/clinvar.GRCh38.vcf.gz]
# Research-grade exploration, not a clinical diagnosis.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

usage() { sed -n '2,6p' "$0"; }

main() {
  INPUT=""; OUTDIR="out"
  GWAS_TBL="$SCRIPT_DIR/data/gwas-catalog.tsv"
  CLINVAR_VCF="$SCRIPT_DIR/data/clinvar.GRCh38.vcf.gz"
  while getopts ":i:o:g:c:h" opt; do
    case "$opt" in
      i) INPUT="$OPTARG" ;;
      o) OUTDIR="$OPTARG" ;;
      g) GWAS_TBL="$OPTARG" ;;
      c) CLINVAR_VCF="$OPTARG" ;;
      h) usage; exit 0 ;;
      :) die "option -$OPTARG requires an argument" ;;
      \?) die "unknown option -$OPTARG" ;;
    esac
  done
  [[ -n "$INPUT" ]] || { usage; exit 1; }
  [[ -f "$INPUT" ]] || die "input not found: $INPUT"
  have bcftools || die "bcftools not on PATH (install htslib/bcftools first)"
  have tabix    || die "tabix not on PATH (part of htslib)"

  mkdir -p "$OUTDIR"
  BASE="$(basename "${INPUT%.vcf.gz}")"; BASE="${BASE%.vcf}"
  log "validate + index"
  validate_index "$INPUT"
  # shellcheck disable=SC1091
  source "$SCRIPT_DIR/lib/trait-xref.sh"
  run_trait_xref
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
