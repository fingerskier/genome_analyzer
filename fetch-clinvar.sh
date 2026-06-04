#!/usr/bin/env bash
# fetch-clinvar.sh — download NCBI ClinVar GRCh38 VCF (+ index) to data/.
# Contigs are bare names (1,2,...,X,Y,MT), matching Sequencing.com VCFs.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CLINVAR_URL="${CLINVAR_URL:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz}"

main() {
  local out="$SCRIPT_DIR/data/clinvar.GRCh38.vcf.gz"
  mkdir -p "$SCRIPT_DIR/data"
  if [[ -s "$out" && "${1:-}" != "-f" && "${1:-}" != "--force" ]]; then
    echo "already present: $out. Use -f to refresh." >&2
    return 0
  fi
  command -v curl  >/dev/null 2>&1 || { echo "ERROR: curl required"        >&2; exit 1; }
  command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix (htslib) required" >&2; exit 1; }
  echo "downloading $CLINVAR_URL" >&2
  curl -fsSL "$CLINVAR_URL" -o "$out"
  if ! curl -fsSL "${CLINVAR_URL}.tbi" -o "${out}.tbi" 2>/dev/null; then
    echo "no published .tbi; building index" >&2
    tabix -p vcf "$out"
  fi
  echo "wrote $out" >&2
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
