#!/usr/bin/env bash
# extract-triage.sh — digest a run's trait-xref hit tables into a compact,
# deterministic triage extract for the vcf-triage skill. Codifies the filter
# thresholds (frequency buckets, review-status weights, GWAS effect cutoff)
# so report generation never re-derives them ad hoc.
#
# Usage: extract-triage.sh <dir> <base>
#        extract-triage.sh <path-to-any-run-artifact>
# Prints the extract to stdout. Missing inputs are reported, never fatal.
set -uo pipefail

usage() { echo "Usage: $(basename "$0") <dir> <base> | <run-artifact-path>" >&2; exit 1; }

if [ $# -eq 2 ]; then
  DIR="$1"; BASE="$2"
elif [ $# -eq 1 ]; then
  DIR="$(dirname "$1")"; f="$(basename "$1")"
  case "$f" in
    *.traits.*)   BASE="${f%%.traits.*}" ;;
    *.vcf.gz)     BASE="${f%.vcf.gz}" ;;
    *.SUMMARY.md) BASE="${f%.SUMMARY.md}" ;;
    *)            BASE="$f" ;;
  esac
else
  usage
fi

CTSV="$DIR/$BASE.traits.clinvar.tsv"
GTSV="$DIR/$BASE.traits.gwas.tsv"
PTSV="$DIR/$BASE.traits.pharmgkb.tsv"
TAB="$(printf '\t')"
tmpd="$(mktemp -d)"; trap 'rm -rf "$tmpd"' EXIT

emit_inputs() {
  echo "== INPUTS =="
  if [ -f "$CTSV" ]; then
    echo "clinvar: found ($(wc -l < "$CTSV" | tr -d ' ') rows)"
  else
    echo "clinvar: MISSING — run xref-traits.sh (needs data/clinvar.GRCh38.vcf.gz via fetch-clinvar.sh)"
  fi
  if [ -f "$GTSV" ]; then
    echo "gwas: found ($(wc -l < "$GTSV" | tr -d ' ') rows)"
  else
    echo "gwas: MISSING — run xref-traits.sh (needs data/gwas-catalog.tsv via fetch-gwas.sh)"
  fi
  sums="$(ls "$DIR"/*.SUMMARY.md 2>/dev/null | grep -v '\.traits\.' || true)"
  if [ -n "$sums" ]; then
    echo "run summaries (read these directly):"
    printf '%s\n' "$sums" | sed 's/^/  /'
  else
    echo "run summaries: none found in $DIR (run analyze.sh)"
  fi
}

# Filled in by later tasks.
emit_clinvar_sections() { :; }
emit_gwas_sections() { :; }

emit_pharmgkb() {
  echo ""
  echo "== PHARMGKB =="
  if [ -f "$PTSV" ]; then
    echo "present ($(wc -l < "$PTSV" | tr -d ' ') rows) — include a drug-gene interaction section"
  else
    echo "absent — PharmGKB cross-reference not yet implemented (open roadmap item); skip that section"
  fi
}

emit_inputs
emit_clinvar_sections
emit_gwas_sections
emit_pharmgkb
