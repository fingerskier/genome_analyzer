#!/usr/bin/env bash
#
# analyze.sh — detect a single-sample VCF's class (snp-indel | sv | cnv) and
# dispatch to the matching analyzer. See docs/superpowers/specs/.
#
# Usage: ./analyze.sh -i sample.vcf.gz [-o out] [-g genes.bed] [-t N]
#                     [--type snp-indel|sv|cnv] [-a none|snpeff|vep] [-b GRCh38] [-r ref.fa]
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

# detect_type <vcf> -> echoes snp-indel|sv|cnv
detect_type() {
  local vcf="$1" hdr recs
  hdr="$(bcftools view -h "$vcf" 2>/dev/null)"
  recs="$(bcftools view -H "$vcf" 2>/dev/null | head -200)"
  if printf '%s' "$hdr" | grep -q '##source=Canvas' \
     || printf '%s' "$hdr" | grep -Eq '##ALT=<ID=CN[0-9]' \
     || printf '%s' "$recs" | grep -q 'SVTYPE=CNV'; then
    printf 'cnv'; return
  fi
  if printf '%s' "$hdr" | grep -Eqi 'Manta|GenerateSVCandidates' \
     || printf '%s' "$hdr" | grep -Eq '##ALT=<ID=(DEL|INS|DUP|INV|BND)' \
     || printf '%s' "$recs" | grep -Eq 'SVTYPE=(DEL|INS|DUP|INV|BND)'; then
    printf 'sv'; return
  fi
  printf 'snp-indel'
}

main() {
  INPUT=""; OUTDIR="out"; GENES_BED=""; THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
  FORCE_TYPE=""; ANNOTATE="none"; BUILD="GRCh38"; REF=""
  while getopts ":i:o:g:t:a:b:r:h-:" opt; do
    case "$opt" in
      i) INPUT="$OPTARG" ;;
      o) OUTDIR="$OPTARG" ;;
      g) GENES_BED="$OPTARG" ;;
      t) THREADS="$OPTARG" ;;
      a) ANNOTATE="$OPTARG" ;;
      b) BUILD="$OPTARG" ;;
      r) REF="$OPTARG" ;;
      h) sed -n '2,8p' "$0"; exit 0 ;;
      -) case "$OPTARG" in
           type) FORCE_TYPE="${!OPTIND}"; OPTIND=$((OPTIND+1)) ;;
           type=*) FORCE_TYPE="${OPTARG#*=}" ;;
           *) die "unknown option --$OPTARG" ;;
         esac ;;
      :) die "option -$OPTARG requires an argument" ;;
      \?) die "unknown option -$OPTARG" ;;
    esac
  done

  [[ -n "$INPUT" ]] || { sed -n '2,8p' "$0"; exit 1; }
  [[ -f "$INPUT" ]] || die "input not found: $INPUT"
  have bcftools || die "bcftools not on PATH (install htslib/bcftools first)"
  have tabix    || die "tabix not on PATH (part of htslib)"

  mkdir -p "$OUTDIR"
  BASE="$(basename "${INPUT%.vcf.gz}")"; BASE="${BASE%.vcf}"

  log "validate + index"
  validate_index "$INPUT"

  local TYPE="${FORCE_TYPE:-$(detect_type "$INPUT")}"
  log "detected: $TYPE"
  case "$TYPE" in
    snp-indel) source "$SCRIPT_DIR/lib/snp-indel.sh"; run_snp_indel ;;
    sv)        source "$SCRIPT_DIR/lib/sv.sh";        run_sv ;;
    cnv)       source "$SCRIPT_DIR/lib/cnv.sh";       run_cnv ;;
    *) die "unknown type: $TYPE (use --type snp-indel|sv|cnv)" ;;
  esac
}

# Allow sourcing for tests without running main.
if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
