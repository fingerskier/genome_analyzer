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

emit_clinvar_sections() {
  echo ""
  echo "== CLINVAR PRIORITY =="
  if [ ! -f "$CTSV" ]; then
    echo "(no clinvar table)"
    echo ""; echo "== CLINVAR PHARMA =="; echo "(no clinvar table)"
    echo ""; echo "== CLINVAR COMMON =="; echo "(no clinvar table)"
    return
  fi
  # One pass tags every row PRIORITY / PHARMA / COMMON / COUNTS; sections are
  # then grepped out of the temp file. Thresholds live here and nowhere else.
  awk -F"$TAB" '
    function bucket(af) {
      if (af == "" || af == "NA" || af !~ /^[0-9.eE+-]+$/) return "unknown"
      if (af+0 < 0.01) return "rare"
      if (af+0 < 0.05) return "low-frequency"
      return "common"
    }
    function weight(rev) {
      if (rev ~ /practice_guideline|reviewed_by_expert_panel/) return 4
      if (rev ~ /multiple_submitters/) return 3
      if (rev ~ /criteria_provided/) return 2
      return 1
    }
    {
      sig=$4; rev=$6; b=bucket($8)
      cond=$5; gsub(/\|/, "; ", cond); if (length(cond) > 100) cond = substr(cond, 1, 100) "..."
      afs=(b == "unknown" ? "AF unknown" : sprintf("AF %s (%s)", $8, $9))
      isPath=(sig ~ /athogenic/ && sig !~ /Conflicting/)
      if (isPath && (b == "rare" || b == "unknown"))
        printf "%d\t%.6f\tPRIORITY\t%s | %s | %s | %s | review: %s | %s\n", \
          weight(rev), (b == "unknown" ? 0.0099 : $8+0), $7, $3, sig, afs, rev, cond
      else if (sig ~ /athogenic/)
        printf "0\t0\tCOMMON\t%s | %s | %s\n", $7, sig, b
      if (sig ~ /drug_response/ && rev ~ /practice_guideline|reviewed_by_expert_panel/)
        printf "0\t0\tPHARMA\t%s | %s | %s | %s\n", $7, $3, b, cond
      n[(sig ~ /drug_response/) ? "drug" : (sig ~ /risk_factor/) ? "risk" : "path"]++
      total++
    }
    END {
      printf "0\t0\tCOUNTS\ttotal %d rows: drug_response %d, risk_factor %d, pathogenic-class/other %d\n", \
        total, n["drug"], n["risk"], n["path"]
    }
  ' "$CTSV" > "$tmpd/clinvar.tagged"

  prio="$(grep "${TAB}PRIORITY${TAB}" "$tmpd/clinvar.tagged" | sort -t"$TAB" -k1,1nr -k2,2g | cut -f4-)"
  if [ -n "$prio" ]; then printf '%s\n' "$prio"
  else echo "(none — no rare or unknown-frequency pathogenic-class hits)"; fi

  echo ""
  echo "== CLINVAR PHARMA =="
  ph="$(grep "${TAB}PHARMA${TAB}" "$tmpd/clinvar.tagged" | cut -f4-)"
  if [ -n "$ph" ]; then printf '%s\n' "$ph"
  else echo "(none — no expert-panel drug-response hits)"; fi

  echo ""
  echo "== CLINVAR COMMON =="
  grep "${TAB}COUNTS${TAB}" "$tmpd/clinvar.tagged" | cut -f4-
  co="$(grep "${TAB}COMMON${TAB}" "$tmpd/clinvar.tagged" | cut -f4-)"
  if [ -n "$co" ]; then printf '%s\n' "$co"
  else echo "(no common pathogenic-class rows to defuse)"; fi
}

emit_gwas_sections() {
  echo ""
  echo "== GWAS HEADLINE =="
  if [ ! -f "$GTSV" ]; then
    echo "(no gwas table)"
    echo ""; echo "== GWAS NOTABLE =="; echo "(no gwas table)"
    return
  fi
  awk -F"$TAB" '
    { t++; if ($5 == "2") h++; if ($10 == "ambiguous") a++ }
    END { printf "total %d associations; %d homozygous risk-allele; %d strand-ambiguous (not scored)\n", t, h, a }
  ' "$GTSV"

  echo ""
  echo "== GWAS NOTABLE =="
  echo "(flag ok or strand_flipped, OR >= 2.0, dedup by rsid keeping strongest p, cap 40)"
  # p-values like 1E-619 underflow doubles; compare exponent + log10(mantissa).
  nb="$(awk -F"$TAB" '
    function pstrength(p,  a, n, v) {
      n = split(toupper(p), a, "E")
      if (n > 1) return a[2]+0 + (a[1]+0 > 0 ? log(a[1]+0)/log(10) : 0)
      v = p+0
      return (v > 0 ? log(v)/log(10) : -99999)
    }
    ($10 == "ok" || $10 == "strand_flipped") && $7 ~ /^[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/ && $7+0 >= 2.0 {
      s = pstrength($6)
      if (!($1 in best) || s < bp[$1]) { best[$1] = $0; bp[$1] = s }
    }
    END { for (r in best) print best[r] }
  ' "$GTSV" | sort -t"$TAB" -k7,7gr | awk -F"$TAB" '
    NR <= 40 {
      raf = ($11 == "" || $11 == "NA" ? "unknown" : sprintf("%.0f%%", $11*100))
      printf "%s | %s | genotype %s | risk %s x%s | p %s | OR %s | %s | RAF %s\n", \
        $2, $8, $3, $4, $5, $6, $7, $1, raf
    }
    END { if (NR > 40) printf "... %d more notable hits not shown\n", NR-40 }
  ')"
  if [ -n "$nb" ]; then printf '%s\n' "$nb"
  else echo "(none — no unambiguous hits with OR >= 2.0)"; fi
}

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
