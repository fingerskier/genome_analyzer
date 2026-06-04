# common.sh — shared helpers for the genome analyzers. Sourced, never executed.
# Targets macOS bash 3.2: no printf '%(...)T', no associative arrays, no mapfile.

log()  { printf '\033[1;36m[%s]\033[0m %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
die()  { printf '\033[1;31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }
have() { command -v "$1" >/dev/null 2>&1; }

# validate_index <vcf>: confirm readable bgzipped VCF; create .tbi if no index.
validate_index() {
  local vcf="$1"
  if ! bcftools view "$vcf" -h >/dev/null 2>&1; then
    die "$vcf is not a readable VCF (bgzipped, not plain gzip?). Try: zcat in.vcf.gz | bgzip > out.vcf.gz"
  fi
  if [[ ! -f "${vcf}.tbi" && ! -f "${vcf}.csi" ]]; then
    log "  no index found, creating .tbi"
    tabix -p vcf "$vcf"
  fi
}
