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

# gene_overlap <events.bed> <genes.bed> <out.tsv>
# events.bed: chrom \t start0 \t end \t event_id  (0-based start, BED convention)
# Writes: event_id \t chrom \t start \t end \t n_genes \t gene1,gene2,...(+N more)
gene_overlap() {
  local events="$1" genes="$2" out="$3"
  have bedtools || die "gene_overlap called without bedtools"
  bedtools intersect -a "$events" -b "$genes" -wa -wb \
    | awk -F'\t' '
        { id=$4; key=id; chrom[id]=$1; s[id]=$2; e[id]=$3; g[id]=g[id] (g[id]?",":"") $8; n[id]++ }
        END{
          for(id in n){
            split(g[id], arr, ",")
            shown=""; cap=25
            for(i=1;i<=n[id] && i<=cap;i++) shown=shown (i>1?",":"") arr[i]
            if(n[id]>cap) shown=shown ",(+" (n[id]-cap) " more)"
            printf "%s\t%s\t%s\t%s\t%d\t%s\n", id, chrom[id], s[id], e[id], n[id], shown
          }
        }' | sort > "$out"
}
