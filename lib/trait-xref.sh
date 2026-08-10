# lib/trait-xref.sh — trait cross-reference engine. Sourced by xref-traits.sh.
# Expects globals (for run_trait_xref): INPUT OUTDIR BASE GWAS_TBL CLINVAR_VCF.
# Targets macOS bash 3.2: no associative arrays, no mapfile.

: "${SCRIPT_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"

# Single source of truth for allele dosage, reused by risk_dosage() (unit-tested)
# and by the GWAS join in run_trait_xref(). dosage() returns "<n><SUBSEP><flag>"
# where n is 0/1/2 or "NA". A het genotype carries exactly one copy of whichever
# of {ref,alt} the (possibly strand-flipped) risk allele equals; a hom-alt
# genotype carries 2 copies of alt, 0 of ref. Palindromic SNPs (A/T, C/G) cannot
# be strand-resolved, so they are flagged "ambiguous" and not scored.
DOSAGE_FN='
function comp(b){ return b=="A"?"T":b=="T"?"A":b=="C"?"G":b=="G"?"C":"N" }
function dosage(ref,alt,gt,risk,   rc){
  if(length(ref)==1 && length(alt)==1 && comp(ref)==alt){ return "NA" SUBSEP "ambiguous" }
  if(risk==ref || risk==alt){ return (gt=="hom" ? (risk==alt?2:0) : 1) SUBSEP "ok" }
  if(length(ref)==1 && length(alt)==1 && length(risk)==1){
    rc=comp(risk)
    if(rc==ref || rc==alt){ return (gt=="hom" ? (rc==alt?2:0) : 1) SUBSEP "strand_flipped" }
  }
  return "NA" SUBSEP "allele_mismatch"
}'

# risk_dosage <ref> <alt> <het|hom> <risk> -> "<dosage>\t<flag>"
risk_dosage() {
  awk -v A="$1" -v B="$2" -v G="$3" -v R="$4" "$DOSAGE_FN"'
    BEGIN{ split(dosage(A,B,G,R), p, SUBSEP); print p[1]"\t"p[2]; exit }'
}

# Population-frequency buckets for triage context — single source of truth,
# shared (like DOSAGE_FN) by freq_bucket() and the summary rendering.
# bucket() returns "<bucket><SUBSEP><display>": common >=5%, low-frequency
# 1-5%, rare <1%, unknown for absent/non-numeric/impossible values.
FREQ_FN='
function bucket(af,   pct){
  if(af !~ /^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ || af+0>1) return "unknown" SUBSEP "unknown"
  pct=(af+0)*100
  if(af+0>=0.05) return "common" SUBSEP sprintf("common (~%.0f%%)", pct)
  if(af+0>=0.01) return "low-frequency" SUBSEP sprintf("low-frequency (~%.0f%%)", pct)
  return "rare" SUBSEP "rare (<1%)"
}'

# freq_bucket <af> -> "<bucket>\t<display>"
freq_bucket() {
  awk -v A="$1" "$FREQ_FN"'
    BEGIN{ split(bucket(A), p, SUBSEP); print p[1]"\t"p[2]; exit }'
}

# run_trait_xref: cross-reference the sample VCF ($INPUT) against the GWAS Catalog
# ($GWAS_TBL) and ClinVar ($CLINVAR_VCF). Writes three artifacts under $OUTDIR and
# prints the SUMMARY path. Degrades gracefully when a database file is absent.
run_trait_xref() {
  local SUMMARY="$OUTDIR/${BASE}.traits.SUMMARY.md"
  local GTSV="$OUTDIR/${BASE}.traits.gwas.tsv"
  local CTSV="$OUTDIR/${BASE}.traits.clinvar.tsv"
  local tmp; tmp="$(mktemp -d)"
  # shellcheck disable=SC2064  # intentional: expand $tmp now so the trap works
  # even on bash 3.2 where RETURN traps leak into the caller's scope.
  trap "rm -rf $(printf %q "$tmp")" RETURN

  # This tool is single-sample by design; the query below reads the first sample
  # only. Warn (don't fail) if the VCF carries more, so partial results aren't a
  # silent surprise.
  local nsamp; nsamp="$(bcftools query -l "$INPUT" 2>/dev/null | grep -c . || true)"
  if [[ "${nsamp:-0}" -gt 1 ]]; then
    log "warning: $nsamp samples present; only the first is cross-referenced"
  fi

  # 1) One normalize+query pass -> carried biallelic calls.
  #    carried.tsv: chrom \t pos \t id \t ref \t alt \t zyg(het|hom)
  log "extracting carried variants (norm -m- + query)"
  bcftools norm -m- "$INPUT" 2>/dev/null \
    | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        { gt=$6; m=gsub(/1/,"1",gt); if(m==0) next; print $1,$2,$3,$4,$5,(m>=2?"hom":"het") }' \
    > "$tmp/carried.tsv"

  # 2) GWAS join (by rsID, p < 5e-8). Load significant rows keyed by rsid, then
  #    stream the sample's carried calls and score dosage with the shared dosage().
  local gwas_note
  if [[ ! -s "$GWAS_TBL" ]]; then
    gwas_note="_skipped: no GWAS table at \`$GWAS_TBL\`. Run \`./fetch-gwas.sh\`._"
    : > "$GTSV"
  else
    log "GWAS Catalog cross-reference (p < 5e-8)"
    awk -F'\t' "$DOSAGE_FN"'
      BEGIN{OFS="\t"; THRESH=5e-8}
      FNR==NR{
        if($4 ~ /[0-9]/ && ($4+0)<THRESH){
          rec=$2 SUBSEP $3 SUBSEP $4 SUBSEP $5 SUBSEP $6 SUBSEP $7
          g[$1]=($1 in g)? g[$1] "\x1f" rec : rec
        }
        next
      }
      ($3 in g){
        rs=$3; ref=$4; alt=$5; z=$6
        geno=(z=="hom"? alt"/"alt : ref"/"alt)
        n=split(g[rs], aa, "\x1f")
        for(i=1;i<=n;i++){
          split(aa[i], f, SUBSEP)
          split(dosage(ref,alt,z,f[1]), dd, SUBSEP)
          dose=dd[1]; flag=dd[2]
          if(flag=="allele_mismatch") continue
          if(flag=="ambiguous" || dose+0>=1)
            print rs, f[2], geno, f[1], dose, f[3], f[4], f[5], f[6], flag
        }
      }' "$GWAS_TBL" "$tmp/carried.tsv" \
      | LC_ALL=C sort -t"$(printf '\t')" -k2,2 > "$GTSV"
    local ghits homhits ambn
    ghits="$(awk -F'\t'  '$10!="ambiguous"' "$GTSV" | grep -c . || true)"
    homhits="$(awk -F'\t' '$5==2' "$GTSV" | grep -c . || true)"
    ambn="$(awk -F'\t'   '$10=="ambiguous"' "$GTSV" | grep -c . || true)"
    gwas_note="$ghits trait association(s) carried ($homhits homozygous; $ambn strand-ambiguous, not scored) — see \`$(basename "$GTSV")\`."
  fi

  # 3) ClinVar join (position+allele, pathogenic-class only).
  local clinvar_note
  if [[ ! -s "$CLINVAR_VCF" ]]; then
    clinvar_note="_skipped: no ClinVar VCF at \`$CLINVAR_VCF\`. Run \`./fetch-clinvar.sh\`._"
    : > "$CTSV"
  else
    log "ClinVar cross-reference (pathogenic-class)"
    bcftools query \
      -i 'CLNSIG ~ "Pathogenic" || CLNSIG ~ "Likely_pathogenic" || CLNSIG ~ "risk_factor" || CLNSIG ~ "drug_response"' \
      -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\n' \
      "$CLINVAR_VCF" 2>/dev/null \
      | awk -F'\t' 'BEGIN{OFS="\t"}
          FNR==NR{
            k=$1":"$2":"$3":"$4
            g=$8; gsub(/:[0-9]+/,"",g); gsub(/\|/,", ",g)
            sig[k]=$5; dn[k]=$6; rev[k]=$7; gene[k]=g; next
          }
          { k=$1":"$2":"$4":"$5
            if(k in sig){
              cond=dn[k]; gsub(/_/," ",cond)
              print $1":"$2, $5, ($6=="hom"?"homozygous":"heterozygous"), sig[k], cond, rev[k], gene[k]
            }
          }' - "$tmp/carried.tsv" > "$CTSV"
    local chits
    chits="$(grep -c . "$CTSV" || true)"
    clinvar_note="$chits pathogenic-class variant(s) carried — see \`$(basename "$CTSV")\`."
  fi

  # 4) Build SUMMARY tables (capped for readability; full data in the TSVs).
  local clinvar_rows gwas_rows
  # CLNDN/CLNSIG pack multiple values with '|', which would collide with the
  # markdown column delimiter — collapse to "; " and cap the condition list at 3.
  clinvar_rows="$(awk -F'\t' 'NR<=50{
      sig=$4; gene=$7; gsub(/\|/,"; ",sig); gsub(/\|/,"; ",gene)
      n=split($5, cc, /\|/); cond=""
      for(i=1;i<=n && i<=3;i++) cond=cond (i>1?"; ":"") cc[i]
      if(n>3) cond=cond " (+" (n-3) " more)"
      printf "| %s | %s | %s | %s | %s |\n", cond, gene, $3, sig, $6
    }' "$CTSV")"
  : "${clinvar_rows:=| _none carried_ |  |  |  |  |}"
  gwas_rows="$(awk -F'\t' '$10!="ambiguous"' "$GTSV" \
    | LC_ALL=C sort -t"$(printf '\t')" -k6,6g \
    | awk -F'\t' 'NR<=50{hl=($5==2?"**":""); printf "| %s | %s | %s%s%s | %s | %s | %s | %s |\n", $2,$8,hl,$3,hl,$4,$5,$6,$7}')"
  : "${gwas_rows:=| _none carried_ |  |  |  |  |  |  |}"

  log "writing $SUMMARY"
  cat > "$SUMMARY" <<EOF
# Trait Cross-Reference — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by xref-traits.sh. **Research-grade exploration, not a diagnosis.**_

## ClinVar — pathogenic-class variants you carry
$clinvar_note

| Condition | Gene | Your call | Significance | Review status |
|---|---|---|---|---|
$clinvar_rows

## GWAS Catalog — trait associations you carry (p < 5e-8)
$gwas_note

| Trait | Gene | Your genotype | Risk allele | Copies | p-value | OR/beta |
|---|---|---|---|---|---|---|
$gwas_rows

(**Bold** genotype = homozygous for the risk allele. "Copies" is how many of your
two alleles match the catalogued risk allele.)

## How to read this
- These are **statistical associations and clinical annotations from public
  databases**, not a diagnosis or medical advice. Carrying a risk allele does not
  mean you have or will get a condition.
- **Strand-ambiguous SNPs** (A/T or C/G) cannot be reliably oriented from this VCF
  alone, so they are listed in the table only with an \`ambiguous\` flag and are not
  scored for copies.
- ClinVar **review status** indicates how well-supported a classification is
  (more submitters / expert panels = stronger).
- This input is a genome VCF with reference blocks, so a known risk site you are
  *not* listed at was genuinely called homozygous-reference (0 copies), not missing.

## Future work
- Lower the GWAS p-value threshold or include sub-significant / VUS entries.
- Polygenic risk scores; haplotype (multi-SNP) associations; ancestry weighting.

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$GTSV")\` | Full GWAS hit table (rsid, trait, genotype, risk allele, copies, p, OR/beta, gene, pmid, flag) |
| \`$(basename "$CTSV")\` | Full ClinVar pathogenic-class hit table |
EOF

  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
