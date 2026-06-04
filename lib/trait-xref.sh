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
