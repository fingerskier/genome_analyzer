# test_17_vcf_triage_extract.sh — vcf-triage extract: CLI, sections, thresholds, degradation.
XT="$REPO/.claude/skills/vcf-triage/extract-triage.sh"
XDIR="$TMP/triage"; mkdir -p "$XDIR"

# --- CLI ---
out="$(bash "$XT" 2>&1)"; rc=$?
assert_eq "1" "$rc" "triage: no args exits 1"
assert_contains "$out" "Usage" "triage: prints usage on no args"

# --- Missing everything: INPUTS notes each absent file, exit 0 ---
out="$(bash "$XT" "$XDIR" nothere 2>&1)"; rc=$?
assert_eq "0" "$rc" "triage: missing inputs still exits 0"
assert_contains "$out" "== INPUTS ==" "triage: INPUTS section present"
assert_contains "$out" "clinvar: MISSING" "triage: notes missing clinvar table"
assert_contains "$out" "gwas: MISSING" "triage: notes missing gwas table"
assert_contains "$out" "xref-traits.sh" "triage: missing note names the producing script"
assert_contains "$out" "== PHARMGKB ==" "triage: PHARMGKB section present"
assert_contains "$out" "absent" "triage: pharmgkb absent note"

# --- Single-artifact path form derives dir/base ---
printf 'x\n' > "$XDIR/mybase.traits.clinvar.tsv"
out="$(bash "$XT" "$XDIR/mybase.traits.clinvar.tsv" 2>&1)"
assert_contains "$out" "clinvar: found" "triage: artifact-path form finds sibling table"
rm "$XDIR/mybase.traits.clinvar.tsv"

# --- ClinVar triage: build a synthetic table covering every branch ---
CT="$XDIR/synth.traits.clinvar.tsv"
{
  printf '2:100\tA\theterozygous\tPathogenic\tRareCondition\tcriteria_provided,_multiple_submitters,_no_conflicts\tGENE_RARE\t0.0002\tExAC\n'
  printf '3:200\tG\theterozygous\tLikely_pathogenic\tUnknownCondition\tcriteria_provided,_single_submitter\tGENE_UNK\tNA\tNA\n'
  printf '4:300\tC\thomozygous\tPathogenic\tCommonCondition\tno_assertion_criteria_provided\tGENE_COMMON\t0.45\t1000G\n'
  printf '6:600\tT\theterozygous\tConflicting_classifications_of_pathogenicity|risk_factor\tConflictCondition\tcriteria_provided,_conflicting_classifications\tGENE_CONFL\t0.004\tExAC\n'
  printf '1:50\tT\theterozygous\tdrug_response\tdrugX response - Toxicity\treviewed_by_expert_panel\tGENE_DRUG\t0.09\tExAC\n'
  printf '1:60\tT\theterozygous\tdrug_response\tdrugY response\tno_assertion_criteria_provided\tGENE_WEAKDRUG\t0.2\t1000G\n'
  printf '5:400\tA\thomozygous\trisk_factor\tSomeSusceptibility\tno_assertion_criteria_provided\tGENE_RISK\t0.5\t1000G\n'
} > "$CT"
out="$(bash "$XT" "$XDIR" synth 2>&1)"; rc=$?
assert_eq "0" "$rc" "clinvar: exits 0"
prio="$(printf '%s\n' "$out" | sed -n '/== CLINVAR PRIORITY ==/,/== CLINVAR PHARMA ==/p')"
assert_contains "$prio" "GENE_RARE" "clinvar: rare pathogenic in PRIORITY"
assert_contains "$prio" "GENE_UNK" "clinvar: unknown-AF pathogenic in PRIORITY"
case "$prio" in *GENE_COMMON*) assert_eq "no" "yes" "clinvar: common pathogenic NOT in PRIORITY";; *) assert_eq "ok" "ok" "clinvar: common pathogenic NOT in PRIORITY";; esac
case "$prio" in *GENE_CONFL*) assert_eq "no" "yes" "clinvar: conflicting NOT in PRIORITY";; *) assert_eq "ok" "ok" "clinvar: conflicting NOT in PRIORITY";; esac
# weight 3 (GENE_RARE) sorts above weight 2 (GENE_UNK)
first="$(printf '%s\n' "$prio" | grep -n 'GENE_RARE\|GENE_UNK' | head -1)"
assert_contains "$first" "GENE_RARE" "clinvar: higher review weight sorts first in PRIORITY"
pharma="$(printf '%s\n' "$out" | sed -n '/== CLINVAR PHARMA ==/,/== CLINVAR COMMON ==/p')"
assert_contains "$pharma" "GENE_DRUG" "clinvar: expert-panel drug_response in PHARMA"
case "$pharma" in *GENE_WEAKDRUG*) assert_eq "no" "yes" "clinvar: weak drug_response NOT in PHARMA";; *) assert_eq "ok" "ok" "clinvar: weak drug_response NOT in PHARMA";; esac
common="$(printf '%s\n' "$out" | sed -n '/== CLINVAR COMMON ==/,/== GWAS HEADLINE ==/p')"
assert_contains "$common" "GENE_COMMON" "clinvar: common pathogenic in COMMON defuse list"
assert_contains "$common" "GENE_CONFL" "clinvar: conflicting row in COMMON defuse list"
assert_contains "$common" "total 7" "clinvar: COMMON carries total row count"

# Missing clinvar table -> sections still emitted with placeholder
out="$(bash "$XT" "$XDIR" nothere 2>&1)"
assert_contains "$out" "== CLINVAR PRIORITY ==" "clinvar degrade: PRIORITY header present without table"
assert_contains "$out" "(no clinvar table)" "clinvar degrade: placeholder note"
