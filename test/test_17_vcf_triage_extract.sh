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

# --- GWAS triage ---
GT="$XDIR/synth.traits.gwas.tsv"
{
  printf 'rs1\tDiseaseOne\tA/G\tA\t1\t1E-20\t2.4\tGENE1\t111\tok\t0.3\n'
  printf 'rs2\tMolecularTrait\tC/T\tT\t2\t1E-300\t0.08\tGENE2\t222\tok\t0.5\n'
  printf 'rs3\tDiseaseAmb\tA/T\tA\tNA\t1E-30\t3.0\tGENE3\t333\tambiguous\tNA\n'
  printf 'rs4\tDupWeak\tG/T\tT\t1\t1E-10\t2.2\tGENE4\t444\tok\t0.2\n'
  printf 'rs4\tDupStrong\tG/T\tT\t1\t1E-50\t2.6\tGENE4\t445\tok\t0.2\n'
  printf 'rs5\tDiseaseFlip\tG/T\tT\t2\t1E-40\t5.0\tGENE5\t555\tstrand_flipped\t0.15\n'
} > "$GT"
out="$(bash "$XT" "$XDIR" synth 2>&1)"
head_s="$(printf '%s\n' "$out" | sed -n '/== GWAS HEADLINE ==/,/== GWAS NOTABLE ==/p')"
assert_contains "$head_s" "total 6" "gwas: headline total count"
assert_contains "$head_s" "2 homozygous" "gwas: headline homozygous count"
assert_contains "$head_s" "1 strand-ambiguous" "gwas: headline ambiguous count"
notable="$(printf '%s\n' "$out" | sed -n '/== GWAS NOTABLE ==/,/== PHARMGKB ==/p')"
assert_contains "$notable" "DiseaseOne" "gwas: OR 2.4 hit is notable"
assert_contains "$notable" "DiseaseFlip" "gwas: strand_flipped scored hit is notable"
case "$notable" in *MolecularTrait*) assert_eq "no" "yes" "gwas: beta 0.08 NOT notable";; *) assert_eq "ok" "ok" "gwas: beta 0.08 NOT notable";; esac
case "$notable" in *DiseaseAmb*) assert_eq "no" "yes" "gwas: ambiguous NOT notable";; *) assert_eq "ok" "ok" "gwas: ambiguous NOT notable";; esac
assert_contains "$notable" "DupStrong" "gwas: dedup keeps strongest p"
case "$notable" in *DupWeak*) assert_eq "no" "yes" "gwas: dedup drops weaker p";; *) assert_eq "ok" "ok" "gwas: dedup drops weaker p";; esac

# --- Cap at 40 with a "more" note ---
CAPT="$XDIR/cap.traits.gwas.tsv"
: > "$CAPT"
i=1
while [ "$i" -le 45 ]; do
  printf 'rsc%d\tCapTrait%d\tA/G\tG\t1\t1E-9\t2.5\tGENEC\t9%d\tok\t0.1\n' "$i" "$i" "$i" >> "$CAPT"
  i=$((i+1))
done
out="$(bash "$XT" "$XDIR" cap 2>&1)"
shown="$(printf '%s\n' "$out" | grep -c 'CapTrait')"
assert_eq "40" "$shown" "gwas: notable capped at 40 rows"
assert_contains "$out" "5 more" "gwas: cap emits a 'more' note"

# --- Full section order with all inputs present ---
order="$(bash "$XT" "$XDIR" synth 2>&1 | grep '^== ' | tr -d '=' | tr -d ' ' | tr '\n' ',')"
assert_eq "INPUTS,CLINVARPRIORITY,CLINVARPHARMA,CLINVARCOMMON,GWASHEADLINE,GWASNOTABLE,PHARMGKB," "$order" "triage: sections in spec order"

# --- Empty placeholder tables (DB skipped by xref-traits.sh) are not clean negatives ---
EDIR="$TMP/triage_empty"; mkdir -p "$EDIR"
: > "$EDIR/skip.traits.clinvar.tsv"
: > "$EDIR/skip.traits.gwas.tsv"
out="$(bash "$XT" "$EDIR" skip 2>&1)"; rc=$?
assert_eq "0" "$rc" "empty: exits 0"
assert_contains "$out" "clinvar: present but EMPTY" "empty: clinvar empty note in INPUTS"
assert_contains "$out" "gwas: present but EMPTY" "empty: gwas empty note in INPUTS"
assert_contains "$out" "database was likely skipped" "empty: skip explanation present"
case "$out" in *"(none — no rare"*) assert_eq "no" "yes" "empty: PRIORITY does not claim a clean negative";; *) assert_eq "ok" "ok" "empty: PRIORITY does not claim a clean negative";; esac
assert_contains "$out" "(clinvar table empty — analysis may not have run; not a clean negative)" "empty: clinvar sections carry empty-table note"
assert_contains "$out" "(gwas table empty — analysis may not have run; not a clean negative)" "empty: gwas sections carry empty-table note"

# --- Summaries listing is scoped to the selected run ---
printf '# a\n' > "$XDIR/alpha.snp-indel.genome.SUMMARY.md"
printf '# a\n' > "$XDIR/alpha.sv.SUMMARY.md"
printf '# b\n' > "$XDIR/beta.snp-indel.genome.SUMMARY.md"
out="$(bash "$XT" "$XDIR" alpha.snp-indel.genome.pass 2>&1)"
assert_contains "$out" "alpha.sv.SUMMARY.md" "scope: same-run sv summary listed"
case "$out" in *beta.snp-indel*) assert_eq "no" "yes" "scope: other-run summary NOT listed";; *) assert_eq "ok" "ok" "scope: other-run summary NOT listed";; esac
rm "$XDIR/alpha.snp-indel.genome.SUMMARY.md" "$XDIR/alpha.sv.SUMMARY.md" "$XDIR/beta.snp-indel.genome.SUMMARY.md"
