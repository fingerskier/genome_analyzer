# vcf-triage Skill Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A project Claude skill (`.claude/skills/vcf-triage/`) that turns a completed pipeline run into a plain-language triage report, with a deterministic `extract-triage.sh` codifying all filter thresholds.

**Architecture:** A bash/awk extract script digests the (potentially 20 MB+) hit TSVs into a < 200-line text extract with fixed sections; `SKILL.md` instructs Claude to run it, read the small `*.SUMMARY.md` files directly, and write the report following an embedded template. Missing inputs degrade to notes, never failures.

**Tech Stack:** bash 3.2-compatible shell (macOS `/bin/bash`), BWK awk (no GNU extensions), existing `test/lib.sh` assertion harness.

**Spec:** `docs/superpowers/specs/2026-08-14-vcf-triage-skill-design.md`

## Global Constraints

- bash 3.2 compatible (macOS ships 3.2): no `declare -A`, no `${var,,}`, no mapfile.
- BWK awk only: no `asort`, no `gensub`, no `length(array)`.
- Extract script exit code: always 0 except usage error (exit 1). Missing input files are reported in the extract, never fatal.
- Input TSVs are headerless and tab-separated.
  - ClinVar columns: 1 chrom:pos, 2 allele, 3 zygosity, 4 significance, 5 condition (`|`-separated), 6 review_status, 7 gene, 8 pop_af, 9 af_source.
  - GWAS columns: 1 rsid, 2 trait, 3 genotype, 4 risk_allele, 5 copies, 6 p, 7 or_beta, 8 gene, 9 pmid, 10 flag (`ok`/`ambiguous`), 11 raf.
- Thresholds (single source of truth, from the spec): frequency buckets rare < 0.01, low-frequency < 0.05, common ≥ 0.05, unknown = `NA`/empty; review weights practice_guideline/expert_panel 4 > multiple_submitters 3 > other criteria_provided 2 > else 1; GWAS notable = flag `ok` AND numeric or_beta ≥ 2.0 (no protective side — the column mixes ORs and betas).
- Extract section order (exact header strings): `== INPUTS ==`, `== CLINVAR PRIORITY ==`, `== CLINVAR PHARMA ==`, `== CLINVAR COMMON ==`, `== GWAS HEADLINE ==`, `== GWAS NOTABLE ==`, `== PHARMGKB ==`.
- No personal data in any committed file. Tests use synthetic fixtures only.
- Tests are **sourced** by `test/run-tests.sh` (not executed): no `#!` needed, use `$REPO`, `$TMP`, `assert_eq`, `assert_contains` from `test/lib.sh`. Run the suite with `./test/run-tests.sh`.

---

### Task 1: extract-triage.sh scaffold — CLI, INPUTS, PHARMGKB sections

**Files:**
- Create: `.claude/skills/vcf-triage/extract-triage.sh`
- Test: `test/test_17_vcf_triage_extract.sh`

**Interfaces:**
- Produces: `extract-triage.sh <dir> <base>` or `extract-triage.sh <artifact-path>`; prints extract to stdout; defines stub functions `emit_clinvar_sections` and `emit_gwas_sections` that Tasks 2–3 fill in; shell variables `CTSV`, `GTSV`, `PTSV`, `TAB`, `tmpd` available to those functions.

- [ ] **Step 1: Write the failing tests**

Create `test/test_17_vcf_triage_extract.sh`:

```bash
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `./test/run-tests.sh 2>&1 | grep -A 20 test_17`
Expected: FAIL lines (script does not exist; bash reports no such file, rc mismatch).

- [ ] **Step 3: Write the scaffold implementation**

Create `.claude/skills/vcf-triage/extract-triage.sh`:

```bash
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
```

Then: `chmod +x .claude/skills/vcf-triage/extract-triage.sh`

- [ ] **Step 4: Run tests to verify they pass**

Run: `./test/run-tests.sh 2>&1 | grep -A 20 test_17`
Expected: all `triage:` assertions `ok`; suite still `0 failed` overall.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/vcf-triage/extract-triage.sh test/test_17_vcf_triage_extract.sh
git commit -m "feat: vcf-triage extract scaffold — CLI, INPUTS/PHARMGKB sections, graceful degradation"
```

---

### Task 2: ClinVar sections — PRIORITY / PHARMA / COMMON

**Files:**
- Modify: `.claude/skills/vcf-triage/extract-triage.sh` (replace the `emit_clinvar_sections` stub)
- Modify: `test/test_17_vcf_triage_extract.sh` (append)

**Interfaces:**
- Consumes: `CTSV`, `TAB`, `tmpd` from Task 1's scaffold.
- Produces: extract sections `== CLINVAR PRIORITY ==` (rare/unknown pathogenic-class, review-weight then AF order), `== CLINVAR PHARMA ==` (expert-panel drug_response), `== CLINVAR COMMON ==` (counts line + one-liner defuse list).

- [ ] **Step 1: Write the failing tests**

Append to `test/test_17_vcf_triage_extract.sh`:

```bash
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
```

(Note: `== GWAS HEADLINE ==` does not exist until Task 3; until then the COMMON
`sed` range runs to end-of-output, which still contains the COMMON lines — the
assertions pass either way.)

- [ ] **Step 2: Run tests to verify they fail**

Run: `./test/run-tests.sh 2>&1 | grep -A 40 test_17`
Expected: new `clinvar:` assertions FAIL (sections absent — stub emits nothing).

- [ ] **Step 3: Implement — replace the `emit_clinvar_sections` stub**

Replace `emit_clinvar_sections() { :; }` in `extract-triage.sh` with:

```bash
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
      if (af == "" || af == "NA") return "unknown"
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
      cond=$5; if (length(cond) > 100) cond = substr(cond, 1, 100) "..."
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
```

Also extend the degradation assertions from Task 1 by appending to the test file:

```bash
# Missing clinvar table -> sections still emitted with placeholder
out="$(bash "$XT" "$XDIR" nothere 2>&1)"
assert_contains "$out" "== CLINVAR PRIORITY ==" "clinvar degrade: PRIORITY header present without table"
assert_contains "$out" "(no clinvar table)" "clinvar degrade: placeholder note"
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `./test/run-tests.sh 2>&1 | grep -c FAIL`
Expected: `0`

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/vcf-triage/extract-triage.sh test/test_17_vcf_triage_extract.sh
git commit -m "feat: vcf-triage extract — ClinVar PRIORITY/PHARMA/COMMON triage sections"
```

---

### Task 3: GWAS sections — HEADLINE counts and NOTABLE hits

**Files:**
- Modify: `.claude/skills/vcf-triage/extract-triage.sh` (replace the `emit_gwas_sections` stub)
- Modify: `test/test_17_vcf_triage_extract.sh` (append)

**Interfaces:**
- Consumes: `GTSV`, `TAB` from Task 1's scaffold.
- Produces: `== GWAS HEADLINE ==` (total / homozygous / ambiguous counts) and `== GWAS NOTABLE ==` (flag ok, or_beta ≥ 2.0, rsid-deduplicated keeping strongest p, effect-sorted, capped at 40 with a "more" note). Completes the full section order asserted below.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_17_vcf_triage_extract.sh`:

```bash
# --- GWAS triage ---
GT="$XDIR/synth.traits.gwas.tsv"
{
  printf 'rs1\tDiseaseOne\tA/G\tA\t1\t1E-20\t2.4\tGENE1\t111\tok\t0.3\n'
  printf 'rs2\tMolecularTrait\tC/T\tT\t2\t1E-300\t0.08\tGENE2\t222\tok\t0.5\n'
  printf 'rs3\tDiseaseAmb\tA/T\tA\tNA\t1E-30\t3.0\tGENE3\t333\tambiguous\tNA\n'
  printf 'rs4\tDupWeak\tG/T\tT\t1\t1E-10\t2.2\tGENE4\t444\tok\t0.2\n'
  printf 'rs4\tDupStrong\tG/T\tT\t1\t1E-50\t2.6\tGENE4\t445\tok\t0.2\n'
} > "$GT"
out="$(bash "$XT" "$XDIR" synth 2>&1)"
head_s="$(printf '%s\n' "$out" | sed -n '/== GWAS HEADLINE ==/,/== GWAS NOTABLE ==/p')"
assert_contains "$head_s" "total 5" "gwas: headline total count"
assert_contains "$head_s" "1 homozygous" "gwas: headline homozygous count"
assert_contains "$head_s" "1 strand-ambiguous" "gwas: headline ambiguous count"
notable="$(printf '%s\n' "$out" | sed -n '/== GWAS NOTABLE ==/,/== PHARMGKB ==/p')"
assert_contains "$notable" "DiseaseOne" "gwas: OR 2.4 hit is notable"
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
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `./test/run-tests.sh 2>&1 | grep -A 60 test_17`
Expected: `gwas:` assertions FAIL (stub emits nothing).

- [ ] **Step 3: Implement — replace the `emit_gwas_sections` stub**

Replace `emit_gwas_sections() { :; }` in `extract-triage.sh` with:

```bash
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
  echo "== GWAS NOTABLE ==   (flag ok, OR >= 2.0, dedup by rsid keeping strongest p, cap 40)"
  # p-values like 1E-619 underflow doubles; compare exponent + log10(mantissa).
  nb="$(awk -F"$TAB" '
    function pstrength(p,  a, n, v) {
      n = split(toupper(p), a, "E")
      if (n > 1) return a[2]+0 + (a[1]+0 > 0 ? log(a[1]+0)/log(10) : 0)
      v = p+0
      return (v > 0 ? log(v)/log(10) : -99999)
    }
    $10 == "ok" && $7 ~ /^[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/ && $7+0 >= 2.0 {
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
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `./test/run-tests.sh 2>&1 | tail -3`
Expected: `... passed, 0 failed` (suite grows from 104 by all new test_17 assertions).

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/vcf-triage/extract-triage.sh test/test_17_vcf_triage_extract.sh
git commit -m "feat: vcf-triage extract — GWAS headline counts + notable hits (OR>=2, dedup, cap 40)"
```

---

### Task 4: SKILL.md, docs, real-data smoke test

**Files:**
- Create: `.claude/skills/vcf-triage/SKILL.md`
- Modify: `PLAN.md` (tick the Phase 3 vcf-triage line)
- Modify: `README.md` (new section after "Trait cross-reference")

**Interfaces:**
- Consumes: `extract-triage.sh <dir> <base>` contract and section names from Tasks 1–3.
- Produces: the user-facing skill; report contract `<dir>/<base>.report.md`.

- [ ] **Step 1: Write SKILL.md**

Create `.claude/skills/vcf-triage/SKILL.md` with exactly this content:

````markdown
---
name: vcf-triage
description: Use when asked to triage, interpret, or write a plain-language report over genome-analyzer run outputs (analyze.sh / xref-traits.sh results — *.SUMMARY.md and *.traits.*.tsv files). Produces a full report pairing technical detail with layman interpretation, using codified triage thresholds. Never diagnoses.
---

# vcf-triage — plain-language triage report over pipeline outputs

## Workflow

1. **Locate the run.** If the user names a basename, directory, or file, use it.
   Otherwise take the newest `*.traits.SUMMARY.md` under `work/`, then `out/`.
   The run basename is the filename up to `.traits.`. Tell the user which run
   you selected.
2. **Run the extract** (never read `*.traits.gwas.tsv` directly — it can exceed
   20 MB):

   ```bash
   .claude/skills/vcf-triage/extract-triage.sh <dir> <base>
   ```

   Its sections are the complete, already-thresholded evidence base:
   INPUTS, CLINVAR PRIORITY / PHARMA / COMMON, GWAS HEADLINE / NOTABLE,
   PHARMGKB. Trust its selection; do not re-filter the TSVs yourself.
3. **Read the small run summaries directly** — the per-type `*.SUMMARY.md`
   files listed under INPUTS (SNP/indel QC, SV, CNV).
4. **Write the report** to `<dir>/<base>.report.md` following the template
   below. These directories are git-ignored; never commit a report.
5. **Deliver** the report file to the user (render it) and restate the headline
   findings in chat in plain language.

## Hard rules

- **Never diagnose.** Every section keeps the research-grade framing; the
  report opens with the warning box. Nothing here changes medical decisions
  without clinical confirmation.
- Anything potentially actionable points to **clinical confirmation and
  genetic counseling** — that is always the next step, never your
  interpretation alone.
- **Frequency defuses.** Explicitly explain that scary-named common variants
  (CLINVAR COMMON list) are population-normal database noise. Rarity is a
  reason to look closer, not a verdict.
- **Heterozygous + recessive condition = carrier status**, and say so: one
  working copy means the person is healthy; relevance is family planning.
- Weight ClinVar claims by review status: expert panel > multiple submitters >
  single submitter > no assertion criteria. Say when support is weak.
- Frame the GWAS total (often 100k+) as normal for every genome; only the
  NOTABLE rows deserve prose, and only disease-relevant ones at that.
- No pronoun or relationship assumptions about the sample's person; use the
  name/pronouns the user used, otherwise "the sample".
- Missing inputs (per INPUTS/PHARMGKB sections) collapse to one line in the
  report: what is missing and which command produces it.

## Report template

Pair every section: **The science** (metrics, thresholds, review tiers,
p-values) then **In plain terms** (what it means for a person, no jargon).

1. **Header** — sample, run date, databases used, then a blockquote warning:
   research-grade exploration, not a diagnosis or medical advice.
2. **Data quality** — records/SNP/indel counts, Ti/Tv vs expected ~2.0–2.1,
   PASS fraction. Plain terms: is the data trustworthy; ~4–5M variants is
   normal for every human.
3. **Structural variants** — counts by type, sizes, gene overlap. Plain terms:
   thousands of small SVs are ordinary; large events near immune (HLA) regions
   are expected diversity.
4. **Copy-number variants** — losses/gains, burden in Mb, largest events.
   Plain terms: flag that repeat-rich regions (centromeres, acrocentric arms)
   are where callers are least reliable.
5. **Pharmacogenomics** (CLINVAR PHARMA + PHARMGKB when present) — table of
   gene / genotype / drug(s) / frequency. The most practically useful section:
   these change how the body handles specific drugs, matter only if the drug
   is ever proposed, and then warrant confirmatory clinical testing.
6. **ClinVar findings triaged by rarity** — CLINVAR PRIORITY rows first
   (carrier framing for recessives, review-status caveats), then the COMMON
   defuse list summarized in a sentence or two.
7. **GWAS** — HEADLINE counts with the why-a-huge-number-is-normal framing;
   NOTABLE rows worth individual plain-language callouts (odds ratio ≥ 2
   disease associations), each with lifestyle/screening context where honest.
8. **Caveats** — single-variant sensitivity of consumer WGS, strand-ambiguous
   SNPs unscored, repeat-region CNV reliability, database label quality, no
   polygenic scores or star-allele calling.
9. **Bottom line** — action-tier table: reassurance / mention-if-drug-proposed /
   family-planning-only / routine-screening / no action. Close plainly: is
   anything urgent (usually no), does anything suggest a present health
   problem (usually no).
````

- [ ] **Step 2: Update PLAN.md and README.md**

In `PLAN.md`, change:

```
- [ ] Author a local skill `vcf-triage` that codifies our filter thresholds + report format
```

to:

```
- [x] Author a local skill `vcf-triage` that codifies our filter thresholds + report format
      (`.claude/skills/vcf-triage/` — extract-triage.sh + SKILL.md)
```

In `README.md`, insert after the "Trait cross-reference" section, before "## Testing":

```markdown
## Plain-language report (vcf-triage skill)

With [Claude Code](https://claude.com/claude-code), ask for a triage report
("triage the latest run", "interpret these results") and the project-local
`vcf-triage` skill turns a run's outputs into `<base>.report.md`: technical
detail and layman interpretation side by side, ClinVar hits triaged by
population frequency and review status, an expert-panel pharmacogenomics
table, and GWAS hits filtered to unambiguous odds-ratio ≥ 2 associations.
The thresholds live in `.claude/skills/vcf-triage/extract-triage.sh`, so the
selection is deterministic even though the prose is generated. Reports land
next to the run outputs (git-ignored). Same disclaimer as everything else
here: research-grade exploration, not a diagnosis.
```

- [ ] **Step 3: Full suite + shellcheck**

Run: `./test/run-tests.sh 2>&1 | tail -3` — expected `0 failed`.
Run: `shellcheck .claude/skills/vcf-triage/extract-triage.sh` (if installed) — fix any errors (info/style advisories at your judgment).

- [ ] **Step 4: Real-data smoke test**

Run the extract against the existing real run (output stays local, never committed):

```bash
.claude/skills/vcf-triage/extract-triage.sh work "$(basename "$(ls work/*.traits.SUMMARY.md | head -1)" .traits.SUMMARY.md)"
```

Expected: all seven sections; CLINVAR PRIORITY contains a small number of rows (order 2–5); CLINVAR PHARMA lists the expert-panel drug-response hits (order ~13); GWAS NOTABLE ≤ 40 rows, no molecular-trait beta flood; runtime under ~30 s on the 23 MB GWAS table.

- [ ] **Step 5: Commit**

```bash
git add .claude/skills/vcf-triage/SKILL.md PLAN.md README.md
git commit -m "feat: vcf-triage skill — plain-language triage report over run outputs (Phase 3)"
```

---

## Self-review notes

- Spec coverage: SKILL.md workflow/hard-rules/template (Task 4), extract contract + all seven sections and thresholds (Tasks 1–3), all six spec test bullets present (entry/degradation Task 1+2, ClinVar routing Task 2, GWAS filter/dedup/cap Task 3, section order Task 3), privacy (constraints + Task 4 README wording).
- Spec deviation, already amended in the spec: GWAS notable is ≥ 2.0 only (no ≤ 0.5 protective side — beta/OR ambiguity).
- Type consistency: section header strings, `extract-triage.sh <dir> <base>` invocation, and `(no clinvar table)` / `(no gwas table)` placeholders are identical across tasks and tests.
