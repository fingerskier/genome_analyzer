# Population-Frequency Context Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add population allele-frequency context (value + plain-language bucket) to every trait cross-reference hit, harvested from data we already cache — zero new downloads.

**Architecture:** `fetch-gwas.sh` keeps the GWAS Catalog's `RISK ALLELE FREQUENCY` column as an 8th normalized column; `lib/trait-xref.sh` passes it through to the GWAS hit table, extracts `AF_EXAC`/`AF_TGP`/`AF_ESP` from the cached ClinVar VCF for ClinVar hits, and renders both plus frequency buckets in SUMMARY.md. A shared awk fragment (`FREQ_FN`, like the existing `DOSAGE_FN`) is the single source of truth for bucket thresholds.

**Tech Stack:** bash 3.2 (macOS default), awk, bcftools/htslib. Zero-dependency test harness in `test/`.

## Global Constraints

- All scripts must run on **bash 3.2**: no associative arrays, no `mapfile`, no `${var,,}`.
- **No new dependencies, no new downloads, no network at xref time** (privacy: never send carried rsIDs anywhere).
- **Append-only TSV columns**: new columns go at the end; existing column positions must not move.
- Test fixtures are **synthetic only** — fake rsIDs/conditions, never real genomic data.
- Sorts that feed joins use `LC_ALL=C`.
- Spec: `docs/superpowers/specs/2026-08-09-frequency-context-design.md`.
- Run the whole suite with `./test/run-tests.sh` from the repo root; it must end `0 failed`.

---

### Task 1: `normalize_gwas` keeps RISK ALLELE FREQUENCY (col 27)

**Files:**
- Modify: `fetch-gwas.sh` (normalize_gwas, lines 11–27)
- Modify: `test/fixtures/gwas_raw.tsv`
- Create: `test/fixtures/gwas_raw_drift.tsv`
- Test: `test/test_13_gwas_normalize.sh`

**Interfaces:**
- Produces: `data/gwas-catalog.tsv` rows become 8 columns: `rsid  risk  trait  pval  orbeta  gene  pmid  raf`, where `raf` is a plain number in [0,1] or the literal string `NA`. Task 3's engine reads column 8.

- [ ] **Step 1: Update the raw fixture and add a header-drift fixture**

Replace `test/fixtures/gwas_raw.tsv` with (columns are TAB-separated; col 27 header is now `RISK ALLELE FREQUENCY`; two new clean rows exercise `NR` and annotated frequencies):

```
COL1	PUBMEDID	x	x	x	x	x	DISEASE/TRAIT	x	x	x	x	x	x	MAPPED_GENE	x	x	x	x	x	STRONGEST SNP-RISK ALLELE	x	x	x	x	x	RISK ALLELE FREQUENCY	P-VALUE	x	x	OR or BETA
x	111	x	x	x	x	x	Height	x	x	x	x	x	x	GENEH	x	x	x	x	x	rs100-G	x	x	x	x	x	0.23	1e-12	x	x	1.2
x	222	x	x	x	x	x	Multi locus	x	x	x	x	x	x	GENEM	x	x	x	x	x	rs7; rs8-A	x	x	x	x	x	x	2e-12	x	x	1.3
x	333	x	x	x	x	x	Unknown allele	x	x	x	x	x	x	GENEU	x	x	x	x	x	rs300-?	x	x	x	x	x	x	3e-12	x	x	1.1
x	444	x	x	x	x	x	Indel risk	x	x	x	x	x	x	GENEI	x	x	x	x	x	rs999-AT	x	x	x	x	x	x	4e-12	x	x	1.0
x	555	x	x	x	x	x	TraitNR	x	x	x	x	x	x	GENEN	x	x	x	x	x	rs101-A	x	x	x	x	x	NR	2e-12	x	x	1.1
x	666	x	x	x	x	x	TraitEA	x	x	x	x	x	x	GENEE	x	x	x	x	x	rs102-C	x	x	x	x	x	0.23 (EA)	3e-12	x	x	1.0
```

Create `test/fixtures/gwas_raw_drift.tsv` — the OLD header (col 27 is `x`, i.e. the layout drifted) plus one clean row that DOES have a number at col 27; the header check must force raf to NA anyway:

```
COL1	PUBMEDID	x	x	x	x	x	DISEASE/TRAIT	x	x	x	x	x	x	MAPPED_GENE	x	x	x	x	x	STRONGEST SNP-RISK ALLELE	x	x	x	x	x	x	P-VALUE	x	x	OR or BETA
x	111	x	x	x	x	x	Height	x	x	x	x	x	x	GENEH	x	x	x	x	x	rs100-G	x	x	x	x	x	0.23	1e-12	x	x	1.2
```

- [ ] **Step 2: Write the failing tests**

Replace `test/test_13_gwas_normalize.sh` with (needles contain literal TABs):

```bash
# test_13_gwas_normalize.sh — fetch-gwas.sh::normalize_gwas (no network).
source "$REPO/fetch-gwas.sh" --source-only

out="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw.tsv" 2>/dev/null)"
assert_eq "3" "$(printf '%s\n' "$out" | grep -c .)" "normalize: the 3 clean rows survive"
assert_contains "$out" "rs100	G	Height	1e-12	1.2	GENEH	111	0.23" "normalize: raf kept as 8th column"
assert_contains "$out" "rs101	A	TraitNR	2e-12	1.1	GENEN	555	NA" "normalize: NR raf -> NA"
assert_contains "$out" "rs102	C	TraitEA	3e-12	1.0	GENEE	666	NA" "normalize: annotated '0.23 (EA)' raf -> NA"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Multi locus')"    "normalize: multi-SNP row dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Unknown allele')" "normalize: '?' risk allele dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Indel risk')"     "normalize: multi-base risk allele dropped"

# Header drift: col 27 is not RISK ALLELE FREQUENCY -> warn, raf always NA.
drift_out="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw_drift.tsv" 2>/dev/null)"
assert_contains "$drift_out" "rs100	G	Height	1e-12	1.2	GENEH	111	NA" "normalize: header drift -> raf NA, not a wrong column"
drift_err="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw_drift.tsv" 2>&1 >/dev/null)"
assert_contains "$drift_err" "WARNING" "normalize: header drift warns on stderr"
```

- [ ] **Step 3: Run the suite to verify the new asserts fail**

Run: `./test/run-tests.sh 2>&1 | sed -n '/test_13/,/test_14/p'`
Expected: FAILs for the raf column asserts (rows currently end at col 7).

- [ ] **Step 4: Implement in `fetch-gwas.sh`**

Replace the `normalize_gwas` function (and update its doc comment) with:

```bash
# normalize_gwas: GWAS Catalog associations TSV on stdin -> normalized table.
# Columns read (1-based): 2 PUBMEDID, 8 DISEASE/TRAIT, 15 MAPPED_GENE,
# 21 STRONGEST SNP-RISK ALLELE (e.g. "rs1234-A"), 27 RISK ALLELE FREQUENCY,
# 28 P-VALUE, 31 OR or BETA. Col 27 is sanity-checked against the header: on
# drift we warn and emit NA rather than harvest a wrong column.
normalize_gwas() {
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1{
      if($27 != "RISK ALLELE FREQUENCY"){
        badhdr=1
        print "WARNING: header col 27 is not RISK ALLELE FREQUENCY; raf column set to NA" > "/dev/stderr"
      }
      next
    }
    {
      sra=$21
      p=0; for(i=length(sra);i>=1;i--){ if(substr(sra,i,1)=="-"){p=i; break} }
      if(p==0) next
      rsid=substr(sra,1,p-1); risk=substr(sra,p+1)
      if(rsid !~ /^rs[0-9]+$/) next
      if(risk !~ /^[ACGT]$/)   next
      trait=$8; gsub(/[\t\r]/," ",trait)
      raf=$27; gsub(/\r/,"",raf)
      if(badhdr || raf !~ /^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ || raf+0>1) raf="NA"
      print rsid, risk, trait, $28, $31, $15, $2, raf
    }'
}
```

Also update the file-top comment (line 2–4) to mention the 8th `raf` column.

- [ ] **Step 5: Run the suite to verify test_13 passes**

Run: `./test/run-tests.sh 2>&1 | tail -3`
Expected: `0 failed` (test_15/15b GWAS fixtures are separate files, unaffected).

- [ ] **Step 6: Commit**

```bash
git add fetch-gwas.sh test/fixtures/gwas_raw.tsv test/fixtures/gwas_raw_drift.tsv test/test_13_gwas_normalize.sh
git commit -m "feat: normalize_gwas keeps RISK ALLELE FREQUENCY as raf column (NA on drift/non-numeric)"
```

---

### Task 2: `freq_bucket` — shared frequency-bucket function

**Files:**
- Modify: `lib/trait-xref.sh` (add after `risk_dosage`, line 29)
- Create: `test/test_14a_freq_bucket.sh`

**Interfaces:**
- Produces: awk fragment `FREQ_FN` defining `bucket(af)` → `"<bucket>" SUBSEP "<display>"`, and shell wrapper `freq_bucket <af>` → `"<bucket>\t<display>"`. Buckets: `common` (≥0.05), `low-frequency` (≥0.01), `rare` (<0.01), `unknown` (non-numeric or >1). Displays: `common (~40%)`, `low-frequency (~2%)`, `rare (<1%)`, `unknown`. Tasks 3–4 prepend `$FREQ_FN` to summary awk programs.

- [ ] **Step 1: Write the failing test**

Create `test/test_14a_freq_bucket.sh`:

```bash
# test_14a_freq_bucket.sh — population-frequency buckets (lib/trait-xref.sh).
source "$REPO/lib/trait-xref.sh"

# freq_bucket <af> -> "<bucket>\t<display>"
assert_eq "common	common (~40%)"               "$(freq_bucket 0.4)"   "bucket: 40% -> common"
assert_eq "common	common (~5%)"                "$(freq_bucket 0.05)"  "bucket: 5% boundary -> common"
assert_eq "low-frequency	low-frequency (~2%)" "$(freq_bucket 0.02)"  "bucket: 2% -> low-frequency"
assert_eq "low-frequency	low-frequency (~1%)" "$(freq_bucket 0.01)"  "bucket: 1% boundary -> low-frequency"
assert_eq "rare	rare (<1%)"                    "$(freq_bucket 0.001)" "bucket: 0.1% -> rare"
assert_eq "rare	rare (<1%)"                    "$(freq_bucket 5e-3)"  "bucket: scientific notation accepted"
assert_eq "unknown	unknown"                   "$(freq_bucket NA)"    "bucket: NA -> unknown"
assert_eq "unknown	unknown"                   "$(freq_bucket abc)"   "bucket: junk -> unknown"
assert_eq "unknown	unknown"                   "$(freq_bucket "")"    "bucket: empty -> unknown"
assert_eq "unknown	unknown"                   "$(freq_bucket 1.5)"   "bucket: >1 is not a frequency -> unknown"
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh 2>&1 | sed -n '/test_14a/,/test_15/p'`
Expected: FAILs — `freq_bucket: command not found` (empty actual).

- [ ] **Step 3: Implement in `lib/trait-xref.sh`**

Insert after the `risk_dosage()` function (line 29):

```bash
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh 2>&1 | tail -3`
Expected: `0 failed`.

- [ ] **Step 5: Commit**

```bash
git add lib/trait-xref.sh test/test_14a_freq_bucket.sh
git commit -m "feat: freq_bucket shared frequency-bucket function (common/low-frequency/rare/unknown)"
```

---

### Task 3: engine passes GWAS raf through to TSV + summary

**Files:**
- Modify: `lib/trait-xref.sh` (`run_trait_xref`: GWAS join awk ~lines 68–90, `gwas_rows` ~line 138, summary GWAS table + artifacts listing)
- Modify: `test/fixtures/gwas_norm.tsv`
- Test: `test/test_15_xref.sh`

**Interfaces:**
- Consumes: 8-column GWAS table from Task 1 (col 8 = raf, number or `NA`).
- Produces: `*.traits.gwas.tsv` becomes 11 columns: `rsid trait genotype risk_allele copies pval or_beta gene pmid flag raf`. Rows read from a legacy 7-column table get raf `NA` (Task 5 adds the summary note + test for that). Summary GWAS table gains a `Risk-allele freq` column rendering `23%` or `—`.

- [ ] **Step 1: Extend the normalized-GWAS fixture**

Replace `test/fixtures/gwas_norm.tsv` with (adds col 8):

```
rs100	G	Height	1e-12	1.2	GENEH	111	0.23
rs200	T	Coffee consumption	3e-20	1.5	GENEC	222	0.10
rs200	A	Weak association	1e-6	0.9	GENEC	333	0.50
rs300	A	Palindromic trait	1e-9	1.1	GENEP	444	0.30
```

- [ ] **Step 2: Write the failing tests**

Append to `test/test_15_xref.sh` (needles contain literal TABs):

```bash
# Frequency context: raf is appended as GTSV column 11 and rendered as a
# percent in the summary GWAS table.
assert_contains "$gtsv" "GENEH	111	ok	0.23" "gwas tsv: raf appended after flag (col 11)"
assert_contains "$sum" "Risk-allele freq" "gwas summary: frequency column header present"
assert_contains "$sum" "| 23% |"          "gwas summary: raf rendered as rounded percent"
```

- [ ] **Step 3: Run to verify the new asserts fail**

Run: `./test/run-tests.sh 2>&1 | sed -n '/test_15/,/test_15a/p'`
Expected: the three new asserts FAIL; pre-existing asserts still pass.

- [ ] **Step 4: Implement in `lib/trait-xref.sh`**

(a) In the GWAS join awk (inside `run_trait_xref`), change the record-packing
line so raf rides along (legacy 7-col tables fall back to `NA`):

```awk
          rec=$2 SUBSEP $3 SUBSEP $4 SUBSEP $5 SUBSEP $6 SUBSEP $7 SUBSEP (NF>=8?$8:"NA")
```

and the output line to append it after `flag`:

```awk
            print rs, f[2], geno, f[1], dose, f[3], f[4], f[5], f[6], flag, f[7]
```

(b) Replace the `gwas_rows` builder so the summary renders the percent (em
dash for NA/absent):

```bash
  gwas_rows="$(awk -F'\t' '$10!="ambiguous"' "$GTSV" \
    | LC_ALL=C sort -t"$(printf '\t')" -k6,6g \
    | awk -F'\t' 'NR<=50{
        hl=($5==2?"**":"")
        raf=($11=="" || $11=="NA" ? "—" : sprintf("%.0f%%", $11*100))
        printf "| %s | %s | %s%s%s | %s | %s | %s | %s | %s |\n", $2,$8,hl,$3,hl,$4,$5,$6,$7,raf}')"
  : "${gwas_rows:=| _none carried_ |  |  |  |  |  |  |  |}"
```

(c) In the SUMMARY heredoc, widen the GWAS table:

```markdown
| Trait | Gene | Your genotype | Risk allele | Copies | p-value | OR/beta | Risk-allele freq |
|---|---|---|---|---|---|---|---|
```

(d) In the Artifacts table, update the GTSV description to
`(rsid, trait, genotype, risk allele, copies, p, OR/beta, gene, pmid, flag, risk-allele freq)`.

- [ ] **Step 5: Run to verify test_15 passes**

Run: `./test/run-tests.sh 2>&1 | tail -3`
Expected: `0 failed`.

- [ ] **Step 6: Commit**

```bash
git add lib/trait-xref.sh test/fixtures/gwas_norm.tsv test/test_15_xref.sh
git commit -m "feat: GWAS risk-allele frequency passed through to hit table and summary"
```

---

### Task 4: engine extracts ClinVar AF + buckets in summary

**Files:**
- Modify: `lib/trait-xref.sh` (`run_trait_xref`: ClinVar query/join ~lines 104–120, `clinvar_rows` ~line 130, summary ClinVar table, "How to read this" bullets, artifacts listing)
- Modify: `test/fixtures/clinvar.vcf.txt`
- Test: `test/test_15_xref.sh`

**Interfaces:**
- Consumes: `FREQ_FN` from Task 2.
- Produces: `*.traits.clinvar.tsv` becomes 9 columns: `chrom:pos alt zygosity significance condition review gene pop_af af_source`, with `pop_af`/`af_source` = `NA`/`NA` when no cohort reported it (or the VCF lacks the AF tags entirely — Task 5 tests that). `af_source` ∈ `ExAC`,`1000G`,`ESP`,`NA`. Summary ClinVar table gains a `How common?` column using `bucket()` displays.

- [ ] **Step 1: Extend the ClinVar fixture**

Replace `test/fixtures/clinvar.vcf.txt` with (adds the three AF INFO header
lines; the 2:500 Pathogenic record gains `AF_EXAC=0.30`; the 1:200 record
stays AF-less to exercise `unknown`):

```
##fileformat=VCFv4.2
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="significance">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="disease name">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="review status">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="gene">
##INFO=<ID=AF_EXAC,Number=1,Type=Float,Description="allele frequencies from ExAC">
##INFO=<ID=AF_TGP,Number=1,Type=Float,Description="allele frequencies from TGP">
##INFO=<ID=AF_ESP,Number=1,Type=Float,Description="allele frequencies from GO-ESP">
##contig=<ID=1>
##contig=<ID=2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	500	900	G	C	.	.	CLNSIG=Pathogenic;CLNDN=Example_condition;CLNREVSTAT=criteria_provided;GENEINFO=GENEX:999;AF_EXAC=0.30
1	100	901	A	G	.	.	CLNSIG=Benign;CLNDN=Nothing_notable;CLNREVSTAT=criteria_provided;GENEINFO=GENEY:888
1	200	902	C	T	.	.	CLNSIG=Pathogenic;CLNDN=Cond_one|Cond_two|Cond_three|Cond_four;CLNREVSTAT=criteria_provided;GENEINFO=GENEZ:777
```

- [ ] **Step 2: Write the failing tests**

Append to `test/test_15_xref.sh`:

```bash
# ClinVar population frequency: AF_EXAC harvested with source, absent AF ->
# NA/unknown, buckets rendered in the summary with a plain-language bullet.
assert_contains "$ctsv" "ExAC"            "clinvar tsv: AF source recorded for the ExAC-backed hit"
assert_contains "$ctsv" "NA	NA"           "clinvar tsv: AF-less record -> pop_af NA, source NA"
assert_contains "$sum" "How common?"      "clinvar summary: frequency column header present"
assert_contains "$sum" "common (~30%)"    "clinvar summary: AF_EXAC=0.30 -> common bucket with percent"
assert_contains "$sum" "| unknown |"      "clinvar summary: AF-less hit shows unknown"
assert_contains "$sum" "large fraction of the population" "summary: plain-language common-variant caveat present"
```

- [ ] **Step 3: Run to verify the new asserts fail**

Run: `./test/run-tests.sh 2>&1 | sed -n '/test_15/,/test_15a/p'`
Expected: the six new asserts FAIL; pre-existing asserts still pass.

- [ ] **Step 4: Implement in `lib/trait-xref.sh`**

(a) Before the ClinVar `bcftools query`, build the format string with an AF
sanity check (a custom `-c` VCF may lack the tags; `bcftools query` errors on
undefined tags, so degrade to no-AF extraction):

```bash
    log "ClinVar cross-reference (pathogenic-class)"
    local cv_fmt='%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/GENEINFO'
    if bcftools view -h "$CLINVAR_VCF" 2>/dev/null | grep -q 'ID=AF_EXAC'; then
      cv_fmt="$cv_fmt"'\t%INFO/AF_EXAC\t%INFO/AF_TGP\t%INFO/AF_ESP'
    fi
```

and use it: `-f "$cv_fmt"'\n'` (replacing the literal format string).

(b) In the join awk's `FNR==NR` block, harvest the first available cohort AF
(largest cohort first) and record its source:

```awk
          FNR==NR{
            k=$1":"$2":"$3":"$4
            g=$8; gsub(/:[0-9]+/,"",g); gsub(/\|/,", ",g)
            af="NA"; src="NA"
            if(NF>=9 && $9!="." && $9!=""){ af=$9; src="ExAC" }
            else if(NF>=10 && $10!="." && $10!=""){ af=$10; src="1000G" }
            else if(NF>=11 && $11!="." && $11!=""){ af=$11; src="ESP" }
            sig[k]=$5; dn[k]=$6; rev[k]=$7; gene[k]=g; paf[k]=af; psrc[k]=src; next
          }
```

and append both to the output line:

```awk
              print $1":"$2, $5, ($6=="hom"?"homozygous":"heterozygous"), sig[k], cond, rev[k], gene[k], paf[k], psrc[k]
```

(c) Replace the `clinvar_rows` builder — prepend `$FREQ_FN` and render the
bucket display as a new last column:

```bash
  clinvar_rows="$(awk -F'\t' "$FREQ_FN"'NR<=50{
      sig=$4; gene=$7; gsub(/\|/,"; ",sig); gsub(/\|/,"; ",gene)
      n=split($5, cc, /\|/); cond=""
      for(i=1;i<=n && i<=3;i++) cond=cond (i>1?"; ":"") cc[i]
      if(n>3) cond=cond " (+" (n-3) " more)"
      split(bucket($8), bb, SUBSEP)
      printf "| %s | %s | %s | %s | %s | %s |\n", cond, gene, $3, sig, $6, bb[2]
    }' "$CTSV")"
  : "${clinvar_rows:=| _none carried_ |  |  |  |  |  |}"
```

(d) In the SUMMARY heredoc, widen the ClinVar table:

```markdown
| Condition | Gene | Your call | Significance | Review status | How common? |
|---|---|---|---|---|---|
```

(e) Add two bullets to "How to read this" (after the review-status bullet):

```markdown
- **How common?** buckets population frequency: *common* (≥5%), *low-frequency*
  (1–5%), *rare* (<1%). A variant carried by a large fraction of the population is
  almost never seriously harmful on its own, despite alarming condition names —
  rarity is a reason to look closer, not a verdict. ClinVar frequencies come from
  the ExAC / 1000 Genomes / ESP cohorts; GWAS risk-allele frequencies from the
  catalog's reporting study.
- An \`unknown\` frequency means those cohorts did not report the variant — often
  a hint of rarity, but also routine for indels and recently catalogued variants.
```

(f) In the Artifacts table, update the CTSV description to
`(chrom:pos, allele, zygosity, significance, condition, review status, gene, pop AF, AF source)`.

- [ ] **Step 5: Run to verify test_15 passes**

Run: `./test/run-tests.sh 2>&1 | tail -3`
Expected: `0 failed`.

- [ ] **Step 6: Commit**

```bash
git add lib/trait-xref.sh test/fixtures/clinvar.vcf.txt test/test_15_xref.sh
git commit -m "feat: ClinVar hits gain population AF (ExAC/1000G/ESP) + frequency buckets in summary"
```

---

### Task 5: degradation — legacy 7-col GWAS cache and AF-less ClinVar VCF

**Files:**
- Modify: `lib/trait-xref.sh` (`run_trait_xref`: gwas_note, ~line 95)
- Create: `test/fixtures/gwas_norm_legacy.tsv`
- Create: `test/fixtures/clinvar_noaf.vcf.txt`
- Test: create `test/test_15b_xref_freq_degrade.sh`

**Interfaces:**
- Consumes: Task 3's GTSV (col 11 raf) and Task 4's CTSV (cols 8–9 `pop_af`,`af_source`) and header-check degradation.
- Produces: a legacy (pre-raf) GWAS table yields raf `NA` on every row plus a summary refresh note containing `fetch-gwas.sh -f`; an AF-less ClinVar VCF yields `NA	NA` and still exits 0.

- [ ] **Step 1: Create the degradation fixtures**

`test/fixtures/gwas_norm_legacy.tsv` (the pre-feature 7-column shape):

```
rs100	G	Height	1e-12	1.2	GENEH	111
rs200	T	Coffee consumption	3e-20	1.5	GENEC	222
```

`test/fixtures/clinvar_noaf.vcf.txt` (no AF_* header tags at all):

```
##fileformat=VCFv4.2
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="significance">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="disease name">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="review status">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="gene">
##contig=<ID=1>
##contig=<ID=2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	500	910	G	C	.	.	CLNSIG=Pathogenic;CLNDN=NoAF_condition;CLNREVSTAT=criteria_provided;GENEINFO=GENEQ:111
```

- [ ] **Step 2: Write the failing test**

Create `test/test_15b_xref_freq_degrade.sh`:

```bash
# test_15b_xref_freq_degrade.sh — frequency degradation: legacy 7-col GWAS
# cache (pre-raf) and a ClinVar VCF without AF_* tags must not break the run.
( source "$REPO/common.sh"; source "$REPO/lib/trait-xref.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/xref_sample.vcf.txt")"
  CLINVAR_VCF="$(make_fixture "$REPO/test/fixtures/clinvar_noaf.vcf.txt")"
  GWAS_TBL="$REPO/test/fixtures/gwas_norm_legacy.tsv"
  OUTDIR="$TMP/xref_degrade_out"; mkdir -p "$OUTDIR"; BASE="xrefd"
  run_trait_xref )
rc=$?
assert_eq "0" "$rc" "freq degrade: run exits 0"

dsum="$(cat "$TMP/xref_degrade_out/xrefd.traits.SUMMARY.md")"
dgtsv="$(cat "$TMP/xref_degrade_out/xrefd.traits.gwas.tsv")"
dctsv="$(cat "$TMP/xref_degrade_out/xrefd.traits.clinvar.tsv")"

assert_eq "0" "$(printf '%s\n' "$dgtsv" | awk -F'\t' '$11!="NA"' | grep -c .)" "freq degrade: legacy table -> every raf is NA"
assert_contains "$dsum" "fetch-gwas.sh -f" "freq degrade: summary tells the user how to refresh the catalog"
assert_contains "$dctsv" "NoAF condition" "freq degrade: AF-less ClinVar hit still reported"
assert_contains "$dctsv" "NA	NA" "freq degrade: AF-less ClinVar VCF -> pop_af/source NA"
```

- [ ] **Step 3: Run to verify it fails**

Run: `./test/run-tests.sh 2>&1 | sed -n '/test_15b/,/test_16/p'`
Expected: the `fetch-gwas.sh -f` assert FAILs (note not implemented yet); the NA asserts already pass via Tasks 3–4 fallbacks.

- [ ] **Step 4: Implement the refresh note in `lib/trait-xref.sh`**

In `run_trait_xref`, inside the GWAS `else` branch, detect the table width
once (before the join awk):

```bash
    log "GWAS Catalog cross-reference (p < 5e-8)"
    local gwas_cols; gwas_cols="$(head -n 1 "$GWAS_TBL" | awk -F'\t' '{print NF}')"
```

and after `gwas_note=` is assigned, append the refresh hint for legacy tables:

```bash
    if [[ "${gwas_cols:-8}" -lt 8 ]]; then
      gwas_note="$gwas_note _Frequencies unavailable in this cached catalog — refresh with \`./fetch-gwas.sh -f\`._"
    fi
```

- [ ] **Step 5: Run the full suite**

Run: `./test/run-tests.sh 2>&1 | tail -3`
Expected: `0 failed`.

- [ ] **Step 6: Commit**

```bash
git add lib/trait-xref.sh test/fixtures/gwas_norm_legacy.tsv test/fixtures/clinvar_noaf.vcf.txt test/test_15b_xref_freq_degrade.sh
git commit -m "feat: graceful frequency degradation (legacy GWAS cache note, AF-less ClinVar VCF)"
```

---

### Task 6: docs — README, PLAN.md, full-suite verification

**Files:**
- Modify: `README.md` (trait cross-reference section, ~lines 89–124)
- Modify: `PLAN.md` (Phase 2 gnomAD line, line 27)

**Interfaces:**
- Consumes: everything above; documentation only.

- [ ] **Step 1: Update README.md**

In the "Trait cross-reference" section, replace the outputs paragraph
(currently starting "Outputs: `*.traits.SUMMARY.md` …") with:

```markdown
Outputs: `*.traits.SUMMARY.md` (readable report — ClinVar pathogenic hits, then
GWAS trait hits you carry ≥1 risk allele for, homozygous hits highlighted),
`*.traits.gwas.tsv`, and `*.traits.clinvar.tsv` (full tables). Every hit now
carries **population-frequency context**: ClinVar hits get an allele frequency
harvested from the ExAC / 1000 Genomes / ESP fields already inside the ClinVar
VCF, GWAS hits get the catalog's reported risk-allele frequency, and the
summary buckets them in plain language (*common* ≥5%, *low-frequency* 1–5%,
*rare* <1%, *unknown*). Frequency is triage context, not a verdict — a
"pathogenic"-labelled variant that a third of the population carries is almost
always low-impact, while *rare* or *unknown* is a reason to look closer. If you
fetched the GWAS catalog before this feature, refresh it once with
`./fetch-gwas.sh -f` (~68 MB); until then the summary notes that frequencies
are unavailable. If a database file is missing, the run still completes and the
summary tells you which fetch script to run. **Strand-ambiguous SNPs** (A/T,
C/G) can't be oriented from the VCF alone, so they're flagged and not scored —
honest over tidy.
```

- [ ] **Step 2: Update PLAN.md**

Replace line 27 (`- [ ] Add gnomAD frequency annotation → flag common variants (almost all benign) for triage`) with:

```markdown
- [x] Population-frequency context on trait-xref hits (harvested from ClinVar's
      ExAC/1000G/ESP fields + GWAS catalog RAF — local gnomAD is infeasible here:
      ~72 GB for chr1 alone; offline dbSNP ALFA (19.7 GB) is the future-work
      upgrade path if richer coverage is ever needed; never per-variant APIs,
      which would leak carried rsIDs)
```

- [ ] **Step 3: Full-suite + shellcheck verification**

Run: `./test/run-tests.sh 2>&1 | tail -3` — expected `0 failed`.
Run: `shellcheck -s bash fetch-gwas.sh lib/trait-xref.sh 2>/dev/null || true` — no new warnings versus main (shellcheck may not be installed; skip if absent).

- [ ] **Step 4: Commit**

```bash
git add README.md PLAN.md
git commit -m "docs: population-frequency context in README; tick PLAN.md Phase 2 frequency item"
```
