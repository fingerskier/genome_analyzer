# Trait Cross-Reference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cross-reference a SNP/indel PASS VCF against the GWAS Catalog and ClinVar, reporting which trait/risk and pathogenic alleles the sample carries, with zygosity.

**Architecture:** A standalone entry point `xref-traits.sh` validates the input and dispatches to `lib/trait-xref.sh::run_trait_xref`, which does one `bcftools norm -m-` pass to extract carried biallelic calls, then joins them against a normalized GWAS table (by rsID, p < 5e-8) and ClinVar (by position+allele, pathogenic-class). Two one-time downloaders (`fetch-gwas.sh`, `fetch-clinvar.sh`) cache the databases to `data/`. All allele-dosage logic lives in one reused awk function (`DOSAGE_FN`) so it is unit-testable and DRY.

**Tech Stack:** bash 3.2 (macOS), bcftools/tabix (htslib), awk, curl, unzip. No bedtools. Zero-dependency bash test harness in `test/`.

**Spec:** `docs/superpowers/specs/2026-06-04-trait-xref-design.md`

**Conventions to match (read these first):** `common.sh` (`log`/`die`/`have`/`validate_index`), `fetch-genes.sh` (the `--source-only` guard + testable parse function + `main`), `lib/cnv.sh` (globals `INPUT OUTDIR BASE`, summary heredoc, precomputed awk row strings), `test/lib.sh` (`assert_eq`/`assert_contains`/`make_fixture`), `test/run-tests.sh` (each `test_*.sh` sourced with `set +e -uo pipefail`).

---

### Task 1: Allele-dosage engine (`DOSAGE_FN` + `risk_dosage`)

The crux of correctness. One awk function computes copies-of-risk-allele with strand handling; a bash wrapper makes it unit-testable.

**Files:**
- Create: `lib/trait-xref.sh`
- Test: `test/test_14_dosage.sh`

- [ ] **Step 1: Write the failing test**

Create `test/test_14_dosage.sh`:

```bash
# test_14_dosage.sh — allele dosage + strand logic (lib/trait-xref.sh).
source "$REPO/lib/trait-xref.sh"

# risk_dosage <ref> <alt> <het|hom> <risk> -> "<dosage>\t<flag>"
assert_eq "1	ok"             "$(risk_dosage A G het G)" "dosage: het, risk=alt -> 1 copy"
assert_eq "1	ok"             "$(risk_dosage A G het A)" "dosage: het, risk=ref -> 1 copy"
assert_eq "2	ok"             "$(risk_dosage A G hom G)" "dosage: hom-alt, risk=alt -> 2 copies"
assert_eq "0	ok"             "$(risk_dosage A G hom A)" "dosage: hom-alt, risk=ref -> 0 copies"
assert_eq "1	strand_flipped" "$(risk_dosage A G het C)" "dosage: risk=C matches alt G via complement (non-palindromic)"
assert_eq "NA	ambiguous"      "$(risk_dosage A T het A)" "dosage: A/T palindromic SNP -> ambiguous, not scored"
assert_eq "NA	ambiguous"      "$(risk_dosage C G hom C)" "dosage: C/G palindromic SNP -> ambiguous, not scored"
assert_eq "NA	allele_mismatch" "$(risk_dosage A G het T)" "dosage: risk T absent at A/G site (comp=A=ref) -> wait"
```

Note: the last assertion is intentionally wrong and will be corrected in Step 3's verification — replace it with the real expectation once you confirm behavior. The correct expectation: `risk=T` at an `A/G` site — `T` is not ref/alt; `comp(T)=A` which equals ref, and `A/G` is **not** palindromic, so it resolves as `strand_flipped` with dosage `1` (het). Use this corrected line instead:

```bash
assert_eq "1	strand_flipped" "$(risk_dosage A G het T)" "dosage: risk=T resolves to ref A via complement -> 1 copy, flipped"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_14`
Expected: FAIL — `risk_dosage: command not found` (file does not exist yet).

- [ ] **Step 3: Write the implementation**

Create `lib/trait-xref.sh` (engine added in Task 4; this task adds only the dosage primitives):

```bash
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
  if(risk==ref || risk==alt){ return (gt=="hom" ? (risk==alt?2:0) : 1) SUBSEP "ok" }
  if(length(ref)==1 && length(alt)==1 && length(risk)==1){
    if(comp(ref)==alt){ return "NA" SUBSEP "ambiguous" }
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_14`
Expected: all `test_14` assertions `ok`. If the `allele_mismatch` case is hard to construct (every single base has a complement among A/C/G/T), that is expected — the corrected Step 1 line uses `strand_flipped`, and a true `allele_mismatch` arises only at multi-base (indel) sites; no unit assertion is required for it.

- [ ] **Step 5: Commit**

```bash
git add lib/trait-xref.sh test/test_14_dosage.sh
git commit -m "feat: allele-dosage engine for trait xref (DOSAGE_FN + risk_dosage)"
```

---

### Task 2: GWAS Catalog downloader + normalizer (`fetch-gwas.sh`)

**Files:**
- Create: `fetch-gwas.sh`
- Test: `test/test_13_gwas_normalize.sh`
- Fixture: `test/fixtures/gwas_raw.tsv`

- [ ] **Step 1: Write the fixture**

Create `test/fixtures/gwas_raw.tsv` — a 1-row header (only the columns we read need realistic values; pad others with `x`) plus four data rows. Build it so columns land at the real indices: 2=PUBMEDID, 8=DISEASE/TRAIT, 15=MAPPED_GENE, 21=STRONGEST SNP-RISK ALLELE, 28=P-VALUE, 31=OR or BETA. Use a tab-filled template; the simplest reliable way is an awk generator, but a literal file is clearer. Create it with this exact content (tab-separated; 31 columns each):

```bash
gen() { awk -v rsa="$1" -v trait="$2" -v pv="$3" -v ob="$4" -v gene="$5" -v pmid="$6" 'BEGIN{
  OFS="\t"; for(i=1;i<=31;i++) f[i]="x";
  f[2]=pmid; f[8]=trait; f[15]=gene; f[21]=rsa; f[28]=pv; f[31]=ob;
  line=f[1]; for(i=2;i<=31;i++) line=line OFS f[i]; print line }'; }
{
  # header row
  awk 'BEGIN{OFS="\t"; n=31; h[2]="PUBMEDID";h[8]="DISEASE/TRAIT";h[15]="MAPPED_GENE";
    h[21]="STRONGEST SNP-RISK ALLELE";h[28]="P-VALUE";h[31]="OR or BETA";
    l="COL1"; for(i=2;i<=n;i++) l=l OFS (h[i]?h[i]:"x"); print l}'
  gen "rs100-G" "Height"          "1e-12" "1.2" "GENEH" "111"   # clean, kept
  gen "rs7; rs8-A" "Multi locus"  "2e-12" "1.3" "GENEM" "222"   # multi-SNP -> dropped
  gen "rs300-?" "Unknown allele"  "3e-12" "1.1" "GENEU" "333"   # '?' allele -> dropped
  gen "rs999-AT" "Indel risk"     "4e-12" "1.0" "GENEI" "444"   # 2-base allele -> dropped
} > "$PWD/test/fixtures/gwas_raw.tsv"
```

Run that snippet once from the repo root to generate the file, then verify it has 5 lines (`wc -l test/fixtures/gwas_raw.tsv` → 5) and commit the resulting file (not the generator).

- [ ] **Step 2: Write the failing test**

Create `test/test_13_gwas_normalize.sh`:

```bash
# test_13_gwas_normalize.sh — fetch-gwas.sh::normalize_gwas (no network).
source "$REPO/fetch-gwas.sh" --source-only

out="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw.tsv")"
assert_eq "1" "$(printf '%s\n' "$out" | grep -c .)" "normalize: only the 1 clean row survives"
assert_contains "$out" "rs100	G	Height	1e-12	1.2	GENEH	111" "normalize: clean row -> 7 cols, rsid/allele split on last dash"
# Dropped rows must not appear.
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Multi locus')"   "normalize: multi-SNP row dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Unknown allele')" "normalize: '?' risk allele dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Indel risk')"     "normalize: multi-base risk allele dropped"
```

- [ ] **Step 3: Run test to verify it fails**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_13`
Expected: FAIL — `normalize_gwas: command not found`.

- [ ] **Step 4: Write the implementation**

Create `fetch-gwas.sh`:

```bash
#!/usr/bin/env bash
# fetch-gwas.sh — download GWAS Catalog associations, normalize to a compact
# join-ready table at data/gwas-catalog.tsv: rsid risk trait pval orbeta gene pmid.
# Only single-rsID rows with a single-base risk allele are kept (v1 scope).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# "ontology-annotated, alternative, full" associations zip. Override with GWAS_URL.
GWAS_URL="${GWAS_URL:-https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip}"

# normalize_gwas: GWAS Catalog associations TSV on stdin -> normalized table.
# Columns read (1-based): 2 PUBMEDID, 8 DISEASE/TRAIT, 15 MAPPED_GENE,
# 21 STRONGEST SNP-RISK ALLELE (e.g. "rs1234-A"), 28 P-VALUE, 31 OR or BETA.
normalize_gwas() {
  awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1{next}
    {
      sra=$21
      p=0; for(i=length(sra);i>=1;i--){ if(substr(sra,i,1)=="-"){p=i; break} }
      if(p==0) next
      rsid=substr(sra,1,p-1); risk=substr(sra,p+1)
      if(rsid !~ /^rs[0-9]+$/) next       # single rsID only (drops "rs7; rs8")
      if(risk !~ /^[ACGT]$/)   next       # single-base risk allele only (drops "?", indels)
      trait=$8; gsub(/[\t\r]/," ",trait)
      print rsid, risk, trait, $28, $31, $15, $2
    }'
}

main() {
  local out="$SCRIPT_DIR/data/gwas-catalog.tsv"
  mkdir -p "$SCRIPT_DIR/data"
  if [[ -s "$out" && "${1:-}" != "-f" && "${1:-}" != "--force" ]]; then
    echo "already present: $out ($(wc -l < "$out" | tr -d ' ') rows). Use -f to refresh." >&2
    return 0
  fi
  command -v curl  >/dev/null 2>&1 || { echo "ERROR: curl required"  >&2; exit 1; }
  command -v unzip >/dev/null 2>&1 || { echo "ERROR: unzip required" >&2; exit 1; }
  echo "downloading $GWAS_URL" >&2
  local zip; zip="$(mktemp)"
  curl -fsSL "$GWAS_URL" -o "$zip"
  unzip -p "$zip" | normalize_gwas | LC_ALL=C sort -k1,1 > "$out"
  rm -f "$zip"
  echo "wrote $out ($(wc -l < "$out" | tr -d ' ') associations)" >&2
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
```

- [ ] **Step 5: Run test to verify it passes**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_13`
Expected: all `test_13` assertions `ok`.

- [ ] **Step 6: Commit**

```bash
git add fetch-gwas.sh test/test_13_gwas_normalize.sh test/fixtures/gwas_raw.tsv
git commit -m "feat: fetch-gwas.sh downloader + normalize_gwas with unit test"
```

---

### Task 3: ClinVar downloader (`fetch-clinvar.sh`)

Thin network downloader. Smoke-tested for sourceability and config only (no network in tests).

**Files:**
- Create: `fetch-clinvar.sh`
- Test: append to `test/test_13_gwas_normalize.sh` is wrong — create `test/test_15a_fetch_clinvar.sh`

- [ ] **Step 1: Write the failing test**

Create `test/test_15a_fetch_clinvar.sh`:

```bash
# test_15a_fetch_clinvar.sh — fetch-clinvar.sh sources cleanly + sane defaults.
source "$REPO/fetch-clinvar.sh" --source-only
assert_contains "$CLINVAR_URL" "clinvar.vcf.gz" "fetch-clinvar: default URL points at clinvar.vcf.gz"
assert_contains "$CLINVAR_URL" "vcf_GRCh38"     "fetch-clinvar: uses the GRCh38 build"
assert_eq "0" "$(bash -n "$REPO/fetch-clinvar.sh"; echo $?)" "fetch-clinvar: passes bash -n syntax check"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_15a`
Expected: FAIL — file not found / `CLINVAR_URL` unbound.

- [ ] **Step 3: Write the implementation**

Create `fetch-clinvar.sh`:

```bash
#!/usr/bin/env bash
# fetch-clinvar.sh — download NCBI ClinVar GRCh38 VCF (+ index) to data/.
# Contigs are bare names (1,2,...,X,Y,MT), matching Sequencing.com VCFs.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CLINVAR_URL="${CLINVAR_URL:-https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz}"

main() {
  local out="$SCRIPT_DIR/data/clinvar.GRCh38.vcf.gz"
  mkdir -p "$SCRIPT_DIR/data"
  if [[ -s "$out" && "${1:-}" != "-f" && "${1:-}" != "--force" ]]; then
    echo "already present: $out. Use -f to refresh." >&2
    return 0
  fi
  command -v curl  >/dev/null 2>&1 || { echo "ERROR: curl required"        >&2; exit 1; }
  command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix (htslib) required" >&2; exit 1; }
  echo "downloading $CLINVAR_URL" >&2
  curl -fsSL "$CLINVAR_URL" -o "$out"
  if ! curl -fsSL "${CLINVAR_URL}.tbi" -o "${out}.tbi" 2>/dev/null; then
    echo "no published .tbi; building index" >&2
    tabix -p vcf "$out"
  fi
  echo "wrote $out" >&2
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_15a`
Expected: all `test_15a` assertions `ok`.

- [ ] **Step 5: Commit**

```bash
git add fetch-clinvar.sh test/test_15a_fetch_clinvar.sh
git commit -m "feat: fetch-clinvar.sh downloader (GRCh38 VCF + index)"
```

---

### Task 4: Cross-reference engine (`run_trait_xref`)

Adds the main function to `lib/trait-xref.sh`: one norm pass, GWAS join (rsID, p<5e-8), ClinVar join (position+allele, pathogenic-class), SUMMARY emit, and graceful degradation.

**Files:**
- Modify: `lib/trait-xref.sh` (append `run_trait_xref`)
- Test: `test/test_15_xref.sh`
- Fixtures: `test/fixtures/xref_sample.vcf.txt`, `test/fixtures/gwas_norm.tsv`, `test/fixtures/clinvar.vcf.txt`

- [ ] **Step 1: Write the fixtures**

Create `test/fixtures/xref_sample.vcf.txt`:

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="passed">
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	100	rs100	A	G	.	PASS	.	GT	0/1
1	200	rs200	C	T	.	PASS	.	GT	1/1
1	300	rs300	A	T	.	PASS	.	GT	0/1
2	500	.	G	C	.	PASS	.	GT	1/1
```

Create `test/fixtures/gwas_norm.tsv` (already-normalized: rsid risk trait pval orbeta gene pmid — TAB separated):

```
rs100	G	Height	1e-12	1.2	GENEH	111
rs200	T	Coffee consumption	3e-20	1.5	GENEC	222
rs200	A	Weak association	1e-6	0.9	GENEC	333
rs300	A	Palindromic trait	1e-9	1.1	GENEP	444
```

Create `test/fixtures/clinvar.vcf.txt`:

```
##fileformat=VCFv4.2
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="significance">
##INFO=<ID=CLNDN,Number=.,Type=String,Description="disease name">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="review status">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="gene">
##contig=<ID=1>
##contig=<ID=2>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
2	500	900	G	C	.	.	CLNSIG=Pathogenic;CLNDN=Example_condition;CLNREVSTAT=criteria_provided;GENEINFO=GENEX:999
1	100	901	A	G	.	.	CLNSIG=Benign;CLNDN=Nothing_notable;CLNREVSTAT=criteria_provided;GENEINFO=GENEY:888
```

- [ ] **Step 2: Write the failing test**

Create `test/test_15_xref.sh`:

```bash
# test_15_xref.sh — end-to-end trait xref on fixtures (GWAS + ClinVar).
( source "$REPO/common.sh"; source "$REPO/lib/trait-xref.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/xref_sample.vcf.txt")"
  CLINVAR_VCF="$(make_fixture "$REPO/test/fixtures/clinvar.vcf.txt")"
  GWAS_TBL="$REPO/test/fixtures/gwas_norm.tsv"
  OUTDIR="$TMP/xref_out"; mkdir -p "$OUTDIR"; BASE="xref"
  run_trait_xref )

sum="$(cat "$TMP/xref_out/xref.traits.SUMMARY.md")"
gtsv="$(cat "$TMP/xref_out/xref.traits.gwas.tsv")"
ctsv="$(cat "$TMP/xref_out/xref.traits.clinvar.tsv")"

# GWAS: rs100 het risk=G -> 1 copy Height ; rs200 hom risk=T -> 2 copies Coffee.
assert_contains "$gtsv" "Height"            "gwas: rs100 Height association reported"
assert_contains "$gtsv" "Coffee consumption" "gwas: rs200 Coffee association reported"
assert_eq "0" "$(printf '%s\n' "$gtsv" | grep -c 'Weak association')" "gwas: p=1e-6 below 5e-8 threshold excluded"
assert_contains "$sum" "**T/T**"            "gwas summary: homozygous risk genotype highlighted"
assert_contains "$sum" "strand-ambiguous"   "gwas summary: mentions the ambiguous-SNP bucket"
# rs300 is A/T palindromic -> ambiguous, present in tsv with flag, not scored.
assert_contains "$gtsv" "ambiguous"         "gwas tsv: palindromic rs300 flagged ambiguous"

# ClinVar: chr2:500 G>C hom matches Pathogenic; chr1:100 Benign must be excluded.
assert_contains "$ctsv" "Example condition" "clinvar: pathogenic hit reported (underscores->spaces)"
assert_contains "$ctsv" "homozygous"        "clinvar: zygosity reported"
assert_contains "$ctsv" "GENEX"             "clinvar: gene symbol reported"
assert_eq "0" "$(printf '%s\n' "$ctsv" | grep -c 'Nothing notable')" "clinvar: benign variant excluded"
assert_contains "$sum" "not a diagnosis"    "summary: carries research-grade caveat"
```

- [ ] **Step 3: Run test to verify it fails**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_15`
Expected: FAIL — `run_trait_xref: command not found`.

- [ ] **Step 4: Write the implementation**

Append to `lib/trait-xref.sh`:

```bash
# run_trait_xref: cross-reference the sample VCF ($INPUT) against the GWAS Catalog
# ($GWAS_TBL) and ClinVar ($CLINVAR_VCF). Writes three artifacts under $OUTDIR and
# prints the SUMMARY path. Degrades gracefully when a database file is absent.
run_trait_xref() {
  local SUMMARY="$OUTDIR/${BASE}.traits.SUMMARY.md"
  local GTSV="$OUTDIR/${BASE}.traits.gwas.tsv"
  local CTSV="$OUTDIR/${BASE}.traits.clinvar.tsv"
  local tmp; tmp="$(mktemp -d)"

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
      FNR==NR{                                                  # file1: GWAS table
        if($4 ~ /[0-9]/ && ($4+0)<THRESH){
          rec=$2 SUBSEP $3 SUBSEP $4 SUBSEP $5 SUBSEP $6 SUBSEP $7
          g[$1]=($1 in g)? g[$1] "\x1f" rec : rec
        }
        next
      }
      ($3 in g){                                                # file2: carried.tsv
        rs=$3; ref=$4; alt=$5; z=$6
        geno=(z=="hom"? alt"/"alt : ref"/"alt)
        n=split(g[rs], aa, "\x1f")
        for(i=1;i<=n;i++){
          split(aa[i], f, SUBSEP)                              # f1 risk f2 trait f3 pval f4 or f5 gene f6 pmid
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
          FNR==NR{                                              # file1: clinvar path rows
            k=$1":"$2":"$3":"$4
            g=$8; gsub(/:[0-9]+/,"",g); gsub(/\|/,", ",g)
            sig[k]=$5; dn[k]=$6; rev[k]=$7; gene[k]=g; next
          }
          { k=$1":"$2":"$4":"$5                                 # file2: carried.tsv
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
  clinvar_rows="$(awk -F'\t' 'NR<=50{printf "| %s | %s | %s | %s | %s |\n", $5,$7,$3,$4,$6}' "$CTSV")"
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

  rm -rf "$tmp"
  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
```

- [ ] **Step 5: Run test to verify it passes**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_15`
Expected: all `test_15` assertions `ok`. If `**T/T**` is not found, check that rs200's ALT is `T` and GT `1/1` so `geno=alt"/"alt="T/T"` and `$5==2`.

- [ ] **Step 6: Commit**

```bash
git add lib/trait-xref.sh test/test_15_xref.sh test/fixtures/xref_sample.vcf.txt test/fixtures/gwas_norm.tsv test/fixtures/clinvar.vcf.txt
git commit -m "feat: run_trait_xref engine (GWAS + ClinVar join, summary, degradation)"
```

---

### Task 5: Entry point (`xref-traits.sh`) + degradation test

**Files:**
- Create: `xref-traits.sh`
- Test: `test/test_16_xref_entry.sh`

- [ ] **Step 1: Write the failing test**

Create `test/test_16_xref_entry.sh`:

```bash
# test_16_xref_entry.sh — xref-traits.sh CLI + graceful degradation (no DBs).
sample="$(make_fixture "$REPO/test/fixtures/xref_sample.vcf.txt")"

# Missing input -> nonzero exit, usage shown.
out="$("$REPO/xref-traits.sh" 2>&1)"; rc=$?
assert_eq "1" "$rc" "entry: no -i exits 1"
assert_contains "$out" "Usage" "entry: prints usage when -i missing"

# Both DBs absent -> still writes a SUMMARY with skip notes, exit 0.
od="$TMP/entry_out"
"$REPO/xref-traits.sh" -i "$sample" -o "$od" -g "$TMP/nope.tsv" -c "$TMP/nope.vcf.gz" >/dev/null 2>&1
rc=$?
assert_eq "0" "$rc" "entry: missing DBs still exits 0"
sum="$(cat "$od/xref.traits.SUMMARY.md" 2>/dev/null)"
assert_contains "$sum" "fetch-gwas.sh"    "entry: GWAS skip note guides the user to fetch-gwas.sh"
assert_contains "$sum" "fetch-clinvar.sh" "entry: ClinVar skip note guides the user to fetch-clinvar.sh"
```

Note: `BASE` derives from the fixture filename `xref_sample.vcf.gz` → `xref_sample`, so the summary is `xref_sample.traits.SUMMARY.md`. Adjust the two `od/...SUMMARY.md` paths to `xref_sample.traits.SUMMARY.md`:

```bash
sum="$(cat "$od/xref_sample.traits.SUMMARY.md" 2>/dev/null)"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./test/run-tests.sh 2>&1 | grep -A2 test_16`
Expected: FAIL — `xref-traits.sh: No such file or directory`.

- [ ] **Step 3: Write the implementation**

Create `xref-traits.sh` (then `chmod +x`):

```bash
#!/usr/bin/env bash
# xref-traits.sh — cross-reference a SNP/indel PASS VCF against the GWAS Catalog
# and ClinVar, reporting carried trait/risk and pathogenic alleles with zygosity.
# Usage: ./xref-traits.sh -i sample.snp-indel.genome.pass.vcf.gz [-o out]
#                         [-g data/gwas-catalog.tsv] [-c data/clinvar.GRCh38.vcf.gz]
# Research-grade exploration, not a clinical diagnosis.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

usage() { sed -n '2,6p' "$0"; }

main() {
  INPUT=""; OUTDIR="out"
  GWAS_TBL="$SCRIPT_DIR/data/gwas-catalog.tsv"
  CLINVAR_VCF="$SCRIPT_DIR/data/clinvar.GRCh38.vcf.gz"
  while getopts ":i:o:g:c:h" opt; do
    case "$opt" in
      i) INPUT="$OPTARG" ;;
      o) OUTDIR="$OPTARG" ;;
      g) GWAS_TBL="$OPTARG" ;;
      c) CLINVAR_VCF="$OPTARG" ;;
      h) usage; exit 0 ;;
      :) die "option -$OPTARG requires an argument" ;;
      \?) die "unknown option -$OPTARG" ;;
    esac
  done
  [[ -n "$INPUT" ]] || { usage; exit 1; }
  [[ -f "$INPUT" ]] || die "input not found: $INPUT"
  have bcftools || die "bcftools not on PATH (install htslib/bcftools first)"
  have tabix    || die "tabix not on PATH (part of htslib)"

  mkdir -p "$OUTDIR"
  BASE="$(basename "${INPUT%.vcf.gz}")"; BASE="${BASE%.vcf}"
  log "validate + index"
  validate_index "$INPUT"
  # shellcheck disable=SC1091
  source "$SCRIPT_DIR/lib/trait-xref.sh"
  run_trait_xref
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
```

The word "Usage" appears on line 3 of the file, inside the `sed -n '2,6p'` range, so `usage` prints it — satisfying the `assert_contains "$out" "Usage"` check.

- [ ] **Step 4: Make it executable and run the test**

```bash
chmod +x xref-traits.sh
./test/run-tests.sh 2>&1 | grep -A2 test_16
```
Expected: all `test_16` assertions `ok`.

- [ ] **Step 5: Commit**

```bash
git add xref-traits.sh test/test_16_xref_entry.sh
git commit -m "feat: xref-traits.sh entry point + degradation test"
```

---

### Task 6: `.gitignore`, README, and full-suite green

**Files:**
- Modify: `.gitignore`
- Modify: `README.md`

- [ ] **Step 1: Ignore the new personal-data artifacts**

Edit `.gitignore` — after the line `*.events.bed` (line 25), add:

```
*.traits.gwas.tsv
*.traits.clinvar.tsv
```

(`*.traits.SUMMARY.md` is already covered by the existing `*.SUMMARY.md` rule.)

- [ ] **Step 2: Document the feature in README**

In `README.md`, after the "Outputs (in `out/`)" table (ends line 87) and before "## Testing", insert:

````markdown
## Trait cross-reference (GWAS Catalog + ClinVar)

Once you have a SNP/indel PASS VCF (from `analyze.sh`), cross-reference it against
public databases to see which trait/risk and clinically-flagged alleles you carry.

> **This is research-grade exploration, not a diagnosis or medical advice.**
> A "risk allele" is a statistical association, not a verdict. Carrying one does
> not mean you have or will develop a condition.

One-time database setup (cached under git-ignored `data/`):

```bash
./fetch-gwas.sh        # GWAS Catalog associations (~68 MB download) -> data/gwas-catalog.tsv
./fetch-clinvar.sh     # NCBI ClinVar GRCh38 VCF (large) -> data/clinvar.GRCh38.vcf.gz
```

Run the cross-reference:

```bash
./xref-traits.sh -i sample.snp-indel.genome.pass.vcf.gz -o out/
```

| Flag | Meaning |
|---|---|
| `-i FILE` | input SNP/indel PASS `.vcf.gz` (required) |
| `-o DIR` | output directory (default `out`) |
| `-g FILE` | GWAS table (default `data/gwas-catalog.tsv`) |
| `-c FILE` | ClinVar VCF (default `data/clinvar.GRCh38.vcf.gz`) |
| `-h` | help |

Outputs: `*.traits.SUMMARY.md` (readable report — ClinVar pathogenic hits, then
GWAS trait hits you carry ≥1 risk allele for, homozygous hits highlighted),
`*.traits.gwas.tsv`, and `*.traits.clinvar.tsv` (full tables). If a database file
is missing, the run still completes and the summary tells you which fetch script
to run. **Strand-ambiguous SNPs** (A/T, C/G) can't be oriented from the VCF alone,
so they're flagged and not scored — honest over tidy.
````

- [ ] **Step 3: Run the entire suite**

Run: `./test/run-tests.sh`
Expected: final line `N passed, 0 failed` (N = previous 41 + the new assertions).

- [ ] **Step 4: Commit**

```bash
git add .gitignore README.md
git commit -m "docs: document trait cross-reference; ignore traits TSV artifacts"
```

---

## Self-Review

**Spec coverage:**
- GWAS Catalog rsID xref, p<5e-8 → Task 2 (normalize) + Task 4 (join). ✓
- ClinVar position+allele xref, pathogenic-class → Task 3 (fetch) + Task 4 (join). ✓
- Risk-allele dosage with strand/palindrome handling → Task 1 (`DOSAGE_FN`), reused in Task 4. ✓
- One-time downloaders mirroring `fetch-genes.sh` → Tasks 2, 3. ✓
- Separate `xref-traits.sh` entry point → Task 5. ✓
- Readable SUMMARY + two TSVs → Task 4. ✓
- Graceful degradation when DB absent → Task 4 (skip notes) + Task 5 (test). ✓
- `.gitignore` for personal-data artifacts → Task 6. ✓
- README section, research-grade caveats → Task 6 + Task 4 summary footer. ✓
- Non-goals (PRS, multi-SNP, ancestry) → excluded; documented as future work in the summary. ✓

**Type/identifier consistency:**
- Globals `INPUT OUTDIR BASE GWAS_TBL CLINVAR_VCF` set in `xref-traits.sh::main`, consumed in `run_trait_xref` — names match. ✓
- `DOSAGE_FN` defined once, used by `risk_dosage` (Task 1) and the GWAS join awk (Task 4). ✓
- normalized GWAS columns `rsid risk trait pval orbeta gene pmid` produced by `normalize_gwas` (Task 2) and consumed positionally by the join (Task 4: `$1..$7`, packing `$2..$7`). ✓
- `carried.tsv` columns `chrom pos id ref alt zyg` produced once (Task 4 step 1) and consumed by both joins (GWAS keys on `$3` id; ClinVar keys on `$1:$2:$4:$5`). ✓
- Output filenames `${BASE}.traits.{SUMMARY.md,gwas.tsv,clinvar.tsv}` consistent across engine, entry test, README, `.gitignore`. ✓

**Placeholder scan:** No TBD/TODO; every code step has complete content. The one intentionally-wrong test line in Task 1 Step 1 is explicitly corrected in the same step with the real expectation and rationale. ✓
