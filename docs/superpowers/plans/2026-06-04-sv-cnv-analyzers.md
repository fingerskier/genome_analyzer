# SV & CNV Analyzers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an auto-detecting `analyze.sh` dispatcher that routes SNP/indel, Manta SV, and Canvas CNV VCFs to per-class analyzers, each emitting a `SUMMARY.md` with descriptive stats plus (for SV/CNV) gene overlap.

**Architecture:** A thin `analyze.sh` entry point validates/indexes the input, sniffs its class from the VCF header/records, and dispatches to `lib/{snp-indel,sv,cnv}.sh`. Shared helpers (logging, validation, gene overlap, summary scaffold) live in `common.sh`. Gene overlap uses `bedtools` against a GENCODE-basic GRCh38 BED cached by `fetch-genes.sh`; it degrades to stats-only when those are absent. All scripts target macOS bash 3.2.

**Tech Stack:** bash 3.2, bcftools/tabix/bgzip (htslib), bedtools (optional), awk. Zero-dependency plain-bash test harness.

---

## File Structure

| File | Responsibility |
|---|---|
| `analyze.sh` | Entry point: arg parse, validate+index, detect class, dispatch |
| `common.sh` | Sourced helpers: `log`/`die`/`have`, `validate_index`, `gene_overlap`, summary scaffold |
| `lib/snp-indel.sh` | Existing SNP/indel logic (moved from `analyze-vcf.sh`) |
| `lib/sv.sh` | Manta SV analyzer |
| `lib/cnv.sh` | Canvas CNV analyzer |
| `fetch-genes.sh` | One-time GENCODE-basic GRCh38 → `data/genes.GRCh38.bed` (chr-normalized) |
| `analyze-vcf.sh` | Compatibility shim → `exec analyze.sh "$@"` |
| `README.md` | Install + usage docs (rewritten) |
| `test/lib.sh` | Test harness: assertions + `make_fixture` |
| `test/run-tests.sh` | Test runner: sets up `$TMP`, sources every `test_*.sh` |
| `test/test_*.sh` | Per-component tests |
| `test/fixtures/*.vcf.txt` | Plain-text VCF fixtures (bgzipped at test time) |
| `test/fixtures/genes.bed` | Tiny gene BED for overlap tests |

**Conventions used throughout:**
- Every script starts with `#!/usr/bin/env bash` and `set -euo pipefail`.
- `lib/*.sh` and `common.sh` are **sourced**, not executed; they define functions and read globals (`INPUT`, `OUTDIR`, `BASE`, `THREADS`, `GENES_BED`).
- `SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"` resolves sibling files regardless of caller cwd.

---

## Task 0: Test harness

**Files:**
- Create: `test/lib.sh`
- Create: `test/run-tests.sh`

- [ ] **Step 1: Write the harness library**

Create `test/lib.sh`:

```bash
# test/lib.sh — zero-dependency assertions + fixture builder. Sourced by run-tests.sh.
PASS=0; FAIL=0

assert_eq() {  # expected actual message
  if [[ "$1" == "$2" ]]; then PASS=$((PASS+1)); printf '  ok   %s\n' "$3"
  else FAIL=$((FAIL+1)); printf '  FAIL %s\n    expected: [%s]\n    actual:   [%s]\n' "$3" "$1" "$2"; fi
}

assert_contains() {  # haystack needle message
  if printf '%s' "$1" | grep -qF -- "$2"; then PASS=$((PASS+1)); printf '  ok   %s\n' "$3"
  else FAIL=$((FAIL+1)); printf '  FAIL %s\n    missing: [%s]\n' "$3" "$2"; fi
}

# make_fixture <fixtures/name.vcf.txt> -> prints path to a bgzipped+indexed VCF in $TMP
make_fixture() {
  local src="$1" base out
  base="$(basename "${src%.vcf.txt}")"
  out="$TMP/${base}.vcf.gz"
  bgzip -c "$src" > "$out"
  tabix -f -p vcf "$out"
  printf '%s' "$out"
}

finish() { printf '\n%d passed, %d failed\n' "$PASS" "$FAIL"; [[ "$FAIL" -eq 0 ]]; }
```

- [ ] **Step 2: Write the runner**

Create `test/run-tests.sh`:

```bash
#!/usr/bin/env bash
set -uo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$DIR/.." && pwd)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
export REPO TMP
# shellcheck disable=SC1091
source "$DIR/lib.sh"
for t in "$DIR"/test_*.sh; do
  [[ -e "$t" ]] || continue
  printf '\n=== %s ===\n' "$(basename "$t")"
  # shellcheck disable=SC1090
  source "$t"
done
finish
```

- [ ] **Step 3: Make runner executable and run it (empty pass)**

Run: `chmod +x test/run-tests.sh && ./test/run-tests.sh`
Expected: prints `0 passed, 0 failed` and exits 0 (no `test_*.sh` yet).

- [ ] **Step 4: Commit**

```bash
git add test/lib.sh test/run-tests.sh
git commit -m "test: add zero-dependency bash test harness"
```

---

## Task 1: Fixtures

**Files:**
- Create: `test/fixtures/snp-indel.vcf.txt`
- Create: `test/fixtures/sv.vcf.txt`
- Create: `test/fixtures/cnv.vcf.txt`
- Create: `test/fixtures/genes.bed`

These are committed as `.vcf.txt` (not `.vcf.gz`) so no genomic-data gitignore rule hits them and nothing private is committed. Counts are chosen so later tests can assert exact numbers.

- [ ] **Step 1: SNP/indel fixture** — `test/fixtures/snp-indel.vcf.txt`:

```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1>
##contig=<ID=2>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	1000	.	A	G	100	PASS	.	GT	0/1
1	2000	.	C	T	100	PASS	.	GT	1/1
2	3000	.	GA	G	100	PASS	.	GT	0/1
```

- [ ] **Step 2: SV fixture (Manta-like)** — `test/fixtures/sv.vcf.txt`:

```
##fileformat=VCFv4.1
##source=GenerateSVCandidates 1.6.0
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakend">
##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=MinQUAL,Description="QUAL score is less than 20">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	5000	D1	N	<DEL>	200	PASS	SVTYPE=DEL;SVLEN=-500;END=5500	GT	0/1
1	60000	D2	N	<DEL>	200	PASS	SVTYPE=DEL;SVLEN=-20000;END=80000	GT	0/1
2	7000	U1	N	<DUP:TANDEM>	200	PASS	SVTYPE=DUP;SVLEN=1500;END=8500	GT	0/1
2	9000	I1	N	<INS>	200	PASS	SVTYPE=INS;SVLEN=300;END=9001	GT	0/1
3	1000	X1	N	<DEL>	10	MinQUAL	SVTYPE=DEL;SVLEN=-400;END=1400	GT	0/1
1	40000	B1	N	N[5:90000[	200	PASS	SVTYPE=BND;MATEID=B2;EVENT=E1	GT	0/1
5	90000	B2	N	]1:40000]N	200	PASS	SVTYPE=BND;MATEID=B1;EVENT=E1	GT	0/1
```

Known truths: PASS records = 6; PASS non-BND = 4 (DEL 2, DUP 1, INS 1); BND events = 1 (EVENT=E1, 2 records); non-PASS = 1 (MinQUAL).

- [ ] **Step 3: CNV fixture (Canvas-like)** — `test/fixtures/cnv.vcf.txt`:

```
##fileformat=VCFv4.1
##source=Canvas 1.40.0.1613+master
##ALT=<ID=CN0,Description="Copy number 0">
##ALT=<ID=CN1,Description="Copy number 1">
##ALT=<ID=CN3,Description="Copy number 3">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=CNVLEN,Number=1,Type=Integer,Description="Number of reference positions spanned">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=L10kb,Description="Length shorter than 10 kb">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
1	800	Canvas:REF:1:800-900	N	.	20	PASS	END=900	GT:CN	./.:2
1	1000	Canvas:LOSS:1:1000-51000	N	<CN0>	30	PASS	SVTYPE=CNV;CNVLEN=50000;END=51000	GT:CN	1/1:0
2	2000	Canvas:LOSS:2:2000-22000	N	<CN1>	30	PASS	SVTYPE=CNV;CNVLEN=20000;END=22000	GT:CN	0/1:1
3	60000	Canvas:GAIN:3:60000-160000	N	<CN3>	30	PASS	SVTYPE=CNV;CNVLEN=100000;END=160000	GT:CN	0/1:3
4	5000	Canvas:GAIN:4:5000-13000	N	<CN3>	5	L10kb	SVTYPE=CNV;CNVLEN=8000;END=13000	GT:CN	0/1:3
```

Known truths: non-REF events = 4; PASS events = 3 (losses 2, gain 1); CN0 homozygous deletion = 1; non-PASS = 1 (L10kb); total PASS loss bases = 70000, gain bases = 100000.

- [ ] **Step 4: Gene BED fixture** — `test/fixtures/genes.bed` (bare contig names, tab-separated):

```
1	100	100000	GENEA
2	100	100000	GENEB
3	50000	200000	GENEC
5	80000	95000	GENEE
```

- [ ] **Step 5: Verify fixtures bgzip/index cleanly**

Add `test/test_00_fixtures.sh`:

```bash
for f in snp-indel sv cnv; do
  vcf="$(make_fixture "$REPO/test/fixtures/$f.vcf.txt")"
  hdr="$(bcftools view -h "$vcf" 2>&1)"
  assert_contains "$hdr" "##fileformat" "fixture $f is a readable VCF"
done
```

Run: `./test/run-tests.sh`
Expected: 3 ok, `3 passed, 0 failed`.

- [ ] **Step 6: Commit**

```bash
git add test/fixtures test/test_00_fixtures.sh
git commit -m "test: add synthetic SV/CNV/snp-indel and gene fixtures"
```

---

## Task 2: common.sh helpers

**Files:**
- Create: `common.sh`
- Test: `test/test_01_common.sh`

- [ ] **Step 1: Write the failing test** — `test/test_01_common.sh`:

```bash
source "$REPO/common.sh"
vcf="$(make_fixture "$REPO/test/fixtures/snp-indel.vcf.txt")"
rm -f "${vcf}.tbi"          # force re-index path
validate_index "$vcf"
[[ -f "${vcf}.tbi" ]] && r=yes || r=no
assert_eq "yes" "$r" "validate_index creates a .tbi"
have bcftools && hb=yes || hb=no
assert_eq "yes" "$hb" "have detects bcftools on PATH"
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — `common.sh` does not exist / `validate_index: command not found`.

- [ ] **Step 3: Write `common.sh`**

```bash
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh`
Expected: the two new assertions pass.

- [ ] **Step 5: Commit**

```bash
git add common.sh test/test_01_common.sh
git commit -m "feat: add common.sh with log/die/have/validate_index"
```

---

## Task 3: lib/snp-indel.sh (move existing logic)

**Files:**
- Create: `lib/snp-indel.sh`
- Test: `test/test_02_snp_indel.sh`

Move the current `analyze-vcf.sh` stages 2–5 + summary into a sourced `run_snp_indel` function. Behavior is unchanged; it now relies on `common.sh` for `log`/`die`/`have`/`validate_index` and reads globals `INPUT`, `OUTDIR`, `BASE`, `THREADS`, `ANNOTATE`, `BUILD`, `REF`.

- [ ] **Step 1: Write the failing test** — `test/test_02_snp_indel.sh`:

```bash
( source "$REPO/common.sh"
  source "$REPO/lib/snp-indel.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/snp-indel.vcf.txt")"
  OUTDIR="$TMP/snp_out"; mkdir -p "$OUTDIR"
  BASE="snp-indel"; THREADS=1; ANNOTATE=none; BUILD=GRCh38; REF=""
  run_snp_indel
)
sum="$TMP/snp_out/snp-indel.SUMMARY.md"
[[ -s "$sum" ]] && e=yes || e=no
assert_eq "yes" "$e" "snp-indel emits a non-empty SUMMARY.md"
assert_contains "$(cat "$sum")" "SNPs" "snp-indel summary reports SNP count"
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — `lib/snp-indel.sh` not found.

- [ ] **Step 3: Write `lib/snp-indel.sh`**

Wrap the existing stage 2–5 logic in a function. Copy the body from `analyze-vcf.sh:73-158` verbatim except: drop the per-stage `Stage N/5` numbering is optional, keep messages. Full content:

```bash
# lib/snp-indel.sh — SNP/indel analyzer. Sourced by analyze.sh. Expects globals:
# INPUT OUTDIR BASE THREADS ANNOTATE BUILD REF ; helpers from common.sh.

run_snp_indel() {
  local STATS="$OUTDIR/${BASE}.stats.txt"
  local PASS="$OUTDIR/${BASE}.pass.vcf.gz"
  local ANNO="$OUTDIR/${BASE}.annotated.vcf.gz"
  local SUMMARY="$OUTDIR/${BASE}.SUMMARY.md"

  log "stats"
  bcftools stats "$INPUT" > "$STATS"
  local SAMPLES NSAMP RECORDS SNPS INDELS TITV
  SAMPLES="$(bcftools query -l "$INPUT" | paste -sd, -)"
  NSAMP="$(bcftools query -l "$INPUT" | wc -l | tr -d ' ')"
  get_sn() { awk -F'\t' -v k="$1" '$1=="SN" && $3==k {print $4; exit}' "$STATS"; }
  RECORDS="$(get_sn 'number of records:')"
  SNPS="$(get_sn 'number of SNPs:')"
  INDELS="$(get_sn 'number of indels:')"
  TITV="$(awk -F'\t' '$1=="TSTV"{print $5; exit}' "$STATS")"
  : "${RECORDS:=0}"; : "${SNPS:=0}"; : "${INDELS:=0}"; : "${TITV:=NA}"

  local ASSAY TITV_OK
  if   (( SNPS > 1000000 )); then ASSAY="Whole genome (WGS)"; TITV_OK="~2.0-2.1"
  elif (( SNPS > 20000   )); then ASSAY="Whole exome (WES) / large panel"; TITV_OK="~3.0-3.3"
  else                            ASSAY="Targeted panel / genotyping array"; TITV_OK="varies"
  fi
  log "assay looks like: $ASSAY  (SNPs=$SNPS, Ti/Tv=$TITV)"

  log "PASS-only filter"
  bcftools view -f PASS,. "$INPUT" -Oz -o "$PASS"
  tabix -f -p vcf "$PASS"
  local PASS_N
  PASS_N="$(bcftools index -n "$PASS" 2>/dev/null || bcftools stats "$PASS" | awk -F'\t' '$1=="SN" && $3=="number of records:"{print $4}')"

  local FINAL="$PASS" ANNO_NOTE="skipped (run with -a snpeff|vep)"
  case "$ANNOTATE" in
    none) log "annotation: skipped" ;;
    snpeff)
      have snpeff || die "snpeff not on PATH"
      log "SnpEff annotation ($BUILD)"
      snpeff -noStats "$BUILD" "$PASS" | bgzip > "$ANNO"
      tabix -f -p vcf "$ANNO"; FINAL="$ANNO"; ANNO_NOTE="SnpEff / $BUILD" ;;
    vep)
      have vep || die "vep not on PATH"
      log "VEP annotation ($BUILD)"
      local REF_ARG=(); [[ -n "$REF" ]] && REF_ARG=(--fasta "$REF" --hgvs)
      vep -i "$PASS" --cache --offline --assembly "$BUILD" --everything \
          --vcf --compress_output bgzip --fork "$THREADS" \
          "${REF_ARG[@]}" -o "$ANNO"
      tabix -f -p vcf "$ANNO"; FINAL="$ANNO"; ANNO_NOTE="Ensembl VEP / $BUILD" ;;
    *) die "unknown annotation engine: $ANNOTATE (use none|snpeff|vep)" ;;
  esac

  log "writing $SUMMARY"
  cat > "$SUMMARY" <<EOF
# VCF Summary — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by analyze.sh (snp-indel)_

## Overview
| Field | Value |
|---|---|
| Source file | \`$(basename "$INPUT")\` |
| Sample(s) | $NSAMP — \`${SAMPLES}\` |
| Likely assay | **$ASSAY** |
| Total records | $RECORDS |
| SNPs | $SNPS |
| Indels | $INDELS |
| Ti/Tv ratio | **$TITV** (healthy: $TITV_OK) |
| PASS records | $PASS_N |
| Annotation | $ANNO_NOTE |

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$STATS")\` | Full \`bcftools stats\` dump |
| \`$(basename "$PASS")\` | PASS-filtered VCF (+ index) |
$( [[ "$FINAL" == "$ANNO" ]] && echo "| \`$(basename "$ANNO")\` | Annotated VCF (+ index) |" )

## Next step
Hand \`$(basename "$FINAL")\` to the **Genome Analyzer** skill for SNPedia / GWAS
cross-referencing. Associations are research-grade, not clinical findings.
EOF

  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh`
Expected: both new assertions pass.

- [ ] **Step 5: Commit**

```bash
git add lib/snp-indel.sh test/test_02_snp_indel.sh
git commit -m "feat: extract snp-indel analyzer into lib/snp-indel.sh"
```

---

## Task 4: analyze.sh dispatcher + type detection

**Files:**
- Create: `analyze.sh`
- Test: `test/test_03_detect.sh`

- [ ] **Step 1: Write the failing test** — `test/test_03_detect.sh`:

```bash
source "$REPO/common.sh"
source "$REPO/analyze.sh" --source-only   # defines functions without running main
for pair in "snp-indel:snp-indel" "sv:sv" "cnv:cnv"; do
  f="${pair%%:*}"; want="${pair##*:}"
  vcf="$(make_fixture "$REPO/test/fixtures/$f.vcf.txt")"
  got="$(detect_type "$vcf")"
  assert_eq "$want" "$got" "detect_type($f) == $want"
done
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — `analyze.sh` not found.

- [ ] **Step 3: Write `analyze.sh`**

```bash
#!/usr/bin/env bash
#
# analyze.sh — detect a single-sample VCF's class (snp-indel | sv | cnv) and
# dispatch to the matching analyzer. See docs/superpowers/specs/.
#
# Usage: ./analyze.sh -i sample.vcf.gz [-o out] [-g genes.bed] [-t N]
#                     [--type snp-indel|sv|cnv] [-a none|snpeff|vep] [-b GRCh38] [-r ref.fa]
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$SCRIPT_DIR/common.sh"

# detect_type <vcf> -> echoes snp-indel|sv|cnv
detect_type() {
  local vcf="$1" hdr recs
  hdr="$(bcftools view -h "$vcf" 2>/dev/null)"
  recs="$(bcftools view -H "$vcf" 2>/dev/null | head -200)"
  if printf '%s' "$hdr" | grep -q '##source=Canvas' \
     || printf '%s' "$hdr" | grep -Eq '##ALT=<ID=CN[0-9]' \
     || printf '%s' "$recs" | grep -q 'SVTYPE=CNV'; then
    printf 'cnv'; return
  fi
  if printf '%s' "$hdr" | grep -Eqi 'Manta|GenerateSVCandidates' \
     || printf '%s' "$hdr" | grep -Eq '##ALT=<ID=(DEL|INS|DUP|INV|BND)' \
     || printf '%s' "$recs" | grep -Eq 'SVTYPE=(DEL|INS|DUP|INV|BND)'; then
    printf 'sv'; return
  fi
  printf 'snp-indel'
}

main() {
  INPUT=""; OUTDIR="out"; GENES_BED=""; THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)"
  FORCE_TYPE=""; ANNOTATE="none"; BUILD="GRCh38"; REF=""
  while getopts ":i:o:g:t:a:b:r:h-:" opt; do
    case "$opt" in
      i) INPUT="$OPTARG" ;;
      o) OUTDIR="$OPTARG" ;;
      g) GENES_BED="$OPTARG" ;;
      t) THREADS="$OPTARG" ;;
      a) ANNOTATE="$OPTARG" ;;
      b) BUILD="$OPTARG" ;;
      r) REF="$OPTARG" ;;
      h) sed -n '2,8p' "$0"; exit 0 ;;
      -) case "$OPTARG" in
           type) FORCE_TYPE="${!OPTIND}"; OPTIND=$((OPTIND+1)) ;;
           type=*) FORCE_TYPE="${OPTARG#*=}" ;;
           *) die "unknown option --$OPTARG" ;;
         esac ;;
      :) die "option -$OPTARG requires an argument" ;;
      \?) die "unknown option -$OPTARG" ;;
    esac
  done

  [[ -n "$INPUT" ]] || { sed -n '2,8p' "$0"; exit 1; }
  [[ -f "$INPUT" ]] || die "input not found: $INPUT"
  have bcftools || die "bcftools not on PATH (install htslib/bcftools first)"
  have tabix    || die "tabix not on PATH (part of htslib)"

  mkdir -p "$OUTDIR"
  BASE="$(basename "${INPUT%.vcf.gz}")"; BASE="${BASE%.vcf}"

  log "validate + index"
  validate_index "$INPUT"

  local TYPE="${FORCE_TYPE:-$(detect_type "$INPUT")}"
  log "detected: $TYPE"
  case "$TYPE" in
    snp-indel) source "$SCRIPT_DIR/lib/snp-indel.sh"; run_snp_indel ;;
    sv)        source "$SCRIPT_DIR/lib/sv.sh";        run_sv ;;
    cnv)       source "$SCRIPT_DIR/lib/cnv.sh";       run_cnv ;;
    *) die "unknown type: $TYPE (use --type snp-indel|sv|cnv)" ;;
  esac
}

# Allow sourcing for tests without running main.
if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh`
Expected: the three `detect_type` assertions pass. (snp-indel end-to-end already covered in Task 3; sv/cnv `run_*` come next.)

- [ ] **Step 5: Make executable and commit**

```bash
chmod +x analyze.sh
git add analyze.sh test/test_03_detect.sh
git commit -m "feat: add analyze.sh dispatcher with type detection"
```

---

## Task 5: analyze-vcf.sh compatibility shim

**Files:**
- Modify: `analyze-vcf.sh` (replace entire contents)
- Test: `test/test_04_shim.sh`

- [ ] **Step 1: Write the failing test** — `test/test_04_shim.sh`:

```bash
out="$("$REPO/analyze-vcf.sh" -i "$(make_fixture "$REPO/test/fixtures/snp-indel.vcf.txt")" -o "$TMP/shim_out" 2>/dev/null)"
[[ -s "$TMP/shim_out/snp-indel.SUMMARY.md" ]] && e=yes || e=no
assert_eq "yes" "$e" "analyze-vcf.sh shim still produces a summary"
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — old `analyze-vcf.sh` writes to a path derived differently / or passes but we want it to route through analyze.sh. (If it passes incidentally, still replace per Step 3 for single-source-of-truth.)

- [ ] **Step 3: Replace `analyze-vcf.sh` with a shim**

```bash
#!/usr/bin/env bash
# analyze-vcf.sh — compatibility shim. The analyzer now lives in analyze.sh,
# which auto-detects snp-indel / sv / cnv. This forwards all arguments.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/analyze.sh" "$@"
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh`
Expected: shim assertion passes.

- [ ] **Step 5: Commit**

```bash
git add analyze-vcf.sh test/test_04_shim.sh
git commit -m "refactor: make analyze-vcf.sh a shim over analyze.sh"
```

---

## Task 6: lib/sv.sh — stats + summary

**Files:**
- Create: `lib/sv.sh`
- Test: `test/test_05_sv.sh`

Gene overlap is wired in Task 9; this task produces stats + summary with a placeholder gene-overlap line.

- [ ] **Step 1: Write the failing test** — `test/test_05_sv.sh`:

```bash
( source "$REPO/common.sh"; source "$REPO/lib/sv.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/sv.vcf.txt")"
  OUTDIR="$TMP/sv_out"; mkdir -p "$OUTDIR"; BASE="sv"; THREADS=1; GENES_BED=""
  run_sv )
sum="$(cat "$TMP/sv_out/sv.SUMMARY.md")"
assert_contains "$sum" "DEL | 2" "SV summary: 2 PASS deletions"
assert_contains "$sum" "DUP | 1" "SV summary: 1 PASS duplication"
assert_contains "$sum" "INS | 1" "SV summary: 1 PASS insertion"
assert_contains "$sum" "BND | 1" "SV summary: 1 PASS breakend event"
assert_contains "$sum" "MinQUAL" "SV summary: lists non-PASS filter reason"
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — `lib/sv.sh` not found.

- [ ] **Step 3: Write `lib/sv.sh`**

```bash
# lib/sv.sh — Manta structural-variant analyzer. Sourced by analyze.sh.
# Expects globals: INPUT OUTDIR BASE GENES_BED ; helpers from common.sh.

run_sv() {
  local PASS="$OUTDIR/${BASE}.pass.vcf.gz"
  local SUMMARY="$OUTDIR/${BASE}.SUMMARY.md"
  local GENES_TSV="$OUTDIR/${BASE}.genes.tsv"

  log "PASS-only filter"
  bcftools view -f PASS,. "$INPUT" -Oz -o "$PASS"
  tabix -f -p vcf "$PASS"

  # Per-record: SVTYPE \t SVLEN \t EVENT \t ID   (PASS only)
  local Q
  Q="$(bcftools query -f '%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/EVENT\t%ID\n' "$PASS")"

  # Counts by type. BND counted as distinct EVENT (fallback: records/2 rounded up).
  local typecounts
  typecounts="$(printf '%s\n' "$Q" | awk -F'\t' '
    $1=="BND"{ ev=($3=="."?$4:$3); if(!(ev in seen)){seen[ev]=1; bnd++} next }
    { c[$1]++ }
    END{
      printf "DEL %d\n", c["DEL"]+0
      printf "DUP %d\n", c["DUP"]+0
      printf "INS %d\n", c["INS"]+0
      printf "INV %d\n", c["INV"]+0
      printf "BND %d\n", bnd+0
    }')"
  local n_del n_dup n_ins n_inv n_bnd
  n_del="$(awk '$1=="DEL"{print $2}' <<<"$typecounts")"
  n_dup="$(awk '$1=="DUP"{print $2}' <<<"$typecounts")"
  n_ins="$(awk '$1=="INS"{print $2}' <<<"$typecounts")"
  n_inv="$(awk '$1=="INV"{print $2}' <<<"$typecounts")"
  n_bnd="$(awk '$1=="BND"{print $2}' <<<"$typecounts")"

  # Size stats from abs(SVLEN), non-BND, non-empty.
  local maxlen medlen
  maxlen="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="BND" && $2!="."{v=($2<0?-$2:$2); if(v>m)m=v} END{print m+0}')"
  medlen="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="BND" && $2!="."{v=($2<0?-$2:$2); print v}' | sort -n | awk '{a[NR]=$1} END{if(NR==0){print 0} else if(NR%2){print a[(NR+1)/2]} else {print int((a[NR/2]+a[NR/2+1])/2)}}')"

  # Per-chromosome tally (PASS) and the largest non-BND events.
  local perchrom largest
  perchrom="$(bcftools query -f '%CHROM\n' "$PASS" | sort -V | uniq -c | awk '{printf "| %s | %s |\n", $2, $1}')"
  largest="$(bcftools query -f '%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/SVLEN\n' "$PASS" \
    | awk -F'\t' '$3!="BND" && $4!="."{v=($4<0?-$4:$4); print v"\t"$1":"$2"\t"$3}' \
    | sort -rn | head -10 | awk -F'\t' '{printf "| %s | %s | %s |\n", $3, $2, $1}')"

  # Non-PASS filter reason tally.
  local reasons
  reasons="$(bcftools view -H "$INPUT" | awk -F'\t' '$7!="PASS" && $7!="."{print $7}' | sort | uniq -c | awk '{printf "%s (%s), ", $2, $1}' | sed 's/, $//')"
  : "${reasons:=none}"

  # Gene overlap placeholder (wired in Task 9).
  local gene_note="_gene overlap added in a later step_"

  log "writing $SUMMARY"
  cat > "$SUMMARY" <<EOF
# SV Summary — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by analyze.sh (sv / Manta)_

## Events by type (PASS)
| Type | Count |
|---|---|
| DEL | $n_del |
| DUP | $n_dup |
| INS | $n_ins |
| INV | $n_inv |
| BND | $n_bnd |

## Size (non-BND, PASS)
| Metric | bp |
|---|---|
| Median | $medlen |
| Largest | $maxlen |

## Per chromosome (PASS)
| Chrom | Events |
|---|---|
$perchrom

## Largest events (non-BND)
| Type | Location | Size (bp) |
|---|---|---|
$largest

## Filtering
Non-PASS records dropped by: $reasons

## Gene overlap
$gene_note

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$PASS")\` | PASS-filtered SV VCF (+ index) |
| \`$(basename "$GENES_TSV")\` | Event → genes-hit table (when enabled) |
EOF

  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh`
Expected: all five SV assertions pass.

- [ ] **Step 5: Commit**

```bash
git add lib/sv.sh test/test_05_sv.sh
git commit -m "feat: add lib/sv.sh Manta SV stats + summary"
```

---

## Task 7: lib/cnv.sh — stats + summary

**Files:**
- Create: `lib/cnv.sh`
- Test: `test/test_06_cnv.sh`

- [ ] **Step 1: Write the failing test** — `test/test_06_cnv.sh`:

```bash
( source "$REPO/common.sh"; source "$REPO/lib/cnv.sh"
  INPUT="$(make_fixture "$REPO/test/fixtures/cnv.vcf.txt")"
  OUTDIR="$TMP/cnv_out"; mkdir -p "$OUTDIR"; BASE="cnv"; THREADS=1; GENES_BED=""
  run_cnv )
sum="$(cat "$TMP/cnv_out/cnv.SUMMARY.md")"
assert_contains "$sum" "Losses | 2" "CNV summary: 2 PASS losses"
assert_contains "$sum" "Gains | 1" "CNV summary: 1 PASS gain"
assert_contains "$sum" "homozygous deletion" "CNV summary: flags CN0 homozygous deletion"
assert_contains "$sum" "L10kb" "CNV summary: notes filtered calls"
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — `lib/cnv.sh` not found.

- [ ] **Step 3: Write `lib/cnv.sh`**

```bash
# lib/cnv.sh — Canvas copy-number-variant analyzer. Sourced by analyze.sh.
# Expects globals: INPUT OUTDIR BASE GENES_BED ; helpers from common.sh.

run_cnv() {
  local PASS="$OUTDIR/${BASE}.pass.vcf.gz"
  local SUMMARY="$OUTDIR/${BASE}.SUMMARY.md"
  local GENES_TSV="$OUTDIR/${BASE}.genes.tsv"

  # Real events only: Canvas:REF segments carry no SVTYPE and ALT='.'; exclude them.
  # Keep records whose ALT is a symbolic <...>. Then PASS-filter.
  log "drop reference segments + PASS filter"
  bcftools view -f PASS,. "$INPUT" \
    | awk -F'\t' '/^#/ || $5 ~ /^</' \
    | bcftools view -Oz -o "$PASS"
  tabix -f -p vcf "$PASS"

  # Per-event: CN \t CNVLEN
  local Q
  Q="$(bcftools query -f '[%CN]\t%INFO/CNVLEN\n' "$PASS")"

  local n_loss n_gain n_homdel loss_bp gain_bp
  n_loss="$(printf '%s\n' "$Q"  | awk -F'\t' '$1!="." && $1<2' | grep -c . || true)"
  n_gain="$(printf '%s\n' "$Q"  | awk -F'\t' '$1!="." && $1>2' | grep -c . || true)"
  n_homdel="$(printf '%s\n' "$Q"| awk -F'\t' '$1=="0"' | grep -c . || true)"
  loss_bp="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="." && $1<2 && $2!="."{s+=$2} END{print s+0}')"
  gain_bp="$(printf '%s\n' "$Q" | awk -F'\t' '$1!="." && $1>2 && $2!="."{s+=$2} END{print s+0}')"
  local loss_mb gain_mb
  loss_mb="$(awk -v b="$loss_bp" 'BEGIN{printf "%.2f", b/1000000}')"
  gain_mb="$(awk -v b="$gain_bp" 'BEGIN{printf "%.2f", b/1000000}')"

  local reasons
  reasons="$(bcftools view -H "$INPUT" | awk -F'\t' '$5 ~ /^</ && $7!="PASS" && $7!="."{print $7}' | sort | uniq -c | awk '{printf "%s (%s), ", $2, $1}' | sed 's/, $//')"
  : "${reasons:=none}"

  # Per-chromosome tally (PASS events) and the largest CNVs by spanned length.
  local perchrom largest
  perchrom="$(bcftools query -f '%CHROM\n' "$PASS" | sort -V | uniq -c | awk '{printf "| %s | %s |\n", $2, $1}')"
  largest="$(bcftools query -f '%CHROM\t%POS\t[%CN]\t%INFO/CNVLEN\n' "$PASS" \
    | awk -F'\t' '$4!="."{dir=($3<2?"loss":"gain"); print $4"\t"$1":"$2"\t"dir" CN="$3}' \
    | sort -rn | head -10 | awk -F'\t' '{printf "| %s | %s | %s bp |\n", $3, $2, $1}')"

  local homdel_note=""
  (( n_homdel > 0 )) && homdel_note=" (incl. $n_homdel homozygous deletion(s), CN=0)"

  local gene_note="_gene overlap added in a later step_"

  log "writing $SUMMARY"
  cat > "$SUMMARY" <<EOF
# CNV Summary — \`${BASE}\`

_Generated $(date -u '+%Y-%m-%d %H:%M UTC') by analyze.sh (cnv / Canvas)_

## Events (PASS, reference segments excluded)
| Class | Count |
|---|---|
| Losses | $n_loss$homdel_note |
| Gains | $n_gain |

## Genome burden (PASS)
| Direction | Mb |
|---|---|
| Lost | $loss_mb |
| Gained | $gain_mb |

## Per chromosome (PASS)
| Chrom | Events |
|---|---|
$perchrom

## Largest events
| Direction | Location | Size (bp) |
|---|---|---|
$largest

## Filtering
Non-PASS events dropped by: $reasons

## Gene overlap
$gene_note

## Artifacts
| File | What it is |
|---|---|
| \`$(basename "$PASS")\` | PASS-filtered CNV VCF, reference segments removed (+ index) |
| \`$(basename "$GENES_TSV")\` | Event → genes-hit table (when enabled) |
EOF

  log "done: $SUMMARY"
  printf '%s\n' "$SUMMARY"
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh`
Expected: all four CNV assertions pass.

- [ ] **Step 5: Commit**

```bash
git add lib/cnv.sh test/test_06_cnv.sh
git commit -m "feat: add lib/cnv.sh Canvas CNV stats + summary"
```

---

## Task 8: gene_overlap helper in common.sh

**Files:**
- Modify: `common.sh` (append `gene_overlap`)
- Test: `test/test_07_overlap.sh`

`gene_overlap` is the shared bedtools wrapper. Tests **skip** when `bedtools` is absent so the suite stays green on this machine; a dedicated degradation test (Task 9) covers the skip path.

- [ ] **Step 1: Write the failing test** — `test/test_07_overlap.sh`:

```bash
if ! command -v bedtools >/dev/null 2>&1; then
  printf '  skip bedtools not installed; gene_overlap tested manually\n'
else
  source "$REPO/common.sh"
  # events BED: chrom \t start \t end \t id  (1-based POS-1 already applied by caller)
  printf '1\t4999\t5500\tD1\n2\t6999\t8500\tU1\n' > "$TMP/ev.bed"
  gene_overlap "$TMP/ev.bed" "$REPO/test/fixtures/genes.bed" "$TMP/ov.tsv"
  out="$(cat "$TMP/ov.tsv")"
  assert_contains "$out" "GENEA" "overlap: chr1 event hits GENEA"
  assert_contains "$out" "GENEB" "overlap: chr2 event hits GENEB"
fi
```

- [ ] **Step 2: Run to verify it fails or skips**

Run: `./test/run-tests.sh`
Expected: prints the skip line (bedtools not installed here). When bedtools is present it FAILs because `gene_overlap` is undefined.

- [ ] **Step 3: Append `gene_overlap` to `common.sh`**

```bash
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
```

- [ ] **Step 4: Run to verify it passes (or still skips)**

Run: `./test/run-tests.sh`
Expected: on this machine, the skip line prints and suite stays green. On a machine with bedtools, the two overlap assertions pass.

- [ ] **Step 5: Commit**

```bash
git add common.sh test/test_07_overlap.sh
git commit -m "feat: add gene_overlap bedtools wrapper to common.sh"
```

---

## Task 9: Wire gene overlap into sv.sh and cnv.sh + degradation

**Files:**
- Modify: `lib/sv.sh` (replace gene-overlap placeholder block)
- Modify: `lib/cnv.sh` (replace gene-overlap placeholder block)
- Test: `test/test_08_degrade.sh`

Resolves the gene-BED location and the graceful-degradation behavior. The default BED is `data/genes.GRCh38.bed` unless `GENES_BED` is set. BND events contribute **both** endpoints.

- [ ] **Step 1: Write the failing degradation test** — `test/test_08_degrade.sh`:

```bash
# With bedtools forced off PATH, SV/CNV must still emit a summary whose
# gene-overlap section shows the skip note rather than failing.
( PATH="$TMP/nobin"; mkdir -p "$TMP/nobin"
  for b in bcftools tabix bgzip awk sort uniq grep sed head paste wc tr cat mktemp dirname basename getconf; do
    src="$(command -v "$b" 2>/dev/null)"; [[ -n "$src" ]] && ln -sf "$src" "$TMP/nobin/$b"
  done
  export PATH
  source "$REPO/common.sh"; source "$REPO/lib/sv.sh"
  INPUT="$(REPO="$REPO" TMP="$TMP" bash -c 'source "'"$REPO"'/test/lib.sh"; make_fixture "'"$REPO"'/test/fixtures/sv.vcf.txt"')"
  OUTDIR="$TMP/deg_out"; mkdir -p "$OUTDIR"; BASE="sv"; GENES_BED=""
  run_sv ) >/dev/null 2>&1
assert_contains "$(cat "$TMP/deg_out/sv.SUMMARY.md")" "bedtools" "degraded SV summary explains how to enable overlap"
```

(If the subshell PATH juggling proves brittle on the target machine, simpler equivalent: temporarily alias `have(){ return 1; }` before `run_sv`. Keep whichever reliably forces the skip path.)

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — current placeholder text is `_gene overlap added in a later step_`, no `bedtools` note.

- [ ] **Step 3: Replace the SV gene-overlap block in `lib/sv.sh`**

Find in `lib/sv.sh`:

```bash
  # Gene overlap placeholder (wired in Task 9).
  local gene_note="_gene overlap added in a later step_"
```

Replace with:

```bash
  # Gene overlap (graceful: stats stand alone if bedtools / gene BED are absent).
  local gene_note bed="${GENES_BED:-$SCRIPT_DIR/data/genes.GRCh38.bed}"
  if ! have bedtools; then
    gene_note="_skipped: bedtools not found. \`brew install bedtools && ./fetch-genes.sh\`_"
  elif [[ ! -s "$bed" ]]; then
    gene_note="_skipped: no gene BED at \`$bed\`. Run \`./fetch-genes.sh\` (or pass \`-g\`)._"
  else
    # Build events BED. Non-BND: POS-1..END. BND: each endpoint as a 1bp interval
    # (this record's POS, plus the mate position parsed from ALT N[chr:pos[ / ]chr:pos]N).
    local ev="$OUTDIR/${BASE}.events.bed"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%ID\t%ALT\n' "$PASS" \
      | awk -F'\t' 'BEGIN{OFS="\t"}
          $4!="BND" && $3!="."{ print $1, $2-1, $3, $5 }
          $4=="BND"{
            print $1, $2-1, $2, $5
            alt=$6
            if(match(alt, /[][][0-9XYMT]+:[0-9]+/)){
              loc=substr(alt, RSTART, RLENGTH); sub(/^[][]/,"",loc)
              split(loc, m, ":"); print m[1], m[2]-1, m[2], $5
            }
          }' > "$ev"
    gene_overlap "$ev" "$bed" "$GENES_TSV"
    gene_note="$(awk -F'\t' 'END{print NR}' "$GENES_TSV") events overlap genes — see \`$(basename "$GENES_TSV")\`."
  fi
```

- [ ] **Step 4: Replace the CNV gene-overlap block in `lib/cnv.sh`**

Find in `lib/cnv.sh`:

```bash
  local gene_note="_gene overlap added in a later step_"
```

Replace with:

```bash
  local gene_note bed="${GENES_BED:-$SCRIPT_DIR/data/genes.GRCh38.bed}"
  if ! have bedtools; then
    gene_note="_skipped: bedtools not found. \`brew install bedtools && ./fetch-genes.sh\`_"
  elif [[ ! -s "$bed" ]]; then
    gene_note="_skipped: no gene BED at \`$bed\`. Run \`./fetch-genes.sh\` (or pass \`-g\`)._"
  else
    local ev="$OUTDIR/${BASE}.events.bed"
    bcftools query -f '%CHROM\t%POS\t%INFO/END\t%ID\n' "$PASS" \
      | awk -F'\t' 'BEGIN{OFS="\t"} $3!="."{print $1, $2-1, $3, $4}' > "$ev"
    gene_overlap "$ev" "$bed" "$GENES_TSV"
    gene_note="$(awk -F'\t' 'END{print NR}' "$GENES_TSV") events overlap genes — see \`$(basename "$GENES_TSV")\`."
  fi
```

Note: `SCRIPT_DIR` is set by `analyze.sh` before sourcing. For standalone test sourcing, add a guard at the top of both `run_sv`/`run_cnv` is unnecessary because tests set `GENES_BED` or rely on the bedtools-absent path. To be safe, both files should define `SCRIPT_DIR` if unset: add near the top of each lib file (outside the function):

```bash
: "${SCRIPT_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
```

- [ ] **Step 5: Run to verify the degradation test passes**

Run: `./test/run-tests.sh`
Expected: degradation assertion passes; all prior assertions still pass.

- [ ] **Step 6: Commit**

```bash
git add lib/sv.sh lib/cnv.sh test/test_08_degrade.sh
git commit -m "feat: wire gene overlap into sv/cnv with graceful degradation"
```

---

## Task 10: fetch-genes.sh

**Files:**
- Create: `fetch-genes.sh`
- Test: `test/test_09_fetch.sh`

Downloads GENCODE basic GRCh38, reduces to one interval per gene symbol, normalizes `chr1`→`1`, writes `data/genes.GRCh38.bed`. The network fetch is **not** run in tests; the test exercises the parse/normalize step via an injected sample.

- [ ] **Step 1: Write the failing test** — `test/test_09_fetch.sh`:

```bash
source "$REPO/fetch-genes.sh" --source-only
# gtf_to_bed reads GENCODE GTF on stdin, writes chrom\tstart\tend\tgene to stdout, chr-normalized.
sample='chr1	HAVANA	gene	1000	2000	.	+	.	gene_name "GENEA"; gene_type "protein_coding";'
got="$(printf '%s\n' "$sample" | gtf_to_bed)"
assert_eq "1	999	2000	GENEA" "$got" "gtf_to_bed: parses + strips chr + 0-bases start"
```

- [ ] **Step 2: Run to verify it fails**

Run: `./test/run-tests.sh`
Expected: FAIL — `fetch-genes.sh` not found.

- [ ] **Step 3: Write `fetch-genes.sh`**

```bash
#!/usr/bin/env bash
# fetch-genes.sh — download GENCODE basic GRCh38 gene spans, collapse to one
# interval per gene symbol, normalize chr-prefix, cache to data/genes.GRCh38.bed.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# GENCODE basic annotation (release 45, GRCh38). Override with GENCODE_URL env var.
GENCODE_URL="${GENCODE_URL:-https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz}"

# gtf_to_bed: GENCODE GTF on stdin -> 'chrom<TAB>start0<TAB>end<TAB>gene_name', chr stripped.
gtf_to_bed() {
  awk -F'\t' 'BEGIN{OFS="\t"}
    $0 ~ /^#/ {next}
    $3=="gene"{
      name=""
      if(match($9, /gene_name "[^"]+"/)){ name=substr($9, RSTART+11, RLENGTH-12) }
      if(name==""){ next }
      chrom=$1; sub(/^chr/, "", chrom)
      print chrom, $4-1, $5, name
    }'
}

# collapse_genes: BED4 on stdin -> one interval per (gene name, chromosome),
# widest span. Keying on name AND chromosome keeps PAR genes (e.g. SHOX on
# chrX and chrY) from merging into one chimeric interval.
collapse_genes() {
  sort -k4,4 -k1,1 -k2,2n \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        { if($4!=g || $1!=gc){ if(g!=""){print c,s,e,g}; g=$4; gc=$1; c=$1; s=$2; e=$3 }
          else { if($2<s)s=$2; if($3>e)e=$3 } }
        END{ if(g!="") print c,s,e,g }' \
    | sort -k1,1 -k2,2n
}

main() {
  local out="$SCRIPT_DIR/data/genes.GRCh38.bed"
  mkdir -p "$SCRIPT_DIR/data"
  command -v bedtools >/dev/null 2>&1 || echo "note: bedtools not installed; analyzers need it for overlap" >&2
  echo "downloading $GENCODE_URL" >&2
  curl -fsSL "$GENCODE_URL" | gunzip -c | gtf_to_bed | collapse_genes > "$out"
  echo "wrote $out ($(wc -l < "$out" | tr -d ' ') genes)" >&2
}

if [[ "${1:-}" != "--source-only" ]]; then
  main "$@"
fi
```

- [ ] **Step 4: Run to verify it passes**

Run: `./test/run-tests.sh`
Expected: the `gtf_to_bed` assertion passes.

- [ ] **Step 5: Make executable and commit**

```bash
chmod +x fetch-genes.sh
git add fetch-genes.sh test/test_09_fetch.sh
git commit -m "feat: add fetch-genes.sh (GENCODE basic GRCh38, chr-normalized)"
```

---

## Task 11: README + docs

**Files:**
- Modify: `README.md` (rewrite for the dispatcher + SV/CNV + install docs)

- [ ] **Step 1: Rewrite `README.md`**

```markdown
# genome-tools

Personal tooling for interpreting germline sequencing deliverables. A single
entry point, `analyze.sh`, reads a `.vcf.gz`, detects whether it holds **SNP/indel**,
**structural-variant (SV)**, or **copy-number-variant (CNV)** calls, and writes a
readable `SUMMARY.md`.

> **Plain-language primer.** A *SNP/indel* call is a single-letter change or a
> small insertion/deletion. A *structural variant* is a larger rearrangement —
> a big deletion, duplication, or a join between two distant places (a breakend).
> A *copy-number variant* is a stretch of chromosome present in too few or too
> many copies. Sequencing.com delivers all three as separate VCFs.

> **Data never goes in git.** See `.gitignore` and the privacy note in `PLAN.md`.

## Install

Required for every run — **bcftools/htslib** (`bcftools`, `tabix`, `bgzip`):

```bash
# macOS (Homebrew). macOS ships bash 3.2; these scripts are written for it.
brew install bcftools htslib
# Debian/Ubuntu
sudo apt install bcftools tabix
# conda (any OS)
conda install -c bioconda bcftools htslib
```

Optional — **bedtools**, needed only for SV/CNV gene overlap:

```bash
brew install bedtools           # macOS
sudo apt install bedtools       # Debian/Ubuntu
conda install -c bioconda bedtools
```

Optional — **SnpEff** or **Ensembl VEP** for SNP/indel annotation (`-a`).

## One-time gene-model setup (SV/CNV gene overlap)

```bash
./fetch-genes.sh     # downloads GENCODE basic GRCh38 -> data/genes.GRCh38.bed
```

Cached under `data/` (git-ignored). Skip it and SV/CNV runs still produce full
stats — only the "genes hit" table is omitted, with a note on how to enable it.
Use a different gene set any time with `-g /path/to/genes.bed`.

## Usage

```bash
./analyze.sh -i sample.vcf.gz -o out/        # auto-detects type
```

It prints e.g. `detected: cnv` and writes results to `out/`.

| Flag | Meaning |
|---|---|
| `-i FILE` | input `.vcf.gz` (required) |
| `-o DIR` | output directory (default `out`) |
| `-g FILE` | gene BED for overlap (default `data/genes.GRCh38.bed`) |
| `-t N` | threads (annotation) |
| `--type T` | force `snp-indel`\|`sv`\|`cnv` instead of auto-detect |
| `-a ENGINE` | SNP/indel annotation: `none`\|`snpeff`\|`vep` |
| `-b BUILD` | genome build for annotation (default `GRCh38`) |
| `-h` | help |

Examples:

```bash
./analyze.sh -i her.snp-indel.genome.vcf.gz -o out/          # SNP/indel summary
./analyze.sh -i her.sv.vcf.gz  -o out/                       # SV summary + gene overlap
./analyze.sh -i her.cnv.vcf.gz -o out/                       # CNV summary + gene overlap
./analyze.sh -i her.snp-indel.genome.vcf.gz -a snpeff        # + SnpEff annotation
```

`analyze-vcf.sh` still works — it forwards to `analyze.sh`.

## Outputs (in `out/`)

| File | What it is |
|---|---|
| `*.SUMMARY.md` | Human-readable summary |
| `*.pass.vcf.gz` (+index) | PASS-filtered calls |
| `*.genes.tsv` | Event → genes-hit table (SV/CNV, when overlap enabled) |
| `*.stats.txt` | Raw `bcftools stats` dump (SNP/indel) |

## Testing

```bash
./test/run-tests.sh
```

Runs on synthetic fixtures; no private data involved. Gene-overlap tests skip
automatically when `bedtools` is not installed.

## Roadmap
See `PLAN.md`.

## License
Licensed under the [Apache License, Version 2.0](LICENSE).
```

- [ ] **Step 2: Sanity-check the docs match reality**

Run: `./analyze.sh -h`
Expected: usage lines print without error (flags match the table above).

- [ ] **Step 3: Commit**

```bash
git add README.md
git commit -m "docs: rewrite README for analyze.sh dispatcher + SV/CNV + install"
```

---

## Task 12: Full integration on real data (manual)

**Files:** none (verification only)

- [ ] **Step 1: Run the whole suite**

Run: `./test/run-tests.sh`
Expected: `N passed, 0 failed`.

- [ ] **Step 2: Run all three real sample VCFs**

```bash
D=~/Dropbox/My\ Documents/sample
./analyze.sh -i "$D"/*.snp-indel.genome.vcf.gz -o ./work/sample
./analyze.sh -i "$D"/*.sv.vcf.gz  -o ./work/sample
./analyze.sh -i "$D"/*.cnv.vcf.gz -o ./work/sample
```

Expected: each prints the right `detected:` line and writes a `*.SUMMARY.md`.
Open the SV and CNV summaries; confirm counts look sane (SV dominated by DEL;
CNV burden in Mb is plausible). If `bedtools` + gene BED are installed, confirm
the `*.genes.tsv` tables are populated; otherwise confirm the skip note appears.

- [ ] **Step 3: (optional) enable gene overlap and re-run**

```bash
brew install bedtools
./fetch-genes.sh
./analyze.sh -i ~/Dropbox/My\ Documents/sample/*.cnv.vcf.gz -o ./work/sample
```

Expected: CNV summary's gene-overlap section now points to a populated `*.genes.tsv`.

---

## Notes / decisions captured

- All scripts target bash 3.2 (the original line-28 `printf '%(...)T'` bug is why).
- Fixtures are committed as `.vcf.txt` and bgzipped at test time, so no genomic
  data and no `.vcf.gz` (already gitignored) enters the repo.
- BND gene overlap reports genes at both breakend endpoints (per spec decision).
- Gene BED is normalized to bare contig names at fetch time (per spec decision).
- `data/` is already in `.gitignore`, so the cached gene BED stays out of git.
