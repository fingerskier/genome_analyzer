# Trait Cross-Reference — Design

_2026-06-04_

## Problem

`analyze.sh` summarizes a SNP/indel VCF but stops at descriptive stats. The
SUMMARY.md even ends with a "hand this to the Genome Analyzer skill for SNPedia /
GWAS cross-referencing" pointer — but no such tool exists yet. We want to take a
SNP/indel PASS VCF and report **which catalogued traits and clinically-flagged
variants the sample actually carries alleles for**, including zygosity (1 or 2
copies of the relevant allele), not merely "this SNP is associated with X."

This is research-grade exploration, not clinical diagnosis.

## Goals

- Cross-reference a SNP/indel PASS VCF against two public databases:
  - **GWAS Catalog** (EBI) — rsID → trait / risk allele / effect size.
  - **ClinVar** (NCBI) — variant → clinical significance / condition.
- For every hit, compute **risk-/reported-allele dosage** (0/1/2 copies) from the
  sample genotype — the substance of the report.
- Produce a readable `*.traits.SUMMARY.md` plus full machine-readable TSVs.
- One-time downloader scripts mirroring `fetch-genes.sh`, caching to `data/`.
- Degrade gracefully (clear instructions) when a database file is absent.
- Honest handling of strand-ambiguous (palindromic) SNPs — flag, don't guess.

## Non-goals (v1)

- Polygenic risk scores / aggregate genetic burden.
- Haplotype or multi-SNP associations (only single-rsID catalog rows).
- Population/ancestry allele-frequency weighting.
- Clinical diagnosis or medical advice. Output is research-grade.
- Re-deriving anything from FASTQ/BAM.

## Resolved decisions (from brainstorming)

- **Databases:** GWAS Catalog **and** ClinVar. ClinVar crosses the parent
  project's "no clinical interpretation" non-goal deliberately — the user opted
  in — but we still present it as research-grade, never diagnostic.
- **GWAS significance filter:** genome-wide significant only, **p < 5e-8**.
  Documented as the default with a future-work note to expose a threshold flag
  for expansion (lower p-value, include VUS).

## Architecture

```
xref-traits.sh        # entry point: arg parse, validate input, dispatch to engine
lib/trait-xref.sh     # the engine: extract genotypes, join both DBs, dosage, emit
fetch-gwas.sh         # one-time: download GWAS Catalog -> data/gwas-catalog.tsv
fetch-clinvar.sh      # one-time: download clinvar.vcf.gz (+index) -> data/
common.sh             # reuse existing log/die/have/validate_index helpers
data/                 # cached DB files (git-ignored)
```

`xref-traits.sh` is a **separate entry point** from `analyze.sh`, not a new
auto-detected route: it is post-processing on an already-summarized SNP/indel
PASS VCF, a different kind of operation than `analyze.sh`'s class summaries.
All scripts target **bash 3.2** (macOS default), consistent with the rest of the
repo — no associative arrays, no `mapfile`, no `printf '%(...)T'`.

### Dependencies

- Required: `bcftools`, `tabix`/`bgzip` (htslib) — already used everywhere.
- Required tools to fetch: `curl`, `unzip` (GWAS zip), `gzip`.
- **No bedtools** — joins are by rsID or by position+allele via `awk`, not
  interval overlap.

## `fetch-gwas.sh`

Downloads the GWAS Catalog "ontology-annotated, alternative, full" associations
file and normalizes it to a compact, join-ready table.

- **Source:**
  `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip`
  (~68 MB zip → ~710 MB TSV named `gwas-catalog-download-associations-alt-full.tsv`).
- **Relevant columns** (1-based, verified against the live header):
  | Col | Header | Use |
  |---|---|---|
  | 2 | `PUBMEDID` | study reference |
  | 8 | `DISEASE/TRAIT` | trait label |
  | 15 | `MAPPED_GENE` | gene context |
  | 21 | `STRONGEST SNP-RISK ALLELE` | `rs1234-A` → rsid + risk allele |
  | 22 | `SNPS` | rsid (fallback / validation) |
  | 28 | `P-VALUE` | significance filter |
  | 31 | `OR or BETA` | effect size |
- **Normalization** (`awk`, skipping the header):
  - Parse `STRONGEST SNP-RISK ALLELE` on the last `-`: left = rsid, right = risk
    allele. Keep only rows whose **rsid matches `^rs[0-9]+$`** and whose **risk
    allele is a single base in `{A,C,G,T}`**. This drops multi-SNP haplotype rows
    (`rs1; rs2`), interaction rows, and `rs1234-?` unknown-allele rows — exactly
    the v1 non-goals.
  - Emit one row per association:
    `rsid \t risk_allele \t trait \t pval \t or_beta \t gene \t pmid`
    to `data/gwas-catalog.tsv`. **No p-value pre-filtering at fetch** — the full
    normalized table is cached so a future threshold change needs no re-download;
    the engine applies p < 5e-8 at join time.
  - Sort by rsid so the engine can join efficiently.
- Idempotent: skip download if `data/gwas-catalog.tsv` already exists unless
  `-f`/`--force` is passed. Log final row count.

## `fetch-clinvar.sh`

- **Source:** `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz`
  (+ `.tbi`). Contigs are **bare names** (`1`, `2`, … `X`, `Y`, `MT`) — already
  matching the sample VCF, so no contig reconciliation needed.
- Download the `.vcf.gz`; fetch the published `.tbi` if available, else build it
  with `tabix -p vcf`. Cache both to `data/clinvar.GRCh38.vcf.gz(.tbi)`.
- Relevant INFO fields (verified present): `CLNSIG`, `CLNSIGCONF`, `CLNDN`,
  `CLNREVSTAT`, `GENEINFO`.
- Idempotent with `-f`/`--force`. Log record count.

## `lib/trait-xref.sh` — engine (`run_trait_xref`)

Input: a SNP/indel PASS `.vcf.gz` (typically `*.snp-indel.genome.pass.vcf.gz`).

1. **Normalize genotypes.** `bcftools norm -m- <in>` to split multiallelic sites
   so every record is biallelic and each `GT` is `0/0`, `0/1`, or `1/1`. Stream
   the result; do not write a second large VCF unless convenient.
2. **Extract the sample's variant calls** with
   `bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'`, keeping only
   records where `GT` contains a `1` (non-reference). From each record derive the
   sample's **two alleles**: `0/1` → {REF, ALT}; `1/1` → {ALT, ALT}. (Phased `|`
   treated same as `/`.) Build two streams:
   - **rsID stream** for GWAS: `rsid \t allele1 \t allele2` (rsID = `%ID` when it
     matches `^rs`).
   - **position stream** for ClinVar: `chrom \t pos \t ref \t alt \t zygosity`
     (`het` for 0/1, `hom` for 1/1).
3. **GWAS join (by rsID).** Join the rsID stream against `data/gwas-catalog.tsv`.
   For each association with **p < 5e-8** at an rsID the sample carries:
   - Compute **dosage** = count of the sample's two alleles equal to
     `risk_allele`.
   - **Strand handling:** if `risk_allele` is not among the sample's REF/ALT
     alleles, try its complement. If the complement matches *and* the SNP is not
     palindromic (not A/T, not C/G), count the dosage and tag `strand_flipped`.
     If the SNP **is** palindromic (A/T or C/G), set dosage to `NA` and tag
     `ambiguous` — never guess strand on a palindrome.
   - Keep rows with `dosage >= 1` **or** `dosage = NA` (carried, or unresolved).
     Output `rsid, trait, your_genotype, risk_allele, dosage, pval, or_beta,
     gene, pmid, flag` to `*.traits.gwas.tsv`.
4. **ClinVar join (by position + allele).** Pull pathogenic-class records once:
   `bcftools query -i 'CLNSIG ~ "Pathogenic" || CLNSIG ~ "Likely_pathogenic" ||
   CLNSIG ~ "risk_factor" || CLNSIG ~ "drug_response"' -f
   '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNDN\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\n'`
   from `data/clinvar.GRCh38.vcf.gz`, **excluding** rows whose CLNSIG is purely
   `Benign`/`Likely_benign`. Join against the sample's position stream on
   `chrom:pos:ref:alt` (exact match — `bcftools norm` already left-aligned and
   split the sample side; ClinVar is already normalized). Output matches with
   the sample's `zygosity`, the condition (`CLNDN`, underscores → spaces), review
   status (`CLNREVSTAT`, surfaced as a star/"conflicting" note via `CLNSIGCONF`),
   and gene to `*.traits.clinvar.tsv`.
5. **Emit `*.traits.SUMMARY.md`** (see Outputs).

### Genotype semantics note

The input is a genome VCF with reference blocks, so a site the sample is
explicitly called `0/0` at is genuinely homozygous reference — real "0 copies"
signal, not missing data. v1 reports only carried alleles (dosage ≥ 1) plus
ambiguous, but this fact is noted in the summary so the absence of a known risk
variant can be read as informative.

## Outputs (per run, in the chosen `-o` dir; all git-ignored)

| Artifact | Description |
|---|---|
| `*.traits.SUMMARY.md` | Readable report: ClinVar pathogenic hits first (condition + zygosity + review status), then GWAS trait hits the sample carries ≥1 risk allele for, grouped by trait, **homozygous hits highlighted**, ambiguous-strand SNPs listed separately. Counts at top. Caveats footer. |
| `*.traits.gwas.tsv` | Full GWAS hit table (one row per association). |
| `*.traits.clinvar.tsv` | Full ClinVar pathogenic-class hit table. |

`.gitignore` adds `*.traits.tsv`, `*.traits.gwas.tsv`, `*.traits.clinvar.tsv`
(genomic-derived personal data) alongside the existing `*.SUMMARY.md` rule.

## Degradation (graceful)

- If `data/gwas-catalog.tsv` is absent: skip the GWAS section, summary notes
  `[gwas] skipped: run ./fetch-gwas.sh`.
- If `data/clinvar.GRCh38.vcf.gz` is absent: skip the ClinVar section with the
  analogous `./fetch-clinvar.sh` note.
- If both absent: the run still validates the input and writes a SUMMARY.md whose
  body is the two skip notes (so the user always gets actionable guidance).
- `-g`/`-c` flags allow pointing at custom GWAS table / ClinVar VCF paths.

## CLI (`xref-traits.sh`)

```
./xref-traits.sh -i sample.snp-indel.genome.pass.vcf.gz -o work
```

| Flag | Meaning |
|---|---|
| `-i FILE` | input SNP/indel PASS `.vcf.gz` (required) |
| `-o DIR` | output directory (default `out`) |
| `-g FILE` | GWAS table (default `data/gwas-catalog.tsv`) |
| `-c FILE` | ClinVar VCF (default `data/clinvar.GRCh38.vcf.gz`) |
| `-h` | help |

## Testing

Zero-dependency bash harness, consistent with the existing `test/` suite.
Synthetic fixtures only — no private data, fake rsIDs/conditions.

- **fetch-gwas normalization:** feed a tiny fixture mimicking the GWAS TSV header
  + rows (one clean `rs1-A`, one multi-SNP `rs2; rs3`, one `rs4-?`, one
  non-rs); assert only the clean row survives with correct columns.
- **dosage:** sample `1/1` for a `risk=ALT` SNP → dosage 2; `0/1` → 1; risk =
  REF with `0/1` → 1; risk = REF with `1/1` → 0 (dropped).
- **strand:** non-palindromic risk allele matching the complement → counted +
  `strand_flipped`; palindromic A/T risk allele → `NA` + `ambiguous`.
- **p-value filter:** association at p = 1e-6 excluded; p = 1e-9 kept.
- **ClinVar join:** sample variant matching a fixture Pathogenic record at
  chrom:pos:ref:alt → reported with correct zygosity; a Benign record at a
  carried site → not reported.
- **degradation:** missing GWAS table and/or ClinVar VCF → SUMMARY.md written
  with the correct skip note(s), exit 0.
- **summary smoke:** end-to-end on fixtures produces a non-empty SUMMARY.md with
  the expected section headers and a homozygous-hit highlight.

## README

Add a "Trait cross-reference" section: what it is (research-grade, not clinical),
the two one-time `fetch-gwas.sh` / `fetch-clinvar.sh` downloads (sizes, that they
cache to git-ignored `data/`), the `xref-traits.sh` usage and flag table, the
outputs, and a prominent **not-a-diagnosis** caveat plus the palindrome/strand
limitation in plain language.
