# Population-Frequency Context — Design

_2026-08-09_

## Problem

The trait cross-reference (`xref-traits.sh`) reports which ClinVar pathogenic-class
and GWAS trait alleles the sample carries, but gives no sense of **how common each
variant is in the general population**. Frequency is the single best triage signal:
a "pathogenic"-labelled variant carried by a third of the population is almost
certainly low-impact, while a variant absent from population databases deserves
more attention. PLAN.md Phase 2 calls this out as
"flag common variants (almost all benign) for triage".

## Goals

- Every ClinVar hit row gains a **population allele frequency** and a
  plain-language **frequency bucket** (common / low-frequency / rare / unknown).
- Every GWAS hit row gains the catalog's **risk-allele frequency**.
- The SUMMARY.md explains, in plain language, why a common "pathogenic" variant
  is usually not alarming.
- **Zero new downloads, fully offline.** No new fetch scripts, no new `data/` files.
- Degrade gracefully against an old (7-column) cached GWAS table.

## Non-goals

- Genome-wide frequency annotation of the full VCF (the SnpEff/VEP `-a` path).
- gnomAD as a data source (see decision below).
- Ancestry-specific frequency weighting; the numbers are aggregate cohort AFs.
- Re-scoring or filtering hits by frequency — frequency is **context, not a filter**.
  Every row that appears today still appears.

## Resolved decisions

- **Data sources: harvest what we already cache — no new database.** Verified
  sizes (2026-08-09) rule out local copies of the real frequency databases on
  this machine (72 GB free): gnomAD v4.1 joint sites is ~72 GB *for chr1 alone*,
  full dbSNP is 29.5 GB, dbSNP ALFA is 19.7 GB. Instead:
  - **ClinVar hits:** the cached `data/clinvar.GRCh38.vcf.gz` already carries
    `AF_EXAC`, `AF_TGP`, `AF_ESP` INFO fields (23% / 10% / 10% of records).
    Coverage skews toward exactly the variants we need to flag — common variants
    are the ones population cohorts observed. A missing AF is itself weak
    evidence of rarity, and is reported honestly as `unknown`.
  - **GWAS hits:** the GWAS Catalog raw file has a `RISK ALLELE FREQUENCY`
    column (col 27; the layout is already verified through col 31) that
    `normalize_gwas` currently drops. Keep it. Requires a one-time
    `./fetch-gwas.sh -f` re-run (~68 MB, same download as before).
- **No online per-variant API lookups (privacy).** Batch-querying Ensembl/gnomAD
  APIs for the sample's carried rsIDs would transmit which variants this genome
  contains to a third party. This repo treats the data as maximally sensitive;
  frequency context is not worth a genotype side-channel. If richer coverage is
  ever wanted, the offline ALFA download (19.7 GB) is the documented future-work
  path — never an API.
- **Frequency buckets:** `common` ≥ 5%, `low-frequency` 1–5%, `rare` < 1%,
  `unknown` otherwise. Standard clinical-genetics convention (BA1/BS1-adjacent),
  chosen for readability, not classification.
- **Append-only columns.** New TSV columns go at the end of each row so existing
  column positions (and any downstream awk) stay valid.

## Architecture

No new files. Three existing pieces change:

```
fetch-gwas.sh        # normalize_gwas keeps RAF -> 8th column of data/gwas-catalog.tsv
lib/trait-xref.sh    # engine: pass RAF through; extract ClinVar AF_*; buckets in summary
test/                # extended fixtures + new freq_bucket unit test
```

### `fetch-gwas.sh` — `normalize_gwas`

- Emit an 8th column `raf` from raw col 27 (`RISK ALLELE FREQUENCY`).
- **Header sanity check:** on the header line, if col 27 is not
  `RISK ALLELE FREQUENCY`, print a warning to stderr (the layout drifted) and
  emit `NA` for every raf rather than harvesting a wrong column. Cols 28/31 were
  live-verified in the trait-xref work; col 27 is between them but gets checked
  at run time, and the first real `-f` re-fetch during implementation confirms it.
- Validation: keep the value only if it is a plain number in [0, 1]
  (`^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$` and `<= 1`); otherwise emit `NA`.
  The raw column contains `NR`, ranges, and annotated values like `0.23 (EA)` —
  all become `NA` rather than guessed at.
- New cached table shape: `rsid  risk  trait  pval  orbeta  gene  pmid  raf`.

### `lib/trait-xref.sh` — engine

- **`freq_bucket <af>`** — new small function beside `risk_dosage`, unit-testable:
  maps a frequency to `common`/`low-frequency`/`rare`/`unknown` (junk or `NA` →
  `unknown`). Single source of truth for the thresholds, implemented like
  `DOSAGE_FN` (awk fragment shared by the unit wrapper and the summary pass).
- **GWAS join:** pack `$8` (raf) into the per-rsID record; append it to each
  output row. `*.traits.gwas.tsv` becomes 11 columns:
  `rsid trait genotype risk_allele copies pval or_beta gene pmid flag raf`.
  - **Legacy-table degradation:** if the cached GWAS table has only 7 columns
    (pre-this-feature cache), every raf is `NA` and the GWAS section note in the
    summary appends: `frequencies unavailable — refresh with ./fetch-gwas.sh -f`.
    Detected once from the table's first line, not per-row.
- **ClinVar join:** extend the `bcftools query -f` format with
  `%INFO/AF_EXAC`, `%INFO/AF_TGP`, `%INFO/AF_ESP`. Per hit, `pop_af` = first
  non-`.` of EXAC → TGP → ESP (largest cohort first) and `af_source` names it
  (`ExAC`/`1000G`/`ESP`); both `NA` when absent. `*.traits.clinvar.tsv` becomes
  9 columns: `chrom:pos alt zygosity significance condition review gene pop_af af_source`.

### SUMMARY.md changes

- **ClinVar table** gains a `How common?` column: `common (~30%)`,
  `low-frequency (~2%)`, `rare (<1%)`, or `unknown` (bucket + rounded percent).
- **GWAS table** gains a `Risk-allele freq` column (`23%` or `—` for NA).
- **"How to read this"** gains two bullets:
  - A variant that is *common* is carried by a large fraction of the population;
    despite alarming condition names, common variants are almost never seriously
    harmful on their own (natural selection would have removed them). Rarity is
    a reason to look closer, not a verdict.
  - `unknown` means the public cohorts behind ClinVar's frequency fields
    (ExAC / 1000 Genomes / ESP) didn't report it — often a sign of rarity, but
    also common for indels and recently-catalogued variants.
- Artifacts table: update both TSV column listings.

### Unchanged

`xref-traits.sh` flags and CLI, `fetch-clinvar.sh`, `analyze.sh` and all
`lib/{snp-indel,sv,cnv}.sh`, degradation behavior for missing databases.

## Testing

Extend the existing zero-dependency harness; synthetic fixtures only.

- **`test_13_gwas_normalize.sh`:** extend the fixture header/rows through col 27+;
  assert a numeric RAF (`0.23`) survives as col 8, and `NR` / `0.23 (EA)` /
  empty → `NA`.
- **New `freq_bucket` unit test:** `0.4` → common, `0.02` → low-frequency,
  `0.001` → rare, `NA`/`abc`/empty → unknown; boundary values `0.05` → common,
  `0.01` → low-frequency.
- **`test_15_xref.sh` (engine e2e):**
  - GWAS fixture table gains an 8th column; assert raf lands in GTSV col 11 and
    the summary GWAS table shows a percent.
  - ClinVar fixture VCF gains `AF_EXAC` on one record (header + INFO); assert
    CTSV pop_af/af_source and a `common` bucket in the summary for it, and
    `unknown` for a record without AF fields.
  - **Legacy 7-column GWAS table** → all raf `NA`, summary carries the
    `fetch-gwas.sh -f` refresh note, exit 0.
- All existing tests must stay green (column additions are append-only).

## README

In the trait cross-reference section: mention the frequency columns, the
buckets, that ClinVar frequencies come from ExAC/1000G/ESP via the cached
ClinVar VCF, and that upgrading an existing GWAS cache needs one
`./fetch-gwas.sh -f`. Plain-language framing: frequency is context for triage,
not a verdict.

## PLAN.md

Tick the Phase 2 gnomAD line, rewording to what was actually built (population
frequency via harvested sources; gnomAD-scale local databases documented as
infeasible on this hardware, ALFA noted as the offline future-work option).
