# SV & CNV Analyzers — Design

_2026-06-04_

## Problem

`analyze-vcf.sh` summarizes single-sample germline **SNP/indel** VCFs. A typical
Sequencing.com 30x WGS delivery also includes two other VCFs that the current
script cannot meaningfully interpret:

- `*.sv.vcf.gz` — structural variants from **Manta 1.6.0** (DEL/INS/DUP/BND,
  symbolic ALTs, `SVTYPE`/`SVLEN`/`END`).
- `*.cnv.vcf.gz` — copy-number variants from **Canvas 1.40.0** (`<CN0>`/`<CN2>`…
  symbolic ALTs, `CN` copy-number genotype, `SVTYPE=CNV`, `CNVLEN`).

The raw `*.fq.gz` paired FASTQ reads are out of scope — they are unaligned reads,
not variant calls.

We want analyzers that produce the same kind of tidy `SUMMARY.md` for these two
variant classes, including gene-level overlap of the events.

## Goals

- One entry point that auto-detects VCF class (SNP/indel | SV | CNV) and routes.
- Descriptive stats appropriate to each class.
- Gene overlap for SV and CNV events (which genes each event hits).
- No regression in existing SNP/indel behavior.
- Degrade gracefully when optional tooling (bedtools / gene BED) is absent.
- A README that documents installing every required and optional tool, plus
  usage for each route — written for someone new to genetics tooling.

## Non-goals

- Re-calling variants from FASTQ.
- Clinical interpretation / pathogenicity classification.
- Dosage-sensitivity or ClinVar/known-disease-gene scoring (future work).
- Multi-sample / somatic logic.

## Architecture

```
analyze.sh            # entry point: arg parse, validate+index, sniff type, dispatch
common.sh             # sourced helpers: log/die/have, validate+index, gene_overlap(), summary scaffold
lib/snp-indel.sh      # current analyze-vcf.sh logic, moved ~verbatim
lib/sv.sh             # Manta SV analyzer
lib/cnv.sh            # Canvas CNV analyzer
fetch-genes.sh        # one-time downloader: GENCODE basic GRCh38 -> data/genes.GRCh38.bed
data/                 # cached gene BED (git-ignored)
analyze-vcf.sh        # compatibility shim -> exec analyze.sh "$@"
```

All scripts target **bash 3.2** (macOS default) — no `printf '%(...)T'`, no
associative arrays, no `mapfile`. (This is the bug already fixed in the current
script: line 28 used a bash-4.2 time format that broke on macOS bash 3.2.57.)

### Type detection (`analyze.sh`)

Read the header + a sample of records once, decide in this order:

1. **CNV** if header has `##source=Canvas` OR any `##ALT=<ID=CN%` OR records carry
   `SVTYPE=CNV`.
2. **SV** if header has `Manta`/`GenerateSVCandidates` OR `##ALT=<ID=DEL/INS/DUP/INV/BND`
   OR records carry `SVTYPE=` (DEL/INS/DUP/INV/BND) without CNV.
3. **SNP/indel** otherwise.

Detection is logged (`detected: Canvas CNV`) so the user can confirm the route.
A `--type snp-indel|sv|cnv` override is available to force a route.

### Shared flow (`common.sh`)

- `log`/`die`/`have` — same helpers as today.
- `validate_index <vcf>` — the current stage-1 bgzip/readability check + tabix.
- `gene_overlap <pass.bed> <out.tsv>` — `bedtools intersect` of events vs the gene
  BED, producing `event_id \t chrom \t start \t end \t n_genes \t gene1,gene2,...`,
  capping the gene list (e.g. first 25 + `(+N more)`) for very large events.
- `summary_header`/`summary_footer` — shared markdown scaffold so all three
  SUMMARY.md files look consistent.

## SV analyzer (`lib/sv.sh`)

1. Validate + index (shared).
2. PASS filter: keep `FILTER ∈ {PASS, .}` → `*.sv.pass.vcf.gz`.
3. Stats (bcftools query + awk):
   - Counts by `SVTYPE`: DEL / INS / DUP / INV / BND. BND counted as paired
     **events** (dedupe via `MATEID`/`EVENT`), not raw breakend records.
   - Size distribution from `abs(SVLEN)`: median, max, and buckets
     (<1kb, 1–10kb, 10–100kb, 100kb–1Mb, >1Mb). BNDs excluded from size stats
     (no meaningful length).
   - Per-chromosome tally.
   - Largest N events (default 10).
   - Tally of non-PASS FILTER reasons (`MinQUAL`, `NoPairSupport`, `MaxDepth`, …).
4. Gene overlap on PASS events → `*.sv.genes.tsv`. Non-BND events use their
   `POS..END` span. **BND** events report genes at *both* breakend endpoints
   (POS on this record + the mate position parsed from the ALT `N[chr:pos[`
   notation), since the gene at each side is the point of a translocation.
5. Emit `*.sv.SUMMARY.md`.

## CNV analyzer (`lib/cnv.sh`)

1. Validate + index (shared).
2. Drop reference segments: exclude `Canvas:REF` records (ALT `.`, CN=2) — these
   are not events. Keep `GAIN`/`LOSS` with symbolic `<CN%>` ALT.
3. PASS filter on the remaining events.
4. Stats (bcftools query of `INFO/CNVLEN` + `FORMAT/CN`):
   - Losses (CN<2; CN0 = homozygous deletion, flagged) vs gains (CN>2).
   - Size distribution from `CNVLEN`.
   - **Genome burden**: total Mb gained and lost.
   - Per-chromosome tally; largest events.
   - Note on `L10kb` / `q7` / `FailedFT` filtered calls.
5. Gene overlap on PASS events → `*.cnv.genes.tsv`.
6. Emit `*.cnv.SUMMARY.md`.

## Gene model (`fetch-genes.sh`)

- Default: **GENCODE basic** annotation for GRCh38 (broad coverage incl.
  protein-coding + lncRNA), collapsed to one interval per gene symbol.
  Rationale: best fit for a personal-WGS exploration workflow where we want to
  catch any gene an event touches; switch with `-g` if a curated set is wanted.
- Output: `data/genes.GRCh38.bed` (`chrom \t start \t end \t gene_symbol`),
  chromosome naming **normalized to bare contig names** (`chr1` → `1`) on fetch
  so it matches the Manta/Canvas VCFs without per-run reconciliation.
- Override: `analyze.sh -g /path/to/custom.bed`.

## Degradation (graceful)

If `bedtools` is not on PATH, or the gene BED is absent and not provided:
the SV/CNV run completes the **stats** portion and writes SUMMARY.md, with the
gene-overlap section replaced by:

```
[gene overlap] skipped: bedtools not found.
  brew install bedtools && ./fetch-genes.sh
```

The SNP/indel route never requires bedtools.

## Dependencies

- Required (all routes): `bcftools`, `tabix` (htslib) — already used.
- Optional (SV/CNV gene overlap): `bedtools`, `data/genes.GRCh38.bed`.

## README (`README.md`)

Written for someone new to genetics tooling. Covers:

- **What this is** — one paragraph: turns Sequencing.com WGS VCFs into readable
  summaries; brief plain-language note on what SNP/indel vs SV vs CNV mean.
- **Install** — macOS (Homebrew) and Linux (apt/conda) commands for the required
  `bcftools`/`htslib` and the optional `bedtools`; note that macOS bash 3.2 is
  fine (scripts are written for it).
- **One-time gene model setup** — `./fetch-genes.sh` (what it downloads, where it
  caches, that it's only needed for SV/CNV gene overlap).
- **Usage** — `./analyze.sh -i file.vcf.gz -o out/` auto-detects type; the full
  flag list (`-o`, `-g`, `-t`, `--type`, `-h`); one example per route.
- **Outputs** — what each artifact file is.
- **Caveat** — associations are research-grade, not clinical.

## Outputs (per run)

| Artifact | Description |
|---|---|
| `*.{sv,cnv}.SUMMARY.md` | Human-readable summary, same style as SNP/indel |
| `*.{sv,cnv}.pass.vcf.gz` (+index) | PASS-filtered events |
| `*.{sv,cnv}.genes.tsv` | Event → genes-hit table (when overlap enabled) |
| `*.{sv,cnv}.stats.txt` | Raw bcftools stats dump |

## Testing

- Use the real sample SV (490K) and CNV (57K) files as fixtures (small, fast).
- Per analyzer: assert type detection routes correctly; assert counts are
  internally consistent (sum of SVTYPE counts == PASS record count; gains+losses
  == CNV event count); assert SUMMARY.md is generated and non-empty.
- Degradation test: run with `bedtools` masked off PATH → stats still produced,
  overlap section shows the skip note.
- Regression: existing SNP/indel file routes to snp-indel and produces the same
  SUMMARY as today.

## Resolved decisions

- **BND gene overlap** → report genes at *both* breakend endpoints (the gene on
  each side is the substance of a translocation). Folded into the SV analyzer.
- **`chr`-prefix** → normalize the gene BED to bare contig names (`chr1` → `1`)
  at fetch time, matching the Manta/Canvas VCFs. Folded into `fetch-genes.sh`.
