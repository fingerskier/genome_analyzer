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
| `-r FILE` | reference FASTA (lets VEP add HGVS notation) |
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
