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
