# genome-tools

Personal tooling for interpreting germline sequencing deliverables (`.vcf.gz` + `.fq.gz`).
Start from the called variants (fast), keep the FASTQs for optional re-analysis.

> **Data never goes in git.** See `.gitignore` and the privacy note in `PLAN.md`.

## Requirements
Core path needs only **htslib/bcftools** (`bcftools`, `tabix`, `bgzip`).
Optional annotation: **SnpEff** or **Ensembl VEP**. FASTQ checks: **seqkit**.

```bash
# macOS
brew install bcftools htslib seqkit
# Debian/Ubuntu
sudo apt install bcftools tabix
```

## Quick start
```bash
# core: stats + PASS filter + summary  (no heavy deps)
./scripts/analyze-vcf.sh -i data/her.vcf.gz -o out

# with SnpEff annotation
./scripts/analyze-vcf.sh -i data/her.vcf.gz -o out -a snpeff -b GRCh38

# multiple per-chromosome shards? concat first
bcftools concat -Oz -o data/merged.vcf.gz data/chr*.vcf.gz
./scripts/analyze-vcf.sh -i data/merged.vcf.gz
```

Output lands in `out/`: a `*.stats.txt`, a PASS-filtered `*.pass.vcf.gz`, and a readable
`*.SUMMARY.md`. Hand the final VCF to the **Genome Analyzer** Claude Code skill for
interpretation.

## Roadmap
See `PLAN.md`.
