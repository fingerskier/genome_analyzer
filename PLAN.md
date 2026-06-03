# PLAN — genome-tools

A small, personal repo of tooling for interpreting germline sequencing data. Tooling lives
in git; **genomic data never does** (see `.gitignore` and the privacy note below).

## Goal
Turn provider deliverables (`.vcf.gz` + `.fq.gz`) into readable, well-caveated insight, with
Claude Code skills/MCP doing the interpretation heavy-lifting.

## Phase 0 — Scaffold (done)
- [x] Repo layout: `scripts/`, `data/` (gitignored), `out/` (gitignored)
- [x] `analyze-vcf.sh`: validate → stats → assay detection → PASS filter → summary
- [x] `.gitignore` excludes all sequence/variant data
- [ ] Pin tool versions in `env.md` (bcftools, htslib, SnpEff/VEP) for reproducibility

## Phase 1 — Inventory the delivery
- [ ] `bcftools query -l` each VCF → confirm single sample, consistent sample ID
- [ ] Classify each `.vcf.gz`: per-chromosome shard vs SNV/indel/SV split
      - shards → `bcftools concat -Oz -o merged.vcf.gz *.vcf.gz`
- [ ] Run `analyze-vcf.sh` on the merged/primary VCF → first `SUMMARY.md`
- [ ] Verify FASTQ integrity: `seqkit stats *.fq.gz` (read counts, length, %GC)
- [ ] Record provider, reference build (GRCh37 vs GRCh38!), and pipeline if documented

## Phase 2 — Annotation layer
- [ ] Decide engine: **SnpEff** (fast setup, DB download) vs **VEP** (richer, cache is ~25 GB)
- [ ] Confirm build matches the VCF — annotating GRCh38 calls against a GRCh37 DB is silently wrong
- [ ] Add gnomAD frequency annotation → flag common variants (almost all benign) for triage
- [ ] Wire ClinVar + PharmGKB cross-reference

## Phase 3 — Interpretation via Claude
- [ ] Install **Genome Analyzer** skill; feed it the annotated PASS VCF
- [ ] Add **BioMCP**: `claude mcp add biomcp -- uv run --with biomcp-python biomcp run`
- [ ] (optional) **AlphaGenome** via BioMCP for regulatory effect on "uncertain" variants
      (needs `ALPHAGENOME_API_KEY`)
- [ ] Author a local skill `vcf-triage` that codifies our filter thresholds + report format

## Phase 4 — Re-analysis from FASTQ (only if needed)
Trigger this phase only if you want structural variants, a fresh caller, or a newer reference —
the provider VCF is SNV/indel-only.
- [ ] Align: `BWA-MEM2` (or `minimap2`) to GRCh38 with read groups
- [ ] Call: **DeepVariant** (cleaner than provider GATK, easy to run)
- [ ] Structural variants / CNV: **Manta**, **GRIDSS**
- [ ] Consider nf-core/sarek to run the whole germline path with sane defaults
- [ ] Compute budget: WGS alignment needs ~30–40 GB RAM + many CPU-hours + ~100 GB scratch.
      The Nitro 5 won't comfortably do WGS alignment — plan a short cloud burst for this step.

## Phase 5 — Reproducibility & hygiene
- [ ] Containerize the toolchain (Docker/Apptainer) so runs are deterministic
- [ ] `Makefile` or `cross-x` job wrapping the stages for one-command runs
- [ ] Checksums (`sha256`) of source deliverables stored in repo (the hashes, not the data)

## Repo layout
```
genome-tools/
├── README.md
├── PLAN.md
├── .gitignore          # excludes ALL genomic data
├── scripts/
│   └── analyze-vcf.sh
├── data/               # gitignored — raw deliverables live here
└── out/                # gitignored — generated artifacts
```

## Privacy note (read before first commit)
This is another person's genome. Treat it as the most sensitive data you handle.
- Never commit `*.vcf*`, `*.fq*`, `*.bam`, or summaries containing real coordinates/genotypes.
- A private repo is not sufficient protection for the raw data — keep data out of git entirely.
- If you ever share a `SUMMARY.md`, confirm she's comfortable with it; genetic risk info is hers.
- Genetic counseling is the right venue for anything that looks clinically actionable. The tools
  here produce research-grade associations, not diagnoses.
