# vcf-triage Skill — Design

**Date:** 2026-08-14
**Status:** approved (design discussion in-session)
**Fulfills:** PLAN.md Phase 3 — "Author a local skill `vcf-triage` that codifies our filter thresholds + report format"

## Goal

A project Claude skill that turns a completed pipeline run (analyze.sh + xref-traits.sh
outputs) into a full plain-language report: technical detail and layman interpretation
paired per section, triage driven by codified thresholds, heavily caveated, never
diagnostic. Audience is a genetics newcomer; genetic counseling is the standing pointer
for anything that looks actionable.

## Components

### 1. `.claude/skills/vcf-triage/SKILL.md`

Frontmatter `name: vcf-triage`, `description` triggering on: triage / interpret /
report / explain results over genome analysis outputs.

Workflow it instructs:

1. **Locate the run.** Accept an explicit basename or directory argument; otherwise
   use the newest `*.traits.SUMMARY.md` under `work/` then `out/`. Derive the run
   basename from it. Report which run was selected.
2. **Run the extract.** `.claude/skills/vcf-triage/extract-triage.sh <dir> <base>` —
   never read `*.traits.gwas.tsv` directly (it can exceed 20 MB).
3. **Read the small inputs directly:** the per-type `*.SUMMARY.md` files (snp-indel,
   sv, cnv) when present.
4. **Write the report** to `<dir>/<base>.report.md` following the embedded template.
5. **Deliver** the report file to the user (render), restating the headline findings
   in chat.

Hard rules (verbatim in the skill):

- Never present anything as a diagnosis; every section carries the
  research-grade framing.
- Anything flagged as potentially actionable points to clinical confirmation and
  genetic counseling.
- Explicitly defuse scary-named common variants (frequency is the triage tool).
- Heterozygous hits for recessive conditions are always framed as carrier status.
- No pronoun/relationship assumptions about the sample's person; refer to
  "the sample" or the name/pronouns the user has used.

### 2. `.claude/skills/vcf-triage/extract-triage.sh`

Bash 3.2 compatible (macOS), awk-based, no new dependencies. Sources `common.sh`
conventions where useful but must run standalone.

**Invocation:** `extract-triage.sh <dir> <base>` (also accepts a single path to any
run artifact and derives dir/base). Exit 0 always, except usage errors (exit 1).

**Inputs** (all optional; absence is reported, never fatal):

| File | Columns (headerless, tab-separated) |
|---|---|
| `<base>.traits.clinvar.tsv` | chrom:pos, allele, zygosity, significance, condition (\|-sep), review_status, gene, pop_af, af_source |
| `<base>.traits.gwas.tsv` | rsid, trait, genotype, risk_allele, copies, p, or_beta, gene, pmid, flag (ok/ambiguous), raf |
| `<base>.traits.pharmgkb.tsv` | future — presence noted as "PharmGKB section available" when it exists |

**Output:** plain text to stdout, small (target < 200 lines), sections:

```
== INPUTS ==            per-file: found/missing (+ which fetch/xref step enables it)
== CLINVAR PRIORITY ==  pathogenic-class rows, rare (<1%) or unknown AF first,
                        ordered by review-status weight, full row detail
== CLINVAR PHARMA ==    drug_response rows with review status
                        reviewed_by_expert_panel, full row detail
== CLINVAR COMMON ==    remaining pathogenic-class rows: count + one line each
                        (gene, significance, bucket) — the "defuse" list
== GWAS HEADLINE ==     counts: total, homozygous, ambiguous-unscored
== GWAS NOTABLE ==      flag=ok rows with numeric or_beta ≥ 2.0 or ≤ 0.5,
                        deduplicated by rsid (keep strongest p), sorted by
                        effect size, capped at 40 rows with a "N more" note
== PHARMGKB ==          present/absent note (future table)
```

**Codified thresholds** (single source of truth for the skill):

- Frequency buckets: rare < 1 %, low-frequency 1–5 %, common ≥ 5 %, unknown = NA/empty.
- Review-status weight: `reviewed_by_expert_panel` (4) > `criteria_provided,_multiple_submitters,_no_conflicts` (3) > `criteria_provided,_single_submitter` / `criteria_provided,_conflicting_classifications` (2) > `no_assertion_criteria_provided` (1).
- Priority = pathogenic-class AND (rare or unknown) — everything else pathogenic-class
  goes to the defuse list.
- GWAS notable = flag `ok` or `strand_flipped` AND or_beta parses as a number
  AND ≥ 2.0. (Strand-flipped hits are unambiguously resolved and fully scored
  by the pipeline; only `ambiguous` rows are unscored and excluded.)
  (The ≤ 0.5 protective side is deliberately dropped: the catalog's OR/beta
  column mixes odds ratios and betas, and a small beta is indistinguishable
  from a protective OR — including it would flood the report with
  molecular-trait rows. Betas rarely reach 2, so the ≥ 2.0 side stays clean.)

### 3. Report template (embedded in SKILL.md)

Generalized from `work/genome-analysis-report-2026-08-14.md` — structure only, no
data-derived variant names or sample details (same policy as commit eeda86f):

1. Header: sample, run date, databases, research-grade warning box
2. Data quality (Ti/Tv, PASS counts) — science + plain terms
3. Structural variants — science + plain terms
4. Copy-number variants — science + plain terms
5. Pharmacogenomics table (from CLINVAR PHARMA) — the "most useful" section
6. ClinVar pathogenic hits triaged by rarity: priority findings (carrier framing),
   then the defused-common list
7. GWAS: why the big number is normal + NOTABLE rows worth a plain-language callout
8. Caveats (single-copy sensitivity, strand-ambiguity, repeat regions, database
   label quality, no polygenic scores / star alleles)
9. Bottom line: action-tier table (reassurance / mention-if-drug-proposed /
   family-planning-only / routine-screening / none)

Sections whose inputs are missing collapse to one line: what's missing and which
command produces it.

## Testing

`test/test_17_vcf_triage_extract.sh` in the existing test style (synthetic fixtures,
no private data):

- exits 0 with all inputs present; sections appear in order
- rare + unknown pathogenic rows land in PRIORITY; common pathogenic rows in COMMON
- expert-panel drug_response rows land in PHARMA (lower review tiers do not)
- GWAS: ambiguous rows excluded, or_beta 2.4 included, beta 0.08 excluded,
  rsid dedup keeps the stronger p
- missing gwas.tsv / clinvar.tsv → INPUTS notes it, other sections still emit,
  exit 0
- cap: > 40 notable rows → 40 emitted + "more" note

The SKILL.md itself is prose and is exercised manually (run the skill on the
existing real outputs and compare against the hand-written report).

## Privacy

Skill + script + tests contain zero personal data and merge to main. Reports land
next to the run outputs in git-ignored directories.

## Out of scope

- PharmGKB cross-reference itself (separate open TODO 3636); the extract only
  reserves the section.
- Polygenic scores, star-allele calling, BioMCP/AlphaGenome integration.
