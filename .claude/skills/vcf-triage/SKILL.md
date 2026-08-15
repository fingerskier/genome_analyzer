---
name: vcf-triage
description: Use when asked to triage, interpret, explain, or report on genome-analyzer run outputs (analyze.sh / xref-traits.sh results — *.SUMMARY.md and *.traits.*.tsv files). Produces a full report pairing technical detail with layman interpretation, using codified triage thresholds. Never diagnoses.
---

# vcf-triage — plain-language triage report over pipeline outputs

## Workflow

1. **Locate the run.** If the user names a basename, directory, or file, use it.
   Otherwise take the newest `*.traits.SUMMARY.md` under `work/`, then `out/`.
   The run basename is the filename up to `.traits.`. Tell the user which run
   you selected.
2. **Run the extract** (never read `*.traits.gwas.tsv` directly — it can exceed
   20 MB):

   ```bash
   .claude/skills/vcf-triage/extract-triage.sh <dir> <base>
   ```

   Its sections are the complete, already-thresholded evidence base:
   INPUTS, CLINVAR PRIORITY / PHARMA / COMMON, GWAS HEADLINE / NOTABLE,
   PHARMGKB. Trust its selection; do not re-filter the TSVs yourself.
3. **Read the small run summaries directly** — the per-type `*.SUMMARY.md`
   files listed under INPUTS (SNP/indel QC, SV, CNV).
4. **Write the report** to `<dir>/<base>.report.md` following the template
   below. These directories are git-ignored; never commit a report.
5. **Deliver** the report file to the user (render it) and restate the headline
   findings in chat in plain language.

## Hard rules

- **Never diagnose.** Every section keeps the research-grade framing; the
  report opens with the warning box. Nothing here changes medical decisions
  without clinical confirmation.
- Anything potentially actionable points to **clinical confirmation and
  genetic counseling** — that is always the next step, never your
  interpretation alone.
- **Frequency defuses.** Explicitly explain that scary-named common variants
  (CLINVAR COMMON list) are population-normal database noise. Rarity is a
  reason to look closer, not a verdict.
- **Heterozygous + recessive condition = carrier status**, and say so: one
  working copy means the person is healthy; relevance is family planning.
- Weight ClinVar claims by review status: expert panel > multiple submitters >
  single submitter > no assertion criteria. Say when support is weak.
- Frame the GWAS total (often 100k+) as normal for every genome; only the
  NOTABLE rows deserve prose, and only disease-relevant ones at that.
- No pronoun or relationship assumptions about the sample's person; use the
  name/pronouns the user used, otherwise "the sample".
- Missing inputs (per INPUTS/PHARMGKB sections) collapse to one line in the
  report: what is missing and which command produces it.

## Report template

Pair every section: **The science** (metrics, thresholds, review tiers,
p-values) then **In plain terms** (what it means for a person, no jargon).

1. **Header** — sample, run date, databases used, then a blockquote warning:
   research-grade exploration, not a diagnosis or medical advice.
2. **Data quality** — records/SNP/indel counts, Ti/Tv vs expected ~2.0–2.1,
   PASS fraction. Plain terms: is the data trustworthy; ~4–5M variants is
   normal for every human.
3. **Structural variants** — counts by type, sizes, gene overlap. Plain terms:
   thousands of small SVs are ordinary; large events near immune (HLA) regions
   are expected diversity.
4. **Copy-number variants** — losses/gains, burden in Mb, largest events.
   Plain terms: flag that repeat-rich regions (centromeres, acrocentric arms)
   are where callers are least reliable.
5. **Pharmacogenomics** (CLINVAR PHARMA + PHARMGKB when present) — table of
   gene / genotype / drug(s) / frequency. The most practically useful section:
   these change how the body handles specific drugs, matter only if the drug
   is ever proposed, and then warrant confirmatory clinical testing.
6. **ClinVar findings triaged by rarity** — CLINVAR PRIORITY rows first
   (carrier framing for recessives, review-status caveats), then the COMMON
   defuse list summarized in a sentence or two.
7. **GWAS** — HEADLINE counts with the why-a-huge-number-is-normal framing;
   NOTABLE rows worth individual plain-language callouts (odds ratio ≥ 2
   disease associations), each with lifestyle/screening context where honest.
8. **Caveats** — single-variant sensitivity of consumer WGS, strand-ambiguous
   SNPs unscored, repeat-region CNV reliability, database label quality, no
   polygenic scores or star-allele calling.
9. **Bottom line** — action-tier table: reassurance / mention-if-drug-proposed /
   family-planning-only / routine-screening / no action. Close plainly: is
   anything urgent (usually no), does anything suggest a present health
   problem (usually no).
