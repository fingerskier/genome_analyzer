# test_09_fetch.sh — fetch-genes.sh: gtf_to_bed parsing (no network).

source "$REPO/fetch-genes.sh" --source-only
# gtf_to_bed reads GENCODE GTF on stdin, writes chrom\tstart\tend\tgene to stdout, chr-normalized.
sample='chr1	HAVANA	gene	1000	2000	.	+	.	gene_name "GENEA"; gene_type "protein_coding";'
got="$(printf '%s\n' "$sample" | gtf_to_bed)"
assert_eq "1	999	2000	GENEA" "$got" "gtf_to_bed: parses + strips chr + 0-bases start"
