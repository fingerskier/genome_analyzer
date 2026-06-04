# test_09_fetch.sh — fetch-genes.sh: gtf_to_bed parsing (no network).

source "$REPO/fetch-genes.sh" --source-only
# gtf_to_bed reads GENCODE GTF on stdin, writes chrom\tstart\tend\tgene to stdout, chr-normalized.
sample='chr1	HAVANA	gene	1000	2000	.	+	.	gene_name "GENEA"; gene_type "protein_coding";'
got="$(printf '%s\n' "$sample" | gtf_to_bed)"
assert_eq "1	999	2000	GENEA" "$got" "gtf_to_bed: parses + strips chr + 0-bases start"

# collapse_genes: one interval per (gene name, chromosome). PAR genes stay split.
cin="$(printf '1\t500\t900\tGENEA\n1\t100\t400\tGENEA\nX\t60000\t61000\tSHOX\nY\t10000\t11000\tSHOX\n')"
cout="$(printf '%s\n' "$cin" | collapse_genes)"
assert_contains "$cout" "1	100	900	GENEA" "collapse: GENEA merged to widest span"
assert_contains "$cout" "X	60000	61000	SHOX" "collapse: SHOX kept on chrX"
assert_contains "$cout" "Y	10000	11000	SHOX" "collapse: SHOX kept on chrY (not merged across chromosomes)"
