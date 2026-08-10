# test_13_gwas_normalize.sh — fetch-gwas.sh::normalize_gwas (no network).
source "$REPO/fetch-gwas.sh" --source-only

out="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw.tsv" 2>/dev/null)"
assert_eq "3" "$(printf '%s\n' "$out" | grep -c .)" "normalize: the 3 clean rows survive"
assert_contains "$out" "rs100	G	Height	1e-12	1.2	GENEH	111	0.23" "normalize: raf kept as 8th column"
assert_contains "$out" "rs101	A	TraitNR	2e-12	1.1	GENEN	555	NA" "normalize: NR raf -> NA"
assert_contains "$out" "rs102	C	TraitEA	3e-12	1.0	GENEE	666	NA" "normalize: annotated '0.23 (EA)' raf -> NA"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Multi locus')"    "normalize: multi-SNP row dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Unknown allele')" "normalize: '?' risk allele dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Indel risk')"     "normalize: multi-base risk allele dropped"

# Header drift: col 27 is not RISK ALLELE FREQUENCY -> warn, raf always NA.
drift_out="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw_drift.tsv" 2>/dev/null)"
assert_contains "$drift_out" "rs100	G	Height	1e-12	1.2	GENEH	111	NA" "normalize: header drift -> raf NA, not a wrong column"
drift_err="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw_drift.tsv" 2>&1 >/dev/null)"
assert_contains "$drift_err" "WARNING" "normalize: header drift warns on stderr"
