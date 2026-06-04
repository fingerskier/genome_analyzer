# test_13_gwas_normalize.sh — fetch-gwas.sh::normalize_gwas (no network).
source "$REPO/fetch-gwas.sh" --source-only

out="$(normalize_gwas < "$REPO/test/fixtures/gwas_raw.tsv")"
assert_eq "1" "$(printf '%s\n' "$out" | grep -c .)" "normalize: only the 1 clean row survives"
assert_contains "$out" "rs100	G	Height	1e-12	1.2	GENEH	111" "normalize: clean row -> 7 cols, rsid/allele split on last dash"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Multi locus')"   "normalize: multi-SNP row dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Unknown allele')" "normalize: '?' risk allele dropped"
assert_eq "0" "$(printf '%s\n' "$out" | grep -c 'Indel risk')"     "normalize: multi-base risk allele dropped"
