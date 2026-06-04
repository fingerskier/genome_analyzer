# test_15a_fetch_clinvar.sh — fetch-clinvar.sh sources cleanly + sane defaults.
source "$REPO/fetch-clinvar.sh" --source-only
assert_contains "$CLINVAR_URL" "clinvar.vcf.gz" "fetch-clinvar: default URL points at clinvar.vcf.gz"
assert_contains "$CLINVAR_URL" "vcf_GRCh38"     "fetch-clinvar: uses the GRCh38 build"
assert_eq "0" "$(bash -n "$REPO/fetch-clinvar.sh"; echo $?)" "fetch-clinvar: passes bash -n syntax check"
