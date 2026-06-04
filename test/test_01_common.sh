source "$REPO/common.sh"
vcf="$(make_fixture "$REPO/test/fixtures/snp-indel.vcf.txt")"
rm -f "${vcf}.tbi"          # force re-index path
validate_index "$vcf"
[[ -f "${vcf}.tbi" ]] && r=yes || r=no
assert_eq "yes" "$r" "validate_index creates a .tbi"
have bcftools && hb=yes || hb=no
assert_eq "yes" "$hb" "have detects bcftools on PATH"
