for f in snp-indel sv cnv; do
  vcf="$(make_fixture "$REPO/test/fixtures/$f.vcf.txt")"
  hdr="$(bcftools view -h "$vcf" 2>&1)"
  assert_contains "$hdr" "##fileformat" "fixture $f is a readable VCF"
done
