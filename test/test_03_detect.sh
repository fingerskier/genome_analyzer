source "$REPO/common.sh"
source "$REPO/analyze.sh" --source-only   # defines functions without running main
for pair in "snp-indel:snp-indel" "sv:sv" "cnv:cnv"; do
  f="${pair%%:*}"; want="${pair##*:}"
  vcf="$(make_fixture "$REPO/test/fixtures/$f.vcf.txt")"
  got="$(detect_type "$vcf")"
  assert_eq "$want" "$got" "detect_type($f) == $want"
done
