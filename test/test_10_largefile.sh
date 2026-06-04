# test_10_largefile.sh
# Regression: real files hold thousands of events; `head` closing a pipe early
# under set -o pipefail (SIGPIPE) must not abort detection or the analyzers.

big="$TMP/big.sv.vcf.txt"
{
  printf '##fileformat=VCFv4.1\n##source=GenerateSVCandidates 1.6.0\n'
  printf '##ALT=<ID=DEL,Description="Deletion">\n'
  printf '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="t">\n'
  printf '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="l">\n'
  printf '##INFO=<ID=END,Number=1,Type=Integer,Description="e">\n'
  printf '##INFO=<ID=EVENT,Number=1,Type=String,Description="ev">\n'
  printf '##FILTER=<ID=PASS,Description="ok">\n'
  printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
  printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
  awk 'BEGIN{for(i=1;i<=3000;i++){p=i*1000; printf "1\t%d\tD%d\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;SVLEN=-%d;END=%d\tGT\t0/1\n", p, i, 500+i, p+500+i}}'
} > "$big"
bigvcf="$(make_fixture "$big")"
out="$TMP/big_out"; mkdir -p "$out"
"$REPO/analyze.sh" -i "$bigvcf" -o "$out" >/dev/null 2>&1; rc=$?
assert_eq "0" "$rc" "analyze.sh survives a 3000-event SV file under set -e (no SIGPIPE abort)"
[[ -s "$out/big.sv.SUMMARY.md" ]] && w=yes || w=no
assert_eq "yes" "$w" "large SV run wrote a SUMMARY.md"
rows="$(sed -n '/## Largest events/,/## Filtering/p' "$out/big.sv.SUMMARY.md" | grep -c '^| DEL ' || true)"
assert_eq "10" "$rows" "Largest events table capped at 10 rows"
