# test_11_perchrom.sh — verify primary contigs listed individually; alt/decoy collapsed.

pc="$TMP/pc.sv.vcf.txt"
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
  printf '1\t5000\tA\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;SVLEN=-500;END=5500\tGT\t0/1\n'
  printf '1\t9000\tB\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;SVLEN=-300;END=9300\tGT\t0/1\n'
  printf 'X\t1000\tC\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;SVLEN=-200;END=1200\tGT\t0/1\n'
  printf 'Un_JTFH01000116v1_decoy\t100\tD\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;SVLEN=-100;END=200\tGT\t0/1\n'
  printf 'Un_KI270442v1\t100\tE\tN\t<DEL>\t200\tPASS\tSVTYPE=DEL;SVLEN=-100;END=200\tGT\t0/1\n'
} > "$pc"

( source "$REPO/common.sh"; source "$REPO/lib/sv.sh"
  INPUT="$(make_fixture "$pc")"
  OUTDIR="$TMP/pc_out"; mkdir -p "$OUTDIR"; BASE="pc"; GENES_BED=""
  run_sv ) >/dev/null 2>&1

sec="$(sed -n '/## Per chromosome/,/## Largest/p' "$TMP/pc_out/pc.SUMMARY.md")"
assert_contains "$sec" "| 1 | 2 |" "perchrom: primary chr1 shown individually (2 events)"
assert_contains "$sec" "| X | 1 |" "perchrom: primary X shown individually"
assert_contains "$sec" "other contigs (2) | 2" "perchrom: 2 decoy contigs collapsed into one row"
