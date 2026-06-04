# test_14_dosage.sh — allele dosage + strand logic (lib/trait-xref.sh).
source "$REPO/lib/trait-xref.sh"

# risk_dosage <ref> <alt> <het|hom> <risk> -> "<dosage>\t<flag>"
assert_eq "1	ok"             "$(risk_dosage A G het G)" "dosage: het, risk=alt -> 1 copy"
assert_eq "1	ok"             "$(risk_dosage A G het A)" "dosage: het, risk=ref -> 1 copy"
assert_eq "2	ok"             "$(risk_dosage A G hom G)" "dosage: hom-alt, risk=alt -> 2 copies"
assert_eq "0	ok"             "$(risk_dosage A G hom A)" "dosage: hom-alt, risk=ref -> 0 copies"
assert_eq "1	strand_flipped" "$(risk_dosage A G het C)" "dosage: risk=C matches alt G via complement (non-palindromic)"
assert_eq "NA	ambiguous"      "$(risk_dosage A T het A)" "dosage: A/T palindromic SNP -> ambiguous, not scored"
assert_eq "NA	ambiguous"      "$(risk_dosage C G hom C)" "dosage: C/G palindromic SNP -> ambiguous, not scored"
assert_eq "1	strand_flipped" "$(risk_dosage A G het T)" "dosage: risk=T resolves to ref A via complement -> 1 copy, flipped"
