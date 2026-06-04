#!/usr/bin/env bash
# test_04_shim.sh — verify analyze-vcf.sh forwards to analyze.sh
out="$("$REPO/analyze-vcf.sh" -i "$(make_fixture "$REPO/test/fixtures/snp-indel.vcf.txt")" -o "$TMP/shim_out" 2>/dev/null)"
[[ -s "$TMP/shim_out/snp-indel.SUMMARY.md" ]] && e=yes || e=no
assert_eq "yes" "$e" "analyze-vcf.sh shim still produces a summary"
