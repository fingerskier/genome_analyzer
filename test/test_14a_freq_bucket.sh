# test_14a_freq_bucket.sh — population-frequency buckets (lib/trait-xref.sh).
source "$REPO/lib/trait-xref.sh"

# freq_bucket <af> -> "<bucket>\t<display>"
assert_eq "common	common (~40%)"               "$(freq_bucket 0.4)"   "bucket: 40% -> common"
assert_eq "common	common (~5%)"                "$(freq_bucket 0.05)"  "bucket: 5% boundary -> common"
assert_eq "low-frequency	low-frequency (~2%)" "$(freq_bucket 0.02)"  "bucket: 2% -> low-frequency"
assert_eq "low-frequency	low-frequency (~1%)" "$(freq_bucket 0.01)"  "bucket: 1% boundary -> low-frequency"
assert_eq "rare	rare (<1%)"                    "$(freq_bucket 0.001)" "bucket: 0.1% -> rare"
assert_eq "rare	rare (<1%)"                    "$(freq_bucket 5e-3)"  "bucket: scientific notation accepted"
assert_eq "unknown	unknown"                   "$(freq_bucket NA)"    "bucket: NA -> unknown"
assert_eq "unknown	unknown"                   "$(freq_bucket abc)"   "bucket: junk -> unknown"
assert_eq "unknown	unknown"                   "$(freq_bucket "")"    "bucket: empty -> unknown"
assert_eq "unknown	unknown"                   "$(freq_bucket 1.5)"   "bucket: >1 is not a frequency -> unknown"
