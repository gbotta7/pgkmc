#!/bin/bash

set -uo pipefail

# locate the binary regardless of where the script is called from
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PGKMC="$ROOT/src/pgkmc"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

n_fail=0

# expect_fail for tests that should fail with exit code 1 and a specific error message
expect_fail() {
    local desc=$1 pattern=$2; shift 2
    local out rc
    out=$("$@" 2>&1) && rc=0 || rc=$?
    if [ "$rc" -ne 1 ]; then
        echo "FAIL: $desc (expected exit 1, got $rc)"; n_fail=$((n_fail+1))
    elif ! grep -q -- "$pattern" <<<"$out"; then
        echo "FAIL: $desc (message '$pattern' not found)"; n_fail=$((n_fail+1))
    else
        echo "OK: $desc"
    fi
}

### TEST DETECT
# 1) runs to completion (exit 0) and writes a non-empty output that matches expectations on full assemblies
if ! "$PGKMC" detect -t3 -k31 -o "$TMP"/snpmers.txt "$ROOT"/test/data/test*.fa.gz 2> "$TMP"/detect.log; then
    echo "FAIL: detect exited with an error"; cat "$TMP/detect.log"; n_fail=$((n_fail+1));
elif [ ! -s "$TMP/snpmers.txt" ]; then
    echo "FAIL: detect produced an empty output"; n_fail=$((n_fail+1));
else
    echo "OK: detect ran to completion"

    # 2) the SNP-mers are the expected ones (sort: output order may depend on hashing/threads)
    if diff <(LC_ALL=C sort "$TMP/snpmers.txt") <(LC_ALL=C sort "$ROOT/test/results/expected_snpmers_full.txt") > "$TMP/snpmers_full.diff"; then
        echo "OK: detect SNP-mers match expected on full assemblies"
    else
        echo "FAIL: detect SNP-mers differ from expected (first 20 lines of diff):"
        head -n 20 "$TMP/snpmers_full.diff"; n_fail=$((n_fail+1));
    fi
fi

# 3) runs to completion (exit 0) and writes a non-empty output that matches expectations on specified bed regions
BEDS=( "$ROOT"/test/data/test*.bed )
printf '%s\n' "${BEDS[@]}" > "$TMP/bed_list.txt"
if ! "$PGKMC" detect -t3 -k31 -b "$TMP"/bed_list.txt -o "$TMP"/snpmers.txt "$ROOT"/test/data/test*.fa.gz 2> "$TMP"/detect.log; then
    echo "FAIL: detect exited with an error"; cat "$TMP/detect.log"; n_fail=$((n_fail+1));
elif [ ! -s "$TMP/snpmers.txt" ]; then
    echo "FAIL: detect produced an empty output"; n_fail=$((n_fail+1));
else
    echo "OK: detect ran to completion"

    # 4) the SNP-mers are the expected ones (sort: output order may depend on hashing/threads)
    if diff <(LC_ALL=C sort "$TMP/snpmers.txt") <(LC_ALL=C sort "$ROOT/test/results/expected_snpmers_bed.txt") > "$TMP/snpmers_bed.diff"; then
        echo "OK: detect SNP-mers match expected on specified bed regions"
    else
        echo "FAIL: detect SNP-mers differ from expected (first 20 lines of diff):"
        head -n 20 "$TMP/snpmers_bed.diff"; n_fail=$((n_fail+1));
    fi
fi

# 5) fails if the input files are missing
expect_fail "detect: no input files" "Usage: pgkmc detect" "$PGKMC" detect
expect_fail "detect: options but no input files" "Usage: pgkmc detect" "$PGKMC" detect -k 31 -o "$TMP/snpmers.txt"

# 6) fails when k is too long or even
expect_fail "detect: even k" "must be odd and <=31" "$PGKMC" detect -k20 "$ROOT"/test/data/test*.fa.gz
expect_fail "detect: k too large (33)" "must be odd and <=31" "$PGKMC" detect -k33 "$ROOT"/test/data/test*.fa.gz

# 7) fails when prefix is too small
expect_fail "detect: prefix too small" "must not be smaller than" "$PGKMC" detect -p 1 "$ROOT"/test/data/test*.fa.gz

# 8) fails if there are fewer bed files than assemblies
echo "$ROOT/test/a.bed" > "$TMP/one_bed.txt"
expect_fail "detect: fewer BEDs than assemblies" "BED files given for" "$PGKMC" detect -k 31 -b "$TMP/one_bed.txt" "$ROOT"/test/data/test*.fa.gz

# 9) fails for missing files
expect_fail "detect: bed list does not exist" "failed to open" "$PGKMC" detect -k 31 -b "$TMP/nope.txt" "$ROOT"/test/data/test*.fa.gz
expect_fail "detect: input file does not exist" "failed to open" "$PGKMC" detect -k 31 "$TMP/nope.fa"

### TEST COUNT
# 1) runs to completion (exit 0) and writes a non-empty output that matches expectations on full assemblies (test on first 10 SNP-mers and k-mers)
LC_ALL=C sort "$ROOT/test/results/expected_snpmers_full.txt" | head -n 10 > "$TMP/first10_snpmers.txt"
"$PGKMC" count --snp -w --kmers "$TMP/first10_snpmers.txt" -t1 -k31 "$ROOT"/test/data/test1.fa.gz 2>/dev/null > "$TMP/first10_snpmers_cts.tsv" # -t1 so the order of the position lists doesn't depend on threads

# 2) the SNP-mer counts are the expected ones
{ head -n 1 "$TMP/first10_snpmers_cts.tsv"; tail -n +2 "$TMP/first10_snpmers_cts.tsv" | LC_ALL=C sort; } > "$TMP/first10_snpmers_cts_sorted.tsv"
n_diff=$(LC_ALL=C comm -13 "$ROOT/test/results/expected_first10_snpmers_cts.tsv" "$TMP/first10_snpmers_cts_sorted.tsv" | wc -l)
if [ -s "$TMP/first10_snpmers_cts.tsv" ] && [ "$n_diff" -eq 0 ]; then
    echo "OK: count --snp -w on first 10 SNP-mers matches expected counts"
else
    echo "FAIL: count --snp -w on first 10 SNP-mers differs from expected counts:"
    LC_ALL=C comm -13 "$ROOT/test/results/expected_first10_snpmers_cts.tsv" "$TMP/first10_snpmers_cts_sorted.tsv" | head -n 20
    n_fail=$((n_fail+1))
fi

# 3) fails if the input files are missing
expect_fail "count: no input file" "Usage: pgkmc count" "$PGKMC" count
expect_fail "count: options but no input file" "Usage: pgkmc count" "$PGKMC" count -k 31 -o "$TMP/snpmers_cts.tsv"

# 4) fails when passed wrong number of positional arguments
expect_fail "count: multiple input files" "Usage: pgkmc count" "$PGKMC" count -k 31 "$ROOT"/test/data/test*.fa.gz

# 5) fails when k is too long, even, or different than what was used in pgkmc detect (both modes)
expect_fail "count: k too large (33) in --snp mode" "must be odd and <=31" "$PGKMC" count --snp -k33 "$ROOT"/test/data/test1.fa.gz
expect_fail "count: k too large (33)" "must be <=31" "$PGKMC" count -k33 "$ROOT"/test/data/test1.fa.gz
expect_fail "count: even k in --snp mode" "must be odd and <=31" "$PGKMC" count --snp --kmers "$ROOT/test/results/expected_snpmers_full.txt" -k20 "$ROOT"/test/data/test1.fa.gz
expect_fail "count: k is different than detect" "failed to parse SNP-mer" "$PGKMC" count --snp -k21 --kmers "$ROOT/test/results/expected_snpmers_full.txt" "$ROOT"/test/data/test1.fa.gz

# 6) fails in --snp mode without passing --kmers
expect_fail "count: no --kmers specified in --snp mode" "--kmers is required" "$PGKMC" count --snp -k31 "$ROOT"/test/data/test1.fa.gz

# 7) fails when min_qual is out of range
expect_fail "count: min_qual negative" "must be in \[0,93\]" "$PGKMC" count -k31 -q -1 "$ROOT"/test/data/test1.fa.gz
expect_fail "count: min_qual too large" "must be in \[0,93\]" "$PGKMC" count -k31 -q 94 "$ROOT"/test/data/test1.fa.gz

# 8) fails for missing files
expect_fail "count: kmers file does not exist in --snp mode" "failed to open" "$PGKMC" count --snp --kmers "$TMP/nope.txt" -k31 "$ROOT"/test/data/test1.fa.gz
expect_fail "count: bed file does not exist" "failed to open" "$PGKMC" count -k31 -b "$TMP/nope.bed" "$ROOT"/test/data/test1.fa.gz
expect_fail "count: input file does not exist" "failed to open" "$PGKMC" count -k31 "$TMP/nope.fa"


if [ "$n_fail" -ne 0 ]; then
    echo "$n_fail test(s) failed"; exit 1
fi
echo "All tests passed"