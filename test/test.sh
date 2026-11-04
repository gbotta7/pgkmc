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
# 1) runs to completion (exit 0) and writes a non-empty output on full assemblies
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

# 3) runs to completion (exit 0) and writes a non-empty output on specified bed regions
if ! "$PGKMC" detect -t3 -k31 -b "$ROOT"/test/data/bed_list.txt -o "$TMP"/snpmers.txt "$ROOT"/test/data/test*.fa.gz 2> "$TMP"/detect.log; then
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

# 5) fails with an error if the input files are missing
expect_fail "no input files" "Usage: pgkmc detect" "$PGKMC" detect
expect_fail "options but no input files" "Usage: pgkmc detect" "$PGKMC" detect -k 31 -o "$TMP/snpmers.txt"

# 6) fails with an error if k is too long or even
expect_fail "even k" "must be odd and <=31" "$PGKMC" detect -k20 "$ROOT"/test/data/test*.fa.gz
expect_fail "k too large (33)" "must be odd and <=31" "$PGKMC" detect -k33 "$ROOT"/test/data/test*.fa.gz

# 7) fails with a prefix that is too small
expect_fail "p too small" "must not be smaller than" "$PGKMC" detect -p 1 "$ROOT"/test/data/test*.fa.gz

# 8) fails if a bed file does not exist
expect_fail "bed list does not exist" "failed to open" "$PGKMC" detect -k 31 -b "$TMP/nope.txt" "$ROOT"/test/data/test*.fa.gz

# 9) fails if there are fewer bed files than assemblies
echo "$ROOT/test/a.bed" > "$TMP/one_bed.txt"
expect_fail "fewer BEDs than FASTAs" "BED files given for" "$PGKMC" detect -k 31 -b "$TMP/one_bed.txt" "$ROOT"/test/data/test*.fa.gz

if [ "$n_fail" -ne 0 ]; then
    echo "$n_fail test(s) failed"; exit 1
fi
echo "All tests passed"