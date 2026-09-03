#!/usr/bin/env bash
# run_math_tests.sh - build & run the mathematical correctness suites.
#
# Pre-refactor mode (default): compiles each curve's monolithic .c + the test
# file. The monolith's own `main` is renamed via -Dmain so the test file's
# assertions drive the program, and -fcommon merges the tentative-definition
# globals declared in the curve headers between the two translation units.
#
# Usage:  tests/run_math_tests.sh [curve ...]     # curves: bn bls12 kss16
#         MATH_GMP=-L/path/to/gmp                 # override gmp link flags
set -u
cd "$(dirname "$0")/.." || exit 1
ROOT="$PWD"
OUT="$(mktemp -d /tmp/pairingma128_tests.XXXXXX)"
GMP_LDFLAGS="${MATH_GMP:--L/opt/homebrew/lib}"
GMP_INC="-I/opt/homebrew/include"
CC="${CC:-cc}"
CFLAGS="-O2 -fcommon"

run_one() { # $1=name  $2=dir  $3=header  $4=source
    local name="$1" dir="$2" src="$4"
    local obj_exe="$OUT/${name}_math_test"
    local digest_mismatch=0
    if [ -f "$dir/field.c" ]; then
        # ---- post-refactor: link split modules + test TU (main.c excluded) ----
        echo "==== $name: building (split modules + test TU) ===="
        local mod objs=()
        for mod in field curve pairing params util; do
            if [ -f "$dir/$mod.c" ]; then
                "$CC" -O2 $GMP_INC -I"$dir" -c "$dir/$mod.c" -o "$OUT/${name}_mod_${mod}.o" || { echo "$name BUILD FAILED ($mod)"; return 2; }
                objs+=("$OUT/${name}_mod_${mod}.o")
            fi
        done
        "$CC" -O2 $GMP_INC -I"$dir" -I"$ROOT/tests" -c "tests/${name}_math_test.c" -o "$OUT/${name}_test.o" || { echo "$name BUILD FAILED (test)"; return 2; }
        "$CC" "${objs[@]}" "$OUT/${name}_test.o" $GMP_LDFLAGS -lgmp -o "$obj_exe" || { echo "$name BUILD FAILED (link)"; return 2; }
    else
        # ---- pre-refactor: single monolith .c owns main ----
        echo "==== $name: building (monolith + test TU) ===="
        "$CC" $CFLAGS $GMP_INC -Dmain=monolith_main -I"$dir" -c "$dir/$src" -o "$OUT/${name}_mono.o" || { echo "$name BUILD FAILED (monolith)"; return 2; }
        "$CC" $CFLAGS $GMP_INC -I"$dir" -I"$ROOT/tests" -c "tests/${name}_math_test.c" -o "$OUT/${name}_test.o" || { echo "$name BUILD FAILED (test)"; return 2; }
        "$CC" "$OUT/${name}_mono.o" "$OUT/${name}_test.o" $GMP_LDFLAGS -lgmp -o "$obj_exe" || { echo "$name BUILD FAILED (link)"; return 2; }
    fi
    echo "==== $name: running ===="
    "$obj_exe" > "$OUT/${name}.log" 2>&1
    local rc=$?
    echo "==== $name: exit code = $rc ===="
    tail -n 12 "$OUT/${name}.log"
    echo "==== $name: full log: $OUT/${name}.log ===="
    # golden digest parity (deterministic pairing values, e.g. BLS12 fixed vectors)
    if [ -f "$ROOT/tests/golden/${name}_digest.txt" ]; then
        if grep -E '^  z_0|^  z_1' "$OUT/${name}.log" | diff -q - "$ROOT/tests/golden/${name}_digest.txt" > /dev/null; then
            echo "==== $name: golden digest MATCHES ===="
        else
            echo "==== $name: golden digest MISMATCH ===="
            digest_mismatch=1
        fi
    fi
    if [ $rc -ne 0 ] || [ $digest_mismatch -ne 0 ]; then return 1; fi
    return 0
}

overall=0
for c in "${@:-bn bls12 kss16}"; do
    case "$c" in
        bn)    run_one bn    BN    BN12_header.h  EFp12_optimal.c || overall=1 ;;
        bls12) run_one bls12 BLS12 BLS_12.h       main.c          || overall=1 ;;
        kss16) run_one kss16 KSS16 KSS_16.h       main.c          || overall=1 ;;
        *) echo "unknown curve: $c" >&2; exit 2 ;;
    esac
done
exit $overall
