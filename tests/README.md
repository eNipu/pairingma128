# tests/ — pairing mathematical-correctness suites

Test-driven-dev gate for the pairingma128 refactor: each suite asserts
**mathematical identities** that any correct pairing implementation must
satisfy. Verdict = process exit code (0 = all checks passed), not stdout text.

## Files

- `math_assert.h`      — tiny CHECK/report harness
- `bn_math_test.c`     — BN optimal-ate battery (16 checks)
- `bls12_math_test.c`  — BLS12 optimal-ate battery (17 checks)
- `kss16_math_test.c`  — KSS16 batteries incl. all 5 pairing entry points (32)
- `run_math_tests.sh`  — build + run one or all suites
- `golden/bls12_digest.txt` — byte-exact pairing value on fixed vectors

## Run

    bash tests/run_math_tests.sh            # all curves
    bash tests/run_math_tests.sh bn bls12   # subset

GMP is expected under /opt/homebrew (Homebrew). Override:

    MATH_GMP=-L/usr/local/lib bash tests/run_math_tests.sh

The runner builds each curve's split modules (field/curve/pairing/params/util)
plus the test file, links against GMP, runs, and checks exit code + golden
digest parity. Each curve directory also has its own `CMakeLists.txt` with
`<curve>` and `<curve>_math_test` targets.

## What is verified (see also the plan document)

- r prime, r | #E(Fp); points on the curve equation and of exact order r
  (cofactor-cleared in the tests, since the G1 generators return raw points);
- pairing output non-degenerate, in the order-r subgroup (z^r == 1),
  bilinear (e([a]P,[b]Q) == z^(ab)), alternating in both slots, deterministic;
- field Frobenius pi(x) == x^p (top extension level);
- KSS16: Pseudo/Sparse/Optimal flavors produce identical values; Ate/Tate pass
  their own full battery (valid but differently-normalized pairings);
- BLS12 golden digest parity (byte-identical pairing value across runs, and
  pre- vs post-refactor).

## Current state (post-refactor)

All three curves are split into modules. `run_math_tests.sh` links the split
modules (excluding main.c/test.c) with the test TU; the `#include`s in each
`*_math_test.c` point at the new module headers (field/curve/pairing/params).
No assertion logic changed from the pre-refactor baseline.
