# Problem
Each curve implementation (`BN/`, `BLS12/`, `KSS16/`) is a single monolithic `.c` file (3.4k–5.3k lines) plus one header. Each header mixes together: global mutable state (`mpz_t prime`, precomputed basis/frobenius constants), finite-field tower struct definitions (`Fp`/`Fp2`/`Fp6`/`Fp12` for BN & BLS12; `Fp`/`Fp2`/`Fp4`/`Fp8`/`Fp16` for KSS16), elliptic-curve point arithmetic (`EFp*`), Miller-loop/pairing code, parameter generation, and test functions — all as prototypes in one header and definitions in one `.c` file, in that mixed order. There is no build system currently checked in (`build/` contains stale CMake cache files referencing a `CMakeLists.txt` that no longer exists in the repo and isn't tracked by git) and no automated tests.


# Goal
Decouple the finite-field arithmetic layer (Fp and its extension towers) from curve-point arithmetic, pairing/Miller-loop logic, parameter generation, and tests — within each curve directory independently (no cross-curve sharing in this pass, per your choice). Preserve behavior exactly; this is cryptographic code, so correctness must be verifiable, not just "it compiles."

# Target layout (applied identically to BN/, BLS12/, KSS16/)
For each curve directory, split the current `header.h` + `main.c`/`EFp12_optimal.c` into:
* `field.h` / `field.c`: struct definitions for the field tower (`Fp`, `Fp2`, `Fp6`, `Fp12` for BN/BLS12; `Fp`, `Fp2`, `Fp4`, `Fp8`, `Fp16` for KSS16) and all pure field arithmetic (init/set/random/clear/print, add/sub/mul/inv/legendre/sqrt/pow/cmp, frobenius maps, and KSS16's precomputed frobenius constants + `pre_calculate_frob_p*` functions, since those operate purely at the field level). Depends only on the global prime/basis constants declared in `params.h`.
* `curve.h` / `curve.c`: `EFp`/`EFp2`/.../`EFp12|EFp16` struct definitions and point arithmetic (`ECD`, `ECA`, `SCM*`, `rational_point`, `frobenius`, `skew_frobenius`, twist mappings, `generate_G1`/`generate_G2`). Depends on `field.h`.
* `pairing.h` / `pairing.c`: sparse mappings, all Miller-loop variants, `Final_exp*`/`Final_Exp`/`final_exp_hard`, and the top-level pairing entry points (`Opt_ate_pairing`, `Ate_pairing`, `Tate_Pairing`, KSS16's `Pseudo_*`/`Sparse_*` families). Depends on `curve.h` and `field.h`.
* `params.h` / `params.c`: global parameter state (`prime`/`PRIME_P`, `EFp_order`, basis constants, etc.) as `extern` declarations in the header with single definitions in the `.c` file, plus `init_parameters`/`set_parameters`/`clear_parameters`/`generate_*`/`get_*`/`weil`/`print_parameters`/`KSS_16_parameters`/`pre_calc_vector_final_exp`. This becomes the one place that defines global state; `field.h`/`curve.h` only declare `extern` references to the parts they use.
* `test.h` / `test.c`: all `test_*`, `check_num_of_Fp_mul`, `check_Pairing`, `Check_SCM`, `rational_point_check`.
* `util.h` / `util.c`: `timedifference_msec`/`timedifference_usec` and (KSS16 only) the diagnostic SCM/NAF counters currently declared loose in `main.c`.
* `main.c`: reduced to just `int main()` plus any KSS16-specific top-level orchestration currently in `main()`, calling into the modules above.

Key mechanical concern: today's headers contain **tentative definitions** of globals (e.g. `mpz_t prime;`, `struct Fp Fp_basis;`) directly in the header, which only works because each header is included by exactly one `.c` file. Once multiple `.c` files include shared headers, these must become `extern` declarations in headers with exactly one defining `.c` file (`params.c`, or `field.c`/`util.c` for the few globals that live there), or the build will fail with duplicate-symbol/link errors under `-fno-common` (the default on modern gcc/clang).

# Verification strategy (no existing tests/build system)
1. Before touching a curve directory, compile its current single-file form as-is with `cc *.c -I/opt/homebrew/include -L/opt/homebrew/lib -lgmp -o <curve>_before`, run it, and save stdout as a baseline. The existing `test_*`/`check_Pairing` functions print a deterministic pass/fail marker (`"test success"`, `"success."`, `"Success"`, or their failure counterparts) even though the underlying random scalars differ per run — use presence/absence of failure markers plus overall exit code as the regression signal, not full stdout diffing (values are RNG/time-seeded and will legitimately differ between runs).
2. Add a minimal `CMakeLists.txt` (one per curve directory, or one top-level with three executables) so each split module set builds into a single `<curve>` executable, linked against system GMP. This replaces the missing/stale build config and gives durable, repeatable verification going forward.
3. After splitting a curve directory, rebuild and rerun; confirm the same pass/fail markers appear and the binary still runs to completion without crashes/leless obviously-wrong output (e.g. `prime length`, curve equation printout should be identical bit-lengths/values to baseline since these are deterministic, not random).
4. Do curves one at a time in this order: BN → BLS12 → KSS16 (increasing size/complexity), verifying build+behavior parity after each before moving to the next, so a regression is caught early and stays isolated to one directory.

# Non-goals
* No cross-curve shared library (BN/BLS12 field-tower duplication stays as-is, per your choice).
* No algorithmic changes, no renaming of existing public function names/signatures, no performance changes — this is a structural/file-organization refactor only.
* Not fixing the pre-existing benign warnings/typos already in the code (e.g. duplicate `mpz_mpz_mul` in KSS16's unsigned-long declaration list) unless they block compilation.

# Orchestration
**Decision**: Use a single child agent (per your choice) on Warp's default harness with model `deepseek-v4-pro-fireworks`, working sequentially through BN → BLS12 → KSS16 in the current working tree (no worktrees needed since it's one sequential agent, not parallel ones).
**Dependencies and ordering**: Research (done) → plan approval → agent splits BN, builds+verifies → splits BLS12, builds+verifies → splits KSS16, builds+verifies → final full-repo build check → report back with build/test output for review.
**Launch config**: Local execution, single agent, model `deepseek-v4-pro-fireworks` on the default Oz harness, per the attached orchestration config.
**Child agents**:
* `field-decoupler` — refactoring agent: works directly in `/Users/lamda/Documents/GitHub/pairingma128` on the current branch (no worktree needed, single agent, no concurrent writers); for each curve directory in order (BN, BLS12, KSS16), captures a pre-refactor baseline build+run output, performs the header/source split described above, fixes the global-state `extern`/definition issue, adds a `CMakeLists.txt` for that curve, rebuilds, reruns, and compares pass/fail markers against baseline before moving to the next curve; reports final status (files created, build commands, before/after test markers) back as its completion message.
**Merge strategy**: Single agent, no merge needed; changes land directly in the working tree for your review (diff/PR as you prefer afterward).
# Addendum: TDD mathematical-correctness gate (added 2025-09)

The refactor must be gated by a **mathematical correctness suite**, not just
pass/fail markers. The suites now exist under `tests/` and all pass against the
current monolithic code (baseline established):

| curve | suite file                | checks | status (pre-refactor baseline) |
|-------|---------------------------|--------|-------------------------------|
| BN    | tests/bn_math_test.c      | 16     | ALL PASS (exit 0)              |
| BLS12 | tests/bls12_math_test.c   | 17     | ALL PASS (exit 0)              |
| KSS16 | tests/kss16_math_test.c   | 32     | ALL PASS (exit 0)              |

Run: `tests/run_math_tests.sh [bn|bls12|kss16]` (GMP assumed under
`/opt/homebrew`; override with `MATH_GMP=`). Pre-refactor it compiles each
monolith `.c` with `-Dmain=monolith_main` plus `-fcommon` (the curve headers
hold tentative global definitions; -fcommon merges them across the two TUs).
Verdicts come from **exit codes**, not stdout markers.

## What the suites assert (per curve)

- r is prime and divides #E(Fp); fixed/cleared points satisfy their curve
  equation (`y^2 = x^3 + Ax + B` for BN/BLS12, `y^2 = x^3 + a_x x` for KSS16;
  note: BN stores the curve as `x^3 + Ax - B`, B=4).
- Test points are order r: `[r]P = ∞`, `[r]Q = ∞`. (G1 points are
  cofactor-cleared inside the tests: the original BLS12/BN G1 generators return
  raw rational points and #E(Fp) has a large cofactor there.)
- Pairing value z = e(P,Q): non-degenerate (z ≠ 1), in the order-r subgroup of
  GT (z^r = 1), bilinear (e([a]P,[b]Q) = z^{ab}), alternating in each slot
  (e(-P,Q)e(P,Q) = 1 and e(P,-Q)e(P,Q) = 1), deterministic on identical inputs.
- Field Frobenius identity: π(x) = x^p and π²(x) = π(π(x)).
- KSS16 additionally: all reduced optimal-ate flavors produce the *identical*
  value e(P,Q) (Pseudo_Sparse == Sparse == Optimal) and Ate/Tate pass the full
  self-consistency battery in their own (differently normalized) convention.

## Golden-value parity

BLS12 parameters and test points are deterministic, so its pairing value is
byte-identical run to run. Golden digests live in `tests/golden/` and the runner
fails if a run diverges — this catches subtle pre/post-refactor drift that
identity checks could still miss. (BN points are RNG-generated; KSS16 uses fixed
points; both are stable run-to-run via the identities.)

## Mandatory gate for the refactor (supersedes "pass/fail markers")

1. Before touching a curve directory, run its suite (must be ALL PASS).
2. After splitting a curve, build the **split modules + the test file** (do NOT
   link the reduced `main.c` into the test executable — it now owns `main`;
   rename it with `-Dmain=` or exclude it, as the runner does today). Update
   only the `#include`s in the test file to the new module headers
   (field.h/curve.h/pairing.h/params.h/util.h as needed); the assertions are
   unchanged. Drop `-fcommon` once the headers use extern declarations.
3. Require: exit 0, "ALL PASS", and (BLS12) golden digest MATCHES. Then move to
   the next curve (BN → BLS12 → KSS16, as before).
4. `nm`-symbol parity check (optional but cheap): the union of global symbols
   over the split objects must equal the monolith object's global symbols, to
   catch dropped/misplaced functions and duplicate or missing globals.

## Facts the tests encode (do not "fix" during a structural refactor)

- BN curve equation: y² = x³ + Ax − B (generator subtracts B); BLS12 adds B.
- BLS12 `Tate_pairing` / `Miller_algo_for_tate` are **empty stubs** (not tested).
- BLS12/BN G1 generators do not cofactor-clear; pairing identities hold for raw
  points but the tests use cofactor-cleared G1 for the textbook order claims.
- KSS16 `Fp16_set_ui` writes *every* leaf (not scalar embedding), so the field
  identity is derived in tests, never assumed via set_ui(1).
- KSS16 Ate/Tate drive Miller lines from G2 (literature-ate normalization):
  their values legitimately differ from the optimal-ate family; each is a valid
  pairing in its own right (batteries pass).

# Status: REFACTOR COMPLETE (all three curves)

Executed BN → BLS12 → KSS16 sequentially. Each curve now has the target layout
(field/curve/pairing/params/test/util .h/.c + reduced main.c + CMakeLists.txt).
The old monolith headers (BN12_header.h, BLS_12.h, KSS_16.h) and the BN
monolith source (EFp12_optimal.c) are removed; BLS12/KSS16 main.c are reduced
to `int main()` + orchestration.

## Verification results (TDD gate, all green)

| curve | checks | status | symbol parity | CMake build |
|-------|--------|--------|---------------|-------------|
| BN    | 16/16  | ALL PASS | 309 == 309 | OK (bn + bn_math_test) |
| BLS12 | 17/17  | ALL PASS + golden digest MATCHES | 310 == 310 | OK (bls12 + bls12_math_test) |
| KSS16 | 32/32  | ALL PASS | 387 == 387 | OK (kss16 + kss16_math_test) |

- Symbol parity: `nm` global-symbol union over split objects == monolith object,
  for all three curves (no dropped/duplicated functions or globals).
- Reduced `main()` executables rebuilt from split modules and smoke-tested:
  BN/BLS12 print "test success", KSS16 prints "Success" (bilinearity check),
  matching the original run markers.
- BLS12 pairing value on fixed vectors is byte-identical to the pre-refactor
  golden digest.

## Implementation notes (deviations/details from the original plan)

- The `extern`/definition split for globals: scalar curve params (prime/order/
  trace/etc.) live in params.c; field-level state (basis, frobenius constants,
  ZERO/ONE, KSS16's p_*/p2_*…p8_* tower constants) lives in field.c; diagnostic
  counters live in util.c — each declared extern in its own header. No duplicate
  symbols under -fno-common.
- Headers were extended (not just copied) to declare functions that the
  monolith defined but its header omitted, so they remain global across TUs:
  BN: Fp_mul_basis, Fp2_isCNR, Fp6_isCNR, EFp12_G1_SCM, EFp12_G2_SCM_2div;
  BLS12: Fp_mul_basis, Fp2_isCNR, Fp6_isCNR, EFp_skew_frobenius,
  EFp12_G2_SCM_normal, ff_ltt_vtt, f_ltp_vtp, Miller_algo_for_tate,
  Tate_pairing, test_tate_pairing, test_ate_pairing;
  KSS16: (all tower functions already declared).
- Pre-existing dead declarations (Fp_frobenius_1..10, Fp12_G3_pow,
  EFp12_ECD/ECA_G2optimal, EFp_skew_frobenius_2, Miller_algo_for_ate,
  Ate_pairing, Check_SCM, rational_point_check, etc.) were preserved verbatim
  in the appropriate headers.
- Pre-existing benign quirks preserved as-is: BLS12/KSS16 duplicate
  `mpz_mpz_mul` in the counter declaration, KSS16's `pre_calculate_frob_p8();;`
  double semicolon, BLS12's empty Tate/Miller stubs.
- KSS16 params.c depends on curve.h and util.h (its parameter search uses EFp
  points and timedifference_msec); this is the only cross-layer dependency
  beyond field<-curve<-pairing.
- `tests/run_math_tests.sh` detects split vs monolith by presence of field.c,
  builds split modules + the test TU, and enforces BLS12 golden-digest parity.
