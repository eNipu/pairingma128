# pairingma128

Reference C implementation of **optimal-Ate bilinear pairings at the 128-bit
security level** on three pairing-friendly curves — **KSS-16**, **BN**, and
**BLS-12** — accompanied by a BLS short-signature demo that runs natively and
in the browser via WebAssembly.

This repository is the implementation behind the paper:

> **Efficient Optimal Ate Pairing at 128-Bit Security Level**
> Md. Al-Amin Khandaker, Yuki Nanjo, Loubna Ghammam, Sylvain Duquesne,
> Yasuyuki Nogami, Yuta Kodera — *IndoCrypt 2017*.

- PDF (ePrint): <https://eprint.iacr.org/2017/1174.pdf>
- Springer: <https://link.springer.com/chapter/10.1007/978-3-319-71667-1_10>

```bibtex
@inproceedings{khandaker2017efficient,
  title     = {Efficient Optimal Ate Pairing at 128-Bit Security Level},
  author    = {Khandaker, Md. Al-Amin and Nanjo, Yuki and Ghammam, Loubna and
               Duquesne, Sylvain and Nogami, Yasuyuki and Kodera, Yuta},
  booktitle = {Progress in Cryptology -- INDOCRYPT 2017},
  series    = {Lecture Notes in Computer Science},
  volume    = {10698},
  pages     = {186--205},
  publisher = {Springer},
  year      = {2017},
  doi       = {10.1007/978-3-319-71667-1_10}
}
```

## What the paper is about

Pairing-friendly curves for the 128-bit security level come in two families:
BN and BLS-12 use a **sextic twist**, while KSS-16 is the atypical, less
studied choice with a **quartic twist**. The paper optimises Miller's
algorithm for the optimal-Ate pairing on KSS-16 by deriving efficient *sparse*
multiplication, and implements the result alongside BN and BLS-12 to
experimentally verify Barbulescu et al.'s estimation. The conclusion is that
the derived **pseudo 8-sparse multiplication** makes Miller's algorithm on
KSS-16 the most efficient of the three at the 128-bit level, supporting
Barbulescu–Duquesne's recommendation.

## Curves

| curve | embedding degree | twist    | field tower                       | pairing entry points                     |
|-------|------------------|----------|-----------------------------------|------------------------------------------|
| BN    | 12               | sextic   | Fp → Fp² → Fp⁶ → Fp¹²             | optimal-Ate, Ate, Tate                    |
| BLS12 | 12               | sextic   | Fp → Fp² → Fp⁶ → Fp¹²             | optimal-Ate, Ate, Tate (Tate is a stub)   |
| KSS16 | 16               | quartic  | Fp → Fp² → Fp⁴ → Fp⁸ → Fp¹⁶       | optimal-Ate, Ate, Tate, plus Pseudo/Sparse variants |

Each curve directory is split into the same module layout
(`field` / `curve` / `pairing` / `params` / `test` / `util` `.h` + `.c`, plus a
reduced `main.c`), decoupling finite-field arithmetic from curve arithmetic and
the Miller-loop/pairing logic.

## Repository layout

```
BN/                  BN curve (sextic twist), field/curve/pairing/params/test/util
BLS12/               BLS12 curve (sextic twist), same module layout
KSS16/               KSS16 curve (quartic twist), same module layout
tests/               mathematical-correctness suites + golden digest
examples/bls/        BLS short signatures + aggregation (native CLI & WASM)
docs/                static BLS demo site (index.html, app.js, bls.js, bls.wasm)
Decouple finite field arithmetic in pairingma128.md   refactor design notes
```

## Requirements

- A C compiler (gcc/clang)
- [GMP](https://gmplib.org/) (`brew install gmp` on macOS)
- CMake ≥ 3.10
- (WebAssembly rebuild only) [Emscripten](https://emscripten.org/) and a WASM
  build of GMP

## Build

Each curve directory has its own `CMakeLists.txt` producing a self-check
executable and a math-test executable:

```sh
cmake -S BLS12 -B build/bls12 && cmake --build build/bls12
./build/bls12/bls12          # prints curve parameters, then runs a pairing check
```

The same works for `BN/` (targets `bn`, `bn_math_test`) and `KSS16/` (targets
`kss16`, `kss16_math_test`). `main()` prints the curve parameters and runs the
curve's own pairing self-check (`"test success"` / `"Success"`).

## Mathematical correctness tests

`tests/` gates the implementation with pairing identities, not just
pass/fail markers. Verdicts come from exit codes (0 = all checks passed).

```sh
bash tests/run_math_tests.sh            # all curves
bash tests/run_math_tests.sh bn bls12   # subset
```

| curve | suite                | checks |
|-------|----------------------|--------|
| BN    | `tests/bn_math_test.c`    | 16 |
| BLS12 | `tests/bls12_math_test.c` | 17 (+ byte-exact golden digest) |
| KSS16 | `tests/kss16_math_test.c` | 32 |

Asserted properties include: prime order `r` dividing `#E(Fp)`; test points on
the curve and of exact order `r`; pairing output non-degenerate, in the order-r
subgroup, bilinear, alternating, and deterministic; field Frobenius identity
`π(x) = x^p`; and on KSS16 that the Pseudo/Sparse/Optimal flavours produce the
identical value. The BLS12 golden digest in `tests/golden/` pins the pairing
value byte-for-byte across runs and refactors. See `tests/README.md`.

## BLS signatures demo

`examples/bls/` implements Boneh–Lynn–Shacham short signatures and signature
aggregation on the BLS12 curve, on top of the pairing library:

```
pk = sk·Q ∈ G₂ ,   σ = sk·H(m) ∈ G₁ ,   verify: e(σ, Q) == e(H(m), pk)
```

Two frontends share the same core (`bls.c`): a **native CLI** (`demo.c`) and a
**WebAssembly** build (`wasm_export.c`) that runs the real 128-bit optimal-Ate
pairing in the browser.

```sh
cmake -S examples/bls -B build/bls && cmake --build build/bls
./build/bls/bls_demo         # "BLS demo: ALL PASS", exit 0
```

To rebuild the browser assets (`docs/bls.js`, `docs/bls.wasm`):

```sh
EMSDK=$HOME/emsdk GMP_DIR=$HOME/gmp-6.3.0 ./examples/bls/build_wasm.sh
```

A fixed, order-r G2 generator is hardcoded in `bls_init()` (validated by
`tests/bls12_math_test.c`) to avoid the library's impractically slow
`EFp12_generate_G2()` in the software-GMP WASM build. See `examples/bls/README.md`.

## Live demo

`docs/` is a self-contained static site. It is vendored into the portfolio repo
(`eNipu/eNipu.github.io`) as `public/bls/` and served at
**<https://al-am.in/bls/>** by Cloudflare Workers Builds — no separate Pages
project, credentials, or DNS record. See `examples/bls/DEPLOY.md`.

## Status

All three curves are refactored into decoupled modules, all math suites pass
(BN 16/16, BLS12 17/17 + golden digest, KSS16 32/32), and `nm` symbol-parity
checks confirm the split objects expose the same globals as the original
monoliths. This is research/reference code; it has not been audited and is not
intended for production use.
