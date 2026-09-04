# BLS short signatures (+ aggregation)

Boneh–Lynn–Shacham signatures on the BLS12 curve, built on the refactored
pairing library. Two frontends share the same core (`bls.c`):

- **Native CLI** (`demo.c`) — prints `[ ok ]/[FAIL]` checks, exit code based.
- **WebAssembly** (`wasm_export.c`) — the browser runs the *actual* 128-bit
  Optimal-Ate pairing via the `docs/` GitHub Pages site.

## The math (one equation)

`pk = sk·Q ∈ G₂`, `σ = sk·H(m) ∈ G₁`, verify: `e(σ, Q) == e(H(m), pk)`.

Aggregation is just `σ = σ₁ + σ₂ + …`; verification becomes
`e(σ, Q) == e(H(m₁),pk₁) · e(H(m₂),pk₂) · …`.

## Layout

| file | purpose |
|------|---------|
| `sha256.c/h` | compact SHA-256 for hash-to-curve |
| `bls.c/h`    | keygen / sign / verify / aggregate + hex serialization |
| `wasm_export.c` | Emscripten `EMSCRIPTEN_KEEPALIVE` surface |
| `demo.c`     | native CLI demo (6 TDD checks) |
| `build_wasm.sh` | build `docs/bls.js` + `docs/bls.wasm` |

## Native build

```sh
cmake -S examples/bls -B build/bls && cmake --build build/bls
./build/bls/bls_demo        # expect "BLS demo: ALL PASS", exit 0
```

## WebAssembly build (→ GitHub Pages)

Prereqs: Emscripten (`emsdk`) and a WASM build of GMP (`libgmp.a`).

```sh
# one-time GMP build
emconfigure ./configure --host=none --disable-assembly --disable-shared \
    --enable-static CC_FOR_BUILD=cc CFLAGS="-O2"
emmake make          # produces .libs/libgmp.a

# build the site assets
EMSDK=$HOME/emsdk GMP_DIR=$HOME/gmp-6.3.0 ./examples/bls/build_wasm.sh
# -> docs/bls.js, docs/bls.wasm
```

`docs/` is a self-contained static site (`index.html`, `style.css`, `app.js`,
`bls.js`, `bls.wasm`). Serve it or point GitHub Pages at it.

### Deploy

The demo ships as part of the portfolio site (`eNipu/eNipu.github.io`), which
Cloudflare Workers Builds deploys on every push to `master`. `docs/` is vendored
there as `public/bls/` and served at <https://al-am.in/bls/>.

To publish a rebuilt `docs/`:

```sh
rsync -a docs/ ../eNipu.github.io/public/bls/   # then commit + push there
```

`docs/` uses only relative URLs, so it also works from any other static host
(`python3 -m http.server` in `docs/`, GitHub Pages, etc.). See `DEPLOY.md`.

## Notes

- A fixed, order-r G2 generator is hardcoded in `bls_init()` (standard practice;
  it is validated by `tests/bls12_math_test.c`). This avoids the library's
  `EFp12_generate_G2()`, which does a ~4600-bit scalar multiplication that is
  impractically slow in the software-GMP WASM build.
- `bls.c` targets the BLS12 curve; the BN/KSS16 modules expose the same API, so
  porting is a matter of swapping the header include and the generator.
