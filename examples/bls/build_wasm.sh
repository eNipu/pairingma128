#!/usr/bin/env bash
# Build docs/bls.js + docs/bls.wasm from the BLS12 modules + GMP (WebAssembly).
# Requires Emscripten (emcc) and a WASM build of GMP (libgmp.a).
# Override with EMSDK=... and GMP_DIR=... if needed.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EMSDK="${EMSDK:-$HOME/emsdk}"
GMP_DIR="${GMP_DIR:-$HOME/gmp-6.3.0}"
cd "$ROOT"

source "$EMSDK/emsdk_env.sh" >/dev/null 2>&1

emcc -O2 \
  -Iexamples/bls -IBLS12 -I"$GMP_DIR" \
  BLS12/field.c BLS12/curve.c BLS12/pairing.c BLS12/params.c BLS12/util.c \
  examples/bls/bls.c examples/bls/sha256.c examples/bls/wasm_export.c \
  "$GMP_DIR/.libs/libgmp.a" \
  -s EXPORTED_FUNCTIONS='["_wasm_bls_init","_wasm_bls_keygen","_wasm_bls_sign","_wasm_bls_verify","_wasm_bls_aggregate","_wasm_bls_aggregate_verify","_malloc","_free"]' \
  -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' \
  -s MODULARIZE=1 \
  -s EXPORT_NAME=createBlsModule \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s STACK_SIZE=1048576 \
  -o docs/bls.js

echo "built docs/bls.js + docs/bls.wasm"
ls -la docs/bls.js docs/bls.wasm
