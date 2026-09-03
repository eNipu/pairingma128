/* wasm_export.c - Emscripten export surface (called from docs/app.js).
 * Strings are marshalled by cwrap; each function returns into the bls.c
 * static buffer, which JS copies immediately. */
#include <emscripten.h>
#include <string.h>
#include "bls.h"

EMSCRIPTEN_KEEPALIVE
int wasm_bls_init(void) {
    return bls_init();
}

EMSCRIPTEN_KEEPALIVE
char *wasm_bls_keygen(void) {
    return bls_keygen();
}

EMSCRIPTEN_KEEPALIVE
char *wasm_bls_sign(const char *sk_hex, const char *msg) {
    return bls_sign(sk_hex, (const uint8_t *)msg, strlen(msg));
}

EMSCRIPTEN_KEEPALIVE
int wasm_bls_verify(const char *pk_hex, const char *sig_hex, const char *msg) {
    return bls_verify(pk_hex, sig_hex, (const uint8_t *)msg, strlen(msg));
}

EMSCRIPTEN_KEEPALIVE
char *wasm_bls_aggregate(const char *sigs_hex) {
    return bls_aggregate(sigs_hex);
}

EMSCRIPTEN_KEEPALIVE
int wasm_bls_aggregate_verify(const char *agg_hex, const char *pks_hex,
                              const char *msgs_hex) {
    return bls_aggregate_verify(agg_hex, pks_hex, msgs_hex);
}
