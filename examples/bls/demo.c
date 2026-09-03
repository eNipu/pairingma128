/* demo.c - native CLI demo of BLS signatures on the BLS12 curve.
 * Prints [ ok ]/[FAIL] per check and exits nonzero on any failure. */
#include "bls.h"
#include <stdio.h>
#include <string.h>

static int checks = 0, fails = 0;

#define CHECK(cond, name)                                                    \
    do {                                                                     \
        checks++;                                                            \
        if (cond) { printf("  [ ok ] %s\n", name); }                         \
        else      { fails++; printf("  [FAIL] %s\n", name); }                \
    } while (0)

/* copy a bls_* result (points into a shared static buffer) */
static void take(char *dst, size_t n, const char *src) {
    if (src && strlen(src) < n) strcpy(dst, src);
    else dst[0] = 0;
}

int main(void) {
    char kp1[32768], kp2[32768];
    char sk1[1024], pk1[8192], pk2[8192];
    char sig1[8192], sig2[8192];
    char *sep;
    const char *m1 = "hello pairing";
    const char *m2 = "world pairing";
    const char *m3 = "evil tampered";
    char sigs[32768], pks[32768], msgs[32768], msgs_bad[32768];

    if (bls_init() != 0) { printf("init failed\n"); return 2; }

    /* keygen: "<sk>|<pk>" */
    take(kp1, sizeof(kp1), bls_keygen());
    sep = strchr(kp1, '|'); if (sep) *sep = 0;
    strcpy(sk1, kp1);
    strcpy(pk1, sep + 1);

    take(sig1, sizeof(sig1), bls_sign(sk1, (const uint8_t *)m1, strlen(m1)));
    take(sig2, sizeof(sig2), bls_sign(sk1, (const uint8_t *)m2, strlen(m2)));

    printf("BLS12 short signatures + aggregation (native)\n");
    printf("  sk = %s\n", sk1);

    /* core checks */
    CHECK(bls_verify(pk1, sig1, (const uint8_t *)m1, strlen(m1)) == 1,
          "valid signature verifies");
    CHECK(bls_verify(pk1, sig1, (const uint8_t *)m3, strlen(m3)) == 0,
          "tampered message is rejected");

    /* a different key must not verify */
    take(kp2, sizeof(kp2), bls_keygen());
    sep = strchr(kp2, '|'); if (sep) *sep = 0;
    strcpy(pk2, sep + 1);
    CHECK(bls_verify(pk2, sig1, (const uint8_t *)m1, strlen(m1)) == 0,
          "wrong public key is rejected");

    /* determinism */
    CHECK(bls_verify(pk1, sig1, (const uint8_t *)m1, strlen(m1)) ==
          bls_verify(pk1, sig1, (const uint8_t *)m1, strlen(m1)),
          "verification is deterministic");

    /* aggregation */
    snprintf(sigs, sizeof(sigs), "%s;%s", sig1, sig2);
    char agg[8192];
    take(agg, sizeof(agg), bls_aggregate(sigs));

    snprintf(pks, sizeof(pks), "%s;%s", pk1, pk1);
    snprintf(msgs, sizeof(msgs), "%s\n%s", m1, m2);
    CHECK(bls_aggregate_verify(agg, pks, msgs) == 1,
          "aggregate signature verifies both messages");

    snprintf(msgs_bad, sizeof(msgs_bad), "%s\n%s", m1, m3);
    CHECK(bls_aggregate_verify(agg, pks, msgs_bad) == 0,
          "aggregate with one tampered message is rejected");

    printf("----------------------------------------\n");
    printf("BLS demo: %d checks, %d failure(s)\n", checks, fails);
    printf("BLS demo: %s\n", fails == 0 ? "ALL PASS" : "FAILED");
    return fails == 0 ? 0 : 1;
}
