/* demo.c - native CLI demo of BLS signatures on the BLS12 curve.
 *
 * The story: a hospital lab signs patient reports with short BLS signatures,
 * then folds two reports into a single aggregate signature. A falsified lab
 * value is caught, and a different lab's key cannot pass as the real one.
 *
 * Prints [ ok ]/[FAIL] per check and exits nonzero on any failure. */
#include "bls.h"
#include <stdio.h>
#include <string.h>

static int checks = 0, fails = 0;

#define CHECK(cond, name)                                                    \
    do {                                                                     \
        checks++;                                                            \
        if (cond) { printf("    [ ok ] %s\n", name); }                       \
        else      { fails++; printf("    [FAIL] %s\n", name); }              \
    } while (0)

/* copy a bls_* result (points into a shared static buffer) */
static void take(char *dst, size_t n, const char *src) {
    if (src && strlen(src) < n) strcpy(dst, src);
    else dst[0] = 0;
}

/* print a long hex value readably: full if short, first n chars + "…" + last 8 */
static void show(const char *label, const char *hex, size_t n) {
    size_t len = strlen(hex);
    if (len <= n + 9) { printf("    %s = %s\n", label, hex); return; }
    printf("    %s = %.*s…%s\n", label, (int)n, hex, hex + len - 8);
}

int main(void) {
    char kp1[32768], kp2[32768];
    char sk1[1024], pk1[8192], pk2[8192];
    char sig1[8192], sig2[8192];
    char *sep;

    /* Patient Ana Rivera's lab reports — the exact bytes that get signed. */
    const char *report1 = "PATIENT: Ana Rivera | TEST: fasting glucose | VALUE: 92 mg/dL | RESULT: NORMAL";
    const char *report2 = "PATIENT: Ana Rivera | TEST: LDL cholesterol | VALUE: 108 mg/dL | RESULT: NORMAL";
    /* Someone edits one lab value before the report reaches the pharmacy. */
    const char *forged  = "PATIENT: Ana Rivera | TEST: fasting glucose | VALUE: 192 mg/dL | RESULT: ABNORMAL";

    char sigs[32768], pks[32768], msgs[32768], msgs_bad[32768];

    if (bls_init() != 0) { printf("init failed\n"); return 2; }

    printf("BLS short signatures — a hospital lab signs patient reports\n");
    printf("============================================================\n");
    printf("St. Mary's lab signs each report with a short BLS signature.\n");
    printf("Anyone with the lab's public key can check a report, and two\n");
    printf("reports' signatures add together into one aggregate signature.\n\n");

    /* keygen: "<sk>|<pk>" */
    take(kp1, sizeof(kp1), bls_keygen());
    sep = strchr(kp1, '|'); if (sep) *sep = 0;
    strcpy(sk1, kp1);
    strcpy(pk1, sep + 1);

    printf("[1] The lab generates its signing keypair\n");
    printf("    sk (secret — never leaves the lab's server) = %s\n", sk1);
    show("pk (public — printed on every report)", pk1, 24);
    printf("\n");

    take(sig1, sizeof(sig1), bls_sign(sk1, (const uint8_t *)report1, strlen(report1)));
    take(sig2, sizeof(sig2), bls_sign(sk1, (const uint8_t *)report2, strlen(report2)));

    printf("[2] The lab signs two of Ana Rivera's reports\n");
    printf("    report 1: %s\n", report1);
    printf("    report 2: %s\n\n", report2);

    printf("[3] The pharmacy verifies reports with the lab's public key\n");
    CHECK(bls_verify(pk1, sig1, (const uint8_t *)report1, strlen(report1)) == 1,
          "the genuine report 1 verifies");
    CHECK(bls_verify(pk1, sig1, (const uint8_t *)forged, strlen(forged)) == 0,
          "a falsified value (92 -> 192 mg/dL) is rejected");

    /* a different lab's key must not verify */
    take(kp2, sizeof(kp2), bls_keygen());
    sep = strchr(kp2, '|'); if (sep) *sep = 0;
    strcpy(pk2, sep + 1);
    CHECK(bls_verify(pk2, sig1, (const uint8_t *)report1, strlen(report1)) == 0,
          "a different lab's public key is rejected");

    /* determinism */
    CHECK(bls_verify(pk1, sig1, (const uint8_t *)report1, strlen(report1)) ==
          bls_verify(pk1, sig1, (const uint8_t *)report1, strlen(report1)),
          "the same report always verifies the same way");

    /* aggregation */
    snprintf(sigs, sizeof(sigs), "%s;%s", sig1, sig2);
    char agg[8192];
    take(agg, sizeof(agg), bls_aggregate(sigs));

    printf("\n[4] The insurer checks both reports with one aggregate signature\n");
    snprintf(pks, sizeof(pks), "%s;%s", pk1, pk1);
    snprintf(msgs, sizeof(msgs), "%s\n%s", report1, report2);
    CHECK(bls_aggregate_verify(agg, pks, msgs) == 1,
          "one aggregate signature verifies both reports");

    snprintf(msgs_bad, sizeof(msgs_bad), "%s\n%s", report1, forged);
    CHECK(bls_aggregate_verify(agg, pks, msgs_bad) == 0,
          "aggregate with one falsified report is rejected");

    printf("\n============================================================\n");
    printf("BLS demo: %d checks, %d failure(s)\n", checks, fails);
    printf("BLS demo: %s\n", fails == 0 ? "ALL PASS" : "FAILED");
    return fails == 0 ? 0 : 1;
}
