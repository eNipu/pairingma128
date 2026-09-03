/* bls.c - BLS short signatures on the BLS12 curve, built on the refactored
 * pairing library modules (field/curve/pairing/params). */
#include "bls.h"
#include "sha256.h"

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "field.h"
#include "curve.h"
#include "pairing.h"
#include "params.h"

static struct EFp12 g_Q;   /* fixed G2 generator */
static int g_ready = 0;
static char g_out[65536];
static unsigned long g_counter = 0;

/* ------------------------------------------------------------------ */
/* helpers                                                             */
/* ------------------------------------------------------------------ */

static void random_sk(mpz_t sk) {
    gmp_randstate_t st;
    gmp_randinit_default(st);
    gmp_randseed_ui(st, (unsigned long)time(NULL) + (g_counter++ * 2654435761UL));
    mpz_urandomm(sk, st, EFp_order);
    if (mpz_sgn(sk) == 0) mpz_set_ui(sk, 1);   /* avoid the zero key */
    gmp_randclear(st);
}

/* try-and-increment hash-to-G1: SHA256(msg || ctr) -> x -> y = sqrt(x^3+Ax+B),
 * then cofactor-clear into the order-r subgroup. */
static void hash_to_G1(struct EFp12 *out, const uint8_t *msg, size_t msglen) {
    struct Fp x, y, x2, x3, ax, rhs;
    uint8_t hash[32], ctrbuf[4];
    uint32_t ctr = 0;
    mpz_t t, h1;
    SHA256_CTX sctx;

    Fp_init(&x); Fp_init(&y); Fp_init(&x2); Fp_init(&x3);
    Fp_init(&ax); Fp_init(&rhs);
    mpz_init(t); mpz_init(h1);

    for (;; ctr++) {
        ctrbuf[0] = (uint8_t)(ctr >> 24); ctrbuf[1] = (uint8_t)(ctr >> 16);
        ctrbuf[2] = (uint8_t)(ctr >> 8);  ctrbuf[3] = (uint8_t)ctr;
        sha256_init(&sctx);
        sha256_update(&sctx, msg, msglen);
        sha256_update(&sctx, ctrbuf, 4);
        sha256_final(&sctx, hash);

        mpz_import(t, 32, 1, 1, 0, 0, hash);
        mpz_mod(t, t, prime);
        Fp_set_mpz(&x, t);

        Fp_mul(&x2, &x, &x);
        Fp_mul(&x3, &x2, &x);
        Fp_mul_mpz(&ax, &x, curve_parameter_A);
        Fp_add(&rhs, &x3, &ax);
        Fp_add_mpz(&rhs, &rhs, curve_parameter_B);
        if (Fp_legendre(&rhs) == 1) {
            Fp_sqrt(&y, &rhs);
            break;
        }
    }

    EFp12_init(out);
    EFp12_set_ui(out, 0);
    Fp_set(&out->x.x0.x0.x0, &x);
    Fp_set(&out->y.x0.x0.x0, &y);
    out->flag = 0;

    mpz_tdiv_q(h1, EFp_total, EFp_order);  /* #E(Fp)/r, the G1 cofactor */
    {
        struct EFp12 P;
        EFp12_init(&P);
        EFp12_set(&P, out);
        EFp12_SCM(out, &P, h1);
        EFp12_clear(&P);
    }

    Fp_clear(&x); Fp_clear(&y); Fp_clear(&x2); Fp_clear(&x3);
    Fp_clear(&ax); Fp_clear(&rhs);
    mpz_clear(t); mpz_clear(h1);
}

/* G1 point: base-field x,y */
static void serialize_g1(char *buf, const struct EFp12 *P) {
    gmp_sprintf(buf, "%Zx:%Zx", P->x.x0.x0.x0.x0, P->y.x0.x0.x0.x0);
}

/* G2 point: twist Fp2 coords */
static void serialize_g2(char *buf, const struct EFp12 *Q) {
    gmp_sprintf(buf, "%Zx:%Zx:%Zx:%Zx",
                Q->x.x0.x2.x0.x0, Q->x.x0.x2.x1.x0,
                Q->y.x1.x1.x0.x0, Q->y.x1.x1.x1.x0);
}

static int parse_g1(struct EFp12 *P, const char *s) {
    char x[512], y[512];
    const char *colon = strchr(s, ':');
    size_t n;
    if (!colon) return -1;
    n = (size_t)(colon - s);
    if (n >= sizeof(x)) return -1;
    memcpy(x, s, n); x[n] = 0;
    if (strlen(colon + 1) >= sizeof(y)) return -1;
    strcpy(y, colon + 1);
    EFp12_init(P);
    EFp12_set_ui(P, 0);
    P->flag = 0;
    if (mpz_set_str(P->x.x0.x0.x0.x0, x, 16) != 0) return -1;
    if (mpz_set_str(P->y.x0.x0.x0.x0, y, 16) != 0) return -1;
    return 0;
}

static int parse_g2(struct EFp12 *Q, const char *s) {
    char f[4][512];
    const char *p = s;
    int i;
    for (i = 0; i < 4; i++) {
        const char *colon = strchr(p, ':');
        size_t n = colon ? (size_t)(colon - p) : strlen(p);
        if (n >= sizeof(f[i])) return -1;
        memcpy(f[i], p, n); f[i][n] = 0;
        if (i < 3) { if (!colon) return -1; p = colon + 1; }
    }
    EFp12_init(Q);
    EFp12_set_ui(Q, 0);
    Q->flag = 0;
    if (mpz_set_str(Q->x.x0.x2.x0.x0, f[0], 16) != 0) return -1;
    if (mpz_set_str(Q->x.x0.x2.x1.x0, f[1], 16) != 0) return -1;
    if (mpz_set_str(Q->y.x1.x1.x0.x0, f[2], 16) != 0) return -1;
    if (mpz_set_str(Q->y.x1.x1.x1.x0, f[3], 16) != 0) return -1;
    return 0;
}

/* ------------------------------------------------------------------ */
/* public API                                                          */
/* ------------------------------------------------------------------ */

int bls_init(void) {
    if (!g_ready) {
        init_parameters();
        set_parameters();

        /* Fixed G2 generator (order r, verified by the math-correctness suite).
         * We hardcode it instead of calling EFp12_generate_G2(), which performs
         * a ~4600-bit scalar multiplication that is far too slow in the
         * software-GMP WebAssembly build. Fixed generators are standard
         * practice (e.g. BLS12-381 uses fixed points). */
        EFp12_init(&g_Q);
        EFp12_set_ui(&g_Q, 0);
        g_Q.flag = 0;
        mpz_set_str(g_Q.x.x0.x2.x0.x0,
            "1716096407462707739249584854219035887699085689300518488744116638769194609975441542417090151676720383761847673487573410049924094001307902771", 10);
        mpz_set_str(g_Q.x.x0.x2.x1.x0,
            "1884774744539223134856244661850811554160388946147660758387803399777014198117917660870489415649762044750981220496189685117651039305687807588", 10);
        mpz_set_str(g_Q.y.x1.x1.x0.x0,
            "862911509351712110897605125200443999396611767860495548282350223126731430707014429484782130102632482076556245081709440912377961089052176300", 10);
        mpz_set_str(g_Q.y.x1.x1.x1.x0,
            "14917642869256953479675328417653271711757572421578255369823731891492952725387400103789740069957122839470285312924022325651200303158794609", 10);
        g_ready = 1;
    }
    return 0;
}

char *bls_keygen(void) {
    mpz_t sk;
    struct EFp12 pk;
    char skbuf[512], pkbuf[4096];
    mpz_init(sk);
    random_sk(sk);
    EFp12_init(&pk);
    EFp12_SCM(&pk, &g_Q, sk);
    gmp_sprintf(skbuf, "%Zx", sk);
    serialize_g2(pkbuf, &pk);
    gmp_sprintf(g_out, "%s|%s", skbuf, pkbuf);
    mpz_clear(sk);
    EFp12_clear(&pk);
    return g_out;
}

char *bls_sign(const char *sk_hex, const uint8_t *msg, size_t msglen) {
    mpz_t sk;
    struct EFp12 H, sig;
    mpz_init(sk);
    mpz_set_str(sk, sk_hex, 16);
    hash_to_G1(&H, msg, msglen);
    EFp12_init(&sig);
    EFp12_SCM(&sig, &H, sk);
    serialize_g1(g_out, &sig);
    mpz_clear(sk);
    EFp12_clear(&H);
    EFp12_clear(&sig);
    return g_out;
}

int bls_verify(const char *pk_hex, const char *sig_hex,
               const uint8_t *msg, size_t msglen) {
    struct EFp12 pk, sig, H;
    struct Fp12 e1, e2;
    int ok;

    if (parse_g2(&pk, pk_hex) != 0) return 0;
    if (parse_g1(&sig, sig_hex) != 0) return 0;
    hash_to_G1(&H, msg, msglen);

    Fp12_init(&e1);
    Fp12_init(&e2);
    Opt_ate_pairing(&e1, &g_Q, &sig);   /* e(sig, Q) */
    Opt_ate_pairing(&e2, &pk,  &H);     /* e(H, pk)  */
    ok = (Fp12_cmp(&e1, &e2) == 0);

    Fp12_clear(&e1);
    Fp12_clear(&e2);
    EFp12_clear(&pk);
    EFp12_clear(&sig);
    EFp12_clear(&H);
    return ok;
}

char *bls_aggregate(const char *sigs_hex) {
    struct EFp12 agg, s;
    int first = 1;
    const char *p = sigs_hex;
    EFp12_init(&agg);
    EFp12_init(&s);

    while (*p) {
        const char *semi = strchr(p, ';');
        size_t n = semi ? (size_t)(semi - p) : strlen(p);
        char tok[8192];
        if (n >= sizeof(tok)) n = sizeof(tok) - 1;
        memcpy(tok, p, n); tok[n] = 0;
        if (parse_g1(&s, tok) == 0) {
            if (first) { EFp12_set(&agg, &s); first = 0; }
            else EFp12_ECA(&agg, &agg, &s);
        }
        if (!semi) break;
        p = semi + 1;
    }

    serialize_g1(g_out, &agg);
    EFp12_clear(&agg);
    EFp12_clear(&s);
    return g_out;
}

int bls_aggregate_verify(const char *agg_hex, const char *pks_hex,
                         const char *msgs_hex) {
    struct EFp12 agg, pk, H;
    struct Fp12 e1, e2, prod;
    int first = 1, ok = 0;
    const char *pp = pks_hex, *mp = msgs_hex;

    if (parse_g1(&agg, agg_hex) != 0) return 0;
    Fp12_init(&e1); Fp12_init(&e2); Fp12_init(&prod);
    Opt_ate_pairing(&e1, &g_Q, &agg);   /* e(agg, Q) */

    while (*pp && *mp) {
        const char *semi = strchr(pp, ';');
        const char *nl = strchr(mp, '\n');
        size_t pn = semi ? (size_t)(semi - pp) : strlen(pp);
        size_t mn = nl ? (size_t)(nl - mp) : strlen(mp);
        char pktok[8192], msgtok[8192];
        if (pn >= sizeof(pktok)) pn = sizeof(pktok) - 1;
        if (mn >= sizeof(msgtok)) mn = sizeof(msgtok) - 1;
        memcpy(pktok, pp, pn); pktok[pn] = 0;
        memcpy(msgtok, mp, mn); msgtok[mn] = 0;

        if (parse_g2(&pk, pktok) == 0) {
            hash_to_G1(&H, (const uint8_t *)msgtok, mn);
            Opt_ate_pairing(&e2, &pk, &H);   /* e(H_i, pk_i) */
            if (first) { Fp12_set(&prod, &e2); first = 0; }
            else Fp12_mul(&prod, &prod, &e2);
            EFp12_clear(&pk);
            EFp12_clear(&H);
        }
        if (!semi) break;
        pp = semi + 1;
        if (nl) mp = nl + 1; else break;
    }

    ok = (!first) && (Fp12_cmp(&e1, &prod) == 0);
    Fp12_clear(&e1); Fp12_clear(&e2); Fp12_clear(&prod);
    EFp12_clear(&agg);
    return ok;
}
