/* kss16_math_test.c - mathematical correctness suite for the KSS16
 * optimal-ate pairing (reference: eprint 2017/1174).
 *
 * TDD gate: every CHECK below is a mathematical identity that must hold for a
 * correct implementation. Failure => nonzero exit.
 *
 * Battery runs on Pseudo_Sparse_Optimal_Ate_Pairing with the same fixed
 * G1/G2 points as the code's own check_Pairing (deterministic, so pairing
 * values are run-to-run reproducible). Cross-entry consistency and the
 * Frobenius-map semantics are printed as probes first; promote to CHECK only
 * once their relationship is established.
 */
#include <stdio.h>
#include <gmp.h>
#include "math_assert.h"
#include "field.h"
#include "curve.h"
#include "pairing.h"
#include "params.h"

static void fp16_hex_leaf(FILE *f, struct Fp *p){ gmp_fprintf(f, "%Zx", p->x0); }

/* --- point/order helpers ------------------------------------------------ */
static int is_inf_EFp(const struct EFp *P){ return P->infity == TRUE; }
static int is_inf_EFp16(const struct EFp16 *P){ return P->infity == TRUE; }

/* Full pairing-identity battery for an entry point with signature
 * pair(ANS, X, Y) where X, Y are EFp16 points of order r. All assertions are
 * convention-agnostic (bilinearity is tested per slot). */
static void battery16(const char *tag,
                      void (*pair)(struct Fp16 *, struct EFp16 *, struct EFp16 *),
                      const struct EFp16 *X, const struct EFp16 *Y,
                      mpz_t r, mpz_t a, mpz_t b, mpz_t ab, mpz_t rminus1){
    struct Fp16 z, zr, one, oz, zp, z2, zt, e_mX, e_mY, prod;
    struct EFp16 aX, bY, mX, mY;
    Fp16_init(&z); Fp16_init(&zr); Fp16_init(&one); Fp16_init(&oz);
    Fp16_init(&zp); Fp16_init(&z2); Fp16_init(&zt);
    Fp16_init(&e_mX); Fp16_init(&e_mY); Fp16_init(&prod);
    EFp16_init(&aX); EFp16_init(&bY); EFp16_init(&mX); EFp16_init(&mY);

    pair(&z, (struct EFp16 *)X, (struct EFp16 *)Y);
    Fp16_div(&one, &z, &z);
    Fp16_mul(&oz, &one, &z);
    CHECK(Fp16_cmp(&oz, &z) == 0, "[b16] one*z == z sanity");
    CHECK(Fp16_cmp(&z, &one) != 0, "[b16] non-degenerate (z != 1)");
    Fp16_pow(&zp, &z, r);
    CHECK(Fp16_cmp(&zp, &one) == 0, "[b16] z^r == 1 (subgroup)");
    EFp16_SCM_BIN(&aX, (struct EFp16 *)X, a);
    EFp16_SCM_BIN(&bY, (struct EFp16 *)Y, b);
    pair(&z2, &aX, &bY);
    Fp16_pow(&zt, &z, ab);
    CHECK(Fp16_cmp(&z2, &zt) == 0, "[b16] bilinear in both slots ([a]X,[b]Y)");
    EFp16_SCM_BIN(&mX, (struct EFp16 *)X, rminus1);
    pair(&e_mX, &mX, (struct EFp16 *)Y);
    Fp16_mul(&prod, &e_mX, &z);
    CHECK(Fp16_cmp(&prod, &one) == 0, "[b16] alternating in slot 1");
    EFp16_SCM_BIN(&mY, (struct EFp16 *)Y, rminus1);
    pair(&e_mY, (struct EFp16 *)X, &mY);
    Fp16_mul(&prod, &e_mY, &z);
    CHECK(Fp16_cmp(&prod, &one) == 0, "[b16] alternating in slot 2");
    pair(&zr, (struct EFp16 *)X, (struct EFp16 *)Y);
    CHECK(Fp16_cmp(&zr, &z) == 0, "[b16] deterministic on identical inputs");

    Fp16_clear(&z); Fp16_clear(&zr); Fp16_clear(&one); Fp16_clear(&oz);
    Fp16_clear(&zp); Fp16_clear(&z2); Fp16_clear(&zt);
    Fp16_clear(&e_mX); Fp16_clear(&e_mY); Fp16_clear(&prod);
    EFp16_clear(&aX); EFp16_clear(&bY); EFp16_clear(&mX); EFp16_clear(&mY);
    (void)tag;
}

int main(void){
    mpz_t a, b, ab, rminus1;

    /* same initialization sequence as the original main() */
    mpz_init(X);
    generate_X();
    mpz_init(PRIME_P);
    mpz_init(order_r);
    mpz_init(order_EFp);
    mpz_init(trace_t);
    mpz_init(a_x);
    KSS_16_parameters();
    pre_calculate_frob_p();
    pre_calculate_frob_p2();
    pre_calculate_frob_p3();
    pre_calculate_frob_p4();
    pre_calculate_frob_p5();
    pre_calculate_frob_p6();
    pre_calculate_frob_p7();
    pre_calculate_frob_p8();
    pre_calc_vector_final_exp();

    gmp_printf("-- KSS16: p bits=%d r bits=%d #E(Fp) bits=%d --\n",
               (int)mpz_sizeinbase(PRIME_P,2), (int)mpz_sizeinbase(order_r,2),
               (int)mpz_sizeinbase(order_EFp,2));

    /* ---------- setup sanity ---------- */
    CHECK(mpz_probab_prime_p(order_r, 40) >= 1, "r is (probable) prime");
    {
        mpz_t rem; mpz_init(rem);
        mpz_mod(rem, order_EFp, order_r);
        CHECK(mpz_sgn(rem) == 0, "r divides #E(Fp)");
        mpz_clear(rem);
    }

    /* ---------- fixed points (same constants as check_Pairing) ---------- */
    struct EFp R;          /* fixed point on E(Fp) */
    struct EFp16 Q;        /* fixed point on the twist (G2) */
    mpz_t x, y;
    EFp_init(&R);
    EFp16_init(&Q);
    mpz_inits(x, y, (mpz_ptr)NULL);
    mpz_set_str(x, "585634432000126707887057629201458798521445673169852167601690963969803976546284829414253996667593644070", 10);
    mpz_set_str(y, "354142610898165644039571248952651466071208533406755516950532183093299777682028181848353083070437124285", 10);
    Fp_set_mpz(&R.x, x);
    Fp_set_mpz(&R.y, y);
    {
        mpz_t qx1, qx2, qx3, qx4, qy1, qy2, qy3, qy4;
        mpz_inits(qx1, qx2, qx3, qx4, qy1, qy2, qy3, qy4, (mpz_ptr)NULL);
        mpz_set_str(qx1, "244491741495785299367612042330502507448971522000836606635913541960426621971155507113603357863094662406", 10);
        mpz_set_str(qx2, "218139782565182912336439776658495825407576298090161716969444130503315739115994708740875588793005720529", 10);
        mpz_set_str(qx3, "489284606975189850057489668560782153703917092806842888347849447445917065112045634272024913067442383327", 10);
        mpz_set_str(qx4, "246161136867400494929650830056723287620786052741538082500453356439354042813046547782023695596667138627", 10);
        mpz_set_str(qy1, "48629828925217256502074899893821067838880538260652902679139733711894991914937403050096788961116072841", 10);
        mpz_set_str(qy2, "338343177397929794197126001730091007944155731601968256694856657191699538746843500939799293603025576582", 10);
        mpz_set_str(qy3, "38952231663050097715674495759233974069804273632029475897530255418880509349197207272848037344062136227", 10);
        mpz_set_str(qy4, "250604271284793810746824678011335939079727656816914850518650242122895946469531880085329357599729742073", 10);
        Fp_set_mpz(&Q.x.x0.x1.x0.x0, qx1);
        Fp_set_mpz(&Q.x.x0.x1.x0.x1, qx2);
        Fp_set_mpz(&Q.x.x0.x1.x1.x0, qx3);
        Fp_set_mpz(&Q.x.x0.x1.x1.x1, qx4);
        Fp_set_mpz(&Q.y.x1.x1.x0.x0, qy1);
        Fp_set_mpz(&Q.y.x1.x1.x0.x1, qy2);
        Fp_set_mpz(&Q.y.x1.x1.x1.x0, qy3);
        Fp_set_mpz(&Q.y.x1.x1.x1.x1, qy4);
        mpz_clears(qx1, qx2, qx3, qx4, qy1, qy2, qy3, qy4, (mpz_ptr)NULL);
    }

    /* fixed R on E(Fp): y^2 == x^3 + a_x x */
    {
        struct Fp x2, x3, axx, rhs, y2;
        Fp_init(&x2); Fp_init(&x3); Fp_init(&axx);
        Fp_init(&rhs); Fp_init(&y2);
        Fp_mul(&y2, &R.y, &R.y);
        Fp_sqr(&x2, &R.x);
        Fp_mul(&x3, &x2, &R.x);
        Fp_mul_mpz(&axx, &R.x, a_x);
        Fp_add(&rhs, &x3, &axx);
        CHECK(mpz_cmp(y2.x0, rhs.x0) == 0, "fixed R is on E(Fp): y^2 == x^3 + a_x x");
        Fp_clear(&x2); Fp_clear(&x3); Fp_clear(&axx);
        Fp_clear(&rhs); Fp_clear(&y2);
    }
    CHECK(!is_inf_EFp(&R), "R is not the point at infinity");
    CHECK(!is_inf_EFp16(&Q), "Q is not the point at infinity");

    /* cofactor-clear R into G1 = E(Fp)[r] */
    struct EFp P;
    mpz_t h1;
    EFp_init(&P);
    mpz_init(h1);
    mpz_tdiv_q(h1, order_EFp, order_r);
    EFp_SCM_BIN(&P, &R, h1);
    mpz_clear(h1);
    CHECK(!is_inf_EFp(&P), "P = [h1]R (cofactor-cleared) is not infinity");

    mpz_init(a); mpz_init(b); mpz_init(ab); mpz_init(rminus1);
    mpz_set_ui(a, 31); mpz_set_ui(b, 11); mpz_mul(ab, a, b);
    mpz_sub_ui(rminus1, order_r, 1);

    /* ---------- point orders: [r]P = O, [r]Q = O ---------- */
    {
        struct EFp T1; struct EFp16 T2;
        EFp_init(&T1); EFp16_init(&T2);
        EFp_SCM_BIN(&T1, &P, order_r);
        CHECK(is_inf_EFp(&T1), "[r]P == infinity (P in G1 has order r)");
        EFp16_SCM_BIN(&T2, &Q, order_r);
        CHECK(is_inf_EFp16(&T2), "[r]Q == infinity (Q in G2 has order r)");
        EFp_clear(&T1); EFp16_clear(&T2);
    }

    /* ---------- core pairing identities on Pseudo_Sparse_Optimal_Ate ---------- */
    {
        struct Fp16 z, zr, one, oz, zp, z2, zt, e_mP, e_mQ, prod;
        struct EFp aP, mP;
        struct EFp16 bQ, mQ;
        Fp16_init(&z); Fp16_init(&zr); Fp16_init(&one); Fp16_init(&oz);
        Fp16_init(&zp); Fp16_init(&z2); Fp16_init(&zt);
        Fp16_init(&e_mP); Fp16_init(&e_mQ); Fp16_init(&prod);
        EFp_init(&aP); EFp_init(&mP);
        EFp16_init(&bQ); EFp16_init(&mQ);

        Pseudo_Sparse_Optimal_Ate_Pairing(&z, &P, &Q);

        /* identity derived from z itself (Fp16_set_ui(1) is NOT the field
         * identity in this tower: set_ui fills every leaf) */
        Fp16_div(&one, &z, &z);
        Fp16_mul(&oz, &one, &z);
        CHECK(Fp16_cmp(&oz, &z) == 0, "GT identity sanity (one*z == z)");
        CHECK(Fp16_cmp(&z, &one) != 0, "e(P,Q) != 1 (non-degenerate)");

        Fp16_pow(&zp, &z, order_r);
        CHECK(Fp16_cmp(&zp, &one) == 0, "e(P,Q)^r == 1 (subgroup)");

        /* bilinearity */
        EFp_SCM_BIN(&aP, &P, a);
        EFp16_SCM_BIN(&bQ, &Q, b);
        Pseudo_Sparse_Optimal_Ate_Pairing(&z2, &aP, &bQ);
        Fp16_pow(&zt, &z, ab);
        CHECK(Fp16_cmp(&z2, &zt) == 0, "e([a]P,[b]Q) == e(P,Q)^(ab) (bilinear, a=31 b=11)");

        /* alternating in G1 */
        EFp_SCM_BIN(&mP, &P, rminus1);
        Pseudo_Sparse_Optimal_Ate_Pairing(&e_mP, &mP, &Q);
        Fp16_mul(&prod, &e_mP, &z);
        CHECK(Fp16_cmp(&prod, &one) == 0, "e(-P,Q)*e(P,Q) == 1 (alternating in G1)");

        /* alternating in G2 */
        EFp16_SCM_BIN(&mQ, &Q, rminus1);
        Pseudo_Sparse_Optimal_Ate_Pairing(&e_mQ, &P, &mQ);
        Fp16_mul(&prod, &e_mQ, &z);
        CHECK(Fp16_cmp(&prod, &one) == 0, "e(P,-Q)*e(P,Q) == 1 (alternating in G2)");

        /* reproducibility */
        Pseudo_Sparse_Optimal_Ate_Pairing(&zr, &P, &Q);
        CHECK(Fp16_cmp(&zr, &z) == 0, "pairing is deterministic on identical inputs");

        /* ---- cross-entry consistency (all reduced optimal-ate flavors must
         *      produce the identical value e(P,Q)) ---- */
        {
            struct Fp16 e_sp, e_opt;
            struct EFp16 P16;
            Fp16_init(&e_sp); Fp16_init(&e_opt); EFp16_init(&P16);
            EFp16_set_EFp(&P16, &P);
            Sparse_Optimal_Ate_Pairing(&e_sp, &P, &Q);
            Optimal_Ate_Pairing(&e_opt, &P16, &Q);
            CHECK(Fp16_cmp(&e_sp, &z) == 0,
                  "Sparse_Optimal_Ate == Pseudo_Sparse (same pairing value)");
            CHECK(Fp16_cmp(&e_opt, &z) == 0,
                  "Optimal_Ate == Pseudo_Sparse (same pairing value)");
            Fp16_clear(&e_sp); Fp16_clear(&e_opt); EFp16_clear(&P16);
        }

        /* Ate_Pairing / Tate_Pairing use a different internal normalization
         * (they drive the Miller lines from G2 instead of G1), so their VALUE
         * is not expected to equal e(P,Q); instead verify they are each a
         * valid bilinear alternating pairing in their own convention. */
        {
            struct EFp16 P16;
            EFp16_init(&P16);
            EFp16_set_EFp(&P16, &P);
            battery16("Ate_Pairing",  Ate_Pairing,  &P16, &Q,
                      order_r, a, b, ab, rminus1);
            battery16("Tate_Pairing", Tate_Pairing, &P16, &Q,
                      order_r, a, b, ab, rminus1);
            EFp16_clear(&P16);
        }

        Fp16_clear(&z); Fp16_clear(&zr); Fp16_clear(&one); Fp16_clear(&oz);
        Fp16_clear(&zp); Fp16_clear(&z2); Fp16_clear(&zt);
        Fp16_clear(&e_mP); Fp16_clear(&e_mQ); Fp16_clear(&prod);
        EFp_clear(&aP); EFp_clear(&mP);
        EFp16_clear(&bQ); EFp16_clear(&mQ);
    }

    /* ---------- Frobenius identity: Fp16_frobenius_map(x) == x^p ---------- */
    {
        struct Fp16 xr, fx, xp;
        Fp16_init(&xr); Fp16_init(&fx); Fp16_init(&xp);
        Fp16_random(&xr);
        Fp16_frobenius_map(&fx, &xr);
        Fp16_pow(&xp, &xr, PRIME_P);
        CHECK(Fp16_cmp(&fx, &xp) == 0, "Frobenius: Fp16_frobenius_map(x) == x^p");
        Fp16_clear(&xr); Fp16_clear(&fx); Fp16_clear(&xp);
    }

    mpz_clear(a); mpz_clear(b); mpz_clear(ab); mpz_clear(rminus1);
    mpz_clears(x, y, (mpz_ptr)NULL);
    EFp_clear(&R); EFp_clear(&P); EFp16_clear(&Q);

    mpz_clear(PRIME_P); mpz_clear(order_r); mpz_clear(order_EFp);
    mpz_clear(trace_t); mpz_clear(a_x); mpz_clear(X);

    MTEST_REPORT("KSS16");
    return 1; /* unreachable */
}
