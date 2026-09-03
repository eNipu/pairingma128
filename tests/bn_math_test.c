/* bn_math_test.c - mathematical correctness suite for the BN optimal-ate
 * pairing (reference: eprint 2017/1174).
 *
 * TDD gate: every check below is a *mathematical identity* that must hold for a
 * correct implementation. Failure => nonzero exit.
 *
 * G1/G2 points are generated at runtime (as the original main does), so the
 * pairing VALUES are not deterministic between runs, but every asserted
 * identity holds for any valid choice of points/scalars.
 */
#include <stdio.h>
#include <gmp.h>
#include "math_assert.h"
#include "field.h"
#include "curve.h"
#include "pairing.h"
#include "params.h"

int main(void){
    gmp_randstate_t state;
    mpz_t a, b, ab, rminus1;

    init_parameters();
    set_parameters();

    gmp_printf("-- BN: p bits=%d r bits=%d #E(Fp) bits=%d --\n",
               (int)mpz_sizeinbase(prime,2), (int)mpz_sizeinbase(EFp_order,2),
               (int)mpz_sizeinbase(EFp_total,2));

    /* ---------- setup sanity ---------- */
    CHECK(mpz_probab_prime_p(EFp_order, 40) >= 1,
          "r is (probable) prime");
    {
        mpz_t rem; mpz_init(rem);
        mpz_mod(rem, EFp_total, EFp_order);
        CHECK(mpz_sgn(rem) == 0, "r divides #E(Fp)");
        mpz_clear(rem);
    }

    /* ---------- points ---------- */
    /* R: random point on E(Fp). P: cofactor-cleared into G1 = E(Fp)[r]. */
    struct EFp12 R, P, Q;
    mpz_t h1;
    EFp12_init(&R); EFp12_init(&P); EFp12_init(&Q);
    mpz_init(h1);
    EFp12_generate_G1(&R);
    EFp12_generate_G2(&Q);
    mpz_tdiv_q(h1, EFp_total, EFp_order); /* cofactor of #E(Fp) over r */
    EFp12_SCM(&P, &R, h1);
    mpz_clear(h1);
    EFp12_clear(&R);
    CHECK(P.flag == 0, "P (cofactor-cleared) is not the point at infinity");
    CHECK(Q.flag == 0, "Q is not the point at infinity");

    /* P must lie on E(Fp): y^2 = x^3 + A x + B */
    {
        struct Fp x, y, y2, x2, x3, ax, rhs;
        Fp_init(&x); Fp_init(&y); Fp_init(&y2); Fp_init(&x2); Fp_init(&x3);
        Fp_init(&ax); Fp_init(&rhs);
        Fp_set(&x, &P.x.x0.x0.x0);
        Fp_set(&y, &P.y.x0.x0.x0);
        Fp_mul(&y2, &y, &y);
        Fp_mul(&x2, &x, &x);
        Fp_mul(&x3, &x2, &x);
        Fp_mul_mpz(&ax, &x, curve_parameter_A);
        Fp_add(&rhs, &x3, &ax);
        Fp_sub_mpz(&rhs, &rhs, curve_parameter_B);
        CHECK(Fp_cmp(&y2, &rhs) == 0, "generated P is on E(Fp): y^2 == x^3 + A x - B");        Fp_clear(&x); Fp_clear(&y); Fp_clear(&y2); Fp_clear(&x2); Fp_clear(&x3);
        Fp_clear(&ax); Fp_clear(&rhs);
    }

    mpz_init(a); mpz_init(b); mpz_init(ab); mpz_init(rminus1);
    mpz_set_ui(a, 31); mpz_set_ui(b, 11); mpz_mul(ab, a, b);
    mpz_sub_ui(rminus1, EFp_order, 1);

    /* ---------- point orders: [r]P = O, [r]Q = O ---------- */
    {
        struct EFp12 T;
        EFp12_init(&T);
        EFp12_G1_SCM(&T, &P, EFp_order);
        CHECK(T.flag == 1, "[r]P == point at infinity (P has order r)");
        EFp12_clear(&T);
        EFp12_init(&T);
        EFp12_G2_SCM_2div(&T, &Q, EFp_order);
        CHECK(T.flag == 1, "[r]Q == point at infinity (Q has order r)");
        EFp12_clear(&T);
    }

    /* ---------- core pairing identities ---------- */
    {
        struct Fp12 z, zr, zi, one, oz, zp, z2, zt, e_mP, e_mQ, prod;
        struct EFp12 aP, bQ, mP, mQ;
        Fp12_init(&z); Fp12_init(&zr); Fp12_init(&zi); Fp12_init(&one);
        Fp12_init(&oz); Fp12_init(&zp); Fp12_init(&z2); Fp12_init(&zt);
        Fp12_init(&e_mP); Fp12_init(&e_mQ); Fp12_init(&prod);
        EFp12_init(&aP); EFp12_init(&bQ); EFp12_init(&mP); EFp12_init(&mQ);

        /* reference value */
        Opt_ate_pairing(&z, &Q, &P);

        /* multiplicative identity, derived (no assumption about Fp12_ONE) */
        Fp12_inv(&zi, &z);
        Fp12_mul(&one, &z, &zi);
        Fp12_mul(&oz, &one, &z);
        CHECK(Fp12_cmp(&oz, &z) == 0, "GT identity element sanity (one*z == z)");

        /* non-degeneracy */
        CHECK(Fp12_cmp(&z, &one) != 0, "e(P,Q) != 1 (non-degenerate)");

        /* reduced pairing output lies in the order-r subgroup */
        Fp12_pow(&zp, &z, EFp_order);
        CHECK(Fp12_cmp(&zp, &one) == 0, "e(P,Q)^r == 1 (subgroup)");

        /* bilinearity: e([a]P,[b]Q) == e(P,Q)^{a*b} */
        EFp12_G1_SCM(&aP, &P, a);
        EFp12_G2_SCM_2div(&bQ, &Q, b);
        Opt_ate_pairing(&z2, &bQ, &aP);
        Fp12_pow(&zt, &z, ab);
        CHECK(Fp12_cmp(&z2, &zt) == 0, "e([a]P,[b]Q) == e(P,Q)^(ab) (bilinear, a=31 b=11)");

        /* alternating in arg1 */
        EFp12_G1_SCM(&mP, &P, rminus1);
        Opt_ate_pairing(&e_mP, &Q, &mP);
        Fp12_mul(&prod, &e_mP, &z);
        CHECK(Fp12_cmp(&prod, &one) == 0, "e(-P,Q)*e(P,Q) == 1 (alternating in G1)");

        /* alternating in arg2 */
        EFp12_G2_SCM_2div(&mQ, &Q, rminus1);
        Opt_ate_pairing(&e_mQ, &mQ, &P);
        Fp12_mul(&prod, &e_mQ, &z);
        CHECK(Fp12_cmp(&prod, &one) == 0, "e(P,-Q)*e(P,Q) == 1 (alternating in G2)");

        /* reproducibility */
        Opt_ate_pairing(&zr, &Q, &P);
        CHECK(Fp12_cmp(&zr, &z) == 0, "pairing is deterministic on identical inputs");

        EFp12_clear(&aP); EFp12_clear(&bQ); EFp12_clear(&mP); EFp12_clear(&mQ);
        Fp12_clear(&z); Fp12_clear(&zr); Fp12_clear(&zi); Fp12_clear(&one);
        Fp12_clear(&oz); Fp12_clear(&zp); Fp12_clear(&z2); Fp12_clear(&zt);
        Fp12_clear(&e_mP); Fp12_clear(&e_mQ); Fp12_clear(&prod);
    }

    /* ---------- field Frobenius identity: pi_k(x) == x^{p^k} ---------- */
    {
        struct Fp12 x, fx, xp, f2, fx2;
        Fp12_init(&x); Fp12_init(&fx); Fp12_init(&xp);
        Fp12_init(&f2); Fp12_init(&fx2);
        gmp_randinit_default(state);
        gmp_randseed_ui(state, 42UL);
        Fp12_random(&x, state);
        Fp12_frobenius_1(&fx, &x);
        Fp12_pow(&xp, &x, prime);
        CHECK(Fp12_cmp(&fx, &xp) == 0, "Frobenius: pi(x) == x^p (Fp12)");
        Fp12_frobenius_2(&f2, &x);
        Fp12_frobenius_1(&fx2, &fx);
        CHECK(Fp12_cmp(&f2, &fx2) == 0, "Frobenius: pi2(x) == pi(pi(x)) (Fp12)");
        gmp_randclear(state);
        Fp12_clear(&x); Fp12_clear(&fx); Fp12_clear(&xp);
        Fp12_clear(&f2); Fp12_clear(&fx2);
    }

    mpz_clear(a); mpz_clear(b); mpz_clear(ab); mpz_clear(rminus1);
    EFp12_clear(&P); EFp12_clear(&Q);

    clear_parameters();
    MTEST_REPORT("BN");
    return 1; /* unreachable */
}
