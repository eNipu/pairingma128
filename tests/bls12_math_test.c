/* bls12_math_test.c - mathematical correctness suite for the BLS12 optimal-ate
 * pairing (reference: eprint 2017/1174).
 *
 * TDD gate: every check below is a *mathematical identity* that must hold for a
 * correct implementation. Failure => nonzero exit. The pairing value e(P,Q) on
 * the fixed points below is deterministic and is printed as hex at the end so a
 * post-refactor run can be diffed against a pre-refactor run for exact parity.
 *
 * Compile against the monolithic source (pre-refactor) or the split modules
 * (post-refactor); this file only depends on the public API.
 */
#include <stdio.h>
#include <gmp.h>
#include "math_assert.h"
#include "field.h"
#include "curve.h"
#include "pairing.h"
#include "params.h"

static void fp_hex(FILE *f, mpz_t x){ gmp_fprintf(f, "%Zx", x); }

int main(void){
    gmp_randstate_t state;
    mpz_t a, b, ab, rminus1, tmp;

    init_parameters();
    set_parameters();

    gmp_printf("-- BLS12: p bits=%d r bits=%d #E(Fp) bits=%d --\n",
               (int)mpz_sizeinbase(prime,2), (int)mpz_sizeinbase(EFp_order,2),
               (int)mpz_sizeinbase(EFp_total,2));

    /* ---------- setup sanity ---------- */
    CHECK(mpz_probab_prime_p(EFp_order, 40) == 2 || mpz_probab_prime_p(EFp_order, 40) == 1,
          "r is (probable) prime");
    mpz_init(tmp); mpz_mod(tmp, EFp_total, EFp_order);
    CHECK(mpz_sgn(tmp) == 0, "r divides #E(Fp)");
    mpz_clear(tmp);

    /* ---------- deterministic test points ---------- */
    /* R: fixed point on E(Fp). P: cofactor-cleared into G1 = E(Fp)[r]. */
    struct EFp12 R, P, Q;
    mpz_t h1;
    EFp12_init(&R); EFp12_init(&P); EFp12_init(&Q);
    mpz_init(h1);
    R.flag = 0; P.flag = 0; Q.flag = 0;
    mpz_set_str(R.x.x0.x0.x0.x0, "1577263467895074691751656396777546817676899779239812956811535508811145903113183105315855416956275013763132428190349521477916841537817217749", 10);
    mpz_set_str(R.y.x0.x0.x0.x0, "1482623217062500042620995282298816432047975830544022708334643354696131009901634528980080824925757942965723245158102345395591975011061998174", 10);
    mpz_set_str(Q.x.x0.x2.x0.x0, "1716096407462707739249584854219035887699085689300518488744116638769194609975441542417090151676720383761847673487573410049924094001307902771", 10);
    mpz_set_str(Q.x.x0.x2.x1.x0, "1884774744539223134856244661850811554160388946147660758387803399777014198117917660870489415649762044750981220496189685117651039305687807588", 10);
    mpz_set_str(Q.y.x1.x1.x0.x0, "862911509351712110897605125200443999396611767860495548282350223126731430707014429484782130102632482076556245081709440912377961089052176300", 10);
    mpz_set_str(Q.y.x1.x1.x1.x0, "14917642869256953479675328417653271711757572421578255369823731891492952725387400103789740069957122839470285312924022325651200303158794609", 10);

    /* R must lie on E(Fp): y^2 = x^3 + A x + B  (BLS12: A=0, B=4) */
    {
        struct Fp x, y, y2, x2, x3, ax, rhs;
        Fp_init(&x); Fp_init(&y); Fp_init(&y2); Fp_init(&x2); Fp_init(&x3);
        Fp_init(&ax); Fp_init(&rhs);
        Fp_set(&x, &R.x.x0.x0.x0);
        Fp_set(&y, &R.y.x0.x0.x0);
        Fp_mul(&y2, &y, &y);
        Fp_mul(&x2, &x, &x);
        Fp_mul(&x3, &x2, &x);
        Fp_mul_mpz(&ax, &x, curve_parameter_A);
        Fp_add(&rhs, &x3, &ax);
        Fp_add_mpz(&rhs, &rhs, curve_parameter_B);
        CHECK(Fp_cmp(&y2, &rhs) == 0, "fixed R is on E(Fp)");
        Fp_clear(&x); Fp_clear(&y); Fp_clear(&y2); Fp_clear(&x2); Fp_clear(&x3);
        Fp_clear(&ax); Fp_clear(&rhs);
    }
    CHECK(R.flag == 0, "R is not the point at infinity");
    CHECK(Q.flag == 0, "Q is not the point at infinity");

    /* cofactor-clear R into the order-r subgroup of E(Fp) */
    mpz_tdiv_q(h1, EFp_total, EFp_order); /* h1 = #E(Fp) / r */
    EFp12_SCM(&P, &R, h1);
    mpz_clear(h1);
    CHECK(P.flag == 0, "P = [h1]R (cofactor-cleared) is not the point at infinity");

    mpz_init(a); mpz_init(b); mpz_init(ab); mpz_init(rminus1);
    mpz_set_ui(a, 31); mpz_set_ui(b, 11); mpz_mul(ab, a, b);
    mpz_sub_ui(rminus1, EFp_order, 1);

    /* ---------- point orders: [r]P = O, [r]Q = O ---------- */
    {
        struct EFp12 T;
        EFp12_init(&T);
        EFp12_SCM(&T, &P, EFp_order);
        CHECK(T.flag == 1, "[r]P == point at infinity (P has order r)");
        EFp12_clear(&T);
        EFp12_init(&T);
        EFp12_SCM(&T, &Q, EFp_order);
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
        printf("pairing value digest (run-to-run deterministic on fixed points):\n");
        printf("  z_00: "); fp_hex(stdout, z.x0.x0.x0.x0); printf("\n");
        printf("  z_01: "); fp_hex(stdout, z.x0.x2.x0.x0); printf("\n");
        printf("  z_10: "); fp_hex(stdout, z.x1.x0.x0.x0); printf("\n");

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
        EFp12_SCM(&aP, &P, a);
        EFp12_SCM(&bQ, &Q, b);
        Opt_ate_pairing(&z2, &bQ, &aP);
        Fp12_pow(&zt, &z, ab);
        CHECK(Fp12_cmp(&z2, &zt) == 0, "e([a]P,[b]Q) == e(P,Q)^(ab) (bilinear, a=31 b=11)");

        /* alternating/antisymmetry in arg1: e([r-1]P,Q)*e(P,Q) == 1 */
        EFp12_SCM(&mP, &P, rminus1);
        Opt_ate_pairing(&e_mP, &Q, &mP);
        Fp12_mul(&prod, &e_mP, &z);
        CHECK(Fp12_cmp(&prod, &one) == 0, "e(-P,Q)*e(P,Q) == 1 (alternating in G1)");

        /* alternating in arg2 */
        EFp12_SCM(&mQ, &Q, rminus1);
        Opt_ate_pairing(&e_mQ, &mQ, &P);
        Fp12_mul(&prod, &e_mQ, &z);
        CHECK(Fp12_cmp(&prod, &one) == 0, "e(P,-Q)*e(P,Q) == 1 (alternating in G2)");

        /* reproducibility on identical inputs */
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
    EFp12_clear(&R); EFp12_clear(&P); EFp12_clear(&Q);

    clear_parameters();
    MTEST_REPORT("BLS12");
    return 1; /* unreachable */
}
