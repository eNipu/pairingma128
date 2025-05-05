#include "Fp12.h"

void Fp12_init(struct Fp12 *P) {
    Fp6_init(&P->x0);
    Fp6_init(&P->x1);
}

void Fp12_set(struct Fp12 *P, struct Fp12 *A) {
    Fp6_set(&P->x0, &A->x0);
    Fp6_set(&P->x1, &A->x1);
}

void Fp12_set_ui(struct Fp12 *P, unsigned long int a) {
    Fp6_set_ui(&P->x0, a);
    Fp6_set_ui(&P->x1, a);
}

void Fp12_set_mpz(struct Fp12 *P, mpz_t a) {
    Fp6_set_mpz(&P->x0, a);
    Fp6_set_mpz(&P->x1, a);
}

void Fp12_set_neg(struct Fp12 *P, struct Fp12 *A) {
    Fp6_set_neg(&P->x0, &A->x0);
    Fp6_set_neg(&P->x1, &A->x1);
}

void Fp12_random(struct Fp12 *P, gmp_randstate_t state) {
    Fp6_random(&P->x0, state);
    Fp6_random(&P->x1, state);
}

void Fp12_clear(struct Fp12 *P) {
    Fp6_clear(&P->x0);
    Fp6_clear(&P->x1);
}

void Fp12_printf(struct Fp12 *P, char *name) {
    printf("%s", name);
    printf("(");
    Fp6_printf(&P->x0, "");
    printf(",");
    Fp6_printf(&P->x1, "");
    printf(")");
}

void Fp12_mul(struct Fp12 *ANS, struct Fp12 *A, struct Fp12 *B) {
    struct Fp6 tmp1, tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);

    Fp6_mul(&tmp2, &A->x1, &B->x1);
    Fp6_add(&tmp1, &A->x0, &A->x1);
    Fp6_add(&ANS->x1, &B->x0, &B->x1);
    Fp6_mul(&ANS->x1, &tmp1, &ANS->x1);
    Fp6_mul(&tmp1, &A->x0, &B->x0);

    Fp6_mul_basis(&ANS->x0, &tmp2);
    Fp6_add(&ANS->x0, &ANS->x0, &tmp1);

    Fp6_sub(&ANS->x1, &ANS->x1, &tmp1);
    Fp6_sub(&ANS->x1, &ANS->x1, &tmp2);

    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}

void Fp12_mul_ui(struct Fp12 *ANS, struct Fp12 *A, unsigned long int a) {
    Fp6_mul_ui(&ANS->x0, &A->x0, a);
    Fp6_mul_ui(&ANS->x1, &A->x1, a);
}

void Fp12_mul_mpz(struct Fp12 *ANS, struct Fp12 *A, mpz_t a) {
    Fp6_mul_mpz(&ANS->x0, &A->x0, a);
    Fp6_mul_mpz(&ANS->x1, &A->x1, a);
}

void Fp12_squaring(struct Fp12 *ANS, struct Fp12 *A) {
    struct Fp6 tmp1, tmp2, tmp3;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&tmp3);

    Fp6_add(&tmp1, &A->x0, &A->x1);
    Fp6_mul_basis(&tmp2, &A->x1);
    Fp6_add(&tmp2, &tmp2, &A->x0);
    Fp6_mul(&tmp3, &A->x0, &A->x1);

    Fp6_mul(&ANS->x0, &tmp1, &tmp2);
    Fp6_sub(&ANS->x0, &ANS->x0, &tmp3);
    Fp6_mul_basis(&tmp1, &tmp3);
    Fp6_sub(&ANS->x0, &ANS->x0, &tmp1);

    Fp6_add(&ANS->x1, &tmp3, &tmp3);

    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&tmp3);
}

void Fp12_add(struct Fp12 *ANS, struct Fp12 *A, struct Fp12 *B) {
    Fp6_add(&ANS->x0, &A->x0, &B->x0);
    Fp6_add(&ANS->x1, &A->x1, &B->x1);
}

void Fp12_add_ui(struct Fp12 *ANS, struct Fp12 *A, unsigned long int a) {
    Fp6_add_ui(&ANS->x0, &A->x0, a);
    Fp6_add_ui(&ANS->x1, &A->x1, a);
}

void Fp12_add_mpz(struct Fp12 *ANS, struct Fp12 *A, mpz_t a) {
    Fp6_add_mpz(&ANS->x0, &ANS->x0, a);
    Fp6_add_mpz(&ANS->x1, &ANS->x1, a);
}

void Fp12_sub(struct Fp12 *ANS, struct Fp12 *A, struct Fp12 *B) {
    Fp6_sub(&ANS->x0, &A->x0, &B->x0);
    Fp6_sub(&ANS->x1, &A->x1, &B->x1);
}

void Fp12_sub_ui(struct Fp12 *ANS, struct Fp12 *A, unsigned long int a) {
    Fp6_sub_ui(&ANS->x0, &ANS->x0, a);
    Fp6_sub_ui(&ANS->x1, &ANS->x1, a);
}

void Fp12_sub_mpz(struct Fp12 *ANS, struct Fp12 *A, mpz_t a) {
    Fp6_sub_mpz(&ANS->x0, &ANS->x0, a);
    Fp6_sub_mpz(&ANS->x1, &ANS->x1, a);
}

void Fp12_inv(struct Fp12 *ANS, struct Fp12 *A) {
    struct Fp12 frob, buf;
    Fp12_init(&frob);
    Fp12_init(&buf);

    Fp12_inv_map(&frob, A);
    Fp12_mul(&buf, A, &frob);
    Fp6_inv(&buf.x0, &buf.x0);
    Fp6_mul(&ANS->x0, &frob.x0, &buf.x0);
    Fp6_mul(&ANS->x1, &frob.x1, &buf.x0);

    Fp12_clear(&frob);
    Fp12_clear(&buf);
}

void Fp12_inv_map(struct Fp12 *ANS, struct Fp12 *A) {
    Fp6_set(&ANS->x0, &A->x0);
    Fp6_set_neg(&ANS->x1, &A->x1);
}

int Fp12_legendre(struct Fp12 *A) {
    mpz_t exp;
    mpz_init(exp);
    struct Fp12 buf;
    Fp12_init(&buf);

    mpz_pow_ui(exp, prime, 12);
    mpz_sub_ui(exp, exp, 1);
    mpz_tdiv_q_ui(exp, exp, 2);
    Fp12_pow(&buf, A, exp);

    mpz_clear(exp);
    if (Fp12_cmp_one(&buf) == 0) {
        Fp12_clear(&buf);
        return 1;
    } else if (Fp12_cmp_zero(&buf) == 0) {
        Fp12_clear(&buf);
        return 0;
    } else {
        Fp12_clear(&buf);
        return -1;
    }
}

void Fp12_sqrt(struct Fp12 *ANS, struct Fp12 *A) {
    struct Fp12 x, y, t, k, n, buf;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&buf);
    unsigned long int e, m;
    mpz_t exp, q, z, result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));

    Fp12_random(&n, state);
    while (Fp12_legendre(&n) != -1) {
        Fp12_random(&n, state);
    }
    mpz_pow_ui(q, prime, 12);
    mpz_sub_ui(q, q, 1);
    mpz_mod_ui(result, q, 2);
    e = 0;
    while (mpz_cmp_ui(result, 0) == 0) {
        mpz_tdiv_q_ui(q, q, 2);
        mpz_mod_ui(result, q, 2);
        e++;
    }
    Fp12_pow(&y, &n, q);
    mpz_set_ui(z, e);
    mpz_sub_ui(exp, q, 1);
    mpz_tdiv_q_ui(exp, exp, 2);
    Fp12_pow(&x, A, exp);
    Fp12_mul(&buf, &x, &x);
    Fp12_mul(&k, &buf, A);
    Fp12_mul(&x, &x, A);
    while (Fp12_cmp_one(&k) != 0) {
        m = 1;
        mpz_ui_pow_ui(exp, 2, m);
        Fp12_pow(&buf, &k, exp);
        while (Fp12_cmp_one(&buf) != 0) {
            m++;
            mpz_ui_pow_ui(exp, 2, m);
            Fp12_pow(&buf, &k, exp);
        }
        mpz_sub_ui(exp, z, m);
        mpz_sub_ui(exp, exp, 1);
        mpz_ui_pow_ui(result, 2, mpz_get_ui(exp));
        Fp12_pow(&t, &y, result);
        Fp12_mul(&y, &t, &t);
        mpz_set_ui(z, m);
        Fp12_mul(&x, &x, &t);
        Fp12_mul(&k, &k, &y);
    }
    Fp12_set(ANS, &x);

    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
    Fp12_clear(&x);
    Fp12_clear(&y);
    Fp12_clear(&t);
    Fp12_clear(&k);
    Fp12_clear(&n);
    Fp12_clear(&buf);
}

void Fp12_pow(struct Fp12 *ANS, struct Fp12 *A, mpz_t a) {
    int i, length;
    length = (int)mpz_sizeinbase(a, 2);
    char binary[length];
    mpz_get_str(binary, 2, a);
    struct Fp12 buf;
    Fp12_init(&buf);
    Fp12_set(&buf, A);

    for (i = 1; binary[i] != '\0'; i++) {
        Fp12_squaring(&buf, &buf);
        if (binary[i] == '1') {
            Fp12_mul(&buf, A, &buf);
        }
    }

    Fp12_set(ANS, &buf);
    Fp12_clear(&buf);
}

int Fp12_cmp(struct Fp12 *A, struct Fp12 *B) {
    if (Fp6_cmp(&A->x0, &B->x0) == 0 && Fp6_cmp(&A->x1, &B->x1) == 0) {
        return 0;
    }
    return 1;
}

int Fp12_cmp_ui(struct Fp12 *A, unsigned long int a) {
    if (Fp6_cmp_ui(&A->x0, a) == 0 && Fp6_cmp_ui(&A->x1, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp12_cmp_mpz(struct Fp12 *A, mpz_t a) {
    if (Fp6_cmp_mpz(&A->x0, a) == 0 && Fp6_cmp_mpz(&A->x1, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp12_cmp_zero(struct Fp12 *A) {
    if (Fp6_cmp_zero(&A->x0) == 0 && Fp6_cmp_zero(&A->x1) == 0) {
        return 0;
    }
    return 1;
}

int Fp12_cmp_one(struct Fp12 *A) {
    if (Fp6_cmp_one(&A->x0) == 0 && Fp6_cmp_zero(&A->x1) == 0) {
        return 0;
    }
    return 1;
}

void Fp12_frobenius_1(struct Fp12 *ANS, struct Fp12 *A) {
    struct Fp tmp;
    Fp_init(&tmp);

    Fp_set(&ANS->x0.x0.x0, &A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1, &A->x0.x0.x1);
    Fp_set(&tmp, &A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0, &A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1, &tmp);
    Fp2_mul_mpz(&ANS->x0.x1, &ANS->x0.x1, Fp2_basis_prime_1_div_3_1.x1.x0);

    Fp_set(&ANS->x0.x2.x0, &A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1, &A->x0.x2.x1);
    Fp2_mul_mpz(&ANS->x0.x2, &ANS->x0.x2, Fp2_basis_prime_1_div_3_2.x0.x0);

    Fp_set(&ANS->x1.x0.x0, &A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1, &A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0, &ANS->x1.x0, &Fp2_basis_prime_1_div_6);

    Fp_set(&tmp, &A->x1.x1.x0);
    Fp_set(&ANS->x1.x1.x0, &A->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1, &tmp);
    Fp2_mul_mpz(&ANS->x1.x1, &ANS->x1.x1, Fp2_basis_prime_1_div_3_1.x1.x0);
    Fp2_mul(&ANS->x1.x1, &ANS->x1.x1, &Fp2_basis_prime_1_div_6);

    Fp_set(&ANS->x1.x2.x0, &A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1, &A->x1.x2.x1);
    Fp2_mul_mpz(&ANS->x1.x2, &ANS->x1.x2, Fp2_basis_prime_1_div_3_2.x0.x0);
    Fp2_mul(&ANS->x1.x2, &ANS->x1.x2, &Fp2_basis_prime_1_div_6);

    Fp_clear(&tmp);
}

void Fp12_frobenius_2(struct Fp12 *ANS, struct Fp12 *A) {
    Fp2_set(&ANS->x0.x0, &A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1, &A->x0.x1, Fp2_basis_prime_2_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2, &A->x0.x2, Fp2_basis_prime_2_div_3_2.x0.x0);

    Fp2_mul_mpz(&ANS->x1.x0, &A->x1.x0, Fp2_basis_prime_2_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1, &A->x1.x1, Fp2_basis_prime_2_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1, &ANS->x1.x1, Fp2_basis_prime_2_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2, &A->x1.x2, Fp2_basis_prime_2_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2, &ANS->x1.x2, Fp2_basis_prime_2_div_6.x0.x0);
}

void Fp12_frobenius_3(struct Fp12 *ANS, struct Fp12 *A) {
    struct Fp tmp;
    Fp_init(&tmp);

    Fp_set(&ANS->x0.x0.x0, &A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1, &A->x0.x0.x1);
    Fp_set(&tmp, &A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0, &A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1, &tmp);
    Fp_set_neg(&ANS->x0.x2.x0, &A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1, &A->x0.x2.x1);

    Fp_set(&ANS->x1.x0.x0, &A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1, &A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0, &ANS->x1.x0, &Fp2_basis_prime_3_div_6);

    Fp_set(&tmp, &A->x1.x1.x0);
    Fp_set(&ANS->x1.x1.x0, &A->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1, &tmp);
    Fp2_mul(&ANS->x1.x1, &ANS->x1.x1, &Fp2_basis_prime_3_div_6);

    Fp_set_neg(&ANS->x1.x2.x0, &A->x1.x2.x0);
    Fp_set(&ANS->x1.x2.x1, &A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2, &ANS->x1.x2, &Fp2_basis_prime_3_div_6);

    Fp_clear(&tmp);
}

void Fp12_frobenius_4(struct Fp12 *ANS, struct Fp12 *A) {
    Fp2_set(&ANS->x0.x0, &A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1, &A->x0.x1, Fp2_basis_prime_4_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2, &A->x0.x2, Fp2_basis_prime_4_div_3_2.x0.x0);

    Fp2_mul_mpz(&ANS->x1.x0, &A->x1.x0, Fp2_basis_prime_4_div_6.x0.x0);
    Fp2_set(&ANS->x1.x1, &A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2, &A->x1.x2, Fp2_basis_prime_4_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2, &ANS->x1.x2, Fp2_basis_prime_4_div_6.x0.x0);
}

void Fp12_frobenius_6(struct Fp12 *ANS, struct Fp12 *A) {
    Fp6_set(&ANS->x0, &A->x0);
    Fp6_set_neg(&ANS->x1, &A->x1);
}

void Fp12_frobenius_8(struct Fp12 *ANS, struct Fp12 *A) {
    Fp2_set(&ANS->x0.x0, &A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1, &A->x0.x1, Fp2_basis_prime_8_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2, &A->x0.x2, Fp2_basis_prime_8_div_3_2.x0.x0);

    Fp2_mul_mpz(&ANS->x1.x0, &A->x1.x0, Fp2_basis_prime_8_div_6.x0.x0);
    Fp2_set(&ANS->x1.x1, &A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2, &A->x1.x2, Fp2_basis_prime_8_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2, &ANS->x1.x2, Fp2_basis_prime_8_div_6.x0.x0);
}

void Fp12_frobenius_10(struct Fp12 *ANS, struct Fp12 *A) {
    Fp2_set(&ANS->x0.x0, &A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1, &A->x0.x1, Fp2_basis_prime_10_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2, &A->x0.x2, Fp2_basis_prime_10_div_3_2.x0.x0);

    Fp2_mul_mpz(&ANS->x1.x0, &A->x1.x0, Fp2_basis_prime_10_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1, &A->x1.x1, Fp2_basis_prime_10_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1, &ANS->x1.x1, Fp2_basis_prime_10_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2, &A->x1.x2, Fp2_basis_prime_10_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2, &ANS->x1.x2, Fp2_basis_prime_10_div_6.x0.x0);
}
