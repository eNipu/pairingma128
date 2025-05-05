#include "Fp2.h"

void Fp2_init(struct Fp2 *P) {
    Fp_init(&P->x0);
    Fp_init(&P->x1);
}

void Fp2_set(struct Fp2 *P, struct Fp2 *A) {
    Fp_set(&P->x0, &A->x0);
    Fp_set(&P->x1, &A->x1);
}

void Fp2_set_ui(struct Fp2 *P, unsigned long int a) {
    Fp_set_ui(&P->x0, a);
    Fp_set_ui(&P->x1, a);
}

void Fp2_set_mpz(struct Fp2 *P, mpz_t a) {
    Fp_set_mpz(&P->x0, a);
    Fp_set_mpz(&P->x1, a);
}

void Fp2_set_neg(struct Fp2 *P, struct Fp2 *A) {
    Fp_set_neg(&P->x0, &A->x0);
    Fp_set_neg(&P->x1, &A->x1);
}

void Fp2_random(struct Fp2 *P, gmp_randstate_t state) {
    Fp_random(&P->x0, state);
    Fp_random(&P->x1, state);
}

void Fp2_clear(struct Fp2 *P) {
    Fp_clear(&P->x0);
    Fp_clear(&P->x1);
}

void Fp2_printf(struct Fp2 *P, char *name) {
    printf("%s", name);
    printf("(");
    Fp_printf(&P->x0, "");
    printf(",");
    Fp_printf(&P->x1, "");
    printf(")");
}

void Fp2_mul(struct Fp2 *ANS, struct Fp2 *A, struct Fp2 *B) {
    struct Fp tmp1, tmp2, tmp3, tmp4;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    Fp_init(&tmp4);

    Fp_mul(&tmp1, &A->x0, &B->x0);
    Fp_mul(&tmp2, &A->x1, &B->x1);
    Fp_add(&tmp3, &A->x0, &A->x1);
    Fp_add(&tmp4, &B->x0, &B->x1);

    Fp_mul_basis(&ANS->x0, &tmp2);
    Fp_add(&ANS->x0, &ANS->x0, &tmp1);

    Fp_mul(&ANS->x1, &tmp3, &tmp4);
    Fp_sub(&ANS->x1, &ANS->x1, &tmp1);
    Fp_sub(&ANS->x1, &ANS->x1, &tmp2);

    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
    Fp_clear(&tmp4);
}

void Fp2_mul_ui(struct Fp2 *ANS, struct Fp2 *A, unsigned long int a) {
    Fp_mul_ui(&ANS->x0, &A->x0, a);
    Fp_mul_ui(&ANS->x1, &A->x1, a);
}

void Fp2_mul_mpz(struct Fp2 *ANS, struct Fp2 *A, mpz_t a) {
    Fp_mul_mpz(&ANS->x0, &A->x0, a);
    Fp_mul_mpz(&ANS->x1, &A->x1, a);
}

void Fp2_mul_basis(struct Fp2 *ANS, struct Fp2 *A) {
    struct Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp, &A->x0);

    Fp_sub(&ANS->x0, &tmp, &A->x1);
    Fp_add(&ANS->x1, &tmp, &A->x1);

    Fp_clear(&tmp);
}

void Fp2_squaring(struct Fp2 *ANS, struct Fp2 *A) {
    struct Fp tmp1, tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);

    Fp_add(&tmp1, &A->x0, &A->x1);
    Fp_sub(&tmp2, &A->x0, &A->x1);

    Fp_mul(&ANS->x1, &A->x0, &A->x1);
    Fp_add(&ANS->x1, &ANS->x1, &ANS->x1);

    Fp_mul(&ANS->x0, &tmp1, &tmp2);

    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
}

void Fp2_inv_basis(struct Fp2 *ANS, struct Fp2 *A) {
    struct Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp, A);

    Fp_add(&ANS->x0, &tmp.x0, &tmp.x1);
    Fp_mul_mpz(&ANS->x0, &ANS->x0, Fp2_basis_inv.x0.x0);
    Fp_sub(&ANS->x1, &tmp.x1, &tmp.x0);
    Fp_mul_mpz(&ANS->x1, &ANS->x1, Fp2_basis_inv.x0.x0);

    Fp2_clear(&tmp);
}

void Fp2_add(struct Fp2 *ANS, struct Fp2 *A, struct Fp2 *B) {
    Fp_add(&ANS->x0, &A->x0, &B->x0);
    Fp_add(&ANS->x1, &A->x1, &B->x1);
}

void Fp2_add_ui(struct Fp2 *ANS, struct Fp2 *A, unsigned long int a) {
    Fp_add_ui(&ANS->x0, &A->x0, a);
    Fp_add_ui(&ANS->x1, &A->x1, a);
}

void Fp2_add_mpz(struct Fp2 *ANS, struct Fp2 *A, mpz_t a) {
    Fp_add_mpz(&ANS->x0, &A->x0, a);
    Fp_add_mpz(&ANS->x1, &A->x1, a);
}

void Fp2_sub(struct Fp2 *ANS, struct Fp2 *A, struct Fp2 *B) {
    Fp_sub(&ANS->x0, &A->x0, &B->x0);
    Fp_sub(&ANS->x1, &A->x1, &B->x1);
}

void Fp2_sub_ui(struct Fp2 *ANS, struct Fp2 *A, unsigned long int a) {
    Fp_sub_ui(&ANS->x0, &A->x0, a);
    Fp_sub_ui(&ANS->x1, &A->x1, a);
}

void Fp2_sub_mpz(struct Fp2 *ANS, struct Fp2 *A, mpz_t a) {
    Fp_sub_mpz(&ANS->x0, &A->x0, a);
    Fp_sub_mpz(&ANS->x1, &A->x1, a);
}

void Fp2_inv(struct Fp2 *ANS, struct Fp2 *A) {
    struct Fp2 frob, buf;
    Fp2_init(&frob);
    Fp2_init(&buf);

    Fp2_inv_map(&frob, A);
    Fp2_mul(&buf, A, &frob);
    Fp_inv(&buf.x0, &buf.x0);
    Fp2_mul_mpz(ANS, &frob, buf.x0.x0);

    Fp2_clear(&frob);
    Fp2_clear(&buf);
}

void Fp2_inv_map(struct Fp2 *ANS, struct Fp2 *A) {
    Fp_set(&ANS->x0, &A->x0);
    Fp_set_neg(&ANS->x1, &A->x1);
}

int Fp2_legendre(struct Fp2 *A) {
    struct Fp2 buf;
    Fp2_init(&buf);

    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp, prime, 2);
    mpz_sub_ui(exp, exp, 1);
    mpz_tdiv_q_ui(exp, exp, 2);
    Fp2_pow(&buf, A, exp);

    mpz_clear(exp);
    if (Fp2_cmp_one(&buf) == 0) {
        Fp2_clear(&buf);
        return 1;
    } else if (Fp2_cmp_zero(&buf) == 0) {
        Fp2_clear(&buf);
        return 0;
    } else {
        Fp2_clear(&buf);
        return -1;
    }
}

void Fp2_sqrt(struct Fp2 *ANS, struct Fp2 *A) {
    struct Fp2 x, y, t, k, n, buf;
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&t);
    Fp2_init(&k);
    Fp2_init(&n);
    Fp2_init(&buf);
    unsigned long int e, m;
    mpz_t exp, q, z, result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));

    Fp2_random(&n, state);
    while (Fp2_legendre(&n) != -1) {
        Fp2_random(&n, state);
    }
    mpz_pow_ui(q, prime, 2);
    mpz_sub_ui(q, q, 1);
    mpz_mod_ui(result, q, 2);
    e = 0;
    while (mpz_cmp_ui(result, 0) == 0) {
        mpz_tdiv_q_ui(q, q, 2);
        mpz_mod_ui(result, q, 2);
        e++;
    }
    Fp2_pow(&y, &n, q);
    mpz_set_ui(z, e);
    mpz_sub_ui(exp, q, 1);
    mpz_tdiv_q_ui(exp, exp, 2);
    Fp2_pow(&x, A, exp);
    Fp2_mul(&buf, &x, &x);
    Fp2_mul(&k, &buf, A);
    Fp2_mul(&x, &x, A);
    while (Fp2_cmp_one(&k) != 0) {
        m = 1;
        mpz_ui_pow_ui(exp, 2, m);
        Fp2_pow(&buf, &k, exp);
        while (Fp2_cmp_one(&buf) != 0) {
            m++;
            mpz_ui_pow_ui(exp, 2, m);
            Fp2_pow(&buf, &k, exp);
        }
        mpz_sub_ui(exp, z, m);
        mpz_sub_ui(exp, exp, 1);
        mpz_ui_pow_ui(result, 2, mpz_get_ui(exp));
        Fp2_pow(&t, &y, result);
        Fp2_mul(&y, &t, &t);
        mpz_set_ui(z, m);
        Fp2_mul(&x, &x, &t);
        Fp2_mul(&k, &k, &y);
    }
    Fp2_set(ANS, &x);

    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&t);
    Fp2_clear(&k);
    Fp2_clear(&n);
    Fp2_clear(&buf);
}

void Fp2_pow(struct Fp2 *ANS, struct Fp2 *A, mpz_t a) {
    int i, length;
    length = (int)mpz_sizeinbase(a, 2);
    char binary[length];
    mpz_get_str(binary, 2, a);
    struct Fp2 buf;
    Fp2_init(&buf);
    Fp2_set(&buf, A);

    for (i = 1; binary[i] != '\0'; i++) {
        Fp2_squaring(&buf, &buf);
        if (binary[i] == '1') {
            Fp2_mul(&buf, A, &buf);
        }
    }

    Fp2_set(ANS, &buf);
    Fp2_clear(&buf);
}

int Fp2_cmp(struct Fp2 *A, struct Fp2 *B) {
    if (Fp_cmp(&A->x0, &B->x0) == 0 && Fp_cmp(&A->x1, &B->x1) == 0) {
        return 0;
    }
    return 1;
}

int Fp2_cmp_ui(struct Fp2 *A, unsigned long int a) {
    if (Fp_cmp_ui(&A->x0, a) == 0 && Fp_cmp_ui(&A->x1, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp2_cmp_mpz(struct Fp2 *A, mpz_t a) {
    if (Fp_cmp_mpz(&A->x0, a) == 0 && Fp_cmp_mpz(&A->x1, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp2_cmp_zero(struct Fp2 *A) {
    if (Fp_cmp_zero(&A->x0) == 0 && Fp_cmp_zero(&A->x1) == 0) {
        return 0;
    }
    return 1;
}

int Fp2_cmp_one(struct Fp2 *A) {
    if (Fp_cmp_one(&A->x0) == 0 && Fp_cmp_zero(&A->x1) == 0) {
        return 0;
    }
    return 1;
}

void Fp2_frobenius_1(struct Fp2 *ANS, struct Fp2 *A) {
    Fp_set(&ANS->x0, &A->x0);
    Fp_set_neg(&ANS->x1, &A->x1);
}

void Fp2_frobenius_2(struct Fp2 *ANS, struct Fp2 *A) {
    Fp2_set(ANS, A);
}

void Fp2_frobenius_3(struct Fp2 *ANS, struct Fp2 *A) {
    Fp_set(&ANS->x0, &A->x0);
    Fp_set_neg(&ANS->x1, &A->x1);
}

void Fp2_frobenius_4(struct Fp2 *ANS, struct Fp2 *A) {
    Fp2_set(ANS, A);
}

void Fp2_frobenius_6(struct Fp2 *ANS, struct Fp2 *A) {
    Fp2_set(ANS, A);
}

void Fp2_frobenius_8(struct Fp2 *ANS, struct Fp2 *A) {
    Fp2_set(ANS, A);
}

void Fp2_frobenius_10(struct Fp2 *ANS, struct Fp2 *A) {
    Fp2_set(ANS, A);
}
