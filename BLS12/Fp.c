#include "Fp.h"

void Fp_init(struct Fp *P) {
    mpz_init(P->x0);
}

void Fp_set(struct Fp *P, struct Fp *A) {
    mpz_set(P->x0, A->x0);
}

void Fp_set_ui(struct Fp *P, unsigned long int a) {
    mpz_set_ui(P->x0, a);
}

void Fp_set_mpz(struct Fp *P, mpz_t a) {
    mpz_set(P->x0, a);
}

void Fp_set_neg(struct Fp *P, struct Fp *A) {
    mpz_sub(P->x0, prime, A->x0);
}

void Fp_random(struct Fp *P, gmp_randstate_t state) {
    mpz_urandomm(P->x0, state, prime);
}

void Fp_clear(struct Fp *P) {
    mpz_clear(P->x0);
}

void Fp_printf(struct Fp *P, char *name) {
    printf("%s", name);
    mpz_out_str(stdout, 10, P->x0);
}

void Fp_mul(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_mul(ANS->x0, A->x0, B->x0);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_mul_ui(struct Fp *ANS, struct Fp *A, unsigned long int a) {
    mpz_mul_ui(ANS->x0, A->x0, a);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_mul_mpz(struct Fp *ANS, struct Fp *A, mpz_t a) {
    mpz_mul(ANS->x0, A->x0, a);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_add(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_add(ANS->x0, A->x0, B->x0);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_add_ui(struct Fp *ANS, struct Fp *A, unsigned long int a) {
    mpz_add_ui(ANS->x0, A->x0, a);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_add_mpz(struct Fp *ANS, struct Fp *A, mpz_t a) {
    mpz_add(ANS->x0, A->x0, a);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_sub(struct Fp *ANS, struct Fp *A, struct Fp *B) {
    mpz_sub(ANS->x0, A->x0, B->x0);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_sub_ui(struct Fp *ANS, struct Fp *A, unsigned long int a) {
    mpz_sub_ui(ANS->x0, A->x0, a);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_sub_mpz(struct Fp *ANS, struct Fp *A, mpz_t a) {
    mpz_sub(ANS->x0, A->x0, a);
    mpz_mod(ANS->x0, ANS->x0, prime);
}

void Fp_inv(struct Fp *ANS, struct Fp *A) {
    mpz_invert(ANS->x0, A->x0, prime);
}

int Fp_legendre(struct Fp *A) {
    return mpz_legendre(A->x0, prime);
}

int Fp_isCNR(struct Fp *A) {
    struct Fp buf;
    Fp_init(&buf);
    mpz_t exp;
    mpz_init(exp);

    mpz_sub_ui(exp, prime, 1);
    mpz_tdiv_q_ui(exp, exp, 3);
    Fp_pow(&buf, A, exp);

    mpz_clear(exp);
    if (Fp_cmp_one(&buf) == 0) {
        Fp_clear(&buf);
        return 1;
    } else if (Fp_cmp_zero(&buf) == 0) {
        Fp_clear(&buf);
        return 0;
    } else {
        Fp_clear(&buf);
        return -1;
    }
}

void Fp_sqrt(struct Fp *ANS, struct Fp *A) {
    struct Fp x, y, t, k, n, buf;
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&t);
    Fp_init(&k);
    Fp_init(&n);
    Fp_init(&buf);
    unsigned long int e, m;
    mpz_t exp, q, z, result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));
    Fp_random(&n, state);

    while (Fp_legendre(&n) != -1) {
        Fp_random(&n, state);
    }
    mpz_sub_ui(q, prime, 1);
    mpz_mod_ui(result, q, 2);
    e = 0;
    while (mpz_cmp_ui(result, 0) == 0) {
        mpz_tdiv_q_ui(q, q, 2);
        mpz_mod_ui(result, q, 2);
        e++;
    }
    Fp_pow(&y, &n, q);
    mpz_set_ui(z, e);
    mpz_sub_ui(exp, q, 1);
    mpz_tdiv_q_ui(exp, exp, 2);
    Fp_pow(&x, A, exp);
    Fp_mul(&buf, &x, &x);
    Fp_mul(&k, &buf, A);
    Fp_mul(&x, &x, A);
    while (Fp_cmp_one(&k) != 0) {
        m = 1;
        mpz_ui_pow_ui(exp, 2, m);
        Fp_pow(&buf, &k, exp);
        while (Fp_cmp_one(&buf) != 0) {
            m++;
            mpz_ui_pow_ui(exp, 2, m);
            Fp_pow(&buf, &k, exp);
        }
        mpz_sub_ui(exp, z, m);
        mpz_sub_ui(exp, exp, 1);
        mpz_ui_pow_ui(result, 2, mpz_get_ui(exp));
        Fp_pow(&t, &y, result);
        Fp_mul(&y, &t, &t);
        mpz_set_ui(z, m);
        Fp_mul(&x, &x, &t);
        Fp_mul(&k, &k, &y);
    }
    Fp_set(ANS, &x);

    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&t);
    Fp_clear(&k);
    Fp_clear(&n);
    Fp_clear(&buf);
}

void Fp_pow(struct Fp *ANS, struct Fp *A, mpz_t a) {
    int i, length;
    length = (int)mpz_sizeinbase(a, 2);
    char binary[length];
    mpz_get_str(binary, 2, a);
    struct Fp buf;
    Fp_init(&buf);
    Fp_set(&buf, A);

    for (i = 1; binary[i] != '\0'; i++) {
        Fp_mul(&buf, &buf, &buf);
        if (binary[i] == '1') {
            Fp_mul(&buf, A, &buf);
        }
    }
    Fp_set(ANS, &buf);

    Fp_clear(&buf);
}

int Fp_cmp(struct Fp *A, struct Fp *B) {
    if (mpz_cmp(A->x0, B->x0) == 0) {
        return 0;
    }
    return 1;
}

int Fp_cmp_ui(struct Fp *A, unsigned long int a) {
    if (mpz_cmp_ui(A->x0, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp_cmp_mpz(struct Fp *A, mpz_t a) {
    if (mpz_cmp(A->x0, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp_cmp_zero(struct Fp *A) {
    if (mpz_cmp_ui(A->x0, 0) == 0) {
        return 0;
    }
    return 1;
}

int Fp_cmp_one(struct Fp *A) {
    if (mpz_cmp_ui(A->x0, 1) == 0) {
        return 0;
    }
    return 1;
}

void Fp_frobenius_1(struct Fp *ANS, struct Fp *A) {
    Fp_set(ANS, A);
}

void Fp_frobenius_2(struct Fp *ANS, struct Fp *A) {
    Fp_set(ANS, A);
}

void Fp_frobenius_3(struct Fp *ANS, struct Fp *A) {
    Fp_set(ANS, A);
}

void Fp_frobenius_4(struct Fp *ANS, struct Fp *A) {
    Fp_set(ANS, A);
}

void Fp_frobenius_6(struct Fp *ANS, struct Fp *A) {
    Fp_set(ANS, A);
}

void Fp_frobenius_8(struct Fp *ANS, struct Fp *A) {
    Fp_set(ANS, A);
}

void Fp_frobenius_10(struct Fp *ANS, struct Fp *A) {
    Fp_set(ANS, A);
}
