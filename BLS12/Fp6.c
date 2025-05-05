#include "Fp6.h"

void Fp6_init(struct Fp6 *P) {
    Fp2_init(&P->x0);
    Fp2_init(&P->x1);
    Fp2_init(&P->x2);
}

void Fp6_set(struct Fp6 *P, struct Fp6 *A) {
    Fp2_set(&P->x0, &A->x0);
    Fp2_set(&P->x1, &A->x1);
    Fp2_set(&P->x2, &A->x2);
}

void Fp6_set_ui(struct Fp6 *P, unsigned long int a) {
    Fp2_set_ui(&P->x0, a);
    Fp2_set_ui(&P->x1, a);
    Fp2_set_ui(&P->x2, a);
}

void Fp6_set_mpz(struct Fp6 *P, mpz_t a) {
    Fp2_set_mpz(&P->x0, a);
    Fp2_set_mpz(&P->x1, a);
    Fp2_set_mpz(&P->x2, a);
}

void Fp6_set_neg(struct Fp6 *P, struct Fp6 *A) {
    Fp2_set_neg(&P->x0, &A->x0);
    Fp2_set_neg(&P->x1, &A->x1);
    Fp2_set_neg(&P->x2, &A->x2);
}

void Fp6_random(struct Fp6 *P, gmp_randstate_t state) {
    Fp2_random(&P->x0, state);
    Fp2_random(&P->x1, state);
    Fp2_random(&P->x2, state);
}

void Fp6_clear(struct Fp6 *P) {
    Fp2_clear(&P->x0);
    Fp2_clear(&P->x1);
    Fp2_clear(&P->x2);
}

void Fp6_printf(struct Fp6 *P, char *name) {
    printf("%s", name);
    printf("(");
    Fp2_printf(&P->x0, "");
    printf(",");
    Fp2_printf(&P->x1, "");
    printf(",");
    Fp2_printf(&P->x2, "");
    printf(")");
}

void Fp6_mul(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B) {
    struct Fp2 tmp00, tmp11, tmp22, buf, t0, t1, t2;
    Fp2_init(&tmp00);
    Fp2_init(&tmp11);
    Fp2_init(&tmp22);
    Fp2_init(&buf);
    Fp2_init(&t0);
    Fp2_init(&t1);
    Fp2_init(&t2);

    Fp2_mul(&tmp00, &A->x0, &B->x0);
    Fp2_mul(&tmp11, &A->x1, &B->x1);
    Fp2_mul(&tmp22, &A->x2, &B->x2);

    Fp2_add(&t0, &A->x0, &A->x1);
    Fp2_add(&buf, &B->x0, &B->x1);
    Fp2_mul(&t0, &t0, &buf);

    Fp2_add(&t1, &A->x1, &A->x2);
    Fp2_add(&buf, &B->x1, &B->x2);
    Fp2_mul(&t1, &t1, &buf);

    Fp2_add(&t2, &B->x0, &B->x2);
    Fp2_add(&buf, &A->x0, &A->x2);
    Fp2_mul(&t2, &t2, &buf);

    Fp2_sub(&t1, &t1, &tmp11);
    Fp2_sub(&t1, &t1, &tmp22);
    Fp2_mul_basis(&buf, &t1);
    Fp2_add(&ANS->x0, &tmp00, &buf);

    Fp2_sub(&t0, &t0, &tmp00);
    Fp2_sub(&t0, &t0, &tmp11);
    Fp2_mul_basis(&buf, &tmp22);
    Fp2_add(&ANS->x1, &buf, &t0);

    Fp2_sub(&t2, &t2, &tmp00);
    Fp2_sub(&t2, &t2, &tmp22);
    Fp2_add(&ANS->x2, &tmp11, &t2);

    Fp2_clear(&tmp00);
    Fp2_clear(&tmp11);
    Fp2_clear(&tmp22);
    Fp2_clear(&buf);
    Fp2_clear(&t0);
    Fp2_clear(&t1);
    Fp2_clear(&t2);
}

void Fp6_mul_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int a) {
    Fp2_mul_ui(&ANS->x0, &A->x0, a);
    Fp2_mul_ui(&ANS->x1, &A->x1, a);
    Fp2_mul_ui(&ANS->x2, &A->x2, a);
}

void Fp6_mul_mpz(struct Fp6 *ANS, struct Fp6 *A, mpz_t a) {
    Fp2_mul_mpz(&ANS->x0, &A->x0, a);
    Fp2_mul_mpz(&ANS->x1, &A->x1, a);
    Fp2_mul_mpz(&ANS->x2, &A->x2, a);
}

void Fp6_mul_basis(struct Fp6 *ANS, struct Fp6 *A) {
    struct Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp, A);

    Fp_sub(&ANS->x0.x0, &tmp.x2.x0, &tmp.x2.x1);
    Fp_add(&ANS->x0.x1, &tmp.x2.x0, &tmp.x2.x1);
    Fp_set(&ANS->x1.x0, &tmp.x0.x0);
    Fp_set(&ANS->x1.x1, &tmp.x0.x1);
    Fp_set(&ANS->x2.x0, &tmp.x1.x0);
    Fp_set(&ANS->x2.x1, &tmp.x1.x1);

    Fp6_clear(&tmp);
}

void Fp6_squaring(struct Fp6 *ANS, struct Fp6 *A) {
    struct Fp2 tmp00, tmp12_2, tmp01_2, tmp22, buf;
    Fp2_init(&tmp00);
    Fp2_init(&tmp22);
    Fp2_init(&tmp12_2);
    Fp2_init(&tmp01_2);
    Fp2_init(&buf);

    Fp2_squaring(&tmp00, &A->x0);
    Fp2_squaring(&tmp22, &A->x2);
    Fp2_add(&buf, &A->x1, &A->x1);
    Fp2_mul(&tmp12_2, &buf, &A->x2);
    Fp2_mul(&tmp01_2, &A->x0, &buf);
    Fp2_add(&buf, &A->x0, &A->x1);
    Fp2_add(&buf, &buf, &A->x2);

    Fp2_mul_basis(&ANS->x0, &tmp12_2);
    Fp2_add(&ANS->x0, &ANS->x0, &tmp00);

    Fp2_mul_basis(&ANS->x1, &tmp22);
    Fp2_add(&ANS->x1, &ANS->x1, &tmp01_2);

    Fp2_squaring(&ANS->x2, &buf);
    Fp2_add(&buf, &tmp00, &tmp22);
    Fp2_add(&buf, &buf, &tmp12_2);
    Fp2_add(&buf, &buf, &tmp01_2);
    Fp2_sub(&ANS->x2, &ANS->x2, &buf);

    Fp2_clear(&tmp00);
    Fp2_clear(&tmp22);
    Fp2_clear(&tmp12_2);
    Fp2_clear(&tmp01_2);
    Fp2_clear(&buf);
}

void Fp6_add(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B) {
    Fp2_add(&ANS->x0, &A->x0, &B->x0);
    Fp2_add(&ANS->x1, &A->x1, &B->x1);
    Fp2_add(&ANS->x2, &A->x2, &B->x2);
}

void Fp6_add_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int a) {
    Fp2_add_ui(&ANS->x0, &A->x0, a);
    Fp2_add_ui(&ANS->x1, &A->x1, a);
    Fp2_add_ui(&ANS->x2, &A->x2, a);
}

void Fp6_add_mpz(struct Fp6 *ANS, struct Fp6 *A, mpz_t a) {
    Fp2_add_mpz(&ANS->x0, &A->x0, a);
    Fp2_add_mpz(&ANS->x1, &A->x1, a);
    Fp2_add_mpz(&ANS->x2, &A->x2, a);
}

void Fp6_sub(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B) {
    Fp2_sub(&ANS->x0, &A->x0, &B->x0);
    Fp2_sub(&ANS->x1, &A->x1, &B->x1);
    Fp2_sub(&ANS->x2, &A->x2, &B->x2);
}

void Fp6_sub_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int a) {
    Fp2_sub_ui(&ANS->x0, &A->x0, a);
    Fp2_sub_ui(&ANS->x1, &A->x1, a);
    Fp2_sub_ui(&ANS->x2, &A->x2, a);
}

void Fp6_sub_mpz(struct Fp6 *ANS, struct Fp6 *A, mpz_t a) {
    Fp2_sub_mpz(&ANS->x0, &A->x0, a);
    Fp2_sub_mpz(&ANS->x1, &A->x1, a);
    Fp2_sub_mpz(&ANS->x2, &A->x2, a);
}

void Fp6_inv(struct Fp6 *ANS, struct Fp6 *A) {
    struct Fp6 frob1, frob2, buf1, buf2;
    Fp6_init(&frob1);
    Fp6_init(&frob2);
    Fp6_init(&buf1);
    Fp6_init(&buf2);

    Fp6_inv_map_1(&frob1, A);
    Fp6_inv_map_2(&frob2, A);
    Fp6_mul(&buf1, &frob1, &frob2);
    Fp6_mul(&buf2, &buf1, A);
    Fp2_inv(&buf2.x0, &buf2.x0);
    Fp2_mul(&ANS->x0, &buf1.x0, &buf2.x0);
    Fp2_mul(&ANS->x1, &buf1.x1, &buf2.x0);
    Fp2_mul(&ANS->x2, &buf1.x2, &buf2.x0);

    Fp6_clear(&frob1);
    Fp6_clear(&frob2);
    Fp6_clear(&buf1);
    Fp6_clear(&buf2);
}

void Fp6_inv_map_1(struct Fp6 *ANS, struct Fp6 *A) {
    Fp2_set(&ANS->x0, &A->x0);
    Fp2_mul_mpz(&ANS->x1, &A->x1, inv_CNR1.x0);
    Fp2_mul_mpz(&ANS->x2, &A->x2, inv_CNR2.x0);
}

void Fp6_inv_map_2(struct Fp6 *ANS, struct Fp6 *A) {
    Fp2_set(&ANS->x0, &A->x0);
    Fp2_mul_mpz(&ANS->x1, &A->x1, inv_CNR2.x0);
    Fp2_mul_mpz(&ANS->x2, &A->x2, inv_CNR1.x0);
}

int Fp6_legendre(struct Fp6 *A) {
    mpz_t exp;
    mpz_init(exp);
    struct Fp6 buf;
    Fp6_init(&buf);

    mpz_pow_ui(exp, prime, 6);
    mpz_sub_ui(exp, exp, 1);
    mpz_tdiv_q_ui(exp, exp, 2);
    Fp6_pow(&buf, A, exp);

    mpz_clear(exp);
    if (Fp6_cmp_one(&buf) == 0) {
        Fp6_clear(&buf);
        return 1;
    } else if (Fp6_cmp_zero(&buf) == 0) {
        Fp6_clear(&buf);
        return 0;
    } else {
        Fp6_clear(&buf);
        return -1;
    }
}

void Fp6_sqrt(struct Fp6 *ANS, struct Fp6 *A) {
    struct Fp6 buf1, buf2;
    Fp6_init(&buf1);
    Fp6_init(&buf2);
    mpz_t exp, buf;
    mpz_init(exp);
    mpz_init(buf);

    Fp6_frobenius_4(&buf1, A);
    Fp6_frobenius_2(&buf2, A);
    Fp6_mul(&buf1, &buf1, &buf2);
    Fp6_mul(&buf1, &buf1, A);
    Fp6_set_ui(&buf2, 0);
    Fp2_sqrt(&buf2.x0, &buf1.x0);
    Fp2_inv(&buf2.x0, &buf2.x0);
    Fp2_set(&buf2.x0, &buf2.x0);
    mpz_pow_ui(exp, prime, 8);
    mpz_pow_ui(buf, prime, 4);
    mpz_add(exp, exp, buf);
    mpz_add_ui(exp, exp, 2);
    mpz_tdiv_q_ui(exp, exp, 2);
    Fp6_pow(&buf1, A, exp);
    Fp6_mul(&buf1, &buf1, &buf2);
    Fp6_set(ANS, &buf1);

    mpz_clear(exp);
    mpz_clear(buf);
    Fp6_clear(&buf1);
    Fp6_clear(&buf2);
}

void Fp6_pow(struct Fp6 *ANS, struct Fp6 *A, mpz_t a) {
    int i, length;
    length = (int)mpz_sizeinbase(a, 2);
    char binary[length];
    mpz_get_str(binary, 2, a);
    struct Fp6 buf;
    Fp6_init(&buf);
    Fp6_set(&buf, A);

    for (i = 1; binary[i] != '\0'; i++) {
        Fp6_squaring(&buf, &buf);
        if (binary[i] == '1') {
            Fp6_mul(&buf, A, &buf);
        }
    }

    Fp6_set(ANS, &buf);
    Fp6_clear(&buf);
}

int Fp6_cmp(struct Fp6 *A, struct Fp6 *B) {
    if (Fp2_cmp(&A->x0, &B->x0) == 0 && Fp2_cmp(&A->x1, &B->x1) == 0 && Fp2_cmp(&A->x2, &B->x2) == 0) {
        return 0;
    }
    return 1;
}

int Fp6_cmp_ui(struct Fp6 *A, unsigned long int a) {
    if (Fp2_cmp_ui(&A->x0, a) == 0 && Fp2_cmp_ui(&A->x1, a) == 0 && Fp2_cmp_ui(&A->x2, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp6_cmp_mpz(struct Fp6 *A, mpz_t a) {
    if (Fp2_cmp_mpz(&A->x0, a) == 0 && Fp2_cmp_mpz(&A->x1, a) == 0 && Fp2_cmp_mpz(&A->x2, a) == 0) {
        return 0;
    }
    return 1;
}

int Fp6_cmp_zero(struct Fp6 *A) {
    if (Fp2_cmp_zero(&A->x0) == 0 && Fp2_cmp_zero(&A->x1) == 0 && Fp2_cmp_zero(&A->x2) == 0) {
        return 0;
    }
    return 1;
}

int Fp6_cmp_one(struct Fp6 *A) {
    if (Fp2_cmp_one(&A->x0) == 0 && Fp2_cmp_zero(&A->x1) == 0 && Fp2_cmp_zero(&A->x2) == 0) {
        return 0;
    }
    return 1;
}

void Fp6_frobenius_1(struct Fp6 *ANS, struct Fp6 *A) {
    struct Fp tmp;
    Fp_init(&tmp);

    Fp_set(&ANS->x0.x0, &A->x0.x0);
    Fp_set_neg(&ANS->x0.x1, &A->x0.x1);

    Fp_set(&tmp, &A->x1.x0);
    Fp_set(&ANS->x1.x0, &A->x1.x1);
    Fp_set(&ANS->x1.x1, &tmp);
    Fp2_mul_mpz(&ANS->x1, &ANS->x1, Fp2_basis_prime_1_div_3_1.x1.x0);

    Fp_set(&ANS->x2.x0, &A->x2.x0);
    Fp_set_neg(&ANS->x2.x1, &A->x2.x1);
    Fp2_mul_mpz(&ANS->x2, &ANS->x2, Fp2_basis_prime_1_div_3_2.x0.x0);

    Fp_clear(&tmp);
}

void Fp6_frobenius_2(struct Fp6 *ANS, struct Fp6 *A) {
    Fp2_set(&ANS->x0, &A->x0);
    Fp2_mul_mpz(&ANS->x1, &A->x1, Fp2_basis_prime_2_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2, &A->x2, Fp2_basis_prime_2_div_3_2.x0.x0);
}

void Fp6_frobenius_3(struct Fp6 *ANS, struct Fp6 *A) {
    struct Fp tmp;
    Fp_init(&tmp);

    Fp_set(&ANS->x0.x0, &A->x0.x0);
    Fp_set_neg(&ANS->x0.x1, &A->x0.x1);

    Fp_set(&tmp, &A->x1.x0);
    Fp_set(&ANS->x1.x0, &A->x1.x1);
    Fp_set(&ANS->x1.x1, &tmp);

    Fp_set_neg(&ANS->x2.x0, &A->x2.x0);
    Fp_set(&ANS->x2.x1, &A->x2.x1);

    Fp_clear(&tmp);
}

void Fp6_frobenius_4(struct Fp6 *ANS, struct Fp6 *A) {
    Fp2_set(&ANS->x0, &A->x0);
    Fp2_mul_mpz(&ANS->x1, &A->x1, Fp2_basis_prime_4_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2, &A->x2, Fp2_basis_prime_4_div_3_2.x0.x0);
}

void Fp6_frobenius_6(struct Fp6 *ANS, struct Fp6 *A) {
    Fp6_set(ANS, A);
}

void Fp6_frobenius_8(struct Fp6 *ANS, struct Fp6 *A) {
    Fp2_set(&ANS->x0, &A->x0);
    Fp2_mul_mpz(&ANS->x1, &A->x1, Fp2_basis_prime_8_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2, &A->x2, Fp2_basis_prime_8_div_3_2.x0.x0);
}

void Fp6_frobenius_10(struct Fp6 *ANS, struct Fp6 *A) {
    Fp2_set(&ANS->x0, &A->x0);
    Fp2_mul_mpz(&ANS->x1, &A->x1, Fp2_basis_prime_10_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2, &A->x2, Fp2_basis_prime_10_div_3_2.x0.x0);
}
