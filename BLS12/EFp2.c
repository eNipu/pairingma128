#include "EFp2.h"

void EFp2_init(struct EFp2 *P) {
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    P->flag = 0;
}

void EFp2_set(struct EFp2 *P, struct EFp2 *A) {
    Fp2_set(&P->x, &A->x);
    Fp2_set(&P->y, &A->y);
    P->flag = A->flag;
}

void EFp2_set_ui(struct EFp2 *P, unsigned long int a) {
    Fp2_set_ui(&P->x, a);
    Fp2_set_ui(&P->y, a);
    P->flag = 0;
}

void EFp2_set_mpz(struct EFp2 *P, mpz_t a) {
    Fp2_set_mpz(&P->x, a);
    Fp2_set_mpz(&P->y, a);
    P->flag = 0;
}

void EFp2_set_neg(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_set(&ANS->x, &P->x);
    Fp2_set_neg(&ANS->y, &P->y);
    ANS->flag = P->flag;
}

void EFp2_clear(struct EFp2 *P) {
    Fp2_clear(&P->x);
    Fp2_clear(&P->y);
}

void EFp2_printf(struct EFp2 *P, char *name) {
    printf("%s", name);
    if (P->flag == 0) {
        printf("(");
        Fp2_printf(&P->x, "X");
        printf(",");
        Fp2_printf(&P->y, "Y");
        printf(")");
    } else {
        printf("0");
    }
}

void EFp2_rational_point(struct EFp2 *P) {
    struct Fp2 buf1, buf2, R;
    Fp2_init(&buf1);
    Fp2_init(&buf2);
    Fp2_init(&R);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));

    while (1) {
        Fp2_random(&P->x, state);
        Fp2_mul(&buf1, &P->x, &P->x);
        Fp2_mul(&buf2, &buf1, &P->x);
        Fp2_mul_mpz(&buf1, &P->x, curve_parameter_A);
        Fp2_add(&R, &buf1, &buf2);
        mpz_add(R.x0.x0, R.x0.x0, curve_parameter_B);
        if (Fp2_legendre(&R) == 1) {
            Fp2_sqrt(&P->y, &R);
            break;
        }
    }

    Fp2_clear(&buf1);
    Fp2_clear(&buf2);
    Fp2_clear(&R);
}

void EFp2_ECD(struct EFp2 *ANS, struct EFp2 *P) {
    if (Fp2_cmp_zero(&P->y) == 0) {
        ANS->flag = 1;
        return;
    }

    struct EFp2 Tmp;
    EFp2_init(&Tmp);
    EFp2_set(&Tmp, P);
    struct Fp2 Buf1, Buf2, C;
    Fp2_init(&Buf1);
    Fp2_init(&Buf2);
    Fp2_init(&C);

    Fp2_mul_ui(&Buf1, &Tmp.y, 2);
    Fp2_inv(&Buf1, &Buf1);
    Fp2_mul(&Buf2, &Tmp.x, &Tmp.x);
    Fp2_mul_ui(&Buf2, &Buf2, 3);
    mpz_add(Buf2.x0.x0, Buf2.x0.x0, curve_parameter_A);
    Fp2_mul(&C, &Buf1, &Buf2);
    Fp2_squaring(&Buf1, &C);
    Fp2_mul_ui(&Buf2, &Tmp.x, 2);
    Fp2_sub(&ANS->x, &Buf1, &Buf2);
    Fp2_sub(&Buf1, &Tmp.x, &ANS->x);
    Fp2_mul(&Buf2, &C, &Buf1);
    Fp2_sub(&ANS->y, &Buf2, &Tmp.y);

    Fp2_clear(&Buf1);
    Fp2_clear(&Buf2);
    Fp2_clear(&C);
    EFp2_clear(&Tmp);
}

void EFp2_ECA(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2) {
    if (P1->flag == 1) {
        EFp2_set(ANS, P2);
        return;
    } else if (P2->flag == 1) {
        EFp2_set(ANS, P1);
        return;
    } else if (Fp2_cmp(&P1->x, &P2->x) == 0) {
        if (Fp2_cmp(&P1->y, &P2->y) != 0) {
            ANS->flag = 1;
            return;
        } else {
            EFp2_ECD(ANS, P1);
            return;
        }
    }

    struct EFp2 Tmp1, Tmp2;
    EFp2_init(&Tmp1);
    EFp2_set(&Tmp1, P1);
    EFp2_init(&Tmp2);
    EFp2_set(&Tmp2, P2);
    struct Fp2 Buf1, Buf2, C;
    Fp2_init(&Buf1);
    Fp2_init(&Buf2);
    Fp2_init(&C);

    Fp2_sub(&Buf1, &Tmp2.x, &Tmp1.x);
    Fp2_inv(&Buf1, &Buf1);
    Fp2_sub(&Buf2, &Tmp2.y, &Tmp1.y);
    Fp2_mul(&C, &Buf1, &Buf2);
    Fp2_squaring(&Buf1, &C);
    Fp2_sub(&Buf2, &Buf1, &Tmp1.x);
    Fp2_sub(&ANS->x, &Buf2, &Tmp2.x);
    Fp2_sub(&Buf1, &Tmp1.x, &ANS->x);
    Fp2_mul(&Buf2, &C, &Buf1);
    Fp2_sub(&ANS->y, &Buf2, &Tmp1.y);

    Fp2_clear(&Buf1);
    Fp2_clear(&Buf2);
    Fp2_clear(&C);
    EFp2_clear(&Tmp1);
    EFp2_clear(&Tmp2);
}

void EFp2_SCM(struct EFp2 *ANS, struct EFp2 *P, mpz_t R) {
    if (mpz_cmp_ui(R, 0) == 0) {
        ANS->flag = 1;
        return;
    } else if (mpz_cmp_ui(R, 1) == 0) {
        EFp2_set(ANS, P);
        return;
    }

    struct EFp2 Tmp, next_P;
    EFp2_init(&Tmp);
    EFp2_set(&Tmp, P);
    EFp2_init(&next_P);
    int i, length;
    length = (int)mpz_sizeinbase(R, 2);
    char binary[length];
    mpz_get_str(binary, 2, R);

    EFp2_set(&next_P, &Tmp);
    for (i = 1; binary[i] != '\0'; i++) {
        EFp2_ECD(&next_P, &next_P);
        if (binary[i] == '1') {
            EFp2_ECA(&next_P, &next_P, &Tmp);
        }
    }

    EFp2_set(ANS, &next_P);

    EFp2_clear(&next_P);
    EFp2_clear(&Tmp);
}

void EFp2_frobenius_1(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_frobenius_1(&ANS->x, &P->x);
    Fp2_frobenius_1(&ANS->y, &P->y);
}

void EFp2_frobenius_2(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_frobenius_2(&ANS->x, &P->x);
    Fp2_frobenius_2(&ANS->y, &P->y);
}

void EFp2_frobenius_3(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_frobenius_3(&ANS->x, &P->x);
    Fp2_frobenius_3(&ANS->y, &P->y);
}

void EFp2_frobenius_4(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_frobenius_4(&ANS->x, &P->x);
    Fp2_frobenius_4(&ANS->y, &P->y);
}

void EFp2_frobenius_6(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_frobenius_6(&ANS->x, &P->x);
    Fp2_frobenius_6(&ANS->y, &P->y);
}

void EFp2_frobenius_8(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_frobenius_8(&ANS->x, &P->x);
    Fp2_frobenius_8(&ANS->y, &P->y);
}

void EFp2_frobenius_10(struct EFp2 *ANS, struct EFp2 *P) {
    Fp2_frobenius_10(&ANS->x, &P->x);
    Fp2_frobenius_10(&ANS->y, &P->y);
}

void EFp2_skew_frobenius_1(struct EFp2 *ANS, struct EFp2 *A) {
    Fp2_frobenius_1(&ANS->x, &A->x);
    Fp2_mul(&ANS->x, &ANS->x, &Fp2_basis_inv_prime_1_div_3);
    Fp2_frobenius_1(&ANS->y, &A->y);
    Fp2_mul(&ANS->y, &ANS->y, &Fp2_basis_inv_prime_1_div_2);
}

void EFp2_skew_frobenius_2(struct EFp2 *ANS, struct EFp2 *A) {
    Fp2_frobenius_2(&ANS->x, &A->x);
    Fp2_mul(&ANS->x, &ANS->x, &Fp2_basis_inv_prime_2_div_3);
    Fp2_frobenius_2(&ANS->y, &A->y);
    Fp2_mul(&ANS->y, &ANS->y, &Fp2_basis_inv_prime_2_div_2);
}

void EFp2_skew_frobenius_3(struct EFp2 *ANS, struct EFp2 *A) {
    Fp2_frobenius_3(&ANS->x, &A->x);
    Fp2_mul(&ANS->x, &ANS->x, &Fp2_basis_inv_prime_3_div_3);
    Fp2_frobenius_3(&ANS->y, &A->y);
    Fp2_mul(&ANS->y, &ANS->y, &Fp2_basis_inv_prime_3_div_2);
}

void EFp2_skew_frobenius_4(struct EFp2 *ANS, struct EFp2 *A) {
    Fp2_mul_mpz(&ANS->x, &A->x, Fp2_basis_prime_4_div_3_2.x0.x0);
    Fp2_set(&ANS->y, &A->y);
}

void EFp2_skew_frobenius_10(struct EFp2 *ANS, struct EFp2 *A) {
    Fp2_frobenius_10(&ANS->x, &A->x);
    Fp2_mul(&ANS->x, &ANS->x, &Fp2_basis_inv_prime_10_div_3);
    Fp2_frobenius_10(&ANS->y, &A->y);
    Fp2_mul(&ANS->y, &ANS->y, &Fp2_basis_inv_prime_10_div_2);
}
