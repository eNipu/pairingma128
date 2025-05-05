#include "EFp.h"

void EFp_init(struct EFp *P) {
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->flag = 0;
}

void EFp_set(struct EFp *P, struct EFp *A) {
    Fp_set(&P->x, &A->x);
    Fp_set(&P->y, &A->y);
    P->flag = A->flag;
}

void EFp_set_ui(struct EFp *P, unsigned long int a) {
    Fp_set_ui(&P->x, a);
    Fp_set_ui(&P->y, a);
    P->flag = 0;
}

void EFp_set_mpz(struct EFp *P, mpz_t a) {
    Fp_set_mpz(&P->x, a);
    Fp_set_mpz(&P->y, a);
    P->flag = 0;
}

void EFp_set_neg(struct EFp *P, struct EFp *A) {
    Fp_set(&P->x, &A->x);
    Fp_set_neg(&P->y, &A->y);
    P->flag = A->flag;
}

void EFp_clear(struct EFp *P) {
    Fp_clear(&P->x);
    Fp_clear(&P->y);
}

void EFp_printf(struct EFp *P, char *name) {
    printf("%s", name);
    if (P->flag == 0) {
        printf("(");
        Fp_printf(&P->x, "");
        printf(",");
        Fp_printf(&P->y, "");
        printf(")");
    } else {
        printf("0");
    }
}

void EFp_rational_point(struct EFp *P) {
    struct Fp buf1, buf2, R;
    Fp_init(&buf1);
    Fp_init(&buf2);
    Fp_init(&R);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));

    while (1) {
        Fp_random(&P->x, state);
        Fp_mul(&buf1, &P->x, &P->x);
        Fp_mul(&buf2, &buf1, &P->x);
        Fp_mul_mpz(&buf1, &P->x, curve_parameter_A);
        Fp_add(&R, &buf1, &buf2);
        Fp_add_mpz(&R, &R, curve_parameter_B);
        if (Fp_legendre(&R) == 1) {
            Fp_sqrt(&P->y, &R);
            break;
        }
    }

    Fp_clear(&buf1);
    Fp_clear(&buf2);
    Fp_clear(&R);
}

void EFp_ECD(struct EFp *ANS, struct EFp *P) {
    if (Fp_cmp_zero(&P->y) == 0) {
        ANS->flag = 1;
        return;
    }

    struct EFp Tmp;
    EFp_init(&Tmp);
    EFp_set(&Tmp, P);
    struct Fp Buf1, Buf2, C;
    Fp_init(&Buf1);
    Fp_init(&Buf2);
    Fp_init(&C);

    Fp_mul_ui(&Buf1, &Tmp.y, 2);
    Fp_inv(&Buf1, &Buf1);
    Fp_mul(&Buf2, &Tmp.x, &Tmp.x);
    Fp_mul_ui(&Buf2, &Buf2, 3);
    Fp_add_mpz(&Buf2, &Buf2, curve_parameter_A);
    Fp_mul(&C, &Buf1, &Buf2);
    Fp_mul(&Buf1, &C, &C);
    Fp_mul_ui(&Buf2, &Tmp.x, 2);
    Fp_sub(&ANS->x, &Buf1, &Buf2);
    Fp_sub(&Buf1, &Tmp.x, &ANS->x);
    Fp_mul(&Buf2, &C, &Buf1);
    Fp_sub(&ANS->y, &Buf2, &Tmp.y);

    Fp_clear(&Buf1);
    Fp_clear(&Buf2);
    Fp_clear(&C);
    EFp_clear(&Tmp);
}

void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2) {
    if (P1->flag == 1) {
        EFp_set(ANS, P2);
        return;
    } else if (P2->flag == 1) {
        EFp_set(ANS, P1);
        return;
    } else if (Fp_cmp(&P1->x, &P2->x) == 0) {
        if (Fp_cmp(&P1->y, &P2->y) != 0) {
            ANS->flag = 1;
            return;
        } else {
            EFp_ECD(ANS, P1);
            return;
        }
    }

    struct EFp Tmp1, Tmp2;
    EFp_init(&Tmp1);
    EFp_set(&Tmp1, P1);
    EFp_init(&Tmp2);
    EFp_set(&Tmp2, P2);
    struct Fp Buf1, Buf2, C;
    Fp_init(&Buf1);
    Fp_init(&Buf2);
    Fp_init(&C);

    Fp_sub(&Buf1, &Tmp2.x, &Tmp1.x);
    Fp_inv(&Buf1, &Buf1);
    Fp_sub(&Buf2, &Tmp2.y, &Tmp1.y);
    Fp_mul(&C, &Buf1, &Buf2);
    Fp_mul(&Buf1, &C, &C);
    Fp_sub(&Buf2, &Buf1, &Tmp1.x);
    Fp_sub(&ANS->x, &Buf2, &Tmp2.x);
    Fp_sub(&Buf1, &Tmp1.x, &ANS->x);
    Fp_mul(&Buf2, &C, &Buf1);
    Fp_sub(&ANS->y, &Buf2, &Tmp1.y);

    Fp_clear(&Buf1);
    Fp_clear(&Buf2);
    Fp_clear(&C);
    EFp_clear(&Tmp1);
    EFp_clear(&Tmp2);
}

void EFp_SCM(struct EFp *ANS, struct EFp *P, mpz_t R) {
    if (mpz_cmp_ui(R, 0) == 0) {
        ANS->flag = 1;
        return;
    } else if (mpz_cmp_ui(R, 1) == 0) {
        EFp_set(ANS, P);
        return;
    }

    struct EFp Tmp, next_P;
    EFp_init(&Tmp);
    EFp_set(&Tmp, P);
    EFp_init(&next_P);
    int i, length;
    length = (int)mpz_sizeinbase(R, 2);
    char binary[length];
    mpz_get_str(binary, 2, R);

    EFp_set(&next_P, &Tmp);
    for (i = 1; binary[i] != '\0'; i++) {
        EFp_ECD(&next_P, &next_P);
        if (binary[i] == '1') {
            EFp_ECA(&next_P, &next_P, &Tmp);
        }
    }

    EFp_set(ANS, &next_P);

    EFp_clear(&next_P);
    EFp_clear(&Tmp);
}

void EFp_skew_frobenius(struct EFp *ANS, struct EFp *A) {
    Fp_mul(&ANS->x, &A->x, &inv_CNR1);
    Fp_set_neg(&ANS->y, &A->y);
}
