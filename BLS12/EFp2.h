#ifndef EFP2_H
#define EFP2_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include "Fp2.h"

struct EFp2 {
    struct Fp2 x;
    struct Fp2 y;
    int flag;
};

void EFp2_init(struct EFp2 *P);
void EFp2_set(struct EFp2 *P, struct EFp2 *A);
void EFp2_set_ui(struct EFp2 *P, unsigned long int a);
void EFp2_set_mpz(struct EFp2 *P, mpz_t a);
void EFp2_set_neg(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_clear(struct EFp2 *P);
void EFp2_printf(struct EFp2 *P, char *name);
void EFp2_rational_point(struct EFp2 *P);
void EFp2_ECD(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_ECA(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2);
void EFp2_SCM(struct EFp2 *ANS, struct EFp2 *P, mpz_t R);
void EFp2_frobenius_1(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_frobenius_2(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_frobenius_3(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_frobenius_4(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_frobenius_6(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_frobenius_8(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_frobenius_10(struct EFp2 *ANS, struct EFp2 *P);
void EFp2_skew_frobenius_1(struct EFp2 *ANS, struct EFp2 *A);
void EFp2_skew_frobenius_2(struct EFp2 *ANS, struct EFp2 *A);
void EFp2_skew_frobenius_3(struct EFp2 *ANS, struct EFp2 *A);
void EFp2_skew_frobenius_4(struct EFp2 *ANS, struct EFp2 *A);
void EFp2_skew_frobenius_10(struct EFp2 *ANS, struct EFp2 *A);

#endif // EFP2_H
