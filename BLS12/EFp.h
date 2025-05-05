#ifndef EFP_H
#define EFP_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include "Fp.h"

struct EFp {
    struct Fp x;
    struct Fp y;
    int flag;
};

void EFp_init(struct EFp *P);
void EFp_set(struct EFp *P, struct EFp *A);
void EFp_set_ui(struct EFp *P, unsigned long int a);
void EFp_set_mpz(struct EFp *P, mpz_t a);
void EFp_set_neg(struct EFp *P, struct EFp *A);
void EFp_clear(struct EFp *P);
void EFp_printf(struct EFp *P, char *name);
void EFp_rational_point(struct EFp *P);
void EFp_ECD(struct EFp *ANS, struct EFp *P);
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2);
void EFp_SCM(struct EFp *ANS, struct EFp *P, mpz_t R);
void EFp_skew_frobenius(struct EFp *ANS, struct EFp *A);

#endif // EFP_H
