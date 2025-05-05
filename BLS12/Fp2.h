#ifndef FP2_H
#define FP2_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include "Fp.h"

struct Fp2 {
    struct Fp x0;
    struct Fp x1;
};

void Fp2_init(struct Fp2 *P);
void Fp2_set(struct Fp2 *P, struct Fp2 *A);
void Fp2_set_ui(struct Fp2 *P, unsigned long int a);
void Fp2_set_mpz(struct Fp2 *P, mpz_t a);
void Fp2_set_neg(struct Fp2 *P, struct Fp2 *A);
void Fp2_random(struct Fp2 *P, gmp_randstate_t state);
void Fp2_clear(struct Fp2 *P);
void Fp2_printf(struct Fp2 *P, char *name);
void Fp2_mul(struct Fp2 *ANS, struct Fp2 *A, struct Fp2 *B);
void Fp2_mul_ui(struct Fp2 *ANS, struct Fp2 *A, unsigned long int a);
void Fp2_mul_mpz(struct Fp2 *ANS, struct Fp2 *A, mpz_t a);
void Fp2_mul_basis(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_squaring(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_inv_basis(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_add(struct Fp2 *ANS, struct Fp2 *A, struct Fp2 *B);
void Fp2_add_ui(struct Fp2 *ANS, struct Fp2 *A, unsigned long int a);
void Fp2_add_mpz(struct Fp2 *ANS, struct Fp2 *A, mpz_t a);
void Fp2_sub(struct Fp2 *ANS, struct Fp2 *A, struct Fp2 *B);
void Fp2_sub_ui(struct Fp2 *ANS, struct Fp2 *A, unsigned long int a);
void Fp2_sub_mpz(struct Fp2 *ANS, struct Fp2 *A, mpz_t a);
void Fp2_inv(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_inv_map(struct Fp2 *ANS, struct Fp2 *A);
int Fp2_legendre(struct Fp2 *A);
void Fp2_sqrt(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_pow(struct Fp2 *ANS, struct Fp2 *A, mpz_t a);
int Fp2_cmp(struct Fp2 *A, struct Fp2 *B);
int Fp2_cmp_ui(struct Fp2 *A, unsigned long int a);
int Fp2_cmp_mpz(struct Fp2 *A, mpz_t a);
int Fp2_cmp_zero(struct Fp2 *A);
int Fp2_cmp_one(struct Fp2 *A);
void Fp2_frobenius_1(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_frobenius_2(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_frobenius_3(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_frobenius_4(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_frobenius_6(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_frobenius_8(struct Fp2 *ANS, struct Fp2 *A);
void Fp2_frobenius_10(struct Fp2 *ANS, struct Fp2 *A);

#endif // FP2_H
