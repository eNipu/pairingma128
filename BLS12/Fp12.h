#ifndef FP12_H
#define FP12_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include "Fp6.h"

struct Fp12 {
    struct Fp6 x0;
    struct Fp6 x1;
};

void Fp12_init(struct Fp12 *P);
void Fp12_set(struct Fp12 *P, struct Fp12 *A);
void Fp12_set_ui(struct Fp12 *P, unsigned long int a);
void Fp12_set_mpz(struct Fp12 *P, mpz_t a);
void Fp12_set_neg(struct Fp12 *P, struct Fp12 *A);
void Fp12_random(struct Fp12 *P, gmp_randstate_t state);
void Fp12_clear(struct Fp12 *P);
void Fp12_printf(struct Fp12 *P, char *name);
void Fp12_mul(struct Fp12 *ANS, struct Fp12 *A, struct Fp12 *B);
void Fp12_mul_ui(struct Fp12 *ANS, struct Fp12 *A, unsigned long int a);
void Fp12_mul_mpz(struct Fp12 *ANS, struct Fp12 *A, mpz_t a);
void Fp12_squaring(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_add(struct Fp12 *ANS, struct Fp12 *A, struct Fp12 *B);
void Fp12_add_ui(struct Fp12 *ANS, struct Fp12 *A, unsigned long int a);
void Fp12_add_mpz(struct Fp12 *ANS, struct Fp12 *A, mpz_t a);
void Fp12_sub(struct Fp12 *ANS, struct Fp12 *A, struct Fp12 *B);
void Fp12_sub_ui(struct Fp12 *ANS, struct Fp12 *A, unsigned long int a);
void Fp12_sub_mpz(struct Fp12 *ANS, struct Fp12 *A, mpz_t a);
void Fp12_inv(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_inv_map(struct Fp12 *ANS, struct Fp12 *A);
int Fp12_legendre(struct Fp12 *A);
void Fp12_sqrt(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_pow(struct Fp12 *ANS, struct Fp12 *A, mpz_t a);
int Fp12_cmp(struct Fp12 *A, struct Fp12 *B);
int Fp12_cmp_ui(struct Fp12 *A, unsigned long int a);
int Fp12_cmp_mpz(struct Fp12 *A, mpz_t a);
int Fp12_cmp_zero(struct Fp12 *A);
int Fp12_cmp_one(struct Fp12 *A);
void Fp12_frobenius_1(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_frobenius_2(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_frobenius_3(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_frobenius_4(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_frobenius_6(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_frobenius_8(struct Fp12 *ANS, struct Fp12 *A);
void Fp12_frobenius_10(struct Fp12 *ANS, struct Fp12 *A);

#endif // FP12_H
