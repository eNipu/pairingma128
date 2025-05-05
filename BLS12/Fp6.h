#ifndef FP6_H
#define FP6_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include "Fp2.h"

struct Fp6 {
    struct Fp2 x0;
    struct Fp2 x1;
    struct Fp2 x2;
};

void Fp6_init(struct Fp6 *P);
void Fp6_set(struct Fp6 *P, struct Fp6 *A);
void Fp6_set_ui(struct Fp6 *P, unsigned long int a);
void Fp6_set_mpz(struct Fp6 *P, mpz_t a);
void Fp6_set_neg(struct Fp6 *P, struct Fp6 *A);
void Fp6_random(struct Fp6 *P, gmp_randstate_t state);
void Fp6_clear(struct Fp6 *P);
void Fp6_printf(struct Fp6 *P, char *name);
void Fp6_mul(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_mul_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int a);
void Fp6_mul_mpz(struct Fp6 *ANS, struct Fp6 *A, mpz_t a);
void Fp6_mul_basis(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_squaring(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_add(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_add_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int a);
void Fp6_add_mpz(struct Fp6 *ANS, struct Fp6 *A, mpz_t a);
void Fp6_sub(struct Fp6 *ANS, struct Fp6 *A, struct Fp6 *B);
void Fp6_sub_ui(struct Fp6 *ANS, struct Fp6 *A, unsigned long int a);
void Fp6_sub_mpz(struct Fp6 *ANS, struct Fp6 *A, mpz_t a);
void Fp6_inv(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_inv_map_1(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_inv_map_2(struct Fp6 *ANS, struct Fp6 *A);
int Fp6_legendre(struct Fp6 *A);
void Fp6_sqrt(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_pow(struct Fp6 *ANS, struct Fp6 *A, mpz_t a);
int Fp6_cmp(struct Fp6 *A, struct Fp6 *B);
int Fp6_cmp_ui(struct Fp6 *A, unsigned long int a);
int Fp6_cmp_mpz(struct Fp6 *A, mpz_t a);
int Fp6_cmp_zero(struct Fp6 *A);
int Fp6_cmp_one(struct Fp6 *A);
void Fp6_frobenius_1(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_frobenius_2(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_frobenius_3(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_frobenius_4(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_frobenius_6(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_frobenius_8(struct Fp6 *ANS, struct Fp6 *A);
void Fp6_frobenius_10(struct Fp6 *ANS, struct Fp6 *A);

#endif // FP6_H
