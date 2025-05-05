#ifndef FP_H
#define FP_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

struct Fp {
    mpz_t x0;
};

void Fp_init(struct Fp *P);
void Fp_set(struct Fp *P, struct Fp *A);
void Fp_set_ui(struct Fp *P, unsigned long int a);
void Fp_set_mpz(struct Fp *P, mpz_t a);
void Fp_set_neg(struct Fp *P, struct Fp *A);
void Fp_random(struct Fp *P, gmp_randstate_t state);
void Fp_clear(struct Fp *P);
void Fp_printf(struct Fp *P, char *name);
void Fp_mul(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_mul_ui(struct Fp *ANS, struct Fp *A, unsigned long int a);
void Fp_mul_mpz(struct Fp *ANS, struct Fp *A, mpz_t a);
void Fp_add(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_add_ui(struct Fp *ANS, struct Fp *A, unsigned long int a);
void Fp_add_mpz(struct Fp *ANS, struct Fp *A, mpz_t a);
void Fp_sub(struct Fp *ANS, struct Fp *A, struct Fp *B);
void Fp_sub_ui(struct Fp *ANS, struct Fp *A, unsigned long int a);
void Fp_sub_mpz(struct Fp *ANS, struct Fp *A, mpz_t a);
void Fp_inv(struct Fp *ANS, struct Fp *A);
int Fp_legendre(struct Fp *A);
int Fp_isCNR(struct Fp *A);
void Fp_sqrt(struct Fp *ANS, struct Fp *A);
void Fp_pow(struct Fp *ANS, struct Fp *A, mpz_t a);
int Fp_cmp(struct Fp *A, struct Fp *B);
int Fp_cmp_ui(struct Fp *A, unsigned long int a);
int Fp_cmp_mpz(struct Fp *A, mpz_t a);
int Fp_cmp_zero(struct Fp *A);
int Fp_cmp_one(struct Fp *A);
void Fp_frobenius_1(struct Fp *ANS, struct Fp *A);
void Fp_frobenius_2(struct Fp *ANS, struct Fp *A);
void Fp_frobenius_3(struct Fp *ANS, struct Fp *A);
void Fp_frobenius_4(struct Fp *ANS, struct Fp *A);
void Fp_frobenius_6(struct Fp *ANS, struct Fp *A);
void Fp_frobenius_8(struct Fp *ANS, struct Fp *A);
void Fp_frobenius_10(struct Fp *ANS, struct Fp *A);

#endif // FP_H
