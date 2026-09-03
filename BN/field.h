#ifndef BN_FIELD_H
#define BN_FIELD_H
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

struct Fp{
    mpz_t x0;
};
struct Fp2{
    struct Fp x0;
    struct Fp x1;
};
struct Fp6{
    struct Fp2 x0;
    struct Fp2 x1;
    struct Fp2 x2;
};
struct Fp12{
    struct Fp6 x0;
    struct Fp6 x1;
};

void Fp_init(struct Fp *P);
void Fp_set(struct Fp *P,struct Fp *A);
void Fp_set_ui(struct Fp *P,unsigned long int a);
void Fp_set_mpz(struct Fp *P,mpz_t a);
void Fp_set_neg(struct Fp *P,struct Fp *A);
void Fp_random(struct Fp *P,gmp_randstate_t state);
void Fp_clear(struct Fp *P);
void Fp_printf(struct Fp *P,char *name);
void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B);
void Fp_mul_ui(struct Fp *ANS,struct Fp *A,unsigned long int a);
void Fp_mul_mpz(struct Fp *ANS,struct Fp *A,mpz_t a);
void Fp_add(struct Fp *ANS,struct Fp *A,struct Fp *B);
void Fp_add_ui(struct Fp *ANS,struct Fp *A,unsigned long int a);
void Fp_add_mpz(struct Fp *ANS,struct Fp *A,mpz_t a);
void Fp_sub(struct Fp *ANS,struct Fp *A,struct Fp *B);
void Fp_sub_ui(struct Fp *ANS,struct Fp *A,unsigned long int a);
void Fp_sub_mpz(struct Fp *ANS,struct Fp *A,mpz_t a);
void Fp_inv(struct Fp *ANS,struct Fp *A);
int Fp_legendre(struct Fp *A);
int Fp_isCNR(struct Fp *A);
void Fp_sqrt(struct Fp *ANS,struct Fp *A);
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t a);
int Fp_cmp(struct Fp *A,struct Fp *B);
int Fp_cmp_ui(struct Fp *A,unsigned long int a);
int Fp_cmp_mpz(struct Fp *A,mpz_t a);
int Fp_cmp_zero(struct Fp *A);
int Fp_cmp_one(struct Fp *A);
void Fp_frobenius_1(struct Fp *ANS,struct Fp *A);
void Fp_frobenius_2(struct Fp *ANS,struct Fp *A);
void Fp_frobenius_3(struct Fp *ANS,struct Fp *A);
void Fp_frobenius_4(struct Fp *ANS,struct Fp *A);
void Fp_frobenius_6(struct Fp *ANS,struct Fp *A);
void Fp_frobenius_8(struct Fp *ANS,struct Fp *A);
void Fp_frobenius_10(struct Fp *ANS,struct Fp *A);
void Fp2_init(struct Fp2 *P);
void Fp2_set(struct Fp2 *P,struct Fp2 *A);
void Fp2_set_ui(struct Fp2 *P,unsigned long int a);
void Fp2_set_mpz(struct Fp2 *P,mpz_t a);
void Fp2_set_neg(struct Fp2 *P,struct Fp2 *A);
void Fp2_random(struct Fp2 *P,gmp_randstate_t state);
void Fp2_clear(struct Fp2 *P);
void Fp2_printf(struct Fp2 *P,char *name);
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a);
void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
void Fp2_mul_basis(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_squaring(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_inv_basis(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a);
void Fp2_add_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_sub_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a);
void Fp2_sub_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
void Fp2_inv(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_inv_map(struct Fp2 *ANS,struct Fp2 *A);
int Fp2_legendre(struct Fp2 *A);
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t a);
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B);
int Fp2_cmp_ui(struct Fp2 *A,unsigned long int a);
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t a);
int Fp2_cmp_zero(struct Fp2 *A);
int Fp2_cmp_one(struct Fp2 *A);
void Fp2_frobenius_1(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_frobenius_2(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_frobenius_3(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_frobenius_4(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_frobenius_6(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_frobenius_8(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_frobenius_10(struct Fp2 *ANS,struct Fp2 *A);
void Fp6_init(struct Fp6 *P);
void Fp6_set(struct Fp6 *P,struct Fp6 *A);
void Fp6_set_ui(struct Fp6 *P,unsigned long int a);
void Fp6_set_mpz(struct Fp6 *P,mpz_t a);
void Fp6_set_neg(struct Fp6 *P,struct Fp6 *A);
void Fp6_random(struct Fp6 *P,gmp_randstate_t state);
void Fp6_clear(struct Fp6 *P);
void Fp6_printf(struct Fp6 *P,char *name);
void Fp6_mul(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_mul_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int a);
void Fp6_mul_mpz(struct Fp6 *ANS,struct Fp6 *A,mpz_t a);
void Fp6_mul_basis(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_squaring(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_add(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_add_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int a);
void Fp6_add_mpz(struct Fp6 *ANS,struct Fp6 *A,mpz_t a);
void Fp6_sub(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B);
void Fp6_sub_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int a);
void Fp6_sub_mpz(struct Fp6 *ANS,struct Fp6 *A,mpz_t a);
void Fp6_inv(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_inv_map_1(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_inv_map_2(struct Fp6 *ANS,struct Fp6 *A);
int Fp6_legendre(struct Fp6 *A);
void Fp6_sqrt(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_pow(struct Fp6 *ANS,struct Fp6 *A,mpz_t a);
int Fp6_cmp(struct Fp6 *A,struct Fp6 *B);
int Fp6_cmp_ui(struct Fp6 *A,unsigned long int a);
int Fp6_cmp_mpz(struct Fp6 *A,mpz_t a);
int Fp6_cmp_zero(struct Fp6 *A);
int Fp6_cmp_one(struct Fp6 *A);
void Fp6_frobenius_1(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_frobenius_2(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_frobenius_3(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_frobenius_4(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_frobenius_6(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_frobenius_8(struct Fp6 *ANS,struct Fp6 *A);
void Fp6_frobenius_10(struct Fp6 *ANS,struct Fp6 *A);
void Fp12_init(struct Fp12 *P);
void Fp12_set(struct Fp12 *P,struct Fp12 *A);
void Fp12_set_ui(struct Fp12 *P,unsigned long int a);
void Fp12_set_mpz(struct Fp12 *P,mpz_t a);
void Fp12_set_neg(struct Fp12 *P,struct Fp12 *A);
void Fp12_random(struct Fp12 *P,gmp_randstate_t state);
void Fp12_clear(struct Fp12 *P);
void Fp12_printf(struct Fp12 *P,char *name);
void Fp12_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_mul_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a);
void Fp12_mul_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
void Fp12_squaring(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_add(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_add_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a);
void Fp12_add_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
void Fp12_sub(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void Fp12_sub_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a);
void Fp12_sub_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
void Fp12_inv(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_inv_map(struct Fp12 *ANS,struct Fp12 *A);
int Fp12_legendre(struct Fp12 *A);
void Fp12_sqrt(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_pow(struct Fp12 *ANS,struct Fp12 *A,mpz_t a);
void Fp12_G3_pow(struct Fp12 *ANS,struct Fp12 *A,mpz_t S);
int Fp12_cmp(struct Fp12 *A,struct Fp12 *B);
int Fp12_cmp_ui(struct Fp12 *A,unsigned long int a);
int Fp12_cmp_mpz(struct Fp12 *A,mpz_t a);
int Fp12_cmp_zero(struct Fp12 *A);
int Fp12_cmp_one(struct Fp12 *A);
void Fp12_frobenius_1(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_2(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_3(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_4(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_6(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_8(struct Fp12 *ANS,struct Fp12 *A);
void Fp12_frobenius_10(struct Fp12 *ANS,struct Fp12 *A);
void Fp_mul_basis(struct Fp *ANS,struct Fp *A);
int Fp2_isCNR(struct Fp2 *A);
int Fp6_isCNR(struct Fp6 *A);

#endif
