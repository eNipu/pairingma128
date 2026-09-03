#ifndef KSS16_FIELD_H
#define KSS16_FIELD_H
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#define TRUE 1
#define FALSE 0
#define x_bit 35
#define c1 2

struct Fp{
    mpz_t x0;
};
struct Fp2{
    struct Fp x0,x1;
};
struct Fp4{
    struct Fp2 x0,x1;
};
struct Fp8{
    struct Fp4 x0,x1;
};
struct Fp16{
    struct Fp8 x0,x1;
};

extern struct Fp4 beta,beta_inv;
extern struct Fp C1,C1_INV;
extern struct Fp p_C4, p_M_C4;
extern struct Fp p_C8, p_M_C_C8;
extern struct Fp p_C_C4_C8, p_M_C_C4_C8;
extern struct Fp p_C16, p_M_C_C16;
extern struct Fp p_C4_C16, p_C_C4_C16, p_M_C_C4_C16;
extern struct Fp p_C4_C8_C16, p_C_C4_C8_C16, p_M_C_C_C4_C8_C16;
extern struct Fp p_C8_C16, p_C_C8_C16, p_M_C_C8_C16;
extern mpz_t p8p1dr;
extern struct Fp p3_C4;
extern struct Fp p3_M_C4;
extern struct Fp p3_M_C_C8;
extern struct Fp p3_C8;
extern struct Fp p3_M_C_C4_C8;
extern struct Fp p3_C8_C4;
extern struct Fp p3_M_C_C4_C16;
extern struct Fp p3_C4_C16;
extern struct Fp p3_C16;
extern struct Fp p3_M_C16;
extern struct Fp p3_C_C4_C8_C16;
extern struct Fp p3_M_C_C4_C8_C16;
extern struct Fp p3_M_C_C8_C16;
extern struct Fp p3_C8_C16;
extern struct Fp p5_C4, p5_M_C4;
extern struct Fp p5_C8, p5_M_C8;
extern struct Fp p5_C_C4_C8, p5_M_C_C4_C8;
extern struct Fp p5_C16, p5_M_C_C16;
extern struct Fp p5_C4_C16, p5_C_C4_C16, p5_M_C_C4_C16;
extern struct Fp p5_C4_C8_C16, p5_C_C4_C8_C16, p5_M_C_C_C4_C8_C16;
extern struct Fp p5_C8_C16, p5_C_C8_C16, p5_M_C_C8_C16;
extern struct Fp p7_C4;
extern struct Fp p7_M_C4;
extern struct Fp p7_M_C_C8;
extern struct Fp p7_C8;
extern struct Fp p7_M_C_C4_C8;
extern struct Fp p7_C8_C4;
extern struct Fp p7_M_C_C4_C16;
extern struct Fp p7_C4_C16;
extern struct Fp p7_C16;
extern struct Fp p7_M_C16;
extern struct Fp p7_C_C4_C8_C16;
extern struct Fp p7_M_C_C4_C8_C16;
extern struct Fp p7_M_C_C8_C16;
extern struct Fp p7_C8_C16;
extern struct Fp p4_C4;
extern struct Fp p4_C8;
extern struct Fp p4_C16;
extern struct Fp p4_C4_C8, p4_C8_C16;
extern struct Fp p4_C4_C16, p4_C4_C8_C16;
extern struct Fp p8_C4;
extern struct Fp p8_C8;
extern struct Fp p8_C16;
extern struct Fp p8_C4_C8, p8_C8_C16;
extern struct Fp p8_C4_C16, p8_C4_C8_C16;
extern struct Fp p2_C8, p2_M_C8;
extern struct Fp p2_C16, p2_M_C16, p2_C_C16, p2_M_C_C16;
extern struct Fp p2_C8_C16, p2_C_C8_C16, p2_M_C8_C16, p2_M_C_C8_C16;
extern struct Fp p6_C8, p6_M_C8;
extern struct Fp p6_C16, p6_M_C16, p6_C_C16, p6_M_C_C16;
extern struct Fp p6_C_C8_C16, p6_M_C_C8_C16, p6_C8_C16, p6_M_C8_C16;
extern struct Fp4 z_inv2;
extern struct Fp z_inv2_test;
extern mpz_t x2;
extern mpz_t x3;
extern mpz_t m00;
extern mpz_t m11;
extern mpz_t m22;
extern mpz_t m33;
extern mpz_t m44;
extern mpz_t m55;
extern mpz_t m66;
extern mpz_t m77;
extern mpz_t tmp1;

void Fp4_mul_betainv(struct Fp4 *ANS);
void Fp_init(struct Fp *A);
void Fp_set(struct Fp *ANS,struct Fp *A);
void Fp_set_ui(struct Fp *A,signed long int B);
void Fp_set_mpz(struct Fp *A, mpz_t a);
void Fp_random(struct Fp *A);
void Fp_clear(struct Fp *A);
void Fp_printf(struct Fp *A);
void Fp_add(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a+b mod p
void Fp_add_ui(struct Fp *ans,struct Fp *a,unsigned long int b);//ans=a+b mod p
void Fp_add_mpz(struct Fp *ans,struct Fp *a,mpz_t b);//ans=a+b mod p
void Fp_sub(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a-b mod p
void Fp_sub_ui(struct Fp *ans,struct Fp *a,unsigned long int b);//ans=a+b mod p
void Fp_sqr(struct Fp *ANS,struct Fp *A);
void Fp_mul(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a*b mod p
void Fp_mul_mpz(struct Fp *ans,struct Fp *a,mpz_t b);
void Fp_mul_ui(struct Fp *ans,struct Fp *a,unsigned long int b);//ans=a*b mod p
void Fp_mul_basis(struct Fp *ans,struct Fp *a);
void Fp_div(struct Fp *ans,struct Fp *a,struct Fp *b);//ans=a/b mod p
void Fp_pow(struct Fp *ans,struct Fp *a,mpz_t b);
void Fp_sqrt(struct Fp *ans,struct Fp *a);//x^2=a mod p
int  Fp_cmp_mpz(struct Fp *A,mpz_t B);
void Fp_mul_mpz(struct Fp *ANS,struct Fp *A,mpz_t B);
void Fp_neg(struct Fp *ANS,struct Fp *A);
int  Fp_cmp(struct Fp *A,struct Fp *B);
int Fp_cmp_mpz(struct Fp *A,mpz_t B);
void Fp2_init(struct Fp2 *A);
void Fp2_set(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_set_ui(struct Fp2 *A,signed long int B);
void Fp2_random(struct Fp2 *A);
void Fp2_clear(struct Fp2 *A);
void Fp2_printf(struct Fp2 *A);
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B);
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_sqr(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_sqr_complex(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_mul_i(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B);
void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t B);
void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B);
void Fp2_neg(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_invert(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B);
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B);
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A);//x^2=a mod p
int  Fp2_legendre(struct Fp2 *A);
int  Fp2_cmp(struct Fp2 *A,struct Fp2 *B);
int  Fp2_cmp_mpz(struct Fp2 *A,mpz_t B);
void Fp2_neg(struct Fp2 *ANS,struct Fp2 *A);
void Fp2_frobenius_map(struct Fp2 *ANS, struct Fp2 *A);
void Fp4_init(struct Fp4 *A);
void Fp4_set(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_set_ui(struct Fp4 *A,signed long int B);
void Fp4_random(struct Fp4 *A);
void Fp4_clear(struct Fp4 *A);
void Fp4_printf(struct Fp4 *A);
void Fp4_add(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B);
void Fp4_add_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B);
void Fp4_sub(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B);
void Fp4_sqr(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_sqr_complex(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_mul(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B);
void Fp4_mul_v(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_mul_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B);
void Fp4_mul_mpz(struct Fp4 *ANS,struct Fp4 *A,mpz_t B);
void Fp4_mul_Fp(struct Fp4 *ANS,struct Fp4 *A,struct Fp *B);
void Fp4_neg(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_invert(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_div(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B);
void Fp4_pow(struct Fp4 *ANS,struct Fp4 *A,mpz_t B);
void Fp4_sqrt(struct Fp4 *ANS,struct Fp4 *A);//x^2=a mod p
int  Fp4_legendre(struct Fp4 *A);
int  Fp4_cmp(struct Fp4 *A,struct Fp4 *B);
int  Fp4_cmp_mpz(struct Fp4 *A,mpz_t B);
void Fp4_neg(struct Fp4 *ANS,struct Fp4 *A);
void Fp4_frobenius_map(struct Fp4 *ANS, struct Fp4 *A);
void Fp8_init(struct Fp8 *A);
void Fp8_set(struct Fp8 *ANS,struct Fp8 *A);
void Fp8_set_ui(struct Fp8 *A,signed long int B);
void Fp8_random(struct Fp8 *A);
void Fp8_clear(struct Fp8 *A);
void Fp8_printf(struct Fp8 *A);
void Fp8_add(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
void Fp8_add_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B);
void Fp8_sub(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
void Fp8_mul(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
void Fp8_sqr(struct Fp8 *ANS,struct Fp8 *A);
void Fp8_sqr_complex(struct Fp8 *ANS,struct Fp8 *A);
void Fp8_mul_v(struct Fp8 *ANS,struct Fp8 *A);
void Fp8_mul_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B);
void Fp8_mul_mpz(struct Fp8 *ANS,struct Fp8 *A,mpz_t B);
void Fp8_mul_Fp(struct Fp8 *ANS,struct Fp8 *A,struct Fp *B);
void Fp8_neg(struct Fp8 *ANS,struct Fp8 *A);
void Fp8_invert(struct Fp8 *ANS,struct Fp8 *A);
void Fp8_div(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B);
void Fp8_pow(struct Fp8 *ANS,struct Fp8 *A,mpz_t B);
void Fp8_sqrt(struct Fp8 *ANS,struct Fp8 *A);//x^2=a mod p
int  Fp8_legendre(struct Fp8 *A);
int  Fp8_cmp(struct Fp8 *A,struct Fp8 *B);
int  Fp8_cmp_mpz(struct Fp8 *A,mpz_t B);
void Fp8_frobenius_map(struct Fp8 *ANS, struct Fp8 *A);
void Fp16_init(struct Fp16 *A);
void Fp16_set(struct Fp16 *ANS,struct Fp16 *A);
void Fp16_set_ui(struct Fp16 *A,signed long int B);
void Fp16_random(struct Fp16 *A);
void Fp16_clear(struct Fp16 *A);
void Fp16_printf(struct Fp16 *A);
void Fp16_add(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B);
void Fp16_add_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B);
void Fp16_sub(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B);
void Fp16_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B);
void Fp16_mul_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B);
void Fp16_mul_mpz(struct Fp16 *ANS,struct Fp16 *A,mpz_t B);
void Fp16_mul_Fp(struct Fp16 *ANS,struct Fp16 *A,struct Fp *B);
void Fp16_neg(struct Fp16 *ANS,struct Fp16 *A);
void Fp16_invert(struct Fp16 *ANS,struct Fp16 *A);
void Fp16_div(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B);
void Fp16_pow(struct Fp16 *ANS,struct Fp16 *A,mpz_t B);
void Fp16_sqrt(struct Fp16 *ANS,struct Fp16 *A);//x^2=a mod p
void Fp16_sqr(struct Fp16 *ANS,struct Fp16 *A);
void Fp16_sqr_complex(struct Fp16 *ANS,struct Fp16 *A);
int  Fp16_legendre(struct Fp16 *A);
int  Fp16_cmp(struct Fp16 *A,struct Fp16 *B);
int  Fp16_cmp_mpz(struct Fp16 *A,mpz_t B);
void Fp16_neg(struct Fp16 *ANS,struct Fp16 *A);
void Fp16_frobenius_map(struct Fp16 *ANS, struct Fp16 *A);
void Fp16_frobenius_map_p3(struct Fp16 *ANS, struct Fp16 *A);
void Fp16_frobenius_map_p5(struct Fp16 *ANS, struct Fp16 *A);
void Fp16_frobenius_map_p7(struct Fp16 *ANS, struct Fp16 *A);
void Fp16_frobenius_map_p2(struct Fp16 *ANS, struct Fp16 *A);
void Fp16_frobenius_map_p4(struct Fp16 *ANS, struct Fp16 *A);
void Fp16_frobenius_map_p6(struct Fp16 *ANS, struct Fp16 *A);
void Fp16_frobenius_map_p8(struct Fp16 *ANS, struct Fp16 *A);
void pre_calculate_frob_p();
void pre_calculate_frob_p2();
void pre_calculate_frob_p3();
void pre_calculate_frob_p4();
void pre_calculate_frob_p5();
void pre_calculate_frob_p6();
void pre_calculate_frob_p7();
void pre_calculate_frob_p8();;

#endif
