#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#define x_bit 77
int X_bit_binary[x_bit+1];

mpz_t prime;
mpz_t mother_parameter;
int sign;
mpz_t trace_t;
mpz_t EFp_order;
mpz_t EFp_total;
mpz_t EFp2_total;
mpz_t EFp6_total;
mpz_t EFp12_total;
mpz_t curve_parameter_A;
mpz_t curve_parameter_B;
mpz_t final_exp;
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
struct EFp{
    struct Fp x;
    struct Fp y;
    int flag;
};
struct EFp2{
    struct Fp2 x;
    struct Fp2 y;
    int flag;
};
struct EFp6{
    struct Fp6 x;
    struct Fp6 y;
    int flag;
};
struct EFp12{
    struct Fp12 x;
    struct Fp12 y;
    int flag;
};
struct Fp Fp_basis;
struct Fp2 Fp2_basis;
struct Fp2 Fp2_basis_inv;
struct Fp6 Fp6_basis;
//ZERO
struct Fp Fp_ZERO;
struct Fp2 Fp2_ZERO;
struct Fp6 Fp6_ZERO;
struct Fp12 Fp12_ZERO;
//ONE
struct Fp Fp_ONE;
struct Fp2 Fp2_ONE;
struct Fp6 Fp6_ONE;
struct Fp12 Fp12_ONE;
//precomputed
struct Fp inv_CNR1;
struct Fp inv_CNR2;
struct Fp epsilon_1;
struct Fp epsilon_2;
//frobenius
struct Fp2 Fp2_basis_prime_1_div_3_1;
struct Fp2 Fp2_basis_prime_1_div_3_2;
struct Fp2 Fp2_basis_prime_1_div_6;
struct Fp2 Fp2_basis_prime_2_div_3_1;
struct Fp2 Fp2_basis_prime_2_div_3_2;
struct Fp2 Fp2_basis_prime_2_div_6;
struct Fp2 Fp2_basis_prime_3_div_3_1;
struct Fp2 Fp2_basis_prime_3_div_3_2;
struct Fp2 Fp2_basis_prime_3_div_6;
struct Fp2 Fp2_basis_prime_4_div_3_1;
struct Fp2 Fp2_basis_prime_4_div_3_2;
struct Fp2 Fp2_basis_prime_4_div_6;
struct Fp2 Fp2_basis_prime_8_div_3_1;
struct Fp2 Fp2_basis_prime_8_div_3_2;
struct Fp2 Fp2_basis_prime_8_div_6;
struct Fp2 Fp2_basis_prime_10_div_3_1;
struct Fp2 Fp2_basis_prime_10_div_3_2;
struct Fp2 Fp2_basis_prime_10_div_6;
//skew frobenius
struct Fp2 Fp2_basis_inv_prime_1_div_3;
struct Fp2 Fp2_basis_inv_prime_1_div_2;
struct Fp2 Fp2_basis_inv_prime_2_div_3;
struct Fp2 Fp2_basis_inv_prime_2_div_2;
struct Fp2 Fp2_basis_inv_prime_3_div_3;
struct Fp2 Fp2_basis_inv_prime_3_div_2;
struct Fp2 Fp2_basis_inv_prime_10_div_3;
struct Fp2 Fp2_basis_inv_prime_10_div_2;

/*============================================================================*/
/* functions                                                                  */
/*============================================================================*/
/*-----------------------------------Fp-----------------------------------*/
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

/*-----------------------------------Fp2----------------------------------*/
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

/*-----------------------------------Fp6----------------------------------*/
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

/*-----------------------------------Fp12---------------------------------*/
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
void Fp12_G3_exp_2div(struct Fp12 *ANS,struct Fp12 *A,mpz_t S);				//****
void Fp12_G3_exp_4div(struct Fp12 *ANS,struct Fp12 *A,mpz_t S);				//****
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

/*-----------------------------------EFp----------------------------------*/
void EFp_init(struct EFp *P);
void EFp_set(struct EFp *P,struct EFp *A);
void EFp_set_ui(struct EFp *P,unsigned long int a);
void EFp_set_mpz(struct EFp *P,mpz_t a);
void EFp_set_neg(struct EFp *P,struct EFp *A);
void EFp_clear(struct EFp *P);
void EFp_printf(struct EFp *P,char *name);
void EFp_rational_point(struct EFp *P);
void EFp_ECD(struct EFp *ANS,struct EFp *P);
void EFp_ECA(struct EFp *ANS,struct EFp *P1,struct EFp *P2);
void EFp_SCM(struct EFp *ANS,struct EFp *P,mpz_t R);
void EFp_skew_frobenius_2(struct EFp *ANS,struct EFp *A);					//*****

/*-----------------------------------EFp2---------------------------------*/
void EFp2_init(struct EFp2 *P);
void EFp2_set(struct EFp2 *P,struct EFp2 *A);
void EFp2_set_ui(struct EFp2 *P,unsigned long int a);
void EFp2_set_mpz(struct EFp2 *P,mpz_t a);
void EFp2_set_neg(struct EFp2 *ANS,struct EFp2 *P);
void EFp2_clear(struct EFp2 *P);
void EFp2_printf(struct EFp2 *P,char *name);
void EFp2_rational_point(struct EFp2 *P);
void EFp2_ECD(struct EFp2 *ANS,struct EFp2 *P);
void EFp2_ECA(struct EFp2 *ANS,struct EFp2 *P1,struct EFp2 *P2);
void EFp2_SCM(struct EFp2 *ANS,struct EFp2 *P,mpz_t R);
void EFp2_frobenius_1(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_frobenius_2(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_frobenius_3(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_frobenius_4(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_frobenius_6(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_frobenius_8(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_frobenius_10(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_1(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_2(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_3(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_4(struct EFp2 *ANS,struct EFp2 *A);
void EFp2_skew_frobenius_10(struct EFp2 *ANS,struct EFp2 *A);

/*-----------------------------------EFp6---------------------------------*/
void EFp6_init(struct EFp6 *P);
void EFp6_set(struct EFp6 *P,struct EFp6 *A);
void EFp6_set_ui(struct EFp6 *P,unsigned long int a);
void EFp6_set_mpz(struct EFp6 *P,mpz_t a);
void EFp6_set_neg(struct EFp6 *P,struct EFp6 *A);
void EFp6_clear(struct EFp6 *P);
void EFp6_printf(struct EFp6 *P,char *name);
void EFp6_rational_point(struct EFp6 *P);
void EFp6_ECD(struct EFp6 *ANS,struct EFp6 *P);
void EFp6_ECA(struct EFp6 *ANS,struct EFp6 *P1,struct EFp6 *P2);
void EFp6_SCM(struct EFp6 *ANS,struct EFp6 *P,mpz_t R);
void EFp6_frobenius_1(struct EFp6 *ANS,struct EFp6 *A);
void EFp6_frobenius_2(struct EFp6 *ANS,struct EFp6 *A);
void EFp6_frobenius_3(struct EFp6 *ANS,struct EFp6 *A);
void EFp6_frobenius_4(struct EFp6 *ANS,struct EFp6 *A);
void EFp6_frobenius_6(struct EFp6 *ANS,struct EFp6 *A);
void EFp6_frobenius_8(struct EFp6 *ANS,struct EFp6 *A);
void EFp6_frobenius_10(struct EFp6 *ANS,struct EFp6 *A);

/*-----------------------------------EFp12--------------------------------*/
void EFp12_init(struct EFp12 *P);
void EFp12_set(struct EFp12 *P,struct EFp12 *A);
void EFp12_set_ui(struct EFp12 *P,unsigned long int a);
void EFp12_set_mpz(struct EFp12 *P,mpz_t a);
void EFp12_set_neg(struct EFp12 *P,struct EFp12 *A);
void EFp12_clear(struct EFp12 *P);
void EFp12_printf(struct EFp12 *P,char *name);
void EFp12_rational_point(struct EFp12 *P);
void EFp12_generate_G1(struct EFp12 *P);
void EFp12_generate_G2(struct EFp12 *Q);
void EFp12_ECD(struct EFp12 *ANS,struct EFp12 *P);
void EFp12_ECA(struct EFp12 *ANS,struct EFp12 *P1,struct EFp12 *P2);
void EFp12_ECD_G2optimal(struct EFp12 *ANS,struct EFp12 *P);
void EFp12_ECA_G2optimal(struct EFp12 *ANS,struct EFp12 *P1,struct EFp12 *P2);
void EFp12_SCM(struct EFp12 *ANS,struct EFp12 *P,mpz_t R);
void EFp12_frobenius_1(struct EFp12 *ANS,struct EFp12 *A);
void EFp12_frobenius_2(struct EFp12 *ANS,struct EFp12 *A);
void EFp12_frobenius_3(struct EFp12 *ANS,struct EFp12 *A);
void EFp12_frobenius_4(struct EFp12 *ANS,struct EFp12 *A);
void EFp12_frobenius_6(struct EFp12 *ANS,struct EFp12 *A);
void EFp12_frobenius_8(struct EFp12 *ANS,struct EFp12 *A);
void EFp12_frobenius_10(struct EFp12 *ANS,struct EFp12 *A);

/*---------------------------------twist-----------------------------------*/
void EFp12_to_EFp2(struct EFp2 *ANS,struct EFp12 *P);
void EFp2_to_EFp12(struct EFp12 *ANS,struct EFp2 *P);

/*-------------------------------final exp-----------------------------------*/
void Final_exp_female_researchers_algo(struct Fp12 *ANS,struct Fp12 *A);
void Final_exp_normal(struct Fp12 *ANS,struct Fp12 *A);


/*--------------------------------8 sparse-----------------------------------*/
void Pseudo_8_sparse_mapping(struct EFp *P,struct EFp2 *Q,struct Fp *L);
void Pseudo_8_sparse_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void ff_ltt(struct Fp12 *f,struct EFp2 *T,struct EFp *P,struct Fp *L);
void f_ltq(struct Fp12 *f,struct EFp2 *T,struct EFp2 *Q,struct EFp *P,struct Fp *L);

/*-------------------------------ate pairing---------------------------------*/
void Miller_algo_for_ate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);		//**
void Ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);			//**

/*----------------------------opt-ate pairing--------------------------------*/
void Miller_algo_for_opt_ate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
void Opt_ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);

/*-----------------------------------init---------------------------------*/
void init_parameters();
void set_parameters();
void clear_parameters();
void generate_mother_parameter();
int generate_prime();
int generate_order();
void generate_trace();
void generate_basis();
void get_epsilon();
void get_scalar_of_final_exp();
void weil();
void print_parameters();

/*-----------------------------------time---------------------------------*/
float timedifference_msec(struct timeval t0, struct timeval t1);
float timedifference_usec(struct timeval t0, struct timeval t1);

/*-----------------------------------test---------------------------------*/
void test();										//**(easy)
void test_opt_ate_pairing();
void test_frobenius();
void check_num_of_Fp_mul();
