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

char X_bit_binary[x_bit+1];

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
//=============================//
struct EFp{
    struct Fp x,y;
    int infity;
};
struct EFp2{
    struct Fp2 x,y;
    int infity;
};
struct EFp4{
    struct Fp4 x,y;
    int infity;
};
struct EFp8{
    struct Fp8 x,y;
    int infity;
};
struct EFp16{
    struct Fp16 x,y;
    int infity;
};

mpz_t X;
mpz_t PRIME_P,order_r,trace_t, order_EFp, a_x;
mpz_t tmp_a;

struct Fp4 beta,beta_inv;

void Fp4_mul_betainv(struct Fp4 *ANS);
void EFp4_SCM_ML(struct EFp4 *ANS, struct EFp4 *P,mpz_t j);
void EFp16_SCM_ML(struct EFp16 *ANS, struct EFp16 *P,mpz_t j);


// #pragma mark Fp methods
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
//-----------------------------------------------------------------------------------------
// #pragma mark Fp2 methods
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

//-----------------------------------------------------------------------------------------
// #pragma mark Fp4 methods
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

//-----------------------------------------------------------------------------------------
// #pragma mark Fp8 methods
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

//-----------------------------------------------------------------------------------------
// #pragma mark Fp16 methods
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
//-----------------------------------------------------------------------------------------
// #pragma mark EFp methods
void EFp_init(struct EFp *A);
void EFp_set(struct EFp *A,struct EFp *B);
void EFp_set_infity(struct EFp *A);
void EFp_clear(struct EFp *A);
void EFp_printf(struct EFp *A);
void EFp_SCM_BIN(struct EFp *ANS, struct EFp *P,mpz_t j);
void EFp_SCM_NAF(struct EFp *ANS, struct EFp *P, mpz_t scalar);
void EFp_SCM_WIN(struct EFp *ANS, struct EFp *P, mpz_t scalar);
void EFp_ECD(struct EFp *ANS, struct EFp *P);//ANS=2*P
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2);//ANS=P1+P2
int  EFp_cmp(struct EFp *A,struct EFp *B);
void EFp_random_set(struct EFp *ANS);//random set EFp on curve
void EFp_neg(struct EFp *ANS, struct EFp *A);

//-----------------------------------------------------------------------------------------
// #pragma mark EFp2 methods
void EFp2_init(struct EFp2 *A);
void EFp2_set(struct EFp2 *A,struct EFp2 *B);
void EFp2_set_infity(struct EFp2 *A);
void EFp2_clear(struct EFp2 *A);
void EFp2_printf(struct EFp2 *A);
void EFp2_ECD(struct EFp2 *ANS, struct EFp2 *P);//ANS=2*P
void EFp2_ECA(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2);//ANS=P1+P2
int  EFp2_cmp(struct EFp2 *A,struct EFp2 *B);
void EFp2_set_EFp(struct EFp2 *A,struct EFp *B);
void EFp2_random_set(struct EFp2 *ANS);
void EFp2_SCM_BIN(struct EFp2 *ANS, struct EFp2 *P, mpz_t j);
void EFp2_SCM_NAF(struct EFp2 *ANS, struct EFp2 *P, mpz_t scalar);
void EFp2_SCM_WIN(struct EFp2 *ANS, struct EFp2 *P, mpz_t scalar);
void EFp2_neg(struct EFp2 *ANS, struct EFp2 *A);

//-----------------------------------------------------------------------------------------
// #pragma mark EFp4 methods
void EFp4_init(struct EFp4 *A);
void EFp4_set(struct EFp4 *A,struct EFp4 *B);
void EFp4_set_infity(struct EFp4 *A);
void EFp4_set_EFp(struct EFp4 *ANS,struct EFp *A);
void EFp4_clear(struct EFp4 *A);
void EFp4_printf(struct EFp4 *A);
void EFp4_ECD(struct EFp4 *ANS, struct EFp4 *P);//ANS=2*P
void EFp4_ECA(struct EFp4 *ANS, struct EFp4 *P1, struct EFp4 *P2);//ANS=P1+P2
int  EFp4_cmp(struct EFp4 *A,struct EFp4 *B);
//void EFp4_random_set(struct EFp4 *ANS);
void EFp4_SCM_BIN(struct EFp4 *ANS, struct EFp4 *P, mpz_t j);
//void EFp4_SCM_NAF(struct EFp4 *ANS, struct EFp4 *P, mpz_t scalar);
void EFp4_SCM_WIN(struct EFp4 *ANS, struct EFp4 *P, mpz_t scalar);
void EFp4_neg(struct EFp4 *ANS, struct EFp4 *A);
void EFp4_SCM_BIN_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j);
void EFp4_SCM_BIN_Pseudo_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j);
void EFp4_ECD_Pseudo_Sparse(struct EFp4 *ANS, struct EFp4 *P);
void EFp4_ECD_Sparse(struct EFp4 *ANS, struct EFp4 *P);
//-----------------------------------------------------------------------------------------
// #pragma mark EFp8 methods
void EFp8_init(struct EFp8 *A);
void EFp8_set(struct EFp8 *A,struct EFp8 *B);
void EFp8_set_infity(struct EFp8 *A);
void EFp8_clear(struct EFp8 *A);
void EFp8_printf(struct EFp8 *A);
void EFp8_ECD(struct EFp8 *ANS, struct EFp8 *P);//ANS=2*P
void EFp8_ECA(struct EFp8 *ANS, struct EFp8 *P1, struct EFp8 *P2);//ANS=P1+P2
int  EFp8_cmp(struct EFp8 *A,struct EFp8 *B);
//void EFp8_random_set(struct EFp8 *ANS);
void EFp8_SCM_BIN(struct EFp8 *ANS, struct EFp8 *P, mpz_t j);
//void EFp8_neg(struct EFp8 *ANS, struct EFp8 *A);

//-----------------------------------------------------------------------------------------
// #pragma mark EFp16 methods
void EFp16_init(struct EFp16 *A);
void EFp16_set(struct EFp16 *A,struct EFp16 *B);
void EFp16_set_infity(struct EFp16 *A);
void EFp16_set_EFp(struct EFp16 *A,struct EFp *B);
void EFp16_clear(struct EFp16 *A);
void EFp16_printf(struct EFp16 *A);
void EFp16_ECD(struct EFp16 *ANS, struct EFp16 *P);//ANS=2*P
void EFp16_ECA(struct EFp16 *ANS, struct EFp16 *P1, struct EFp16 *P2);//ANS=P1+P2
int  EFp16_cmp(struct EFp16 *A,struct EFp16 *B);
void EFp16_SCM_BIN(struct EFp16 *ANS, struct EFp16 *P, mpz_t j);
void EFp16_SCM_WIN(struct EFp16 *ANS, struct EFp16 *P, mpz_t scalar);
void EFp16_frobenius_map(struct EFp16 *ANS,struct EFp16 *A);
void EFp16_random_set(struct EFp16 *ANS);
void EFp16_random_set_G2(struct EFp16 *ANS);
void final_exp_hard(struct Fp16 *ANS,struct Fp16 *A);
//-----------------------------------------------------------------------------------------

#pragma mark util methods
void KSS_16_parameters(void);
void EFp16_to_EFp4_map(struct EFp4 *ANS,struct EFp16 *A);
void EFp4_to_EFp16_map(struct EFp16 *ANS,struct EFp4 *A);
void Check_SCM();
void Skew_Frobenius_map(struct EFp4 *ANS, struct EFp4 *Qt);
void generate_X();
void pre_calc_vector_final_exp(void);
void pre_calculate_frob_p();
void pre_calculate_frob_p2();
void pre_calculate_frob_p3();
void pre_calculate_frob_p4();
void pre_calculate_frob_p5();
void pre_calculate_frob_p6();
void pre_calculate_frob_p7();
void pre_calculate_frob_p8();;
void dealloc_constants();
float timedifference_msec(struct timeval t0, struct timeval t1);
void check_Pairing(void);

#pragma mark Pairing methods
void Miller_algo(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q,mpz_t roop);
void Optimal_Miller(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q,mpz_t roop);
void Final_Exp(struct Fp16 *ANS,struct Fp16 *A);
void Tate_Pairing(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q);
void Ate_Pairing(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q);
void Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2);
//void Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2);
void ltt_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q);
void v2t_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q);
void ltp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q);
void vtp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q);
void ADD_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q,struct Fp16 *Qx_neg);
void DBL_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *Q,struct Fp16 *Qx_neg);
#pragma mark Pseudo Sparse
void Pseudo_Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2);
void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2);
void Pseudo_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Pseudo_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Pseudo_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp4 *Q,struct EFp *P,struct Fp *L);
void Pseudo_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp *P,struct Fp *L);
void Pseudo_type1_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B);

#pragma mark Pseudo Sparse
void Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q);
void Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2);

void Sparse_type1_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Sparse_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Sparse_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *Px_neg);
void Sparse_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *Px_neg);

void rational_point_check(struct EFp4 *A);
