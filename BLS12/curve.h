#ifndef BLS12_CURVE_H
#define BLS12_CURVE_H
#include "field.h"

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
void EFp12_to_EFp2(struct EFp2 *ANS,struct EFp12 *P);
void EFp2_to_EFp12(struct EFp12 *ANS,struct EFp2 *P);
void EFp_skew_frobenius(struct EFp *ANS,struct EFp *A);
void EFp12_G2_SCM_normal(struct EFp12 *ANS,struct EFp12 *Q,mpz_t S);

#endif
