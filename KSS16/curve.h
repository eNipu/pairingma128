#ifndef KSS16_CURVE_H
#define KSS16_CURVE_H
#include "field.h"

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

void EFp4_SCM_ML(struct EFp4 *ANS, struct EFp4 *P,mpz_t j);
void EFp16_SCM_ML(struct EFp16 *ANS, struct EFp16 *P,mpz_t j);
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
void EFp4_init(struct EFp4 *A);
void EFp4_set(struct EFp4 *A,struct EFp4 *B);
void EFp4_set_infity(struct EFp4 *A);
void EFp4_set_EFp(struct EFp4 *ANS,struct EFp *A);
void EFp4_clear(struct EFp4 *A);
void EFp4_printf(struct EFp4 *A);
void EFp4_ECD(struct EFp4 *ANS, struct EFp4 *P);//ANS=2*P
void EFp4_ECA(struct EFp4 *ANS, struct EFp4 *P1, struct EFp4 *P2);//ANS=P1+P2
int  EFp4_cmp(struct EFp4 *A,struct EFp4 *B);
void EFp4_SCM_BIN(struct EFp4 *ANS, struct EFp4 *P, mpz_t j);
void EFp4_SCM_WIN(struct EFp4 *ANS, struct EFp4 *P, mpz_t scalar);
void EFp4_neg(struct EFp4 *ANS, struct EFp4 *A);
void EFp4_SCM_BIN_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j);
void EFp4_SCM_BIN_Pseudo_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j);
void EFp4_ECD_Pseudo_Sparse(struct EFp4 *ANS, struct EFp4 *P);
void EFp4_ECD_Sparse(struct EFp4 *ANS, struct EFp4 *P);
void EFp8_init(struct EFp8 *A);
void EFp8_set(struct EFp8 *A,struct EFp8 *B);
void EFp8_set_infity(struct EFp8 *A);
void EFp8_clear(struct EFp8 *A);
void EFp8_printf(struct EFp8 *A);
void EFp8_ECD(struct EFp8 *ANS, struct EFp8 *P);//ANS=2*P
void EFp8_ECA(struct EFp8 *ANS, struct EFp8 *P1, struct EFp8 *P2);//ANS=P1+P2
int  EFp8_cmp(struct EFp8 *A,struct EFp8 *B);
void EFp8_SCM_BIN(struct EFp8 *ANS, struct EFp8 *P, mpz_t j);
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
void EFp16_to_EFp4_map(struct EFp4 *ANS,struct EFp16 *A);
void EFp4_to_EFp16_map(struct EFp16 *ANS,struct EFp4 *A);
void Skew_Frobenius_map(struct EFp4 *ANS, struct EFp4 *Qt);

#endif
