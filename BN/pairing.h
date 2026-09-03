#ifndef BN_PAIRING_H
#define BN_PAIRING_H
#include "curve.h"
#include "params.h"

void Final_exp_1(struct Fp12 *ANS,struct Fp12 *A);
void Final_exp_2(struct Fp12 *ANS,struct Fp12 *A);
void Final_exp_normal(struct Fp12 *ANS,struct Fp12 *A);
void EFp_ECD_return_lambda(struct EFp *ANS,struct Fp *lambda,struct EFp *P);
void EFp_ECA_return_lambda(struct EFp *ANS,struct Fp *lambda,struct EFp *P1,struct EFp *P2);
void Pseudo_8_sparse_mapping(struct EFp *P,struct EFp2 *Q,struct Fp *L);
void Pseudo_8_sparse_mul_2(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void ff_ltt_2(struct Fp12 *f,struct EFp2 *T,struct EFp *P,struct Fp *L);
void f_ltq_2(struct Fp12 *f,struct EFp2 *T,struct EFp2 *Q,struct EFp *P,struct Fp *L);
void Miller_algo_for_opt_ate_2(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
void Opt_ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);

#endif
