#ifndef BLS12_PAIRING_H
#define BLS12_PAIRING_H
#include "curve.h"
#include "params.h"

void Final_exp_female_researchers_algo(struct Fp12 *ANS,struct Fp12 *A);
void Final_exp_normal(struct Fp12 *ANS,struct Fp12 *A);
void Pseudo_8_sparse_mapping(struct EFp *P,struct EFp2 *Q,struct Fp *L);
void Pseudo_8_sparse_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B);
void ff_ltt(struct Fp12 *f,struct EFp2 *T,struct EFp *P,struct Fp *L);
void f_ltq(struct Fp12 *f,struct EFp2 *T,struct EFp2 *Q,struct EFp *P,struct Fp *L);
void Miller_algo_for_ate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);		//**
void Ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);			//**
void Miller_algo_for_opt_ate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
void Opt_ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
void ff_ltt_vtt(struct Fp12 *f,struct EFp12 *T,struct EFp12 *Q);
void f_ltp_vtp(struct Fp12 *f,struct EFp12 *T,struct EFp12 *P,struct EFp12 *Q);
void Miller_algo_for_tate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);
void Tate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P);

#endif
