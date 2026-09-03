#ifndef KSS16_PAIRING_H
#define KSS16_PAIRING_H
#include "curve.h"
#include "params.h"

void final_exp_hard(struct Fp16 *ANS,struct Fp16 *A);
void Miller_algo(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q,mpz_t roop);
void Optimal_Miller(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q,mpz_t roop);
void Final_Exp(struct Fp16 *ANS,struct Fp16 *A);
void Tate_Pairing(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q);
void Ate_Pairing(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q);
void Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2);
void ltt_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q);
void v2t_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q);
void ltp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q);
void vtp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q);
void ADD_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q,struct Fp16 *Qx_neg);
void DBL_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *Q,struct Fp16 *Qx_neg);
void Pseudo_Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2);
void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2);
void Pseudo_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Pseudo_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Pseudo_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp4 *Q,struct EFp *P,struct Fp *L);
void Pseudo_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp *P,struct Fp *L);
void Pseudo_type1_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B);
void Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q);
void Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2);
void Sparse_type1_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Sparse_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop);
void Sparse_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *Px_neg);
void Sparse_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *Px_neg);

#endif
