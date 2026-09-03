#ifndef KSS16_PARAMS_H
#define KSS16_PARAMS_H
#include "field.h"

extern char X_bit_binary[x_bit+1];
extern mpz_t X;
extern mpz_t PRIME_P,order_r,trace_t, order_EFp, a_x;
extern mpz_t tmp_a;

void KSS_16_parameters(void);
void generate_X();
void pre_calc_vector_final_exp(void);
void dealloc_constants();

#endif
