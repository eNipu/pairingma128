#ifndef BLS12_PARAMS_H
#define BLS12_PARAMS_H
#include "field.h"

#define x_bit 77
extern int X_bit_binary[x_bit+1];

extern mpz_t prime;
extern mpz_t mother_parameter;
extern int sign;
extern mpz_t trace_t;
extern mpz_t EFp_order;
extern mpz_t EFp_total;
extern mpz_t EFp2_total;
extern mpz_t EFp6_total;
extern mpz_t EFp12_total;
extern mpz_t curve_parameter_A;
extern mpz_t curve_parameter_B;
extern mpz_t final_exp;
extern struct Fp Fp_basis;
extern struct Fp2 Fp2_basis;
extern struct Fp2 Fp2_basis_inv;
extern struct Fp6 Fp6_basis;
//ZERO
extern struct Fp Fp_ZERO;
extern struct Fp2 Fp2_ZERO;
extern struct Fp6 Fp6_ZERO;
extern struct Fp12 Fp12_ZERO;
//ONE
extern struct Fp Fp_ONE;
extern struct Fp2 Fp2_ONE;
extern struct Fp6 Fp6_ONE;
extern struct Fp12 Fp12_ONE;
//precomputed
extern struct Fp inv_CNR1;
extern struct Fp inv_CNR2;
extern struct Fp epsilon_1;
extern struct Fp epsilon_2;
//frobenius
extern struct Fp2 Fp2_basis_prime_1_div_3_1;
extern struct Fp2 Fp2_basis_prime_1_div_3_2;
extern struct Fp2 Fp2_basis_prime_1_div_6;
extern struct Fp2 Fp2_basis_prime_2_div_3_1;
extern struct Fp2 Fp2_basis_prime_2_div_3_2;
extern struct Fp2 Fp2_basis_prime_2_div_6;
extern struct Fp2 Fp2_basis_prime_3_div_3_1;
extern struct Fp2 Fp2_basis_prime_3_div_3_2;
extern struct Fp2 Fp2_basis_prime_3_div_6;
extern struct Fp2 Fp2_basis_prime_4_div_3_1;
extern struct Fp2 Fp2_basis_prime_4_div_3_2;
extern struct Fp2 Fp2_basis_prime_4_div_6;
extern struct Fp2 Fp2_basis_prime_8_div_3_1;
extern struct Fp2 Fp2_basis_prime_8_div_3_2;
extern struct Fp2 Fp2_basis_prime_8_div_6;
extern struct Fp2 Fp2_basis_prime_10_div_3_1;
extern struct Fp2 Fp2_basis_prime_10_div_3_2;
extern struct Fp2 Fp2_basis_prime_10_div_6;
//skew frobenius
extern struct Fp2 Fp2_basis_inv_prime_1_div_3;
extern struct Fp2 Fp2_basis_inv_prime_1_div_2;
extern struct Fp2 Fp2_basis_inv_prime_2_div_3;
extern struct Fp2 Fp2_basis_inv_prime_2_div_2;
extern struct Fp2 Fp2_basis_inv_prime_3_div_3;
extern struct Fp2 Fp2_basis_inv_prime_3_div_2;
extern struct Fp2 Fp2_basis_inv_prime_10_div_3;
extern struct Fp2 Fp2_basis_inv_prime_10_div_2;

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

#endif
