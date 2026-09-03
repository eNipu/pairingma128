#ifndef BN_TEST_H
#define BN_TEST_H
#include "pairing.h"

void test();
void test_G1_SCM();
void test_G2_SCM();
void test_G3_pow();
void test_tate_pairing();
void test_ate_pairing();
void test_twist();
void test_frobenius();
void test_skew_frobenius();
void check_num_of_Fp_mul();

#endif
