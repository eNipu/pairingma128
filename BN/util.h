#ifndef BN_UTIL_H
#define BN_UTIL_H
#include <sys/time.h>

extern unsigned long int num,mpz_mpz_mul,mpz_ui_mul,Fp_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num;

float timedifference_msec(struct timeval t0, struct timeval t1);
float timedifference_usec(struct timeval t0, struct timeval t1);

#endif
