#include "util.h"

unsigned long int num,mpz_mpz_mul,mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num;

float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}

float timedifference_usec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec);
}
