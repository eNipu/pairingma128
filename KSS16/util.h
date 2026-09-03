#ifndef KSS16_UTIL_H
#define KSS16_UTIL_H
#include <sys/time.h>

extern int count_eca_new, count_ecd_new,count_eca_BIN, count_ecd_BIN,count_eca_BIN_map, count_ecd_BIN_map, count_eca_WIN, count_ecd_WIN, count_eca_WIN_map,count_ecd_WIN_map, count_eca_NAF, count_ecd_NAF,count_eca_NAF_map,count_ecd_NAF_map, count_eca_ML, count_ecd_ML,count_eca_ML_map, count_ecd_ML_map;
extern int count_eca_new_pre, count_ecd_new_pre;
extern int Big_M, Small_m, Big_add, Sqr;
extern long int myNAF[5000];
extern long int naf_index;
extern struct timeval tonaf_t0;
extern struct timeval tonaf_t1;
extern float elapsed_tonaf_t0;
extern float elapsed_tonaf_t1;
extern int fp16_pow;
extern unsigned long int num,mpz_mpz_mul,mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num;

float timedifference_msec(struct timeval t0, struct timeval t1);

#endif
