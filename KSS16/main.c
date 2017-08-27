//
//  main.c
//  Fp16_EffSM
//
//  Created by ***** on 6/7/16.
//  Copyright © 2016 ******. All rights reserved.
//

#include"KSS_16.h"
#include <time.h>

int count_eca_new, count_ecd_new,count_eca_BIN, count_ecd_BIN,count_eca_BIN_map, count_ecd_BIN_map, count_eca_WIN, count_ecd_WIN, count_eca_WIN_map,count_ecd_WIN_map, count_eca_NAF, count_ecd_NAF,count_eca_NAF_map,count_ecd_NAF_map, count_eca_ML, count_ecd_ML,count_eca_ML_map, count_ecd_ML_map;

int count_eca_new_pre, count_ecd_new_pre;

int Big_M, Small_m, Big_add, Sqr;
long int myNAF[5000];
long int naf_index = 0;

struct timeval tonaf_t0;
struct timeval tonaf_t1;
float elapsed_tonaf_t0;
float elapsed_tonaf_t1;

struct Fp C1,C1_INV;
struct Fp p_C4, p_M_C4;
struct Fp p_C8, p_M_C_C8;
struct Fp p_C_C4_C8, p_M_C_C4_C8;
struct Fp p_C16, p_M_C_C16;
struct Fp p_C4_C16, p_C_C4_C16, p_M_C_C4_C16;
struct Fp p_C4_C8_C16, p_C_C4_C8_C16, p_M_C_C_C4_C8_C16;
struct Fp p_C8_C16, p_C_C8_C16, p_M_C_C8_C16;
mpz_t p8p1dr;

struct Fp p3_C4;
struct Fp p3_M_C4;
struct Fp p3_M_C_C8;
struct Fp p3_C8;
struct Fp p3_M_C_C4_C8;
struct Fp p3_C8_C4;
struct Fp p3_M_C_C4_C16;
struct Fp p3_C4_C16;
struct Fp p3_C16;
struct Fp p3_M_C16;
struct Fp p3_C_C4_C8_C16;
struct Fp p3_M_C_C4_C8_C16;
struct Fp p3_M_C_C8_C16;
struct Fp p3_C8_C16;
//
struct Fp p5_C4, p5_M_C4;
struct Fp p5_C8, p5_M_C8;
struct Fp p5_C_C4_C8, p5_M_C_C4_C8;
struct Fp p5_C16, p5_M_C_C16;
struct Fp p5_C4_C16, p5_C_C4_C16, p5_M_C_C4_C16;
struct Fp p5_C4_C8_C16, p5_C_C4_C8_C16, p5_M_C_C_C4_C8_C16;
struct Fp p5_C8_C16, p5_C_C8_C16, p5_M_C_C8_C16;
//

struct Fp p7_C4;
struct Fp p7_M_C4;
struct Fp p7_M_C_C8;
struct Fp p7_C8;
struct Fp p7_M_C_C4_C8;
struct Fp p7_C8_C4;
struct Fp p7_M_C_C4_C16;
struct Fp p7_C4_C16;
struct Fp p7_C16;
struct Fp p7_M_C16;
struct Fp p7_C_C4_C8_C16;
struct Fp p7_M_C_C4_C8_C16;
struct Fp p7_M_C_C8_C16;
struct Fp p7_C8_C16;

//
struct Fp p4_C4;
struct Fp p4_C8;
struct Fp p4_C16;
struct Fp p4_C4_C8, p4_C8_C16;
struct Fp p4_C4_C16, p4_C4_C8_C16;
//
struct Fp p8_C4;
struct Fp p8_C8;
struct Fp p8_C16;
struct Fp p8_C4_C8, p8_C8_C16;
struct Fp p8_C4_C16, p8_C4_C8_C16;
//
struct Fp p2_C8, p2_M_C8;
struct Fp p2_C16, p2_M_C16, p2_C_C16, p2_M_C_C16;
struct Fp p2_C8_C16, p2_C_C8_C16, p2_M_C8_C16, p2_M_C_C8_C16;

struct Fp p6_C8, p6_M_C8;
struct Fp p6_C16, p6_M_C16, p6_C_C16, p6_M_C_C16;
struct Fp p6_C_C8_C16, p6_M_C_C8_C16, p6_C8_C16, p6_M_C8_C16;

struct Fp4 z_inv2;
struct Fp z_inv2_test;

mpz_t x2;
mpz_t x3;
mpz_t m00;
mpz_t m11;
mpz_t m22;
mpz_t m33;
mpz_t m44;
mpz_t m55;
mpz_t m66;
mpz_t m77;
mpz_t tmp1;

int fp16_pow;

unsigned long int num,mpz_mpz_mul,mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num;

int main(void){
    mpz_init(X);
    
    generate_X();
    mpz_init(PRIME_P);
    mpz_init(order_r);
    mpz_init(order_EFp);
    mpz_init(trace_t);
    mpz_init(a_x);
    
    KSS_16_parameters();
    pre_calculate_frob_p();
    pre_calculate_frob_p2();
    pre_calculate_frob_p3();
    pre_calculate_frob_p4();
    pre_calculate_frob_p5();
    pre_calculate_frob_p6();
    pre_calculate_frob_p7();
    pre_calculate_frob_p8();
    pre_calc_vector_final_exp();
    
    gmp_printf("X = %Zd\n",X);
    gmp_printf("p = %Zd\n",PRIME_P);
    gmp_printf("r = %Zd\n",order_r);
    gmp_printf("t = %Zd\n",trace_t);
    gmp_printf("#E(Fp) = %Zd\n",order_EFp);
    
    printf("X = %dbit\n",(int)mpz_sizeinbase(X,2));
    printf("p = %dbit\n",(int)mpz_sizeinbase(PRIME_P,2));
    printf("r = %dbit\n",(int)mpz_sizeinbase(order_r,2));
    printf("t = %dbit\n",(int)mpz_sizeinbase(trace_t,2));
    gmp_printf("y^2 = x^3 + %Zdx\n",a_x);
    

    check_Pairing();
    
    mpz_clear(PRIME_P);
    mpz_clear(order_r);
    mpz_clear(order_EFp);
    mpz_clear(trace_t);
    mpz_clear(a_x);
    dealloc_constants();
    return 0;
}

void dealloc_constants()
{
    Fp_clear(&C1);
    Fp_init(&C1_INV);
    Fp_clear(&p_C16);
    Fp_clear(&p_C8);
    Fp_clear(&p_C4);
    Fp_clear(&p_M_C4);
    Fp_clear(&p_M_C_C8);
    Fp_clear(&p_M_C_C16);
    Fp_clear(&p_C8_C16);
    Fp_clear(&p_C_C8_C16);
    Fp_clear(&p_M_C_C8_C16);
    Fp_clear(&p_C4_C16);
    Fp_clear(&p_C_C4_C8);
    Fp_clear(&p_M_C_C4_C8);
    Fp_clear(&p_C4_C8_C16);
    Fp_clear(&p_M_C_C_C4_C8_C16);
    Fp_clear(&p_C_C4_C8_C16);
}

void pre_calculate_frob_p()
{
    mpz_init(p8p1dr);
    Fp_init(&C1);
    Fp_init(&p_C16);
    Fp_init(&p_C8);
    Fp_init(&p_C4);
    Fp_init(&p_M_C4);
    Fp_init(&p_M_C_C8);
    Fp_init(&p_M_C_C16);
    Fp_init(&p_C8_C16);
    Fp_init(&p_C_C8_C16);
    Fp_init(&p_M_C_C8_C16);
    Fp_init(&p_C4_C16);
    Fp_init(&p_C_C4_C8);
    Fp_init(&p_M_C_C4_C8);
    Fp_init(&p_C4_C8_C16);
    Fp_init(&p_M_C_C_C4_C8_C16);
    Fp_init(&p_C_C4_C8_C16);
    Fp_init(&p_C_C4_C16);
    Fp_init(&p_M_C_C4_C16);
    
    Fp_set_ui(&C1,c1);
    
    mpz_pow_ui(p8p1dr,PRIME_P,8);
    mpz_add_ui(p8p1dr,p8p1dr,1);
    mpz_tdiv_q(p8p1dr,p8p1dr,order_r);
    
    mpz_invert(C1_INV.x0,C1.x0,PRIME_P);
    Fp_inv_num++;
    
    mpz_sub_ui(p_C16.x0,PRIME_P,13);
    mpz_tdiv_q_ui(p_C16.x0,p_C16.x0,16);
    Fp_pow(&p_C16,&C1,p_C16.x0);
    
    mpz_sub_ui(p_C8.x0,PRIME_P,5);
    mpz_tdiv_q_ui(p_C8.x0,p_C8.x0,8);
    Fp_pow(&p_C8,&C1,p_C8.x0);
    
    mpz_sub_ui(p_C4.x0,PRIME_P,1);
    mpz_tdiv_q_ui(p_C4.x0,p_C4.x0,4);
    Fp_pow(&p_C4,&C1,p_C4.x0);
    Fp_neg(&p_M_C4, &p_C4);
    
    mpz_sub(p_M_C_C8.x0,PRIME_P,C1.x0);
    Fp_mul(&p_M_C_C8, &p_C8, &p_M_C_C8);
    
    Fp_mul(&p_M_C_C16, &p_C16, &C1);
    Fp_neg(&p_M_C_C16, &p_M_C_C16);
    
    Fp_mul(&p_C4_C16, &p_C4, &p_C8);
    Fp_mul(&p_C_C4_C8, &p_C4_C16, &C1); // c.c^(p-1)/4.c^(p-5)/8
    mpz_sub(p_M_C_C4_C8.x0,PRIME_P,p_C_C4_C8.x0);
    
    Fp_mul(&p_C8_C16,&p_C16,&p_C8);// c^(p-13)/16.c^(p-5)/8
    Fp_mul(&p_C_C8_C16,&p_C8_C16,&C1);// c.c^(p-13)/16.c^(p-5)/8
    mpz_sub(p_M_C_C8_C16.x0,PRIME_P,p_C_C8_C16.x0);
    
    Fp_mul(&p_C4_C8_C16, &p_C_C8_C16, &p_C4); // c.c^(p-1)/4.c^(p-13)/16.c^(p-5)/8
    
    Fp_mul(&p_C_C4_C8_C16, &p_C4_C8_C16, &C1); // c.c.c^(p-13)/16.c^(p-5)/8
    mpz_sub(p_M_C_C_C4_C8_C16.x0,PRIME_P,p_C_C4_C8_C16.x0);
    
    Fp_mul(&p_C_C4_C16, &p_C4, &p_C16);
    Fp_mul(&p_C_C4_C16, &p_C_C4_C16, &C1);
    Fp_neg(&p_M_C_C4_C16, &p_C_C4_C16);
}

void pre_calculate_frob_p2()
{
    Fp_init(&C1);
    Fp_init(&p2_C16);
    Fp_init(&p2_M_C16);
    Fp_init(&p2_C8);
    Fp_init(&p2_M_C8);
    Fp_init(&p2_C_C16);
    Fp_init(&p2_M_C_C16);
    Fp_init(&p2_C8_C16);
    Fp_init(&p2_M_C8_C16);
    Fp_init(&p2_C_C8_C16);
    Fp_init(&p2_M_C_C8_C16);
    
    
    Fp_set_ui(&C1,c1);
    
    mpz_invert(C1_INV.x0,C1.x0,PRIME_P);
    Fp_inv_num++;
    
    mpz_pow_ui(p2_C16.x0,PRIME_P,2);
    mpz_sub_ui(p2_C16.x0,p2_C16.x0,9);
    mpz_tdiv_q_ui(p2_C16.x0,p2_C16.x0,16);
    Fp_pow(&p2_C16,&C1,p2_C16.x0);
    Fp_neg(&p2_M_C16,&p2_C16);
    
    mpz_pow_ui(p2_C8.x0,PRIME_P,2);
    mpz_sub_ui(p2_C8.x0,p2_C8.x0,1);
    mpz_tdiv_q_ui(p2_C8.x0,p2_C8.x0,8);
    Fp_pow(&p2_C8,&C1,p2_C8.x0);
    Fp_neg(&p2_M_C8,&p2_C8);
    
    Fp_mul(&p2_C_C16,&p2_C16,&C1);//c1*c1^(p^2-9/16)
    Fp_neg(&p2_M_C_C16, &p2_C_C16);
    
    Fp_mul(&p2_C8_C16, &p2_C8, &p2_C16);
    Fp_mul(&p2_C_C8_C16, &p2_C8_C16, &C1);
    Fp_neg(&p2_M_C8_C16, &p2_C8_C16);
    Fp_neg(&p2_M_C_C8_C16, &p2_C_C8_C16);
}
void pre_calculate_frob_p3()
{
    struct Fp TMP;
    Fp_init(&TMP);
    Fp_init(&p3_C4);
    Fp_init(&p3_M_C4);
    Fp_init(&p3_M_C_C8);
    Fp_init(&p3_M_C_C4_C8);
    Fp_init(&p3_C8_C4);
    Fp_init(&p3_M_C_C4_C16);
    Fp_init(&p3_C4_C16);
    Fp_init(&p3_C16);
    Fp_init(&p3_M_C16);
    Fp_init(&p3_C4);
    Fp_init(&p3_C_C4_C8_C16);
    Fp_init(&p3_M_C_C4_C8_C16);
    Fp_init(&p3_M_C_C8_C16);
    Fp_init(&p3_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_3;
    mpz_init(prime_3);
    mpz_pow_ui(prime_3,PRIME_P,3);
    
    mpz_sub_ui(tmp,prime_3,5);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p3_C16,&C1,tmp);
    Fp_neg(&p3_M_C16, &p3_C16);
    
    mpz_sub_ui(tmp,prime_3,5);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p3_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_3,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p3_C4,&C1,tmp);
    Fp_printf(&p3_C4);
    
    Fp_neg(&p3_M_C4, &p3_C4);
    
    Fp_mul(&p3_M_C_C8, &p3_C8, &C1);
    Fp_neg(&p3_M_C_C8, &p3_M_C_C8);
    ;
    Fp_mul(&p3_M_C_C4_C8, &p3_M_C_C8, &p3_C4);
    
    Fp_mul(&p3_C8_C4, &p3_C8, &p3_C4);
    
    Fp_mul(&p3_C8_C16, &p3_C8, &p3_C16);
    Fp_mul(&TMP, &p3_C8_C16, &C1);
    Fp_neg(&p3_M_C_C8_C16, &TMP);
    
    Fp_mul(&p3_C_C4_C8_C16, &TMP, &p3_C4);
    Fp_neg(&p3_M_C_C4_C8_C16, &p3_C_C4_C8_C16);
    
    Fp_mul(&p3_C4_C16, &p3_C4, &p3_C16);
    Fp_mul(&p3_M_C_C4_C16, &p3_C4_C16, &C1);
    Fp_neg(&p3_M_C_C4_C16, &p3_M_C_C4_C16);
    
    Fp_clear(&TMP);
    mpz_clear(tmp);
    mpz_clear(prime_3);
}

void pre_calculate_frob_p4()
{
    Fp_init(&p4_C4);
    Fp_init(&p4_C8);
    Fp_init(&p4_C16);
    Fp_init(&p4_C4_C8);
    Fp_init(&p4_C8_C16);
    Fp_init(&p4_C4_C16);
    Fp_init(&p4_C4_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_4;
    mpz_init(prime_4);
    mpz_pow_ui(prime_4,PRIME_P,4);
    
    mpz_sub_ui(tmp,prime_4,1);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p4_C16,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_4,1);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p4_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_4,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p4_C4,&C1,tmp);
    
    Fp_mul(&p4_C4_C8, &p4_C4, &p4_C8);
    Fp_mul(&p4_C4_C8_C16, &p4_C4_C8, &p4_C16);
    
    Fp_mul(&p4_C8_C16, &p4_C16, &p4_C8);
    Fp_mul(&p4_C4_C16, &p4_C16, &p4_C4);
    
    mpz_clear(tmp);
    mpz_clear(prime_4);
}

void pre_calculate_frob_p8()
{
    Fp_init(&p8_C4);
    Fp_init(&p8_C8);
    Fp_init(&p8_C16);
    Fp_init(&p8_C4_C8);
    Fp_init(&p8_C8_C16);
    Fp_init(&p8_C4_C16);
    Fp_init(&p8_C4_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_8;
    mpz_init(prime_8);
    mpz_pow_ui(prime_8,PRIME_P,8);
    
    mpz_sub_ui(tmp,prime_8,1);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p8_C16,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_8,1);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p8_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_8,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p8_C4,&C1,tmp);
    
    Fp_mul(&p8_C4_C8, &p8_C4, &p8_C8);
    Fp_mul(&p8_C4_C8_C16, &p8_C4_C8, &p8_C16);
    
    Fp_mul(&p8_C8_C16, &p8_C16, &p8_C8);
    Fp_mul(&p8_C4_C16, &p8_C16, &p8_C4);
    
    mpz_clear(tmp);
    mpz_clear(prime_8);
}


void pre_calculate_frob_p5()
{
    Fp_init(&p5_C16);
    Fp_init(&p5_C8);
    Fp_init(&p5_C4);
    Fp_init(&p5_M_C4);
    Fp_init(&p5_M_C8);
    Fp_init(&p5_M_C_C16);
    Fp_init(&p5_C8_C16);
    Fp_init(&p5_C_C8_C16);
    Fp_init(&p5_M_C_C8_C16);
    Fp_init(&p5_C4_C16);
    Fp_init(&p5_C_C4_C8);
    Fp_init(&p5_M_C_C4_C8);
    Fp_init(&p5_C4_C8_C16);
    Fp_init(&p5_M_C_C_C4_C8_C16);
    Fp_init(&p5_C_C4_C8_C16);
    Fp_init(&p5_C_C4_C16);
    Fp_init(&p5_M_C_C4_C16);
    
    mpz_t prime_5;
    mpz_init(prime_5);
    mpz_pow_ui(prime_5,PRIME_P,5);
    
    mpz_sub_ui(p5_C16.x0,prime_5,13);
    mpz_tdiv_q_ui(p5_C16.x0,p5_C16.x0,16);
    Fp_pow(&p5_C16,&C1,p5_C16.x0);
    
    mpz_sub_ui(p5_C8.x0,prime_5,5);
    mpz_tdiv_q_ui(p5_C8.x0,p5_C8.x0,8);
    Fp_pow(&p5_C8,&C1,p5_C8.x0);
    
    mpz_sub_ui(p5_C4.x0,prime_5,1);
    mpz_tdiv_q_ui(p5_C4.x0,p5_C4.x0,4);
    Fp_pow(&p5_C4,&C1,p5_C4.x0);
    Fp_neg(&p5_M_C4, &p5_C4);
    
    mpz_sub(p5_M_C8.x0,PRIME_P,C1.x0);
    Fp_mul(&p5_M_C8, &p5_C8, &p5_M_C8);
    
    Fp_mul(&p5_M_C_C16, &p5_C16, &C1);
    Fp_neg(&p5_M_C_C16, &p5_M_C_C16);
    
    Fp_mul(&p5_C4_C16, &p5_C4, &p5_C8);
    Fp_mul(&p5_C_C4_C8, &p5_C4_C16, &C1); // c.c^(p-1)/4.c^(p-5)/8
    mpz_sub(p5_M_C_C4_C8.x0,PRIME_P,p5_C_C4_C8.x0);
    
    Fp_mul(&p5_C8_C16,&p5_C16,&p5_C8);// c^(p-13)/16.c^(p-5)/8
    Fp_mul(&p5_C_C8_C16,&p5_C8_C16,&C1);// c.c^(p-13)/16.c^(p-5)/8
    mpz_sub(p5_M_C_C8_C16.x0,PRIME_P,p5_C_C8_C16.x0);
    
    Fp_mul(&p5_C4_C8_C16, &p5_C_C8_C16, &p5_C4); // c.c^(p-1)/4.c^(p-13)/16.c^(p-5)/8
    
    Fp_mul(&p5_C_C4_C8_C16, &p5_C4_C8_C16, &C1); // c.c.c^(p-13)/16.c^(p-5)/8
    Fp_neg(&p5_M_C_C_C4_C8_C16, &p5_C_C4_C8_C16);
    
    Fp_mul(&p5_C_C4_C16, &p5_C4, &p5_C16);
    Fp_mul(&p5_C_C4_C16, &p5_C_C4_C16, &C1);
    Fp_neg(&p5_M_C_C4_C16, &p5_C_C4_C16);
}
void pre_calculate_frob_p6()
{
    Fp_init(&C1);
    Fp_init(&p6_C16);
    Fp_init(&p6_M_C16);
    Fp_init(&p6_C8);
    Fp_init(&p6_M_C8);
    Fp_init(&p6_C_C16);
    Fp_init(&p6_M_C_C16);
    Fp_init(&p6_C_C8_C16);
    Fp_init(&p6_M_C_C8_C16);
    Fp_init(&p6_C8_C16);
    Fp_init(&p6_M_C8_C16);
    
    Fp_set_ui(&C1,c1);
    
    mpz_pow_ui(p6_C16.x0,PRIME_P,6);
    mpz_sub_ui(p6_C16.x0,p6_C16.x0,9);
    mpz_tdiv_q_ui(p6_C16.x0,p6_C16.x0,16);
    Fp_pow(&p6_C16,&C1,p6_C16.x0);
    Fp_neg(&p6_M_C16,&p6_C16);
    
    mpz_pow_ui(p6_C8.x0,PRIME_P,6);
    mpz_sub_ui(p6_C8.x0,p6_C8.x0,1);
    mpz_tdiv_q_ui(p6_C8.x0,p6_C8.x0,8);
    Fp_pow(&p6_C8,&C1,p6_C8.x0);
    Fp_neg(&p6_M_C8,&p6_C8);
    
    Fp_mul(&p6_C_C16,&p6_C16,&C1);//c1*c1^(p^2-9/16)
    Fp_neg(&p6_M_C_C16, &p6_C_C16);
    
    Fp_mul(&p6_C8_C16, &p6_C8, &p6_C16);
    Fp_mul(&p6_C_C8_C16, &p6_C8_C16, &C1);
    Fp_neg(&p6_M_C8_C16, &p6_C8_C16);
    Fp_neg(&p6_M_C_C8_C16, &p6_C_C8_C16);
}
void pre_calculate_frob_p7()
{
    struct Fp TMP;
    Fp_init(&TMP);
    Fp_init(&p7_C4);
    Fp_init(&p7_M_C4);
    Fp_init(&p7_M_C_C8);
    Fp_init(&p7_M_C_C4_C8);
    Fp_init(&p7_C8_C4);
    Fp_init(&p7_M_C_C4_C16);
    Fp_init(&p7_C4_C16);
    Fp_init(&p7_C16);
    Fp_init(&p7_M_C16);
    Fp_init(&p7_C4);
    Fp_init(&p7_C_C4_C8_C16);
    Fp_init(&p7_M_C_C4_C8_C16);
    Fp_init(&p7_M_C_C8_C16);
    Fp_init(&p7_C8_C16);
    
    mpz_t tmp;
    mpz_init(tmp);
    
    mpz_t prime_7;
    mpz_init(prime_7);
    mpz_pow_ui(prime_7,PRIME_P,7);
    
    mpz_sub_ui(tmp,prime_7,5);
    mpz_tdiv_q_ui(tmp,tmp,16);
    Fp_pow(&p7_C16,&C1,tmp);
    Fp_neg(&p7_M_C16, &p7_C16);
    
    mpz_sub_ui(tmp,prime_7,5);
    mpz_tdiv_q_ui(tmp,tmp,8);
    Fp_pow(&p7_C8,&C1,tmp);
    
    mpz_sub_ui(tmp,prime_7,1);
    mpz_tdiv_q_ui(tmp,tmp,4);
    Fp_pow(&p7_C4,&C1,tmp);
    Fp_printf(&p7_C4);
    
    Fp_neg(&p7_M_C4, &p7_C4);
    
    Fp_mul(&p7_M_C_C8, &p7_C8, &C1);
    Fp_neg(&p7_M_C_C8, &p7_M_C_C8);
    ;
    Fp_mul(&p7_M_C_C4_C8, &p7_M_C_C8, &p7_C4);
    
    Fp_mul(&p7_C8_C4, &p7_C8, &p7_C4);
    
    Fp_mul(&p7_C8_C16, &p7_C8, &p7_C16);
    Fp_mul(&TMP, &p7_C8_C16, &C1);
    Fp_neg(&p7_M_C_C8_C16, &TMP);
    
    Fp_mul(&p7_C_C4_C8_C16, &TMP, &p7_C4);
    Fp_neg(&p7_M_C_C4_C8_C16, &p7_C_C4_C8_C16);
    
    Fp_mul(&p7_C4_C16, &p7_C4, &p7_C16);
    Fp_mul(&p7_M_C_C4_C16, &p7_C4_C16, &C1);
    Fp_neg(&p7_M_C_C4_C16, &p7_M_C_C4_C16);
    
    Fp_clear(&TMP);
    mpz_clear(tmp);
    mpz_clear(prime_7);
}

void generate_X(){
    //c1 = 2
    // 2^ -2^32-2^18+2^8+1
    X_bit_binary[35]=1;
    X_bit_binary[32]=-1;
    X_bit_binary[18]=-1;
    X_bit_binary[8]=1;
    X_bit_binary[0]=1;
    //2^49+2^26+2^15-2^7-1
    //    X_bit_binary[49]=1;
    //    X_bit_binary[26]=1;
    //    X_bit_binary[15]=1;
    //    X_bit_binary[7]=-1;
    //    X_bit_binary[0]=-1;
    
    mpz_t tmp,set_2;
    mpz_init(tmp);
    mpz_init(set_2);
    mpz_set_ui(set_2,2);
    
    int i;
    for(i=x_bit;i>=0;i--){
        printf("%d",X_bit_binary[i]);
        if(X_bit_binary[i]==1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_add(X,X,tmp);
        }else if(X_bit_binary[i]==-1){
            mpz_pow_ui(tmp,set_2,i);
            mpz_sub(X,X,tmp);
        }
    }
    printf("\n");
    mpz_out_str(stdout,10,X);
    printf("\n");
    return;
}
//mpz_mpz_mul,mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num;
//-----------------------------------------------------------------------------------------
// #pragma mark Fp methods
void Fp_init(struct Fp *A){
    mpz_init(A->x0);
}
void Fp_set(struct Fp *ANS,struct Fp *E){
    mpz_set(ANS->x0,E->x0);
}
void Fp_set_ui(struct Fp *A,signed long int B){
    mpz_set_ui(A->x0,B);
}
void Fp_set_mpz(struct Fp *A, mpz_t a)
{
    mpz_set(A->x0,a);
}
void Fp_random(struct Fp *A){
    mpz_random(A->x0,10);
    mpz_mod(A->x0,A->x0,PRIME_P);
}
void Fp_clear(struct Fp *A){
    mpz_clear(A->x0);
}
void Fp_printf(struct Fp *A){
    gmp_printf("%Zd\n",A->x0);
}
void Fp_add(struct Fp *ANS,struct Fp *A,struct Fp *B){
    mpz_mpz_add++;
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_add_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
    mpz_ui_add++;
    mpz_add_ui(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_add_mpz(struct Fp *ANS,struct Fp *A,mpz_t B){
    mpz_add(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_sub(struct Fp *ANS,struct Fp *A,struct Fp *B){
    mpz_mpz_add++;
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_sub_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
    mpz_ui_add++;
    mpz_sub_ui(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    
}
void Fp_sqr(struct Fp *ANS,struct Fp *A){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_mul(tmp.x0,A->x0,A->x0);//A^2: GMP mpz_mul take care of squaring
    Sqr++;
    mpz_mod(tmp.x0,tmp.x0,PRIME_P);
    
    Fp_set(ANS,&tmp);
    Fp_clear(&tmp);
}
void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        Fp_mpz_sqr++;
    }else{
        mpz_mpz_mul++;
    }
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_mul_mpz(struct Fp *ANS,struct Fp *A,mpz_t B){
    mpz_mpz_mul++;
    mpz_mul(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_mul_ui(struct Fp *ANS,struct Fp *A,unsigned long int B){
    mpz_ui_mul++;
    mpz_mul_ui(ANS->x0,A->x0,B);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_mul_basis(struct Fp *ANS,struct Fp *A){
    basis_mul_num++;
    mpz_mul_ui(ANS->x0,A->x0,c1);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
}
void Fp_div(struct Fp *ANS,struct Fp *A,struct Fp *B){
    struct Fp tmp;
    Fp_init(&tmp);
    
    mpz_invert(tmp.x0,B->x0,PRIME_P);
    Fp_inv_num++;
    mpz_mul(ANS->x0,A->x0,tmp.x0);
    mpz_mod(ANS->x0,ANS->x0,PRIME_P);
    
    Fp_clear(&tmp);
}
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t j){
    int i,length;
    length=(int)mpz_sizeinbase(j,2);
    char binary[length];
    mpz_get_str(binary,2,j);
    struct Fp buf;
    Fp_init(&buf);
    Fp_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp_mul(&buf,&buf,&buf);
        if(binary[i]=='1'){
            Fp_mul(&buf,A,&buf);
        }
    }
    Fp_set(ANS,&buf);
    
    Fp_clear(&buf);
}
void Fp_sqrt(struct Fp *ANS,struct Fp *A){
    struct Fp n_tmp,y_tmp,x_tmp,b_tmp,t_tmp,tmp_Fp;
    Fp_init(&n_tmp);
    Fp_init(&y_tmp);
    Fp_init(&x_tmp);
    Fp_init(&b_tmp);
    Fp_init(&t_tmp);
    Fp_init(&tmp_Fp);
    
    Fp_set(&n_tmp,A);
    
    mpz_t tmp_mpz,q_tmp,e_tmp,r_tmp,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q_tmp);
    mpz_init(e_tmp);
    mpz_init(r_tmp);
    mpz_init(set_1);
    mpz_init(set_2);
    
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(mpz_legendre(n_tmp.x0,PRIME_P)!=-1){
        Fp_add_ui(&n_tmp,&n_tmp,1);
    }
    
    mpz_set(q_tmp,PRIME_P);
    mpz_sub_ui(q_tmp,q_tmp,1);
    mpz_set_ui(e_tmp,0);
    
    while(mpz_odd_p(q_tmp)==0){
        mpz_add_ui(e_tmp,e_tmp,1);
        mpz_div_ui(q_tmp,q_tmp,2);
    }
    
    Fp_pow(&y_tmp,&n_tmp,q_tmp);
    
    mpz_set(r_tmp,e_tmp);
    
    mpz_sub_ui(tmp_mpz,q_tmp,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    Fp_pow(&x_tmp,A,tmp_mpz);
    Fp_pow(&tmp_Fp,&x_tmp,set_2);
    Fp_mul(&b_tmp,&tmp_Fp,A);
    Fp_mul(&x_tmp,&x_tmp,A);
    
    int m;
    
    while(Fp_cmp_mpz(&b_tmp,set_1)==1){
        m=-1;
        while(Fp_cmp_mpz(&tmp_Fp,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp_pow(&tmp_Fp,&b_tmp,tmp_mpz);
        }
        //        gmp_printf("%Zd\n",tmp_Fp.x0);
        mpz_sub_ui(tmp_mpz,r_tmp,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        Fp_pow(&t_tmp,&y_tmp,tmp_mpz);
        Fp_pow(&y_tmp,&t_tmp,set_2);
        mpz_set_ui(r_tmp,m);
        Fp_mul(&x_tmp,&x_tmp,&t_tmp);
        Fp_mul(&b_tmp,&b_tmp,&y_tmp);
        Fp_set(&tmp_Fp, &b_tmp);
    }
    
    
    
    Fp_set(ANS,&x_tmp);
    
    Fp_clear(&n_tmp);
    Fp_clear(&y_tmp);
    Fp_clear(&x_tmp);
    Fp_clear(&b_tmp);
    Fp_clear(&t_tmp);
    Fp_clear(&tmp_Fp);
    mpz_clear(tmp_mpz);
    mpz_clear(q_tmp);
    mpz_clear(e_tmp);
    mpz_clear(r_tmp);
    mpz_clear(set_1);
}
void Fp_neg(struct Fp *ANS,struct Fp *A){
    mpz_sub(ANS->x0,PRIME_P,A->x0);
}
int Fp_cmp(struct Fp *A,struct Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        return 0;
    }
    return 1;
}
int Fp_cmp_mpz(struct Fp *A,mpz_t B){
    if(mpz_cmp(A->x0,B)==0){
        return 0;
    }
    return 1;
}


// #pragma mark Fp2 Methods
void Fp2_init(struct Fp2 *A){
    Fp_init(&A->x0);
    Fp_init(&A->x1);
}
void Fp2_set(struct Fp2 *ANS,struct Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set(&ANS->x1,&A->x1);
}
void Fp2_set_ui(struct Fp2 *A,signed long int B){
    Fp_set_ui(&A->x0,B);
    Fp_set_ui(&A->x1,B);
}
void Fp2_random(struct Fp2 *A){
    Fp_random(&A->x0);
    Fp_random(&A->x1);
}
void Fp2_clear(struct Fp2 *A){
    Fp_clear(&A->x0);
    Fp_clear(&A->x1);
}
void Fp2_printf(struct Fp2 *A){
    gmp_printf("%Zd,%Zd\n",A->x0.x0,A->x1.x0);
}
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
    Fp_add_ui(&ANS->x0,&A->x0,B);
    Fp_add_ui(&ANS->x1,&A->x1,B);
}
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
}
void Fp2_sqr(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp tmp1,tmp2,tmp3;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    
    Fp_add(&tmp1,&A->x0,&A->x1);
    Fp_mul_basis(&tmp2,&A->x1);
    Fp_add(&tmp2,&tmp2,&A->x0);
    Fp_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp_mul(&ANS->x0,&tmp1,&tmp2);
    Fp_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp_mul_basis(&tmp1,&tmp3);
    Fp_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
}

void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp tmp1,tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    
    //set
    Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp_mul_basis(&ANS->x0,&tmp2);//b*d*v
    Fp_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
}
void Fp2_mul_i(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp,&A->x0);
    
    Fp_mul_basis(&ANS->x0,&A->x1);
    Fp_set(&ANS->x1,&tmp);
    
    Fp_clear(&tmp);
}
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int B){
    Fp_mul_ui(&ANS->x0,&A->x0,B);
    Fp_mul_ui(&ANS->x1,&A->x1,B);
}
void Fp2_mul_Fp(struct Fp2 *ANS,struct Fp2 *A,struct Fp *B){
    Fp_mul(&ANS->x0,&ANS->x0,B);
    Fp_mul(&ANS->x1,&ANS->x1,B);
}
void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t B){
    Fp_mul_mpz(&ANS->x0,&A->x0,B);
    Fp_mul_mpz(&ANS->x1,&A->x1,B);
}
void Fp2_invert(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 frob,buf;
    Fp2_init(&frob);
    Fp2_init(&buf);
    
    //Fp2_inv_map(&frob,A);
    Fp_set(&frob.x0,&A->x0);
    Fp_neg(&frob.x1,&A->x1);
    Fp2_mul(&buf,A,&frob);
    mpz_invert(buf.x0.x0,buf.x0.x0,PRIME_P);
    Fp_inv_num++;
    Fp2_mul_mpz(ANS,&frob,buf.x0.x0);
    
    Fp2_clear(&frob);
    Fp2_clear(&buf);
}
void Fp2_div(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp2 tmp;
    Fp2_init(&tmp);
    
    Fp2_invert(&tmp,B);
    Fp2_mul(ANS,A,&tmp);
    
    Fp2_clear(&tmp);
}
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t B){
    int i,length;
    length=(int)mpz_sizeinbase(B,2);
    char binary[length];
    mpz_get_str(binary,2,B);
    struct Fp2 buf;
    Fp2_init(&buf);
    Fp2_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        //Fp2_mul(&buf,&buf,&buf);
        Fp2_sqr(&buf,&buf);
        if(binary[i]=='1'){
            Fp2_mul(&buf,A,&buf);
        }
    }
    
    Fp2_set(ANS,&buf);
    Fp2_clear(&buf);
}
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 n,y,x,b,t,tmp_Fp2;
    Fp2_init(&n);
    Fp2_init(&y);
    Fp2_init(&x);
    Fp2_init(&b);
    Fp2_init(&t);
    Fp2_init(&tmp_Fp2);
    Fp2_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    
    while(Fp2_legendre(&n)!=-1){
        Fp2_random(&n);
    }
    
    mpz_pow_ui(q,PRIME_P,2);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    
    Fp2_pow(&y,&n,q);
    mpz_set(r,e);
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    Fp2_pow(&x,A,tmp_mpz);
    Fp2_pow(&tmp_Fp2,&x,set_2);
    Fp2_mul(&b,&tmp_Fp2,A);
    Fp2_mul(&x,&x,A);
    
    int m;
    
    while(Fp2_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp2_set(&tmp_Fp2,&b);
        while(Fp2_cmp_mpz(&tmp_Fp2,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp2_pow(&tmp_Fp2,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        Fp2_pow(&t,&y,tmp_mpz);
        Fp2_pow(&y,&t,set_2);
        mpz_set_ui(r,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&b,&b,&y);
    }
    
    Fp2_set(ANS,&x);
    
    Fp2_clear(&n);
    Fp2_clear(&y);
    Fp2_clear(&x);
    Fp2_clear(&b);
    Fp2_clear(&t);
    Fp2_clear(&tmp_Fp2);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t B){
    if(Fp_cmp_mpz(&A->x0,B)==0 && Fp_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}
int Fp2_legendre(struct Fp2 *a){
    mpz_t i;
    struct Fp2 tmp;
    Fp2_init(&tmp);
    mpz_init(i);
    
    mpz_pow_ui(i,PRIME_P,2);
    mpz_sub_ui(i,i,1);
    mpz_div_ui(i,i,2);
    
    Fp2_pow(&tmp,a,i);
    
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    
    if((Fp2_cmp_mpz(&tmp,cmp))==0){
        Fp2_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp2_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
void Fp2_neg(struct Fp2 *ans,struct Fp2 *a){
    Fp_neg(&ans->x0,&a->x0);
    Fp_neg(&ans->x1,&a->x1);
}

void Fp2_frobenius_map(struct Fp2 *ANS, struct Fp2 *A){
    struct Fp2 t_ans;
    Fp2_init(&t_ans);
    
    Fp_set(&t_ans.x0,&A->x0);
    if (mpz_cmp_ui(A->x1.x0,0)==0) {
        Fp_set(&t_ans.x1,&A->x1);
    }
    else{
        mpz_sub(t_ans.x1.x0,PRIME_P,A->x1.x0);
    }
    
    
    Fp2_set(ANS,&t_ans);
    
    Fp2_clear(&t_ans);
}

#pragma mark Fp4 methods

void Fp4_init(struct Fp4 *A){
    Fp2_init(&A->x0);
    Fp2_init(&A->x1);
}
void Fp4_set(struct Fp4 *ANS,struct Fp4 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_set(&ANS->x1,&A->x1);
}
void Fp4_set_ui(struct Fp4 *A,signed long int B){
    Fp2_set_ui(&A->x0,B);
    Fp2_set_ui(&A->x1,B);
}
void Fp4_random(struct Fp4 *A){
    Fp2_random(&A->x0);
    Fp2_random(&A->x1);
}
void Fp4_clear(struct Fp4 *A){
    Fp2_clear(&A->x0);
    Fp2_clear(&A->x1);
}
void Fp4_printf(struct Fp4 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0,A->x0.x1.x0,A->x1.x0.x0,A->x1.x1.x0);
}
void Fp4_add(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    Fp2_add(&ANS->x0,&A->x0,&B->x0);
    Fp2_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp4_add_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B){
    Fp2_add_ui(&ANS->x0,&A->x0,B);
    Fp2_add_ui(&ANS->x1,&A->x1,B);
}
void Fp4_sub(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    Fp2_sub(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp4_sqr(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp2 tmp1,tmp2,tmp3;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&tmp3);
    
    Fp2_add(&tmp1,&A->x0,&A->x1);
    Fp2_mul_i(&tmp2,&A->x1);
    Fp2_add(&tmp2,&tmp2,&A->x0);
    Fp2_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp2_mul(&ANS->x0,&tmp1,&tmp2);
    Fp2_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp2_mul_i(&tmp1,&tmp3);
    Fp2_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp2_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&tmp3);
}

void Fp4_mul(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp2 tmp1,tmp2;
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    
    //set
    Fp2_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp2_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp2_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp2_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp2_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp2_mul_i(&ANS->x0,&tmp2);//b*d*v
    Fp2_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp2_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp2_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
}
void Fp4_mul_v(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set(&tmp,A);
    
    Fp_mul_basis(&ANS->x0.x0,&tmp.x1.x1);
    Fp_set(&ANS->x0.x1,&tmp.x1.x0);
    Fp_set(&ANS->x1.x0,&tmp.x0.x0);
    Fp_set(&ANS->x1.x1,&tmp.x0.x1);
    
    Fp4_clear(&tmp);
}
void Fp4_mul_ui(struct Fp4 *ANS,struct Fp4 *A,unsigned long int B){
    Fp2_mul_ui(&ANS->x0,&A->x0,B);
    Fp2_mul_ui(&ANS->x1,&A->x1,B);
}
void Fp4_mul_mpz(struct Fp4 *ANS,struct Fp4 *A,mpz_t B){
    Fp2_mul_mpz(&ANS->x0,&A->x0,B);
    Fp2_mul_mpz(&ANS->x1,&A->x1,B);
}
void Fp4_mul_Fp(struct Fp4 *ANS,struct Fp4 *A,struct Fp *B){
    Fp2_mul_Fp(&ANS->x0,&A->x0,B);
    Fp2_mul_Fp(&ANS->x1,&A->x1,B);
}
void Fp4_invert(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 frob,buf;
    Fp4_init(&frob);
    Fp4_init(&buf);
    
    Fp2_set(&frob.x0,&A->x0);
    Fp2_neg(&frob.x1,&A->x1);
    Fp4_mul(&buf,A,&frob);
    Fp2_invert(&buf.x0,&buf.x0);
    Fp2_mul(&ANS->x0,&frob.x0,&buf.x0);
    Fp2_mul(&ANS->x1,&frob.x1,&buf.x0);
    
    Fp4_clear(&frob);
    Fp4_clear(&buf);
}
void Fp4_div(struct Fp4 *ANS,struct Fp4 *A,struct Fp4 *B){
    struct Fp4 tmp,t_ans;
    Fp4_init(&tmp);
    Fp4_init(&t_ans);
    
    Fp4_invert(&tmp,B);
    Fp4_mul(&t_ans,A,&tmp);
    
    Fp4_set(ANS,&t_ans);
    
    Fp4_clear(&tmp);
    Fp4_clear(&t_ans);
}
void Fp4_pow(struct Fp4 *ANS,struct Fp4 *A,mpz_t B){
    int i,length;
    length=(int)mpz_sizeinbase(B,2);
    char binary[length];
    mpz_get_str(binary,2,B);
    struct Fp4 buf;
    Fp4_init(&buf);
    Fp4_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        //Fp2_mul(&buf,&buf,&buf);
        Fp4_sqr(&buf,&buf);
        if(binary[i]=='1'){
            Fp4_mul(&buf,A,&buf);
        }
    }
    
    Fp4_set(ANS,&buf);
    Fp4_clear(&buf);
}
void Fp4_sqrt(struct Fp4 *ANS,struct Fp4 *A){
    struct Fp4 n,y,x,b,t,tmp_Fp4;
    Fp4_init(&n);
    Fp4_init(&y);
    Fp4_init(&x);
    Fp4_init(&b);
    Fp4_init(&t);
    Fp4_init(&tmp_Fp4);
    Fp4_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp4_legendre(&n)!=-1){
        Fp4_random(&n);
    }
    mpz_pow_ui(q,PRIME_P,12);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp4_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp4_pow(&x,A,tmp_mpz);
    Fp4_pow(&tmp_Fp4,&x,set_2);
    Fp4_mul(&b,&tmp_Fp4,A);
    Fp4_mul(&x,&x,A);
    
    int m;
    
    while(Fp4_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp4_set(&tmp_Fp4,&b);
        while(Fp4_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp4_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp4_pow(&t,&y,tmp_mpz);
        Fp4_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp4_mul(&x,&x,&t);
        Fp4_mul(&b,&b,&y);
    }
    
    Fp4_set(ANS,&x);
    
    Fp4_clear(&n);
    Fp4_clear(&y);
    Fp4_clear(&x);
    Fp4_clear(&b);
    Fp4_clear(&t);
    Fp4_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp4_legendre(struct Fp4 *a){
    mpz_t i,cmp;
    struct Fp4 tmp;
    Fp4_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,PRIME_P,4);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp4_pow(&tmp,a,i);
    
    if((Fp4_cmp_mpz(&tmp,cmp))==0){
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp4_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp4_cmp(struct Fp4 *A,struct Fp4 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp4_cmp_mpz(struct Fp4 *A,mpz_t B){
    if(Fp2_cmp_mpz(&A->x0,B)==0 && Fp2_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}
void Fp4_neg(struct Fp4 *ans,struct Fp4 *a){
    Fp2_neg(&ans->x0,&a->x0);
    Fp2_neg(&ans->x1,&a->x1);
}

void Fp4_mul_betainv(struct Fp4 *ANS)
{
    struct Fp4 tmp;
    Fp4_init(&tmp);
    Fp4_set_ui(&tmp, 0);
    mpz_set(tmp.x1.x1.x0,a_x);
    mpz_mul(tmp.x1.x1.x0,tmp.x1.x1.x0,C1_INV.x0);
    
    Fp4_set(ANS, &tmp);
    Fp4_clear(&tmp);
}

void Fp4_frobenius_map(struct Fp4 *ANS, struct Fp4 *A){
    struct Fp4 t_ans;
    Fp4_init(&t_ans);
    
    Fp2_frobenius_map(&t_ans.x0,&A->x0);
    Fp2_frobenius_map(&t_ans.x1,&A->x1);
    Fp2_mul_Fp(&t_ans.x1,&t_ans.x1,&p_C4);
    
    Fp4_set(ANS,&t_ans);
    
    Fp4_clear(&t_ans);
}


// #pragma mark Fp8 methods
void Fp8_init(struct Fp8 *A){
    Fp4_init(&A->x0);
    Fp4_init(&A->x1);
}
void Fp8_set(struct Fp8 *ANS,struct Fp8 *A){
    Fp4_set(&ANS->x0,&A->x0);
    Fp4_set(&ANS->x1,&A->x1);
}
void Fp8_set_ui(struct Fp8 *A,signed long int B){
    Fp4_set_ui(&A->x0,B);
    Fp4_set_ui(&A->x1,B);
}
void Fp8_random(struct Fp8 *A){
    Fp4_random(&A->x0);
    Fp4_random(&A->x1);
}
void Fp8_clear(struct Fp8 *A){
    Fp4_clear(&A->x0);
    Fp4_clear(&A->x1);
}
void Fp8_printf(struct Fp8 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x0.x0.x0.x0,A->x0.x0.x1.x0,A->x0.x1.x0.x0,A->x0.x1.x1.x0);
    gmp_printf("(%Zd,%Zd,%Zd,%Zd\n",A->x1.x0.x0.x0,A->x1.x0.x1.x0,A->x1.x1.x0.x0,A->x1.x1.x1.x0);
}
void Fp8_add(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    Fp4_add(&ANS->x0,&A->x0,&B->x0);
    Fp4_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp8_add_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    Fp4_add_ui(&ANS->x0,&A->x0,B);
    Fp4_add_ui(&ANS->x1,&A->x1,B);
}
void Fp8_sub(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    Fp4_sub(&ANS->x0,&A->x0,&B->x0);
    Fp4_sub(&ANS->x1,&A->x1,&B->x1);
}

void Fp8_sqr(struct Fp8 *ANS,struct Fp8 *A){
    
    struct Fp4 tmp1,tmp2,tmp3;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    
    Fp4_add(&tmp1,&A->x0,&A->x1);
    Fp4_mul_v(&tmp2,&A->x1);
    Fp4_add(&tmp2,&tmp2,&A->x0);
    Fp4_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp4_mul(&ANS->x0,&tmp1,&tmp2);
    Fp4_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp4_mul_v(&tmp1,&tmp3);
    Fp4_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp4_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    
}

void Fp8_mul(struct Fp8 *ANS,struct Fp8 *A, struct Fp8 *B){
    struct Fp4 tmp1,tmp2;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    
    //set
    Fp4_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp4_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp4_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp4_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp4_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp4_mul_v(&ANS->x0,&tmp2);//b*d*v
    Fp4_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp4_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp4_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
}
void Fp8_mul_v(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    Fp8_set(&tmp,A);
    
    Fp_mul_basis(&ANS->x0.x0.x0,&tmp.x1.x1.x1);
    Fp_set(&ANS->x0.x0.x1,&tmp.x1.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&tmp.x1.x0.x0);
    Fp_set(&ANS->x0.x1.x1,&tmp.x1.x0.x1);
    Fp_set(&ANS->x1.x0.x0,&tmp.x0.x0.x0);
    Fp_set(&ANS->x1.x0.x1,&tmp.x0.x0.x1);
    Fp_set(&ANS->x1.x1.x0,&tmp.x0.x1.x0);
    Fp_set(&ANS->x1.x1.x1,&tmp.x0.x1.x1);
    
    Fp8_clear(&tmp);
}
void Fp8_mul_ui(struct Fp8 *ANS,struct Fp8 *A,unsigned long int B){
    Fp4_mul_ui(&ANS->x0,&A->x0,B);
    Fp4_mul_ui(&ANS->x1,&A->x1,B);
}
void Fp8_mul_Fp(struct Fp8 *ANS,struct Fp8 *A,struct Fp *B){
    Fp4_mul_Fp(&ANS->x0,&A->x0,B);
    Fp4_mul_Fp(&ANS->x1,&A->x1,B);
}
void Fp8_mul_mpz(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    Fp4_mul_mpz(&ANS->x0,&A->x0,B);
    Fp4_mul_mpz(&ANS->x1,&A->x1,B);
}
void Fp8_invert(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 frob,buf;
    Fp8_init(&frob);
    Fp8_init(&buf);
    
    Fp4_set(&frob.x0,&A->x0);
    Fp4_neg(&frob.x1,&A->x1);
    Fp8_mul(&buf,A,&frob);
    Fp4_invert(&buf.x0,&buf.x0);
    Fp4_mul(&ANS->x0,&frob.x0,&buf.x0);
    Fp4_mul(&ANS->x1,&frob.x1,&buf.x0);
    
    Fp8_clear(&frob);
    Fp8_clear(&buf);
}
void Fp8_div(struct Fp8 *ANS,struct Fp8 *A,struct Fp8 *B){
    struct Fp8 tmp;
    Fp8_init(&tmp);
    
    Fp8_invert(&tmp,B);
    Fp8_mul(ANS,A,&tmp);
    
    Fp8_clear(&tmp);
}
void Fp8_pow(struct Fp8 *ANS,struct Fp8 *A,mpz_t B){
    int i,length;
    length=(int)mpz_sizeinbase(B,2);
    char binary[length];
    mpz_get_str(binary,2,B);
    struct Fp8 buf;
    Fp8_init(&buf);
    Fp8_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        //Fp2_mul(&buf,&buf,&buf);
        Fp8_sqr(&buf,&buf);
        if(binary[i]=='1'){
            Fp8_mul(&buf,A,&buf);
        }
    }
    
    Fp8_set(ANS,&buf);
    Fp8_clear(&buf);
}
void Fp8_sqrt(struct Fp8 *ANS,struct Fp8 *A){
    struct Fp8 n,y,x,b,t,tmp_Fp4;
    Fp8_init(&n);
    Fp8_init(&y);
    Fp8_init(&x);
    Fp8_init(&b);
    Fp8_init(&t);
    Fp8_init(&tmp_Fp4);
    Fp8_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp8_legendre(&n)!=-1){
        Fp8_random(&n);
    }
    mpz_pow_ui(q,PRIME_P,12);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp8_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp8_pow(&x,A,tmp_mpz);
    Fp8_pow(&tmp_Fp4,&x,set_2);
    Fp8_mul(&b,&tmp_Fp4,A);
    Fp8_mul(&x,&x,A);
    
    int m;
    
    while(Fp8_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp8_set(&tmp_Fp4,&b);
        while(Fp8_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp8_pow(&tmp_Fp4,&b,tmp_mpz);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp8_pow(&t,&y,tmp_mpz);
        Fp8_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp8_mul(&x,&x,&t);
        Fp8_mul(&b,&b,&y);
    }
    
    Fp8_set(ANS,&x);
    
    Fp8_clear(&n);
    Fp8_clear(&y);
    Fp8_clear(&x);
    Fp8_clear(&b);
    Fp8_clear(&t);
    Fp8_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp8_legendre(struct Fp8 *a){
    mpz_t i,cmp;
    struct Fp8 tmp;
    Fp8_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,PRIME_P,8);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp8_pow(&tmp,a,i);
    
    if((Fp8_cmp_mpz(&tmp,cmp))==0){
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp8_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp8_cmp(struct Fp8 *A,struct Fp8 *B){
    if(Fp4_cmp(&A->x0,&B->x0)==0 && Fp4_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp8_cmp_mpz(struct Fp8 *A,mpz_t B){
    if(Fp4_cmp_mpz(&A->x0,B)==0 && Fp4_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}
void Fp8_neg(struct Fp8 *ans,struct Fp8 *a){
    Fp4_neg(&ans->x0,&a->x0);
    Fp4_neg(&ans->x1,&a->x1);
}
void Fp8_frobenius_map(struct Fp8 *ANS, struct Fp8 *A){
    struct Fp8 tmp_ans;
    struct Fp4 ans_tmp4;
    Fp4_init(&ans_tmp4);
    Fp8_init(&tmp_ans);
    
    Fp4_frobenius_map(&tmp_ans.x0,&A->x0);
    Fp4_frobenius_map(&tmp_ans.x1,&A->x1);
    Fp4_mul_Fp(&tmp_ans.x1,&tmp_ans.x1,&p_C8);
    
    Fp2_mul_i(&tmp_ans.x1.x0, &tmp_ans.x1.x0);
    Fp2_mul_i(&tmp_ans.x1.x1, &tmp_ans.x1.x1);
    
    Fp8_set(ANS,&tmp_ans);
    
    //    Fp_clear(&pm5d8);
    Fp4_clear(&ans_tmp4);
    //    Fp_clear(&set_c1);
    Fp8_clear(&tmp_ans);
}

// #pragma mark Fp8 methods
void Fp16_init(struct Fp16 *A){
    Fp8_init(&A->x0);
    Fp8_init(&A->x1);
}
void Fp16_set(struct Fp16 *ANS,struct Fp16 *A){
    Fp8_set(&ANS->x0,&A->x0);
    Fp8_set(&ANS->x1,&A->x1);
}
void Fp16_set_ui(struct Fp16 *A,signed long int B){
    Fp8_set_ui(&A->x0,B);
    Fp8_set_ui(&A->x1,B);
}
void Fp16_random(struct Fp16 *A){
    Fp8_random(&A->x0);
    Fp8_random(&A->x1);
}
void Fp16_clear(struct Fp16 *A){
    Fp8_clear(&A->x0);
    Fp8_clear(&A->x1);
}
void Fp16_printf(struct Fp16 *A){
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x0.x0.x0.x0.x0,A->x0.x0.x0.x1.x0,A->x0.x0.x1.x0.x0,A->x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x0.x1.x0.x0.x0,A->x0.x1.x0.x1.x0,A->x0.x1.x1.x0.x0,A->x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x1.x0.x0.x0.x0,A->x1.x0.x0.x1.x0,A->x1.x0.x1.x0.x0,A->x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x1.x1.x0.x0.x0,A->x1.x1.x0.x1.x0,A->x1.x1.x1.x0.x0,A->x1.x1.x1.x1.x0);
}
void Fp16_add(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    Fp8_add(&ANS->x0,&A->x0,&B->x0);
    Fp8_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp16_add_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B){
    Fp8_add_ui(&ANS->x0,&A->x0,B);
    Fp8_add_ui(&ANS->x1,&A->x1,B);
}
void Fp16_sub(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    Fp8_sub(&ANS->x0,&A->x0,&B->x0);
    Fp8_sub(&ANS->x1,&A->x1,&B->x1);
}
void Fp16_sqr(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp8 tmp1,tmp2,tmp3;
    Fp8_init(&tmp1);
    Fp8_init(&tmp2);
    Fp8_init(&tmp3);
    
    Fp8_add(&tmp1,&A->x0,&A->x1);
    Fp8_mul_v(&tmp2,&A->x1);
    Fp8_add(&tmp2,&tmp2,&A->x0);
    Fp8_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp8_mul(&ANS->x0,&tmp1,&tmp2);
    Fp8_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp8_mul_v(&tmp1,&tmp3);
    Fp8_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp8_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp8_clear(&tmp1);
    Fp8_clear(&tmp2);
    Fp8_clear(&tmp3);
}
void Fp16_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp8 tmp1,tmp2;
    Fp8_init(&tmp1);
    Fp8_init(&tmp2);
    
    //set
    Fp8_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp8_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp8_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp8_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp8_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp8_mul_v(&ANS->x0,&tmp2);//b*d*v
    Fp8_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp8_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp8_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp8_clear(&tmp1);
    Fp8_clear(&tmp2);
}
void Fp16_mul_v(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_v(&tmp.x0,&A->x1);
    Fp8_set(&tmp.x1,&A->x0);
    
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void Fp16_mul_ui(struct Fp16 *ANS,struct Fp16 *A,unsigned long int B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp8_mul_ui(&tmp.x0,&A->x0,B);
    Fp8_mul_ui(&tmp.x1,&A->x1,B);
    
    Fp16_set(ANS,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_mul_Fp(struct Fp16 *ANS,struct Fp16 *A,struct Fp *B){
    Fp8_mul_Fp(&ANS->x0,&A->x0,B);
    Fp8_mul_Fp(&ANS->x1,&A->x1,B);
}
void Fp16_mul_mpz(struct Fp16 *ANS,struct Fp16 *A,mpz_t B){
    Fp8_mul_mpz(&ANS->x0,&A->x0,B);
    Fp8_mul_mpz(&ANS->x1,&A->x1,B);
}
void Fp16_invert(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 frob,buf;
    Fp16_init(&frob);
    Fp16_init(&buf);
    
    Fp8_set(&frob.x0,&A->x0);
    Fp8_neg(&frob.x1,&A->x1);
    Fp16_mul(&buf,A,&frob);
    Fp8_invert(&buf.x0,&buf.x0);
    Fp8_mul(&ANS->x0,&frob.x0,&buf.x0);
    Fp8_mul(&ANS->x1,&frob.x1,&buf.x0);
    
    Fp16_clear(&frob);
    Fp16_clear(&buf);
}
void Fp16_div(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp16_invert(&tmp,B);
    Fp16_mul(ANS,A,&tmp);
    
    Fp16_clear(&tmp);
}
void Fp16_pow(struct Fp16 *ANS,struct Fp16 *A,mpz_t B){
    int i,length;
    length=(int)mpz_sizeinbase(B,2);
    char binary[length];
    mpz_get_str(binary,2,B);
    struct Fp16 buf;
    Fp16_init(&buf);
    Fp16_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        //Fp2_mul(&buf,&buf,&buf);
        Fp16_sqr(&buf,&buf);
        if(binary[i]=='1'){
            Fp16_mul(&buf,A,&buf);
        }
    }
    
    Fp16_set(ANS,&buf);
    Fp16_clear(&buf);
}
void Fp16_sqrt(struct Fp16 *ANS,struct Fp16 *A){
    struct Fp16 n,y,x,b,t,tmp_Fp4;
    Fp16_init(&n);
    Fp16_init(&y);
    Fp16_init(&x);
    Fp16_init(&b);
    Fp16_init(&t);
    Fp16_init(&tmp_Fp4);
    Fp16_set(&n,A);
    
    mpz_t tmp_mpz,q,e,r,set_1,set_2;
    mpz_init(tmp_mpz);
    mpz_init(q);
    mpz_init(e);
    mpz_init(r);
    mpz_init(set_1);
    mpz_init(set_2);
    mpz_set_ui(set_1,1);
    mpz_set_ui(set_2,2);
    
    while(Fp16_legendre(&n)!=-1){
        Fp16_random(&n);
    }
    mpz_pow_ui(q,PRIME_P,16);
    mpz_sub_ui(q,q,1);
    mpz_set_ui(e,0);
    while(mpz_odd_p(q)==0){
        mpz_add_ui(e,e,1);
        mpz_div_ui(q,q,2);
    }
    Fp16_pow(&y,&n,q);
    
    mpz_set(r,e);
    
    mpz_sub_ui(tmp_mpz,q,1);
    mpz_div_ui(tmp_mpz,tmp_mpz,2);
    
    Fp16_pow(&x,A,tmp_mpz);
    Fp16_pow(&tmp_Fp4,&x,set_2);
    Fp16_mul(&b,&tmp_Fp4,A);
    Fp16_mul(&x,&x,A);
    
    int m;
    
    while(Fp16_cmp_mpz(&b,set_1)==1){
        m=-1;
        Fp16_set(&tmp_Fp4,&b);
        while(Fp16_cmp_mpz(&tmp_Fp4,set_1)==1){
            m++;
            mpz_pow_ui(tmp_mpz,set_2,m);
            Fp16_pow(&tmp_Fp4,&b,tmp_mpz);
            // Fp16_printf(&tmp_Fp4);
        }
        mpz_sub_ui(tmp_mpz,r,m);
        mpz_sub_ui(tmp_mpz,tmp_mpz,1);
        mpz_powm(tmp_mpz,set_2,tmp_mpz,PRIME_P);
        // gmp_printf("%Zd,%Zd,%d\n",tmp_mpz,r,m);
        Fp16_pow(&t,&y,tmp_mpz);
        Fp16_pow(&y,&t,set_2);
        // gmp_printf("%Zd,%Zd,\n",y.x0.x0.x0,y.x0.x1.x0);
        mpz_set_ui(r,m);
        Fp16_mul(&x,&x,&t);
        Fp16_mul(&b,&b,&y);
        
    }
    
    Fp16_set(ANS,&x);
    
    Fp16_clear(&n);
    Fp16_clear(&y);
    Fp16_clear(&x);
    Fp16_clear(&b);
    Fp16_clear(&t);
    Fp16_clear(&tmp_Fp4);
    mpz_clear(tmp_mpz);
    mpz_clear(q);
    mpz_clear(e);
    mpz_clear(r);
    mpz_clear(set_1);
}
int Fp16_legendre(struct Fp16 *a){
    mpz_t i,cmp;
    struct Fp16 tmp;
    Fp16_init(&tmp);
    mpz_init(i);
    mpz_init(cmp);
    mpz_set_ui(cmp,1);
    mpz_pow_ui(i,PRIME_P,16);
    mpz_sub_ui(i,i,1);
    mpz_tdiv_q_ui(i,i,2);
    Fp16_pow(&tmp,a,i);
    
    if((Fp16_cmp_mpz(&tmp,cmp))==0){
        Fp16_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return 1;
    }else{
        Fp16_clear(&tmp);
        mpz_clear(i);
        mpz_clear(cmp);
        return -1;
    }
}
int Fp16_cmp(struct Fp16 *A,struct Fp16 *B){
    if(Fp8_cmp(&A->x0,&B->x0)==0 && Fp8_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp16_cmp_mpz(struct Fp16 *A,mpz_t B){
    if(Fp8_cmp_mpz(&A->x0,B)==0 && Fp8_cmp_mpz(&A->x1,B)==0){
        return 0;
    }
    return 1;
}
void Fp16_neg(struct Fp16 *ans,struct Fp16 *a){
    Fp8_neg(&ans->x0,&a->x0);
    Fp8_neg(&ans->x1,&a->x1);
}

void Fp16_frobenius_map(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    if (mpz_cmp_ui(tmp_ans.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0, &p_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1, &p_M_C4);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x1, &p_M_C_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x0, &p_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x1, &p_M_C_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x0, &p_C4_C16);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x1.x0, &p_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x1.x1, &p_M_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x0.x1, &p_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x0.x0, &p_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x1.x1, &p_M_C_C_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x1.x0, &p_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x0.x0, &p_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x0.x1, &p_M_C_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}

void Fp16_frobenius_map_p2(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    Fp_neg(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0);
    Fp_neg(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x0, &p2_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x1, &p2_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x0, &p2_M_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x1, &p2_M_C8);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x0.x1, &p2_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x0.x0, &p2_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x1.x1, &p2_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x1.x0, &p2_M_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x0.x1, &p2_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x0.x0, &p2_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x1.x1, &p2_M_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x1.x0, &p2_M_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}


void Fp16_frobenius_map_p3(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0); //a0
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);//-a1
    if (mpz_cmp_ui(TMP.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p3_C4);//c4 a2
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p3_M_C4);//-c4 a3
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x1, &p3_M_C_C8); //-cc8 a5
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x0, &p3_C8);//c8 a4
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x1, &p3_M_C_C4_C8); // -c8c4c a7
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x0, &p3_C8_C4); //c8c4
    
    //from 9-16
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x1.x1, &p3_M_C_C4_C16);//-cc4c16
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x1.x0, &p3_C4_C16);//c4c16
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x0.x0, &p3_C16);//c16
    
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x0.x1, &p3_M_C16);//-c16
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x1.x0, &p3_C_C4_C8_C16);//cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x1.x1, &p3_M_C_C4_C8_C16);//-cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x0.x1, &p3_M_C_C8_C16);//-cc8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x0.x0, &p3_C8_C16);//c8c16
    
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}

void Fp16_frobenius_map_p4(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p4_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p4_C4);
    
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x0, &p4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x1, &p4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x0, &p4_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x1, &p4_C4_C8);
    
    //from 9-11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x0.x0, &p4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x0.x1, &p4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x1.x0, &p4_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x1.x1, &p4_C4_C16);
    //12-16
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x0.x0, &p4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x0.x1, &p4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x1.x0, &p4_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x1.x1, &p4_C4_C8_C16);
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}

void Fp16_frobenius_map_p8(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p8_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p8_C4);
    
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x0, &p8_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x1, &p8_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x0, &p8_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x1, &p8_C4_C8);
    
    //from 9-11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x0.x0, &p8_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x0.x1, &p8_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x1.x0, &p8_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x1.x1, &p8_C4_C16);
    //12-16
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x0.x0, &p8_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x0.x1, &p8_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x1.x0, &p8_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x1.x1, &p8_C4_C8_C16);
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}


void Fp16_frobenius_map_p5(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    if (mpz_cmp_ui(tmp_ans.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0, &p5_C4);
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1, &p5_M_C4);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x1, &p5_M_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x0, &p5_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x1, &p5_M_C_C4_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x0, &p5_C4_C16);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x1.x0, &p5_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x1.x1, &p5_M_C_C4_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x0.x1, &p5_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x0.x0, &p5_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x1.x1, &p5_M_C_C_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x1.x0, &p5_C4_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x0.x0, &p5_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x0.x1, &p5_M_C_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}


void Fp16_frobenius_map_p6(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp_ans;
    Fp16_init(&tmp_ans);
    Fp16_set(&tmp_ans, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //a0-a3
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &tmp_ans.x0.x0.x0.x0);
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &tmp_ans.x0.x0.x0.x1);
    Fp_neg(&ans_tmp8.x0.x0.x1.x0, &tmp_ans.x0.x0.x1.x0);
    Fp_neg(&ans_tmp8.x0.x0.x1.x1, &tmp_ans.x0.x0.x1.x1);
    
    //a4-a7
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &tmp_ans.x0.x1.x0.x0, &p6_C8);
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &tmp_ans.x0.x1.x0.x1, &p6_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &tmp_ans.x0.x1.x1.x0, &p6_M_C8);
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &tmp_ans.x0.x1.x1.x1, &p6_M_C8);
    
    //from a8-a11
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &tmp_ans.x1.x0.x0.x1, &p6_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &tmp_ans.x1.x0.x0.x0, &p6_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &tmp_ans.x1.x0.x1.x1, &p6_M_C_C16);
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &tmp_ans.x1.x0.x1.x0, &p6_M_C16);
    
    //a12-a15
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &tmp_ans.x1.x1.x0.x1, &p6_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &tmp_ans.x1.x1.x0.x0, &p6_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &tmp_ans.x1.x1.x1.x1, &p6_M_C_C8_C16);
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &tmp_ans.x1.x1.x1.x0, &p6_M_C8_C16);
    
    Fp8_set(&tmp_ans.x0, &ans_tmp8.x0);
    Fp8_set(&tmp_ans.x1, &ans_tmp8.x1);
    Fp16_set(ANS,&tmp_ans);
    
    Fp16_clear(&tmp_ans);
}

void Fp16_frobenius_map_p7(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 TMP;
    Fp16_init(&TMP);
    Fp16_set(&TMP, A);
    struct Fp16 ans_tmp8;
    Fp16_init(&ans_tmp8);
    
    //1-4
    Fp_set(&ans_tmp8.x0.x0.x0.x0, &TMP.x0.x0.x0.x0); //a0
    Fp_set(&ans_tmp8.x0.x0.x0.x1, &TMP.x0.x0.x0.x1);//-a1
    if (mpz_cmp_ui(TMP.x0.x0.x0.x1.x0, 0) != 0) {
        Fp_neg(&ans_tmp8.x0.x0.x0.x1, &ans_tmp8.x0.x0.x0.x1);
    }
    Fp_mul(&ans_tmp8.x0.x0.x1.x0, &TMP.x0.x0.x1.x0, &p7_C4);//c4 a2
    Fp_mul(&ans_tmp8.x0.x0.x1.x1, &TMP.x0.x0.x1.x1, &p7_M_C4);//-c4 a3
    //5-8
    Fp_mul(&ans_tmp8.x0.x1.x0.x0, &TMP.x0.x1.x0.x1, &p7_M_C_C8); //-cc8 a5
    Fp_mul(&ans_tmp8.x0.x1.x0.x1, &TMP.x0.x1.x0.x0, &p7_C8);//c8 a4
    Fp_mul(&ans_tmp8.x0.x1.x1.x0, &TMP.x0.x1.x1.x1, &p7_M_C_C4_C8); // -c8c4c a7
    Fp_mul(&ans_tmp8.x0.x1.x1.x1, &TMP.x0.x1.x1.x0, &p7_C8_C4); //c8c4
    //from 9-16
    Fp_mul(&ans_tmp8.x1.x0.x0.x0, &TMP.x1.x0.x1.x1, &p7_M_C_C4_C16);//-cc4c16
    Fp_mul(&ans_tmp8.x1.x0.x0.x1, &TMP.x1.x0.x1.x0, &p7_C4_C16);//c4c16
    Fp_mul(&ans_tmp8.x1.x0.x1.x0, &TMP.x1.x0.x0.x0, &p7_C16);//c16
    Fp_mul(&ans_tmp8.x1.x0.x1.x1, &TMP.x1.x0.x0.x1, &p7_M_C16);//-c16
    Fp_mul(&ans_tmp8.x1.x1.x0.x0, &TMP.x1.x1.x1.x0, &p7_C_C4_C8_C16);//cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x0.x1, &TMP.x1.x1.x1.x1, &p7_M_C_C4_C8_C16);//-cc4c8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x0, &TMP.x1.x1.x0.x1, &p7_M_C_C8_C16);//-cc8c16
    Fp_mul(&ans_tmp8.x1.x1.x1.x1, &TMP.x1.x1.x0.x0, &p7_C8_C16);//c8c16
    
    
    Fp8_set(&TMP.x0, &ans_tmp8.x0);
    Fp8_set(&TMP.x1, &ans_tmp8.x1);
    Fp16_set(ANS,&TMP);
    
    Fp16_clear(&TMP);
    Fp16_clear(&ans_tmp8);
}


//-----------------------------------------------------------------------------------------
// #pragma mark EFp2 methods
void EFp2_init(struct EFp2 *A){
    Fp2_init(&A->x);
    Fp2_init(&A->y);
    A->infity=FALSE;
}
void EFp2_set(struct EFp2 *A,struct EFp2 *B){
    Fp2_set(&A->x,&B->x);
    Fp2_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp2_set_infity(struct EFp2 *A){
    Fp2_set_ui(&A->x,0);
    Fp2_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp2_clear(struct EFp2 *A){
    Fp2_clear(&A->x);
    Fp2_clear(&A->y);
}
void EFp2_printf(struct EFp2 *A){
    gmp_printf("(%Zd,%Zd)(%Zd,%Zd)\n",A->x.x0.x0,A->x.x1.x0,A->y.x0.x0,A->y.x1.x0);
}
void EFp2_SCM_BIN(struct EFp2 *ANS,struct EFp2 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp2 Q,R;
    EFp2_init(&Q);
    EFp2_set(&Q,P);
    EFp2_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp2_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp2_ECA(&Q,&Q,P);
        }
    }
    EFp2_set(ANS,&Q);
    
    EFp2_clear(&Q);
    EFp2_clear(&R);
    return;
}
void EFp2_ECD(struct EFp2 *ANS, struct EFp2 *P){
    if(P->infity==TRUE){
        EFp2_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp2_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp2_set_infity(ANS);
        return;
    }
    
    struct Fp2 x,y,lambda,tmp;
    struct EFp2 t_ans;
    Fp2_init(&x);
    Fp2_init(&lambda);
    Fp2_init(&tmp);
    Fp2_init(&y);
    EFp2_init(&t_ans);
    
    Fp2_sqr(&x,&P->x);
    Fp2_add(&tmp,&x,&x);
    Fp2_add(&x,&tmp,&x);
    Fp_add_mpz(&x.x0,&x.x0,a_x);
    Fp2_add(&y,&P->y,&P->y);
    Fp2_div(&lambda,&x,&y);
    Fp2_sqr(&tmp,&lambda);
    Fp2_add(&x,&P->x,&P->x);
    Fp2_sub(&x,&tmp,&x);
    Fp2_sub(&tmp,&P->x,&x);
    Fp2_set(&t_ans.x,&x);
    Fp2_mul(&tmp,&tmp,&lambda);
    Fp2_sub(&t_ans.y,&tmp,&P->y);
    
    EFp2_set(ANS,&t_ans);
    
    Fp2_clear(&x);
    Fp2_clear(&lambda);
    Fp2_clear(&y);
    Fp2_clear(&tmp);
    EFp2_clear(&t_ans);
}
void EFp2_ECA(struct EFp2 *ANS, struct EFp2 *P1, struct EFp2 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp2_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp2_set(ANS,P2);
        return;
    }
    else if(Fp2_cmp(&P1->x,&P2->x)==0&&Fp2_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp2_set_infity(ANS);
        return;
    }
    else if(EFp2_cmp(P1,P2)==0){ // P=Q
        EFp2_ECD(ANS,P1);
        return;
    }
    
    struct Fp2 x,y,lambda,tmp;
    struct EFp2 t_ans;
    
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&lambda);
    Fp2_init(&tmp);
    EFp2_init(&t_ans);
    
    Fp2_sub(&x,&P2->x,&P1->x);
    Fp2_sub(&y,&P2->y,&P1->y);
    Fp2_div(&lambda,&y,&x);
    Fp2_sqr(&tmp,&lambda);
    Fp2_add(&x,&P1->x,&P2->x);
    Fp2_sub(&x,&tmp,&x);
    Fp2_sub(&tmp,&P1->x,&x);
    Fp2_set(&t_ans.x,&x);
    Fp2_mul(&tmp,&tmp,&lambda);
    Fp2_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp2_set(ANS,&t_ans);
    
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&lambda);
    Fp2_clear(&tmp);
    EFp2_clear(&t_ans);
}
int EFp2_cmp(struct EFp2 *A,struct EFp2 *B){
    if(Fp2_cmp(&A->x,&B->x)==0 && Fp2_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp2_random_set(struct EFp2 *ANS){
    struct EFp2 P;
    EFp2_init(&P);
    
    struct Fp2 x,a,tmp_fp;
    Fp2_init(&a);
    Fp2_init(&x);
    Fp2_init(&tmp_fp);
    
    mpz_t t2,p2,p22,tmp,r2;
    mpz_t set_3;
    mpz_init(set_3);
    mpz_set_ui(set_3,3);
    
    mpz_init(t2);
    mpz_init(p2);
    mpz_init(p22);
    mpz_init(tmp);
    mpz_init(r2);
    
    mpz_pow_ui(t2,trace_t,2);
    mpz_mul_ui(p2,PRIME_P,2);
    mpz_sub(tmp,t2,p2);
    
    mpz_pow_ui(p22,PRIME_P,2);
    mpz_add_ui(p22,p22,1);
    mpz_sub(p22,p22,tmp);
    
    do{
        Fp2_random(&x);
        //    Fp2_printf(&x);
        Fp2_pow(&a,&x,set_3);
        //    Fp2_printf(&a);
        Fp2_mul_mpz(&tmp_fp, &x, a_x);
        Fp2_add(&a, &a, &tmp_fp);
        //    mpz_add(a.x0.x0,a.x0.x0,tmp_fp.x0.x0);
        // Fp2_printf(&a);
        // printf("before crash %d\n",Fp2_legendre(&a));
        // printf("%d\n",Fp2_legendre(&a));
    }while(Fp2_legendre(&a)!=1);
    // Fp2_printf(&a);
    Fp2_sqrt(&P.y,&a);
    Fp2_set(&P.x,&x);
    
    mpz_t r12_div_r2;
    mpz_init(r12_div_r2);
    mpz_div(r12_div_r2,p22,order_r);
    mpz_div(r12_div_r2,r12_div_r2,order_r);
    
    printf("before scm bit leg====================\n");
    EFp2_SCM_BIN(ANS,&P,r12_div_r2);
    // EFp2_SCM_BIN(ANS,&P,p22);
    printf("after scm bit leg====================\n");
    EFp2_clear(&P);
    Fp2_clear(&a);
    Fp2_clear(&x);
    Fp2_clear(&tmp_fp);
    mpz_clear(tmp);
    mpz_clear(t2);
    mpz_clear(p22);
    mpz_clear(p2);
}

//-----------------------------------------------------------------------------------------
// #pragma mark EFp4 methods

void EFp4_init(struct EFp4 *A){
    Fp4_init(&A->x);
    Fp4_init(&A->y);
    A->infity=FALSE;
}
void EFp4_set(struct EFp4 *A,struct EFp4 *B){
    Fp4_set(&A->x,&B->x);
    Fp4_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp4_set_infity(struct EFp4 *A){
    Fp4_set_ui(&A->x,0);
    Fp4_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp4_set_EFp(struct EFp4 *ANS,struct EFp *A){
    Fp4_set_ui(&ANS->x,0);
    Fp4_set_ui(&ANS->y,0);
    
    Fp_set(&ANS->x.x0.x0,&A->x);
    Fp_set(&ANS->y.x0.x0,&A->y);
    ANS->infity=A->infity;
}
void EFp4_clear(struct EFp4 *A){
    Fp4_clear(&A->x);
    Fp4_clear(&A->y);
}
void EFp4_printf(struct EFp4 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd)\n",A->x.x0.x0.x0,A->x.x0.x1.x0,A->x.x1.x0.x0, A->x.x1.x1.x0);
    gmp_printf("(%Zd,%Zd,%Zd,%Zd)\n",A->y.x0.x0.x0,A->y.x0.x1.x0,A->y.x1.x0.x0, A->y.x1.x1.x0);
}
void EFp4_SCM_BIN(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp4_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp4_ECA(&Q,&Q,P);
        }
    }
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}
void EFp4_ECD(struct EFp4 *ANS, struct EFp4 *P){
    if(P->infity==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_infity(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_sqr(&x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    //    Fp4_mul_betainv(&a_betainv);
    //    Fp4_add(&x, &x, &a_betainv);
    
    Fp_add_mpz(&x.x0.x0,&x.x0.x0,a_x);
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_sqr(&tmp,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
void EFp4_ECD_Sparse(struct EFp4 *ANS, struct EFp4 *P){
    if(P->infity==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_infity(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_sqr(&x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    //Fp4_mul_betainv(&a_betainv);
    //Fp4_add(&x, &x, &a_betainv);
    
    // Fp_add_mpz(&x.x0.x0,&x.x0.x0,a_x);
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_sqr(&tmp,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
void EFp4_ECD_Pseudo_Sparse(struct EFp4 *ANS, struct EFp4 *P){
    if(P->infity==TRUE){
        EFp4_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp4_set_infity(ANS);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp,a_betainv;
    struct EFp4 t_ans;
    Fp4_init(&x);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    Fp4_init(&y);
    Fp4_init(&a_betainv);
    
    EFp4_init(&t_ans);
    
    
    Fp4_sqr(&x,&P->x);
    Fp4_add(&tmp,&x,&x);
    Fp4_add(&x,&tmp,&x);
    
    Fp4_mul_betainv(&a_betainv);
    Fp4_mul(&a_betainv, &a_betainv, &z_inv2);
    Fp4_add(&x, &x, &a_betainv);
    
    Fp4_add(&y,&P->y,&P->y);
    Fp4_div(&lambda,&x,&y);
    Fp4_sqr(&tmp,&lambda);
    Fp4_add(&x,&P->x,&P->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&lambda);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}

void EFp4_ECA(struct EFp4 *ANS, struct EFp4 *P1, struct EFp4 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp4_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp4_set(ANS,P2);
        return;
    }
    else if(Fp4_cmp(&P1->x,&P2->x)==0&&Fp4_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp4_set_infity(ANS);
        return;
    }
    else if(EFp4_cmp(P1,P2)==0){ // P=Q
        EFp4_ECD(ANS,P1);
        return;
    }
    
    struct Fp4 x,y,lambda,tmp;
    struct EFp4 t_ans;
    
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&lambda);
    Fp4_init(&tmp);
    EFp4_init(&t_ans);
    
    Fp4_sub(&x,&P2->x,&P1->x);
    Fp4_sub(&y,&P2->y,&P1->y);
    Fp4_div(&lambda,&y,&x);
    Fp4_sqr(&tmp,&lambda);
    Fp4_add(&x,&P1->x,&P2->x);
    Fp4_sub(&x,&tmp,&x);
    Fp4_sub(&tmp,&P1->x,&x);
    Fp4_set(&t_ans.x,&x);
    Fp4_mul(&tmp,&tmp,&lambda);
    Fp4_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp4_set(ANS,&t_ans);
    
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&lambda);
    Fp4_clear(&tmp);
    EFp4_clear(&t_ans);
}
int EFp4_cmp(struct EFp4 *A,struct EFp4 *B){
    if(Fp4_cmp(&A->x,&B->x)==0 && Fp4_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
//-----------------------------------------------------------------------------------------
// #pragma mark EFp8 methods

void EFp8_init(struct EFp8 *A){
    Fp8_init(&A->x);
    Fp8_init(&A->y);
    A->infity=FALSE;
}
void EFp8_set(struct EFp8 *A,struct EFp8 *B){
    Fp8_set(&A->x,&B->x);
    Fp8_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp8_set_infity(struct EFp8 *A){
    Fp8_set_ui(&A->x,0);
    Fp8_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp8_clear(struct EFp8 *A){
    Fp8_clear(&A->x);
    Fp8_clear(&A->y);
}
void EFp8_printf(struct EFp8 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd",A->x.x0.x0.x0.x0,A->x.x0.x0.x1.x0,A->x.x0.x1.x0.x0,A->x.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->x.x1.x0.x0.x0,A->x.x1.x0.x1.x0,A->x.x1.x1.x0.x0,A->x.x1.x1.x1.x0);
    
    gmp_printf("(%Zd,%Zd,%Zd,%Zd",A->y.x0.x0.x0.x0,A->y.x0.x0.x1.x0,A->y.x0.x1.x0.x0,A->y.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->y.x1.x0.x0.x0,A->y.x1.x0.x1.x0,A->y.x1.x1.x0.x0,A->y.x1.x1.x1.x0);
    
}
void EFp8_SCM_BIN(struct EFp8 *ANS,struct EFp8 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp8 Q,R;
    EFp8_init(&Q);
    EFp8_set(&Q,P);
    EFp8_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp8_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp8_ECA(&Q,&Q,P);
        }
    }
    EFp8_set(ANS,&Q);
    
    EFp8_clear(&Q);
    EFp8_clear(&R);
    return;
}
void EFp8_ECD(struct EFp8 *ANS, struct EFp8 *P){
    if(P->infity==TRUE){
        EFp8_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp8_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp8_set_infity(ANS);
        return;
    }
    
    struct Fp8 x,y,lambda,tmp;
    struct EFp8 t_ans;
    Fp8_init(&x);
    Fp8_init(&lambda);
    Fp8_init(&tmp);
    Fp8_init(&y);
    EFp8_init(&t_ans);
    
    Fp8_sqr(&x,&P->x);
    Fp8_add(&tmp,&x,&x);
    Fp8_add(&x,&tmp,&x);
    Fp8_add(&y,&P->y,&P->y);
    Fp8_div(&lambda,&x,&y);
    Fp8_sqr(&tmp,&lambda);
    Fp8_add(&x,&P->x,&P->x);
    Fp8_sub(&x,&tmp,&x);
    Fp8_sub(&tmp,&P->x,&x);
    Fp8_set(&t_ans.x,&x);
    Fp8_mul(&tmp,&tmp,&lambda);
    Fp8_sub(&t_ans.y,&tmp,&P->y);
    
    EFp8_set(ANS,&t_ans);
    
    Fp8_clear(&x);
    Fp8_clear(&lambda);
    Fp8_clear(&y);
    Fp8_clear(&tmp);
    EFp8_clear(&t_ans);
}
void EFp8_ECA(struct EFp8 *ANS, struct EFp8 *P1, struct EFp8 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp8_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp8_set(ANS,P2);
        return;
    }
    else if(Fp8_cmp(&P1->x,&P2->x)==0&&Fp8_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp8_set_infity(ANS);
        return;
    }
    else if(EFp8_cmp(P1,P2)==0){ // P=Q
        EFp8_ECD(ANS,P1);
        return;
    }
    
    struct Fp8 x,y,lambda,tmp;
    struct EFp8 t_ans;
    
    Fp8_init(&x);
    Fp8_init(&y);
    Fp8_init(&lambda);
    Fp8_init(&tmp);
    EFp8_init(&t_ans);
    
    Fp8_sub(&x,&P2->x,&P1->x);
    Fp8_sub(&y,&P2->y,&P1->y);
    Fp8_div(&lambda,&y,&x);
    Fp8_sqr(&tmp,&lambda);
    Fp8_add(&x,&P1->x,&P2->x);
    Fp8_sub(&x,&tmp,&x);
    Fp8_sub(&tmp,&P1->x,&x);
    Fp8_set(&t_ans.x,&x);
    Fp8_mul(&tmp,&tmp,&lambda);
    Fp8_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp8_set(ANS,&t_ans);
    
    Fp8_clear(&x);
    Fp8_clear(&y);
    Fp8_clear(&lambda);
    Fp8_clear(&tmp);
    EFp8_clear(&t_ans);
}
int EFp8_cmp(struct EFp8 *A,struct EFp8 *B){
    if(Fp8_cmp(&A->x,&B->x)==0 && Fp8_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
//-----------------------------------------------------------------------------------------
// #pragma mark EFp16 methods

void EFp16_init(struct EFp16 *A){
    Fp16_init(&A->x);
    Fp16_init(&A->y);
    A->infity=FALSE;
}
void EFp16_set(struct EFp16 *A,struct EFp16 *B){
    Fp16_set(&A->x,&B->x);
    Fp16_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp16_set_infity(struct EFp16 *A){
    Fp16_set_ui(&A->x,0);
    Fp16_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp16_set_EFp(struct EFp16 *A,struct EFp *B){
    Fp16_set_ui(&A->x,0);
    Fp16_set_ui(&A->y,0);
    
    Fp_set(&A->x.x0.x0.x0.x0,&B->x);
    Fp_set(&A->y.x0.x0.x0.x0,&B->y);
    A->infity=B->infity;
}
void EFp16_clear(struct EFp16 *A){
    Fp16_clear(&A->x);
    Fp16_clear(&A->y);
}
void EFp16_printf(struct EFp16 *A){
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,\n",A->x.x0.x0.x0.x0.x0,A->x.x0.x0.x0.x1.x0,A->x.x0.x0.x1.x0.x0,A->x.x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x.x0.x1.x0.x0.x0,A->x.x0.x1.x0.x1.x0,A->x.x0.x1.x1.x0.x0,A->x.x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->x.x1.x0.x0.x0.x0,A->x.x1.x0.x0.x1.x0,A->x.x1.x0.x1.x0.x0,A->x.x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->x.x1.x1.x0.x0.x0,A->x.x1.x1.x0.x1.x0,A->x.x1.x1.x1.x0.x0,A->x.x1.x1.x1.x1.x0);
    
    gmp_printf("(%Zd,%Zd,%Zd,%Zd,\n",A->y.x0.x0.x0.x0.x0,A->y.x0.x0.x0.x1.x0,A->y.x0.x0.x1.x0.x0,A->y.x0.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->y.x0.x1.x0.x0.x0,A->y.x0.x1.x0.x1.x0,A->y.x0.x1.x1.x0.x0,A->y.x0.x1.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd,\n",A->y.x1.x0.x0.x0.x0,A->y.x1.x0.x0.x1.x0,A->y.x1.x0.x1.x0.x0,A->y.x1.x0.x1.x1.x0);
    gmp_printf("%Zd,%Zd,%Zd,%Zd)\n",A->y.x1.x1.x0.x0.x0,A->y.x1.x1.x0.x1.x0,A->y.x1.x1.x1.x0.x0,A->y.x1.x1.x1.x1.x0);
}
void EFp16_SCM_BIN(struct EFp16 *ANS,struct EFp16 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp16 Q,R;
    EFp16_init(&Q);
    EFp16_set(&Q,P);
    EFp16_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp16_ECD(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp16_ECA(&Q,&Q,P);
        }
    }
    EFp16_set(ANS,&Q);
    
    EFp16_clear(&Q);
    EFp16_clear(&R);
    return;
}
void EFp16_ECD(struct EFp16 *ANS, struct EFp16 *P){
    if(P->infity==TRUE){
        EFp16_set(ANS,P);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp16_cmp_mpz(&P->y,cmp)==0){//P.y==0
        EFp16_set_infity(ANS);
        return;
    }
    
    struct Fp16 x,y,lambda,tmp;
    struct EFp16 t_ans;
    Fp16_init(&x);
    Fp16_init(&lambda);
    Fp16_init(&tmp);
    Fp16_init(&y);
    EFp16_init(&t_ans);
    
    Fp16_sqr(&x,&P->x);
    Fp16_add(&tmp,&x,&x);
    Fp16_add(&x,&tmp,&x);
    Fp_add_mpz(&x.x0.x0.x0.x0,&x.x0.x0.x0.x0,a_x);
    Fp16_add(&y,&P->y,&P->y);
    Fp16_div(&lambda,&x,&y);
    Fp16_sqr(&tmp,&lambda);
    Fp16_add(&x,&P->x,&P->x);
    Fp16_sub(&x,&tmp,&x);
    Fp16_sub(&tmp,&P->x,&x);
    Fp16_set(&t_ans.x,&x);
    Fp16_mul(&tmp,&tmp,&lambda);
    Fp16_sub(&t_ans.y,&tmp,&P->y);
    
    EFp16_set(ANS,&t_ans);
    
    Fp16_clear(&x);
    Fp16_clear(&lambda);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&t_ans);
}
void EFp16_ECA(struct EFp16 *ANS, struct EFp16 *P1, struct EFp16 *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp16_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp16_set(ANS,P2);
        return;
    }
    else if(Fp16_cmp(&P1->x,&P2->x)==0&&Fp16_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp16_set_infity(ANS);
        return;
    }
    else if(EFp16_cmp(P1,P2)==0){ // P=Q
        EFp16_ECD(ANS,P1);
        return;
    }
    
    struct Fp16 x,y,lambda,tmp;
    struct EFp16 t_ans;
    
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&lambda);
    Fp16_init(&tmp);
    EFp16_init(&t_ans);
    
    Fp16_sub(&x,&P2->x,&P1->x);
    Fp16_sub(&y,&P2->y,&P1->y);
    Fp16_div(&lambda,&y,&x);
    Fp16_sqr(&tmp,&lambda);
    Fp16_add(&x,&P1->x,&P2->x);
    Fp16_sub(&x,&tmp,&x);
    Fp16_sub(&tmp,&P1->x,&x);
    Fp16_set(&t_ans.x,&x);
    Fp16_mul(&tmp,&tmp,&lambda);
    Fp16_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp16_set(ANS,&t_ans);
    
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&lambda);
    Fp16_clear(&tmp);
    EFp16_clear(&t_ans);
}
int EFp16_cmp(struct EFp16 *A,struct EFp16 *B){
    if(Fp16_cmp(&A->x,&B->x)==0 && Fp16_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp16_random_set(struct EFp16 *ANS){
    struct EFp16 P,ans_temp;
    EFp16_init(&P);
    EFp16_init(&ans_temp);
    EFp16_init(&P);
    
    struct Fp16 x,a,tmp16;
    Fp16_init(&a);
    Fp16_init(&x);
    Fp16_init(&tmp16);
    
    //t16=a^16+b^16=((t^2-2p)^2-2p^2)^2-2p^4)^2-2p^8
    mpz_t t2,p2,p22,p4,p8,tmp1,tmp2,t16;
    mpz_init(t2);
    mpz_init(p2);
    mpz_init(p22);
    mpz_init(p4);
    mpz_init(p8);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(t16);
    
    mpz_pow_ui(tmp1,trace_t,2);//t^2
    mpz_mul_ui(p2,PRIME_P,2);//2p
    mpz_sub(t2,tmp1,p2); //t2=t^2-2p
    mpz_pow_ui(t2,t2,2);//t2=(t^2-2p)^2
    
    mpz_pow_ui(tmp1,PRIME_P,2); //p^2
    mpz_mul_ui(p22,tmp1,2);//2p^2
    mpz_sub(tmp1,t2,p22); // (t^2-2p)^2-2p^2
    mpz_pow_ui(tmp2,tmp1,2);//tmp2=((t^2-2p)^2-2p^2)^2
    
    mpz_pow_ui(tmp1,PRIME_P,4); //p^4
    mpz_mul_ui(p4,tmp1,2);//2p^4
    mpz_sub(tmp1,tmp2,p4); //(((t^2-2p)^2-2p^2)^2-2p^4)
    mpz_pow_ui(tmp2,tmp1,2);//tmp2=(((t^2-2p)^2-2p^2)^2-2p^4)^2
    
    mpz_pow_ui(tmp1,PRIME_P,8); //p^8
    mpz_mul_ui(p8,tmp1,2);//2p^8
    mpz_sub(t16,tmp2,p8);
    
    
    mpz_t r16_div_r2,sEFp_16;
    mpz_init(r16_div_r2);
    mpz_init(sEFp_16);
    mpz_pow_ui(tmp1,PRIME_P,16);
    mpz_add_ui(tmp1,tmp1,1);
    mpz_sub(sEFp_16,tmp1,t16);
    
    mpz_pow_ui(tmp1,order_r,2);
    
    //    printf("r^2 divisible %d\n",(int)mpz_divisible_p(sEFp_16,tmp1));
    mpz_tdiv_q(r16_div_r2,sEFp_16,order_r);
    mpz_tdiv_q(r16_div_r2,r16_div_r2,order_r);
    do{
        Fp16_random(&x);
        Fp16_sqr(&a,&x);
        Fp16_mul(&a,&a,&x);//x^3
        Fp16_mul_mpz(&tmp16,&x, a_x); //ax
        Fp16_add(&a, &a, &tmp16);//x^3+ax
    }while(Fp16_legendre(&a)!=1);
    Fp16_sqrt(&P.y,&a);
    Fp16_set(&P.x,&x);
    EFp16_SCM_BIN(ANS,&P,r16_div_r2);//R
    
    
    EFp16_SCM_BIN(&ans_temp, ANS, order_r);//T
    if (ans_temp.infity == TRUE)
    {
        printf("Check Successful \n");
    }
    
    EFp16_clear(&P);
    Fp16_clear(&a);
    Fp16_clear(&x);
    mpz_clear(t2);
    mpz_clear(p2);
    mpz_clear(p22);
    mpz_clear(p4);
    mpz_clear(p8);
    mpz_clear(tmp1);
    mpz_clear(tmp2);
    mpz_clear(t16);
}

void EFp16_random_set_G2(struct EFp16 *ANS){
    struct EFp16 P,P_frobenius,tmp_EFp16;
    EFp16_init(&P);
    EFp16_init(&P_frobenius);
    EFp16_init(&tmp_EFp16);
    
    EFp16_random_set(&P);
    
    EFp16_frobenius_map(&P_frobenius,&P);
    Fp16_neg(&tmp_EFp16.y,&P.y);
    Fp16_set(&tmp_EFp16.x,&P.x);
    
    EFp16_ECA(&tmp_EFp16,&tmp_EFp16,&P_frobenius);
    //    EFp16_printf(&tmp_EFp16);
    EFp16_set(ANS,&tmp_EFp16);
    
    EFp16_clear(&P);
    EFp16_clear(&P_frobenius);
    EFp16_clear(&tmp_EFp16);
}

//--------------------------
// #pragma mark EFp methods
void EFp_init(struct EFp *A){
    Fp_init(&A->x);
    Fp_init(&A->y);
    A->infity=FALSE;
}
void EFp_set(struct EFp *A,struct EFp *B){
    Fp_set(&A->x,&B->x);
    Fp_set(&A->y,&B->y);
    A->infity=B->infity;
}
void EFp_set_infity(struct EFp *A){
    Fp_set_ui(&A->x,0);
    Fp_set_ui(&A->y,0);
    A->infity=TRUE;
}
void EFp_clear(struct EFp *A){
    Fp_clear(&A->x);
    Fp_clear(&A->y);
}
void EFp_printf(struct EFp *A){
    gmp_printf("(%Zd,%Zd)\n",A->x.x0,A->y.x0);
}
void EFp_SCM_BIN(struct EFp *ANS, struct EFp *P,mpz_t j){
    int i;
    int r;//bit数
    r= (int)mpz_sizeinbase(j,2);
    
    struct EFp Q;
    EFp_init(&Q);
    EFp_set(&Q,P);
    
    for(i=r-2;i>=0;i--){
        if(mpz_tstbit(j,i)==1){
            EFp_ECD(&Q,&Q);
            EFp_ECA(&Q,&Q,P);
        }else{
            EFp_ECD(&Q,&Q);
        }
    }
    
    EFp_set(ANS,&Q);
    EFp_clear(&Q);
    return;
}
void EFp_ECD(struct EFp *ANS, struct EFp *P){
    if(P->infity==TRUE){
        EFp_set(ANS,P);
        return;
    }
    if(mpz_sgn(P->y.x0)==0){//P.y==0
        EFp_set_infity(ANS);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_ans;
    Fp_init(&x);
    Fp_init(&lambda);
    Fp_init(&tmp);
    Fp_init(&y);
    EFp_init(&t_ans);
    
    Fp_mul(&x,&P->x,&P->x);
    Fp_add(&tmp,&x,&x);
    Fp_add(&x,&tmp,&x);//3x^2+a
    //    gmp_printf("tmem A = %Zd\n",tmp_a);
    Fp_add_mpz(&x,&x,tmp_a);
    
    Fp_add(&y,&P->y,&P->y);//2y
    
    Fp_div(&lambda,&x,&y);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P->x,&P->x);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P->x,&x);
    
    
    Fp_set(&t_ans.x,&x);
    
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_ans.y,&tmp,&P->y);
    
    EFp_set(ANS,&t_ans);
    
    Fp_clear(&x);
    Fp_clear(&lambda);
    Fp_clear(&y);
    Fp_clear(&tmp);
    EFp_clear(&t_ans);
}
void EFp_ECA(struct EFp *ANS, struct EFp *P1, struct EFp *P2){
    if(P2->infity==TRUE){//if P2==inf
        EFp_set(ANS,P1);
        return;
    }
    else if(P1->infity==TRUE){//if P1==inf
        EFp_set(ANS,P2);
        return;
    }
    else if(Fp_cmp(&P1->x,&P2->x)==0 && Fp_cmp(&P1->y,&P2->y)==1){ //P1.x==P2.x&&P1.y!=P2.y
        EFp_set_infity(ANS);
        return;
    }
    else if(EFp_cmp(P1,P2)==0){ // P=Q
        EFp_ECD(ANS,P1);
        return;
    }
    
    struct Fp x,y,lambda,tmp;
    struct EFp t_ans;
    
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&lambda);
    Fp_init(&tmp);
    EFp_init(&t_ans);
    
    Fp_sub(&x,&P2->x,&P1->x);
    Fp_sub(&y,&P2->y,&P1->y);
    Fp_div(&lambda,&y,&x);
    Fp_mul(&tmp,&lambda,&lambda);
    Fp_add(&x,&P1->x,&P2->x);
    Fp_sub(&x,&tmp,&x);
    Fp_sub(&tmp,&P1->x,&x);
    Fp_set(&t_ans.x,&x);
    Fp_mul(&tmp,&tmp,&lambda);
    Fp_sub(&t_ans.y,&tmp,&P1->y);
    
    EFp_set(ANS,&t_ans);
    
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&lambda);
    Fp_clear(&tmp);
    EFp_clear(&t_ans);
}
int EFp_cmp(struct EFp *A,struct EFp *B){
    if(Fp_cmp(&A->x,&B->x)==0 && Fp_cmp(&A->y,&B->y)==0){
        return 0;
    }
    return 1;
}
void EFp_neg(struct EFp *ANS, struct EFp *A){
    struct EFp tmp;
    EFp_init(&tmp);
    Fp_neg(&tmp.y,&A->y);
    Fp_set(&tmp.x,&A->x);
    EFp_set(ANS,&tmp);
    EFp_clear(&tmp);
}
void EFp_random_set(struct EFp *ANS){
    struct Fp a,x,tmp;
    Fp_init(&a);
    Fp_init(&x);
    Fp_init(&tmp);
    
    struct EFp P,Q;
    EFp_init(&P);
    EFp_init(&Q);
    
    do{
        Fp_random(&x);
        Fp_mul(&a,&x,&x);
        Fp_mul(&a,&a,&x);
        Fp_mul_mpz(&tmp, &x, a_x);
        Fp_add(&a, &a, &tmp);
    }while(mpz_legendre(a.x0,PRIME_P)!=1);
    Q.infity=0;
    Fp_sqrt(&P.y,&a);
    Fp_set(&P.x,&x);
    EFp_set(ANS,&P);
    
    
    Fp_clear(&a);
    Fp_clear(&x);
    EFp_clear(&P);
}

void EFp16_to_EFp4_map(struct EFp4 *ANS,struct EFp16 *A){
    Fp4_set_ui(&ANS->x,0);
    Fp4_set_ui(&ANS->y,0);
    Fp4_set(&ANS->x,&A->x.x0.x1);
    Fp4_set(&ANS->y,&A->y.x1.x1);
    ANS->infity=A->infity;
}

void EFp4_to_EFp16_map(struct EFp16 *ANS,struct EFp4 *A){
    Fp16_set_ui(&ANS->x,0);
    Fp16_set_ui(&ANS->y,0);
    Fp4_set(&ANS->x.x0.x1,&A->x);
    Fp4_set(&ANS->y.x1.x1,&A->y);
    ANS->infity=A->infity;
}

//-----------------------------------------------------------------------------------------
void Fp16_frobenius_map_old(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp16_pow(&tmp,A,PRIME_P);
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
void EFp16_frobenius_map(struct EFp16 *ANS,struct EFp16 *A){
    struct EFp16 tmp;
    EFp16_init(&tmp);
    
    Fp16_frobenius_map(&tmp.x,&A->x);
    Fp16_frobenius_map(&tmp.y,&A->y);
    
    EFp16_set(ANS,&tmp);
    EFp16_clear(&tmp);
}


void Skew_Frobenius_map(struct EFp4 *ANS, struct EFp4 *Qt)
{
    struct EFp4 tmp_ans;
    EFp4_init(&tmp_ans);
    
    struct Fp4 Qt_x, Qt_y;
    Fp4_init(&Qt_x);
    Fp4_init(&Qt_y);
    Fp4_set(&tmp_ans.x, &Qt->x);
    Fp4_set(&tmp_ans.y, &Qt->y);
    
    
    Fp_mul(&Qt_x.x0.x0, &tmp_ans.x.x0.x1,&p_M_C_C8);
    Fp_mul(&Qt_x.x0.x1,&tmp_ans.x.x0.x0,&p_C8);
    Fp_mul(&Qt_x.x1.x0, &tmp_ans.x.x1.x1, &p_M_C_C4_C8);
    Fp_mul(&Qt_x.x1.x1, &tmp_ans.x.x1.x0, &p_C4_C16);
    
    Fp_mul(&Qt_y.x0.x0, &tmp_ans.y.x1.x1,&p_M_C_C_C4_C8_C16);
    Fp_mul(&Qt_y.x0.x1,&tmp_ans.y.x1.x0,&p_C4_C8_C16);
    Fp_mul(&Qt_y.x1.x0, &tmp_ans.y.x0.x0, &p_C_C8_C16);
    Fp_mul(&Qt_y.x1.x1, &tmp_ans.y.x0.x1, &p_M_C_C8_C16);
    
    
    Fp4_set(&tmp_ans.x, &Qt_x);
    Fp4_set(&tmp_ans.y, &Qt_y);
    
    EFp4_set(ANS,&tmp_ans);
    EFp4_clear(&tmp_ans);
    Fp4_clear(&Qt_x);
    Fp4_clear(&Qt_y);
}

void KSS_16_parameters(void){
    
    mpz_t tmp1,tmp2,two;
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(two);
    
    
    //set p,r
    mpz_t p_tmp,r_tmp,t_tmp;
    mpz_t xpow2,xpow4,xpow5,xpow6,xpow8,xpow9,xpow10;
    //    mpz_t tmp1,tmp2;
    
    mpz_init(p_tmp);
    mpz_init(r_tmp);
    mpz_init(t_tmp);
    mpz_init(xpow2);
    mpz_init(xpow4);
    mpz_init(xpow5);
    mpz_init(xpow6);
    mpz_init(xpow8);
    mpz_init(xpow9);
    mpz_init(xpow10);
    mpz_init(tmp1);
    mpz_init(tmp2);
    
    mpz_mul(xpow2,X,X);
    mpz_mul(xpow4,xpow2,xpow2);
    mpz_mul(xpow5,xpow4,X);
    mpz_mul(xpow6,xpow5,X);
    mpz_mul(xpow8,xpow6,xpow2);
    mpz_mul(xpow9,xpow8,X);
    mpz_mul(xpow10,xpow9,X);
    
    //t=1/35(2x^5+41x+35)
    mpz_mul_ui(tmp1,X,41);
    mpz_add_ui(tmp1,tmp1,35);
    mpz_mul_ui(tmp2,xpow5,2);
    mpz_add(t_tmp,tmp1,tmp2);
    
    mpz_div_ui(trace_t,t_tmp,35);
    
    //r=x^8+48x^4+625
    mpz_mul_ui(tmp1,xpow4,48);
    mpz_add_ui(r_tmp,xpow8,625);
    mpz_add(tmp2,tmp1,r_tmp);
    mpz_tdiv_q_ui(order_r,tmp2,61250);
    //     mpz_tdiv_q_ui(order_r,tmp2,49);
    //     mpz_set(order_r,tmp2);
    //    gmp_printf ("order = %Zd\n",order_r);
    // mpz_set(r,r_tmp);
    
    //p=1/980(x^10+2x^9+5x^8+48x^6+152x^5+240x^4+625x^2+2398x+3125)
    mpz_mul_ui(tmp1,xpow9,2);
    mpz_add(p_tmp,tmp1,xpow10);
    mpz_mul_ui(tmp1,xpow8,5);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow6,48);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow5,152);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow4,240);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,xpow2,625);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_mul_ui(tmp1,X,2398);
    mpz_add(p_tmp,tmp1,p_tmp);
    mpz_add_ui(p_tmp,p_tmp,3125);
    
    mpz_div_ui(PRIME_P,p_tmp,980);
    
    mpz_add_ui(order_EFp,PRIME_P,1);
    mpz_sub(order_EFp,order_EFp,trace_t);
    
    if(mpz_probab_prime_p(PRIME_P,25)==0){
        gmp_printf("p:%Zd\n",PRIME_P);
        printf("not  prime number!\n");
        exit(0);
    }
    
    
    struct EFp P,ANS;
    int legendre;
    struct Fp rhs,tmp_ax,x;
    mpz_init(tmp_a);
    Fp_init(&rhs);
    Fp_init(&tmp_ax);
    EFp_init(&P);
    EFp_init(&ANS);
    Fp_init(&x);
    mpz_init(tmp_a);
    mpz_set_ui(tmp_a,0);
    
    for(;;){
        mpz_add_ui(tmp_a,tmp_a,1);
        Fp_set_ui(&x,1);
        legendre=0;
        while(legendre !=1){
            mpz_powm_ui(rhs.x0,x.x0,3,PRIME_P);
            //            gmp_printf("tmp %Zd\n",tmp_a);
            mpz_mul(tmp_ax.x0,x.x0,tmp_a);
            Fp_add(&rhs, &rhs, &tmp_ax);
            if((legendre = mpz_legendre(rhs.x0,PRIME_P))==1){
                //gmp_printf("a in while = %Zd\n",rhs.x0);
                Fp_printf(&rhs);
                Fp_sqrt(&P.y,&rhs);
                Fp_set(&P.x,&x);
                EFp_SCM_BIN(&ANS,&P,order_EFp);
                //                printf("SCM  ==\n");
                //                EFp_printf(&ANS);
                if(ANS.infity == TRUE){
                    mpz_set(a_x,tmp_a);
                    // mpz_clear(tmp_a);
                    Fp_clear(&rhs);
                    Fp_clear(&x);
                    EFp_clear(&P);
                    EFp_clear(&ANS);
                    return;
                }
            }
            Fp_add_ui(&x,&x,1);
        }
    }
    return;
}


//-----------------------------------------------------------------------------------------
void Miller_algo(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    struct Fp16 l_sum,v_sum;
    Fp16_init(&l_sum);
    Fp16_init(&v_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    Fp_set_ui(&v_sum.x0.x0.x0.x0,1);
    
    struct EFp16 T;
    EFp16_init(&T);
    EFp16_set(&T,P);
    
    
    struct Fp16 ltt,ltp,v2t,vtp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    Fp16_init(&v2t);
    Fp16_init(&vtp);
    
    int i;
    struct Fp16 tmp1;
    Fp16_init(&tmp1);
    //    Fp16_init(&lambda);
    int r_bit;//bit数
    
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        Fp16_sqr(&l_sum,&l_sum);
        Fp16_sqr(&v_sum,&v_sum);
        
        ltt_q(&ltt,&T,Q);
        Fp16_mul(&l_sum,&l_sum,&ltt);
        
        EFp16_ECD(&T,&T);
        v2t_q(&v2t,&T,Q);
        Fp16_mul(&v_sum,&v_sum,&v2t);
        
        if(mpz_tstbit(loop,i)==1){
            ltp_q(&ltp,&T,P,Q);
            Fp16_mul(&l_sum,&l_sum,&ltp);
            
            EFp16_ECA(&T,&T,P);
            vtp_q(&vtp,&T,Q);
            Fp16_mul(&v_sum,&v_sum,&vtp);
        }
    }
    
    
    // EFp16_printf(&T);
    Fp16_div(ANS,&l_sum,&v_sum);
    
    Fp16_clear(&l_sum);
    Fp16_clear(&v_sum);
    EFp16_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&v2t);
    Fp16_clear(&vtp);
    Fp16_clear(&tmp1);
}
void Optimal_Miller(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    
    struct EFp16 T,EFp16_tmp;
    EFp16_init(&T);
    EFp16_init(&EFp16_tmp);
    
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct Fp16 Px_neg;
    Fp16_init(&Px_neg);
    Fp16_neg(&Px_neg,&P->x);//TODO Why neg Px?
    
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    struct EFp16 Q_neg;
    EFp16_init(&Q_neg);
    Fp16_neg(&Q_neg.y,&Q->y);
    Fp16_set(&Q_neg.x,&Q->x);
    
    if(X_bit_binary[x_bit]==-1){
        EFp16_set(&T,&Q_neg);
    }else{
        EFp16_set(&T,Q);
    }
    int i;
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_sqr(&l_sum,&l_sum);
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_sqr(&l_sum,&l_sum);
                
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_sqr(&l_sum,&l_sum);
                
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    //  EFp16_SCM_BIN(&EFp_tmp,Q,prime);
    struct EFp4 Q_bar,EFp4_tmp;
    EFp4_init(&Q_bar);
    EFp4_init(&EFp4_tmp);
    EFp16_to_EFp4_map(&Q_bar, Q);
    Skew_Frobenius_map(&EFp4_tmp,&Q_bar);
    EFp4_to_EFp16_map(&EFp16_tmp, &EFp4_tmp);
    //  EFp16_frobenius_map(&EFp_tmp, Q);
    
    
    ltp_q(&ltp,&T,&EFp16_tmp,P);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    ltt_q(&ltt,Q,P);
    
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp16_clear(&T);
    EFp16_clear(&EFp16_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&Q_bar);
    EFp4_clear(&EFp4_tmp);
}

void pre_calc_vector_final_exp(void)
{
    printf("\n\n VECTOR PRE_CALC \n");
    struct timeval t0;
    struct timeval t1;
    gettimeofday(&t0, 0);
    
    
    mpz_t B, A;
    mpz_init(B);
    mpz_init(A);
    
    mpz_add_ui(B,X,1);
    mpz_mul(B,B,B);
    mpz_add_ui(B,B,4);
    
    
    mpz_pow_ui(A,X,3);
    mpz_mul(A,A,B);
    mpz_add_ui(A,A,56);
    
    
    
    mpz_t x2,x3,x4;
    mpz_inits(x2,x3,x4,(mpz_ptr) NULL);
    mpz_mul(x2,X,X);
    mpz_mul(x3,x2,X);
    mpz_mul(x4,x2,x2);
    
    mpz_t m00, m11, m22, m33, m44, m55, m66, m77, tmp1;
    mpz_inits(m00, m11, m22, m33, m44, m55, m66, m77,tmp1, (mpz_ptr) NULL);
    
    //m00:=2*u^3*A+55*u^2*B;
    mpz_mul(m00,x3,A);
    mpz_mul_ui(m00,m00,2);
    mpz_mul_ui(tmp1,x2,55);
    mpz_mul(tmp1,tmp1,B);
    mpz_add(m00,m00,tmp1);
    
    //m11:=-4*u^2*A-75*u*B;
    mpz_mul(m11,x2,A);
    mpz_mul_si(m11,m11,-4);
    mpz_mul(tmp1,X,B);
    mpz_mul_ui(tmp1,tmp1,75);
    mpz_sub(m11,m11,tmp1);
    
    //m22:=-2*u*A-125*B;
    mpz_mul(m22,X,A);
    mpz_mul_si(m22,m22,-2);
    mpz_mul_ui(tmp1,B,125);
    mpz_sub(m22,m22,tmp1);
    
    //m33:=-u^4*A-24*u^3*B+196;
    mpz_mul(tmp1,x4,A);
    mpz_neg(m33,tmp1);
    mpz_mul(tmp1,x3,B);
    mpz_mul_si(tmp1,tmp1,-24);
    mpz_add(m33,m33,tmp1);
    mpz_add_ui(m33,m33,196);
    //     gmp_printf("m33 = %Zd\n",m33);
    
    //m44:=u^3*A+10*u^2*B;
    mpz_mul(m44,x3,A);
    mpz_mul(tmp1,x2,B);
    mpz_mul_ui(tmp1,tmp1,10);
    mpz_add(m44,m44,tmp1);
    
    //m55:=3*u^2*A+100*u*B;
    mpz_mul(m55,x2,A);
    mpz_mul_ui(m55,m55,3);
    mpz_mul(tmp1,X,B);
    mpz_mul_ui(tmp1,tmp1,100);
    mpz_add(m55,m55,tmp1);
    
    //m66:=-11*u*A-250*B;
    mpz_mul(m66,X,A);
    mpz_mul_si(m66,m66,-11);
    mpz_mul_si(tmp1,B,-250);
    mpz_add(m66,m66,tmp1);
    //    gmp_printf("m66 = %Zd\n",m66);
    //m77:=7*A;
    mpz_mul_ui(m77,A,7);
    
    gettimeofday(&t1, 0);
    double elapsed = timedifference_msec(t0, t1);
    printf("FINAL PREC CALC VECTOR ms: %f [ms]\n", elapsed);
}

void final_exp_hard(struct Fp16 *ANS,struct Fp16 *fd)
{
    mpz_t Xp1,temp_exp;
    mpz_init(Xp1);
    mpz_init(temp_exp);
    mpz_set(Xp1,X);
    mpz_add_ui(Xp1,Xp1,1);
    
    struct Fp16 t,t0, t1, t2, t3, t4, t5, t6,t7, t8, t9, t10, t11, t12, tmp16, tmp16_1;
    Fp16_init(&t0);
    Fp16_init(&t1);
    Fp16_init(&t2);
    Fp16_init(&t3);
    Fp16_init(&t4);
    Fp16_init(&t5);
    Fp16_init(&t6);
    Fp16_init(&t7);
    Fp16_init(&t8);
    Fp16_init(&t9);
    Fp16_init(&t10);
    Fp16_init(&t11);
    Fp16_init(&t12);
    Fp16_init(&tmp16);
    Fp16_init(&tmp16_1);
    Fp16_init(&t);
    
    
    struct Fp16 t00, t01, t02, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,t28, t29,t30, t31, t32,t33, t37, s0, s1, s2, s3;
    Fp16_init(&t00);
    Fp16_init(&t01);
    Fp16_init(&t02);
    Fp16_init(&t13);
    Fp16_init(&t14);
    Fp16_init(&t15);
    Fp16_init(&t16);
    Fp16_init(&t17);
    Fp16_init(&t18);
    Fp16_init(&t19);
    Fp16_init(&t20);
    Fp16_init(&t21);
    Fp16_init(&t22);
    Fp16_init(&t23);
    Fp16_init(&t24);
    Fp16_init(&t25);
    Fp16_init(&t26);
    Fp16_init(&t27);
    Fp16_init(&t28);
    Fp16_init(&t29);
    Fp16_init(&t30);
    Fp16_init(&t31);
    Fp16_init(&t32);
    Fp16_init(&t33);
    Fp16_init(&t37);
    Fp16_init(&s0);
    Fp16_init(&s1);
    Fp16_init(&s2);
    Fp16_init(&s3);
    
    
    Fp16_sqr(&t0, fd);
    Fp16_sqr(&t1, &t0);
    Fp16_pow(&t2, fd, Xp1);
    
    Fp16_pow(&t3, &t2, Xp1);
    Fp16_mul(&t4, &t3, &t1);
    
    Fp16_pow(&t5,&t4,X);
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t6, &t4,temp_exp);
    mpz_set_ui(temp_exp,8);
    Fp16_pow(&t7, &t1,temp_exp);
    mpz_set_ui(temp_exp,2);
    Fp16_pow(&t8, &t7, temp_exp);
    
    Fp8_set(&tmp16_1.x0, &t1.x0);
    Fp8_neg(&tmp16_1.x1, &t1.x1);
    
    Fp16_mul(&t9,&t7,&tmp16_1);
    Fp16_sqr(&t10, &t9);
    Fp16_pow(&t11, &t5, X);
    Fp16_pow(&t12, &t11, X);
    Fp16_mul(&t01, &t12, &t10);
    
    //    Fp16_pow(&tmp16, fd, A);
    //    printf("IF t01 eq M^A %d\n",Fp16_cmp(&t01, &tmp16));
    
    Fp16_pow(&t14, &t01, X);
    Fp16_sqr(&tmp16, &t14);
    
    //Fp16_invert(&t13, &tmp16);
    Fp8_set(&t13.x0, &tmp16.x0);
    Fp8_neg(&t13.x1, &tmp16.x1);
    
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t00, &t6, temp_exp);
    Fp16_pow(&t15, &t00, temp_exp);
    
    //Fp16_invert(&tmp16, &t15);
    Fp8_set(&tmp16.x0, &t15.x0);
    Fp8_neg(&tmp16.x1, &t15.x1);
    Fp16_mul(&t0, &t13, &tmp16);
    
    Fp16_sqr(&t16, &t0);
    mpz_set_ui(temp_exp,4);
    Fp16_pow(&t17, &t13, temp_exp);
    Fp16_mul(&t18, &t17, &t14);
    
    //t2:=t16*t18;
    Fp16_mul(&t2, &t16, &t18);
    
    //t19:=t14^(u);
    Fp16_pow(&t19, &t14, X);
    //t20:=t19^(u);
    Fp16_pow(&t20, &t19, X);
    //t21:=t20^(u);
    Fp16_pow(&t21, &t20, X);
    //t22:=t19^2;
    Fp16_sqr(&t22, &t19);
    
    //t23:=t5^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t23, &t5, temp_exp);
    //t24:=t23^5;
    Fp16_pow(&t24, &t23, temp_exp);
    //t25:=t24^3;
    mpz_set_ui(temp_exp,3);
    Fp16_pow(&t25, &t24, temp_exp);
    //t26:=t24*t25;
    Fp16_mul(&t26, &t24, &t25);
    //t27:=t22^2;
    Fp16_sqr(&t27, &t22);
    //t37:=(t27*t25)^(-1);
    Fp16_mul(&tmp16, &t27, &t25);
    
    // Fp16_invert(&t37, &tmp16);
    Fp8_set(&t37.x0, &tmp16.x0);
    Fp8_neg(&t37.x1, &tmp16.x1);
    
    //    mpz_set_str(temp_exp,"88805286319009225549297960754086960732617029553721085722102881418823661400",10);
    //    Fp16_pow(&tmp16, fd, temp_exp);
    //    Fp16_invert(&tmp16, &tmp16);
    //    Fp16_printf(&tmp16);
    //    printf("IF t01 eq M^A %d\n",Fp16_cmp(&t37, &tmp16));
    //    Fp16_printf(&t0);
    
    //t28:=t27*t19^(-1);
    Fp8_set(&tmp16.x0, &t19.x0);
    Fp8_neg(&tmp16.x1, &t19.x1);
    Fp16_mul(&t28, &t27, &tmp16);
    
    //t3:=t28*t26;
    Fp16_mul(&t3, &t28, &t26);
    
    
    //t29:=t11^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t29, &t11, temp_exp);
    //t30:=t29^2;
    Fp16_sqr(&t30, &t29);
    //t4:=t20*t30;
    Fp16_mul(&t4, &t20, &t30);
    
    //s0:=t20^2;
    Fp16_sqr(&s0, &t20);
    //s1:=t30^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&s1, &t30, temp_exp);
    //s2:=s1*t29;
    Fp16_mul(&s2, &s1, &t29);
    //s3:=s0*s2;
    Fp16_mul(&s3, &s0, &s2);
    
    //t31:=t12^(24);
    mpz_set_ui(temp_exp,24);
    Fp16_pow(&t31, &t12, temp_exp);
    
    //t5:=t21^(-1)*t31^(-1);
    Fp8_set(&tmp16_1.x0, &t31.x0);
    Fp8_neg(&tmp16_1.x1, &t31.x1);
    
    Fp8_set(&tmp16.x0, &t21.x0);
    Fp8_neg(&tmp16.x1, &t21.x1);
    
    Fp16_mul(&t5, &tmp16, &tmp16_1);
    
    //t6:=t8^3*t1;
    mpz_set_ui(temp_exp,3);
    Fp16_pow(&tmp16, &t8, temp_exp);
    Fp16_mul(&t6, &tmp16,&t1);
    
    //t7:=t5*t6;
    Fp16_mul(&t7, &t5, &t6);
    
    //t8:=t01^7;
    mpz_set_ui(temp_exp,7);
    Fp16_pow(&t8, &t01, temp_exp);
    
    //    Fp16_pow(&tmp16, fd, m77);
    //    printf("IF t8 eq M^m77 %d\n",Fp16_cmp(&t8, &tmp16));
    
    //t32:=t37^p*t7^(p^3);
    Fp16_frobenius_map(&tmp16, &t37);
    Fp16_frobenius_map_p3(&tmp16_1, &t7);
    Fp16_mul(&t32, &tmp16_1, &tmp16);
    
    //t32:=t32*t3^(p^5);
    Fp16_frobenius_map_p5(&tmp16, &t3);
    Fp16_mul(&t32, &t32, &tmp16);
    
    //t32:=t32*t8^(p^7);
    Fp16_frobenius_map_p7(&tmp16, &t8);
    Fp16_mul(&t32, &t32, &tmp16);
    
    //t33:=t0^(p^2)*t2^(p^6);
    Fp16_frobenius_map_p2(&tmp16, &t0);
    Fp16_frobenius_map_p6(&tmp16_1, &t2);
    Fp16_mul(&t33, &tmp16, &tmp16_1);
    
    // t7=M^m33 and t2=M^m66.
    
    //t:=t32*t33;
    Fp16_mul(&t, &t33, &t32);
    
    //t:=t*t4^(p^4);
    Fp16_frobenius_map_p4(&tmp16, &t4);
    Fp16_mul(&t, &t, &tmp16);
    
    //t:=t*s3;
    Fp16_mul(&t, &t, &s3);
    
    Fp16_set(ANS, &t);
    
    
    mpz_clears(m00, m11, m22, m33, m44, m55, m66, m77,tmp1, (mpz_ptr) NULL);
    mpz_clears(x2,x3,(mpz_ptr) NULL);
    
    Fp16_clear(&t0);
    Fp16_clear(&t1);
    Fp16_clear(&t2);
    Fp16_clear(&t3);
    Fp16_clear(&t4);
    Fp16_clear(&t5);
    Fp16_clear(&t6);
    Fp16_clear(&t7);
    Fp16_clear(&t8);
    Fp16_clear(&t9);
    Fp16_clear(&t10);
    Fp16_clear(&t11);
    Fp16_clear(&t12);
    Fp16_clear(&tmp16);
    Fp16_clear(&tmp16_1);
    Fp16_clear(&t);
    Fp16_clear(&t00);
    Fp16_clear(&t01);
    Fp16_clear(&t02);
    Fp16_clear(&t13);
    Fp16_clear(&t14);
    Fp16_clear(&t15);
    Fp16_clear(&t16);
    Fp16_clear(&t17);
    Fp16_clear(&t18);
    Fp16_clear(&t19);
    Fp16_clear(&t20);
    Fp16_clear(&t21);
    Fp16_clear(&t22);
    Fp16_clear(&t23);
    Fp16_clear(&t24);
    Fp16_clear(&t25);
    Fp16_clear(&t26);
    Fp16_clear(&t27);
    Fp16_clear(&t28);
    Fp16_clear(&t29);
    Fp16_clear(&t30);
    Fp16_clear(&t31);
    Fp16_clear(&t32);
    Fp16_clear(&t33);
    Fp16_clear(&t37);
    Fp16_clear(&s0);
    Fp16_clear(&s1);
    Fp16_clear(&s2);
    Fp16_clear(&s3);
}
void Final_Exp(struct Fp16 *ANS,struct Fp16 *A){
    
    struct Fp16 M,A_p8;
    Fp16_init(&M);
    Fp16_init(&A_p8);
    
    printf("\n\n f^p^8-1 \n");
    struct timeval t0;
    struct timeval t1;
    gettimeofday(&t0, 0);
    Fp16_frobenius_map_p8(&A_p8,A);
    Fp16_div(&M, &A_p8, A);
    gettimeofday(&t1, 0);
    double elapsed_0 = timedifference_msec(t0, t1);
    printf("FINAL EXP EASY PART: %f [ms]\n", elapsed_0);
    
    
    printf("\n\n M^P^8+1 \n");
    gettimeofday(&t0, 0);
    final_exp_hard(ANS, &M);
    gettimeofday(&t1, 0);
    double elapsed = timedifference_msec(t0, t1);
    printf("FINAL EXP HARD PART ms: %f [ms]\n", elapsed);
    
    Fp16_clear(&M);
    Fp16_clear(&A_p8);
}
void Tate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    Miller_algo(&t_ans,G2,G1,order_r);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
}
void Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,trace_t,1);
    
    Miller_algo(&t_ans,G2,G1,tm1);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
}
void Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 Miller_X, t_ans;
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    Optimal_Miller(&Miller_X,G1,G2,X);
    
    Final_Exp(&t_ans,&Miller_X);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}
void ADD_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q,struct Fp16 *Qx_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_sub(&A,&P->x,&T->x);//xt-xp
    Fp16_sub(&B,&P->y,&T->y);//yt-yp
    Fp16_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_add(&D,&T->x,&P->x);
    Fp16_mul(&tmp1,&C,&C);
    Fp16_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp16_mul(&tmp2,&C,&T->x);
    Fp16_sub(&E,&tmp2,&T->y);
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x);
    Fp16_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp16_set(&l_tmp,&Q->y);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Qx_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
}
//TODO *Q = *P, *Qx_neg = * Px_neg
void DBL_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct Fp16 *Px_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_add(&A,&T->y,&T->y);//2y
    Fp16_mul(&B,&T->x,&T->x);//x^2
    Fp16_mul_ui(&B,&B,3);//3x^2
    Fp_add_mpz(&B.x0.x0.x0.x0,&B.x0.x0.x0.x0,a_x);//lambda=3x^2+a
    Fp16_div(&C,&B,&A);//lambda=3x^2+a/2y
    
    Fp16_add(&D,&T->x,&T->x); //D=2x
    Fp16_mul(&tmp1,&C,&C);// lamda^2
    Fp16_sub(&x3_tmp.x,&tmp1,&D); //x3.x=lamda^2-x2t
    
    Fp16_mul(&tmp2,&C,&T->x);//xt.lamda
    Fp16_sub(&E,&tmp2,&T->y);//xt.lamda-yt
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x); //x3.lamda
    Fp16_sub(&x3_tmp.y,&E,&tmp3); // x3.y = xt.lamda-yt - x3.lamda
    
    Fp16_set(&l_tmp,&P->y);
    // Fp_set_ui(&l_tmp.x0.x0.x0,1);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Px_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    if(T->infity==TRUE){
        EFp16_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp16_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp16_set_infity(T_ANS);
        return;
    }
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
    mpz_clear(cmp);
}
//-------------------
void ltt_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,lambda,ltt;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&lambda);
    Fp16_init(&ltt);
    
    Fp16_mul(&tmp1,&T->x,&T->x);//xt^2
    Fp16_add(&tmp2,&tmp1,&tmp1);
    Fp16_add(&tmp1,&tmp2,&tmp1);//3xt^2
    Fp_add_mpz(&tmp1.x0.x0.x0.x0,&tmp1.x0.x0.x0.x0,a_x);//TODO
    Fp16_add(&tmp2,&T->y,&T->y);//2yt
    
    Fp16_div(&lambda,&tmp1,&tmp2);//lambda=3xt^2+a/2yt
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=xq-xt
    Fp16_mul(&tmp3,&tmp3,&lambda);//tmp3=lambda(xq-xt)
    
    Fp16_sub(&ltt,&Q->y,&T->y);//yq-yt
    Fp16_sub(&ltt,&ltt,&tmp3);//ltt=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltt);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&lambda);
    Fp16_clear(&ltt);
}
void v2t_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 v2t;
    Fp16_init(&v2t);
    
    Fp16_sub(&v2t,&Q->x,&T->x);//v2t=xq-xt
    Fp16_set(ANS,&v2t);
    
    Fp16_clear(&v2t);
}
void ltp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    if((Fp16_cmp(&T->x,&P->x))==0&&(Fp16_cmp(&T->y,&P->y))!=0){//xt==xp&&yt!=yp
        Fp16_sub(&ltp,&Q->x,&T->x);
        Fp16_set(ANS,&ltp);
        
        return;
    }
    
    Fp16_sub(&tmp1,&T->x,&P->x);//xt-xp
    Fp16_sub(&tmp2,&T->y,&P->y);//yt-yp
    Fp16_div(&lambda,&tmp2,&tmp1);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=(xq-xt)
    Fp16_mul(&tmp4,&tmp3,&lambda);//tmp4=lambda(xq-xt)
    
    Fp16_sub(&ltp,&Q->y,&T->y);//ltp=yq-yt
    Fp16_sub(&ltp,&ltp,&tmp4);//ltp=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
}
void vtp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 vtp;
    Fp16_init(&vtp);
    if(T->infity==1){//if T is infity
        Fp16_set_ui(ANS,0);
        Fp_set_ui(&ANS->x0.x0.x0.x0,1);
        return;
    }
    
    Fp16_sub(&vtp,&Q->x,&T->x);
    Fp16_set(ANS,&vtp);
    
    Fp16_clear(&vtp);
}

#pragma mark Pseudo Sparse
//P=Q_map, Q =P_map
void Pseudo_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp4 *Q,struct EFp *P,struct Fp *L){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
    
    Fp4_sub(&A,&Q->x,&T->x);//xt-xp
    Fp4_sub(&B,&Q->y,&T->y);//yt-yp
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&Q->x);
    Fp4_mul(&tmp1,&C,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul_mpz(&l_tmp.x1.x1,&E,L->x0);
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T,&x3_tmp);
    
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
}
void Pseudo_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp *P,struct Fp *L){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltt;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltt);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
    mpz_t ac_inv;
    mpz_init(ac_inv);
    
    Fp4_add(&A,&T->y,&T->y);//2yt
    Fp4_mul(&B,&T->x,&T->x);//xt^2
    Fp4_mul_ui(&B,&B,3);//3xt^2
    Fp4_mul_betainv(&tmp);
    Fp4_mul_mpz(&tmp, &tmp, z_inv2_test.x0);
    Fp4_add(&B, &B, &tmp);
    
    Fp4_div(&C,&B,&A);//
    
    Fp4_add(&D,&T->x,&T->x);//D=2xt
    Fp4_mul(&tmp1,&C,&C);//C^2
    Fp4_sub(&x3_tmp.x,&tmp1,&D);//x3.x = C^2-D
    
    Fp4_mul(&tmp2,&C,&T->x); // C*xt
    Fp4_sub(&E,&tmp2,&T->y); //E=C*xt-yt
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);//C*x3.x
    Fp4_sub(&x3_tmp.y,&E,&tmp3);//x3.y = E-C*x3.x
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul_mpz(&l_tmp.x1.x1,&E,L->x0);//l_tmp=
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    Fp16_set(l_ANS,&l_tmp);
    
    EFp4_set(T,&x3_tmp);
    
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltt);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
}
void Pseudo_type1_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    
    struct Fp16 a,b;
    Fp16_init(&a);
    Fp16_init(&b);
    
    struct Fp4 a0,a1,a2,a3,b2,b3,C_0,C_1,C_2,C_3,T0,T1,T2,T3,T4,a2a3,tmp2,tmp3;
    Fp4_init(&a0);
    Fp4_init(&a1);
    Fp4_init(&a2);
    Fp4_init(&a3);
    Fp4_init(&b2);
    Fp4_init(&b3);
    
    Fp4_init(&C_0);
    Fp4_init(&C_1);
    Fp4_init(&C_2);
    Fp4_init(&C_3);
    
    Fp4_init(&T0);
    Fp4_init(&T1);
    Fp4_init(&T2);
    Fp4_init(&T3);
    Fp4_init(&T4);
    
    Fp4_init(&a2a3);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    
    
    Fp16_set(&a, A);
    Fp16_set(&b, B);
    
    Fp4_set(&a0, &a.x0.x0);
    Fp4_set(&a1, &a.x0.x1);
    Fp4_set(&a2, &a.x1.x0);
    Fp4_set(&a3, &a.x1.x1);
    
    Fp4_set(&b2, &b.x1.x0);
    Fp4_set(&b3, &b.x1.x1);
    
    
    Fp4_add(&a2a3, &a2, &a3);   /**< (a2+a3) */
    Fp4_add(&T4, &b2, &b3);     /**< t4=(b2+b3) */
    
    Fp4_mul(&T1, &a2, &b2);     /**< t1=(a2*b2) */
    Fp4_mul(&T0, &a3, &b3);     /**< t0=(a3*b3) */
    //    Fp4_mul_v(&T0, &tmp2);/**< t0=(a3*b3)*beta */
    
    Fp4_mul(&tmp2, &a2a3, &T4);
    Fp4_sub(&tmp3, &tmp2, &T1);
    Fp4_sub(&C_0, &tmp3, &T0);
    Fp4_mul_v(&tmp3, &C_0);
    Fp4_add(&C_0, &tmp3, &a0);
    
    Fp4_mul_v(&tmp2, &T0);
    Fp4_add(&C_1, &T1, &tmp2);
    Fp4_add(&C_1, &C_1, &a1);
    
    Fp4_mul(&T3, &a0, &b2);
    Fp4_mul(&T2, &a1, &b3);
    Fp4_mul_v(&tmp2, &T2);
    Fp4_add(&C_2, &T3, &tmp2);
    Fp4_add(&C_2, &C_2, &a2);
    
    
    Fp4_add(&tmp2, &a0, &a1);
    Fp4_mul(&tmp3, &tmp2, &T4);
    Fp4_sub(&tmp2, &tmp3, &T3);
    Fp4_sub(&C_3, &tmp2, &T2);
    Fp4_add(&C_3, &C_3, &a3);
    
    Fp4_set(&a.x0.x0,&C_0);
    Fp4_set(&a.x0.x1,&C_1);
    Fp4_set(&a.x1.x0,&C_2);
    Fp4_set(&a.x1.x1,&C_3);
    
    Fp16_set(ANS, &a);
    
    
    Fp16_clear(&a);
    Fp16_clear(&b);
    Fp4_clear(&a0);
    Fp4_clear(&a1);
    Fp4_clear(&a2);
    Fp4_clear(&a3);
    Fp4_clear(&b2);
    Fp4_clear(&b3);
    Fp4_clear(&C_0);
    Fp4_clear(&C_1);
    Fp4_clear(&C_2);
    Fp4_clear(&C_3);
    Fp4_clear(&T0);
    Fp4_clear(&T1);
    Fp4_clear(&T2);
    Fp4_clear(&T3);
    Fp4_clear(&T4);
}
void Pseudo_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){//Q:G2,P:G1
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,Q_map,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&Q_map);
    EFp4_init(&EFp4_tmp);
    
    struct EFp P_map;
    EFp_init(&P_map);
    
    struct Fp L,xy,xy_2,y_inv,tmp,y_tmp;
    Fp_init(&L);
    Fp_init(&xy);
    Fp_init(&xy_2);
    Fp_init(&y_inv);
    Fp_init(&tmp);
    Fp_init(&y_tmp);
    Fp_init(&z_inv2_test);
    
    mpz_invert(y_inv.x0,P->y.x0.x0.x0,PRIME_P);//yp^-1
    Fp_inv_num++;
    Fp_mul(&xy,&P->x.x0.x0,&y_inv);//xp.yp^-1
    
    Fp_sqr(&xy_2,&xy);//xy_2 = xp^2.yp^-2
    Fp_mul_mpz(&P_map.x,&P->x.x0.x0,xy_2.x0);//P.x= xp^3.yp^-2
    Fp_set(&P_map.y,&P_map.x);
    
    Fp_mul(&y_tmp,&xy_2,&xy);// xp^2.yp^-2 * xp.yp^-1 = xp^3.yp^-3
    Fp4_mul_mpz(&Q_map.y,&Q->y,y_tmp.x0); //Q_map.y = yQ'.xp^3.yp^-3
    Fp4_mul_mpz(&Q_map.x,&Q->x,xy_2.x0); //Q_map.x = xQ'.xp^2.yp^-2
    
    mpz_invert(L.x0,P_map.y.x0,PRIME_P); // L =yp_bar^-1
    Fp_inv_num++;
    Fp_set(&z_inv2_test,&xy_2);
    Fp_sqr(&z_inv2_test,&z_inv2_test);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q_map.y);
    Fp4_set(&Q_neg.x,&Q_map.x);
    
    EFp4_set(&T,&Q_map);
    
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_sqr(&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                break;
            case 1:
                Fp16_sqr(&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&P_map,&L);
                Pseudo_type1_ADD_LINE(&ltp,&T,&Q_map,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_sqr(&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&P_map,&L);
                Pseudo_type1_ADD_LINE(&ltp,&T,&Q_neg,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    Skew_Frobenius_map(&EFp4_tmp, &Q_map);
    Pseudo_type1_ADD_LINE(&ltp,&T,&EFp4_tmp,&P_map,&L);
    Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    Fp16_frobenius_map_p3(&tmp_f, &l_sum);
    
    Pseudo_type1_DBL_LINE(&ltt,&Q_map,&P_map,&L);
    Pseudo_type1_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    EFp4_clear(&Q_neg);
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp_clear(&P_map);
    EFp4_clear(&Q_map);
    EFp4_clear(&EFp4_tmp);
    Fp_clear(&L);
    Fp_clear(&xy);
    Fp_clear(&xy_2);
    Fp_clear(&y_inv);
    Fp_clear(&tmp);
    Fp_clear(&y_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&tmp_f);
}

void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct EFp4 G1_EFp4,G2_EFp4;
    EFp4_init(&G1_EFp4);
    EFp4_init(&G2_EFp4);
    
    struct Fp16 ltp,Miller_X,t_ans;
    Fp16_init(&ltp);
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    struct timeval t0;
    struct timeval t1;
    
    EFp4_set_EFp(&G1_EFp4,G1);
    EFp16_to_EFp4_map(&G2_EFp4,G2);
    printf("\n\nPseudo Sparse Miller algo\n");
    
    mpz_mpz_mul=0;
    mpz_ui_mul=0;
    Fp_mpz_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    gettimeofday(&t0, 0);
    Pseudo_type1_Optimal_Miller(&Miller_X,&G1_EFp4,&G2_EFp4,X);
    gettimeofday(&t1, 0);
    double elapsed = timedifference_msec(t0,t1);
    printf("Miller Time : %f [ms]\n", elapsed);
    printf("Miller loop cost\nmpz_mpz_mul:%ld,mpz_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
    
    mpz_mpz_mul=0;
    mpz_ui_mul=0;
    Fp_mpz_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    Final_Exp(&t_ans,&Miller_X);
    printf("Final exp cost\nmpz_mpz_mul:%ld,mpz_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
    
    Fp16_set(ANS, &t_ans);
    
    Fp16_clear(&ltp);
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}

void Sparse_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    Fp4_sub(&A,&P->x,&T->x);//xt-xp
    Fp4_sub(&B,&P->y,&T->y);//yt-yp
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&P->x);
    Fp4_sqr(&tmp1,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
}

void Sparse_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    
    Fp4_add(&A,&T->y,&T->y);//xt-xp
    Fp4_sqr(&B,&T->x);
    Fp4_mul_ui(&B,&B,3);
    struct Fp4 ac_inv;
    Fp4_init(&ac_inv);
    Fp4_mul_betainv(&ac_inv);
    Fp4_add(&B,&B,&ac_inv);
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&T->x);
    Fp4_sqr(&tmp1,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    if(T->infity==TRUE){
        EFp4_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp4_set_infity(T_ANS);
        return;
    }
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
    mpz_clear(cmp);
}

void Sparse_type1_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,EFp_tmp;
    EFp4_init(&T);
    EFp4_init(&EFp_tmp);
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    Fp4_neg(&Px_neg,&P->x);
    
    EFp4_set(&T,Q);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    int r_bit;
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        if(mpz_tstbit(loop,i)==1){
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
            Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
            Fp16_mul(&l_sum,&l_sum,&ltp);
        }else{
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
        }
    }
    Fp16_set(ANS,&l_sum);
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp4_clear(&EFp_tmp);
    Fp4_clear(&Px_neg);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
}

void Sparse_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&EFp4_tmp);
    
    mpz_t p3;
    mpz_init(p3);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    
    Fp4_neg(&Px_neg,&P->x);
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q->y);
    Fp4_set(&Q_neg.x,&Q->x);
    
    if(X_bit_binary[x_bit]==-1){
        EFp4_set(&T,&Q_neg);
    }else{
        EFp4_set(&T,Q);
    }
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Sparse_type1_ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    Skew_Frobenius_map(&EFp4_tmp,Q);
    //    EFp4_SCM_BIN(&EFp4_tmp, &Q_map, prime);
    
    Sparse_type1_ADD_LINE(&ltp,&T,&T,&EFp4_tmp,P,&Px_neg);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    
    //ltt_q(&ltt,Q,P);
    Sparse_type1_DBL_LINE(&ltt,&T,Q,P,&Px_neg);
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&EFp4_tmp);
    
}
void EFp4_SCM_BIN_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp4_ECD_Sparse(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp4_ECA(&Q,&Q,P);
        }
    }
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}
void EFp4_SCM_BIN_Pseudo_Sparse(struct EFp4 *ANS,struct EFp4 *P,mpz_t j){
    int i,length;
    length= (int)mpz_sizeinbase(j,2);
    char j_binary[length];
    mpz_get_str(j_binary,2,j);
    struct EFp4 Q,R;
    EFp4_init(&Q);
    EFp4_set(&Q,P);
    EFp4_init(&R);
    for(i=1;j_binary[i]!='\0';i++){
        EFp4_ECD_Pseudo_Sparse(&Q,&Q);
        if(j_binary[i]=='1'){
            EFp4_ECA(&Q,&Q,P);
        }
    }
    EFp4_set(ANS,&Q);
    
    EFp4_clear(&Q);
    EFp4_clear(&R);
    return;
}



void Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp4 *G1,struct EFp4 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,trace_t,1);
    
    Sparse_type1_Miller(&t_ans,G1,G2,tm1);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
    mpz_clear(tm1);
}


void Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct EFp4 G1_EFp4,G2_EFp4;
    EFp4_init(&G1_EFp4);
    EFp4_init(&G2_EFp4);
    
    struct Fp16 ltp,Miller_X,t_ans;
    Fp16_init(&ltp);
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    EFp4_set_EFp(&G1_EFp4,G1);
    EFp16_to_EFp4_map(&G2_EFp4,G2);
    
    Sparse_type1_Optimal_Miller(&Miller_X,&G1_EFp4,&G2_EFp4,X);
    Final_Exp(&t_ans,&Miller_X);
    Fp16_set(ANS, &t_ans);
    
    Fp16_clear(&ltp);
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}

void check_Pairing(void){
    struct EFp P_EFp, R_EFp;
    EFp_init(&P_EFp);
    EFp_init(&R_EFp);
    
    
    struct EFp16 P_Fp16,Q_Fp16,Q_G2,R_Fp16,S_Fp16;
    EFp16_init(&P_Fp16);
    EFp16_init(&Q_Fp16);
    EFp16_init(&R_Fp16);
    EFp16_init(&S_Fp16);
    EFp16_init(&Q_G2);
    
    struct Fp16 ans_Fp16,tmp1_Fp16,tmp2_Fp16,tmp3_Fp16;
    Fp16_init(&ans_Fp16);
    Fp16_init(&tmp1_Fp16);
    Fp16_init(&tmp2_Fp16);
    Fp16_init(&tmp3_Fp16);
    
    struct EFp4 P_EFp4, Q_EFp4, R_EFp4, S_EFp4;
    EFp4_init(&P_EFp4);
    EFp4_init(&Q_EFp4);
    EFp4_init(&R_EFp4);
    EFp4_init(&S_EFp4);
    
    mpz_t a,b,ab;
    mpz_init(a);
    mpz_init(b);
    mpz_init(ab);
    
    mpz_set_ui(a,31);
    mpz_set_ui(b,11);
    mpz_mul(ab,a,b);
    
    mpz_t x,y;
    mpz_inits(x,y,(mpz_ptr) NULL);
    mpz_set_str(x,"585634432000126707887057629201458798521445673169852167601690963969803976546284829414253996667593644070", 10);
    mpz_set_str(y,"354142610898165644039571248952651466071208533406755516950532183093299777682028181848353083070437124285", 10);
    Fp_set_mpz(&P_EFp.x,x);
    Fp_set_mpz(&P_EFp.y,y);
    EFp16_set_EFp(&P_Fp16,&P_EFp);
    mpz_clears(x,y,(mpz_ptr) NULL);
    
    mpz_t x1, x2, x3, x4, y1, y2, y3, y4;
    mpz_inits(x1,x2,x3,x4,y1,y2,y3,y4,(mpz_ptr) NULL);
    
    mpz_set_str(x1,"244491741495785299367612042330502507448971522000836606635913541960426621971155507113603357863094662406", 10);
    mpz_set_str(x2,"218139782565182912336439776658495825407576298090161716969444130503315739115994708740875588793005720529", 10);
    mpz_set_str(x3,"489284606975189850057489668560782153703917092806842888347849447445917065112045634272024913067442383327", 10);
    mpz_set_str(x4,"246161136867400494929650830056723287620786052741538082500453356439354042813046547782023695596667138627", 10);
    mpz_set_str(y1,"48629828925217256502074899893821067838880538260652902679139733711894991914937403050096788961116072841", 10);
    mpz_set_str(y2,"338343177397929794197126001730091007944155731601968256694856657191699538746843500939799293603025576582", 10);
    mpz_set_str(y3,"38952231663050097715674495759233974069804273632029475897530255418880509349197207272848037344062136227", 10);
    mpz_set_str(y4,"250604271284793810746824678011335939079727656816914850518650242122895946469531880085329357599729742073", 10);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x0.x0,x1);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x0.x1,x2);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x1.x0,x3);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x1.x1,x4);
    
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x0.x0,y1);
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x0.x1,y2);
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x1.x0,y3);
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x1.x1,y4);
    
    printf("\n\nPseudo Sparse Optimal Ate Pairing \n");
    Big_M = 0; Small_m =0, Big_add=0, Sqr=0, fp16_pow =0;
    Pseudo_Sparse_Optimal_Ate_Pairing(&tmp1_Fp16,&P_EFp,&Q_Fp16);
    printf("Fp SQR in Fp16 mul = Fp16 pow %d, %d, %d, %d, %d\n",Big_M, Sqr, Small_m,Big_add,fp16_pow);
    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
 
    
    EFp_SCM_BIN(&R_EFp,&P_EFp,a);
    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
    Pseudo_Sparse_Optimal_Ate_Pairing(&tmp2_Fp16,&R_EFp,&S_Fp16);

    printf("Bilinearity Check: \n");
    if (Fp16_cmp(&tmp2_Fp16, &tmp1_Fp16) == 0) {
        printf("Success \n");
    }
    else{
        printf("Failure\n");
    }
    
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(ab);
    
    EFp_clear(&P_EFp);
    EFp_clear(&R_EFp);
    
    EFp16_clear(&P_Fp16);
    EFp16_clear(&Q_Fp16);
    EFp16_clear(&R_Fp16);
    EFp16_clear(&S_Fp16);
    
    Fp16_clear(&ans_Fp16);
    Fp16_clear(&tmp1_Fp16);
    Fp16_clear(&tmp2_Fp16);
    Fp16_clear(&tmp3_Fp16);
    
    EFp4_clear(&P_EFp4);
    EFp4_clear(&Q_EFp4);
    EFp4_clear(&R_EFp4);
    EFp4_clear(&S_EFp4);
}

float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}

