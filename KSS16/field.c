#include "field.h"
#include "params.h"
#include "util.h"

struct Fp4 beta,beta_inv;
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
struct Fp p5_C4, p5_M_C4;
struct Fp p5_C8, p5_M_C8;
struct Fp p5_C_C4_C8, p5_M_C_C4_C8;
struct Fp p5_C16, p5_M_C_C16;
struct Fp p5_C4_C16, p5_C_C4_C16, p5_M_C_C4_C16;
struct Fp p5_C4_C8_C16, p5_C_C4_C8_C16, p5_M_C_C_C4_C8_C16;
struct Fp p5_C8_C16, p5_C_C8_C16, p5_M_C_C8_C16;
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
struct Fp p4_C4;
struct Fp p4_C8;
struct Fp p4_C16;
struct Fp p4_C4_C8, p4_C8_C16;
struct Fp p4_C4_C16, p4_C4_C8_C16;
struct Fp p8_C4;
struct Fp p8_C8;
struct Fp p8_C16;
struct Fp p8_C4_C8, p8_C8_C16;
struct Fp p8_C4_C16, p8_C4_C8_C16;
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

void Fp16_frobenius_map_old(struct Fp16 *ANS, struct Fp16 *A){
    struct Fp16 tmp;
    Fp16_init(&tmp);
    
    Fp16_pow(&tmp,A,PRIME_P);
    Fp16_set(ANS,&tmp);
    Fp16_clear(&tmp);
}
