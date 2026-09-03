#include "field.h"
#include "params.h"
#include "test.h"

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
