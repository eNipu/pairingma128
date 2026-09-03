#include "params.h"

int X_bit_binary[x_bit+1];
int X_bit_binary_for_opt_ate[x_bit_for_opt_ate+1];
mpz_t prime;
mpz_t mother_parameter;
int sign;
mpz_t trace_t;
mpz_t EFp_order;
mpz_t EFp_total;
mpz_t EFp2_total;
mpz_t EFp6_total;
mpz_t EFp12_total;
mpz_t curve_parameter_A;
mpz_t curve_parameter_B;
mpz_t final_exp;
struct Fp Fp_basis;
struct Fp2 Fp2_basis;
struct Fp2 Fp2_basis_inv;
struct Fp6 Fp6_basis;
struct Fp Fp_ZERO;
struct Fp2 Fp2_ZERO;
struct Fp6 Fp6_ZERO;
struct Fp12 Fp12_ZERO;
struct Fp Fp_ONE;
struct Fp2 Fp2_ONE;
struct Fp6 Fp6_ONE;
struct Fp12 Fp12_ONE;
struct Fp inv_CNR1;
struct Fp inv_CNR2;
struct Fp epsilon_1;
struct Fp epsilon_2;
struct Fp2 Fp2_basis_prime_1_div_3_1;
struct Fp2 Fp2_basis_prime_1_div_3_2;
struct Fp2 Fp2_basis_prime_1_div_6;
struct Fp2 Fp2_basis_prime_2_div_3_1;
struct Fp2 Fp2_basis_prime_2_div_3_2;
struct Fp2 Fp2_basis_prime_2_div_6;
struct Fp2 Fp2_basis_prime_3_div_3_1;
struct Fp2 Fp2_basis_prime_3_div_3_2;
struct Fp2 Fp2_basis_prime_3_div_6;
struct Fp2 Fp2_basis_prime_4_div_3_1;
struct Fp2 Fp2_basis_prime_4_div_3_2;
struct Fp2 Fp2_basis_prime_4_div_6;
struct Fp2 Fp2_basis_prime_8_div_3_1;
struct Fp2 Fp2_basis_prime_8_div_3_2;
struct Fp2 Fp2_basis_prime_8_div_6;
struct Fp2 Fp2_basis_prime_10_div_3_1;
struct Fp2 Fp2_basis_prime_10_div_3_2;
struct Fp2 Fp2_basis_prime_10_div_6;
struct Fp2 Fp2_basis_inv_prime_1_div_3;
struct Fp2 Fp2_basis_inv_prime_1_div_2;
struct Fp2 Fp2_basis_inv_prime_2_div_3;
struct Fp2 Fp2_basis_inv_prime_2_div_2;
struct Fp2 Fp2_basis_inv_prime_3_div_3;
struct Fp2 Fp2_basis_inv_prime_3_div_2;
struct Fp2 Fp2_basis_inv_prime_10_div_3;
struct Fp2 Fp2_basis_inv_prime_10_div_2;

void init_parameters(){
    //parameters
    mpz_init(prime);
    mpz_init(mother_parameter);
    mpz_init(trace_t);
    mpz_init(EFp_order);
    mpz_init(EFp_total);
    mpz_init(EFp2_total);
    mpz_init(EFp6_total);
    mpz_init(EFp12_total);
    mpz_init(curve_parameter_A);
    mpz_init(curve_parameter_B);
    
    //basis
    Fp_init(&Fp_basis);
    Fp2_init(&Fp2_basis_inv);
    Fp2_init(&Fp2_basis);
    Fp6_init(&Fp6_basis);
    
    //precomputed
    Fp_init(&inv_CNR1);
    Fp_init(&inv_CNR2);
    Fp_init(&epsilon_1);
    Fp_init(&epsilon_2);
    //frobenius
    Fp2_init(&Fp2_basis_prime_1_div_3_1);
    Fp2_init(&Fp2_basis_prime_1_div_3_2);
    Fp2_init(&Fp2_basis_prime_1_div_6);
    Fp2_init(&Fp2_basis_prime_2_div_3_1);
    Fp2_init(&Fp2_basis_prime_2_div_3_2);
    Fp2_init(&Fp2_basis_prime_2_div_6);
    Fp2_init(&Fp2_basis_prime_3_div_3_1);
    Fp2_init(&Fp2_basis_prime_3_div_3_2);
    Fp2_init(&Fp2_basis_prime_3_div_6);
    Fp2_init(&Fp2_basis_prime_4_div_3_1);
    Fp2_init(&Fp2_basis_prime_4_div_3_2);
    Fp2_init(&Fp2_basis_prime_4_div_6);
    Fp2_init(&Fp2_basis_prime_8_div_3_1);
    Fp2_init(&Fp2_basis_prime_8_div_3_2);
    Fp2_init(&Fp2_basis_prime_8_div_6);
    Fp2_init(&Fp2_basis_prime_10_div_3_1);
    Fp2_init(&Fp2_basis_prime_10_div_3_2);
    Fp2_init(&Fp2_basis_prime_10_div_6);
    //skew_frobenius
    Fp2_init(&Fp2_basis_inv_prime_1_div_3);
    Fp2_init(&Fp2_basis_inv_prime_1_div_2);
    Fp2_init(&Fp2_basis_inv_prime_2_div_3);
    Fp2_init(&Fp2_basis_inv_prime_2_div_2);
    Fp2_init(&Fp2_basis_inv_prime_3_div_3);
    Fp2_init(&Fp2_basis_inv_prime_3_div_2);
    Fp2_init(&Fp2_basis_inv_prime_10_div_3);
    Fp2_init(&Fp2_basis_inv_prime_10_div_2);
    
    mpz_init(final_exp);
    memset(X_bit_binary,0,sizeof(X_bit_binary));
    memset(X_bit_binary_for_opt_ate,0,sizeof(X_bit_binary_for_opt_ate));
    
    //ZERO
    Fp_init(&Fp_ZERO);
    Fp2_init(&Fp2_ZERO);
    Fp6_init(&Fp6_ZERO);
    Fp12_init(&Fp12_ZERO);
    //set ZERO
    Fp_set_ui(&Fp_ZERO,0);
    Fp2_set_ui(&Fp2_ZERO,0);
    Fp6_set_ui(&Fp6_ZERO,0);
    Fp12_set_ui(&Fp12_ZERO,0);
}

void clear_parameters(){
    //parameters
    mpz_clear(prime);
    mpz_clear(mother_parameter);
    mpz_clear(trace_t);
    mpz_clear(EFp_order);
    mpz_clear(EFp_total);
    mpz_clear(EFp2_total);
    mpz_clear(EFp6_total);
    mpz_clear(EFp12_total);
    mpz_clear(curve_parameter_A);
    mpz_clear(curve_parameter_B);
    
    //basis
    Fp_clear(&Fp_basis);
    Fp2_clear(&Fp2_basis_inv);
    Fp2_clear(&Fp2_basis);
    Fp6_clear(&Fp6_basis);
    //precomputed
    Fp_clear(&inv_CNR1);
    Fp_clear(&inv_CNR2);
    Fp_clear(&epsilon_1);
    Fp_clear(&epsilon_2);
    //frobenius
    Fp2_clear(&Fp2_basis_prime_1_div_3_1);
    Fp2_clear(&Fp2_basis_prime_1_div_3_2);
    Fp2_clear(&Fp2_basis_prime_1_div_6);
    Fp2_clear(&Fp2_basis_prime_2_div_3_1);
    Fp2_clear(&Fp2_basis_prime_2_div_3_2);
    Fp2_clear(&Fp2_basis_prime_2_div_6);
    Fp2_clear(&Fp2_basis_prime_3_div_3_1);
    Fp2_clear(&Fp2_basis_prime_3_div_3_2);
    Fp2_clear(&Fp2_basis_prime_3_div_6);
    Fp2_clear(&Fp2_basis_prime_4_div_3_1);
    Fp2_clear(&Fp2_basis_prime_4_div_3_2);
    Fp2_clear(&Fp2_basis_prime_4_div_6);
    Fp2_clear(&Fp2_basis_prime_8_div_3_1);
    Fp2_clear(&Fp2_basis_prime_8_div_3_2);
    Fp2_clear(&Fp2_basis_prime_8_div_6);
    Fp2_clear(&Fp2_basis_prime_10_div_3_1);
    Fp2_clear(&Fp2_basis_prime_10_div_3_2);
    Fp2_clear(&Fp2_basis_prime_10_div_6);
    //skew_frobenius
    Fp2_clear(&Fp2_basis_inv_prime_1_div_3);
    Fp2_clear(&Fp2_basis_inv_prime_1_div_2);
    Fp2_clear(&Fp2_basis_inv_prime_2_div_3);
    Fp2_clear(&Fp2_basis_inv_prime_2_div_2);
    Fp2_clear(&Fp2_basis_inv_prime_3_div_3);
    Fp2_clear(&Fp2_basis_inv_prime_3_div_2);
    Fp2_clear(&Fp2_basis_inv_prime_10_div_3);
    Fp2_clear(&Fp2_basis_inv_prime_10_div_2);
    
    //ZERO
    Fp_clear(&Fp_ZERO);
    Fp2_clear(&Fp2_ZERO);
    Fp6_clear(&Fp6_ZERO);
    Fp12_clear(&Fp12_ZERO);
}

void set_parameters(){
    mpz_t result;
    mpz_init(result);
    
    //set curve parameter
    mpz_set_ui(curve_parameter_A,0);
    mpz_set_ui(curve_parameter_B,4);
    
    //generate mother_parameter
    generate_mother_parameter();
    sign=1;				//sign of mother_parameter
    
    //generate prime,order,trace
    if(generate_prime()==1 && generate_order()==1){
        generate_trace();
    }else{
        printf("This mother parameter cannot use.\n");
        clear_parameters();
        exit(1);
    }
    
    weil();				//total rational point
    get_epsilon();			//calculate 1^(1/3)
    generate_basis();			//set basis
    get_scalar_of_final_exp();	//get value
    
    mpz_clear(result);
}

void generate_mother_parameter(){
    int i;
    mpz_t buf,set_2;
    mpz_init(buf);
    mpz_init(set_2);
    mpz_set_ui(set_2,2);
    
    //X_bit_binary
    X_bit_binary[114]=1;
    X_bit_binary[101]=1;
    X_bit_binary[14]=-1;
    X_bit_binary[0]=-1;
    //X_bit_binary_for_opt_ate
    X_bit_binary_for_opt_ate[116]=1;
    X_bit_binary_for_opt_ate[115]=1;
    X_bit_binary_for_opt_ate[103]=1;
    X_bit_binary_for_opt_ate[102]=1;
    X_bit_binary_for_opt_ate[16]=-1;
    X_bit_binary_for_opt_ate[15]=-1;
    X_bit_binary_for_opt_ate[2]=-1;
    
    //mother_parameter
    mpz_set_ui(mother_parameter,0);
    for(i=x_bit; i>=0; i--){
        if(X_bit_binary[i]==1){
            mpz_pow_ui(buf,set_2,i);
            mpz_add(mother_parameter,mother_parameter,buf);
        }else if(X_bit_binary[i]==-1){
            mpz_pow_ui(buf,set_2,i);
            mpz_sub(mother_parameter,mother_parameter,buf);
        }
    }
    
    mpz_clear(buf);
    mpz_clear(set_2);
}

int generate_prime(){
    mpz_t buf1,buf2,base;
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(base);
    mpz_set(base,mother_parameter);
    
    //prime
    mpz_pow_ui(buf1,base,4);
    mpz_mul_ui(buf1,buf1,36);
    mpz_set(buf2,buf1);
    mpz_pow_ui(buf1,base,3);
    mpz_mul_ui(buf1,buf1,36);
    if(sign==1){
        mpz_add(buf2,buf2,buf1);
    }else{
        mpz_sub(buf2,buf2,buf1);
    }
    mpz_pow_ui(buf1,base,2);
    mpz_mul_ui(buf1,buf1,24);
    mpz_add(buf2,buf2,buf1);
    mpz_mul_ui(buf1,base,6);
    if(sign==1){
        mpz_add(buf2,buf2,buf1);
    }else{
        mpz_sub(buf2,buf2,buf1);
    }
    mpz_add_ui(buf2,buf2,1);
    
    //isprime
    if(mpz_probab_prime_p(buf2,25)==0){
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(base);
        return 0;
    }else{
        mpz_set(prime,buf2);
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(base);
        return 1;
    }
}

int generate_order(){
    mpz_t buf1,buf2,base;
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(base);
    mpz_set(base,mother_parameter);
    
    //order
    mpz_pow_ui(buf1,base,4);
    mpz_mul_ui(buf1,buf1,36);
    mpz_set(buf2,buf1);
    mpz_pow_ui(buf1,base,3);
    mpz_mul_ui(buf1,buf1,36);
    if(sign==1){
        mpz_add(buf2,buf2,buf1);
    }else{
        mpz_sub(buf2,buf2,buf1);
    }
    mpz_pow_ui(buf1,base,2);
    mpz_mul_ui(buf1,buf1,18);
    mpz_add(buf2,buf2,buf1);
    mpz_mul_ui(buf1,base,6);
    if(sign==1){
        mpz_add(buf2,buf2,buf1);
    }else{
        mpz_sub(buf2,buf2,buf1);
    }
    mpz_add_ui(buf2,buf2,1);
    
    //isprime
    if(mpz_probab_prime_p(buf2,25)==0){
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(base);
        return 0;
    }else{
        mpz_set(EFp_order,buf2);
        mpz_clear(buf1);
        mpz_clear(buf2);
        mpz_clear(base);
        return 1;
    }
}

void generate_trace(){
    mpz_t buf,base;
    mpz_init(buf);
    mpz_init(base);
    
    mpz_set(base,mother_parameter);
    mpz_pow_ui(buf,base,2);
    mpz_mul_ui(buf,buf,6);
    mpz_add_ui(trace_t,buf,1);
    
    mpz_clear(buf);
    mpz_clear(base);
}

void generate_basis(){
    unsigned long int i;
    
    //Fp_basis
    i=0;
    while(mpz_cmp_ui(prime,i)>0){
        Fp_set_ui(&Fp_basis,i);
        if(Fp_legendre(&Fp_basis)==1){
            break;
        }
        i++;
    }
    
    //Fp2_basis
    Fp2_set_ui(&Fp2_basis,0);
    Fp_set_ui(&Fp2_basis.x1,1);
    i=0;
    while(mpz_cmp_ui(prime,i)>0){
        Fp_set_ui(&Fp2_basis.x0,i);
        if(Fp2_legendre(&Fp2_basis)==-1){
            break;
        }
        i++;
    }
    //Fp2_basis_inv
    Fp2_inv(&Fp2_basis_inv,&Fp2_basis);
    
    //Fp12_basis
    Fp6_set_ui(&Fp6_basis,0);
    Fp_set_ui(&Fp6_basis.x1.x0,1);
    i=0;
    while(mpz_cmp_ui(prime,i)>0){
        Fp_set_ui(&Fp6_basis.x0.x0,i);
        if(Fp6_legendre(&Fp6_basis)==-1){
            break;
        }
        i++;
    }
    
}

void get_epsilon(){
    struct Fp inv,root,buf;
    Fp_init(&inv);
    Fp_init(&root);
    Fp_init(&buf);
    mpz_t exp;
    mpz_init(exp);
    
    Fp_set_ui(&buf,2);
    Fp_inv(&inv,&buf);
    mpz_sub_ui(buf.x0,prime,3);
    Fp_sqrt(&root,&buf);
    Fp_sub_ui(&buf,&root,1);
    Fp_mul(&inv_CNR1,&buf,&inv);
    Fp_mul(&inv_CNR2,&inv_CNR1,&inv_CNR1);
    
    mpz_clear(exp);
    Fp_clear(&inv);
    Fp_clear(&root);
    Fp_clear(&buf);
}

void get_scalar_of_final_exp(){
    struct Fp2 Buf;
    Fp2_init(&Buf);
    
    mpz_t exp,buf,p2,p3,p4,p6,p8,p10;
    mpz_init(exp);
    mpz_init(buf);
    mpz_init(p2);
    mpz_init(p3);
    mpz_init(p4);
    mpz_init(p6);
    mpz_init(p8);
    mpz_init(p10);
    
    mpz_mul(p2,prime,prime);
    mpz_mul(p3,p2,prime);
    mpz_mul(p4,p3,prime);
    mpz_mul(p6,p4,p2);
    mpz_mul(p8,p6,p2);
    mpz_mul(p10,p8,p2);
    
    //frobenius_1
    mpz_sub_ui(exp,prime,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&Fp2_basis_prime_1_div_3_1,&Fp2_basis,exp);
    Fp2_mul(&Fp2_basis_prime_1_div_3_2,&Fp2_basis_prime_1_div_3_1,&Fp2_basis_prime_1_div_3_1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_prime_1_div_6,&Fp2_basis,exp);
    //Fp2_printf(&Fp2_basis_prime_1_div_3_1,""); printf("\n");
    //Fp2_printf(&Fp2_basis_prime_1_div_3_2,""); printf("\n");
    //Fp2_printf(&Fp2_basis_prime_1_div_6,""); printf("\n");
    
    
    //frobenius_2
    mpz_sub_ui(exp,p2,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&Fp2_basis_prime_2_div_3_1,&Fp2_basis,exp);
    Fp2_mul(&Fp2_basis_prime_2_div_3_2,&Fp2_basis_prime_2_div_3_1,&Fp2_basis_prime_2_div_3_1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_prime_2_div_6,&Fp2_basis,exp);
    /*Fp2_printf(&Fp2_basis_prime_2_div_3_1,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_2_div_3_2,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_2_div_6,""); printf("\n");
     */
    
    //frobenius_3
    mpz_sub_ui(exp,p3,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&Fp2_basis_prime_3_div_3_1,&Fp2_basis,exp);
    Fp2_mul(&Fp2_basis_prime_3_div_3_2,&Fp2_basis_prime_3_div_3_1,&Fp2_basis_prime_3_div_3_1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_prime_3_div_6,&Fp2_basis,exp);
    /*Fp2_printf(&Fp2_basis_prime_3_div_3_1,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_3_div_3_2,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_3_div_6,""); printf("\n");
     */
    
    //frobenius_4
    mpz_sub_ui(exp,p4,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&Fp2_basis_prime_4_div_3_1,&Fp2_basis,exp);
    Fp2_mul(&Fp2_basis_prime_4_div_3_2,&Fp2_basis_prime_4_div_3_1,&Fp2_basis_prime_4_div_3_1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_prime_4_div_6,&Fp2_basis,exp);
    /*Fp2_printf(&Fp2_basis_prime_4_div_3_1,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_4_div_3_2,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_4_div_6,""); printf("\n");
     */
    
    //frobenius_8
    mpz_sub_ui(exp,p8,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&Fp2_basis_prime_8_div_3_1,&Fp2_basis,exp);
    Fp2_mul(&Fp2_basis_prime_8_div_3_2,&Fp2_basis_prime_8_div_3_1,&Fp2_basis_prime_8_div_3_1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_prime_8_div_6,&Fp2_basis,exp);
    /*Fp2_printf(&Fp2_basis_prime_8_div_3_1,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_8_div_3_2,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_8_div_6,""); printf("\n");
     */
    
    //frobenius_10
    mpz_sub_ui(exp,p10,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&Fp2_basis_prime_10_div_3_1,&Fp2_basis,exp);
    Fp2_mul(&Fp2_basis_prime_10_div_3_2,&Fp2_basis_prime_10_div_3_1,&Fp2_basis_prime_10_div_3_1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_prime_10_div_6,&Fp2_basis,exp);
    /*Fp2_printf(&Fp2_basis_prime_10_div_3_1,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_10_div_3_2,""); printf("\n");
     Fp2_printf(&Fp2_basis_prime_10_div_6,""); printf("\n");
     */
    
    //skew_frobenius_1
    Fp2_inv(&Fp2_basis_inv_prime_1_div_3,&Fp2_basis_prime_1_div_3_1);
    mpz_sub_ui(exp,prime,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_inv_prime_1_div_2,&Fp2_basis,exp);
    Fp2_inv(&Fp2_basis_inv_prime_1_div_2,&Fp2_basis_inv_prime_1_div_2);
    
    //skew_frobenius_2
    Fp2_inv(&Fp2_basis_inv_prime_2_div_3,&Fp2_basis_prime_2_div_3_1);
    mpz_sub_ui(exp,p2,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_inv_prime_2_div_2,&Fp2_basis,exp);
    Fp2_inv(&Fp2_basis_inv_prime_2_div_2,&Fp2_basis_inv_prime_2_div_2);
    
    //skew_frobenius_3
    Fp2_inv(&Fp2_basis_inv_prime_3_div_3,&Fp2_basis_prime_3_div_3_1);
    mpz_sub_ui(exp,p3,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_inv_prime_3_div_2,&Fp2_basis,exp);
    Fp2_inv(&Fp2_basis_inv_prime_3_div_2,&Fp2_basis_inv_prime_3_div_2);
    
    //skew_frobenius_10
    Fp2_inv(&Fp2_basis_inv_prime_10_div_3,&Fp2_basis_prime_10_div_3_1);
    mpz_sub_ui(exp,p10,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&Fp2_basis_inv_prime_10_div_2,&Fp2_basis,exp);
    Fp2_inv(&Fp2_basis_inv_prime_10_div_2,&Fp2_basis_inv_prime_10_div_2);
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp2_clear(&Buf);
    mpz_clear(p2);
    mpz_clear(p3);
    mpz_clear(p4);
    mpz_clear(p6);
    mpz_clear(p8);
    mpz_clear(p10);
}

void weil(){
    mpz_t t2,t6,t12,p_exp_2,p_exp_6,buf;
    mpz_init(t2);
    mpz_init(t6);
    mpz_init(t12);
    mpz_init(p_exp_2);
    mpz_init(p_exp_6);
    mpz_init(buf);
    
    //EFp_total
    mpz_add_ui(buf,prime,1);
    mpz_sub(EFp_total,buf,trace_t);
    
    //t2←α^2+β^2
    mpz_pow_ui(t2,trace_t,2);
    mpz_mul_ui(buf,prime,2);
    mpz_sub(t2,t2,buf);
    //EFp2_total
    mpz_pow_ui(p_exp_2,prime,2);
    mpz_sub(buf,p_exp_2,t2);
    mpz_add_ui(EFp2_total,buf,1);
    
    //α^6+β^6
    mpz_pow_ui(t6,t2,3);
    mpz_mul(buf,t2,p_exp_2);
    mpz_mul_ui(buf,buf,3);
    mpz_sub(t6,t6,buf);
    //EFp6_total
    mpz_pow_ui(p_exp_6,p_exp_2,3);
    mpz_sub(buf,p_exp_6,t6);
    mpz_add_ui(EFp6_total,buf,1);
    
    //α^12+β^12
    mpz_pow_ui(t12,t6,2);
    mpz_mul_ui(buf,p_exp_6,2);
    mpz_sub(t12,t12,buf);
    //EFp12_total
    mpz_pow_ui(buf,p_exp_6,2);
    mpz_sub(buf,buf,t12);
    mpz_add_ui(EFp12_total,buf,1);
    
    mpz_clear(t2);
    mpz_clear(t6);
    mpz_clear(t12);
    mpz_clear(p_exp_2);
    mpz_clear(p_exp_6);
    mpz_clear(buf);
}

void print_parameters(){
    printf("====================================================================================\n");
    printf("prime length:%dbit\n",(int)mpz_sizeinbase(prime,2));
    printf("E:y^2=x^3-");
    if(mpz_cmp_ui(curve_parameter_A,0)!=0){
        mpz_out_str(stdout,10,curve_parameter_A);
        printf("x+");
    }
    if(mpz_cmp_ui(curve_parameter_B,0)!=0){
        mpz_out_str(stdout,10,curve_parameter_B);
    }
    printf("(mod");
    mpz_out_str(stdout,10,prime);
    printf(")\n");
    printf("mother parameter:");
    mpz_out_str(stdout,10,mother_parameter);
    printf("\n");
    printf("Fp2:f(x)=x^2+");
    Fp_printf(&Fp_basis,"");
    printf("\n");
    printf("Fp6:f(x)=x^3-");
    Fp2_printf(&Fp2_basis,"");
    printf("\n");
    printf("Fp12:f(x)=x^2-");
    Fp6_printf(&Fp6_basis,"");
    printf("\n");
    
    printf("EFp_order:");
    mpz_out_str(stdout,10,EFp_order);
    printf("\n");
    printf("trace_t=");
    mpz_out_str(stdout,10,trace_t);
    printf("\n");
    
}
