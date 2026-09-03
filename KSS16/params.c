#include "params.h"

char X_bit_binary[x_bit+1];
mpz_t X;
mpz_t PRIME_P,order_r,trace_t, order_EFp, a_x;
mpz_t tmp_a;

#include "curve.h"
#include "util.h"

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
