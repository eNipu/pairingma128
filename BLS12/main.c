#include "BLS_12.h"

unsigned long int num,mpz_mpz_mul,mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num;
/*============================================================================*/
/* main                                                                       */
/*============================================================================*/
int main(){
    init_parameters();
    set_parameters();
    print_parameters();
    
    test_opt_ate_pairing();
    
    clear_parameters();
    return 0;
}
/*============================================================================*/
/* Fp                                                                         */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void Fp_init(struct Fp *P){
    mpz_init(P->x0);
}
/*---------------------------------set----------------------------------*/
void Fp_set(struct Fp *P,struct Fp *A){
    mpz_set(P->x0,A->x0);
}
void Fp_set_ui(struct Fp *P,unsigned long int a){
    mpz_set_ui(P->x0,a);
}
void Fp_set_mpz(struct Fp *P,mpz_t a){
    mpz_set(P->x0,a);
}
void Fp_set_neg(struct Fp *P,struct Fp *A){
    mpz_sub(P->x0,prime,A->x0);
}
/*---------------------------------random--------------------------------*/
void Fp_random(struct Fp *P,gmp_randstate_t state){
    mpz_urandomm(P->x0,state,prime);
}
/*---------------------------------clear---------------------------------*/
void Fp_clear(struct Fp *P){
    mpz_clear(P->x0);
}
/*---------------------------------print---------------------------------*/
void Fp_printf(struct Fp *P,char *name){
    printf("%s",name);
    mpz_out_str(stdout,10,P->x0);
}
/*---------------------------vector calculation--------------------------*/
void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        Fp_mpz_sqr++;
    }else{
        mpz_mpz_mul++;
    }
    mpz_mul(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_mul_ui(struct Fp *ANS,struct Fp *A,unsigned long int a){
    mpz_ui_mul++;
    mpz_mul_ui(ANS->x0,A->x0,a);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_mul_mpz(struct Fp *ANS,struct Fp *A,mpz_t a){
    mpz_mpz_mul++;
    mpz_mul(ANS->x0,A->x0,a);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_mul_basis(struct Fp *ANS,struct Fp *A){
    basis_mul_num++;
    mpz_sub(ANS->x0,prime,A->x0);
}
void Fp_add(struct Fp *ANS,struct Fp *A,struct Fp *B){
    mpz_mpz_add++;
    mpz_add(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_add_ui(struct Fp *ANS,struct Fp *A,unsigned long int a){
    mpz_ui_add++;
    mpz_add_ui(ANS->x0,A->x0,a);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_add_mpz(struct Fp *ANS,struct Fp *A,mpz_t a){
    mpz_add(ANS->x0,A->x0,a);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_sub(struct Fp *ANS,struct Fp *A,struct Fp *B){
    mpz_mpz_add++;
    mpz_sub(ANS->x0,A->x0,B->x0);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_sub_ui(struct Fp *ANS,struct Fp *A,unsigned long int a){
    mpz_ui_add++;
    mpz_sub_ui(ANS->x0,A->x0,a);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
void Fp_sub_mpz(struct Fp *ANS,struct Fp *A,mpz_t a){
    mpz_sub(ANS->x0,A->x0,a);
    mpz_mod(ANS->x0,ANS->x0,prime);
}
/*------------------------------inverse---------------------------------*/
void Fp_inv(struct Fp *ANS,struct Fp *A){
    Fp_inv_num++;
    mpz_invert(ANS->x0,A->x0,prime);
}
/*------------------------------legendre-------------------------------*/
int Fp_legendre(struct Fp *A){
    return mpz_legendre(A->x0,prime);
}
int Fp_isCNR(struct Fp *A){
    struct Fp buf;
    Fp_init(&buf);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_sub_ui(exp,prime,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp_pow(&buf,A,exp);
    
    mpz_clear(exp);
    if(Fp_cmp_one(&buf)==0){
        Fp_clear(&buf);
        return 1;
    }else if(Fp_cmp_zero(&buf)==0){
        Fp_clear(&buf);
        return 0;
    }else{
        Fp_clear(&buf);
        return -1;
    }
}
/*---------------------------------sqr----------------------------------*/
void Fp_sqrt(struct Fp *ANS,struct Fp *A){
    struct Fp x,y,t,k,n,buf;
    Fp_init(&x);
    Fp_init(&y);
    Fp_init(&t);
    Fp_init(&k);
    Fp_init(&n);
    Fp_init(&buf);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp_random(&n,state);
    
    while(Fp_legendre(&n)!=-1){
        Fp_random(&n,state);
    }
    mpz_sub_ui(q,prime,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp_pow(&x,A,exp);
    Fp_mul(&buf,&x,&x);
    Fp_mul(&k,&buf,A);
    Fp_mul(&x,&x,A);
    while(Fp_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp_pow(&buf,&k,exp);
        while(Fp_cmp_one(&buf)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp_pow(&buf,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp_pow(&t,&y,result);
        Fp_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp_mul(&x,&x,&t);
        Fp_mul(&k,&k,&y);
    }
    Fp_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp_clear(&x);
    Fp_clear(&y);
    Fp_clear(&t);
    Fp_clear(&k);
    Fp_clear(&n);
    Fp_clear(&buf);
}
/*---------------------------------pow----------------------------------*/
void Fp_pow(struct Fp *ANS,struct Fp *A,mpz_t a){
    int i,length;
    length=(int)mpz_sizeinbase(a,2);
    char binary[length];
    mpz_get_str(binary,2,a);
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
/*---------------------------------cmp----------------------------------*/
int Fp_cmp(struct Fp *A,struct Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        return 0;
    }
    return 1;
}
int Fp_cmp_ui(struct Fp *A,unsigned long int a){
    if(mpz_cmp_ui(A->x0,a)==0){
        return 0;
    }
    return 1;
}
int Fp_cmp_mpz(struct Fp *A,mpz_t a){
    if(mpz_cmp(A->x0,a)==0){
        return 0;
    }
    return 1;
}
int Fp_cmp_zero(struct Fp *A){
    if(mpz_cmp_ui(A->x0,0)==0){
        return 0;
    }
    return 1;
}
int Fp_cmp_one(struct Fp *A){
    if(mpz_cmp_ui(A->x0,1)==0){
        return 0;
    }
    return 1;
}
/*============================================================================*/
/* Fp2                                                                        */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void Fp2_init(struct Fp2 *P){
    Fp_init(&P->x0);
    Fp_init(&P->x1);
}
/*---------------------------------set----------------------------------*/
void Fp2_set(struct Fp2 *P,struct Fp2 *A){
    Fp_set(&P->x0,&A->x0);
    Fp_set(&P->x1,&A->x1);
}
void Fp2_set_ui(struct Fp2 *P,unsigned long int a){
    Fp_set_ui(&P->x0,a);
    Fp_set_ui(&P->x1,a);
}
void Fp2_set_mpz(struct Fp2 *P,mpz_t a){
    Fp_set_mpz(&P->x0,a);
    Fp_set_mpz(&P->x1,a);
}
void Fp2_set_neg(struct Fp2 *P,struct Fp2 *A){
    Fp_set_neg(&P->x0,&A->x0);
    Fp_set_neg(&P->x1,&A->x1);
}
/*---------------------------------random--------------------------------*/
void Fp2_random(struct Fp2 *P,gmp_randstate_t state){
    Fp_random(&P->x0,state);
    Fp_random(&P->x1,state);
}
/*---------------------------------clear---------------------------------*/
void Fp2_clear(struct Fp2 *P){
    Fp_clear(&P->x0);
    Fp_clear(&P->x1);
}
/*---------------------------------print---------------------------------*/
void Fp2_printf(struct Fp2 *P,char *name){
    printf("%s",name);
    printf("(");
    Fp_printf(&P->x0,"");
    printf(",");
    Fp_printf(&P->x1,"");
    printf(")");
}
/*---------------------------vector calculation--------------------------*/
void Fp2_mul(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    struct Fp tmp1,tmp2,tmp3,tmp4;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    Fp_init(&tmp3);
    Fp_init(&tmp4);
    
    //set
    Fp_mul(&tmp1,&A->x0,&B->x0);//a*c
    Fp_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp_add(&tmp3,&A->x0,&A->x1);//a+b
    Fp_add(&tmp4,&B->x0,&B->x1);//c+d
    //x0
    Fp_mul_basis(&ANS->x0,&tmp2);//b*d*v
    Fp_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp_mul(&ANS->x1,&tmp3,&tmp4);//(a+b)(c+d)
    Fp_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
    Fp_clear(&tmp3);
    Fp_clear(&tmp4);
}
void Fp2_mul_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a){
    Fp_mul_ui(&ANS->x0,&A->x0,a);
    Fp_mul_ui(&ANS->x1,&A->x1,a);
}
void Fp2_mul_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a){
    Fp_mul_mpz(&ANS->x0,&A->x0,a);
    Fp_mul_mpz(&ANS->x1,&A->x1,a);
}
void Fp2_mul_basis(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp tmp;
    Fp_init(&tmp);
    Fp_set(&tmp,&A->x0);
    
    Fp_sub(&ANS->x0,&tmp,&A->x1);
    Fp_add(&ANS->x1,&tmp,&A->x1);
    
    Fp_clear(&tmp);
}
void Fp2_squaring(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp tmp1,tmp2;
    Fp_init(&tmp1);
    Fp_init(&tmp2);
    
    Fp_add(&tmp1,&A->x0,&A->x1);
    Fp_sub(&tmp2,&A->x0,&A->x1);
    //x1
    Fp_mul(&ANS->x1,&A->x0,&A->x1);
    Fp_add(&ANS->x1,&ANS->x1,&ANS->x1);
    //x0
    Fp_mul(&ANS->x0,&tmp1,&tmp2);
    
    Fp_clear(&tmp1);
    Fp_clear(&tmp2);
}
void Fp2_inv_basis(struct Fp2 *ANS,struct Fp2 *A){
    //Fp2_mul(ANS,A,&Fp2_basis_inv);
    struct Fp2 tmp;
    Fp2_init(&tmp);
    Fp2_set(&tmp,A);
    
    Fp_add(&ANS->x0,&tmp.x0,&tmp.x1);
    Fp_mul_mpz(&ANS->x0,&ANS->x0,Fp2_basis_inv.x0.x0);
    Fp_sub(&ANS->x1,&tmp.x1,&tmp.x0);
    Fp_mul_mpz(&ANS->x1,&ANS->x1,Fp2_basis_inv.x0.x0);
    
    Fp2_clear(&tmp);
}
void Fp2_add(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    Fp_add(&ANS->x0,&A->x0,&B->x0);
    Fp_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp2_add_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a){
    Fp_add_ui(&ANS->x0,&A->x0,a);
    Fp_add_ui(&ANS->x1,&A->x1,a);
}
void Fp2_add_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a){
    Fp_add_mpz(&ANS->x0,&A->x0,a);
    Fp_add_mpz(&ANS->x1,&A->x1,a);
}
void Fp2_sub(struct Fp2 *ANS,struct Fp2 *A,struct Fp2 *B){
    Fp_sub(&ANS->x0,&A->x0,&B->x0);
    Fp_sub(&ANS->x1,&A->x1,&B->x1);
}
void Fp2_sub_ui(struct Fp2 *ANS,struct Fp2 *A,unsigned long int a){
    Fp_sub_ui(&ANS->x0,&A->x0,a);
    Fp_sub_ui(&ANS->x1,&A->x1,a);
}
void Fp2_sub_mpz(struct Fp2 *ANS,struct Fp2 *A,mpz_t a){
    Fp_sub_mpz(&ANS->x0,&A->x0,a);
    Fp_sub_mpz(&ANS->x1,&A->x1,a);
}
/*------------------------------inverse---------------------------------*/
void Fp2_inv(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 frob,buf;
    Fp2_init(&frob);
    Fp2_init(&buf);
    
    Fp2_inv_map(&frob,A);
    Fp2_mul(&buf,A,&frob);
    Fp_inv(&buf.x0,&buf.x0);
    Fp2_mul_mpz(ANS,&frob,buf.x0.x0);
    
    Fp2_clear(&frob);
    Fp2_clear(&buf);
}
void Fp2_inv_map(struct Fp2 *ANS,struct Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
}
/*------------------------------legendre-------------------------------*/
int Fp2_legendre(struct Fp2 *A){
    struct Fp2 buf;
    Fp2_init(&buf);
    
    mpz_t exp;
    mpz_init(exp);
    mpz_pow_ui(exp,prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&buf,A,exp);
    
    mpz_clear(exp);
    if(Fp2_cmp_one(&buf)==0){
        Fp2_clear(&buf);
        return 1;
    }else if(Fp2_cmp_zero(&buf)==0){
        Fp2_clear(&buf);
        return 0;
    }else{
        Fp2_clear(&buf);
        return -1;
    }
}
int Fp2_isCNR(struct Fp2 *A){
    struct Fp2 buf;
    Fp2_init(&buf);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,2);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp2_pow(&buf,A,exp);
    
    mpz_clear(exp);
    if(Fp2_cmp_one(&buf)==0){
        Fp2_clear(&buf);
        return 1;
    }else if(Fp2_cmp(&buf,&Fp2_ZERO)==0){
        Fp2_clear(&buf);
        return 0;
    }else{
        Fp2_clear(&buf);
        return -1;
    }
}
/*---------------------------------sqr----------------------------------*/
void Fp2_sqrt(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 x,y,t,k,n,buf;
    Fp2_init(&x);
    Fp2_init(&y);
    Fp2_init(&t);
    Fp2_init(&k);
    Fp2_init(&n);
    Fp2_init(&buf);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp2_random(&n,state);
    while(Fp2_legendre(&n)!=-1){
        Fp2_random(&n,state);
    }
    mpz_pow_ui(q,prime,2);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp2_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp2_pow(&x,A,exp);
    Fp2_mul(&buf,&x,&x);
    Fp2_mul(&k,&buf,A);
    Fp2_mul(&x,&x,A);
    while(Fp2_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp2_pow(&buf,&k,exp);
        while(Fp2_cmp_one(&buf)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp2_pow(&buf,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp2_pow(&t,&y,result);
        Fp2_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp2_mul(&x,&x,&t);
        Fp2_mul(&k,&k,&y);
    }
    Fp2_set(ANS,&x);
    
    mpz_clear(exp);
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(result);
    Fp2_clear(&x);
    Fp2_clear(&y);
    Fp2_clear(&t);
    Fp2_clear(&k);
    Fp2_clear(&n);
    Fp2_clear(&buf);
}
/*---------------------------------pow----------------------------------*/
void Fp2_pow(struct Fp2 *ANS,struct Fp2 *A,mpz_t a){
    int i,length;
    length=(int)mpz_sizeinbase(a,2);
    char binary[length];
    mpz_get_str(binary,2,a);
    struct Fp2 buf;
    Fp2_init(&buf);
    Fp2_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp2_squaring(&buf,&buf);
        if(binary[i]=='1'){
            Fp2_mul(&buf,A,&buf);
        }
    }
    
    Fp2_set(ANS,&buf);
    Fp2_clear(&buf);
}
/*---------------------------------cmp----------------------------------*/
int Fp2_cmp(struct Fp2 *A,struct Fp2 *B){
    if(Fp_cmp(&A->x0,&B->x0)==0 && Fp_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp2_cmp_ui(struct Fp2 *A,unsigned long int a){
    if(Fp_cmp_ui(&A->x0,a)==0 && Fp_cmp_ui(&A->x1,a)==0){
        return 0;
    }
    return 1;
}
int Fp2_cmp_mpz(struct Fp2 *A,mpz_t a){
    if(Fp_cmp_mpz(&A->x0,a)==0 && Fp_cmp_mpz(&A->x1,a)==0){
        return 0;
    }
    return 1;
}
int Fp2_cmp_zero(struct Fp2 *A){
    if(Fp_cmp_zero(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
int Fp2_cmp_one(struct Fp2 *A){
    if(Fp_cmp_one(&A->x0)==0 && Fp_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
/*--------------------------------frobenius--------------------------------*/
void Fp2_frobenius_1(struct Fp2 *ANS,struct Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);;
}
void Fp2_frobenius_2(struct Fp2 *ANS,struct Fp2 *A){
    Fp2_set(ANS,A);
}
void Fp2_frobenius_3(struct Fp2 *ANS,struct Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
}
void Fp2_frobenius_4(struct Fp2 *ANS,struct Fp2 *A){
    Fp2_set(ANS,A);
}
void Fp2_frobenius_6(struct Fp2 *ANS,struct Fp2 *A){
    Fp2_set(ANS,A);
}
void Fp2_frobenius_8(struct Fp2 *ANS,struct Fp2 *A){
    Fp2_set(ANS,A);
}
void Fp2_frobenius_10(struct Fp2 *ANS,struct Fp2 *A){
    Fp2_set(ANS,A);
}
/*============================================================================*/
/* Fp6                                                                        */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void Fp6_init(struct Fp6 *P){
    Fp2_init(&P->x0);
    Fp2_init(&P->x1);
    Fp2_init(&P->x2);
}
/*---------------------------------set----------------------------------*/
void Fp6_set(struct Fp6 *P,struct Fp6 *A){
    Fp2_set(&P->x0,&A->x0);
    Fp2_set(&P->x1,&A->x1);
    Fp2_set(&P->x2,&A->x2);
}
void Fp6_set_ui(struct Fp6 *P,unsigned long int a){
    Fp2_set_ui(&P->x0,a);
    Fp2_set_ui(&P->x1,a);
    Fp2_set_ui(&P->x2,a);
}
void Fp6_set_mpz(struct Fp6 *P,mpz_t a){
    Fp2_set_mpz(&P->x0,a);
    Fp2_set_mpz(&P->x1,a);
    Fp2_set_mpz(&P->x2,a);
}
void Fp6_set_neg(struct Fp6 *P,struct Fp6 *A){
    Fp2_set_neg(&P->x0,&A->x0);
    Fp2_set_neg(&P->x1,&A->x1);
    Fp2_set_neg(&P->x2,&A->x2);
}
/*---------------------------------random--------------------------------*/
void Fp6_random(struct Fp6 *P,gmp_randstate_t state){
    Fp2_random(&P->x0,state);
    Fp2_random(&P->x1,state);
    Fp2_random(&P->x2,state);
}
/*---------------------------------clear---------------------------------*/
void Fp6_clear(struct Fp6 *P){
    Fp2_clear(&P->x0);
    Fp2_clear(&P->x1);
    Fp2_clear(&P->x2);
}
/*---------------------------------print---------------------------------*/
void Fp6_printf(struct Fp6 *P,char *name){
    printf("%s",name);
    printf("(");
    Fp2_printf(&P->x0,"");
    printf(",");
    Fp2_printf(&P->x1,"");
    printf(",");
    Fp2_printf(&P->x2,"");
    printf(")");
}
/*---------------------------vector calculation--------------------------*/
void Fp6_mul(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
    struct Fp2 tmp00,tmp11,tmp22,buf,t0,t1,t2;
    Fp2_init(&tmp00);
    Fp2_init(&tmp11);
    Fp2_init(&tmp22);
    Fp2_init(&buf);
    Fp2_init(&t0);
    Fp2_init(&t1);
    Fp2_init(&t2);
    
    //set
    Fp2_mul(&tmp00,&A->x0,&B->x0);//x0*y0
    Fp2_mul(&tmp11,&A->x1,&B->x1);//x1*y1
    Fp2_mul(&tmp22,&A->x2,&B->x2);//x2*y2
    
    Fp2_add(&t0,&A->x0,&A->x1);//x0+x1
    Fp2_add(&buf,&B->x0,&B->x1);//y0+y1
    Fp2_mul(&t0,&t0,&buf);//(x0+x1)(y0+y1)
    
    Fp2_add(&t1,&A->x1,&A->x2);//x1+x2
    Fp2_add(&buf,&B->x1,&B->x2);//y1+y2
    Fp2_mul(&t1,&t1,&buf);//(x1+x2)(y1+y2)
    
    Fp2_add(&t2,&B->x0,&B->x2);//y2+y0
    Fp2_add(&buf,&A->x0,&A->x2);//x2+x0
    Fp2_mul(&t2,&t2,&buf);//(x2+x0)(y2+y0)
    //x0
    Fp2_sub(&t1,&t1,&tmp11);
    Fp2_sub(&t1,&t1,&tmp22);//(x1+x2)(y1+y2)-x1y1-x2y2
    Fp2_mul_basis(&buf,&t1);
    Fp2_add(&ANS->x0,&tmp00,&buf);
    //x1
    Fp2_sub(&t0,&t0,&tmp00);
    Fp2_sub(&t0,&t0,&tmp11);
    Fp2_mul_basis(&buf,&tmp22);
    Fp2_add(&ANS->x1,&buf,&t0);
    //x2
    Fp2_sub(&t2,&t2,&tmp00);
    Fp2_sub(&t2,&t2,&tmp22);
    Fp2_add(&ANS->x2,&tmp11,&t2);
    
    //clear
    Fp2_clear(&tmp00);
    Fp2_clear(&tmp11);
    Fp2_clear(&tmp22);
    Fp2_clear(&buf);
    Fp2_clear(&t0);
    Fp2_clear(&t1);
    Fp2_clear(&t2);
}
void Fp6_mul_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int a){
    Fp2_mul_ui(&ANS->x0,&A->x0,a);
    Fp2_mul_ui(&ANS->x1,&A->x1,a);
    Fp2_mul_ui(&ANS->x2,&A->x2,a);
}
void Fp6_mul_mpz(struct Fp6 *ANS,struct Fp6 *A,mpz_t a){
    Fp2_mul_mpz(&ANS->x0,&A->x0,a);
    Fp2_mul_mpz(&ANS->x1,&A->x1,a);
    Fp2_mul_mpz(&ANS->x2,&A->x2,a);
}
void Fp6_mul_basis(struct Fp6 *ANS,struct Fp6 *A){
    struct Fp6 tmp;
    Fp6_init(&tmp);
    Fp6_set(&tmp,A);
    
    Fp_sub(&ANS->x0.x0,&tmp.x2.x0,&tmp.x2.x1);
    Fp_add(&ANS->x0.x1,&tmp.x2.x0,&tmp.x2.x1);
    Fp_set(&ANS->x1.x0,&tmp.x0.x0);
    Fp_set(&ANS->x1.x1,&tmp.x0.x1);
    Fp_set(&ANS->x2.x0,&tmp.x1.x0);
    Fp_set(&ANS->x2.x1,&tmp.x1.x1);
    
    Fp6_clear(&tmp);
}
void Fp6_squaring(struct Fp6 *ANS,struct Fp6 *A){
    struct Fp2 tmp00,tmp12_2,tmp01_2,tmp22,buf;
    Fp2_init(&tmp00);
    Fp2_init(&tmp22);
    Fp2_init(&tmp12_2);
    Fp2_init(&tmp01_2);
    Fp2_init(&buf);
    
    Fp2_squaring(&tmp00,&A->x0);		//x0^2
    Fp2_squaring(&tmp22,&A->x2);		//x2^2
    Fp2_add(&buf,&A->x1,&A->x1);		//2x1
    Fp2_mul(&tmp12_2,&buf,&A->x2);	//2x1x2
    Fp2_mul(&tmp01_2,&A->x0,&buf);	//2x0x1
    Fp2_add(&buf,&A->x0,&A->x1);		//x0+x1+x2
    Fp2_add(&buf,&buf,&A->x2);
    
    //x0
    Fp2_mul_basis(&ANS->x0,&tmp12_2);
    Fp2_add(&ANS->x0,&ANS->x0,&tmp00);
    //x1
    Fp2_mul_basis(&ANS->x1,&tmp22);
    Fp2_add(&ANS->x1,&ANS->x1,&tmp01_2);
    //x2
    Fp2_squaring(&ANS->x2,&buf);
    Fp2_add(&buf,&tmp00,&tmp22);
    Fp2_add(&buf,&buf,&tmp12_2);
    Fp2_add(&buf,&buf,&tmp01_2);
    Fp2_sub(&ANS->x2,&ANS->x2,&buf);
    
    Fp2_clear(&tmp00);
    Fp2_clear(&tmp22);
    Fp2_clear(&tmp12_2);
    Fp2_clear(&tmp01_2);
    Fp2_clear(&buf);
}
void Fp6_add(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
    Fp2_add(&ANS->x0,&A->x0,&B->x0);
    Fp2_add(&ANS->x1,&A->x1,&B->x1);
    Fp2_add(&ANS->x2,&A->x2,&B->x2);
}
void Fp6_add_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int a){
    Fp2_add_ui(&ANS->x0,&A->x0,a);
    Fp2_add_ui(&ANS->x1,&A->x1,a);
    Fp2_add_ui(&ANS->x2,&A->x2,a);
}
void Fp6_add_mpz(struct Fp6 *ANS,struct Fp6 *A,mpz_t a){
    Fp2_add_mpz(&ANS->x0,&A->x0,a);
    Fp2_add_mpz(&ANS->x1,&A->x1,a);
    Fp2_add_mpz(&ANS->x2,&A->x2,a);
}
void Fp6_sub(struct Fp6 *ANS,struct Fp6 *A,struct Fp6 *B){
    Fp2_sub(&ANS->x0,&A->x0,&B->x0);
    Fp2_sub(&ANS->x1,&A->x1,&B->x1);
    Fp2_sub(&ANS->x2,&A->x2,&B->x2);
}
void Fp6_sub_ui(struct Fp6 *ANS,struct Fp6 *A,unsigned long int a){
    Fp2_sub_ui(&ANS->x0,&A->x0,a);
    Fp2_sub_ui(&ANS->x1,&A->x1,a);
    Fp2_sub_ui(&ANS->x2,&A->x2,a);
}
void Fp6_sub_mpz(struct Fp6 *ANS,struct Fp6 *A,mpz_t a){
    Fp2_sub_mpz(&ANS->x0,&A->x0,a);
    Fp2_sub_mpz(&ANS->x1,&A->x1,a);
    Fp2_sub_mpz(&ANS->x2,&A->x2,a);
}
/*-------------------------------inverse--------------------------------*/
void Fp6_inv(struct Fp6 *ANS,struct Fp6 *A){
    struct Fp6 frob1,frob2,buf1, buf2;
    Fp6_init(&frob1);
    Fp6_init(&frob2);
    Fp6_init(&buf1);
    Fp6_init(&buf2);
    
    Fp6_inv_map_1(&frob1,A);
    Fp6_inv_map_2(&frob2,A);
    Fp6_mul(&buf1,&frob1,&frob2);
    Fp6_mul(&buf2,&buf1,A);
    Fp2_inv(&buf2.x0,&buf2.x0);
    Fp2_mul(&ANS->x0,&buf1.x0,&buf2.x0);
    Fp2_mul(&ANS->x1,&buf1.x1,&buf2.x0);
    Fp2_mul(&ANS->x2,&buf1.x2,&buf2.x0);
    
    Fp6_clear(&frob1);
    Fp6_clear(&frob2);
    Fp6_clear(&buf1);
    Fp6_clear(&buf2);
}
void Fp6_inv_map_1(struct Fp6 *ANS,struct Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_mul_mpz(&ANS->x1,&A->x1,inv_CNR1.x0);
    Fp2_mul_mpz(&ANS->x2,&A->x2,inv_CNR2.x0);
}
void Fp6_inv_map_2(struct Fp6 *ANS,struct Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_mul_mpz(&ANS->x1,&A->x1,inv_CNR2.x0);
    Fp2_mul_mpz(&ANS->x2,&A->x2,inv_CNR1.x0);
}
/*-------------------------------legendre-------------------------------*/
int Fp6_legendre(struct Fp6 *A){
    mpz_t exp;		mpz_init(exp);
    struct Fp6 buf;
    Fp6_init(&buf);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&buf,A,exp);
    
    mpz_clear(exp);
    if(Fp6_cmp_one(&buf)==0){
        Fp6_clear(&buf);
        return 1;
    }else if(Fp6_cmp_zero(&buf)==0){
        Fp6_clear(&buf);
        return 0;
    }else{
        Fp6_clear(&buf);
        return -1;
    }
}
int Fp6_isCNR(struct Fp6 *A){
    struct Fp6 buf;
    Fp6_init(&buf);
    mpz_t exp;
    mpz_init(exp);
    
    mpz_pow_ui(exp,prime,6);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,3);
    Fp6_pow(&buf,A,exp);
    
    mpz_clear(exp);
    if(Fp6_cmp_one(&buf)==0){
        Fp6_clear(&buf);
        return 1;
    }else if(Fp6_cmp_zero(&buf)==0){
        Fp6_clear(&buf);
        return 0;
    }else{
        Fp6_clear(&buf);
        return -1;
    }
}
/*---------------------------------sqr----------------------------------*/
void Fp6_sqrt(struct Fp6 *ANS,struct Fp6 *A){
    struct Fp6 buf1,buf2;
    Fp6_init(&buf1);
    Fp6_init(&buf2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp6_frobenius_4(&buf1,A);
    Fp6_frobenius_2(&buf2,A);
    Fp6_mul(&buf1,&buf1,&buf2);
    Fp6_mul(&buf1,&buf1,A);
    Fp6_set_ui(&buf2,0);
    Fp2_sqrt(&buf2.x0,&buf1.x0);
    Fp2_inv(&buf2.x0,&buf2.x0);
    Fp2_set(&buf2.x0,&buf2.x0);
    mpz_pow_ui(exp,prime,8);
    mpz_pow_ui(buf,prime,4);
    mpz_add(exp,exp,buf);
    mpz_add_ui(exp,exp,2);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp6_pow(&buf1,A,exp);
    Fp6_mul(&buf1,&buf1,&buf2);
    Fp6_set(ANS,&buf1);
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp6_clear(&buf1);
    Fp6_clear(&buf2);
}
/*---------------------------------pow----------------------------------*/
void Fp6_pow(struct Fp6 *ANS,struct Fp6 *A,mpz_t a){
    int i,length;
    length=(int)mpz_sizeinbase(a,2);
    char binary[length];
    mpz_get_str(binary,2,a);
    struct Fp6 buf;
    Fp6_init(&buf);
    Fp6_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp6_squaring(&buf,&buf);
        if(binary[i]=='1'){
            Fp6_mul(&buf,A,&buf);
        }
    }
    
    Fp6_set(ANS,&buf);
    Fp6_clear(&buf);
}
/*---------------------------------cmp----------------------------------*/
int Fp6_cmp(struct Fp6 *A,struct Fp6 *B){
    if(Fp2_cmp(&A->x0,&B->x0)==0 && Fp2_cmp(&A->x1,&B->x1)==0 && Fp2_cmp(&A->x2,&B->x2)==0){
        return 0;
    }
    return 1;
}
int Fp6_cmp_ui(struct Fp6 *A,unsigned long int a){
    if(Fp2_cmp_ui(&A->x0,a)==0 && Fp2_cmp_ui(&A->x1,a)==0 && Fp2_cmp_ui(&A->x2,a)==0){
        return 0;
    }
    return 1;
}
int Fp6_cmp_mpz(struct Fp6 *A,mpz_t a){
    if(Fp2_cmp_mpz(&A->x0,a)==0 && Fp2_cmp_mpz(&A->x1,a)==0 && Fp2_cmp_mpz(&A->x2,a)==0){
        return 0;
    }
    return 1;
}
int Fp6_cmp_zero(struct Fp6 *A){
    if(Fp2_cmp_zero(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}
int Fp6_cmp_one(struct Fp6 *A){
    if(Fp2_cmp_one(&A->x0)==0 && Fp2_cmp_zero(&A->x1)==0 && Fp2_cmp_zero(&A->x2)==0){
        return 0;
    }
    return 1;
}
/*--------------------------------frobenius--------------------------------*/
void Fp6_frobenius_1(struct Fp6 *ANS,struct Fp6 *A){
    struct Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0,&A->x0.x0);
    Fp_set_neg(&ANS->x0.x1,&A->x0.x1);
    //x1
    //Fp_set(&ANS->x1.x0,&A->x1.x0);
    //Fp_set_neg(&ANS->x1.x1,&A->x1.x1);
    //Fp2_mul(&ANS->x1,&ANS->x1,&Fp2_basis_prime_1_div_3_1);	//can be efficient
    Fp_set(&tmp,&A->x1.x0);
    Fp_set(&ANS->x1.x0,&A->x1.x1);
    Fp_set(&ANS->x1.x1,&tmp);
    Fp2_mul_mpz(&ANS->x1,&ANS->x1,Fp2_basis_prime_1_div_3_1.x1.x0);
    //x2
    Fp_set(&ANS->x2.x0,&A->x2.x0);
    Fp_set_neg(&ANS->x2.x1,&A->x2.x1);
    Fp2_mul_mpz(&ANS->x2,&ANS->x2,Fp2_basis_prime_1_div_3_2.x0.x0);
    
    Fp_clear(&tmp);
}
void Fp6_frobenius_2(struct Fp6 *ANS,struct Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_mul_mpz(&ANS->x1,&A->x1,Fp2_basis_prime_2_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2,&A->x2,Fp2_basis_prime_2_div_3_2.x0.x0);
}
void Fp6_frobenius_3(struct Fp6 *ANS,struct Fp6 *A){
    struct Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0,&A->x0.x0);
    Fp_set_neg(&ANS->x0.x1,&A->x0.x1);
    //x1
    Fp_set(&tmp,&A->x1.x0);
    Fp_set(&ANS->x1.x0,&A->x1.x1);
    Fp_set(&ANS->x1.x1,&tmp);
    //x2
    Fp_set_neg(&ANS->x2.x0,&A->x2.x0);
    Fp_set(&ANS->x2.x1,&A->x2.x1);
    
    Fp_clear(&tmp);
}
void Fp6_frobenius_4(struct Fp6 *ANS,struct Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_mul_mpz(&ANS->x1,&A->x1,Fp2_basis_prime_4_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2,&A->x2,Fp2_basis_prime_4_div_3_2.x0.x0);
}
void Fp6_frobenius_6(struct Fp6 *ANS,struct Fp6 *A){
    Fp6_set(ANS,A);
}
void Fp6_frobenius_8(struct Fp6 *ANS,struct Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_mul_mpz(&ANS->x1,&A->x1,Fp2_basis_prime_8_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2,&A->x2,Fp2_basis_prime_8_div_3_2.x0.x0);
}
void Fp6_frobenius_10(struct Fp6 *ANS,struct Fp6 *A){
    Fp2_set(&ANS->x0,&A->x0);
    Fp2_mul_mpz(&ANS->x1,&A->x1,Fp2_basis_prime_10_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x2,&A->x2,Fp2_basis_prime_10_div_3_2.x0.x0);
}
/*============================================================================*/
/* Fp12                                                                       */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void Fp12_init(struct Fp12 *P){
    Fp6_init(&P->x0);
    Fp6_init(&P->x1);
}
/*---------------------------------set----------------------------------*/
void Fp12_set(struct Fp12 *P,struct Fp12 *A){
    Fp6_set(&P->x0,&A->x0);
    Fp6_set(&P->x1,&A->x1);
}
void Fp12_set_ui(struct Fp12 *P,unsigned long int a){
    Fp6_set_ui(&P->x0,a);
    Fp6_set_ui(&P->x1,a);
}
void Fp12_set_mpz(struct Fp12 *P,mpz_t a){
    Fp6_set_mpz(&P->x0,a);
    Fp6_set_mpz(&P->x1,a);
}
void Fp12_set_neg(struct Fp12 *P,struct Fp12 *A){
    Fp6_set_neg(&P->x0,&A->x0);
    Fp6_set_neg(&P->x1,&A->x1);
}
/*---------------------------------random--------------------------------*/
void Fp12_random(struct Fp12 *P,gmp_randstate_t state){
    Fp6_random(&P->x0,state);
    Fp6_random(&P->x1,state);
}
/*---------------------------------clear---------------------------------*/
void Fp12_clear(struct Fp12 *P){
    Fp6_clear(&P->x0);
    Fp6_clear(&P->x1);
}
/*---------------------------------print---------------------------------*/
void Fp12_printf(struct Fp12 *P,char *name){
    printf("%s",name);
    printf("(");
    Fp6_printf(&P->x0,"");
    printf(",");
    Fp6_printf(&P->x1,"");
    printf(")");
}
/*---------------------------vector calculation--------------------------*/
void Fp12_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
    struct Fp6 tmp1,tmp2;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    
    //set
    Fp6_mul(&tmp2,&A->x1,&B->x1);//b*d
    Fp6_add(&tmp1,&A->x0,&A->x1);//a+b
    Fp6_add(&ANS->x1,&B->x0,&B->x1);//c+d
    Fp6_mul(&ANS->x1,&tmp1,&ANS->x1);//(a+b)(c+d)
    Fp6_mul(&tmp1,&A->x0,&B->x0);//a*c
    //x0
    Fp6_mul_basis(&ANS->x0,&tmp2);//b*d*v
    Fp6_add(&ANS->x0,&ANS->x0,&tmp1);//a*c+b*d*v
    //x1
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp1);
    Fp6_sub(&ANS->x1,&ANS->x1,&tmp2);
    
    //clear
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
}
void Fp12_mul_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a){
    Fp6_mul_ui(&ANS->x0,&A->x0,a);
    Fp6_mul_ui(&ANS->x1,&A->x1,a);
}
void Fp12_mul_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a){
    Fp6_mul_mpz(&ANS->x0,&A->x0,a);
    Fp6_mul_mpz(&ANS->x1,&A->x1,a);
}
void Fp12_squaring(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp6 tmp1,tmp2,tmp3;
    Fp6_init(&tmp1);
    Fp6_init(&tmp2);
    Fp6_init(&tmp3);
    
    Fp6_add(&tmp1,&A->x0,&A->x1);
    Fp6_mul_basis(&tmp2,&A->x1);
    Fp6_add(&tmp2,&tmp2,&A->x0);
    Fp6_mul(&tmp3,&A->x0,&A->x1);
    
    //x0
    Fp6_mul(&ANS->x0,&tmp1,&tmp2);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp3);
    Fp6_mul_basis(&tmp1,&tmp3);
    Fp6_sub(&ANS->x0,&ANS->x0,&tmp1);
    //x1
    Fp6_add(&ANS->x1,&tmp3,&tmp3);
    
    Fp6_clear(&tmp1);
    Fp6_clear(&tmp2);
    Fp6_clear(&tmp3);
}
void Fp12_add(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
    Fp6_add(&ANS->x0,&A->x0,&B->x0);
    Fp6_add(&ANS->x1,&A->x1,&B->x1);
}
void Fp12_add_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a){
    Fp6_add_ui(&ANS->x0,&A->x0,a);
    Fp6_add_ui(&ANS->x1,&A->x1,a);
}
void Fp12_add_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a){
    Fp6_add_mpz(&ANS->x0,&ANS->x0,a);
    Fp6_add_mpz(&ANS->x1,&ANS->x1,a);
}
void Fp12_sub(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
    Fp6_sub(&ANS->x0,&A->x0,&B->x0);
    Fp6_sub(&ANS->x1,&A->x1,&B->x1);
}
void Fp12_sub_ui(struct Fp12 *ANS,struct Fp12 *A,unsigned long int a){
    Fp6_sub_ui(&ANS->x0,&ANS->x0,a);
    Fp6_sub_ui(&ANS->x1,&ANS->x1,a);
}
void Fp12_sub_mpz(struct Fp12 *ANS,struct Fp12 *A,mpz_t a){
    Fp6_sub_mpz(&ANS->x0,&ANS->x0,a);
    Fp6_sub_mpz(&ANS->x1,&ANS->x1,a);
}
/*------------------------------inverse---------------------------------*/
void Fp12_inv(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp12 frob,buf;
    Fp12_init(&frob);
    Fp12_init(&buf);
    
    Fp12_inv_map(&frob,A);
    Fp12_mul(&buf,A,&frob);
    Fp6_inv(&buf.x0,&buf.x0);
    Fp6_mul(&ANS->x0,&frob.x0,&buf.x0);
    Fp6_mul(&ANS->x1,&frob.x1,&buf.x0);
    
    Fp12_clear(&frob);
    Fp12_clear(&buf);
}
void Fp12_inv_map(struct Fp12 *ANS,struct Fp12 *A){
    Fp6_set(&ANS->x0,&A->x0);
    Fp6_set_neg(&ANS->x1,&A->x1);
}
/*-------------------------------legendre-------------------------------*/
int Fp12_legendre(struct Fp12 *A){
    mpz_t exp;
    mpz_init(exp);
    struct Fp12 buf;
    Fp12_init(&buf);
    
    mpz_pow_ui(exp,prime,12);
    mpz_sub_ui(exp,exp,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&buf,A,exp);
    
    mpz_clear(exp);
    if(Fp12_cmp_one(&buf)==0){
        Fp12_clear(&buf);
        return 1;
    }else if(Fp12_cmp_zero(&buf)==0){
        Fp12_clear(&buf);
        return 0;
    }else{
        Fp12_clear(&buf);
        return -1;
    }
}
/*---------------------------------sqr----------------------------------*/
void Fp12_sqrt(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp12 x,y,t,k,n,buf;
    Fp12_init(&x);
    Fp12_init(&y);
    Fp12_init(&t);
    Fp12_init(&k);
    Fp12_init(&n);
    Fp12_init(&buf);
    unsigned long int e,m;
    mpz_t exp,q,z,result;
    mpz_init(exp);
    mpz_init(q);
    mpz_init(z);
    mpz_init(result);
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    Fp12_random(&n,state);
    while(Fp12_legendre(&n)!=-1){
        Fp12_random(&n,state);
    }
    mpz_pow_ui(q,prime,12);
    mpz_sub_ui(q,q,1);
    mpz_mod_ui(result,q,2);
    e=0;
    while(mpz_cmp_ui(result,0)==0){
        mpz_tdiv_q_ui(q,q,2);
        mpz_mod_ui(result,q,2);
        e++;
    }
    Fp12_pow(&y,&n,q);
    mpz_set_ui(z,e);
    mpz_sub_ui(exp,q,1);
    mpz_tdiv_q_ui(exp,exp,2);
    Fp12_pow(&x,A,exp);
    Fp12_mul(&buf,&x,&x);
    Fp12_mul(&k,&buf,A);
    Fp12_mul(&x,&x,A);
    while(Fp12_cmp_one(&k)!=0){
        m=1;
        mpz_ui_pow_ui(exp,2,m);
        Fp12_pow(&buf,&k,exp);
        while(Fp12_cmp_one(&buf)!=0){
            m++;
            mpz_ui_pow_ui(exp,2,m);
            Fp12_pow(&buf,&k,exp);
        }
        mpz_sub_ui(exp,z,m);
        mpz_sub_ui(exp,exp,1);
        mpz_ui_pow_ui(result,2,mpz_get_ui(exp));
        Fp12_pow(&t,&y,result);
        Fp12_mul(&y,&t,&t);
        mpz_set_ui(z,m);
        Fp12_mul(&x,&x,&t);
        Fp12_mul(&k,&k,&y);
    }
    Fp12_set(ANS,&x);
    
    mpz_clear(q);
    mpz_clear(z);
    mpz_clear(exp);
    mpz_clear(result);
    Fp12_clear(&x);
    Fp12_clear(&y);
    Fp12_clear(&t);
    Fp12_clear(&k);
    Fp12_clear(&n);
    Fp12_clear(&buf);
}
/*---------------------------------pow----------------------------------*/
void Fp12_pow(struct Fp12 *ANS,struct Fp12 *A,mpz_t a){
    int i,length;
    length=(int)mpz_sizeinbase(a,2);
    char binary[length];
    mpz_get_str(binary,2,a);
    struct Fp12 buf;
    Fp12_init(&buf);
    Fp12_set(&buf,A);
    
    for(i=1; binary[i]!='\0'; i++){
        Fp12_squaring(&buf,&buf);
        if(binary[i]=='1'){
            Fp12_mul(&buf,A,&buf);
        }
    }
    
    Fp12_set(ANS,&buf);
    Fp12_clear(&buf);
}


/*---------------------------------cmp----------------------------------*/
int Fp12_cmp(struct Fp12 *A,struct Fp12 *B){
    if(Fp6_cmp(&A->x0,&B->x0)==0 && Fp6_cmp(&A->x1,&B->x1)==0){
        return 0;
    }
    return 1;
}
int Fp12_cmp_ui(struct Fp12 *A,unsigned long int a){
    if(Fp6_cmp_ui(&A->x0,a)==0 && Fp6_cmp_ui(&A->x1,a)==0){
        return 0;
    }
    return 1;
}
int Fp12_cmp_mpz(struct Fp12 *A,mpz_t a){
    if(Fp6_cmp_mpz(&A->x0,a)==0 && Fp6_cmp_mpz(&A->x1,a)==0){
        return 0;
    }
    return 1;
}
int Fp12_cmp_zero(struct Fp12 *A){
    if(Fp6_cmp_zero(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
int Fp12_cmp_one(struct Fp12 *A){
    if(Fp6_cmp_one(&A->x0)==0 && Fp6_cmp_zero(&A->x1)==0){
        return 0;
    }
    return 1;
}
/*-------------------------------frobenius--------------------------------*/
void Fp12_frobenius_1(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    //Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x0);
    //Fp_set_neg(&ANS->x0.x1.x1,&A->x0.x1.x1);
    //Fp2_mul(&ANS->x0.x1,&ANS->x0.x1,&Fp2_basis_prime_1_div_3_1);	//can be efficient
    Fp_set(&tmp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp);
    Fp2_mul_mpz(&ANS->x0.x1,&ANS->x0.x1,Fp2_basis_prime_1_div_3_1.x1.x0);
    
    Fp_set(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set_neg(&ANS->x0.x2.x1,&A->x0.x2.x1);
    Fp2_mul_mpz(&ANS->x0.x2,&ANS->x0.x2,Fp2_basis_prime_1_div_3_2.x0.x0);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&Fp2_basis_prime_1_div_6);
    //Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x0);
    //Fp_set_neg(&ANS->x1.x1.x1,&A->x1.x1.x1);
    //Fp2_mul(&ANS->x1.x1,&ANS->x1.x1,&Fp2_basis_prime_1_div_3_1);	//can be efficient
    Fp_set(&tmp,&A->x1.x1.x0);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,Fp2_basis_prime_1_div_3_1.x1.x0);
    
    Fp2_mul(&ANS->x1.x1,&ANS->x1.x1,&Fp2_basis_prime_1_div_6);
    Fp_set(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set_neg(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,Fp2_basis_prime_1_div_3_2.x0.x0);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&Fp2_basis_prime_1_div_6);
    
    Fp_clear(&tmp);
}
void Fp12_frobenius_2(struct Fp12 *ANS,struct Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,Fp2_basis_prime_2_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,Fp2_basis_prime_2_div_3_2.x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,Fp2_basis_prime_2_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1,&A->x1.x1,Fp2_basis_prime_2_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,Fp2_basis_prime_2_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,Fp2_basis_prime_2_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,Fp2_basis_prime_2_div_6.x0.x0);
}
void Fp12_frobenius_3(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp tmp;
    Fp_init(&tmp);
    
    //x0
    Fp_set(&ANS->x0.x0.x0,&A->x0.x0.x0);
    Fp_set_neg(&ANS->x0.x0.x1,&A->x0.x0.x1);
    Fp_set(&tmp,&A->x0.x1.x0);
    Fp_set(&ANS->x0.x1.x0,&A->x0.x1.x1);
    Fp_set(&ANS->x0.x1.x1,&tmp);
    Fp_set_neg(&ANS->x0.x2.x0,&A->x0.x2.x0);
    Fp_set(&ANS->x0.x2.x1,&A->x0.x2.x1);
    //x1
    Fp_set(&ANS->x1.x0.x0,&A->x1.x0.x0);
    Fp_set_neg(&ANS->x1.x0.x1,&A->x1.x0.x1);
    Fp2_mul(&ANS->x1.x0,&ANS->x1.x0,&Fp2_basis_prime_3_div_6);
    Fp_set(&tmp,&A->x1.x1.x0);
    Fp_set(&ANS->x1.x1.x0,&A->x1.x1.x1);
    Fp_set(&ANS->x1.x1.x1,&tmp);
    Fp2_mul(&ANS->x1.x1,&ANS->x1.x1,&Fp2_basis_prime_3_div_6);
    Fp_set_neg(&ANS->x1.x2.x0,&A->x1.x2.x0);
    Fp_set(&ANS->x1.x2.x1,&A->x1.x2.x1);
    Fp2_mul(&ANS->x1.x2,&ANS->x1.x2,&Fp2_basis_prime_3_div_6);
    
    Fp_clear(&tmp);
}
void Fp12_frobenius_4(struct Fp12 *ANS,struct Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,Fp2_basis_prime_4_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,Fp2_basis_prime_4_div_3_2.x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,Fp2_basis_prime_4_div_6.x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,Fp2_basis_prime_4_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,Fp2_basis_prime_4_div_6.x0.x0);
}
void Fp12_frobenius_6(struct Fp12 *ANS,struct Fp12 *A){
    //x0
    Fp6_set(&ANS->x0,&A->x0);
    //x1
    Fp6_set_neg(&ANS->x1,&A->x1);
}
void Fp12_frobenius_8(struct Fp12 *ANS,struct Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,Fp2_basis_prime_8_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,Fp2_basis_prime_8_div_3_2.x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,Fp2_basis_prime_8_div_6.x0.x0);
    Fp2_set(&ANS->x1.x1,&A->x1.x1);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,Fp2_basis_prime_8_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,Fp2_basis_prime_8_div_6.x0.x0);
}
void Fp12_frobenius_10(struct Fp12 *ANS,struct Fp12 *A){
    //x0
    Fp2_set(&ANS->x0.x0,&A->x0.x0);
    Fp2_mul_mpz(&ANS->x0.x1,&A->x0.x1,Fp2_basis_prime_10_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x0.x2,&A->x0.x2,Fp2_basis_prime_10_div_3_2.x0.x0);
    //x1
    Fp2_mul_mpz(&ANS->x1.x0,&A->x1.x0,Fp2_basis_prime_10_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1,&A->x1.x1,Fp2_basis_prime_10_div_3_1.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x1,&ANS->x1.x1,Fp2_basis_prime_10_div_6.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2,&A->x1.x2,Fp2_basis_prime_10_div_3_2.x0.x0);
    Fp2_mul_mpz(&ANS->x1.x2,&ANS->x1.x2,Fp2_basis_prime_10_div_6.x0.x0);
}
/*============================================================================*/
/* EFp                                                                        */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void EFp_init(struct EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->flag=0;
}
/*---------------------------------set----------------------------------*/
void EFp_set(struct EFp *P,struct EFp *A){
    Fp_set(&P->x,&A->x);
    Fp_set(&P->y,&A->y);
    P->flag=A->flag;
}
void EFp_set_ui(struct EFp *P,unsigned long int a){
    Fp_set_ui(&P->x,a);
    Fp_set_ui(&P->y,a);
    P->flag=0;
}
void EFp_set_mpz(struct EFp *P,mpz_t a){
    Fp_set_mpz(&P->x,a);
    Fp_set_mpz(&P->y,a);
    P->flag=0;
}
void EFp_set_neg(struct EFp *P,struct EFp *A){
    Fp_set(&P->x,&A->x);
    Fp_set_neg(&P->y,&A->y);
    P->flag=A->flag;
}
/*---------------------------------clear--------------------------------*/
void EFp_clear(struct EFp *P){
    Fp_clear(&P->x);
    Fp_clear(&P->y);
}
/*---------------------------------print--------------------------------*/
void EFp_printf(struct EFp *P,char *name){
    printf("%s",name);
    if(P->flag==0){
        printf("(");
        Fp_printf(&P->x,"");
        printf(",");
        Fp_printf(&P->y,"");
        printf(")");
    }else{
        printf("0");
    }
}
/*-----------------------------rational point---------------------------*/
void EFp_rational_point(struct EFp *P){
    struct Fp buf1,buf2,R;
    Fp_init(&buf1);
    Fp_init(&buf2);
    Fp_init(&R);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp_random(&P->x,state);
        Fp_mul(&buf1,&P->x,&P->x);
        Fp_mul(&buf2,&buf1,&P->x);
        Fp_mul_mpz(&buf1,&P->x,curve_parameter_A);
        Fp_add(&R,&buf1,&buf2);
        Fp_add_mpz(&R,&R,curve_parameter_B);
        if(Fp_legendre(&R)==1){
            Fp_sqrt(&P->y,&R);
            break;
        }
    }
    
    Fp_clear(&buf1);
    Fp_clear(&buf2);
    Fp_clear(&R);
}
/*---------------------------------SCM----------------------------------*/
void EFp_ECD(struct EFp *ANS,struct EFp *P){
    if(Fp_cmp_zero(&P->y)==0){
        ANS->flag=1;
        return;
    }
    
    struct EFp Tmp;
    EFp_init(&Tmp);
    EFp_set(&Tmp,P);
    struct Fp Buf1,Buf2,C;
    Fp_init(&Buf1);
    Fp_init(&Buf2);
    Fp_init(&C);
    
    Fp_mul_ui(&Buf1,&Tmp.y,2);
    Fp_inv(&Buf1,&Buf1);
    Fp_mul(&Buf2,&Tmp.x,&Tmp.x);
    Fp_mul_ui(&Buf2,&Buf2,3);
    Fp_add_mpz(&Buf2,&Buf2,curve_parameter_A);
    Fp_mul(&C,&Buf1,&Buf2);
    Fp_mul(&Buf1,&C,&C);
    Fp_mul_ui(&Buf2,&Tmp.x,2);
    Fp_sub(&ANS->x,&Buf1,&Buf2);
    Fp_sub(&Buf1,&Tmp.x,&ANS->x);
    Fp_mul(&Buf2,&C,&Buf1);
    Fp_sub(&ANS->y,&Buf2,&Tmp.y);
    
    //clear
    Fp_clear(&Buf1);
    Fp_clear(&Buf2);
    Fp_clear(&C);
    EFp_clear(&Tmp);
}
void EFp_ECA(struct EFp *ANS,struct EFp *P1,struct EFp *P2){
    if(P1->flag==1){
        EFp_set(ANS,P2);
        return;
    }else if(P2->flag==1){
        EFp_set(ANS,P1);
        return;
    }else if(Fp_cmp(&P1->x,&P2->x)==0){
        if(Fp_cmp(&P1->y,&P2->y)!=0){
            ANS->flag=1;
            return;
        }else{
            EFp_ECD(ANS,P1);
            return;
        }
    }
    
    struct EFp Tmp1,Tmp2;
    EFp_init(&Tmp1);
    EFp_set(&Tmp1,P1);
    EFp_init(&Tmp2);
    EFp_set(&Tmp2,P2);
    struct Fp Buf1,Buf2,C;
    Fp_init(&Buf1);
    Fp_init(&Buf2);
    Fp_init(&C);
    
    Fp_sub(&Buf1,&Tmp2.x,&Tmp1.x);
    Fp_inv(&Buf1,&Buf1);
    Fp_sub(&Buf2,&Tmp2.y,&Tmp1.y);
    Fp_mul(&C,&Buf1,&Buf2);
    Fp_mul(&Buf1,&C,&C);
    Fp_sub(&Buf2,&Buf1,&Tmp1.x);
    Fp_sub(&ANS->x,&Buf2,&Tmp2.x);
    Fp_sub(&Buf1,&Tmp1.x,&ANS->x);
    Fp_mul(&Buf2,&C,&Buf1);
    Fp_sub(&ANS->y,&Buf2,&Tmp1.y);
    
    //clear
    Fp_clear(&Buf1);
    Fp_clear(&Buf2);
    Fp_clear(&C);
    EFp_clear(&Tmp1);
    EFp_clear(&Tmp2);
}
void EFp_SCM(struct EFp *ANS,struct EFp *P,mpz_t R){
    if(mpz_cmp_ui(R,0)==0){
        ANS->flag=1;
        return;
    }else if(mpz_cmp_ui(R,1)==0){
        EFp_set(ANS,P);
        return;
    }
    
    struct EFp Tmp,next_P;
    EFp_init(&Tmp);
    EFp_set(&Tmp,P);
    EFp_init(&next_P);
    int i,length;
    length=(int)mpz_sizeinbase(R,2);
    char binary[length];
    mpz_get_str(binary,2,R);
    
    EFp_set(&next_P,&Tmp);
    for(i=1; binary[i]!='\0'; i++){
        EFp_ECD(&next_P,&next_P);
        if(binary[i]=='1'){
            EFp_ECA(&next_P,&next_P,&Tmp);
        }
    }
    
    EFp_set(ANS,&next_P);
    
    EFp_clear(&next_P);
    EFp_clear(&Tmp);
}
/*----------------------------skew frobenius----------------------------*/
void EFp_skew_frobenius(struct EFp *ANS,struct EFp *A){
    Fp_mul(&ANS->x,&A->x,&inv_CNR1);
    Fp_set_neg(&ANS->y,&A->y);
}

/*============================================================================*/
/* EFp2                                                                       */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void EFp2_init(struct EFp2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    P->flag=0;
}
/*---------------------------------set----------------------------------*/
void EFp2_set(struct EFp2 *P,struct EFp2 *A){
    Fp2_set(&P->x,&A->x);
    Fp2_set(&P->y,&A->y);
    P->flag=A->flag;
}
void EFp2_set_ui(struct EFp2 *P,unsigned long int a){
    Fp2_set_ui(&P->x,a);
    Fp2_set_ui(&P->y,a);
    P->flag=0;
}
void EFp2_set_mpz(struct EFp2 *P,mpz_t a){
    Fp2_set_mpz(&P->x,a);
    Fp2_set_mpz(&P->y,a);
    P->flag=0;
}
void EFp2_set_neg(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_set(&ANS->x,&P->x);
    Fp2_set_neg(&ANS->y,&P->y);
    ANS->flag=P->flag;
}
/*---------------------------------clear--------------------------------*/
void EFp2_clear(struct EFp2 *P){
    Fp2_clear(&P->x);
    Fp2_clear(&P->y);
}
/*---------------------------------print--------------------------------*/
void EFp2_printf(struct EFp2 *P,char *name){
    printf("%s",name);
    if(P->flag==0){
        printf("(");
        Fp2_printf(&P->x,"X");
        printf(",");
        Fp2_printf(&P->y,"Y");
        printf(")");
    }else{
        printf("0");
    }
}
/*-----------------------------rational point---------------------------*/
void EFp2_rational_point(struct EFp2 *P){
    struct Fp2 buf1,buf2,R;
    Fp2_init(&buf1);
    Fp2_init(&buf2);
    Fp2_init(&R);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp2_random(&P->x,state);
        Fp2_mul(&buf1,&P->x,&P->x);
        Fp2_mul(&buf2,&buf1,&P->x);
        Fp2_mul_mpz(&buf1,&P->x,curve_parameter_A);
        Fp2_add(&R,&buf1,&buf2);
        mpz_add(R.x0.x0,R.x0.x0,curve_parameter_B);
        if(Fp2_legendre(&R)==1){
            Fp2_sqrt(&P->y,&R);
            break;
        }
    }
    
    Fp2_clear(&buf1);
    Fp2_clear(&buf2);
    Fp2_clear(&R);
}
/*---------------------------------SCM----------------------------------*/
void EFp2_ECD(struct EFp2 *ANS,struct EFp2 *P){
    if(Fp2_cmp_zero(&P->y)==0){
        ANS->flag=1;
        return;
    }
    
    struct EFp2 Tmp;
    EFp2_init(&Tmp);
    EFp2_set(&Tmp,P);
    struct Fp2 Buf1,Buf2,C;
    Fp2_init(&Buf1);
    Fp2_init(&Buf2);
    Fp2_init(&C);
    
    Fp2_mul_ui(&Buf1,&Tmp.y,2);
    
    Fp2_inv(&Buf1,&Buf1);
    Fp2_mul(&Buf2,&Tmp.x,&Tmp.x);
    Fp2_mul_ui(&Buf2,&Buf2,3);
    mpz_add(Buf2.x0.x0,Buf2.x0.x0,curve_parameter_A);
    Fp2_mul(&C,&Buf1,&Buf2);
    
    Fp2_squaring(&Buf1,&C);
    Fp2_mul_ui(&Buf2,&Tmp.x,2);
    Fp2_sub(&ANS->x,&Buf1,&Buf2);
    
    Fp2_sub(&Buf1,&Tmp.x,&ANS->x);
    Fp2_mul(&Buf2,&C,&Buf1);
    Fp2_sub(&ANS->y,&Buf2,&Tmp.y);
    
    Fp2_clear(&Buf1);
    Fp2_clear(&Buf2);
    Fp2_clear(&C);
    EFp2_clear(&Tmp);
}
void EFp2_ECA(struct EFp2 *ANS,struct EFp2 *P1,struct EFp2 *P2){
    if(P1->flag==1){
        EFp2_set(ANS,P2);
        return;
    }else if(P2->flag==1){
        EFp2_set(ANS,P1);
        return;
    }else if(Fp2_cmp(&P1->x,&P2->x)==0){
        if(Fp2_cmp(&P1->y,&P2->y)!=0){
            ANS->flag=1;
            return;
        }else{
            EFp2_ECD(ANS,P1);
            return;
        }
    }
    
    struct EFp2 Tmp1,Tmp2;
    EFp2_init(&Tmp1);
    EFp2_set(&Tmp1,P1);
    EFp2_init(&Tmp2);
    EFp2_set(&Tmp2,P2);
    struct Fp2 Buf1,Buf2,C;
    Fp2_init(&Buf1);
    Fp2_init(&Buf2);
    Fp2_init(&C);
    
    Fp2_sub(&Buf1,&Tmp2.x,&Tmp1.x);
    Fp2_inv(&Buf1,&Buf1);
    Fp2_sub(&Buf2,&Tmp2.y,&Tmp1.y);
    Fp2_mul(&C,&Buf1,&Buf2);
    Fp2_squaring(&Buf1,&C);
    Fp2_sub(&Buf2,&Buf1,&Tmp1.x);
    Fp2_sub(&ANS->x,&Buf2,&Tmp2.x);
    Fp2_sub(&Buf1,&Tmp1.x,&ANS->x);
    Fp2_mul(&Buf2,&C,&Buf1);
    Fp2_sub(&ANS->y,&Buf2,&Tmp1.y);
    
    //clear
    Fp2_clear(&Buf1);
    Fp2_clear(&Buf2);
    Fp2_clear(&C);
    EFp2_clear(&Tmp1);
    EFp2_clear(&Tmp2);
}
void EFp2_SCM(struct EFp2 *ANS,struct EFp2 *P,mpz_t R){
    if(mpz_cmp_ui(R,0)==0){
        ANS->flag=1;
        return;
    }else if(mpz_cmp_ui(R,1)==0){
        EFp2_set(ANS,P);
        return;
    }
    
    struct EFp2 Tmp,next_P;
    EFp2_init(&Tmp);
    EFp2_set(&Tmp,P);
    EFp2_init(&next_P);
    int i,length;
    length=(int)mpz_sizeinbase(R,2);
    char binary[length];
    mpz_get_str(binary,2,R);
    
    EFp2_set(&next_P,&Tmp);
    for(i=1; binary[i]!='\0'; i++){
        EFp2_ECD(&next_P,&next_P);
        if(binary[i]=='1'){
            EFp2_ECA(&next_P,&next_P,&Tmp);
        }
    }
    EFp2_set(ANS,&next_P);
    
    EFp2_clear(&next_P);
    EFp2_clear(&Tmp);
}
/*-----------------------------frobenius--------------------------------*/
void EFp2_frobenius_1(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_frobenius_1(&ANS->x,&P->x);
    Fp2_frobenius_1(&ANS->y,&P->y);
}
void EFp2_frobenius_2(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_frobenius_2(&ANS->x,&P->x);
    Fp2_frobenius_2(&ANS->y,&P->y);
}
void EFp2_frobenius_3(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_frobenius_3(&ANS->x,&P->x);
    Fp2_frobenius_3(&ANS->y,&P->y);
}
void EFp2_frobenius_4(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_frobenius_4(&ANS->x,&P->x);
    Fp2_frobenius_4(&ANS->y,&P->y);
}
void EFp2_frobenius_6(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_frobenius_6(&ANS->x,&P->x);
    Fp2_frobenius_6(&ANS->y,&P->y);
}
void EFp2_frobenius_8(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_frobenius_8(&ANS->x,&P->x);
    Fp2_frobenius_8(&ANS->y,&P->y);
}
void EFp2_frobenius_10(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_frobenius_10(&ANS->x,&P->x);
    Fp2_frobenius_10(&ANS->y,&P->y);
}
/*-----------------------------skew frobenius--------------------------------*/
void EFp2_skew_frobenius_1(struct EFp2 *ANS,struct EFp2 *A){
    //x
    Fp2_frobenius_1(&ANS->x,&A->x);
    Fp2_mul(&ANS->x,&ANS->x,&Fp2_basis_inv_prime_1_div_3);
    //y
    Fp2_frobenius_1(&ANS->y,&A->y);
    Fp2_mul(&ANS->y,&ANS->y,&Fp2_basis_inv_prime_1_div_2);
}
void EFp2_skew_frobenius_2(struct EFp2 *ANS,struct EFp2 *A){
    //x
    Fp2_frobenius_2(&ANS->x,&A->x);
    Fp2_mul(&ANS->x,&ANS->x,&Fp2_basis_inv_prime_2_div_3);
    //y
    Fp2_frobenius_2(&ANS->y,&A->y);
    Fp2_mul(&ANS->y,&ANS->y,&Fp2_basis_inv_prime_2_div_2);
}
void EFp2_skew_frobenius_3(struct EFp2 *ANS,struct EFp2 *A){
    //x
    Fp2_frobenius_3(&ANS->x,&A->x);
    Fp2_mul(&ANS->x,&ANS->x,&Fp2_basis_inv_prime_3_div_3);
    //y
    Fp2_frobenius_3(&ANS->y,&A->y);
    Fp2_mul(&ANS->y,&ANS->y,&Fp2_basis_inv_prime_3_div_2);
}
void EFp2_skew_frobenius_4(struct EFp2 *ANS,struct EFp2 *A){
    //x
    Fp2_mul_mpz(&ANS->x,&A->x,Fp2_basis_prime_4_div_3_2.x0.x0);
    //y
    Fp2_set(&ANS->y,&A->y);
}
void EFp2_skew_frobenius_10(struct EFp2 *ANS,struct EFp2 *A){
    //x
    Fp2_frobenius_10(&ANS->x,&A->x);
    Fp2_mul(&ANS->x,&ANS->x,&Fp2_basis_inv_prime_10_div_3);
    //y
    Fp2_frobenius_10(&ANS->y,&A->y);
    Fp2_mul(&ANS->y,&ANS->y,&Fp2_basis_inv_prime_10_div_2);
}
/*============================================================================*/
/* EFp6                                                                       */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void EFp6_init(struct EFp6 *P){
    Fp6_init(&P->x);
    Fp6_init(&P->y);
    P->flag=0;
}
void EFp6_set(struct EFp6 *P,struct EFp6 *A){
    Fp6_set(&P->x,&A->x);
    Fp6_set(&P->y,&A->y);
    P->flag=A->flag;
}
void EFp6_set_ui(struct EFp6 *P,unsigned long int a){
    Fp6_set_ui(&P->x,a);
    Fp6_set_ui(&P->y,a);
    P->flag=0;
}
void EFp6_set_mpz(struct EFp6 *P,mpz_t a){
    Fp6_set_mpz(&P->x,a);
    Fp6_set_mpz(&P->y,a);
    P->flag=0;
}
void EFp6_set_neg(struct EFp6 *P,struct EFp6 *A){
    Fp6_set(&P->x,&A->x);
    Fp6_set_neg(&P->y,&A->y);
    P->flag=A->flag;
}
/*---------------------------------clear--------------------------------*/
void EFp6_clear(struct EFp6 *P){
    Fp6_clear(&P->x);
    Fp6_clear(&P->y);
}
/*---------------------------------print--------------------------------*/
void EFp6_printf(struct EFp6 *P,char *name){
    printf("%s",name);
    if(P->flag==0){
        printf("(");
        Fp6_printf(&P->x,"X");
        printf(",");
        Fp6_printf(&P->y,"Y");
        printf(")");
    }else{
        printf("0");
    }
}
/*-----------------------------rational point---------------------------*/
void EFp6_rational_point(struct EFp6 *P){
    struct Fp6 buf1,buf2,R;
    Fp6_init(&buf1);
    Fp6_init(&buf2);
    Fp6_init(&R);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp6_random(&P->x,state);
        Fp6_mul(&buf1,&P->x,&P->x);
        Fp6_mul(&buf2,&buf1,&P->x);
        Fp6_mul_mpz(&buf1,&P->x,curve_parameter_A);
        Fp6_add(&R,&buf1,&buf2);
        mpz_add(R.x0.x0.x0,R.x0.x0.x0,curve_parameter_B);
        if(Fp6_legendre(&R)==1){
            Fp6_sqrt(&P->y,&R);
            break;
        }
    }
    
    Fp6_clear(&buf1);
    Fp6_clear(&buf2);
    Fp6_clear(&R);
}
/*---------------------------------SCM----------------------------------*/
void EFp6_ECD(struct EFp6 *ANS,struct EFp6 *P){
    if(Fp6_cmp_zero(&P->y)==0){
        ANS->flag=1;
        return;
    }
    
    struct EFp6 Tmp;
    EFp6_init(&Tmp);
    EFp6_set(&Tmp,P);
    struct Fp6 Buf1,Buf2,C;
    Fp6_init(&Buf1);
    Fp6_init(&Buf2);
    Fp6_init(&C);
    
    Fp6_mul_ui(&Buf1,&Tmp.y,2);
    
    Fp6_inv(&Buf1,&Buf1);
    Fp6_mul(&Buf2,&Tmp.x,&Tmp.x);
    Fp6_mul_ui(&Buf2,&Buf2,3);
    mpz_add(Buf2.x0.x0.x0,Buf2.x0.x0.x0,curve_parameter_A);
    Fp6_mul(&C,&Buf1,&Buf2);
    Fp6_mul(&Buf1,&C,&C);
    Fp6_mul_ui(&Buf2,&Tmp.x,2);
    Fp6_sub(&ANS->x,&Buf1,&Buf2);
    Fp6_sub(&Buf1,&Tmp.x,&ANS->x);
    Fp6_mul(&Buf2,&C,&Buf1);
    Fp6_sub(&ANS->y,&Buf2,&Tmp.y);
    
    Fp6_clear(&Buf1);
    Fp6_clear(&Buf2);
    Fp6_clear(&C);
    EFp6_clear(&Tmp);
}
void EFp6_ECA(struct EFp6 *ANS,struct EFp6 *P1,struct EFp6 *P2){
    if(P1->flag==1){
        EFp6_set(ANS,P2);
        return;
    }else if(P2->flag==1){
        EFp6_set(ANS,P1);
        return;
    }else if(Fp6_cmp(&P1->x,&P2->x)==0){
        if(Fp6_cmp(&P1->y,&P2->y)!=0){
            ANS->flag=1;
            return;
        }else{
            EFp6_ECD(ANS,P1);
            return;
        }
    }
    
    struct EFp6 Tmp1,Tmp2;
    EFp6_init(&Tmp1);
    EFp6_set(&Tmp1,P1);
    EFp6_init(&Tmp2);
    EFp6_set(&Tmp2,P2);
    struct Fp6 Buf1,Buf2,C;
    Fp6_init(&Buf1);
    Fp6_init(&Buf2);
    Fp6_init(&C);
    
    Fp6_sub(&Buf1,&Tmp2.x,&Tmp1.x);
    Fp6_inv(&Buf1,&Buf1);
    Fp6_sub(&Buf2,&Tmp2.y,&Tmp1.y);
    Fp6_mul(&C,&Buf1,&Buf2);
    Fp6_mul(&Buf1,&C,&C);
    Fp6_sub(&Buf2,&Buf1,&Tmp1.x);
    Fp6_sub(&ANS->x,&Buf2,&Tmp2.x);
    Fp6_sub(&Buf1,&Tmp1.x,&ANS->x);
    Fp6_mul(&Buf2,&C,&Buf1);
    Fp6_sub(&ANS->y,&Buf2,&Tmp1.y);
    
    //clear
    Fp6_clear(&Buf1);
    Fp6_clear(&Buf2);
    Fp6_clear(&C);
    EFp6_clear(&Tmp1);
    EFp6_clear(&Tmp2);
}
void EFp6_SCM(struct EFp6 *ANS,struct EFp6 *P,mpz_t R){
    if(mpz_cmp_ui(R,0)==0){
        ANS->flag=1;
        return;
    }else if(mpz_cmp_ui(R,1)==0){
        EFp6_set(ANS,P);
        return;
    }
    
    struct EFp6 Tmp,next_P;
    EFp6_init(&Tmp);
    EFp6_set(&Tmp,P);
    EFp6_init(&next_P);
    int i,length;
    length=(int)mpz_sizeinbase(R,2);
    char binary[length];
    mpz_get_str(binary,2,R);
    mpz_t order,buf;
    mpz_init(order);
    mpz_init(buf);
    
    EFp6_set(&next_P,&Tmp);
    for(i=1; binary[i]!='\0'; i++){
        EFp6_ECD(&next_P,&next_P);
        if(binary[i]=='1'){
            EFp6_ECA(&next_P,&next_P,&Tmp);
        }
    }
    
    EFp6_set(ANS,&next_P);
    
    EFp6_clear(&next_P);
    EFp6_clear(&Tmp);
}
/*-----------------------------frobenius--------------------------------*/
void EFp6_frobenius_1(struct EFp6 *ANS,struct EFp6 *P){
    Fp6_frobenius_1(&ANS->x,&P->x);
    Fp6_frobenius_1(&ANS->y,&P->y);
}
void EFp6_frobenius_2(struct EFp6 *ANS,struct EFp6 *P){
    Fp6_frobenius_2(&ANS->x,&P->x);
    Fp6_frobenius_2(&ANS->y,&P->y);
}
void EFp6_frobenius_3(struct EFp6 *ANS,struct EFp6 *P){
    Fp6_frobenius_3(&ANS->x,&P->x);
    Fp6_frobenius_3(&ANS->y,&P->y);
}
void EFp6_frobenius_4(struct EFp6 *ANS,struct EFp6 *P){
    Fp6_frobenius_4(&ANS->x,&P->x);
    Fp6_frobenius_4(&ANS->y,&P->y);
}
void EFp6_frobenius_6(struct EFp6 *ANS,struct EFp6 *P){
    Fp6_frobenius_6(&ANS->x,&P->x);
    Fp6_frobenius_6(&ANS->y,&P->y);
}
void EFp6_frobenius_8(struct EFp6 *ANS,struct EFp6 *P){
    Fp6_frobenius_8(&ANS->x,&P->x);
    Fp6_frobenius_8(&ANS->y,&P->y);
}
void EFp6_frobenius_10(struct EFp6 *ANS,struct EFp6 *P){
    Fp6_frobenius_10(&ANS->x,&P->x);
    Fp6_frobenius_10(&ANS->y,&P->y);
}
/*============================================================================*/
/* EFp12                                                                      */
/*============================================================================*/
/*---------------------------------init---------------------------------*/
void EFp12_init(struct EFp12 *P){
    Fp12_init(&P->x);
    Fp12_init(&P->y);
    P->flag=0;
}
/*---------------------------------set----------------------------------*/
void EFp12_set(struct EFp12 *P,struct EFp12 *A){
    Fp12_set(&P->x,&A->x);
    Fp12_set(&P->y,&A->y);
    P->flag=A->flag;
}
void EFp12_set_ui(struct EFp12 *P,unsigned long int a){
    Fp12_set_ui(&P->x,a);
    Fp12_set_ui(&P->y,a);
    P->flag=0;
}
void EFp12_set_mpz(struct EFp12 *P,mpz_t a){
    Fp12_set_mpz(&P->x,a);
    Fp12_set_mpz(&P->y,a);
    P->flag=0;
}
void EFp12_set_neg(struct EFp12 *P,struct EFp12 *A){
    Fp12_set(&P->x,&A->x);
    Fp12_set_neg(&P->y,&A->y);
    P->flag=A->flag;
}
/*---------------------------------clear---------------------------------*/
void EFp12_clear(struct EFp12 *P){
    Fp12_clear(&P->x);
    Fp12_clear(&P->y);
}
/*---------------------------------print---------------------------------*/
void EFp12_printf(struct EFp12 *P,char *name){
    printf("%s",name);
    if(P->flag==0){
        printf("(");
        Fp12_printf(&P->x,"X");
        printf("\n");
        Fp12_printf(&P->y,"Y");
        printf(")");
    }else{
        printf("0");
    }
}
/*-----------------------------rational point---------------------------*/
void EFp12_rational_point(struct EFp12 *P){
    struct Fp12 buf1,buf2,R;
    Fp12_init(&buf1);
    Fp12_init(&buf2);
    Fp12_init(&R);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    while(1){
        Fp12_random(&P->x,state);
        Fp12_squaring(&buf1,&P->x);
        Fp12_mul(&buf2,&buf1,&P->x);
        Fp12_mul_mpz(&buf1,&P->x,curve_parameter_A);
        Fp12_add(&R,&buf1,&buf2);
        mpz_add(R.x0.x0.x0.x0,R.x0.x0.x0.x0,curve_parameter_B);
        if(Fp12_legendre(&R)==1){
            Fp12_sqrt(&P->y,&R);
            break;
        }
    }
    
    Fp12_clear(&buf1);
    Fp12_clear(&buf2);
    Fp12_clear(&R);
}
void EFp12_generate_G1(struct EFp12 *P){
    struct EFp g1;
    EFp_init(&g1);
    
    EFp_rational_point(&g1);
    EFp12_set_ui(P,0);
    Fp_set(&P->x.x0.x0.x0,&g1.x);
    Fp_set(&P->y.x0.x0.x0,&g1.y);
    P->flag=g1.flag;
    
    EFp_clear(&g1);
}
void EFp12_generate_G2(struct EFp12 *Q){
    struct EFp12 random_P,P,frobenius_P;
    EFp12_init(&random_P);
    EFp12_init(&P);
    EFp12_init(&frobenius_P);
    mpz_t exp;
    mpz_init(exp);
    
    EFp12_rational_point(&random_P);
    mpz_pow_ui(exp,EFp_order,2);
    mpz_tdiv_q(exp,EFp12_total,exp);
    EFp12_SCM(&P,&random_P,exp);
    EFp12_frobenius_1(&frobenius_P,&P);
    EFp12_set_neg(&P,&P);
    EFp12_ECA(Q,&P,&frobenius_P);
    
    mpz_clear(exp);
    EFp12_clear(&random_P);
    EFp12_clear(&P);
    EFp12_clear(&frobenius_P);
}
/*---------------------------------SCM----------------------------------*/
void EFp12_ECD(struct EFp12 *ANS,struct EFp12 *P){
    if(Fp12_cmp_zero(&P->y)==0){
        ANS->flag=1;
        return;
    }
    
    struct EFp12 Tmp;
    EFp12_init(&Tmp);
    EFp12_set(&Tmp,P);
    struct Fp12 Buf1,Buf2,C;
    Fp12_init(&Buf1);
    Fp12_init(&Buf2);
    Fp12_init(&C);
    
    Fp12_mul_ui(&Buf1,&Tmp.y,2);
    Fp12_inv(&Buf1,&Buf1);
    //Fp12_mul(&Buf2,&Tmp.x,&Tmp.x);
    Fp12_squaring(&Buf2,&Tmp.x);
    Fp12_mul_ui(&Buf2,&Buf2,3);
    mpz_add(Buf2.x0.x0.x0.x0,Buf2.x0.x0.x0.x0,curve_parameter_A);
    Fp12_mul(&C,&Buf1,&Buf2);
    //Fp12_mul(&Buf1,&C,&C);
    Fp12_squaring(&Buf1,&C);
    Fp12_mul_ui(&Buf2,&Tmp.x,2);
    Fp12_sub(&ANS->x,&Buf1,&Buf2);
    Fp12_sub(&Buf1,&Tmp.x,&ANS->x);
    Fp12_mul(&Buf2,&C,&Buf1);
    Fp12_sub(&ANS->y,&Buf2,&Tmp.y);
    
    Fp12_clear(&Buf1);
    Fp12_clear(&Buf2);
    Fp12_clear(&C);
    EFp12_clear(&Tmp);
}
void EFp12_ECA(struct EFp12 *ANS,struct EFp12 *P1,struct EFp12 *P2){
    if(P1->flag==1){
        EFp12_set(ANS,P2);
        return;
    }else if(P2->flag==1){
        EFp12_set(ANS,P1);
        return;
    }else if(Fp12_cmp(&P1->x,&P2->x)==0){
        if(Fp12_cmp(&P1->y,&P2->y)!=0){
            ANS->flag=1;
            return;
        }else{
            EFp12_ECD(ANS,P1);
            return;
        }
    }
    
    struct EFp12 Tmp1,Tmp2;
    EFp12_init(&Tmp1);
    EFp12_set(&Tmp1,P1);
    EFp12_init(&Tmp2);
    EFp12_set(&Tmp2,P2);
    struct Fp12 Buf1,Buf2,C;
    Fp12_init(&Buf1);
    Fp12_init(&Buf2);
    Fp12_init(&C);
    
    Fp12_sub(&Buf1,&Tmp2.x,&Tmp1.x);
    Fp12_inv(&Buf1,&Buf1);
    Fp12_sub(&Buf2,&Tmp2.y,&Tmp1.y);
    Fp12_mul(&C,&Buf1,&Buf2);
    //Fp12_mul(&Buf1,&C,&C);
    Fp12_squaring(&Buf1,&C);
    Fp12_sub(&Buf2,&Buf1,&Tmp1.x);
    Fp12_sub(&ANS->x,&Buf2,&Tmp2.x);
    Fp12_sub(&Buf1,&Tmp1.x,&ANS->x);
    Fp12_mul(&Buf2,&C,&Buf1);
    Fp12_sub(&ANS->y,&Buf2,&Tmp1.y);
    
    //clear
    Fp12_clear(&Buf1);
    Fp12_clear(&Buf2);
    Fp12_clear(&C);
    EFp12_clear(&Tmp1);
    EFp12_clear(&Tmp2);
}
void EFp12_SCM(struct EFp12 *ANS,struct EFp12 *P,mpz_t R){
    if(mpz_cmp_ui(R,0)==0){
        ANS->flag=1;
        return;
    }else if(mpz_cmp_ui(R,1)==0){
        EFp12_set(ANS,P);
        return;
    }
    
    struct EFp12 Tmp,next_P;
    EFp12_init(&Tmp);
    EFp12_set(&Tmp,P);
    EFp12_init(&next_P);
    int i,length;
    length=(int)mpz_sizeinbase(R,2);
    char binary[length];
    mpz_get_str(binary,2,R);
    
    EFp12_set(&next_P,&Tmp);
    for(i=1; binary[i]!='\0'; i++){
        EFp12_ECD(&next_P,&next_P);
        if(binary[i]=='1'){
            EFp12_ECA(&next_P,&next_P,&Tmp);
        }
    }
    EFp12_set(ANS,&next_P);
    
    EFp12_clear(&next_P);
    EFp12_clear(&Tmp);
}

/*--------------------------------G2 SCM--------------------------------*/
void EFp12_G2_SCM_normal(struct EFp12 *ANS,struct EFp12 *Q,mpz_t S){
    struct EFp2 twisted_Q;
    EFp2_init(&twisted_Q);
    
    EFp12_to_EFp2(&twisted_Q,Q);
    EFp2_SCM(&twisted_Q,&twisted_Q,S);
    EFp2_to_EFp12(ANS,&twisted_Q);
    
    EFp2_clear(&twisted_Q);
}
/*-----------------------------frobenius--------------------------------*/
void EFp12_frobenius_1(struct EFp12 *ANS,struct EFp12 *P){
    Fp12_frobenius_1(&ANS->x,&P->x);
    Fp12_frobenius_1(&ANS->y,&P->y);
}
void EFp12_frobenius_2(struct EFp12 *ANS,struct EFp12 *P){
    Fp12_frobenius_2(&ANS->x,&P->x);
    Fp12_frobenius_2(&ANS->y,&P->y);
}
void EFp12_frobenius_3(struct EFp12 *ANS,struct EFp12 *P){
    Fp12_frobenius_3(&ANS->x,&P->x);
    Fp12_frobenius_3(&ANS->y,&P->y);
}
void EFp12_frobenius_4(struct EFp12 *ANS,struct EFp12 *P){
    Fp12_frobenius_4(&ANS->x,&P->x);
    Fp12_frobenius_4(&ANS->y,&P->y);
}
void EFp12_frobenius_6(struct EFp12 *ANS,struct EFp12 *P){
    Fp12_frobenius_6(&ANS->x,&P->x);
    Fp12_frobenius_6(&ANS->y,&P->y);
}
void EFp12_frobenius_8(struct EFp12 *ANS,struct EFp12 *P){
    Fp12_frobenius_8(&ANS->x,&P->x);
    Fp12_frobenius_8(&ANS->y,&P->y);
}
void EFp12_frobenius_10(struct EFp12 *ANS,struct EFp12 *P){
    Fp12_frobenius_10(&ANS->x,&P->x);
    Fp12_frobenius_10(&ANS->y,&P->y);
}
/*============================================================================*/
/* sextic twist                                                               */
/*============================================================================*/
void EFp12_to_EFp2(struct EFp2 *ANS,struct EFp12 *P){
    Fp2_set_ui(&ANS->x,0);
    Fp2_set(&ANS->x,&P->x.x0.x2);
    Fp2_mul_basis(&ANS->x,&ANS->x);
    Fp2_set_ui(&ANS->y,0);
    Fp2_set(&ANS->y,&P->y.x1.x1);
    Fp2_mul_basis(&ANS->y,&ANS->y);
    ANS->flag=P->flag;
}
void EFp2_to_EFp12(struct EFp12 *ANS,struct EFp2 *P){
    Fp12_set_ui(&ANS->x,0);
    Fp2_set(&ANS->x.x0.x2,&P->x);
    Fp2_inv_basis(&ANS->x.x0.x2,&ANS->x.x0.x2);
    Fp12_set_ui(&ANS->y,0);
    Fp2_set(&ANS->y.x1.x1,&P->y);
    Fp2_inv_basis(&ANS->y.x1.x1,&ANS->y.x1.x1);
    ANS->flag=P->flag;
}
/*============================================================================*/
/* final exp                                                                  */
/*============================================================================*/
void Final_exp_female_researchers_algo(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp12 tmp,t0,t1,t2,t3,t4,t5, test;
    Fp12_init(&tmp);
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t5);
    Fp12_init(&t4);
    Fp12_init(&test);
    mpz_t positive_X,positive_X2;
    mpz_init(positive_X);
    mpz_init(positive_X2);
    
    mpz_neg(positive_X,mother_parameter);
    //mpz_add(positive_X2,positive_X,positive_X);
    mpz_set_str(positive_X2,"75557863162960075030528",10);
    //gmp_printf("X=%Zd\n positive_X =%Zd\n",X, positive_X2);
    
    //f←f^(p^6)*f^-1
    Fp12_frobenius_6(&t0,A);//f^(p^6)
    Fp12_inv(&t1,A);//f^-1
    Fp12_mul(&tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_2(&t0,&tmp);//f^(p^2)
    Fp12_mul(&tmp,&t0,&tmp);//f^(p^2)*f
    
    Fp12_squaring(&t0, &tmp);
    Fp12_pow(&t1, &t0, positive_X);
    Fp12_frobenius_6(&t1, &t1);
    
    Fp12_pow(&t2,&t1,positive_X2);//t2:=t1^(u2);
    Fp12_frobenius_6(&t2,&t2);
    Fp12_frobenius_6(&t3,&tmp);//t3:=f^(-1);
    Fp12_mul(&t1,&t3,&t1);//t1:=t3*t1;
    Fp12_frobenius_6(&t1,&t1);//t1:=t1^(-1);
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    Fp12_pow(&t2,&t1,positive_X);//t2:=t1^(u);
    Fp12_frobenius_6(&t2,&t2);
    Fp12_pow(&t3,&t2,positive_X);//t3:=t2^(u);
    Fp12_frobenius_6(&t3,&t3);
    Fp12_frobenius_6(&t1,&t1);//t1:=t1^(-1);
    
    Fp12_mul(&t3,&t1,&t3);//t3:=t1*t3;
    Fp12_frobenius_6(&t1,&t1);//t1:=t1^(-1);
    Fp12_frobenius_3(&t1,&t1);//t1:=t1^(p^3);
    Fp12_frobenius_2(&t2,&t2);//t2:=t2^(p^2);
    
    
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    Fp12_pow(&t2,&t3,positive_X);//t2:=t3^(u);
    Fp12_frobenius_6(&t2,&t2);
    Fp12_mul(&t2,&t2,&t0);//t2:=t2*t0;
    Fp12_mul(&t2,&t2,&tmp);//t2:=t2*f;
    Fp12_mul(&t1,&t1,&t2);//t1:=t1*t2;
    
    Fp12_frobenius_1(&t2,&t3);//t2:=t3^p;
    Fp12_mul(ANS,&t1,&t2);//t1:=t1*t2;
    
    Fp12_clear(&tmp);
    Fp12_clear(&t0);
    Fp12_clear(&t1);
    Fp12_clear(&t2);
    Fp12_clear(&t3);
    mpz_clear(positive_X);
    mpz_clear(positive_X2);
}
void Final_exp_normal(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp12 Tmp,Buf1,Buf2;
    Fp12_init(&Tmp);
    Fp12_set(&Tmp,A);
    Fp12_init(&Buf1);
    Fp12_init(&Buf2);
    mpz_t exp,buf;
    mpz_init(exp);
    mpz_init(buf);
    
    Fp12_frobenius_6(&Buf1,&Tmp);
    Fp12_inv(&Buf2,&Tmp);
    Fp12_mul(&Tmp,&Buf1,&Buf2);
    
    Fp12_frobenius_2(&Buf1,&Tmp);
    Fp12_mul(&Tmp,&Buf1,&Tmp);
    
    mpz_pow_ui(exp,prime,4);
    mpz_pow_ui(buf,prime,2);
    mpz_sub(exp,exp,buf);
    mpz_add_ui(exp,exp,1);
    mpz_tdiv_q(exp,exp,EFp_order);
    Fp12_pow(ANS,&Tmp,exp);
    
    
    mpz_clear(exp);
    mpz_clear(buf);
    Fp12_clear(&Tmp);
    Fp12_clear(&Buf1);
    Fp12_clear(&Buf2);
}
/*============================================================================*/
/* Tate pairing                                                               */
/*============================================================================*/

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void ff_ltt_vtt(struct Fp12 *f,struct EFp12 *T,struct EFp12 *Q){
    
}
void f_ltp_vtp(struct Fp12 *f,struct EFp12 *T,struct EFp12 *P,struct EFp12 *Q){
    
}
void Miller_algo_for_tate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P){
    
}
void Tate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P){
    
}

/*============================================================================*/
/* Pseudo 8-sparse                                                            */
/*============================================================================*/
void Pseudo_8_sparse_mapping(struct EFp *P,struct EFp2 *Q,struct Fp *L){
    struct EFp2 Tmp_Q;
    EFp2_init(&Tmp_Q);
    struct EFp Tmp_P;
    EFp_init(&Tmp_P);
    struct Fp A,B,C,D,c;
    Fp_init(&A);
    Fp_init(&B);
    Fp_init(&C);
    Fp_init(&D);
    Fp_init(&c);
    
    EFp_set(&Tmp_P,P);
    EFp2_set(&Tmp_Q,Q);
    
    Fp_mul(&A,&Tmp_P.x,&Tmp_P.y);
    Fp_inv(&A,&A);
    Fp_mul(&B,&Tmp_P.x,&Tmp_P.x);
    Fp_mul(&B,&B,&A);
    Fp_mul(&C,&Tmp_P.y,&A);
    Fp_mul(&D,&B,&B);
    
    Fp2_mul_mpz(&Q->x,&Tmp_Q.x,D.x0);
    Fp_mul(&c,&B,&D);
    Fp2_mul_mpz(&Q->y,&Tmp_Q.y,c.x0);
    
    Fp_mul(&P->x,&D,&Tmp_P.x);
    Fp_set(&P->y,&P->x);
    
    Fp_mul(L,&C,&Tmp_P.y);
    Fp_mul(L,L,L);
    Fp_mul(L,L,&C);
    
    
    EFp2_clear(&Tmp_Q);
    EFp_clear(&Tmp_P);
    Fp_clear(&A);
    Fp_clear(&B);
    Fp_clear(&C);
    Fp_clear(&D);
    Fp_clear(&c);
}
void Pseudo_8_sparse_mul(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
    //A= f0 + f1γ^2 + f2γ^4 + f3γ　+ f4γ^3 + f5γ^5
    //B= 1  +                        aγ^3 +  bγ^5
    // x0.x0  x0.x1   x0.x2  x1.x0  x1.x1   x1.x2
    struct Fp12 ans;
    Fp12_init(&ans);
    struct Fp2 tmp0,tmp1,tmp2,tmp3;
    Fp2_init(&tmp0);
    Fp2_init(&tmp1);
    Fp2_init(&tmp2);
    Fp2_init(&tmp3);
    
    Fp2_mul(&tmp0,&A->x0.x0,&B->x1.x1);		//tmp0←a*f0
    Fp2_mul(&tmp1,&A->x0.x1,&B->x1.x2);		//tmp1←b*f1
    Fp2_add(&tmp2,&A->x0.x0,&A->x0.x1);		//tmp2←f0+f1
    Fp2_add(&tmp3,&B->x1.x1,&B->x1.x2);		//tmp3←a+b
    Fp2_mul(&tmp2,&tmp2,&tmp3);			//tmp2←tmp2*tmp3
    Fp2_sub(&tmp2,&tmp2,&tmp0);			//tmp2←tmp2-tmp0
    Fp2_sub(&tmp2,&tmp2,&tmp1);			//tmp2←tmp2-tmp1
    
    Fp2_add(&ans.x1.x2,&tmp2,&A->x1.x2);	//ans[γ^5]←tmp2+f5
    Fp2_add(&ans.x1.x1,&tmp0,&A->x1.x1);	//ans[γ^3]←tmp0+f4
    Fp2_mul(&tmp2,&A->x0.x2,&B->x1.x2);		//tmp2←b*f2
    Fp2_mul_basis(&tmp2,&tmp2);			//tmp2←tmp2*α
    Fp2_add(&ans.x1.x1,&ans.x1.x1,&tmp2);	//ans[γ^3]←ans[γ^3]+tmp2
    Fp2_mul(&tmp0,&A->x0.x2,&B->x1.x1);		//tmp0←a*f2
    Fp2_add(&tmp0,&tmp0,&tmp1);			//tmp0←tmp0+tmp1
    Fp2_mul_basis(&tmp0,&tmp0);			//tmp0←tmp0*α
    Fp2_add(&ans.x1.x0,&tmp0,&A->x1.x0);	//ans[γ]←tmp0+f3
    
    Fp2_mul(&tmp0,&A->x1.x0,&B->x1.x1);		//tmp0←a*f3
    Fp2_mul(&tmp1,&A->x1.x1,&B->x1.x2);		//tmp1←b*f4
    Fp2_add(&tmp2,&A->x1.x0,&A->x1.x1);		//tmp2←f3+f4
    Fp2_mul(&tmp2,&tmp2,&tmp3);			//tmp2←tmp2+tmp3
    Fp2_sub(&tmp2,&tmp2,&tmp0);			//tmp2←tmp2-tmp0
    Fp2_sub(&tmp2,&tmp2,&tmp1);			//tmp2←tmp2-tmp1
    
    Fp2_mul_basis(&tmp2,&tmp2);			//tmp2←tmp2*α
    Fp2_add(&ans.x0.x0,&tmp2,&A->x0.x0);	//ans[1]←tmp2+f0
    
    Fp2_mul(&tmp2,&A->x1.x2,&B->x1.x1);		//tmp2←a*f5
    Fp2_add(&tmp2,&tmp1,&tmp2);			//tmp2←tmp1+tmp2
    Fp2_mul_basis(&tmp2,&tmp2);			//tmp2←tmp2*α
    Fp2_add(&ans.x0.x1,&tmp2,&A->x0.x1);	//ans[γ^2]←tmp2+f1
    Fp2_mul(&tmp3,&A->x1.x2,&B->x1.x2);		//tmp3←b*f5
    Fp2_mul_basis(&tmp3,&tmp3);			//tmp3←tmp3*α
    
    Fp2_add(&tmp0,&tmp0,&tmp3);			//tmp0←tmp0+tmp3
    Fp2_add(&ans.x0.x2,&tmp0,&A->x0.x2);	//ans[γ^2]←tmp0+f2
    
    Fp12_set(ANS,&ans);
    
    Fp12_clear(&ans);
    Fp2_clear(&tmp0);
    Fp2_clear(&tmp1);
    Fp2_clear(&tmp2);
    Fp2_clear(&tmp3);
}
void ff_ltt(struct Fp12 *f,struct EFp2 *T,struct EFp *P,struct Fp *L){
    struct EFp2 Tmp_T;
    EFp2_init(&Tmp_T);
    struct Fp12 ff,ltt;
    Fp12_init(&ff);
    Fp12_init(&ltt);
    struct Fp2 A,B,C,D,E;
    Fp2_init(&A);
    Fp2_init(&B);
    Fp2_init(&C);
    Fp2_init(&D);
    Fp2_init(&E);
    EFp2_set(&Tmp_T,T);
    
    Fp12_squaring(&ff,f);
    
    //ltt
    Fp2_add(&A,&Tmp_T.y,&Tmp_T.y);		//A=1/(2*T.y)
    Fp2_inv(&A,&A);
    Fp2_squaring(&B,&Tmp_T.x);			//B=3(T.x)^2
    Fp2_mul_ui(&B,&B,3);
    Fp2_mul(&C,&A,&B);				//C=A*B
    Fp2_add(&D,&Tmp_T.x,&Tmp_T.x);		//D=2T.x
    Fp2_squaring(&T->x,&C);				//next_T.x=C^2-D
    Fp2_sub(&T->x,&T->x,&D);
    Fp2_mul(&E,&C,&Tmp_T.x);			//E=C*T.x-T.y
    Fp2_sub(&E,&E,&Tmp_T.y);
    Fp2_mul(&T->y,&C,&T->x);			//next_T.y=E-C*next_T.x
    Fp2_sub(&T->y,&E,&T->y);
    
    //set ltt
    Fp_set_ui(&ltt.x0.x0.x0,1);
    Fp2_set_neg(&ltt.x1.x2,&C);
    Fp2_inv_basis(&ltt.x1.x2,&ltt.x1.x2);
    Fp2_mul_mpz(&ltt.x1.x1,&E,L->x0);
    Fp2_inv_basis(&ltt.x1.x1,&ltt.x1.x1);
    
    Pseudo_8_sparse_mul(f,&ff,&ltt);
    
    EFp2_clear(&Tmp_T);
    Fp2_clear(&A);
    Fp2_clear(&B);
    Fp2_clear(&C);
    Fp2_clear(&D);
    Fp2_clear(&E);
    Fp12_clear(&ff);
    Fp12_clear(&ltt);
}
void f_ltq(struct Fp12 *f,struct EFp2 *T,struct EFp2 *Q,struct EFp *P,struct Fp *L){
    struct EFp2 Tmp_T;
    EFp2_init(&Tmp_T);
    struct Fp12 ltq;
    Fp12_init(&ltq);
    struct Fp2 A,B,C,D,E;
    Fp2_init(&A);
    Fp2_init(&B);
    Fp2_init(&C);
    Fp2_init(&D);
    Fp2_init(&E);
    EFp2_set(&Tmp_T,T);
    
    //ltq
    Fp2_sub(&A,&Q->x,&Tmp_T.x);		//A=(Q->x-T.x)^-1
    Fp2_inv(&A,&A);
    Fp2_sub(&B,&Q->y,&Tmp_T.y);		//B=(Q->y-T.y)
    Fp2_mul(&C,&A,&B);			//C=A*B
    Fp2_add(&D,&Tmp_T.x,&Q->x);		//D=Q->x+T.x
    Fp2_squaring(&T->x,&C);			//next_T.x=C^2-D
    Fp2_sub(&T->x,&T->x,&D);
    Fp2_mul(&E,&C,&Tmp_T.x);		//E=C*T.x-T.y
    Fp2_sub(&E,&E,&Tmp_T.y);
    Fp2_mul(&T->y,&C,&T->x);		//next_T.y=E-C*next_T.x
    Fp2_sub(&T->y,&E,&T->y);
    
    //set ltq
    Fp_set_ui(&ltq.x0.x0.x0,1);
    Fp2_set_neg(&ltq.x1.x2,&C);
    Fp2_inv_basis(&ltq.x1.x2,&ltq.x1.x2);
    Fp2_mul_mpz(&ltq.x1.x1,&E,L->x0);
    Fp2_inv_basis(&ltq.x1.x1,&ltq.x1.x1);
    
    Pseudo_8_sparse_mul(f,f,&ltq);
    
    EFp2_clear(&Tmp_T);
    Fp12_clear(&ltq);
    Fp2_clear(&A);
    Fp2_clear(&B);
    Fp2_clear(&C);
    Fp2_clear(&D);
}


/*============================================================================*/
/* Opt-ate pairing                                                            */
/*============================================================================*/
void Miller_algo_for_opt_ate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P){
    struct EFp12 Buf;
    EFp12_init(&Buf);
    struct EFp2 T;
    EFp2_init(&T);
    struct EFp2 mapped_Q,mapped_Q_neg,mapped_Q1,mapped_Q2_neg;
    EFp2_init(&mapped_Q);
    EFp2_init(&mapped_Q_neg);
    EFp2_init(&mapped_Q1);
    EFp2_init(&mapped_Q2_neg);
    struct EFp mapped_P;
    EFp_init(&mapped_P);
    struct Fp12 f;
    Fp12_init(&f);
    struct Fp L;
    Fp_init(&L);
    int i;
    
    //set
    Fp_set(&mapped_P.x,&P->x.x0.x0.x0);	//set P
    Fp_set(&mapped_P.y,&P->y.x0.x0.x0);
    mapped_P.flag=P->flag;
    EFp12_to_EFp2(&mapped_Q,Q);//set mapped_Q
    mapped_Q.flag=Q->flag;
    Pseudo_8_sparse_mapping(&mapped_P,&mapped_Q,&L);
    
    EFp2_set_neg(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    
    EFp2_set(&T,&mapped_Q_neg);		//set T
    Fp12_set_ui(&f,0);			//set f
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=x_bit-1; i>=0; i--){
        switch(X_bit_binary[i]){
            case 0:
                ff_ltt(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt(&f,&T,&mapped_P,&L);
                f_ltq(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
        
    }
    
    Fp12_set(ANS,&f);
    
    EFp12_clear(&Buf);
    Fp12_clear(&f);
    EFp2_clear(&T);
    EFp2_clear(&mapped_Q);
    EFp2_clear(&mapped_Q_neg);
    EFp2_clear(&mapped_Q1);
    EFp2_clear(&mapped_Q2_neg);
    EFp_clear(&mapped_P);
    Fp_clear(&L);
}
void Opt_ate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P){
    struct timeval t0,t1;
    mpz_mpz_mul=0;
    mpz_ui_mul=0;
    Fp_mpz_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    gettimeofday(&t0,NULL);
    Miller_algo_for_opt_ate(ANS,Q,P);
    gettimeofday(&t1,NULL);
    printf("loop time :%.2f[ms]\n",timedifference_msec(t0,t1));
    printf("Miller loop cost\nmpz_mpz_mul:%ld,mpz_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
    //Fp12_printf(ANS,"test2:"); printf("\n");
    
    mpz_mpz_mul=0;
    mpz_ui_mul=0;
    Fp_mpz_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    gettimeofday(&t0,NULL);
    Final_exp_female_researchers_algo(ANS,ANS);
    gettimeofday(&t1,NULL);
    printf("final exp time:%.2f[ms]\n",timedifference_msec(t0,t1));
    printf("final exp cost\nmpz_mpz_mul:%ld,mul_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
}
/*============================================================================*/
/*  init                                                                      */
/*============================================================================*/
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
/*============================================================================*/
/*  clear                                                                     */
/*============================================================================*/
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
/*============================================================================*/
/*  set                                                                       */
/*============================================================================*/
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
    if(generate_prime()==1){
        generate_order();
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
    X_bit_binary[77]=-1;
    X_bit_binary[50]=1;
    X_bit_binary[33]=1;
    
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
    mpz_t result,buf1,buf2,modtest;
    mpz_init(result);
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(modtest);
    
    mpz_sub_ui(result,mother_parameter,1);
    mpz_pow_ui(result,result,2);
    
    mpz_pow_ui(buf1,mother_parameter,4);
    mpz_pow_ui(buf2,mother_parameter,2);
    mpz_sub(buf1,buf1,buf2);
    mpz_add_ui(buf1,buf1,1);
    
    mpz_mul(result,result,buf1);
    
    //check div3
    mpz_mod_ui(modtest,result,3);
    if(mpz_cmp_ui(modtest,0)!=0){
        mpz_init(result);
        mpz_init(buf1);
        mpz_init(buf2);
        mpz_init(modtest);
        printf("cannot devided by 3\n");
        return 0;
    }
    
    mpz_tdiv_q_ui(result,result,3);
    mpz_add(result,result,mother_parameter);
    
    //isprime
    if(mpz_probab_prime_p(result,25)==0){
        mpz_init(result);
        mpz_init(buf1);
        mpz_init(buf2);
        mpz_init(modtest);
        printf("not prime\n");
        return 0;
    }
    
    mpz_set(prime,result);
    
    mpz_init(result);
    mpz_init(buf1);
    mpz_init(buf2);
    mpz_init(modtest);
    return 1;
}
int generate_order(){
    mpz_t buf1,buf2;
    mpz_init(buf1);
    mpz_init(buf2);
    
    mpz_pow_ui(buf1,mother_parameter,4);
    mpz_pow_ui(buf2,mother_parameter,2);
    mpz_sub(EFp_order,buf1,buf2);
    mpz_add_ui(EFp_order,EFp_order,1);
    
    mpz_clear(buf1);
    mpz_clear(buf2);
    return 1;
}
void generate_trace(){
    mpz_add_ui(trace_t,mother_parameter,1);
}
void generate_basis(){
    //Fp_basis
    Fp_set_ui(&Fp_basis,1);
    //Fp2_basis
    Fp2_set_ui(&Fp2_basis,1);
    //Fp2_basis_inv
    Fp2_inv(&Fp2_basis_inv,&Fp2_basis);
    //Fp6_basis
    Fp6_set_ui(&Fp6_basis,0);
    Fp_set_ui(&Fp6_basis.x1.x0,1);
    
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

/*============================================================================*/
/*  print                                                                     */
/*============================================================================*/
void print_parameters(){
    printf("====================================================================================\n");
    printf("prime length:%dbit\n",(int)mpz_sizeinbase(prime,2));
    printf("E:y^2=x^3+");
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
/*============================================================================*/
/* time                                                                       */
/*============================================================================*/
float timedifference_msec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) * 1000.0f + (t1.tv_usec - t0.tv_usec) / 1000.0f;
}
float timedifference_usec(struct timeval t0, struct timeval t1){
    return (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec);
}
/*============================================================================*/
/* test                                                                       */
/*============================================================================*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void test_tate_pairing(){
    printf("====================================================================================\n");
    printf("Tate pairing\n\n");
}
void test_ate_pairing(){
    printf("====================================================================================\n");
    printf("Ate pairing\n\n");
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void test_opt_ate_pairing(){
    printf("====================================================================================\n");
    printf("Opt-Ate pairing\n\n");
    struct EFp12 P,Q,S1_P,S2_P,S1_Q,S2_Q;
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&S1_P);
    EFp12_init(&S2_P);
    EFp12_init(&S1_Q);
    EFp12_init(&S2_Q);
    struct Fp12 Z,Test1,Test2,Test3;
    Fp12_init(&Z);
    Fp12_init(&Test1);
    Fp12_init(&Test2);
    Fp12_init(&Test3);
    mpz_t S1,S2,S12;
    mpz_init(S1);
    mpz_init(S2);
    mpz_init(S12);
    
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    printf("*********scalar*********\n");
    mpz_urandomm(S1,state,EFp_order);	//S1
    printf("S1:");
    gmp_printf("%Zd",S1);
    printf("\n");
    mpz_urandomm(S2,state,EFp_order);	//S2
    printf("S2:");
    gmp_printf("%Zd",S2);
    printf("\n");
    mpz_mul(S12,S1,S2);			//S12
    mpz_mod(S12,S12,EFp_order);
    printf("S12:");
    gmp_printf("%Zd",S12);
    printf("\n\n");
    
    printf("*********G1 & G2*********\n");
    EFp12_generate_G1(&P);			//P
    //mpz_set_str(P.x.x0.x0.x0.x0,"1577263467895074691751656396777546817676899779239812956811535508811145903113183105315855416956275013763132428190349521477916841537817217749",10);
    //mpz_set_str(P.y.x0.x0.x0.x0,"1482623217062500042620995282298816432047975830544022708334643354696131009901634528980080824925757942965723245158102345395591975011061998174",10);
    EFp12_printf(&P,"P:");
    printf("\n");
    EFp12_generate_G2(&Q);			//Q
    //mpz_set_str(Q.x.x0.x2.x0.x0,"1716096407462707739249584854219035887699085689300518488744116638769194609975441542417090151676720383761847673487573410049924094001307902771",10);
    //mpz_set_str(Q.x.x0.x2.x1.x0,"1884774744539223134856244661850811554160388946147660758387803399777014198117917660870489415649762044750981220496189685117651039305687807588",10);
    //mpz_set_str(Q.y.x1.x1.x0.x0,"862911509351712110897605125200443999396611767860495548282350223126731430707014429484782130102632482076556245081709440912377961089052176300",10);
    //mpz_set_str(Q.y.x1.x1.x1.x0,"14917642869256953479675328417653271711757572421578255369823731891492952725387400103789740069957122839470285312924022325651200303158794609",10);
    EFp12_printf(&Q,"Q:");
    printf("\n\n");
    
    printf("*********calculate [S1]P,[S2]P,[S1]Q,[S2]Q*********\n");
    EFp12_SCM(&S1_P,&P,S1);		//S1_P
    EFp12_SCM(&S2_P,&P,S2);		//S2_P
    EFp12_SCM(&S1_Q,&Q,S1);		//S1_Q
    EFp12_SCM(&S2_Q,&Q,S2);		//S2_Q
    printf("\n");
    
    
    printf("-------------------------------opt-ate--------------------------------\n");
    Opt_ate_pairing(&Z,&Q,&P);
    //linary test
    printf("*********linearity test**********\n");
    Fp12_pow(&Test1,&Z,S12);					//Test1 Z^S12
    Opt_ate_pairing(&Test2,&S2_Q,&S1_P);			//Test2 S1_P,S2_Q
    Fp12_printf(&Test2,""); printf("\n\n");
    Opt_ate_pairing(&Test3,&S1_Q,&S2_P);			//Test3 S2_Q,S1_P
    Fp12_printf(&Test3,""); printf("\n\n");
    printf("*********************************\n");
    if(Fp12_cmp(&Test1,&Fp12_ZERO)!=0 && Fp12_cmp(&Test1,&Fp12_ONE)!=0 && Fp12_cmp(&Test1,&Test2)==0 && Fp12_cmp(&Test2,&Test3)==0 && Fp12_cmp(&Test3,&Test1)==0){
        printf("test success\n\n");
    }else{
        printf("test failed\n\n");
    }
    
    mpz_clear(S1);
    mpz_clear(S2);
    mpz_clear(S12);
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&S1_P);
    EFp12_clear(&S2_P);
    EFp12_clear(&S1_Q);
    EFp12_clear(&S2_Q);
    Fp12_clear(&Z);
    Fp12_clear(&Test1);
    Fp12_clear(&Test2);
    Fp12_clear(&Test3);
}

void test_frobenius(){
    printf("====================================================================================\n");
    printf("frobenius mapping\n");
    struct Fp12 P,Q;
    Fp12_init(&P);
    Fp12_init(&Q);
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    
    mpz_t exp;
    mpz_init(exp);
    
    Fp12_random(&P,state);
    
    Fp12_frobenius_1(&Q,&P);
    Fp12_printf(&Q,""); printf("\n");
    mpz_pow_ui(exp,prime,1);
    Fp12_pow(&P,&P,exp);
    Fp12_printf(&P,""); printf("\n");
    if(Fp12_cmp(&P,&Q)==0){
        printf("success.\n");
    }
    
    mpz_clear(exp);	
    Fp12_clear(&P);
    Fp12_clear(&Q);
}
void check_num_of_Fp_mul(){
    struct Fp2 P;
    Fp2_init(&P);
    struct Fp6 Q;
    Fp6_init(&Q);
    struct Fp12 R,Buf;
    Fp12_init(&R);
    Fp12_init(&Buf);
    
    gmp_randstate_t state;
    gmp_randinit_default (state);
    gmp_randseed_ui(state,(unsigned long)time(NULL));
    Fp2_random(&P,state);
    Fp6_random(&Q,state);
    Fp12_random(&R,state);
    
    Fp2_clear(&P);
    Fp6_clear(&Q);
    Fp12_clear(&R);
    Fp12_clear(&Buf);
}
