#include "field.h"
#include "params.h"
#include "util.h"

void Fp_init(struct Fp *P){
    mpz_init(P->x0);
}

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

void Fp_random(struct Fp *P,gmp_randstate_t state){
    mpz_urandomm(P->x0,state,prime);
}

void Fp_clear(struct Fp *P){
    mpz_clear(P->x0);
}

void Fp_printf(struct Fp *P,char *name){
    printf("%s",name);
    mpz_out_str(stdout,10,P->x0);
}

void Fp_mul(struct Fp *ANS,struct Fp *A,struct Fp *B){
    if(mpz_cmp(A->x0,B->x0)==0){
        Fp_sqr++;
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

void Fp_inv(struct Fp *ANS,struct Fp *A){
    Fp_inv_num++;
    mpz_invert(ANS->x0,A->x0,prime);
}

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

void Fp2_init(struct Fp2 *P){
    Fp_init(&P->x0);
    Fp_init(&P->x1);
}

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

void Fp2_random(struct Fp2 *P,gmp_randstate_t state){
    Fp_random(&P->x0,state);
    Fp_random(&P->x1,state);
}

void Fp2_clear(struct Fp2 *P){
    Fp_clear(&P->x0);
    Fp_clear(&P->x1);
}

void Fp2_printf(struct Fp2 *P,char *name){
    printf("%s",name);
    printf("(");
    Fp_printf(&P->x0,"");
    printf(",");
    Fp_printf(&P->x1,"");
    printf(")");
}

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

void Fp2_inv(struct Fp2 *ANS,struct Fp2 *A){
    struct Fp2 frob,buf;
    Fp2_init(&frob);
    Fp2_init(&buf);
    
    Fp2_inv_map(&frob,A);
    Fp2_mul(&buf,A,&frob);
    //Fp_inv(&buf.x0,&buf.x0);
    Fp_inv(&buf.x0,&buf.x0);
    Fp2_mul_mpz(ANS,&frob,buf.x0.x0);
    
    Fp2_clear(&frob);
    Fp2_clear(&buf);
}

void Fp2_inv_map(struct Fp2 *ANS,struct Fp2 *A){
    Fp_set(&ANS->x0,&A->x0);
    Fp_set_neg(&ANS->x1,&A->x1);
}

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

void Fp2_frobenius_1(struct Fp2 *ANS,struct Fp2 *A){
    //x0
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

void Fp6_init(struct Fp6 *P){
    Fp2_init(&P->x0);
    Fp2_init(&P->x1);
    Fp2_init(&P->x2);
}

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

void Fp6_random(struct Fp6 *P,gmp_randstate_t state){
    Fp2_random(&P->x0,state);
    Fp2_random(&P->x1,state);
    Fp2_random(&P->x2,state);
}

void Fp6_clear(struct Fp6 *P){
    Fp2_clear(&P->x0);
    Fp2_clear(&P->x1);
    Fp2_clear(&P->x2);
}

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

void Fp12_init(struct Fp12 *P){
    Fp6_init(&P->x0);
    Fp6_init(&P->x1);
}

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

void Fp12_random(struct Fp12 *P,gmp_randstate_t state){
    Fp6_random(&P->x0,state);
    Fp6_random(&P->x1,state);
}

void Fp12_clear(struct Fp12 *P){
    Fp6_clear(&P->x0);
    Fp6_clear(&P->x1);
}

void Fp12_printf(struct Fp12 *P,char *name){
    printf("%s",name);
    printf("(");
    Fp6_printf(&P->x0,"");
    printf(",");
    Fp6_printf(&P->x1,"");
    printf(")");
}

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
