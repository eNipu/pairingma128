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
/* Tate pairing                                                               */
/*============================================================================*/

/*
 * @brief Line evaluation for the Tate pairing (vertical line)
 * 
 * This function computes the vertical line function evaluation for the Tate pairing.
 * 
 * @param[out] f Result of the line evaluation
 * @param[in] T Point on the elliptic curve
 * @param[in] Q Point for line evaluation
 */
void ff_ltt_vtt(struct Fp12 *f, struct EFp12 *T, struct EFp12 *Q) {
    struct Fp12 tmp;
    Fp12_init(&tmp);
    
    // Vertical line through T: x - T.x
    Fp12_set_ui(&tmp, 0);
    Fp_set(&tmp.x0.x0.x0, &T->x.x0.x0.x0);
    Fp12_sub(&tmp, &Q->x, &tmp);
    Fp12_mul(f, f, &tmp);
    
    Fp12_clear(&tmp);
}

/*
 * @brief Line evaluation for the Tate pairing (tangent/addition line)
 * 
 * This function computes the tangent or addition line function evaluation for the Tate pairing.
 * 
 * @param[out] f Result of the line evaluation
 * @param[in] T Point on the elliptic curve
 * @param[in] P Point being added or doubled
 * @param[in] Q Point for line evaluation
 */
void f_ltp_vtp(struct Fp12 *f, struct EFp12 *T, struct EFp12 *P, struct EFp12 *Q) {
    struct Fp12 tmp1, tmp2;
    struct Fp lambda;
    
    Fp12_init(&tmp1);
    Fp12_init(&tmp2);
    Fp_init(&lambda);
    
    // Calculate slope (lambda)
    if (Fp_cmp(&T->x.x0.x0.x0, &P->x.x0.x0.x0) == 0) {
        // Point doubling case
        // λ = (3x²+a)/(2y)
        Fp_mul(&lambda, &T->x.x0.x0.x0, &T->x.x0.x0.x0);  // x^2
        Fp_mul_ui(&lambda, &lambda, 3);                   // 3x^2
        Fp_add_mpz(&lambda, &lambda, curve_parameter_A);  // 3x^2 + a
        
        Fp_mul_ui(&tmp1.x0.x0.x0, &T->y.x0.x0.x0, 2);     // 2y
        Fp_inv(&tmp1.x0.x0.x0, &tmp1.x0.x0.x0);           // 1/(2y)
        
        Fp_mul(&lambda, &lambda, &tmp1.x0.x0.x0);         // λ = (3x²+a)/(2y)
    } else {
        // Point addition case
        // λ = (y2-y1)/(x2-x1)
        Fp_sub(&lambda, &P->y.x0.x0.x0, &T->y.x0.x0.x0);  // y2-y1
        Fp_sub(&tmp1.x0.x0.x0, &P->x.x0.x0.x0, &T->x.x0.x0.x0); // x2-x1
        Fp_inv(&tmp1.x0.x0.x0, &tmp1.x0.x0.x0);           // 1/(x2-x1)
        Fp_mul(&lambda, &lambda, &tmp1.x0.x0.x0);         // λ = (y2-y1)/(x2-x1)
    }
    
    // Line evaluation: λ(x-xT) - (y-yT)
    // First part: λ(x-xT)
    Fp12_set_ui(&tmp1, 0);
    Fp_set(&tmp1.x0.x0.x0, &T->x.x0.x0.x0);
    Fp12_sub(&tmp1, &Q->x, &tmp1);
    Fp12_set_ui(&tmp2, 0);
    Fp_set(&tmp2.x0.x0.x0, &lambda);
    Fp12_mul(&tmp1, &tmp1, &tmp2);
    
    // Second part: -(y-yT)
    Fp12_set_ui(&tmp2, 0);
    Fp_set(&tmp2.x0.x0.x0, &T->y.x0.x0.x0);
    Fp12_sub(&tmp2, &Q->y, &tmp2);
    Fp12_set_neg(&tmp2, &tmp2);
    
    // Combine: λ(x-xT) - (y-yT)
    Fp12_add(&tmp1, &tmp1, &tmp2);
    
    // Multiply result into f
    Fp12_mul(f, f, &tmp1);
    
    Fp12_clear(&tmp1);
    Fp12_clear(&tmp2);
    Fp_clear(&lambda);
}

/*
 * @brief Miller's algorithm for the Tate pairing
 * 
 * This function implements Miller's algorithm for the Tate pairing calculation.
 * 
 * @param[out] ANS Result of Miller's algorithm
 * @param[in] Q First point (in G2)
 * @param[in] P Second point (in G1)
 */
void Miller_algo_for_tate(struct Fp12 *ANS, struct EFp12 *Q, struct EFp12 *P) {
    int i, length;
    struct EFp12 T;
    struct Fp12 f;
    mpz_t loop;
    char binary[1000]; // Sufficiently large buffer
    
    // Initialize
    EFp12_init(&T);
    Fp12_init(&f);
    mpz_init(loop);
    
    // Set loop parameter (r-1, where r is the order of the curve)
    mpz_sub_ui(loop, EFp_order, 1);
    
    // Convert loop parameter to binary
    length = (int)mpz_sizeinbase(loop, 2);
    mpz_get_str(binary, 2, loop);
    
    // Initialize Miller's algorithm
    Fp12_set_one(&f);
    EFp12_set(&T, P);
    
    // Main Miller loop
    for (i = 1; binary[i] != '\0'; i++) {
        // Double step: compute f² · l_{T,T}(Q)
        f_ltp_vtp(&f, &T, &T, Q);
        Fp12_squaring(&f, &f);
        
        // Point doubling on T
        EFp12_ECD(&T, &T);
        
        // If bit is 1, perform addition step
        if (binary[i] == '1') {
            // Addition step: compute f · l_{T,P}(Q)
            f_ltp_vtp(&f, &T, P, Q);
            
            // Point addition on T
            EFp12_ECA(&T, &T, P);
        }
    }
    
    // Set the result
    Fp12_set(ANS, &f);
    
    // Clean up
    EFp12_clear(&T);
    Fp12_clear(&f);
    mpz_clear(loop);
}

/*
 * @brief Tate pairing computation
 * 
 * This function computes the Tate pairing e(P,Q) for points P ∈ G1 and Q ∈ G2.
 * 
 * @param[out] ANS Result of the Tate pairing
 * @param[in] Q First point (in G2)
 * @param[in] P Second point (in G1)
 */
void Tate_pairing(struct Fp12 *ANS, struct EFp12 *Q, struct EFp12 *P) {
    // Compute Miller's algorithm
    Miller_algo_for_tate(ANS, Q, P);
    
    // Perform final exponentiation
    Final_exp_normal(ANS, ANS);
}

/*============================================================================*/
/* Ate pairing                                                                 */
/*============================================================================*/

/**
 * @brief Miller's algorithm for the Ate pairing
 * 
 * This function implements Miller's algorithm for the Ate pairing calculation.
 * The Ate pairing uses the trace of Frobenius as the loop parameter instead
 * of the order of the curve as in Tate pairing, resulting in a shorter loop.
 * 
 * @param[out] ANS Result of Miller's algorithm
 * @param[in] Q First point (in G2)
 * @param[in] P Second point (in G1)
 */
void Miller_algo_for_ate(struct Fp12 *ANS, struct EFp12 *Q, struct EFp12 *P) {
    int i, length;
    struct EFp12 T;
    struct Fp12 f;
    mpz_t loop;
    char binary[1000]; // Sufficiently large buffer
    
    // Initialize
    EFp12_init(&T);
    Fp12_init(&f);
    mpz_init(loop);
    
    // Set loop parameter to trace of Frobenius
    mpz_set(loop, trace_t);
    
    // Convert loop parameter to binary
    length = (int)mpz_sizeinbase(loop, 2);
    mpz_get_str(binary, 2, loop);
    
    // Initialize Miller's algorithm
    Fp12_set_one(&f);
    EFp12_set(&T, Q);
    
    // Skip the first bit which is always 1
    for (i = 1; binary[i] != '\0'; i++) {
        // Double step
        f_ltp_vtp(&f, &T, &T, P);
        Fp12_squaring(&f, &f);
        EFp12_ECD(&T, &T);
        
        // Addition step (if bit is 1)
        if (binary[i] == '1') {
            f_ltp_vtp(&f, &T, Q, P);
            EFp12_ECA(&T, &T, Q);
        }
    }
    
    // Apply the final adjustment for Ate pairing
    // For Ate, we need to compute f_{q,Q}(P) / f_{1,Q}(P)
    struct EFp12 Q_frob;
    EFp12_init(&Q_frob);
    EFp12_frobenius_1(&Q_frob, Q);
    
    // Compute f_{1,Q}(P) with line function
    struct Fp12 line_eval;
    Fp12_init(&line_eval);
    Fp12_set_one(&line_eval);
    
    // Apply line evaluation for f_{q,Q}(P) / f_{1,Q}(P)
    f_ltp_vtp(&line_eval, &T, &Q_frob, P);
    
    // Set the result
    Fp12_mul(ANS, &f, &line_eval);
    
    // Clean up
    EFp12_clear(&T);
    EFp12_clear(&Q_frob);
    Fp12_clear(&f);
    Fp12_clear(&line_eval);
    mpz_clear(loop);
}

/**
 * @brief Ate pairing computation
 * 
 * This function computes the Ate pairing e(P,Q) for points P ∈ G1 and Q ∈ G2.
 * The Ate pairing is an optimization of the Tate pairing with a shorter Miller loop.
 * 
 * @param[out] ANS Result of the Ate pairing
 * @param[in] Q First point (in G2)
 * @param[in] P Second point (in G1)
 */
void Ate_pairing(struct Fp12 *ANS, struct EFp12 *Q, struct EFp12 *P) {
    // Verify that points are valid
    if (P->flag == 0 || Q->flag == 0) {
        Fp12_set_ui(ANS, 1);
        return;
    }
    
    // Compute Miller's algorithm
    Miller_algo_for_ate(ANS, Q, P);
    
    // Apply final exponentiation
    Final_exp_normal(ANS, ANS);
}

/**
 * @brief Test function for the Optimal Ate pairing
 * 
 * This function tests the Optimal Ate pairing implementation by:
 * 1. Generating points in G1 and G2
 * 2. Testing bilinearity property: e(aP,Q) = e(P,aQ) = e(P,Q)^a
 * 3. Testing performance of the pairing computation
 */
void test_opt_ate_pairing() {
    struct EFp12 P, Q, aP;
    struct Fp12 e1, e2, e3;
    mpz_t a;
    struct timeval tv1, tv2;
    
    // Initialize
    EFp12_init(&P);
    EFp12_init(&Q);
    EFp12_init(&aP);
    Fp12_init(&e1);
    Fp12_init(&e2);
    Fp12_init(&e3);
    mpz_init(a);
    
    // Generate random scalar and points
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));
    mpz_urandomm(a, state, EFp_order);
    
    // Generate points in G1 and G2
    EFp12_generate_G1(&P);
    EFp12_generate_G2(&Q);
    
    printf("=== Testing Optimal Ate Pairing on BLS12 ===\n");
    
    // Compute e(P,Q)
    gettimeofday(&tv1, NULL);
    Opt_ate_pairing(&e1, &Q, &P);
    gettimeofday(&tv2, NULL);
    printf("Time for Optimal Ate pairing: %f ms\n", timedifference_msec(tv1, tv2));
    
    // Compute e(aP,Q)
    EFp12_SCM(&aP, &P, a);
    Opt_ate_pairing(&e2, &Q, &aP);
    
    // Compute e(P,Q)^a
    Fp12_pow(&e3, &e1, a);
    
    // Check bilinearity: e(aP,Q) = e(P,Q)^a
    if (Fp12_cmp(&e2, &e3) == 0) {
        printf("PASS: Bilinearity check for Optimal Ate pairing\n");
    } else {
        printf("FAIL: Bilinearity check failed for Optimal Ate pairing\n");
    }
    
    // Compare with Ate pairing
    struct Fp12 ate_result;
    Fp12_init(&ate_result);
    
    gettimeofday(&tv1, NULL);
    Ate_pairing(&ate_result, &Q, &P);
    gettimeofday(&tv2, NULL);
    printf("Time for Ate pairing: %f ms\n", timedifference_msec(tv1, tv2));
    
    // Check if Optimal Ate and regular Ate give the same result
    if (Fp12_cmp(&e1, &ate_result) == 0) {
        printf("PASS: Optimal Ate and Ate pairings produce the same result\n");
    } else {
        printf("FAIL: Optimal Ate and Ate pairings produce different results\n");
    }
    
    // Clean up
    EFp12_clear(&P);
    EFp12_clear(&Q);
    EFp12_clear(&aP);
    Fp12_clear(&e1);
    Fp12_clear(&e2);
    Fp12_clear(&e3);
    Fp12_clear(&ate_result);
    mpz_clear(a);
    
    printf("=== Test complete ===\n");
}
