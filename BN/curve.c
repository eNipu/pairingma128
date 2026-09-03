#include "curve.h"
#include "params.h"

void EFp_init(struct EFp *P){
    Fp_init(&P->x);
    Fp_init(&P->y);
    P->flag=0;
}

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

void EFp_clear(struct EFp *P){
    Fp_clear(&P->x);
    Fp_clear(&P->y);
}

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
        Fp_sub_mpz(&R,&R,curve_parameter_B);
        if(Fp_legendre(&R)==1){
            Fp_sqrt(&P->y,&R);
            break;
        }
    }
    
    Fp_clear(&buf1);
    Fp_clear(&buf2);
    Fp_clear(&R);
}

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

void EFp_skew_frobenius(struct EFp *ANS,struct EFp *A){
    Fp_mul(&ANS->x,&A->x,&inv_CNR1);
    Fp_set_neg(&ANS->y,&A->y);
}

void EFp2_init(struct EFp2 *P){
    Fp2_init(&P->x);
    Fp2_init(&P->y);
    P->flag=0;
}

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

void EFp2_clear(struct EFp2 *P){
    Fp2_clear(&P->x);
    Fp2_clear(&P->y);
}

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
        mpz_sub(R.x0.x0,R.x0.x0,curve_parameter_B);
        if(Fp2_legendre(&R)==1){
            Fp2_sqrt(&P->y,&R);
            break;
        }
    }
    
    Fp2_clear(&buf1);
    Fp2_clear(&buf2);
    Fp2_clear(&R);
}

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

void EFp2_skew_frobenius_10(struct EFp2 *ANS,struct EFp2 *A){
    //x
    Fp2_frobenius_10(&ANS->x,&A->x);
    Fp2_mul(&ANS->x,&ANS->x,&Fp2_basis_inv_prime_10_div_3);
    //y
    Fp2_frobenius_10(&ANS->y,&A->y);
    Fp2_mul(&ANS->y,&ANS->y,&Fp2_basis_inv_prime_10_div_2);
}

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

void EFp6_clear(struct EFp6 *P){
    Fp6_clear(&P->x);
    Fp6_clear(&P->y);
}

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
        mpz_sub(R.x0.x0.x0,R.x0.x0.x0,curve_parameter_B);
        if(Fp6_legendre(&R)==1){
            Fp6_sqrt(&P->y,&R);
            break;
        }
    }
    
    Fp6_clear(&buf1);
    Fp6_clear(&buf2);
    Fp6_clear(&R);
}

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

void EFp12_init(struct EFp12 *P){
    Fp12_init(&P->x);
    Fp12_init(&P->y);
    P->flag=0;
}

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

void EFp12_clear(struct EFp12 *P){
    Fp12_clear(&P->x);
    Fp12_clear(&P->y);
}

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
        mpz_sub(R.x0.x0.x0.x0,R.x0.x0.x0.x0,curve_parameter_B);
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

void EFp12_G1_SCM(struct EFp12 *ANS,struct EFp12 *P,mpz_t S){
    //struct timeval t0,t1;
    struct EFp tmp_P,skew_P,add_R,next_P,test;
    EFp_init(&tmp_P);
    EFp_init(&skew_P);
    EFp_init(&add_R);
    EFp_init(&next_P);
    EFp_init(&test);
    mpz_t s1,s2,s3,s4,s5,V1,V2,buf,A,B;
    mpz_init(s1);
    mpz_init(s2);
    mpz_init(s3);
    mpz_init(s4);
    mpz_init(s5);
    mpz_init(V1);
    mpz_init(V2);
    mpz_init(buf);
    mpz_init(A);
    mpz_init(B);
    int i,length_A,length_B,loop_length;
    
    //gettimeofday(&t0,NULL);
    //set
    Fp_set(&tmp_P.x,&P->x.x0.x0.x0);	//set tmp_P
    Fp_set(&tmp_P.y,&P->y.x0.x0.x0);
    tmp_P.flag=P->flag;
    EFp_skew_frobenius(&skew_P,&tmp_P);	//set skew_P
    
    EFp_SCM(&next_P,&skew_P,EFp_order);
    
    //set V1
    mpz_mul(V1,mother_parameter,mother_parameter);
    mpz_mul_ui(V1,V1,6);
    mpz_mul_ui(buf,mother_parameter,4);
    mpz_add(V1,V1,buf);
    mpz_add_ui(V1,V1,1);
    //set V2
    mpz_add(V2,mother_parameter,mother_parameter);
    mpz_add_ui(V2,V2,1);
    mpz_tdiv_qr(s1,s2,S,V1);	//s1,s2
    mpz_mul(buf,V2,s1);		//s3,s4
    mpz_tdiv_qr(s3,s4,buf,V1);
    mpz_mul(s5,V2,s3);		//s5
    mpz_add(A,s4,s5);			//A
    mpz_mod(A,A,EFp_order);
    mpz_sub(B,s2,s5);			//B
    mpz_mod(B,B,EFp_order);
    
    //binary A
    length_A=(int)mpz_sizeinbase(A,2);
    printf("G1 length_A:%d\n",length_A);
    //binary B
    length_B=(int)mpz_sizeinbase(B,2);
    if(length_B>400){
        printf("neg patern\n");
        mpz_sub(B,EFp_order,B);
        length_B=(int)mpz_sizeinbase(B,2);
        EFp_set_neg(&tmp_P,&tmp_P);
    }
    printf("G1 length_B:%d\n",length_B);
    
    //add_R
    EFp_ECA(&add_R,&tmp_P,&skew_P);	//set add_R
    
    //set binary
    if(length_A>length_B){
        loop_length=length_A;
    }else{
        loop_length=length_B;
    }
    char binary_A[loop_length+1];
    char binary_B[loop_length+1];
    if(length_A>length_B){
        char binary_buf[length_B+1];
        mpz_get_str(binary_buf,2,B);
        mpz_get_str(binary_A,2,A);
        memset(binary_B,'0',sizeof(binary_B));
        memmove(binary_B+length_A-length_B,binary_buf,sizeof(binary_buf));
        //set next_P
        EFp_set(&next_P,&skew_P);
    }else if(length_A<length_B){
        loop_length=length_B;
        mpz_get_str(binary_B,2,B);
        char binary_buf[length_A+1];
        mpz_get_str(binary_buf,2,A);
        memset(binary_A,'0',sizeof(binary_A));
        memmove(binary_A+length_B-length_A,binary_buf,sizeof(binary_buf));
        //set next_P
        EFp_set(&next_P,&tmp_P);
    }else{
        mpz_get_str(binary_A,2,A);
        mpz_get_str(binary_B,2,B);
        //set next_P
        EFp_set(&next_P,&add_R);
    }
    //gettimeofday(&t1,NULL);
    //printf("set time:%.2f[ms]\n",timedifference_msec(t0,t1));
    
    //gettimeofday(&t0,NULL);
    //SCM
    for(i=1; i<loop_length; i++){
        EFp_ECD(&next_P,&next_P);
        if(binary_A[i] == '0' && binary_B[i] == '1'){
            EFp_ECA(&next_P,&next_P,&tmp_P);
        }else if(binary_A[i] == '1' && binary_B[i] == '0'){
            EFp_ECA(&next_P,&next_P,&skew_P);
        }else if(binary_A[i] == '1' && binary_B[i] == '1'){
            EFp_ECA(&next_P,&next_P,&add_R);
        }
    }
    //gettimeofday(&t1,NULL);
    //printf("loop time:%.2f[ms]\n",timedifference_msec(t0,t1));
    
    Fp12_set_ui(&ANS->x,0);
    Fp12_set_ui(&ANS->y,0);
    Fp_set(&ANS->x.x0.x0.x0,&next_P.x);
    Fp_set(&ANS->y.x0.x0.x0,&next_P.y);
    ANS->flag=next_P.flag;
    
    mpz_clear(s1);
    mpz_clear(s2);
    mpz_clear(s3);
    mpz_clear(s4);
    mpz_clear(s5);
    mpz_clear(V1);
    mpz_clear(V2);
    mpz_clear(buf);
    mpz_clear(A);
    mpz_clear(B);
    EFp_clear(&test);
    EFp_clear(&tmp_P);
    EFp_clear(&skew_P);
    EFp_clear(&add_R);
    EFp_clear(&next_P);
}

void EFp12_G2_SCM_2div(struct EFp12 *ANS,struct EFp12 *Q,mpz_t S){
    //struct timeval t0,t1;
    struct EFp12 Buf;
    EFp12_init(&Buf);
    struct EFp2 twisted_Q,frobenius_Q,add_R,next_Q;
    EFp2_init(&twisted_Q);
    EFp2_init(&frobenius_Q);
    EFp2_init(&add_R);
    EFp2_init(&next_Q);
    mpz_t A,B,buf;
    mpz_init(A);
    mpz_init(B);
    mpz_init(buf);
    int i,length_A,length_B,loop_length;
    
    //gettimeofday(&t0,NULL);
    //set
    EFp12_to_EFp2(&twisted_Q,Q);				//twisted_Q
    EFp12_frobenius_1(&Buf,Q);
    EFp12_to_EFp2(&frobenius_Q,&Buf);			//frobenius_Q
    EFp2_ECA(&add_R,&twisted_Q,&frobenius_Q);		//add_R
    
    //A,B
    mpz_sub_ui(buf,trace_t,1);
    mpz_tdiv_qr(A,B,S,buf);
    
    //binary A
    length_A=(int)mpz_sizeinbase(A,2);
    printf("G2 length_A:%d\n",length_A);
    //binary B
    length_B=(int)mpz_sizeinbase(B,2);
    printf("G2 length_B:%d\n",length_B);
    printf("\n");
    
    //set binary
    if(length_A>length_B){
        loop_length=length_A;
    }else{
        loop_length=length_B;
    }
    char binary_A[loop_length+1];
    char binary_B[loop_length+1];
    if(length_A>length_B){
        char binary_buf[length_B+1];
        mpz_get_str(binary_buf,2,B);
        mpz_get_str(binary_A,2,A);
        memset(binary_B,'0',sizeof(binary_B));
        memmove(binary_B+length_A-length_B,binary_buf,sizeof(binary_buf));
        //set next_Q
        EFp2_set(&next_Q,&frobenius_Q);
    }else if(length_A<length_B){
        mpz_get_str(binary_B,2,B);
        char binary_buf[length_A+1];
        mpz_get_str(binary_buf,2,A);
        memset(binary_A,'0',sizeof(binary_A));
        memmove(binary_A+length_B-length_A,binary_buf,sizeof(binary_buf));
        //set next_Q
        EFp2_set(&next_Q,&twisted_Q);
    }else{
        mpz_get_str(binary_A,2,A);
        mpz_get_str(binary_B,2,B);
        //set next_Q
        EFp2_set(&next_Q,&add_R);
    }
    //gettimeofday(&t1,NULL);
    //printf("set time:%.2f[ms]\n",timedifference_msec(t0,t1));
    
    //gettimeofday(&t0,NULL);
    //SCM
    for(i=1; i<loop_length; i++){
        EFp2_ECD(&next_Q,&next_Q);
        if(binary_A[i] == '0' && binary_B[i] == '1'){
            EFp2_ECA(&next_Q,&next_Q,&twisted_Q);
        }else if(binary_A[i] == '1' && binary_B[i] == '0'){
            EFp2_ECA(&next_Q,&next_Q,&frobenius_Q);
        }else if(binary_A[i] == '1' && binary_B[i] == '1'){
            EFp2_ECA(&next_Q,&next_Q,&add_R);
        }
    }
    //gettimeofday(&t1,NULL);
    //printf("set time:%.2f[ms]\n",timedifference_msec(t0,t1));
    
    EFp2_to_EFp12(ANS,&next_Q);
    ANS->flag=next_Q.flag;
    
    mpz_clear(A);
    mpz_clear(B);
    mpz_clear(buf);
    EFp2_clear(&twisted_Q);
    EFp2_clear(&frobenius_Q);
    EFp2_clear(&add_R);
    EFp2_clear(&next_Q);
    EFp12_clear(&Buf);
}

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

void EFp2_to_EFp2(struct EFp2 *ANS,struct EFp2 *P){
    Fp2_inv_basis(&ANS->x,&P->x);
    Fp2_inv_basis(&ANS->y,&ANS->y);
}
