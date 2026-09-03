#include "curve.h"
#include "params.h"
#include "util.h"

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
