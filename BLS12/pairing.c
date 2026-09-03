#include "pairing.h"
#include "params.h"
#include "util.h"

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

void ff_ltt_vtt(struct Fp12 *f,struct EFp12 *T,struct EFp12 *Q){
    
}

void f_ltp_vtp(struct Fp12 *f,struct EFp12 *T,struct EFp12 *P,struct EFp12 *Q){
    
}

void Miller_algo_for_tate(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P){
    
}

void Tate_pairing(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P){
    
}

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
