#include "pairing.h"
#include "params.h"
#include "util.h"

void Final_exp_1(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp12 Tmp,t0,t1,t2,t3,t4;
    Fp12_init(&Tmp);
    Fp12_init(&t0);
    Fp12_init(&t1);
    Fp12_init(&t2);
    Fp12_init(&t3);
    Fp12_init(&t4);
    
    Fp12_set(&Tmp,A);
    
    //f←f^(p^6)*f^-1
    Fp12_frobenius_6(&t0,&Tmp);//f^(p^6)
    Fp12_inv(&t1,&Tmp);//f^-1
    Fp12_mul(&Tmp,&t0,&t1);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_2(&t0,&Tmp);//f^(p^2)
    Fp12_mul(&Tmp,&t0,&Tmp);//f^(p^2)*f
    
    Fp12_pow(&t0,&Tmp,mother_parameter);	//t0←f^(-u)
    Fp12_frobenius_6(&t0,&t0);
    Fp12_squaring(&t0,&t0);				//t0←t0^2
    Fp12_squaring(&t1,&t0);				//t1←t0^2
    Fp12_mul(&t1,&t0,&t1);				//t1←t0*t1
    Fp12_pow(&t2,&t1,mother_parameter);		//t2←t1^(-u)
    Fp12_frobenius_6(&t2,&t2);
    Fp12_frobenius_6(&t3,&t1);			//t3←t1^-1
    Fp12_mul(&t1,&t2,&t3);				//t1←t2*t3
    Fp12_squaring(&t3,&t2);				//t3←t2^2
    Fp12_pow(&t4,&t3,mother_parameter);		//t4←t3^(-u)
    Fp12_frobenius_6(&t4,&t4);
    Fp12_frobenius_6(&t4,&t4);			//t4←t4^(-1)
    Fp12_mul(&t4,&t4,&t1);				//t4←t4*t1
    Fp12_mul(&t3,&t4,&t0);				//t3←t4*t0
    Fp12_mul(&t0,&t2,&t4);				//t0←t2*t4
    Fp12_mul(&t0,&t0,&Tmp);				//t0←t0*f
    Fp12_frobenius_1(&t2,&t3);			//t2←t3^p
    Fp12_mul(&t0,&t2,&t0);				//t0←t2*t0
    Fp12_frobenius_2(&t2,&t4);			//t2←t4^(p^2)
    Fp12_mul(&t0,&t2,&t0);				//t0←t2*t0
    Fp12_frobenius_6(&t2,&Tmp);			//t2←f^(-1)
    Fp12_mul(&t2,&t2,&t3);				//t2←t2*t3
    Fp12_frobenius_3(&t2,&t2);			//t2←t2^(p^3)
    Fp12_mul(ANS,&t2,&t0);				//t0←t2*t0
    
    Fp12_clear(&Tmp);
    Fp12_clear(&t0);
    Fp12_clear(&t1);
    Fp12_clear(&t2);
    Fp12_clear(&t3);
    Fp12_clear(&t4);
}

void Final_exp_2(struct Fp12 *ANS,struct Fp12 *A){
    struct Fp12 Tmp,Buf1,Buf2,Buf3,a,b,c;
    Fp12_init(&Tmp);
    Fp12_set(&Tmp,A);
    Fp12_init(&Buf1);
    Fp12_init(&Buf2);
    Fp12_init(&Buf3);
    Fp12_init(&a);
    Fp12_init(&b);
    Fp12_init(&c);
    mpz_t exp;
    mpz_init(exp);
    //f←f^(p^6)*f^-1
    Fp12_frobenius_6(&Buf1,&Tmp);//f^(p^6)
    Fp12_inv(&Buf2,&Tmp);//f^-1
    Fp12_mul(&Tmp,&Buf1,&Buf2);//f^(p^6)*f^-1
    
    //f←f^(p^2)*f
    Fp12_frobenius_2(&Buf1,&Tmp);//f^(p^2)
    Fp12_mul(&Tmp,&Buf1,&Tmp);//f^(p^2)*f
    //Fp12_pow(ANS,&Tmp,final_exp);
    
    //a←(f^(f^6))*(6x+5)
    Fp12_frobenius_6(&Buf1,&Tmp);//f^(p^6)
    mpz_mul_ui(exp,mother_parameter,6);
    mpz_add_ui(exp,exp,5);
    Fp12_pow(&a,&Buf1,exp);
    
    //b←a^p
    Fp12_pow(&b,&a,prime);
    //b←a*b
    Fp12_mul(&b,&a,&b);
    
    //compute f^p,f^p^2,f^p^3
    Fp12_frobenius_1(&Buf1,&Tmp);//f^p
    Fp12_frobenius_2(&Buf2,&Tmp);//f^p^2
    Fp12_frobenius_3(&Buf3,&Tmp);//f^p^3
    
    Fp12_squaring(&c,&Buf1);//b*(f^p)^2*f^p^2
    Fp12_mul(&c,&c,&b);
    Fp12_mul(&c,&c,&Buf2);
    
    Fp12_mul(&Buf1,&Buf1,&Tmp);//f^p^3*(c^6)^x^2*c*b*(f^p*f)^9*a*f^4
    Fp12_squaring(&Buf2,&Buf1);//(f^p*f)^9
    Fp12_squaring(&Buf2,&Buf2);
    Fp12_squaring(&Buf2,&Buf2);
    Fp12_mul(&Buf1,&Buf1,&Buf2);
    Fp12_mul(&Buf1,&Buf1,&a);
    Fp12_mul(&Buf1,&Buf1,&b);
    Fp12_mul(&Buf1,&Buf1,&c);//c*b*(f^p*f)^9*a
    
    Fp12_squaring(&Buf2,&c);
    Fp12_squaring(&Buf2,&Buf2);
    Fp12_mul(&Buf2,&Buf2,&c);
    Fp12_mul(&c,&Buf2,&c);//(c^6)
    mpz_mul(exp,mother_parameter,mother_parameter);
    Fp12_pow(&c,&c,exp);//(c^6)^x^2
    Fp12_mul(&Buf1,&Buf1,&c);//(c^6)^x^2*c*b*(f^p*f)^9*a
    
    Fp12_squaring(&Buf2,&Tmp);
    Fp12_squaring(&Buf2,&Buf2);//f^4
    Fp12_mul(&Buf1,&Buf1,&Buf2);//(c^6)^x^2*c*b*(f^p*f)^9*a*f^4
    
    Fp12_mul(ANS,&Buf1,&Buf3);//f^p^3*(c^6)^x^2*c*b*(f^p*f)^9*a*f^4
    
    mpz_clear(exp);
    Fp12_clear(&a);
    Fp12_clear(&b);
    Fp12_clear(&c);
    Fp12_clear(&Tmp);
    Fp12_clear(&Buf1);
    Fp12_clear(&Buf2);
    Fp12_clear(&Buf3);
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

void EFp_ECD_return_lambda(struct EFp *ANS,struct Fp *lambda,struct EFp *P){
    if(Fp_cmp_zero(&P->y)==0){
        ANS->flag=1;
        return;
    }
    
    struct EFp Tmp;
    EFp_init(&Tmp);
    EFp_set(&Tmp,P);
    struct Fp Buf1,Buf2;
    Fp_init(&Buf1);
    Fp_init(&Buf2);
    
    Fp_mul_ui(&Buf1,&Tmp.y,2);
    Fp_inv(&Buf1,&Buf1);
    Fp_mul(&Buf2,&Tmp.x,&Tmp.x);
    Fp_mul_ui(&Buf2,&Buf2,3);
    Fp_add_mpz(&Buf2,&Buf2,curve_parameter_A);
    Fp_mul(lambda,&Buf1,&Buf2);
    Fp_mul(&Buf1,lambda,lambda);
    Fp_mul_ui(&Buf2,&Tmp.x,2);
    Fp_sub(&ANS->x,&Buf1,&Buf2);
    Fp_sub(&Buf1,&Tmp.x,&ANS->x);
    Fp_mul(&Buf2,lambda,&Buf1);
    Fp_sub(&ANS->y,&Buf2,&Tmp.y);
    
    //clear
    Fp_clear(&Buf1);
    Fp_clear(&Buf2);
    EFp_clear(&Tmp);
}

void EFp_ECA_return_lambda(struct EFp *ANS,struct Fp *lambda,struct EFp *P1,struct EFp *P2){
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
    struct Fp Buf1,Buf2;
    Fp_init(&Buf1);
    Fp_init(&Buf2);
    
    Fp_sub(&Buf1,&Tmp2.x,&Tmp1.x);
    Fp_inv(&Buf1,&Buf1);
    Fp_sub(&Buf2,&Tmp2.y,&Tmp1.y);
    Fp_mul(lambda,&Buf1,&Buf2);
    Fp_mul(&Buf1,lambda,lambda);
    Fp_sub(&Buf2,&Buf1,&Tmp1.x);
    Fp_sub(&ANS->x,&Buf2,&Tmp2.x);
    Fp_sub(&Buf1,&Tmp1.x,&ANS->x);
    Fp_mul(&Buf2,lambda,&Buf1);
    Fp_sub(&ANS->y,&Buf2,&Tmp1.y);
    
    //clear
    Fp_clear(&Buf1);
    Fp_clear(&Buf2);
    EFp_clear(&Tmp1);
    EFp_clear(&Tmp2);
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

void Pseudo_8_sparse_mul_2(struct Fp12 *ANS,struct Fp12 *A,struct Fp12 *B){
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

void ff_ltt_2(struct Fp12 *f,struct EFp2 *T,struct EFp *P,struct Fp *L){
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
    
    Pseudo_8_sparse_mul_2(f,&ff,&ltt);
    
    EFp2_clear(&Tmp_T);
    Fp2_clear(&A);
    Fp2_clear(&B);
    Fp2_clear(&C);
    Fp2_clear(&D);
    Fp2_clear(&E);
    Fp12_clear(&ff);
    Fp12_clear(&ltt);
}

void f_ltq_2(struct Fp12 *f,struct EFp2 *T,struct EFp2 *Q,struct EFp *P,struct Fp *L){
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
    
    Pseudo_8_sparse_mul_2(f,f,&ltq);
    
    EFp2_clear(&Tmp_T);
    Fp12_clear(&ltq);
    Fp2_clear(&A);
    Fp2_clear(&B);
    Fp2_clear(&C);
    Fp2_clear(&D);
}

void Miller_algo_for_opt_ate_2(struct Fp12 *ANS,struct EFp12 *Q,struct EFp12 *P){
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
    
    EFp2_set(&mapped_Q_neg,&mapped_Q);//set mapped_Q_neg
    Fp2_set_neg(&mapped_Q_neg.y,&mapped_Q_neg.y);
    EFp2_set(&T,&mapped_Q);		//set T
    Fp12_set_ui(&f,0);			//set f
    Fp_set_ui(&f.x0.x0.x0,1);
    //miller
    for(i=x_bit_for_opt_ate-1; i>=0; i--){
        switch(X_bit_binary_for_opt_ate[i]){
            case 0:
                ff_ltt_2(&f,&T,&mapped_P,&L);
                break;
            case 1:
                ff_ltt_2(&f,&T,&mapped_P,&L);
                f_ltq_2(&f,&T,&mapped_Q,&mapped_P,&L);
                break;
            case -1:
                ff_ltt_2(&f,&T,&mapped_P,&L);
                f_ltq_2(&f,&T,&mapped_Q_neg,&mapped_P,&L);
                break;
            default:
                break;
        }
        
    }
    
    EFp2_skew_frobenius_1(&mapped_Q1,&mapped_Q);//Q^p
    EFp2_skew_frobenius_2(&mapped_Q2_neg,&mapped_Q);//Q^(p^2)
    EFp2_set_neg(&mapped_Q2_neg,&mapped_Q2_neg);
    f_ltq_2(&f,&T,&mapped_Q1,&mapped_P,&L);
    f_ltq_2(&f,&T,&mapped_Q2_neg,&mapped_P,&L);
    
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
    Fp_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    gettimeofday(&t0,NULL);
    Miller_algo_for_opt_ate_2(ANS,Q,P);
    gettimeofday(&t1,NULL);
    //printf("<CASE 2>\n1+(α+1)^-1*a'[βγ]+(α+1)^-1*b'[β^2γ]\n");
    printf("loop time 2 :%.2f[ms]\n",timedifference_msec(t0,t1));
    printf("Miller loop cost\nmpz_mpz_mul:%ld,mpz_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
    //Fp12_printf(ANS,"test2:"); printf("\n");
    
    mpz_mpz_mul=0;
    mpz_ui_mul=0;
    Fp_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    gettimeofday(&t0,NULL);
    Final_exp_1(ANS,ANS);
    gettimeofday(&t1,NULL);
    printf("final exp time:%.2f[ms]\n",timedifference_msec(t0,t1));
    printf("Final exp cost\nmpz_mpz_mul:%ld,mpz_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
}
