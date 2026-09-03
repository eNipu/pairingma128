#include "pairing.h"
#include "params.h"
#include "util.h"

void Miller_algo(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    struct Fp16 l_sum,v_sum;
    Fp16_init(&l_sum);
    Fp16_init(&v_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    Fp_set_ui(&v_sum.x0.x0.x0.x0,1);
    
    struct EFp16 T;
    EFp16_init(&T);
    EFp16_set(&T,P);
    
    
    struct Fp16 ltt,ltp,v2t,vtp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    Fp16_init(&v2t);
    Fp16_init(&vtp);
    
    int i;
    struct Fp16 tmp1;
    Fp16_init(&tmp1);
    //    Fp16_init(&lambda);
    int r_bit;//bit数
    
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        Fp16_sqr(&l_sum,&l_sum);
        Fp16_sqr(&v_sum,&v_sum);
        
        ltt_q(&ltt,&T,Q);
        Fp16_mul(&l_sum,&l_sum,&ltt);
        
        EFp16_ECD(&T,&T);
        v2t_q(&v2t,&T,Q);
        Fp16_mul(&v_sum,&v_sum,&v2t);
        
        if(mpz_tstbit(loop,i)==1){
            ltp_q(&ltp,&T,P,Q);
            Fp16_mul(&l_sum,&l_sum,&ltp);
            
            EFp16_ECA(&T,&T,P);
            vtp_q(&vtp,&T,Q);
            Fp16_mul(&v_sum,&v_sum,&vtp);
        }
    }
    
    
    // EFp16_printf(&T);
    Fp16_div(ANS,&l_sum,&v_sum);
    
    Fp16_clear(&l_sum);
    Fp16_clear(&v_sum);
    EFp16_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&v2t);
    Fp16_clear(&vtp);
    Fp16_clear(&tmp1);
}

void Optimal_Miller(struct Fp16 *ANS,struct EFp16 *P,struct EFp16 *Q, mpz_t loop){
    
    struct EFp16 T,EFp16_tmp;
    EFp16_init(&T);
    EFp16_init(&EFp16_tmp);
    
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct Fp16 Px_neg;
    Fp16_init(&Px_neg);
    Fp16_neg(&Px_neg,&P->x);//TODO Why neg Px?
    
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    struct EFp16 Q_neg;
    EFp16_init(&Q_neg);
    Fp16_neg(&Q_neg.y,&Q->y);
    Fp16_set(&Q_neg.x,&Q->x);
    
    if(X_bit_binary[x_bit]==-1){
        EFp16_set(&T,&Q_neg);
    }else{
        EFp16_set(&T,Q);
    }
    int i;
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_sqr(&l_sum,&l_sum);
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_sqr(&l_sum,&l_sum);
                
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_sqr(&l_sum,&l_sum);
                
                DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    //  EFp16_SCM_BIN(&EFp_tmp,Q,prime);
    struct EFp4 Q_bar,EFp4_tmp;
    EFp4_init(&Q_bar);
    EFp4_init(&EFp4_tmp);
    EFp16_to_EFp4_map(&Q_bar, Q);
    Skew_Frobenius_map(&EFp4_tmp,&Q_bar);
    EFp4_to_EFp16_map(&EFp16_tmp, &EFp4_tmp);
    //  EFp16_frobenius_map(&EFp_tmp, Q);
    
    
    ltp_q(&ltp,&T,&EFp16_tmp,P);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    ltt_q(&ltt,Q,P);
    
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp16_clear(&T);
    EFp16_clear(&EFp16_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&Q_bar);
    EFp4_clear(&EFp4_tmp);
}

void final_exp_hard(struct Fp16 *ANS,struct Fp16 *fd)
{
    mpz_t Xp1,temp_exp;
    mpz_init(Xp1);
    mpz_init(temp_exp);
    mpz_set(Xp1,X);
    mpz_add_ui(Xp1,Xp1,1);
    
    struct Fp16 t,t0, t1, t2, t3, t4, t5, t6,t7, t8, t9, t10, t11, t12, tmp16, tmp16_1;
    Fp16_init(&t0);
    Fp16_init(&t1);
    Fp16_init(&t2);
    Fp16_init(&t3);
    Fp16_init(&t4);
    Fp16_init(&t5);
    Fp16_init(&t6);
    Fp16_init(&t7);
    Fp16_init(&t8);
    Fp16_init(&t9);
    Fp16_init(&t10);
    Fp16_init(&t11);
    Fp16_init(&t12);
    Fp16_init(&tmp16);
    Fp16_init(&tmp16_1);
    Fp16_init(&t);
    
    
    struct Fp16 t00, t01, t02, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23, t24, t25, t26, t27,t28, t29,t30, t31, t32,t33, t37, s0, s1, s2, s3;
    Fp16_init(&t00);
    Fp16_init(&t01);
    Fp16_init(&t02);
    Fp16_init(&t13);
    Fp16_init(&t14);
    Fp16_init(&t15);
    Fp16_init(&t16);
    Fp16_init(&t17);
    Fp16_init(&t18);
    Fp16_init(&t19);
    Fp16_init(&t20);
    Fp16_init(&t21);
    Fp16_init(&t22);
    Fp16_init(&t23);
    Fp16_init(&t24);
    Fp16_init(&t25);
    Fp16_init(&t26);
    Fp16_init(&t27);
    Fp16_init(&t28);
    Fp16_init(&t29);
    Fp16_init(&t30);
    Fp16_init(&t31);
    Fp16_init(&t32);
    Fp16_init(&t33);
    Fp16_init(&t37);
    Fp16_init(&s0);
    Fp16_init(&s1);
    Fp16_init(&s2);
    Fp16_init(&s3);
    
    
    Fp16_sqr(&t0, fd);
    Fp16_sqr(&t1, &t0);
    Fp16_pow(&t2, fd, Xp1);
    
    Fp16_pow(&t3, &t2, Xp1);
    Fp16_mul(&t4, &t3, &t1);
    
    Fp16_pow(&t5,&t4,X);
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t6, &t4,temp_exp);
    mpz_set_ui(temp_exp,8);
    Fp16_pow(&t7, &t1,temp_exp);
    mpz_set_ui(temp_exp,2);
    Fp16_pow(&t8, &t7, temp_exp);
    
    Fp8_set(&tmp16_1.x0, &t1.x0);
    Fp8_neg(&tmp16_1.x1, &t1.x1);
    
    Fp16_mul(&t9,&t7,&tmp16_1);
    Fp16_sqr(&t10, &t9);
    Fp16_pow(&t11, &t5, X);
    Fp16_pow(&t12, &t11, X);
    Fp16_mul(&t01, &t12, &t10);
    
    //    Fp16_pow(&tmp16, fd, A);
    //    printf("IF t01 eq M^A %d\n",Fp16_cmp(&t01, &tmp16));
    
    Fp16_pow(&t14, &t01, X);
    Fp16_sqr(&tmp16, &t14);
    
    //Fp16_invert(&t13, &tmp16);
    Fp8_set(&t13.x0, &tmp16.x0);
    Fp8_neg(&t13.x1, &tmp16.x1);
    
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t00, &t6, temp_exp);
    Fp16_pow(&t15, &t00, temp_exp);
    
    //Fp16_invert(&tmp16, &t15);
    Fp8_set(&tmp16.x0, &t15.x0);
    Fp8_neg(&tmp16.x1, &t15.x1);
    Fp16_mul(&t0, &t13, &tmp16);
    
    Fp16_sqr(&t16, &t0);
    mpz_set_ui(temp_exp,4);
    Fp16_pow(&t17, &t13, temp_exp);
    Fp16_mul(&t18, &t17, &t14);
    
    //t2:=t16*t18;
    Fp16_mul(&t2, &t16, &t18);
    
    //t19:=t14^(u);
    Fp16_pow(&t19, &t14, X);
    //t20:=t19^(u);
    Fp16_pow(&t20, &t19, X);
    //t21:=t20^(u);
    Fp16_pow(&t21, &t20, X);
    //t22:=t19^2;
    Fp16_sqr(&t22, &t19);
    
    //t23:=t5^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t23, &t5, temp_exp);
    //t24:=t23^5;
    Fp16_pow(&t24, &t23, temp_exp);
    //t25:=t24^3;
    mpz_set_ui(temp_exp,3);
    Fp16_pow(&t25, &t24, temp_exp);
    //t26:=t24*t25;
    Fp16_mul(&t26, &t24, &t25);
    //t27:=t22^2;
    Fp16_sqr(&t27, &t22);
    //t37:=(t27*t25)^(-1);
    Fp16_mul(&tmp16, &t27, &t25);
    
    // Fp16_invert(&t37, &tmp16);
    Fp8_set(&t37.x0, &tmp16.x0);
    Fp8_neg(&t37.x1, &tmp16.x1);
    
    //    mpz_set_str(temp_exp,"88805286319009225549297960754086960732617029553721085722102881418823661400",10);
    //    Fp16_pow(&tmp16, fd, temp_exp);
    //    Fp16_invert(&tmp16, &tmp16);
    //    Fp16_printf(&tmp16);
    //    printf("IF t01 eq M^A %d\n",Fp16_cmp(&t37, &tmp16));
    //    Fp16_printf(&t0);
    
    //t28:=t27*t19^(-1);
    Fp8_set(&tmp16.x0, &t19.x0);
    Fp8_neg(&tmp16.x1, &t19.x1);
    Fp16_mul(&t28, &t27, &tmp16);
    
    //t3:=t28*t26;
    Fp16_mul(&t3, &t28, &t26);
    
    
    //t29:=t11^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&t29, &t11, temp_exp);
    //t30:=t29^2;
    Fp16_sqr(&t30, &t29);
    //t4:=t20*t30;
    Fp16_mul(&t4, &t20, &t30);
    
    //s0:=t20^2;
    Fp16_sqr(&s0, &t20);
    //s1:=t30^5;
    mpz_set_ui(temp_exp,5);
    Fp16_pow(&s1, &t30, temp_exp);
    //s2:=s1*t29;
    Fp16_mul(&s2, &s1, &t29);
    //s3:=s0*s2;
    Fp16_mul(&s3, &s0, &s2);
    
    //t31:=t12^(24);
    mpz_set_ui(temp_exp,24);
    Fp16_pow(&t31, &t12, temp_exp);
    
    //t5:=t21^(-1)*t31^(-1);
    Fp8_set(&tmp16_1.x0, &t31.x0);
    Fp8_neg(&tmp16_1.x1, &t31.x1);
    
    Fp8_set(&tmp16.x0, &t21.x0);
    Fp8_neg(&tmp16.x1, &t21.x1);
    
    Fp16_mul(&t5, &tmp16, &tmp16_1);
    
    //t6:=t8^3*t1;
    mpz_set_ui(temp_exp,3);
    Fp16_pow(&tmp16, &t8, temp_exp);
    Fp16_mul(&t6, &tmp16,&t1);
    
    //t7:=t5*t6;
    Fp16_mul(&t7, &t5, &t6);
    
    //t8:=t01^7;
    mpz_set_ui(temp_exp,7);
    Fp16_pow(&t8, &t01, temp_exp);
    
    //    Fp16_pow(&tmp16, fd, m77);
    //    printf("IF t8 eq M^m77 %d\n",Fp16_cmp(&t8, &tmp16));
    
    //t32:=t37^p*t7^(p^3);
    Fp16_frobenius_map(&tmp16, &t37);
    Fp16_frobenius_map_p3(&tmp16_1, &t7);
    Fp16_mul(&t32, &tmp16_1, &tmp16);
    
    //t32:=t32*t3^(p^5);
    Fp16_frobenius_map_p5(&tmp16, &t3);
    Fp16_mul(&t32, &t32, &tmp16);
    
    //t32:=t32*t8^(p^7);
    Fp16_frobenius_map_p7(&tmp16, &t8);
    Fp16_mul(&t32, &t32, &tmp16);
    
    //t33:=t0^(p^2)*t2^(p^6);
    Fp16_frobenius_map_p2(&tmp16, &t0);
    Fp16_frobenius_map_p6(&tmp16_1, &t2);
    Fp16_mul(&t33, &tmp16, &tmp16_1);
    
    // t7=M^m33 and t2=M^m66.
    
    //t:=t32*t33;
    Fp16_mul(&t, &t33, &t32);
    
    //t:=t*t4^(p^4);
    Fp16_frobenius_map_p4(&tmp16, &t4);
    Fp16_mul(&t, &t, &tmp16);
    
    //t:=t*s3;
    Fp16_mul(&t, &t, &s3);
    
    Fp16_set(ANS, &t);
    
    
    mpz_clears(m00, m11, m22, m33, m44, m55, m66, m77,tmp1, (mpz_ptr) NULL);
    mpz_clears(x2,x3,(mpz_ptr) NULL);
    
    Fp16_clear(&t0);
    Fp16_clear(&t1);
    Fp16_clear(&t2);
    Fp16_clear(&t3);
    Fp16_clear(&t4);
    Fp16_clear(&t5);
    Fp16_clear(&t6);
    Fp16_clear(&t7);
    Fp16_clear(&t8);
    Fp16_clear(&t9);
    Fp16_clear(&t10);
    Fp16_clear(&t11);
    Fp16_clear(&t12);
    Fp16_clear(&tmp16);
    Fp16_clear(&tmp16_1);
    Fp16_clear(&t);
    Fp16_clear(&t00);
    Fp16_clear(&t01);
    Fp16_clear(&t02);
    Fp16_clear(&t13);
    Fp16_clear(&t14);
    Fp16_clear(&t15);
    Fp16_clear(&t16);
    Fp16_clear(&t17);
    Fp16_clear(&t18);
    Fp16_clear(&t19);
    Fp16_clear(&t20);
    Fp16_clear(&t21);
    Fp16_clear(&t22);
    Fp16_clear(&t23);
    Fp16_clear(&t24);
    Fp16_clear(&t25);
    Fp16_clear(&t26);
    Fp16_clear(&t27);
    Fp16_clear(&t28);
    Fp16_clear(&t29);
    Fp16_clear(&t30);
    Fp16_clear(&t31);
    Fp16_clear(&t32);
    Fp16_clear(&t33);
    Fp16_clear(&t37);
    Fp16_clear(&s0);
    Fp16_clear(&s1);
    Fp16_clear(&s2);
    Fp16_clear(&s3);
}

void Final_Exp(struct Fp16 *ANS,struct Fp16 *A){
    
    struct Fp16 M,A_p8;
    Fp16_init(&M);
    Fp16_init(&A_p8);
    
    printf("\n\n f^p^8-1 \n");
    struct timeval t0;
    struct timeval t1;
    gettimeofday(&t0, 0);
    Fp16_frobenius_map_p8(&A_p8,A);
    Fp16_div(&M, &A_p8, A);
    gettimeofday(&t1, 0);
    double elapsed_0 = timedifference_msec(t0, t1);
    printf("FINAL EXP EASY PART: %f [ms]\n", elapsed_0);
    
    
    printf("\n\n M^P^8+1 \n");
    gettimeofday(&t0, 0);
    final_exp_hard(ANS, &M);
    gettimeofday(&t1, 0);
    double elapsed = timedifference_msec(t0, t1);
    printf("FINAL EXP HARD PART ms: %f [ms]\n", elapsed);
    
    Fp16_clear(&M);
    Fp16_clear(&A_p8);
}

void Tate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    Miller_algo(&t_ans,G2,G1,order_r);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
}

void Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,trace_t,1);
    
    Miller_algo(&t_ans,G2,G1,tm1);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
}

void Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp16 *G1,struct EFp16 *G2){
    struct Fp16 Miller_X, t_ans;
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    Optimal_Miller(&Miller_X,G1,G2,X);
    
    Final_Exp(&t_ans,&Miller_X);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}

void ADD_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q,struct Fp16 *Qx_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_sub(&A,&P->x,&T->x);//xt-xp
    Fp16_sub(&B,&P->y,&T->y);//yt-yp
    Fp16_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_add(&D,&T->x,&P->x);
    Fp16_mul(&tmp1,&C,&C);
    Fp16_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp16_mul(&tmp2,&C,&T->x);
    Fp16_sub(&E,&tmp2,&T->y);
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x);
    Fp16_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp16_set(&l_tmp,&Q->y);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Qx_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
}

void DBL_LINE(struct Fp16 *l_ANS,struct EFp16 *T_ANS,struct EFp16 *T,struct EFp16 *P,struct Fp16 *Px_neg){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp16 x,y,tmp;
    Fp16_init(&x);
    Fp16_init(&y);
    Fp16_init(&tmp);
    
    struct EFp16 x3_tmp;
    EFp16_init(&x3_tmp);
    struct Fp16 A,B,C,D,E,F;
    Fp16_init(&A);
    Fp16_init(&B);
    Fp16_init(&C);
    Fp16_init(&D);
    Fp16_init(&E);
    Fp16_init(&F);
    
    
    Fp16_add(&A,&T->y,&T->y);//2y
    Fp16_mul(&B,&T->x,&T->x);//x^2
    Fp16_mul_ui(&B,&B,3);//3x^2
    Fp_add_mpz(&B.x0.x0.x0.x0,&B.x0.x0.x0.x0,a_x);//lambda=3x^2+a
    Fp16_div(&C,&B,&A);//lambda=3x^2+a/2y
    
    Fp16_add(&D,&T->x,&T->x); //D=2x
    Fp16_mul(&tmp1,&C,&C);// lamda^2
    Fp16_sub(&x3_tmp.x,&tmp1,&D); //x3.x=lamda^2-x2t
    
    Fp16_mul(&tmp2,&C,&T->x);//xt.lamda
    Fp16_sub(&E,&tmp2,&T->y);//xt.lamda-yt
    
    Fp16_mul(&tmp3,&C,&x3_tmp.x); //x3.lamda
    Fp16_sub(&x3_tmp.y,&E,&tmp3); // x3.y = xt.lamda-yt - x3.lamda
    
    Fp16_set(&l_tmp,&P->y);
    // Fp_set_ui(&l_tmp.x0.x0.x0,1);
    
    Fp16_add(&l_tmp,&l_tmp,&E);
    
    Fp16_mul(&F,&C,Px_neg);
    Fp16_add(&l_tmp,&l_tmp,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp16_set(T_ANS,&x3_tmp);
    
    if(T->infity==TRUE){
        EFp16_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp16_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp16_set_infity(T_ANS);
        return;
    }
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp16_clear(&x);
    Fp16_clear(&y);
    Fp16_clear(&tmp);
    EFp16_clear(&x3_tmp);
    Fp16_clear(&A);
    Fp16_clear(&B);
    Fp16_clear(&C);
    Fp16_clear(&D);
    Fp16_clear(&E);
    Fp16_clear(&F);
    mpz_clear(cmp);
}

void ltt_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,lambda,ltt;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&lambda);
    Fp16_init(&ltt);
    
    Fp16_mul(&tmp1,&T->x,&T->x);//xt^2
    Fp16_add(&tmp2,&tmp1,&tmp1);
    Fp16_add(&tmp1,&tmp2,&tmp1);//3xt^2
    Fp_add_mpz(&tmp1.x0.x0.x0.x0,&tmp1.x0.x0.x0.x0,a_x);//TODO
    Fp16_add(&tmp2,&T->y,&T->y);//2yt
    
    Fp16_div(&lambda,&tmp1,&tmp2);//lambda=3xt^2+a/2yt
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=xq-xt
    Fp16_mul(&tmp3,&tmp3,&lambda);//tmp3=lambda(xq-xt)
    
    Fp16_sub(&ltt,&Q->y,&T->y);//yq-yt
    Fp16_sub(&ltt,&ltt,&tmp3);//ltt=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltt);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&lambda);
    Fp16_clear(&ltt);
}

void v2t_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 v2t;
    Fp16_init(&v2t);
    
    Fp16_sub(&v2t,&Q->x,&T->x);//v2t=xq-xt
    Fp16_set(ANS,&v2t);
    
    Fp16_clear(&v2t);
}

void ltp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *P,struct EFp16 *Q){
    struct Fp16 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp16_init(&tmp1);
    Fp16_init(&tmp2);
    Fp16_init(&tmp3);
    Fp16_init(&tmp4);
    Fp16_init(&lambda);
    Fp16_init(&ltp);
    
    if((Fp16_cmp(&T->x,&P->x))==0&&(Fp16_cmp(&T->y,&P->y))!=0){//xt==xp&&yt!=yp
        Fp16_sub(&ltp,&Q->x,&T->x);
        Fp16_set(ANS,&ltp);
        
        return;
    }
    
    Fp16_sub(&tmp1,&T->x,&P->x);//xt-xp
    Fp16_sub(&tmp2,&T->y,&P->y);//yt-yp
    Fp16_div(&lambda,&tmp2,&tmp1);//lambda=(yt-tp)/(xt-xp)
    
    Fp16_sub(&tmp3,&Q->x,&T->x);//tmp3=(xq-xt)
    Fp16_mul(&tmp4,&tmp3,&lambda);//tmp4=lambda(xq-xt)
    
    Fp16_sub(&ltp,&Q->y,&T->y);//ltp=yq-yt
    Fp16_sub(&ltp,&ltp,&tmp4);//ltp=yq-yt-lambda(xq-xt)
    
    Fp16_set(ANS,&ltp);
    
    Fp16_clear(&tmp1);
    Fp16_clear(&tmp2);
    Fp16_clear(&tmp3);
    Fp16_clear(&tmp4);
    Fp16_clear(&lambda);
    Fp16_clear(&ltp);
}

void vtp_q(struct Fp16 *ANS,struct EFp16 *T,struct EFp16 *Q){
    struct Fp16 vtp;
    Fp16_init(&vtp);
    if(T->infity==1){//if T is infity
        Fp16_set_ui(ANS,0);
        Fp_set_ui(&ANS->x0.x0.x0.x0,1);
        return;
    }
    
    Fp16_sub(&vtp,&Q->x,&T->x);
    Fp16_set(ANS,&vtp);
    
    Fp16_clear(&vtp);
}

void Pseudo_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp4 *Q,struct EFp *P,struct Fp *L){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
    
    Fp4_sub(&A,&Q->x,&T->x);//xt-xp
    Fp4_sub(&B,&Q->y,&T->y);//yt-yp
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&Q->x);
    Fp4_mul(&tmp1,&C,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul_mpz(&l_tmp.x1.x1,&E,L->x0);
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T,&x3_tmp);
    
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
}

void Pseudo_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T,struct EFp *P,struct Fp *L){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltt;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltt);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    
    mpz_t ac_inv;
    mpz_init(ac_inv);
    
    Fp4_add(&A,&T->y,&T->y);//2yt
    Fp4_mul(&B,&T->x,&T->x);//xt^2
    Fp4_mul_ui(&B,&B,3);//3xt^2
    Fp4_mul_betainv(&tmp);
    Fp4_mul_mpz(&tmp, &tmp, z_inv2_test.x0);
    Fp4_add(&B, &B, &tmp);
    
    Fp4_div(&C,&B,&A);//
    
    Fp4_add(&D,&T->x,&T->x);//D=2xt
    Fp4_mul(&tmp1,&C,&C);//C^2
    Fp4_sub(&x3_tmp.x,&tmp1,&D);//x3.x = C^2-D
    
    Fp4_mul(&tmp2,&C,&T->x); // C*xt
    Fp4_sub(&E,&tmp2,&T->y); //E=C*xt-yt
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);//C*x3.x
    Fp4_sub(&x3_tmp.y,&E,&tmp3);//x3.y = E-C*x3.x
    
    Fp_set_ui(&l_tmp.x0.x0.x0.x0,1);
    
    Fp4_mul_mpz(&l_tmp.x1.x1,&E,L->x0);//l_tmp=
    
    Fp4_neg(&l_tmp.x1.x0,&C);
    Fp16_set(l_ANS,&l_tmp);
    
    EFp4_set(T,&x3_tmp);
    
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltt);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
}

void Pseudo_type1_mul(struct Fp16 *ANS,struct Fp16 *A,struct Fp16 *B){
    
    struct Fp16 a,b;
    Fp16_init(&a);
    Fp16_init(&b);
    
    struct Fp4 a0,a1,a2,a3,b2,b3,C_0,C_1,C_2,C_3,T0,T1,T2,T3,T4,a2a3,tmp2,tmp3;
    Fp4_init(&a0);
    Fp4_init(&a1);
    Fp4_init(&a2);
    Fp4_init(&a3);
    Fp4_init(&b2);
    Fp4_init(&b3);
    
    Fp4_init(&C_0);
    Fp4_init(&C_1);
    Fp4_init(&C_2);
    Fp4_init(&C_3);
    
    Fp4_init(&T0);
    Fp4_init(&T1);
    Fp4_init(&T2);
    Fp4_init(&T3);
    Fp4_init(&T4);
    
    Fp4_init(&a2a3);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    
    
    Fp16_set(&a, A);
    Fp16_set(&b, B);
    
    Fp4_set(&a0, &a.x0.x0);
    Fp4_set(&a1, &a.x0.x1);
    Fp4_set(&a2, &a.x1.x0);
    Fp4_set(&a3, &a.x1.x1);
    
    Fp4_set(&b2, &b.x1.x0);
    Fp4_set(&b3, &b.x1.x1);
    
    
    Fp4_add(&a2a3, &a2, &a3);   /**< (a2+a3) */
    Fp4_add(&T4, &b2, &b3);     /**< t4=(b2+b3) */
    
    Fp4_mul(&T1, &a2, &b2);     /**< t1=(a2*b2) */
    Fp4_mul(&T0, &a3, &b3);     /**< t0=(a3*b3) */
    //    Fp4_mul_v(&T0, &tmp2);/**< t0=(a3*b3)*beta */
    
    Fp4_mul(&tmp2, &a2a3, &T4);
    Fp4_sub(&tmp3, &tmp2, &T1);
    Fp4_sub(&C_0, &tmp3, &T0);
    Fp4_mul_v(&tmp3, &C_0);
    Fp4_add(&C_0, &tmp3, &a0);
    
    Fp4_mul_v(&tmp2, &T0);
    Fp4_add(&C_1, &T1, &tmp2);
    Fp4_add(&C_1, &C_1, &a1);
    
    Fp4_mul(&T3, &a0, &b2);
    Fp4_mul(&T2, &a1, &b3);
    Fp4_mul_v(&tmp2, &T2);
    Fp4_add(&C_2, &T3, &tmp2);
    Fp4_add(&C_2, &C_2, &a2);
    
    
    Fp4_add(&tmp2, &a0, &a1);
    Fp4_mul(&tmp3, &tmp2, &T4);
    Fp4_sub(&tmp2, &tmp3, &T3);
    Fp4_sub(&C_3, &tmp2, &T2);
    Fp4_add(&C_3, &C_3, &a3);
    
    Fp4_set(&a.x0.x0,&C_0);
    Fp4_set(&a.x0.x1,&C_1);
    Fp4_set(&a.x1.x0,&C_2);
    Fp4_set(&a.x1.x1,&C_3);
    
    Fp16_set(ANS, &a);
    
    
    Fp16_clear(&a);
    Fp16_clear(&b);
    Fp4_clear(&a0);
    Fp4_clear(&a1);
    Fp4_clear(&a2);
    Fp4_clear(&a3);
    Fp4_clear(&b2);
    Fp4_clear(&b3);
    Fp4_clear(&C_0);
    Fp4_clear(&C_1);
    Fp4_clear(&C_2);
    Fp4_clear(&C_3);
    Fp4_clear(&T0);
    Fp4_clear(&T1);
    Fp4_clear(&T2);
    Fp4_clear(&T3);
    Fp4_clear(&T4);
}

void Pseudo_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){//Q:G2,P:G1
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,Q_map,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&Q_map);
    EFp4_init(&EFp4_tmp);
    
    struct EFp P_map;
    EFp_init(&P_map);
    
    struct Fp L,xy,xy_2,y_inv,tmp,y_tmp;
    Fp_init(&L);
    Fp_init(&xy);
    Fp_init(&xy_2);
    Fp_init(&y_inv);
    Fp_init(&tmp);
    Fp_init(&y_tmp);
    Fp_init(&z_inv2_test);
    
    mpz_invert(y_inv.x0,P->y.x0.x0.x0,PRIME_P);//yp^-1
    Fp_inv_num++;
    Fp_mul(&xy,&P->x.x0.x0,&y_inv);//xp.yp^-1
    
    Fp_sqr(&xy_2,&xy);//xy_2 = xp^2.yp^-2
    Fp_mul_mpz(&P_map.x,&P->x.x0.x0,xy_2.x0);//P.x= xp^3.yp^-2
    Fp_set(&P_map.y,&P_map.x);
    
    Fp_mul(&y_tmp,&xy_2,&xy);// xp^2.yp^-2 * xp.yp^-1 = xp^3.yp^-3
    Fp4_mul_mpz(&Q_map.y,&Q->y,y_tmp.x0); //Q_map.y = yQ'.xp^3.yp^-3
    Fp4_mul_mpz(&Q_map.x,&Q->x,xy_2.x0); //Q_map.x = xQ'.xp^2.yp^-2
    
    mpz_invert(L.x0,P_map.y.x0,PRIME_P); // L =yp_bar^-1
    Fp_inv_num++;
    Fp_set(&z_inv2_test,&xy_2);
    Fp_sqr(&z_inv2_test,&z_inv2_test);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q_map.y);
    Fp4_set(&Q_neg.x,&Q_map.x);
    
    EFp4_set(&T,&Q_map);
    
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_sqr(&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                break;
            case 1:
                Fp16_sqr(&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&P_map,&L);
                Pseudo_type1_ADD_LINE(&ltp,&T,&Q_map,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_sqr(&l_sum,&l_sum);
                Pseudo_type1_DBL_LINE(&ltt,&T,&P_map,&L);
                Pseudo_type1_ADD_LINE(&ltp,&T,&Q_neg,&P_map,&L);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltt);
                Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    Skew_Frobenius_map(&EFp4_tmp, &Q_map);
    Pseudo_type1_ADD_LINE(&ltp,&T,&EFp4_tmp,&P_map,&L);
    Pseudo_type1_mul(&l_sum,&l_sum,&ltp);
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    Fp16_frobenius_map_p3(&tmp_f, &l_sum);
    
    Pseudo_type1_DBL_LINE(&ltt,&Q_map,&P_map,&L);
    Pseudo_type1_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    EFp4_clear(&Q_neg);
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp_clear(&P_map);
    EFp4_clear(&Q_map);
    EFp4_clear(&EFp4_tmp);
    Fp_clear(&L);
    Fp_clear(&xy);
    Fp_clear(&xy_2);
    Fp_clear(&y_inv);
    Fp_clear(&tmp);
    Fp_clear(&y_tmp);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    Fp16_clear(&tmp_f);
}

void Pseudo_Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct EFp4 G1_EFp4,G2_EFp4;
    EFp4_init(&G1_EFp4);
    EFp4_init(&G2_EFp4);
    
    struct Fp16 ltp,Miller_X,t_ans;
    Fp16_init(&ltp);
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    struct timeval t0;
    struct timeval t1;
    
    EFp4_set_EFp(&G1_EFp4,G1);
    EFp16_to_EFp4_map(&G2_EFp4,G2);
    printf("\n\nPseudo Sparse Miller algo\n");
    
    mpz_mpz_mul=0;
    mpz_ui_mul=0;
    Fp_mpz_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    gettimeofday(&t0, 0);
    Pseudo_type1_Optimal_Miller(&Miller_X,&G1_EFp4,&G2_EFp4,X);
    gettimeofday(&t1, 0);
    double elapsed = timedifference_msec(t0,t1);
    printf("Miller Time : %f [ms]\n", elapsed);
    printf("Miller loop cost\nmpz_mpz_mul:%ld,mpz_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
    
    mpz_mpz_mul=0;
    mpz_ui_mul=0;
    Fp_mpz_sqr=0;
    mpz_mpz_add=0;
    mpz_ui_add=0;
    basis_mul_num=0;
    Fp_inv_num=0;
    Final_Exp(&t_ans,&Miller_X);
    printf("Final exp cost\nmpz_mpz_mul:%ld,mpz_ui_mul:%ld,Fp_sqr:%ld,mpz_add:%ld,mpz_add_ui:%ld,basis_mul:%ld,Fp_inv:%ld\n\n",mpz_mpz_mul,mpz_ui_mul,Fp_mpz_sqr,mpz_mpz_add,mpz_ui_add,basis_mul_num,Fp_inv_num);
    
    Fp16_set(ANS, &t_ans);
    
    Fp16_clear(&ltp);
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}

void Sparse_type1_ADD_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *P,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    Fp4_sub(&A,&P->x,&T->x);//xt-xp
    Fp4_sub(&B,&P->y,&T->y);//yt-yp
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&P->x);
    Fp4_sqr(&tmp1,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
}

void Sparse_type1_DBL_LINE(struct Fp16 *l_ANS,struct EFp4 *T_ANS,struct EFp4 *T,struct EFp4 *Q,struct Fp4 *Qx_neg){
    struct Fp4 tmp1,tmp2,tmp3,tmp4,lambda,ltp;
    Fp4_init(&tmp1);
    Fp4_init(&tmp2);
    Fp4_init(&tmp3);
    Fp4_init(&tmp4);
    Fp4_init(&lambda);
    Fp4_init(&ltp);
    
    struct Fp16 l_tmp;
    Fp16_init(&l_tmp);
    
    struct Fp4 x,y,tmp;
    Fp4_init(&x);
    Fp4_init(&y);
    Fp4_init(&tmp);
    
    struct EFp4 x3_tmp;
    EFp4_init(&x3_tmp);
    struct Fp4 A,B,C,D,E,F;
    Fp4_init(&A);
    Fp4_init(&B);
    Fp4_init(&C);
    Fp4_init(&D);
    Fp4_init(&E);
    Fp4_init(&F);
    
    
    Fp4_add(&A,&T->y,&T->y);//xt-xp
    Fp4_sqr(&B,&T->x);
    Fp4_mul_ui(&B,&B,3);
    struct Fp4 ac_inv;
    Fp4_init(&ac_inv);
    Fp4_mul_betainv(&ac_inv);
    Fp4_add(&B,&B,&ac_inv);
    Fp4_div(&C,&B,&A);//lambda=(yt-tp)/(xt-xp)
    
    Fp4_add(&D,&T->x,&T->x);
    Fp4_sqr(&tmp1,&C);
    Fp4_sub(&x3_tmp.x,&tmp1,&D);
    
    Fp4_mul(&tmp2,&C,&T->x);
    Fp4_sub(&E,&tmp2,&T->y);
    
    Fp4_mul(&tmp3,&C,&x3_tmp.x);
    Fp4_sub(&x3_tmp.y,&E,&tmp3);
    
    Fp4_set(&l_tmp.x0.x0,&Q->y);
    
    Fp4_set(&l_tmp.x1.x1,&E);
    
    Fp4_mul(&F,&C,Qx_neg);
    Fp4_set(&l_tmp.x1.x0,&F);
    
    Fp16_set(l_ANS,&l_tmp);
    EFp4_set(T_ANS,&x3_tmp);
    
    if(T->infity==TRUE){
        EFp4_set(T_ANS,T);
        return;
    }
    mpz_t cmp;
    mpz_init(cmp);
    mpz_set_ui(cmp,0);
    if(Fp4_cmp_mpz(&T->y,cmp)==0){//P.y==0
        EFp4_set_infity(T_ANS);
        return;
    }
    Fp4_clear(&tmp1);
    Fp4_clear(&tmp2);
    Fp4_clear(&tmp3);
    Fp4_clear(&tmp4);
    Fp4_clear(&lambda);
    Fp4_clear(&ltp);
    Fp16_clear(&l_tmp);
    Fp4_clear(&x);
    Fp4_clear(&y);
    Fp4_clear(&tmp);
    EFp4_clear(&x3_tmp);
    Fp4_clear(&A);
    Fp4_clear(&B);
    Fp4_clear(&C);
    Fp4_clear(&D);
    Fp4_clear(&E);
    Fp4_clear(&F);
    mpz_clear(cmp);
}

void Sparse_type1_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,EFp_tmp;
    EFp4_init(&T);
    EFp4_init(&EFp_tmp);
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    Fp4_neg(&Px_neg,&P->x);
    
    EFp4_set(&T,Q);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    int r_bit;
    r_bit= (int)mpz_sizeinbase(loop,2);
    
    for(i=r_bit-2;i>=0;i--){
        if(mpz_tstbit(loop,i)==1){
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
            Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
            Fp16_mul(&l_sum,&l_sum,&ltp);
        }else{
            Fp16_mul(&l_sum,&l_sum,&l_sum);
            Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
            Fp16_mul(&l_sum,&l_sum,&ltt);
        }
    }
    Fp16_set(ANS,&l_sum);
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    EFp4_clear(&EFp_tmp);
    Fp4_clear(&Px_neg);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
}

void Sparse_type1_Optimal_Miller(struct Fp16 *ANS,struct EFp4 *P,struct EFp4 *Q,mpz_t loop){
    struct Fp16 l_sum;
    Fp16_init(&l_sum);
    Fp_set_ui(&l_sum.x0.x0.x0.x0,1);
    
    struct EFp4 T,EFp4_tmp;
    EFp4_init(&T);
    EFp4_init(&EFp4_tmp);
    
    mpz_t p3;
    mpz_init(p3);
    
    struct Fp16 ltt,ltp;
    Fp16_init(&ltt);
    Fp16_init(&ltp);
    
    int i;
    
    struct Fp4 Px_neg;
    Fp4_init(&Px_neg);
    
    Fp4_neg(&Px_neg,&P->x);
    
    struct EFp4 Q_neg;
    EFp4_init(&Q_neg);
    Fp4_neg(&Q_neg.y,&Q->y);
    Fp4_set(&Q_neg.x,&Q->x);
    
    if(X_bit_binary[x_bit]==-1){
        EFp4_set(&T,&Q_neg);
    }else{
        EFp4_set(&T,Q);
    }
    for(i=x_bit-1;i>=0;i--){
        switch (X_bit_binary[i]){
            case 0:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Fp16_mul(&l_sum,&l_sum,&ltt);
                break;
                
            case 1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Sparse_type1_ADD_LINE(&ltp,&T,&T,Q,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
            case -1:
                Fp16_mul(&l_sum,&l_sum,&l_sum);
                
                Sparse_type1_DBL_LINE(&ltt,&T,&T,P,&Px_neg);
                Sparse_type1_ADD_LINE(&ltp,&T,&T,&Q_neg,P,&Px_neg);
                
                Fp16_mul(&l_sum,&l_sum,&ltt);
                Fp16_mul(&l_sum,&l_sum,&ltp);
                break;
        }
    }
    
    Skew_Frobenius_map(&EFp4_tmp,Q);
    //    EFp4_SCM_BIN(&EFp4_tmp, &Q_map, prime);
    
    Sparse_type1_ADD_LINE(&ltp,&T,&T,&EFp4_tmp,P,&Px_neg);
    Fp16_mul(&l_sum,&l_sum,&ltp);
    
    
    struct Fp16 tmp_f;
    Fp16_init(&tmp_f);
    
    
    Fp16_frobenius_map(&tmp_f, &l_sum);
    //    Fp16_pow(&l_sum, &tmp_f,prime);
    Fp16_frobenius_map(&l_sum, &tmp_f);
    //    Fp16_pow(&tmp_f, &l_sum,prime);
    Fp16_frobenius_map(&tmp_f, &l_sum);
    
    
    //ltt_q(&ltt,Q,P);
    Sparse_type1_DBL_LINE(&ltt,&T,Q,P,&Px_neg);
    
    Fp16_mul(&l_sum,&tmp_f,&ltt);
    Fp16_set(ANS,&l_sum);
    
    
    Fp16_clear(&l_sum);
    EFp4_clear(&T);
    Fp16_clear(&ltt);
    Fp16_clear(&ltp);
    EFp4_clear(&EFp4_tmp);
    
}

void Sparse_Ate_Pairing(struct Fp16 *ANS,struct EFp4 *G1,struct EFp4 *G2){
    struct Fp16 t_ans;
    Fp16_init(&t_ans);
    
    mpz_t tm1;
    mpz_init(tm1);
    mpz_sub_ui(tm1,trace_t,1);
    
    Sparse_type1_Miller(&t_ans,G1,G2,tm1);
    Final_Exp(&t_ans,&t_ans);
    Fp16_set(ANS,&t_ans);
    
    Fp16_clear(&t_ans);
    mpz_clear(tm1);
}

void Sparse_Optimal_Ate_Pairing(struct Fp16 *ANS,struct EFp *G1,struct EFp16 *G2){
    struct EFp4 G1_EFp4,G2_EFp4;
    EFp4_init(&G1_EFp4);
    EFp4_init(&G2_EFp4);
    
    struct Fp16 ltp,Miller_X,t_ans;
    Fp16_init(&ltp);
    Fp16_init(&Miller_X);
    Fp16_init(&t_ans);
    
    EFp4_set_EFp(&G1_EFp4,G1);
    EFp16_to_EFp4_map(&G2_EFp4,G2);
    
    Sparse_type1_Optimal_Miller(&Miller_X,&G1_EFp4,&G2_EFp4,X);
    Final_Exp(&t_ans,&Miller_X);
    Fp16_set(ANS, &t_ans);
    
    Fp16_clear(&ltp);
    Fp16_clear(&Miller_X);
    Fp16_clear(&t_ans);
}
