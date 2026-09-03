#include "test.h"
#include "params.h"
#include "util.h"

void test_tate_pairing(){
    printf("====================================================================================\n");
    printf("Tate pairing\n\n");
}

void test_ate_pairing(){
    printf("====================================================================================\n");
    printf("Ate pairing\n\n");
}

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
