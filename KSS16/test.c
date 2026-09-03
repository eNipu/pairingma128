#include "test.h"
#include "params.h"
#include "util.h"

void check_Pairing(void){
    struct EFp P_EFp, R_EFp;
    EFp_init(&P_EFp);
    EFp_init(&R_EFp);
    
    
    struct EFp16 P_Fp16,Q_Fp16,Q_G2,R_Fp16,S_Fp16;
    EFp16_init(&P_Fp16);
    EFp16_init(&Q_Fp16);
    EFp16_init(&R_Fp16);
    EFp16_init(&S_Fp16);
    EFp16_init(&Q_G2);
    
    struct Fp16 ans_Fp16,tmp1_Fp16,tmp2_Fp16,tmp3_Fp16;
    Fp16_init(&ans_Fp16);
    Fp16_init(&tmp1_Fp16);
    Fp16_init(&tmp2_Fp16);
    Fp16_init(&tmp3_Fp16);
    
    struct EFp4 P_EFp4, Q_EFp4, R_EFp4, S_EFp4;
    EFp4_init(&P_EFp4);
    EFp4_init(&Q_EFp4);
    EFp4_init(&R_EFp4);
    EFp4_init(&S_EFp4);
    
    mpz_t a,b,ab;
    mpz_init(a);
    mpz_init(b);
    mpz_init(ab);
    
    mpz_set_ui(a,31);
    mpz_set_ui(b,11);
    mpz_mul(ab,a,b);
    
    mpz_t x,y;
    mpz_inits(x,y,(mpz_ptr) NULL);
    mpz_set_str(x,"585634432000126707887057629201458798521445673169852167601690963969803976546284829414253996667593644070", 10);
    mpz_set_str(y,"354142610898165644039571248952651466071208533406755516950532183093299777682028181848353083070437124285", 10);
    Fp_set_mpz(&P_EFp.x,x);
    Fp_set_mpz(&P_EFp.y,y);
    EFp16_set_EFp(&P_Fp16,&P_EFp);
    mpz_clears(x,y,(mpz_ptr) NULL);
    
    mpz_t x1, x2, x3, x4, y1, y2, y3, y4;
    mpz_inits(x1,x2,x3,x4,y1,y2,y3,y4,(mpz_ptr) NULL);
    
    mpz_set_str(x1,"244491741495785299367612042330502507448971522000836606635913541960426621971155507113603357863094662406", 10);
    mpz_set_str(x2,"218139782565182912336439776658495825407576298090161716969444130503315739115994708740875588793005720529", 10);
    mpz_set_str(x3,"489284606975189850057489668560782153703917092806842888347849447445917065112045634272024913067442383327", 10);
    mpz_set_str(x4,"246161136867400494929650830056723287620786052741538082500453356439354042813046547782023695596667138627", 10);
    mpz_set_str(y1,"48629828925217256502074899893821067838880538260652902679139733711894991914937403050096788961116072841", 10);
    mpz_set_str(y2,"338343177397929794197126001730091007944155731601968256694856657191699538746843500939799293603025576582", 10);
    mpz_set_str(y3,"38952231663050097715674495759233974069804273632029475897530255418880509349197207272848037344062136227", 10);
    mpz_set_str(y4,"250604271284793810746824678011335939079727656816914850518650242122895946469531880085329357599729742073", 10);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x0.x0,x1);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x0.x1,x2);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x1.x0,x3);
    Fp_set_mpz(&Q_Fp16.x.x0.x1.x1.x1,x4);
    
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x0.x0,y1);
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x0.x1,y2);
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x1.x0,y3);
    Fp_set_mpz(&Q_Fp16.y.x1.x1.x1.x1,y4);
    
    printf("\n\nPseudo Sparse Optimal Ate Pairing \n");
    Big_M = 0; Small_m =0, Big_add=0, Sqr=0, fp16_pow =0;
    Pseudo_Sparse_Optimal_Ate_Pairing(&tmp1_Fp16,&P_EFp,&Q_Fp16);
    printf("Fp SQR in Fp16 mul = Fp16 pow %d, %d, %d, %d, %d\n",Big_M, Sqr, Small_m,Big_add,fp16_pow);
    Fp16_pow(&tmp1_Fp16,&tmp1_Fp16,ab);
 
    
    EFp_SCM_BIN(&R_EFp,&P_EFp,a);
    EFp16_SCM_BIN(&S_Fp16,&Q_Fp16,b);
    Pseudo_Sparse_Optimal_Ate_Pairing(&tmp2_Fp16,&R_EFp,&S_Fp16);

    printf("Bilinearity Check: \n");
    if (Fp16_cmp(&tmp2_Fp16, &tmp1_Fp16) == 0) {
        printf("Success \n");
    }
    else{
        printf("Failure\n");
    }
    
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(ab);
    
    EFp_clear(&P_EFp);
    EFp_clear(&R_EFp);
    
    EFp16_clear(&P_Fp16);
    EFp16_clear(&Q_Fp16);
    EFp16_clear(&R_Fp16);
    EFp16_clear(&S_Fp16);
    
    Fp16_clear(&ans_Fp16);
    Fp16_clear(&tmp1_Fp16);
    Fp16_clear(&tmp2_Fp16);
    Fp16_clear(&tmp3_Fp16);
    
    EFp4_clear(&P_EFp4);
    EFp4_clear(&Q_EFp4);
    EFp4_clear(&R_EFp4);
    EFp4_clear(&S_EFp4);
}
