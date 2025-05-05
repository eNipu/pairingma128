#include "pairing.h"
#include "BLS12/BLS_12.h"
#include "BN/BN12_header.h"
#include "KSS16/KSS_16.h"

/**
 * @brief Initialize curve parameters
 * 
 * @param params Pointer to curve parameters structure
 * @param type Type of curve to initialize
 * @return 0 on success, error code otherwise
 */
int pairing_init(CurveParams *params, CurveType type) {
    if (!params) {
        return -1;
    }
    
    mpz_init(params->prime);
    mpz_init(params->order);
    mpz_init(params->trace);
    mpz_init(params->parameter);
    
    params->type = type;
    params->curve_specific = NULL;
    
    switch (type) {
        case CURVE_BLS12:
            // Call BLS12 initialization
            init_parameters();
            set_parameters();
            
            // Copy parameters
            mpz_set(params->prime, prime);
            mpz_set(params->order, EFp_order);
            mpz_set(params->trace, trace_t);
            mpz_set(params->parameter, mother_parameter);
            break;
            
        case CURVE_BN:
            // Similar initialization for BN curves
            init_parameters();
            set_parameters();
            
            // Copy parameters
            mpz_set(params->prime, prime);
            mpz_set(params->order, EFp_order);
            mpz_set(params->trace, trace_t);
            mpz_set(params->parameter, mother_parameter);
            break;
            
        case CURVE_KSS16:
            // KSS16 initialization
            KSS_16_parameters();
            
            // Copy parameters
            mpz_set(params->prime, PRIME_P);
            mpz_set(params->order, order_r);
            mpz_set(params->trace, trace_t);
            mpz_set(params->parameter, X);
            break;
            
        default:
            return -1;
    }
    
    return 0;
}

/**
 * @brief Clean up curve parameters
 * 
 * @param params Pointer to curve parameters structure
 */
void pairing_clear(CurveParams *params) {
    if (!params) return;
    
    mpz_clear(params->prime);
    mpz_clear(params->order);
    mpz_clear(params->trace);
    mpz_clear(params->parameter);
    
    switch (params->type) {
        case CURVE_BLS12:
            clear_parameters();
            break;
            
        case CURVE_BN:
            clear_parameters();
            break;
            
        case CURVE_KSS16:
            dealloc_constants();
            break;
    }
    
    params->curve_specific = NULL;
}

/**
 * @brief Initialize a point in G1
 * 
 * @param point Pointer to G1Point structure
 * @param params Curve parameters
 */
void g1_init(G1Point *point, CurveParams *params) {
    if (!point || !params) return;
    
    point->type = params->type;
    point->infinity = 0;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct EFp12 *p = (struct EFp12 *)malloc(sizeof(struct EFp12));
            if (p) {
                EFp12_init(p);
                point->point_data = p;
            }
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *p = (struct EFp12 *)malloc(sizeof(struct EFp12));
            if (p) {
                EFp12_init(p);
                point->point_data = p;
            }
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp *p = (struct EFp *)malloc(sizeof(struct EFp));
            if (p) {
                EFp_init(p);
                point->point_data = p;
            }
            break;
        }
    }
}

/**
 * @brief Initialize a point in G2
 * 
 * @param point Pointer to G2Point structure
 * @param params Curve parameters
 */
void g2_init(G2Point *point, CurveParams *params) {
    if (!point || !params) return;
    
    point->type = params->type;
    point->infinity = 0;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct EFp12 *p = (struct EFp12 *)malloc(sizeof(struct EFp12));
            if (p) {
                EFp12_init(p);
                point->point_data = p;
            }
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *p = (struct EFp12 *)malloc(sizeof(struct EFp12));
            if (p) {
                EFp12_init(p);
                point->point_data = p;
            }
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp16 *p = (struct EFp16 *)malloc(sizeof(struct EFp16));
            if (p) {
                EFp16_init(p);
                point->point_data = p;
            }
            break;
        }
    }
}

/**
 * @brief Initialize an element in GT
 * 
 * @param element Pointer to GTElement structure
 * @param params Curve parameters
 */
void gt_init(GTElement *element, CurveParams *params) {
    if (!element || !params) return;
    
    element->type = params->type;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *e = (struct Fp12 *)malloc(sizeof(struct Fp12));
            if (e) {
                Fp12_init(e);
                element->element_data = e;
            }
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *e = (struct Fp12 *)malloc(sizeof(struct Fp12));
            if (e) {
                Fp12_init(e);
                element->element_data = e;
            }
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *e = (struct Fp16 *)malloc(sizeof(struct Fp16));
            if (e) {
                Fp16_init(e);
                element->element_data = e;
            }
            break;
        }
    }
}

/**
 * @brief Generate a random point in G1
 * 
 * @param point Pointer to G1Point structure
 * @param params Curve parameters
 */
void g1_random(G1Point *point, CurveParams *params) {
    if (!point || !params || !point->point_data) return;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_generate_G1(p);
            point->infinity = p->flag;
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_generate_G1(p);
            point->infinity = p->flag;
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp *p = (struct EFp *)point->point_data;
            EFp_random_set(p);
            point->infinity = p->infity;
            break;
        }
    }
}

/**
 * @brief Generate a random point in G2
 * 
 * @param point Pointer to G2Point structure
 * @param params Curve parameters
 */
void g2_random(G2Point *point, CurveParams *params) {
    if (!point || !params || !point->point_data) return;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_generate_G2(p);
            point->infinity = p->flag;
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_generate_G2(p);
            point->infinity = p->flag;
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp16 *p = (struct EFp16 *)point->point_data;
            EFp16_random_set_G2(p);
            point->infinity = p->infity;
            break;
        }
    }
}

/**
 * @brief Perform scalar multiplication in G1
 * 
 * @param result Output point P*s
 * @param point Input point P
 * @param scalar Scalar value s
 * @param params Curve parameters
 */
void g1_scalar_mul(G1Point *result, G1Point *point, mpz_t scalar, CurveParams *params) {
    if (!result || !point || !params || !result->point_data || !point->point_data) return;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct EFp12 *res = (struct EFp12 *)result->point_data;
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_SCM(res, p, scalar);
            result->infinity = res->flag;
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *res = (struct EFp12 *)result->point_data;
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_SCM(res, p, scalar);
            result->infinity = res->flag;
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp *res = (struct EFp *)result->point_data;
            struct EFp *p = (struct EFp *)point->point_data;
            EFp_SCM_BIN(res, p, scalar);
            result->infinity = res->infity;
            break;
        }
    }
}

/**
 * @brief Perform scalar multiplication in G2
 * 
 * @param result Output point P*s
 * @param point Input point P
 * @param scalar Scalar value s
 * @param params Curve parameters
 */
void g2_scalar_mul(G2Point *result, G2Point *point, mpz_t scalar, CurveParams *params) {
    if (!result || !point || !params || !result->point_data || !point->point_data) return;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct EFp12 *res = (struct EFp12 *)result->point_data;
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_SCM(res, p, scalar);
            result->infinity = res->flag;
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *res = (struct EFp12 *)result->point_data;
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_SCM(res, p, scalar);
            result->infinity = res->flag;
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp16 *res = (struct EFp16 *)result->point_data;
            struct EFp16 *p = (struct EFp16 *)point->point_data;
            EFp16_SCM_BIN(res, p, scalar);
            result->infinity = res->infity;
            break;
        }
    }
}

/**
 * @brief Compute Tate pairing e(P,Q)
 * 
 * @param result Output GT element
 * @param p G1 point P
 * @param q G2 point Q
 * @param params Curve parameters
 */
void tate_pairing(GTElement *result, G1Point *p, G2Point *q, CurveParams *params) {
    if (!result || !p || !q || !params || !result->element_data || !p->point_data || !q->point_data) return;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *res = (struct Fp12 *)result->element_data;
            struct EFp12 *p_point = (struct EFp12 *)p->point_data;
            struct EFp12 *q_point = (struct EFp12 *)q->point_data;
            Tate_pairing(res, q_point, p_point);
            break;
        }
        
        case CURVE_BN: {
            // BN curves typically use Ate pairing directly instead of Tate
            // This is a placeholder for BN Tate pairing if needed
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *res = (struct Fp16 *)result->element_data;
            struct EFp *p_point = (struct EFp *)p->point_data;
            struct EFp16 *q_point = (struct EFp16 *)q->point_data;
            Tate_Pairing(res, p_point, q_point);
            break;
        }
    }
}

/**
 * @brief Compute Ate pairing e(P,Q)
 * 
 * @param result Output GT element
 * @param p G1 point P
 * @param q G2 point Q
 * @param params Curve parameters
 */
void ate_pairing(GTElement *result, G1Point *p, G2Point *q, CurveParams *params) {
    if (!result || !p || !q || !params || !result->element_data || !p->point_data || !q->point_data) return;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *res = (struct Fp12 *)result->element_data;
            struct EFp12 *p_point = (struct EFp12 *)p->point_data;
            struct EFp12 *q_point = (struct EFp12 *)q->point_data;
            Ate_pairing(res, q_point, p_point);
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *res = (struct Fp12 *)result->element_data;
            struct EFp12 *p_point = (struct EFp12 *)p->point_data;
            struct EFp12 *q_point = (struct EFp12 *)q->point_data;
            // Use BN's Ate pairing function
            Opt_ate_pairing(res, q_point, p_point);
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *res = (struct Fp16 *)result->element_data;
            struct EFp *p_point = (struct EFp *)p->point_data;
            struct EFp16 *q_point = (struct EFp16 *)q->point_data;
            Ate_Pairing(res, p_point, q_point);
            break;
        }
    }
}

/**
 * @brief Compute optimal Ate pairing e(P,Q)
 * 
 * @param result Output GT element
 * @param p G1 point P
 * @param q G2 point Q
 * @param params Curve parameters
 */
void optimal_ate_pairing(GTElement *result, G1Point *p, G2Point *q, CurveParams *params) {
    if (!result || !p || !q || !params || !result->element_data || !p->point_data || !q->point_data) return;
    
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *res = (struct Fp12 *)result->element_data;
            struct EFp12 *p_point = (struct EFp12 *)p->point_data;
            struct EFp12 *q_point = (struct EFp12 *)q->point_data;
            Opt_ate_pairing(res, q_point, p_point);
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *res = (struct Fp12 *)result->element_data;
            struct EFp12 *p_point = (struct EFp12 *)p->point_data;
            struct EFp12 *q_point = (struct EFp12 *)q->point_data;
            Opt_ate_pairing(res, q_point, p_point);
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *res = (struct Fp16 *)result->element_data;
            struct EFp *p_point = (struct EFp *)p->point_data;
            struct EFp16 *q_point = (struct EFp16 *)q->point_data;
            Optimal_Ate_Pairing(res, p_point, q_point);
            break;
        }
    }
}

/**
 * @brief Test pairing bilinearity
 * 
 * This function tests whether e(aP,Q) = e(P,aQ) = e(P,Q)^a
 * 
 * @param params Curve parameters
 * @return 1 if test passes, 0 otherwise
 */
int test_pairing_bilinearity(CurveParams *params) {
    if (!params) return 0;
    
    G1Point p, ap;
    G2Point q, aq;
    GTElement e1, e2, e3, e4;
    mpz_t a;
    int result = 1;
    
    // Initialize
    g1_init(&p, params);
    g1_init(&ap, params);
    g2_init(&q, params);
    g2_init(&aq, params);
    gt_init(&e1, params);
    gt_init(&e2, params);
    gt_init(&e3, params);
    gt_init(&e4, params);
    mpz_init(a);
    
    // Generate random scalar and points
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, (unsigned long)time(NULL));
    mpz_urandomm(a, state, params->order);
    
    // Generate random points
    g1_random(&p, params);
    g2_random(&q, params);
    
    // Compute e(P,Q)
    optimal_ate_pairing(&e1, &p, &q, params);
    
    // Compute e(aP,Q)
    g1_scalar_mul(&ap, &p, a, params);
    optimal_ate_pairing(&e2, &ap, &q, params);
    
    // Compute e(P,aQ)
    g2_scalar_mul(&aq, &q, a, params);
    optimal_ate_pairing(&e3, &p, &aq, params);
    
    // Check bilinearity
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *res1 = (struct Fp12 *)e1.element_data;
            struct Fp12 *res2 = (struct Fp12 *)e2.element_data;
            struct Fp12 *res3 = (struct Fp12 *)e3.element_data;
            struct Fp12 power;
            
            Fp12_init(&power);
            Fp12_pow(&power, res1, a);
            
            if (Fp12_cmp(&power, res2) != 0 || Fp12_cmp(&power, res3) != 0) {
                result = 0;
            }
            
            Fp12_clear(&power);
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *res1 = (struct Fp12 *)e1.element_data;
            struct Fp12 *res2 = (struct Fp12 *)e2.element_data;
            struct Fp12 *res3 = (struct Fp12 *)e3.element_data;
            struct Fp12 power;
            
            Fp12_init(&power);
            Fp12_pow(&power, res1, a);
            
            if (Fp12_cmp(&power, res2) != 0 || Fp12_cmp(&power, res3) != 0) {
                result = 0;
            }
            
            Fp12_clear(&power);
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *res1 = (struct Fp16 *)e1.element_data;
            struct Fp16 *res2 = (struct Fp16 *)e2.element_data;
            struct Fp16 *res3 = (struct Fp16 *)e3.element_data;
            struct Fp16 power;
            
            Fp16_init(&power);
            Fp16_pow(&power, res1, a);
            
            if (Fp16_cmp(&power, res2) != 0 || Fp16_cmp(&power, res3) != 0) {
                result = 0;
            }
            
            Fp16_clear(&power);
            break;
        }
    }
    
    // Clean up
    g1_clear(&p);
    g1_clear(&ap);
    g2_clear(&q);
    g2_clear(&aq);
    gt_clear(&e1);
    gt_clear(&e2);
    gt_clear(&e3);
    gt_clear(&e4);
    mpz_clear(a);
    gmp_randclear(state);
    
    return result;
}

/**
 * @brief Benchmark pairing operations
 * 
 * @param params Curve parameters
 */
void benchmark_pairing(CurveParams *params) {
    if (!params) return;
    
    G1Point p;
    G2Point q;
    GTElement e1, e2, e3;
    struct timeval tv1, tv2;
    
    // Initialize
    g1_init(&p, params);
    g2_init(&q, params);
    gt_init(&e1, params);
    gt_init(&e2, params);
    gt_init(&e3, params);
    
    // Generate random points
    g1_random(&p, params);
    g2_random(&q, params);
    
    // Benchmark Tate pairing
    gettimeofday(&tv1, NULL);
    tate_pairing(&e1, &p, &q, params);
    gettimeofday(&tv2, NULL);
    printf("Tate pairing time: %f ms\n", 
           (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec - tv1.tv_usec) / 1000.0);
    
    // Benchmark Ate pairing
    gettimeofday(&tv1, NULL);
    ate_pairing(&e2, &p, &q, params);
    gettimeofday(&tv2, NULL);
    printf("Ate pairing time: %f ms\n", 
           (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec - tv1.tv_usec) / 1000.0);
    
    // Benchmark Optimal Ate pairing
    gettimeofday(&tv1, NULL);
    optimal_ate_pairing(&e3, &p, &q, params);
    gettimeofday(&tv2, NULL);
    printf("Optimal Ate pairing time: %f ms\n", 
           (tv2.tv_sec - tv1.tv_sec) * 1000.0 + (tv2.tv_usec - tv1.tv_usec) / 1000.0);
    
    // Clean up
    g1_clear(&p);
    g2_clear(&q);
    gt_clear(&e1);
    gt_clear(&e2);
    gt_clear(&e3);
}

/**
 * @brief Clean up a G1 point
 * 
 * @param point Pointer to G1Point structure
 */
void g1_clear(G1Point *point) {
    if (!point || !point->point_data) return;
    
    switch (point->type) {
        case CURVE_BLS12: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_clear(p);
            free(p);
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_clear(p);
            free(p);
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp *p = (struct EFp *)point->point_data;
            EFp_clear(p);
            free(p);
            break;
        }
    }
    
    point->point_data = NULL;
}

/**
 * @brief Clean up a G2 point
 * 
 * @param point Pointer to G2Point structure
 */
void g2_clear(G2Point *point) {
    if (!point || !point->point_data) return;
    
    switch (point->type) {
        case CURVE_BLS12: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_clear(p);
            free(p);
            break;
        }
        
        case CURVE_BN: {
            struct EFp12 *p = (struct EFp12 *)point->point_data;
            EFp12_clear(p);
            free(p);
            break;
        }
        
        case CURVE_KSS16: {
            struct EFp16 *p = (struct EFp16 *)point->point_data;
            EFp16_clear(p);
            free(p);
            break;
        }
    }
    
    point->point_data = NULL;
}

/**
 * @brief Clean up a GT element
 * 
 * @param element Pointer to GTElement structure
 */
void gt_clear(GTElement *element) {
    if (!element || !element->element_data) return;
    
    switch (element->type) {
        case CURVE_BLS12: {
            struct Fp12 *e = (struct Fp12 *)element->element_data;
            Fp12_clear(e);
            free(e);
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *e = (struct Fp12 *)element->element_data;
            Fp12_clear(e);
            free(e);
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *e = (struct Fp16 *)element->element_data;
            Fp16_clear(e);
            free(e);
            break;
        }
    }
    
    element->element_data = NULL;
}