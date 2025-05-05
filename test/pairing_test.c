#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include "../pairing.h"

// Colors for test output
#define GREEN "\033[0;32m"
#define RED "\033[0;31m"
#define BLUE "\033[0;34m"
#define RESET "\033[0m"

// Test status tracking
int tests_run = 0;
int tests_passed = 0;

/**
 * @brief Test helper function to print test result
 */
void test_result(const char *test_name, int passed) {
    tests_run++;
    if (passed) {
        tests_passed++;
        printf("%s[PASS]%s %s\n", GREEN, RESET, test_name);
    } else {
        printf("%s[FAIL]%s %s\n", RED, RESET, test_name);
    }
}

/**
 * @brief Test bilinearity property of pairings
 * 
 * This test verifies that e(aP,Q) = e(P,aQ) = e(P,Q)^a
 */
void test_bilinearity(CurveParams *params) {
    printf("%sTesting pairing bilinearity for %s curve...%s\n", 
           BLUE, 
           params->type == CURVE_BLS12 ? "BLS12" : 
           params->type == CURVE_BN ? "BN" : "KSS16",
           RESET);

    G1Point p, ap;
    G2Point q, aq;
    GTElement e1, e2, e3;
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
    
    // Check bilinearity based on curve type
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
                printf("  - e(aP,Q) or e(P,aQ) does not match e(P,Q)^a\n");
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
                printf("  - e(aP,Q) or e(P,aQ) does not match e(P,Q)^a\n");
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
                printf("  - e(aP,Q) or e(P,aQ) does not match e(P,Q)^a\n");
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
    mpz_clear(a);
    gmp_randclear(state);
    
    test_result("Pairing bilinearity", result);
}

/**
 * @brief Test non-degeneracy property of pairings
 * 
 * This test verifies that e(P,Q) != 1 for non-identity points P and Q
 */
void test_non_degeneracy(CurveParams *params) {
    printf("%sTesting non-degeneracy for %s curve...%s\n", 
           BLUE, 
           params->type == CURVE_BLS12 ? "BLS12" : 
           params->type == CURVE_BN ? "BN" : "KSS16",
           RESET);

    G1Point p;
    G2Point q;
    GTElement e;
    int result = 1;
    
    // Initialize
    g1_init(&p, params);
    g2_init(&q, params);
    gt_init(&e, params);
    
    // Generate random points
    g1_random(&p, params);
    g2_random(&q, params);
    
    // Compute e(P,Q)
    optimal_ate_pairing(&e, &p, &q, params);
    
    // Check non-degeneracy based on curve type
    int is_one = 0;
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *res = (struct Fp12 *)e.element_data;
            is_one = (Fp12_cmp_one(res) == 0);
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *res = (struct Fp12 *)e.element_data;
            is_one = (Fp12_cmp_one(res) == 0);
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *res = (struct Fp16 *)e.element_data;
            is_one = (Fp16_cmp_one(res) == 0);
            break;
        }
    }
    
    if (is_one) {
        result = 0;
        printf("  - e(P,Q) = 1, which violates non-degeneracy\n");
    }
    
    // Clean up
    g1_clear(&p);
    g2_clear(&q);
    gt_clear(&e);
    
    test_result("Pairing non-degeneracy", result);
}

/**
 * @brief Test identity element property: e(O,Q) = e(P,O) = 1
 */
void test_identity_property(CurveParams *params) {
    printf("%sTesting identity property for %s curve...%s\n", 
           BLUE, 
           params->type == CURVE_BLS12 ? "BLS12" : 
           params->type == CURVE_BN ? "BN" : "KSS16",
           RESET);

    G1Point p, p_inf;
    G2Point q, q_inf;
    GTElement e1, e2;
    int result = 1;
    
    // Initialize
    g1_init(&p, params);
    g1_init(&p_inf, params);
    g2_init(&q, params);
    g2_init(&q_inf, params);
    gt_init(&e1, params);
    gt_init(&e2, params);
    
    // Generate random points and set infinity points
    g1_random(&p, params);
    g2_random(&q, params);
    p_inf.infinity = 1;
    q_inf.infinity = 1;
    
    // Test e(O,Q) = 1
    optimal_ate_pairing(&e1, &p_inf, &q, params);
    
    // Test e(P,O) = 1
    optimal_ate_pairing(&e2, &p, &q_inf, params);
    
    // Check results based on curve type
    int is_one_e1 = 0, is_one_e2 = 0;
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *res1 = (struct Fp12 *)e1.element_data;
            struct Fp12 *res2 = (struct Fp12 *)e2.element_data;
            is_one_e1 = (Fp12_cmp_one(res1) == 0);
            is_one_e2 = (Fp12_cmp_one(res2) == 0);
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *res1 = (struct Fp12 *)e1.element_data;
            struct Fp12 *res2 = (struct Fp12 *)e2.element_data;
            is_one_e1 = (Fp12_cmp_one(res1) == 0);
            is_one_e2 = (Fp12_cmp_one(res2) == 0);
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *res1 = (struct Fp16 *)e1.element_data;
            struct Fp16 *res2 = (struct Fp16 *)e2.element_data;
            is_one_e1 = (Fp16_cmp_one(res1) == 0);
            is_one_e2 = (Fp16_cmp_one(res2) == 0);
            break;
        }
    }
    
    if (!is_one_e1) {
        result = 0;
        printf("  - e(O,Q) != 1\n");
    }
    
    if (!is_one_e2) {
        result = 0;
        printf("  - e(P,O) != 1\n");
    }
    
    // Clean up
    g1_clear(&p);
    g1_clear(&p_inf);
    g2_clear(&q);
    g2_clear(&q_inf);
    gt_clear(&e1);
    gt_clear(&e2);
    
    test_result("Pairing identity property", result);
}

/**
 * @brief Compare different pairing implementations for consistency
 * 
 * This test verifies that Tate, Ate, and Optimal Ate pairings return
 * consistent results for the same input points.
 */
void test_pairing_consistency(CurveParams *params) {
    printf("%sTesting pairing consistency for %s curve...%s\n", 
           BLUE, 
           params->type == CURVE_BLS12 ? "BLS12" : 
           params->type == CURVE_BN ? "BN" : "KSS16",
           RESET);

    G1Point p;
    G2Point q;
    GTElement e_tate, e_ate, e_opt_ate;
    int result = 1;
    
    // Initialize
    g1_init(&p, params);
    g2_init(&q, params);
    gt_init(&e_tate, params);
    gt_init(&e_ate, params);
    gt_init(&e_opt_ate, params);
    
    // Generate random points
    g1_random(&p, params);
    g2_random(&q, params);
    
    // Compute pairings
    tate_pairing(&e_tate, &p, &q, params);
    ate_pairing(&e_ate, &p, &q, params);
    optimal_ate_pairing(&e_opt_ate, &p, &q, params);
    
    // Check consistency based on curve type
    int tate_ate_equal = 0, tate_opt_equal = 0, ate_opt_equal = 0;
    switch (params->type) {
        case CURVE_BLS12: {
            struct Fp12 *res_tate = (struct Fp12 *)e_tate.element_data;
            struct Fp12 *res_ate = (struct Fp12 *)e_ate.element_data;
            struct Fp12 *res_opt = (struct Fp12 *)e_opt_ate.element_data;
            
            tate_ate_equal = (Fp12_cmp(res_tate, res_ate) == 0);
            tate_opt_equal = (Fp12_cmp(res_tate, res_opt) == 0);
            ate_opt_equal = (Fp12_cmp(res_ate, res_opt) == 0);
            break;
        }
        
        case CURVE_BN: {
            struct Fp12 *res_tate = (struct Fp12 *)e_tate.element_data;
            struct Fp12 *res_ate = (struct Fp12 *)e_ate.element_data;
            struct Fp12 *res_opt = (struct Fp12 *)e_opt_ate.element_data;
            
            tate_ate_equal = (Fp12_cmp(res_tate, res_ate) == 0);
            tate_opt_equal = (Fp12_cmp(res_tate, res_opt) == 0);
            ate_opt_equal = (Fp12_cmp(res_ate, res_opt) == 0);
            break;
        }
        
        case CURVE_KSS16: {
            struct Fp16 *res_tate = (struct Fp16 *)e_tate.element_data;
            struct Fp16 *res_ate = (struct Fp16 *)e_ate.element_data;
            struct Fp16 *res_opt = (struct Fp16 *)e_opt_ate.element_data;
            
            tate_ate_equal = (Fp16_cmp(res_tate, res_ate) == 0);
            tate_opt_equal = (Fp16_cmp(res_tate, res_opt) == 0);
            ate_opt_equal = (Fp16_cmp(res_ate, res_opt) == 0);
            break;
        }
    }
    
    if (!tate_ate_equal || !tate_opt_equal || !ate_opt_equal) {
        result = 0;
        if (!tate_ate_equal) printf("  - Tate and Ate pairings don't match\n");
        if (!tate_opt_equal) printf("  - Tate and Optimal Ate pairings don't match\n");
        if (!ate_opt_equal) printf("  - Ate and Optimal Ate pairings don't match\n");
    }
    
    // Clean up
    g1_clear(&p);
    g2_clear(&q);
    gt_clear(&e_tate);
    gt_clear(&e_ate);
    gt_clear(&e_opt_ate);
    
    test_result("Pairing consistency", result);
}

/**
 * @brief Run all tests for the specified curve
 */
void test_curve(CurveType curve_type) {
    CurveParams params;
    
    if (pairing_init(&params, curve_type) != 0) {
        printf("%s[ERROR]%s Failed to initialize %s curve\n", RED, RESET,
               curve_type == CURVE_BLS12 ? "BLS12" : 
               curve_type == CURVE_BN ? "BN" : "KSS16");
        return;
    }
    
    printf("\n%s=== Testing %s Curve ===%s\n", BLUE, 
           curve_type == CURVE_BLS12 ? "BLS12" : 
           curve_type == CURVE_BN ? "BN" : "KSS16",
           RESET);
    
    test_bilinearity(&params);
    test_non_degeneracy(&params);
    test_identity_property(&params);
    test_pairing_consistency(&params);
    
    pairing_clear(&params);
}

/**
 * @brief Main function
 */
int main() {
    printf("\n%s=== Elliptic Curve Pairing Tests ===%s\n", BLUE, RESET);
    
    // Test all curve types
    test_curve(CURVE_BLS12);
    test_curve(CURVE_BN);
    test_curve(CURVE_KSS16);
    
    // Print summary
    printf("\n%s=== Test Summary ===%s\n", BLUE, RESET);
    printf("Tests run: %d\n", tests_run);
    printf("Tests passed: %s%d%s\n", (tests_passed == tests_run) ? GREEN : RED, tests_passed, RESET);
    printf("Success rate: %.2f%%\n", (tests_run > 0) ? (float)tests_passed / tests_run * 100 : 0);
    
    return (tests_passed == tests_run) ? 0 : 1;
}