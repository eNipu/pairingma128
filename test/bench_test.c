#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include "../pairing.h"

// Colors for output
#define GREEN "\033[0;32m"
#define YELLOW "\033[0;33m"
#define BLUE "\033[0;34m"
#define RESET "\033[0m"

// Number of iterations for benchmarks
#define ITERATIONS 10

/**
 * @brief Calculate time difference in milliseconds
 */
float timedifference_msec(struct timeval start, struct timeval end) {
    return (end.tv_sec - start.tv_sec) * 1000.0f + (end.tv_usec - start.tv_usec) / 1000.0f;
}

/**
 * @brief Benchmark G1 scalar multiplication operations
 */
void benchmark_g1_scalar_mul(CurveParams *params) {
    G1Point p, result;
    mpz_t scalar;
    struct timeval start, end;
    float total_time = 0.0;
    
    // Initialize
    g1_init(&p, params);
    g1_init(&result, params);
    mpz_init(scalar);
    
    printf("Benchmarking G1 scalar multiplication...\n");
    
    // Generate random point and scalar
    g1_random(&p, params);
    mpz_urandomm(scalar, NULL, params->order);
    
    for (int i = 0; i < ITERATIONS; i++) {
        // Time the operation
        gettimeofday(&start, NULL);
        g1_scalar_mul(&result, &p, scalar, params);
        gettimeofday(&end, NULL);
        
        float time_ms = timedifference_msec(start, end);
        total_time += time_ms;
        printf("  Iteration %d: %f ms\n", i + 1, time_ms);
    }
    
    printf("%sAverage G1 scalar multiplication time: %f ms%s\n\n", 
           GREEN, total_time / ITERATIONS, RESET);
    
    // Clean up
    g1_clear(&p);
    g1_clear(&result);
    mpz_clear(scalar);
}

/**
 * @brief Benchmark G2 scalar multiplication operations
 */
void benchmark_g2_scalar_mul(CurveParams *params) {
    G2Point q, result;
    mpz_t scalar;
    struct timeval start, end;
    float total_time = 0.0;
    
    // Initialize
    g2_init(&q, params);
    g2_init(&result, params);
    mpz_init(scalar);
    
    printf("Benchmarking G2 scalar multiplication...\n");
    
    // Generate random point and scalar
    g2_random(&q, params);
    mpz_urandomm(scalar, NULL, params->order);
    
    for (int i = 0; i < ITERATIONS; i++) {
        // Time the operation
        gettimeofday(&start, NULL);
        g2_scalar_mul(&result, &q, scalar, params);
        gettimeofday(&end, NULL);
        
        float time_ms = timedifference_msec(start, end);
        total_time += time_ms;
        printf("  Iteration %d: %f ms\n", i + 1, time_ms);
    }
    
    printf("%sAverage G2 scalar multiplication time: %f ms%s\n\n", 
           GREEN, total_time / ITERATIONS, RESET);
    
    // Clean up
    g2_clear(&q);
    g2_clear(&result);
    mpz_clear(scalar);
}

/**
 * @brief Benchmark pairing operations
 */
void benchmark_pairings(CurveParams *params) {
    G1Point p;
    G2Point q;
    GTElement e_tate, e_ate, e_opt_ate;
    struct timeval start, end;
    float tate_total = 0.0, ate_total = 0.0, opt_ate_total = 0.0;
    
    // Initialize
    g1_init(&p, params);
    g2_init(&q, params);
    gt_init(&e_tate, params);
    gt_init(&e_ate, params);
    gt_init(&e_opt_ate, params);
    
    // Generate random points
    g1_random(&p, params);
    g2_random(&q, params);
    
    // Benchmark Tate pairing
    printf("Benchmarking Tate pairing...\n");
    for (int i = 0; i < ITERATIONS; i++) {
        gettimeofday(&start, NULL);
        tate_pairing(&e_tate, &p, &q, params);
        gettimeofday(&end, NULL);
        
        float time_ms = timedifference_msec(start, end);
        tate_total += time_ms;
        printf("  Iteration %d: %f ms\n", i + 1, time_ms);
    }
    
    printf("%sAverage Tate pairing time: %f ms%s\n\n", 
           GREEN, tate_total / ITERATIONS, RESET);
    
    // Benchmark Ate pairing
    printf("Benchmarking Ate pairing...\n");
    for (int i = 0; i < ITERATIONS; i++) {
        gettimeofday(&start, NULL);
        ate_pairing(&e_ate, &p, &q, params);
        gettimeofday(&end, NULL);
        
        float time_ms = timedifference_msec(start, end);
        ate_total += time_ms;
        printf("  Iteration %d: %f ms\n", i + 1, time_ms);
    }
    
    printf("%sAverage Ate pairing time: %f ms%s\n\n", 
           GREEN, ate_total / ITERATIONS, RESET);
    
    // Benchmark Optimal Ate pairing
    printf("Benchmarking Optimal Ate pairing...\n");
    for (int i = 0; i < ITERATIONS; i++) {
        gettimeofday(&start, NULL);
        optimal_ate_pairing(&e_opt_ate, &p, &q, params);
        gettimeofday(&end, NULL);
        
        float time_ms = timedifference_msec(start, end);
        opt_ate_total += time_ms;
        printf("  Iteration %d: %f ms\n", i + 1, time_ms);
    }
    
    printf("%sAverage Optimal Ate pairing time: %f ms%s\n\n", 
           GREEN, opt_ate_total / ITERATIONS, RESET);
    
    // Clean up
    g1_clear(&p);
    g2_clear(&q);
    gt_clear(&e_tate);
    gt_clear(&e_ate);
    gt_clear(&e_opt_ate);
}

/**
 * @brief Benchmark curve operations
 */
void benchmark_curve(CurveType curve_type) {
    CurveParams params;
    
    if (pairing_init(&params, curve_type) != 0) {
        printf("%s[ERROR]%s Failed to initialize curve\n", "\033[0;31m", RESET);
        return;
    }
    
    printf("\n%s=== Benchmarking %s Curve ===%s\n\n", BLUE, 
           curve_type == CURVE_BLS12 ? "BLS12" : 
           curve_type == CURVE_BN ? "BN" : "KSS16",
           RESET);
    
    benchmark_g1_scalar_mul(&params);
    benchmark_g2_scalar_mul(&params);
    benchmark_pairings(&params);
    
    pairing_clear(&params);
}

/**
 * @brief Main function
 */
int main() {
    printf("\n%s=== Elliptic Curve Pairing Benchmarks ===%s\n", BLUE, RESET);
    
    benchmark_curve(CURVE_BLS12);
    benchmark_curve(CURVE_BN);
    benchmark_curve(CURVE_KSS16);
    
    return 0;
}