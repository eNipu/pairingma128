/**
 * @file pairing.h
 * @brief Common interface for pairing operations across different curves
 * 
 * This file provides a unified interface for bilinear pairing operations
 * across different elliptic curve implementations: BLS12, BN, and KSS16.
 */

#ifndef PAIRING_H
#define PAIRING_H

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <sys/time.h>

/**
 * @brief Enumeration of supported curve types
 */
typedef enum {
    CURVE_BLS12,  /**< BLS12 curve with embedding degree k=12 */
    CURVE_BN,     /**< Barreto-Naehrig curve with embedding degree k=12 */
    CURVE_KSS16   /**< KSS curve with embedding degree k=16 */
} CurveType;

/**
 * @brief Structure to hold curve-specific parameters
 */
typedef struct {
    CurveType type;        /**< Type of the curve */
    mpz_t prime;           /**< Base field characteristic */
    mpz_t order;           /**< Order of the curve */
    mpz_t trace;           /**< Trace of Frobenius */
    mpz_t parameter;       /**< Curve defining parameter */
    void *curve_specific;  /**< Pointer to curve-specific data */
} CurveParams;

/**
 * @brief Generic structure for a point on G1
 */
typedef struct {
    CurveType type;        /**< Type of the curve */
    void *point_data;      /**< Pointer to curve-specific point data */
    int infinity;          /**< Flag for point at infinity */
} G1Point;

/**
 * @brief Generic structure for a point on G2
 */
typedef struct {
    CurveType type;        /**< Type of the curve */
    void *point_data;      /**< Pointer to curve-specific point data */
    int infinity;          /**< Flag for point at infinity */
} G2Point;

/**
 * @brief Generic structure for an element in GT
 */
typedef struct {
    CurveType type;        /**< Type of the curve */
    void *element_data;    /**< Pointer to curve-specific element data */
} GTElement;

/**
 * @brief Initialize curve parameters
 * 
 * @param params Pointer to curve parameters structure
 * @param type Type of curve to initialize
 * @return 0 on success, error code otherwise
 */
int pairing_init(CurveParams *params, CurveType type);

/**
 * @brief Clean up curve parameters
 * 
 * @param params Pointer to curve parameters structure
 */
void pairing_clear(CurveParams *params);

/**
 * @brief Initialize a point in G1
 * 
 * @param point Pointer to G1Point structure
 * @param params Curve parameters
 */
void g1_init(G1Point *point, CurveParams *params);

/**
 * @brief Initialize a point in G2
 * 
 * @param point Pointer to G2Point structure
 * @param params Curve parameters
 */
void g2_init(G2Point *point, CurveParams *params);

/**
 * @brief Initialize an element in GT
 * 
 * @param element Pointer to GTElement structure
 * @param params Curve parameters
 */
void gt_init(GTElement *element, CurveParams *params);

/**
 * @brief Generate a random point in G1
 * 
 * @param point Pointer to G1Point structure
 * @param params Curve parameters
 */
void g1_random(G1Point *point, CurveParams *params);

/**
 * @brief Generate a random point in G2
 * 
 * @param point Pointer to G2Point structure
 * @param params Curve parameters
 */
void g2_random(G2Point *point, CurveParams *params);

/**
 * @brief Perform scalar multiplication in G1
 * 
 * @param result Output point P*s
 * @param point Input point P
 * @param scalar Scalar value s
 * @param params Curve parameters
 */
void g1_scalar_mul(G1Point *result, G1Point *point, mpz_t scalar, CurveParams *params);

/**
 * @brief Perform scalar multiplication in G2
 * 
 * @param result Output point P*s
 * @param point Input point P
 * @param scalar Scalar value s
 * @param params Curve parameters
 */
void g2_scalar_mul(G2Point *result, G2Point *point, mpz_t scalar, CurveParams *params);

/**
 * @brief Compute Tate pairing e(P,Q)
 * 
 * @param result Output GT element
 * @param p G1 point P
 * @param q G2 point Q
 * @param params Curve parameters
 */
void tate_pairing(GTElement *result, G1Point *p, G2Point *q, CurveParams *params);

/**
 * @brief Compute Ate pairing e(P,Q)
 * 
 * @param result Output GT element
 * @param p G1 point P
 * @param q G2 point Q
 * @param params Curve parameters
 */
void ate_pairing(GTElement *result, G1Point *p, G2Point *q, CurveParams *params);

/**
 * @brief Compute optimal Ate pairing e(P,Q)
 * 
 * @param result Output GT element
 * @param p G1 point P
 * @param q G2 point Q
 * @param params Curve parameters
 */
void optimal_ate_pairing(GTElement *result, G1Point *p, G2Point *q, CurveParams *params);

/**
 * @brief Test pairing bilinearity
 * 
 * This function tests whether e(aP,Q) = e(P,aQ) = e(P,Q)^a
 * 
 * @param params Curve parameters
 * @return 1 if test passes, 0 otherwise
 */
int test_pairing_bilinearity(CurveParams *params);

/**
 * @brief Benchmark pairing operations
 * 
 * @param params Curve parameters
 */
void benchmark_pairing(CurveParams *params);

/**
 * @brief Clean up a G1 point
 * 
 * @param point Pointer to G1Point structure
 */
void g1_clear(G1Point *point);

/**
 * @brief Clean up a G2 point
 * 
 * @param point Pointer to G2Point structure
 */
void g2_clear(G2Point *point);

/**
 * @brief Clean up a GT element
 * 
 * @param element Pointer to GTElement structure
 */
void gt_clear(GTElement *element);

#endif /* PAIRING_H */