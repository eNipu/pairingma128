/* bls.h - BLS short signatures (Boneh-Lynn-Shacham) on the BLS12 curve.
 *
 * Serialization conventions (all hex strings):
 *   - secret key:  "<sk>"                      (one big integer)
 *   - G1 point:    "<x>:<y>"                    (signatures, hashes)
 *   - G2 point:    "<x0>:<x1>:<y0>:<y1>"       (public keys; twist coords)
 *   - keygen:      "<sk>|<pk>"
 *   - sig list:    "<sig1>;<sig2>;..."
 *   - pk list:     "<pk1>;<pk2>;..."
 *   - msg list:    "<msg1>\n<msg2>\n..."       (newline-separated)
 */
#ifndef BLS_H
#define BLS_H

#include <stddef.h>
#include <stdint.h>

int  bls_init(void);
char *bls_keygen(void);
char *bls_sign(const char *sk_hex, const uint8_t *msg, size_t msglen);
int  bls_verify(const char *pk_hex, const char *sig_hex, const uint8_t *msg, size_t msglen);
char *bls_aggregate(const char *sigs_hex);
int  bls_aggregate_verify(const char *agg_hex, const char *pks_hex, const char *msgs_hex);

#endif
