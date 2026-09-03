/* math_assert.h - tiny assertion harness for pairing math correctness tests.
 * Each curve test includes this and calls CHECK(cond, "name").
 * A check failure increments a counter; the test main prints a report and
 * returns nonzero if any check failed (so CI / the runner script sees it).
 * The pass/fail verdict is driven by exit code, NOT by stdout markers.
 */
#ifndef MATH_ASSERT_H
#define MATH_ASSERT_H

#include <stdio.h>

static int mtest_checks = 0;
static int mtest_fails  = 0;

#define CHECK(cond, name)                                                    \
    do {                                                                     \
        mtest_checks++;                                                      \
        if (cond) {                                                          \
            printf("  [ ok ] %s\n", name);                                   \
        } else {                                                             \
            mtest_fails++;                                                   \
            printf("  [FAIL] %s\n", name);                                   \
        }                                                                    \
    } while (0)

#define MTEST_REPORT(curve)                                                  \
    do {                                                                     \
        printf("==========================================\n");              \
        printf("%s: %d checks, %d failure(s)\n", curve, mtest_checks,        \
               mtest_fails);                                                 \
        printf("%s: %s\n", curve, mtest_fails == 0 ? "ALL PASS" : "FAILED"); \
        return mtest_fails == 0 ? 0 : 1;                                     \
    } while (0)

#endif /* MATH_ASSERT_H */
