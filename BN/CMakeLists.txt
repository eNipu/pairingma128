cmake_minimum_required(VERSION 3.10)
project(bn_pairing C)

find_library(GMP_LIBRARY gmp REQUIRED)
find_path(GMP_INCLUDE_DIR gmp.h REQUIRED)
set(BN_DIR ${CMAKE_CURRENT_LIST_DIR})

# main executable (original main(), split across modules)
add_executable(bn
    ${BN_DIR}/field.c
    ${BN_DIR}/curve.c
    ${BN_DIR}/pairing.c
    ${BN_DIR}/params.c
    ${BN_DIR}/util.c
    ${BN_DIR}/test.c
    ${BN_DIR}/main.c
)
target_include_directories(bn PRIVATE ${BN_DIR} ${GMP_INCLUDE_DIR})
target_link_libraries(bn PRIVATE ${GMP_LIBRARY})

# mathematical-correctness suite (TDD gate; excludes main.c and test.c)
add_executable(bn_math_test
    ${BN_DIR}/field.c
    ${BN_DIR}/curve.c
    ${BN_DIR}/pairing.c
    ${BN_DIR}/params.c
    ${BN_DIR}/util.c
    ${CMAKE_CURRENT_LIST_DIR}/../tests/bn_math_test.c
)
target_include_directories(bn_math_test PRIVATE ${BN_DIR} ${CMAKE_CURRENT_LIST_DIR}/../tests ${GMP_INCLUDE_DIR})
target_link_libraries(bn_math_test PRIVATE ${GMP_LIBRARY})
