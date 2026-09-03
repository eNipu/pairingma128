cmake_minimum_required(VERSION 3.10)
project(bls12_pairing C)

find_library(GMP_LIBRARY gmp REQUIRED)
find_path(GMP_INCLUDE_DIR gmp.h REQUIRED)
set(BLS12_DIR ${CMAKE_CURRENT_LIST_DIR})

add_executable(bls12
    ${BLS12_DIR}/field.c
    ${BLS12_DIR}/curve.c
    ${BLS12_DIR}/pairing.c
    ${BLS12_DIR}/params.c
    ${BLS12_DIR}/util.c
    ${BLS12_DIR}/test.c
    ${BLS12_DIR}/main.c
)
target_include_directories(bls12 PRIVATE ${BLS12_DIR} ${GMP_INCLUDE_DIR})
target_link_libraries(bls12 PRIVATE ${GMP_LIBRARY})

add_executable(bls12_math_test
    ${BLS12_DIR}/field.c
    ${BLS12_DIR}/curve.c
    ${BLS12_DIR}/pairing.c
    ${BLS12_DIR}/params.c
    ${BLS12_DIR}/util.c
    ${CMAKE_CURRENT_LIST_DIR}/../tests/bls12_math_test.c
)
target_include_directories(bls12_math_test PRIVATE ${BLS12_DIR} ${CMAKE_CURRENT_LIST_DIR}/../tests ${GMP_INCLUDE_DIR})
target_link_libraries(bls12_math_test PRIVATE ${GMP_LIBRARY})
