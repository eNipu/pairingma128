cmake_minimum_required(VERSION 3.10)
project(bls_demo C)

find_library(GMP_LIBRARY gmp REQUIRED)
find_path(GMP_INCLUDE_DIR gmp.h REQUIRED)

set(BLS12_DIR ${CMAKE_CURRENT_LIST_DIR}/../../BLS12)
set(BLS_DIR ${CMAKE_CURRENT_LIST_DIR})

add_executable(bls_demo
    ${BLS_DIR}/bls.c
    ${BLS_DIR}/sha256.c
    ${BLS_DIR}/demo.c
    ${BLS12_DIR}/field.c
    ${BLS12_DIR}/curve.c
    ${BLS12_DIR}/pairing.c
    ${BLS12_DIR}/params.c
    ${BLS12_DIR}/util.c
)
target_include_directories(bls_demo PRIVATE ${BLS_DIR} ${BLS12_DIR} ${GMP_INCLUDE_DIR})
target_link_libraries(bls_demo PRIVATE ${GMP_LIBRARY})
