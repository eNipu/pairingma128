cmake_minimum_required(VERSION 3.10)
project(kss16_pairing C)

find_library(GMP_LIBRARY gmp REQUIRED)
find_path(GMP_INCLUDE_DIR gmp.h REQUIRED)
set(KSS16_DIR ${CMAKE_CURRENT_LIST_DIR})

add_executable(kss16
    ${KSS16_DIR}/field.c
    ${KSS16_DIR}/curve.c
    ${KSS16_DIR}/pairing.c
    ${KSS16_DIR}/params.c
    ${KSS16_DIR}/util.c
    ${KSS16_DIR}/test.c
    ${KSS16_DIR}/main.c
)
target_include_directories(kss16 PRIVATE ${KSS16_DIR} ${GMP_INCLUDE_DIR})
target_link_libraries(kss16 PRIVATE ${GMP_LIBRARY})

add_executable(kss16_math_test
    ${KSS16_DIR}/field.c
    ${KSS16_DIR}/curve.c
    ${KSS16_DIR}/pairing.c
    ${KSS16_DIR}/params.c
    ${KSS16_DIR}/util.c
    ${CMAKE_CURRENT_LIST_DIR}/../tests/kss16_math_test.c
)
target_include_directories(kss16_math_test PRIVATE ${KSS16_DIR} ${CMAKE_CURRENT_LIST_DIR}/../tests ${GMP_INCLUDE_DIR})
target_link_libraries(kss16_math_test PRIVATE ${GMP_LIBRARY})
