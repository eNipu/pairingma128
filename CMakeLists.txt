cmake_minimum_required(VERSION 3.12)
project(pairingma128 C)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Set compiler flags
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -O2")
endif()

# Find GMP
find_library(GMP_LIBRARIES NAMES gmp)
find_path(GMP_INCLUDE_DIRS NAMES gmp.h)

if(NOT GMP_LIBRARIES OR NOT GMP_INCLUDE_DIRS)
    message(FATAL_ERROR "GMP library not found. Please install libgmp-dev.")
endif()

include_directories(${GMP_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

# Define parameter initialization functions
add_definitions(-DUSE_INIT_PARAMS)

# Create parameter initialization source file
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/param_init.c" 
"#include \"BLS12/BLS_12.h\"
#include \"BN/BN12_header.h\"
#include \"KSS16/KSS_16.h\"

void initialize_all_params() {
    // BLS12 and BN curve initialization
    init_parameters();
    set_parameters();
    
    // KSS16 initialization if header exists
    #ifdef KSS_16_H
    KSS_16_parameters();
    #endif
}

void clear_all_params() {
    // BLS12 and BN cleanup
    clear_parameters();
    
    // KSS16 cleanup if header exists
    #ifdef KSS_16_H
    dealloc_constants();
    #endif
}
")

# Add parameter initialization to sources
set(PARAM_INIT "${CMAKE_CURRENT_BINARY_DIR}/param_init.c")

# BLS12 library
file(GLOB BLS12_SOURCES "BLS12/*.c")
if(BLS12_SOURCES)
    list(FILTER BLS12_SOURCES EXCLUDE REGEX "BLS12/main.c$")
    add_library(bls12_lib ${BLS12_SOURCES})
    target_include_directories(bls12_lib PUBLIC BLS12)
    target_link_libraries(bls12_lib ${GMP_LIBRARIES})
else()
    message(STATUS "BLS12 source files not found, using minimal stub.")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/bls12_stub.c" "/* Empty stub file */")
    add_library(bls12_lib "${CMAKE_CURRENT_BINARY_DIR}/bls12_stub.c")
endif()

# BN library
file(GLOB BN_SOURCES "BN/*.c")
if(BN_SOURCES)
    list(FILTER BN_SOURCES EXCLUDE REGEX "BN/main.c$")
    add_library(bn_lib ${BN_SOURCES})
    target_include_directories(bn_lib PUBLIC BN)
    target_link_libraries(bn_lib ${GMP_LIBRARIES})
else()
    message(STATUS "BN source files not found, using minimal stub.")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/bn_stub.c" "/* Empty stub file */")
    add_library(bn_lib "${CMAKE_CURRENT_BINARY_DIR}/bn_stub.c")
endif()

# KSS16 library
file(GLOB KSS16_SOURCES "KSS16/*.c")
if(KSS16_SOURCES)
    list(FILTER KSS16_SOURCES EXCLUDE REGEX "KSS16/main.c$")
    add_library(kss16_lib ${KSS16_SOURCES})
    target_include_directories(kss16_lib PUBLIC KSS16)
    target_link_libraries(kss16_lib ${GMP_LIBRARIES})
else()
    message(STATUS "KSS16 source files not found, using minimal stub.")
    file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/kss16_stub.c" "/* Empty stub file */")
    add_library(kss16_lib "${CMAKE_CURRENT_BINARY_DIR}/kss16_stub.c")
endif()

# Pairing common library
add_library(pairing_lib pairing.c ${PARAM_INIT})
target_include_directories(pairing_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(pairing_lib bls12_lib bn_lib kss16_lib ${GMP_LIBRARIES})

# Build the individual main executables if they exist
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/BLS12/main.c")
    add_executable(bls12_main BLS12/main.c)
    target_link_libraries(bls12_main bls12_lib ${GMP_LIBRARIES})
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/BN/main.c")
    add_executable(bn_main BN/main.c)
    target_link_libraries(bn_main bn_lib ${GMP_LIBRARIES})
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/KSS16/main.c" AND NOT "${KSS16_SOURCES}" MATCHES "KSS16/main.c")
    add_executable(kss16_main KSS16/main.c)
    target_link_libraries(kss16_main kss16_lib ${GMP_LIBRARIES})
endif()

# Create test initialization file
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/test_init.c" 
"#include \"pairing.h\"
#include \"BLS12/BLS_12.h\"
#include \"BN/BN12_header.h\"
#include \"KSS16/KSS_16.h\"

extern void initialize_all_params();
extern void clear_all_params();

void test_init() {
    initialize_all_params();
}

void test_cleanup() {
    clear_all_params();
}
")

# Test executables
add_executable(pairing_test test/pairing_test.c "${CMAKE_CURRENT_BINARY_DIR}/test_init.c")
target_link_libraries(pairing_test pairing_lib ${GMP_LIBRARIES})

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/bench_test.c")
    add_executable(bench_test test/bench_test.c "${CMAKE_CURRENT_BINARY_DIR}/test_init.c")
    target_link_libraries(bench_test pairing_lib ${GMP_LIBRARIES})
endif()

# Create necessary directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test)

# Add a custom target to run all tests
add_custom_target(check
    COMMAND pairing_test
    DEPENDS pairing_test
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Running tests..."
)

# Add a custom target to run benchmarks
if(TARGET bench_test)
    add_custom_target(benchmark
        COMMAND bench_test
        DEPENDS bench_test
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Running benchmarks..."
    )
endif()

# Install targets
install(TARGETS bls12_lib bn_lib kss16_lib pairing_lib
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(FILES pairing.h BLS12/BLS_12.h BN/BN12_header.h KSS16/KSS_16.h 
        DESTINATION include/pairingma128)
