cmake_minimum_required(VERSION 3.10)
project(pairingma128)

set(CMAKE_C_STANDARD 99)

# Include directories
include_directories(BLS12 BN KSS16)

# Library source files (excluding main.c files)
set(BLS12_LIB_SRC
    BLS12/Fp.c
    BLS12/Fp2.c
    BLS12/Fp6.c
    BLS12/Fp12.c
    BLS12/EFp.c
    BLS12/EFp2.c
)

# Create static libraries
add_library(bls12_lib STATIC ${BLS12_LIB_SRC})

# Executables (just main.c files)
add_executable(bls12 BLS12/main.c)
add_executable(kss16 KSS16/main.c)

# Link libraries to executables
target_link_libraries(bls12 bls12_lib gmp)
target_link_libraries(kss16 gmp)

# Link gmp to the static libraries
target_link_libraries(bls12_lib gmp)
