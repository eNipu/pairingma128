cmake_minimum_required(VERSION 3.10)
project(pairingma128)

set(CMAKE_C_STANDARD 99)

# Include directories
include_directories(BLS12 BN KSS16)

# Source files
set(BLS12_SRC
    BLS12/main.c
    BLS12/Fp.c
    BLS12/Fp2.c
    BLS12/Fp6.c
    BLS12/Fp12.c
    BLS12/EFp.c
    BLS12/EFp2.c
    BLS12/EFp6.c
    BLS12/EFp12.c
    BLS12/pairing.c
)

set(BN_SRC
    BN/main.c
    BN/Fp.c
    BN/Fp2.c
    BN/Fp6.c
    BN/Fp12.c
    BN/EFp.c
    BN/EFp2.c
    BN/EFp6.c
    BN/EFp12.c
    BN/pairing.c
)

set(KSS16_SRC
    KSS16/main.c
    KSS16/Fp.c
    KSS16/Fp2.c
    KSS16/Fp6.c
    KSS16/Fp12.c
    KSS16/EFp.c
    KSS16/EFp2.c
    KSS16/EFp6.c
    KSS16/EFp12.c
    KSS16/pairing.c
)

# Executables
add_executable(bls12 ${BLS12_SRC})
add_executable(bn ${BN_SRC})
add_executable(kss16 ${KSS16_SRC})

# Link libraries
target_link_libraries(bls12 gmp)
target_link_libraries(bn gmp)
target_link_libraries(kss16 gmp)
