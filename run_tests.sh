#!/bin/bash
set -e

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=== Pairing Library Test Runner ===${NC}"

# Step 1: Create build directory if it doesn't exist
echo -e "${BLUE}Step 1: Preparing build environment...${NC}"
mkdir -p build
cd build

# Step 2: Configure the build with CMake
echo -e "${BLUE}Step 2: Configuring with CMake...${NC}"
cmake .. -DCMAKE_BUILD_TYPE=Debug

# Step 3: Build the project
echo -e "${BLUE}Step 3: Building the project...${NC}"
make -j$(nproc 2>/dev/null || echo 2)

# Step 4: Run the pairing tests
echo -e "${BLUE}Step 4: Running pairing tests...${NC}"
if [ -f ./pairing_test ]; then
    echo "Running pairing tests with parameter initialization..."
    ./pairing_test
    TEST_RESULT=$?
    if [ $TEST_RESULT -eq 0 ]; then
        echo -e "${GREEN}All pairing tests passed!${NC}"
    else
        echo -e "${RED}Some pairing tests failed!${NC}"
        exit $TEST_RESULT
    fi
else
    echo -e "${RED}Error: pairing_test executable not found${NC}"
    exit 1
fi

# Step 5: Run benchmark tests if available
echo -e "${BLUE}Step 5: Running benchmarks...${NC}"
if [ -f ./bench_test ]; then
    echo "Running benchmark tests..."
    ./bench_test
else
    echo -e "${BLUE}Benchmark executable not found, skipping benchmarks${NC}"
fi

echo -e "${GREEN}All tests completed successfully!${NC}"
exit 0