#!/bin/bash
# Run unit tests for pic-mpi
#
# Usage: ./run_tests.sh [build_dir]
#   build_dir: path to pic-mpi build directory (default: ../build)

# Don't use set -e, it causes issues with arithmetic and test failures

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${1:-$SCRIPT_DIR/../build}"
TEST_DIR="$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "========================================"
echo "PIC-MPI Unit Tests"
echo "========================================"
echo "Build directory: $BUILD_DIR"
echo "Test directory: $TEST_DIR"
echo ""

# Check if build directory exists
if [ ! -d "$BUILD_DIR" ]; then
    echo -e "${RED}Error: Build directory not found: $BUILD_DIR${NC}"
    echo "Please build pic-mpi first with: cmake -B build && cmake --build build"
    exit 1
fi

# Create test build directory
TEST_BUILD_DIR="$BUILD_DIR/test_rma"
mkdir -p "$TEST_BUILD_DIR"

# Find the library
if [ -f "$BUILD_DIR/libpic-mpi.a" ]; then
    LIB_FILE="$BUILD_DIR/libpic-mpi.a"
else
    echo -e "${RED}Error: Cannot find libpic-mpi.a in $BUILD_DIR${NC}"
    exit 1
fi

# Find dependency libraries
PIC_LIB=""
if [ -f "$BUILD_DIR/_deps/pic-build/libpic.a" ]; then
    PIC_LIB="$BUILD_DIR/_deps/pic-build/libpic.a"
fi

# Find module files
MOD_PATH="$BUILD_DIR/modules"
if [ ! -d "$MOD_PATH" ]; then
    MOD_PATH="$BUILD_DIR"
fi

# Also need pic modules
PIC_MOD_PATH="$BUILD_DIR/_deps/pic-build/modules"

echo "Library: $LIB_FILE"
echo "PIC Library: $PIC_LIB"
echo "Modules: $MOD_PATH"
echo ""

# Compiler settings - must use mpifort for mpi_f08 module
FC="mpifort"
FFLAGS="-O2 -g -I$MOD_PATH -I$PIC_MOD_PATH"
LIBS="$LIB_FILE $PIC_LIB"

# Verify mpifort is available
if ! command -v $FC &> /dev/null; then
    echo -e "${RED}Error: mpifort not found in PATH${NC}"
    echo "Please load your MPI module or set PATH"
    exit 1
fi

echo "Compiler: $FC ($(which $FC))"

# Tests to run (test_name:num_procs)
TESTS=(
    "test_rma_win_allocate:2"
    "test_rma_get_put:2"
    "test_rma_nonblocking:4"
    "test_rma_accumulate:4"
    "test_collectives:4"
    "test_communicators:4"
    "test_point_to_point:2"
    "test_darrays:4"
)

PASSED=0
FAILED=0
SKIPPED=0

compile_test() {
    local test_name=$1
    local src_file="$TEST_DIR/${test_name}.f90"
    local exe_file="$TEST_BUILD_DIR/$test_name"

    echo -n "Compiling $test_name... "

    if [ ! -f "$src_file" ]; then
        echo -e "${YELLOW}SKIPPED (source not found)${NC}"
        return 1
    fi

    if $FC $FFLAGS -o "$exe_file" "$src_file" $LIBS 2>/dev/null; then
        echo -e "${GREEN}OK${NC}"
        return 0
    else
        echo -e "${RED}FAILED${NC}"
        # Show error details
        $FC $FFLAGS -o "$exe_file" "$src_file" $LIBS 2>&1 | head -20
        return 1
    fi
}

run_test() {
    local test_name=$1
    local nprocs=$2
    local exe_file="$TEST_BUILD_DIR/$test_name"

    echo ""
    echo "----------------------------------------"
    echo "Running: $test_name (np=$nprocs)"
    echo "----------------------------------------"

    if [ ! -x "$exe_file" ]; then
        echo -e "${YELLOW}SKIPPED (not compiled)${NC}"
        SKIPPED=$((SKIPPED + 1))
        return
    fi

    # Run the test
    if mpirun --oversubscribe -np $nprocs "$exe_file" 2>&1; then
        PASSED=$((PASSED + 1))
    else
        echo -e "${RED}Test execution failed${NC}"
        FAILED=$((FAILED + 1))
    fi
}

# Compile all tests
echo "Compiling tests..."
echo ""
for test_spec in "${TESTS[@]}"; do
    test_name="${test_spec%%:*}"
    compile_test "$test_name" || true
done

echo ""
echo "========================================"
echo "Running tests..."
echo "========================================"

# Run all tests
for test_spec in "${TESTS[@]}"; do
    test_name="${test_spec%%:*}"
    nprocs="${test_spec##*:}"
    run_test "$test_name" "$nprocs"
done

echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo -e "Passed:  ${GREEN}$PASSED${NC}"
echo -e "Failed:  ${RED}$FAILED${NC}"
echo -e "Skipped: ${YELLOW}$SKIPPED${NC}"
echo "========================================"

if [ $FAILED -gt 0 ]; then
    exit 1
fi
exit 0
