#!/usr/bin/env bash
# Test script for teloclip Docker image

# Provide VERSION build arg when building test image
# Convert to PEP 440: v0.3.2-6-gd5ce6fc-dirty → 0.3.2.post6+gd5ce6fc.dirty
# docker build -t teloclip:test --build-arg VERSION=$(git describe --tags --always --dirty | sed 's/^v//' | sed 's/-\([0-9]\+\)-/.post\1+/' | sed 's/-dirty/.dirty/') .
# ./scripts/test-docker.sh

set -euo pipefail

# Configuration
IMAGE_NAME="${IMAGE_NAME:-teloclip}"
TAG="${TAG:-test}"
TEST_DATA_DIR="tests/data"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test counter
TESTS_PASSED=0
TESTS_FAILED=0

echo -e "${GREEN}=== Testing teloclip Docker image ===${NC}"
echo "Image: ${IMAGE_NAME}:${TAG}"
echo ""

# Helper function to run tests
run_test() {
    local test_name="$1"
    local test_cmd="$2"

    echo -e "${YELLOW}▶ Test: ${test_name}${NC}"

    if eval "${test_cmd}"; then
        echo -e "${GREEN}✓ Passed${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "${RED}✗ Failed${NC}"
        ((TESTS_FAILED++))
    fi
    echo ""
}

# Test 1: Image exists
run_test "Image exists" \
    "docker image inspect ${IMAGE_NAME}:${TAG} > /dev/null 2>&1"

# Test 2: Version command
run_test "Version command" \
    "docker run --rm ${IMAGE_NAME}:${TAG} --version"

# Test 3: Help command
run_test "Help command" \
    "docker run --rm ${IMAGE_NAME}:${TAG} --help > /dev/null"

# Test 4: Filter subcommand help
run_test "Filter subcommand help" \
    "docker run --rm ${IMAGE_NAME}:${TAG} filter --help > /dev/null"

# Test 5: Extract subcommand help
run_test "Extract subcommand help" \
    "docker run --rm ${IMAGE_NAME}:${TAG} extract --help > /dev/null"

# Test 6: Extend subcommand help
run_test "Extend subcommand help" \
    "docker run --rm ${IMAGE_NAME}:${TAG} extend --help > /dev/null"

# Test 7: Run with test data (if available)
if [ -d "${TEST_DATA_DIR}" ] && [ -f "${TEST_DATA_DIR}/test.sam" ] && [ -f "${TEST_DATA_DIR}/test.fna.fai" ]; then
    run_test "Process test SAM file" \
        "docker run --rm -v \$(pwd)/${TEST_DATA_DIR}:/data ${IMAGE_NAME}:${TAG} filter --ref-idx /data/test.fna.fai /data/test.sam > /tmp/teloclip_test_output.sam"

    # Check output was created
    if [ -s /tmp/teloclip_test_output.sam ]; then
        echo -e "${GREEN}✓ Output file created${NC}"
        rm -f /tmp/teloclip_test_output.sam
    else
        echo -e "${YELLOW}⚠ Output file is empty (may be expected if no overhangs found)${NC}"
    fi
    echo ""
else
    echo -e "${YELLOW}⚠ Skipping data processing test (test data not found)${NC}"
    echo ""
fi

# Test 8: Check image size
echo -e "${YELLOW}▶ Image size check${NC}"
IMAGE_SIZE=$(docker image inspect ${IMAGE_NAME}:${TAG} --format='{{.Size}}' | awk '{print $1/1024/1024}')
echo "Image size: ${IMAGE_SIZE} MB"
if (( $(echo "${IMAGE_SIZE} < 150" | bc -l) )); then
    echo -e "${GREEN}✓ Image size is acceptable (<150MB)${NC}"
    ((TESTS_PASSED++))
else
    echo -e "${YELLOW}⚠ Image size is larger than expected (>150MB)${NC}"
fi
echo ""

# Summary
echo -e "${GREEN}=== Test Summary ===${NC}"
echo "Tests passed: ${TESTS_PASSED}"
echo "Tests failed: ${TESTS_FAILED}"
echo ""

if [ ${TESTS_FAILED} -eq 0 ]; then
    echo -e "${GREEN}✓ All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}✗ Some tests failed${NC}"
    exit 1
fi
