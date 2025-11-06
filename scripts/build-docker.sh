#!/usr/bin/env bash
# Build script for teloclip Docker image
# Supports multi-architecture builds with buildx

set -euo pipefail

# Configuration
IMAGE_NAME="${IMAGE_NAME:-teloclip}"
TAG="${TAG:-test}"
PLATFORMS="${PLATFORMS:-linux/amd64,linux/arm64}"
PUSH="${PUSH:-false}"

# Get version from git tags and convert to PEP 440 format
# Example: v0.3.2-6-gd5ce6fc-dirty -> 0.3.2.post6+gd5ce6fc.dirty
if [ -z "${VERSION:-}" ]; then
    RAW_VERSION=$(git describe --tags --always --dirty 2>/dev/null || echo '0.0.0')
    # Convert git describe output to PEP 440
    VERSION=$(echo "$RAW_VERSION" | \
        sed 's/^v//' | \
        sed 's/-\([0-9]\+\)-g/.post\1+g/' | \
        sed 's/-dirty/.dirty/' | \
        sed 's/-/+/')
fi

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Building teloclip Docker image ===${NC}"
echo "Image: ${IMAGE_NAME}:${TAG}"
echo "Version: ${VERSION}"
echo "Platforms: ${PLATFORMS}"
echo "Push to registry: ${PUSH}"
echo ""

# Check if buildx is available
if ! docker buildx version &> /dev/null; then
    echo -e "${RED}Error: docker buildx not found${NC}"
    echo "Please install Docker Buildx or use Docker Desktop"
    exit 1
fi

# Create builder instance if it doesn't exist
if ! docker buildx inspect teloclip-builder &> /dev/null; then
    echo -e "${YELLOW}Creating buildx builder instance...${NC}"
    docker buildx create --name teloclip-builder --use
fi

# Use the builder
docker buildx use teloclip-builder

# Build command
BUILD_CMD="docker buildx build"
BUILD_CMD+=" --platform ${PLATFORMS}"
BUILD_CMD+=" --tag ${IMAGE_NAME}:${TAG}"
BUILD_CMD+=" --build-arg VERSION=${VERSION}"

# Add push flag if requested
if [ "${PUSH}" = "true" ]; then
    BUILD_CMD+=" --push"
elif [[ "${PLATFORMS}" == *","* ]]; then
    # Multi-platform build without push requires output to docker-container
    echo -e "${YELLOW}Note: Multi-platform build cannot be loaded to local Docker${NC}"
    echo -e "${YELLOW}Building to build cache only. Use PUSH=true to push to registry.${NC}"
    # Remove --load, buildx will cache the images
else
    BUILD_CMD+=" --load"
fi

# Add build context
BUILD_CMD+=" ."

echo -e "${GREEN}Running: ${BUILD_CMD}${NC}"
echo ""

# Execute build
eval "${BUILD_CMD}"

if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}✓ Build successful!${NC}"
    echo ""
    echo "To test the image, run:"
    echo "  docker run --rm ${IMAGE_NAME}:${TAG} --version"
    echo "  docker run --rm ${IMAGE_NAME}:${TAG} --help"
    echo ""
    echo "To run tests:"
    echo "  ./scripts/test-docker.sh"
else
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
fi
