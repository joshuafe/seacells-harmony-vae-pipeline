#!/bin/bash
# Docker Helper Script for SEACells Pipeline
# Makes it easy to build and run the pipeline in Docker

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
IMAGE_NAME="seacells-vae-pipeline"
CONTAINER_NAME="seacells-pipeline-run"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo -e "${BLUE}================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================${NC}"
}

print_success() {
    echo -e "${GREEN}✓ $1${NC}"
}

print_error() {
    echo -e "${RED}✗ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}→ $1${NC}"
}

show_help() {
    cat << EOF
Docker Helper Script for SEACells Pipeline

Usage: $0 [COMMAND] [OPTIONS]

Commands:
    build       Build the Docker image
    run         Run the pipeline in Docker
    shell       Open an interactive shell in the container
    clean       Remove Docker image and containers
    logs        Show logs from the last run
    help        Show this help message

Run Options:
    --skip-seacells       Skip SEACells and Harmony step
    --skip-clustering     Skip clustering analysis step
    --n-trials N          Number of Optuna trials (default: 50)

Examples:
    # Build the image
    $0 build

    # Run full pipeline
    $0 run --n-trials 50

    # Run only VAE (skip earlier steps)
    $0 run --skip-seacells --skip-clustering --n-trials 100

    # Open interactive shell for debugging
    $0 shell

    # View logs
    $0 logs

EOF
}

build_image() {
    print_header "Building Docker Image"

    if ! command -v docker &> /dev/null; then
        print_error "Docker is not installed or not in PATH"
        exit 1
    fi

    print_info "Building image: ${IMAGE_NAME}"
    print_info "This may take 10-20 minutes on first build..."

    cd "$SCRIPT_DIR"
    docker build -t ${IMAGE_NAME} .

    if [ $? -eq 0 ]; then
        print_success "Docker image built successfully"
        print_info "Image name: ${IMAGE_NAME}"
    else
        print_error "Docker build failed"
        exit 1
    fi
}

run_pipeline() {
    print_header "Running Pipeline in Docker"

    # Check if image exists
    if ! docker image inspect ${IMAGE_NAME} &> /dev/null; then
        print_error "Docker image '${IMAGE_NAME}' not found"
        print_info "Run: $0 build"
        exit 1
    fi

    # Prepare volume mounts
    DATA_DIR="${SCRIPT_DIR}/data"
    SEACELLS_DIR="${SCRIPT_DIR}/seacells_output"
    VAE_DIR="${SCRIPT_DIR}/metacell_vae_output"

    # Create output directories if they don't exist
    mkdir -p "${SEACELLS_DIR}"
    mkdir -p "${VAE_DIR}"

    # Remove old container if exists
    docker rm -f ${CONTAINER_NAME} 2>/dev/null || true

    print_info "Starting container: ${CONTAINER_NAME}"
    print_info "Data directory: ${DATA_DIR}"
    print_info "SEACells output: ${SEACELLS_DIR}"
    print_info "VAE output: ${VAE_DIR}"
    print_info "Pipeline arguments: $@"
    echo ""

    # Run the container
    docker run --name ${CONTAINER_NAME} \
        -v "${DATA_DIR}:/app/data:ro" \
        -v "${SEACELLS_DIR}:/app/seacells_output" \
        -v "${VAE_DIR}:/app/metacell_vae_output" \
        ${IMAGE_NAME} "$@"

    if [ $? -eq 0 ]; then
        print_success "Pipeline completed successfully"
    else
        print_error "Pipeline failed. Check logs with: $0 logs"
        exit 1
    fi
}

open_shell() {
    print_header "Opening Interactive Shell"

    if ! docker image inspect ${IMAGE_NAME} &> /dev/null; then
        print_error "Docker image '${IMAGE_NAME}' not found"
        print_info "Run: $0 build"
        exit 1
    fi

    DATA_DIR="${SCRIPT_DIR}/data"
    SEACELLS_DIR="${SCRIPT_DIR}/seacells_output"
    VAE_DIR="${SCRIPT_DIR}/metacell_vae_output"

    mkdir -p "${SEACELLS_DIR}"
    mkdir -p "${VAE_DIR}"

    print_info "Opening bash shell in container..."
    print_info "Environment is activated and ready to go"
    print_info "Try: python scripts/seacells_batch_correction.py --help"
    echo ""

    docker run -it --rm \
        -v "${DATA_DIR}:/app/data:ro" \
        -v "${SEACELLS_DIR}:/app/seacells_output" \
        -v "${VAE_DIR}:/app/metacell_vae_output" \
        --entrypoint /bin/bash \
        ${IMAGE_NAME}
}

show_logs() {
    print_header "Container Logs"

    if docker container inspect ${CONTAINER_NAME} &> /dev/null; then
        docker logs ${CONTAINER_NAME}
    else
        print_error "Container '${CONTAINER_NAME}' not found"
        print_info "Run the pipeline first with: $0 run"
        exit 1
    fi
}

clean_docker() {
    print_header "Cleaning Docker Resources"

    print_info "Removing container: ${CONTAINER_NAME}"
    docker rm -f ${CONTAINER_NAME} 2>/dev/null || print_info "Container not found"

    print_info "Removing image: ${IMAGE_NAME}"
    docker rmi -f ${IMAGE_NAME} 2>/dev/null || print_info "Image not found"

    print_success "Cleanup complete"
}

# Main script logic
COMMAND="${1:-help}"

case "$COMMAND" in
    build)
        build_image
        ;;
    run)
        shift  # Remove 'run' from arguments
        run_pipeline "$@"
        ;;
    shell)
        open_shell
        ;;
    logs)
        show_logs
        ;;
    clean)
        clean_docker
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        print_error "Unknown command: $COMMAND"
        echo ""
        show_help
        exit 1
        ;;
esac
