# Makefile for gene-summariser project
.PHONY: help install install-dev test lint format type-check clean docker-build docker-run

# Configuration
PYTHON := python3
PIP := $(PYTHON) -m pip
PYTEST := $(PYTHON) -m pytest
RUFF := $(PYTHON) -m ruff
MYPY := $(PYTHON) -m mypy
DOCKER := docker
DOCKER_IMAGE := gene-summariser
DOCKER_TAG := latest

# Colors for output
GREEN := \033[0;32m
YELLOW := \033[1;33m
NC := \033[0m # No Color

# Default target
.DEFAULT_GOAL := help

## help: Show this help message
help:
	@echo "$(GREEN)Gene Model Summariser - Available Commands$(NC)"
	@echo ""
	@echo "$(YELLOW)Setup:$(NC)"
	@echo "  make install          - Install package"
	@echo "  make install-dev      - Install with development dependencies"
	@echo ""
	@echo "$(YELLOW)Testing:$(NC)"
	@echo "  make test             - Run tests"
	@echo "  make test-cov         - Run tests with coverage"
	@echo "  make test-watch       - Run tests in watch mode"
	@echo ""
	@echo "$(YELLOW)Code Quality:$(NC)"
	@echo "  make lint             - Run linter"
	@echo "  make lint-fix         - Fix linting issues"
	@echo "  make format           - Format code"
	@echo "  make format-check     - Check code formatting"
	@echo "  make type-check       - Run type checker"
	@echo "  make quality          - Run all quality checks"
	@echo ""
	@echo "$(YELLOW)Docker:$(NC)"
	@echo "  make docker-build     - Build Docker image"
	@echo "  make docker-run       - Run Docker container"
	@echo "  make docker-test      - Test Docker container"
	@echo "  make docker-run-test  - Run with test data (one command)"
	@echo "  make docker-clean     - Remove Docker image"
	@echo ""
	@echo "$(YELLOW)Cleanup:$(NC)"
	@echo "  make clean            - Remove generated files"
	@echo "  make clean-all        - Remove all build artifacts"

## install: Install the package
install:
	@echo "$(GREEN)Installing gene-summariser...$(NC)"
	$(PIP) install -e .

## install-dev: Install with development dependencies
install-dev:
	@echo "$(GREEN)Installing gene-summariser with dev dependencies...$(NC)"
	$(PIP) install -e ".[dev]"

## test: Run tests
test:
	@echo "$(GREEN)Running tests...$(NC)"
	$(PYTEST)

## test-cov: Run tests with coverage
test-cov:
	@echo "$(GREEN)Running tests with coverage...$(NC)"
	$(PYTEST) --cov=gene_summariser --cov-report=html --cov-report=term-missing
	@echo "$(GREEN)Coverage report: htmlcov/index.html$(NC)"

## test-watch: Run tests in watch mode
test-watch:
	@echo "$(GREEN)Running tests in watch mode...$(NC)"
	$(PYTEST) -f

## lint: Run linter
lint:
	@echo "$(GREEN)Running linter...$(NC)"
	$(RUFF) check src/ tests/

## lint-fix: Fix linting issues
lint-fix:
	@echo "$(GREEN)Fixing linting issues...$(NC)"
	$(RUFF) check --fix src/ tests/

## format: Format code
format:
	@echo "$(GREEN)Formatting code...$(NC)"
	$(RUFF) format src/ tests/

## format-check: Check code formatting
format-check:
	@echo "$(GREEN)Checking code formatting...$(NC)"
	$(RUFF) format --check src/ tests/

## type-check: Run type checker
type-check:
	@echo "$(GREEN)Running type checker...$(NC)"
	$(MYPY) src/

## quality: Run all quality checks
quality: lint format-check type-check
	@echo "$(GREEN)All quality checks passed!$(NC)"

## docker-build: Build Docker image
docker-build:
	@echo "$(GREEN)Building Docker image...$(NC)"
	$(DOCKER) build -t $(DOCKER_IMAGE):$(DOCKER_TAG) .
	@echo "$(GREEN)Docker image built: $(DOCKER_IMAGE):$(DOCKER_TAG)$(NC)"

## docker-run: Run Docker container
docker-run:
	@echo "$(GREEN)Running Docker container...$(NC)"
	@echo "$(YELLOW)Note: Mount your data directory with -v$(NC)"
	$(DOCKER) run --rm $(DOCKER_IMAGE):$(DOCKER_TAG) --help

## docker-test: Test Docker container
docker-test:
	@echo "$(GREEN)Testing Docker container...$(NC)"
	$(DOCKER) run --rm $(DOCKER_IMAGE):$(DOCKER_TAG) --help
	@echo "$(GREEN)Docker container test passed!$(NC)"

## docker-run-test: Run Docker container with test data
docker-run-test:
	@echo "$(GREEN)Running Docker container with test data...$(NC)"
ifeq ($(OS),Windows_NT)
	$(DOCKER) run --rm \
		-v "%cd%/tests/fixtures:/data" \
		-v "%cd%/results:/app/results" \
		$(DOCKER_IMAGE):$(DOCKER_TAG) \
		--gff /data/models.gff3 \
		--fasta /data/testfasta.fasta \
		--outdir /app/results
else
	$(DOCKER) run --rm \
		-v "$(PWD)/tests/fixtures:/data" \
		-v "$(PWD)/results:/app/results" \
		$(DOCKER_IMAGE):$(DOCKER_TAG) \
		--gff /data/models.gff3 \
		--fasta /data/testfasta.fasta \
		--outdir /app/results
endif
	@echo "$(GREEN)Results saved to ./results/$(NC)"

## docker-clean: Remove Docker image
docker-clean:
	@echo "$(GREEN)Removing Docker image...$(NC)"
	$(DOCKER) rmi $(DOCKER_IMAGE):$(DOCKER_TAG) || true

## clean: Remove generated files
clean:
	@echo "$(GREEN)Cleaning generated files...$(NC)"
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf .ruff_cache/
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type f -name "*.db" -delete

## clean-all: Remove all build artifacts
clean-all: clean
	@echo "$(GREEN)Cleaning all build artifacts...$(NC)"
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf results/
	rm -rf qc_report.*
	rm -rf *.log

# Development workflow shortcuts
.PHONY: dev check run

## dev: Setup development environment
dev: install-dev
	@echo "$(GREEN)Development environment ready!$(NC)"

## check: Run all checks before commit
check: quality test
	@echo "$(GREEN)All checks passed! Ready to commit.$(NC)"

## run: Run the tool (example)
run:
	@echo "$(GREEN)Running gene-summariser...$(NC)"
	@echo "$(YELLOW)Example: gene-summariser --gff data/example.gff3$(NC)"