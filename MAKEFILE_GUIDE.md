# Makefile Usage Guide

## What is the Makefile?

The Makefile provides **convenient shortcuts** for common tasks. Instead of remembering long commands, you just type `make <command>`.

## Quick Start

### View All Available Commands

```bash
make help
```

**Output:**

```
Gene Model Summariser - Available Commands

Setup:
  make install          - Install package
  make install-dev      - Install with development dependencies

Testing:
  make test             - Run tests
  make test-cov         - Run tests with coverage
  make test-watch       - Run tests in watch mode

Code Quality:
  make lint             - Run linter
  make lint-fix         - Fix linting issues
  make format           - Format code
  make format-check     - Check code formatting
  make type-check       - Run type checker
  make quality          - Run all quality checks

Docker:
  make docker-build     - Build Docker image
  make docker-run       - Run Docker container
  make docker-test      - Test Docker container
  make docker-clean     - Remove Docker image

Cleanup:
  make clean            - Remove generated files
  make clean-all        - Remove all build artifacts
```

---

## Common Workflows

### Workflow 1: Initial Setup

```bash
# 1. Install with development dependencies
make install-dev

# Output:
# Installing gene-summariser with dev dependencies...
# Successfully installed gene-summariser
```

---

### Workflow 2: Development Cycle

```bash
# 1. Make code changes
vim src/gene_summariser/qc.py

# 2. Run tests
make test

# Output:
# Running tests...
# ======================== test session starts =========================
# tests/test_qc.py ..................................... [100%]
# ========================= 30 passed in 2.5s ==========================

# 3. Check code quality
make quality

# Output:
# Running linter...
# All done! ✨
# Checking code formatting...
# All done! ✨
# Running type checker...
# Success: no issues found
# All quality checks passed!

# 4. If issues found, fix them
make lint-fix
make format
```

---

### Workflow 3: Before Committing

```bash
# Run all checks to ensure everything is good
make check

# This runs:
# - make quality (lint, format, type-check)
# - make test

# Output:
# Running linter...
# Checking code formatting...
# Running type checker...
# Running tests...
# All checks passed! Ready to commit.
```

---

### Workflow 4: Docker Development

```bash
# 1. Build Docker image
make docker-build

# Output:
# Building Docker image...
# Successfully built abc123def456
# Docker image built: gene-summariser:latest

# 2. Test it works
make docker-test

# Output:
# Testing Docker container...
# usage: gene-summariser [-h] -g GFF [-f FASTA] ...
# Docker container test passed!

# 3. Run with your data (manual)
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

---

## Detailed Command Reference

### Setup Commands

#### `make install`

**What it does:**

- Installs the package in normal mode
- Installs runtime dependencies only

**When to use:**

- When you just want to use the tool (not develop it)

**Example:**

```bash
make install

# Equivalent to:
pip install -e .
```

---

#### `make install-dev`

**What it does:**

- Installs package with development dependencies
- Includes pytest, mypy, ruff, etc.

**When to use:**

- When starting development
- After cloning the repository
- When dependencies change

**Example:**

```bash
make install-dev

# Equivalent to:
pip install -e ".[dev]"

# Installs:
# - gene-summariser (your package)
# - pytest, pytest-cov (testing)
# - mypy (type checking)
# - ruff (linting/formatting)
# - pandas-stubs (type stubs)
```

---

### Testing Commands

#### `make test`

**What it does:**

- Runs all tests
- Shows pass/fail results

**When to use:**

- After making changes
- Before committing
- To verify everything works

**Example:**

```bash
make test

# Output:
# Running tests...
# ====================== test session starts =======================
# tests/test_models.py .............. [ 40%]
# tests/test_qc.py ................... [ 80%]
# tests/test_integration.py .....    [100%]
# ====================== 30 passed in 3.2s =========================

# Equivalent to:
pytest
```

---

#### `make test-cov`

**What it does:**

- Runs tests with coverage analysis
- Generates HTML coverage report
- Shows which lines aren't tested

**When to use:**

- To check test coverage
- To find untested code
- Before submission

**Example:**

```bash
make test-cov

# Output:
# Running tests with coverage...
# ====================== test session starts =======================
# tests/test_models.py .............. [ 40%]
# tests/test_qc.py ................... [ 80%]
# tests/test_integration.py .....    [100%]
#
# ---------- coverage: platform linux, python 3.11 ----------
# Name                            Stmts   Miss  Cover
# ---------------------------------------------------
# src/gene_summariser/__init__.py     8      0   100%
# src/gene_summariser/qc.py         150      5    97%
# src/gene_summariser/parser.py     120      8    93%
# ---------------------------------------------------
# TOTAL                             650     25    96%
#
# Coverage report: htmlcov/index.html

# Then open the HTML report:
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
start htmlcov/index.html  # Windows

# Equivalent to:
pytest --cov=gene_summariser --cov-report=html --cov-report=term-missing
```

---

### Code Quality Commands

#### `make lint`

**What it does:**

- Checks code for style issues
- Reports errors without fixing

**When to use:**

- To see what needs fixing
- Before committing
- In CI pipeline

**Example:**

```bash
make lint

# Output (if issues found):
# Running linter...
# src/gene_summariser/qc.py:45:1: E501 Line too long (120 > 100)
# src/gene_summariser/qc.py:67:1: W291 Trailing whitespace
# Found 2 errors

# Output (if clean):
# Running linter...
# All done! ✨ 0 errors found.

# Equivalent to:
ruff check src/ tests/
```

---

#### `make lint-fix`

**What it does:**

- Automatically fixes linting issues
- Modifies files in place

**When to use:**

- When lint reports errors
- To automatically fix style issues

**Example:**

```bash
make lint-fix

# Output:
# Fixing linting issues...
# Fixed 2 errors in src/gene_summariser/qc.py
# All done! ✨

# Equivalent to:
ruff check --fix src/ tests/
```

---

#### `make format`

**What it does:**

- Formats code according to style guide
- Makes code consistent

**When to use:**

- Before committing
- After writing new code
- To clean up formatting

**Example:**

```bash
make format

# Output:
# Formatting code...
# Reformatted src/gene_summariser/qc.py
# Reformatted src/gene_summariser/parser.py
# 2 files reformatted, 5 files left unchanged

# Equivalent to:
ruff format src/ tests/
```

---

#### `make type-check`

**What it does:**

- Checks type annotations
- Finds type errors

**When to use:**

- After adding type hints
- Before committing
- To catch type bugs

**Example:**

```bash
make type-check

# Output (if issues):
# Running type checker...
# src/gene_summariser/qc.py:67: error: Incompatible return value
# Found 1 error in 1 file

# Output (if clean):
# Running type checker...
# Success: no issues found in 7 source files

# Equivalent to:
mypy src/
```

---

#### `make quality`

**What it does:**

- Runs ALL quality checks
  - Linting
  - Format checking
  - Type checking

**When to use:**

- Before committing
- Before creating PR
- To ensure code quality

**Example:**

```bash
make quality

# Output:
# Running linter...
# All done! ✨
# Checking code formatting...
# All done! ✨
# Running type checker...
# Success: no issues found
# All quality checks passed!

# Equivalent to:
ruff check src/ tests/
ruff format --check src/ tests/
mypy src/
```

---

### Docker Commands

#### `make docker-build`

**What it does:**

- Builds Docker image
- Tags it as `gene-summariser:latest`

**When to use:**

- After changing Dockerfile
- Before running Docker container
- For deployment

**Example:**

```bash
make docker-build

# Output:
# Building Docker image...
# [+] Building 45.2s (12/12) FINISHED
# => [internal] load build definition
# => [internal] load .dockerignore
# => [internal] load metadata
# => [1/7] FROM docker.io/library/python:3.11-slim
# => [2/7] WORKDIR /app
# => [3/7] RUN apt-get update
# => [4/7] COPY pyproject.toml /app/
# => [5/7] RUN pip install
# => [6/7] COPY src/ /app/src/
# => [7/7] RUN pip install -e .
# => exporting to image
# => => naming to docker.io/library/gene-summariser:latest
# Docker image built: gene-summariser:latest

# Equivalent to:
docker build -t gene-summariser:latest .
```

---

#### `make docker-test`

**What it does:**

- Tests that Docker container works
- Runs `--help` command

**When to use:**

- After building Docker image
- To verify container works
- Before deploying

**Example:**

```bash
make docker-test

# Output:
# Testing Docker container...
# usage: gene-summariser [-h] -g GFF [-f FASTA] ...
#
# This program takes a gff file and performs necessary QC checks
#
# positional arguments:
#   -g GFF, --gff GFF     please write the path of the gffFile here.
#   ...
#
# Docker container test passed!

# Equivalent to:
docker run --rm gene-summariser:latest --help
```

---

#### `make docker-run`

**What it does:**

- Shows how to run Docker container
- Displays help message

**When to use:**

- As a reminder of Docker usage
- To verify container starts

**Example:**

```bash
make docker-run

# Output:
# Running Docker container...
# Note: Mount your data directory with -v
# usage: gene-summariser [-h] -g GFF ...

# For actual data processing, use:
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

---

#### `make docker-clean`

**What it does:**

- Removes Docker image
- Frees up disk space

**When to use:**

- When rebuilding from scratch
- To clean up old images
- When disk space is low

**Example:**

```bash
make docker-clean

# Output:
# Removing Docker image...
# Untagged: gene-summariser:latest
# Deleted: sha256:abc123...

# Equivalent to:
docker rmi gene-summariser:latest
```

---

### Cleanup Commands

#### `make clean`

**What it does:**

- Removes generated files
  - Test cache
  - Coverage reports
  - Python cache
  - Database files

**When to use:**

- When tests behave strangely
- Before fresh test run
- To clean up workspace

**Example:**

```bash
make clean

# Output:
# Cleaning generated files...
# Removed:
# - htmlcov/
# - .coverage
# - .pytest_cache/
# - __pycache__ directories
# - *.pyc files

# Files removed:
rm -rf htmlcov/
rm -rf .coverage
rm -rf .pytest_cache/
find . -name "__pycache__" -exec rm -rf {} +
find . -name "*.pyc" -delete
```

---

#### `make clean-all`

**What it does:**

- Removes ALL build artifacts
  - Everything from `make clean`
  - Build directories
  - Distribution files
  - Results directories

**When to use:**

- When starting fresh
- Before archiving project
- When build seems corrupted

**Example:**

```bash
make clean-all

# Output:
# Cleaning all build artifacts...
# Removed:
# - build/, dist/, *.egg-info/
# - results/, qc_report.*, *.log
# - All test and cache files

# Files removed:
# (everything from make clean, plus:)
rm -rf build/
rm -rf dist/
rm -rf *.egg-info/
rm -rf results/
```

---

## Practical Examples

### Example 1: Daily Development

```bash
# Morning: Start work
cd gene-summariser-group-f
source venv/bin/activate

# Make some changes
vim src/gene_summariser/qc.py

# Test your changes
make test

# Check code quality
make lint

# If linting errors, auto-fix
make lint-fix

# Format code
make format

# Final check before committing
make check

# If all passes, commit
git add .
git commit -m "Add new QC check"
git push
```

---

### Example 2: Before Creating PR

```bash
# Run complete quality check
make quality

# Run tests with coverage
make test-cov

# Build and test Docker
make docker-build
make docker-test

# If everything passes
git push origin feature/my-feature
# Then create PR on GitHub
```

---

### Example 3: After Cloning Repository

```bash
# Clone
git clone https://github.com/your-org/gene-summariser.git
cd gene-summariser

# Setup
make install-dev

# Verify everything works
make test

# Try the CLI
gene-summariser --help

# Build Docker
make docker-build
```

---

### Example 4: Fixing CI Failures

```bash
# CI failed on linting
# Run locally to see errors
make lint

# Auto-fix
make lint-fix

# Run all quality checks
make quality

# If passes, push
git add .
git commit -m "Fix linting issues"
git push
```

---

## Troubleshooting

### "make: command not found"

**Problem**: Make is not installed

**Solution**:

```bash
# macOS
xcode-select --install

# Ubuntu/Debian
sudo apt-get install build-essential

# Windows
# Use Git Bash or WSL
```

---

### "No rule to make target"

**Problem**: Typo in command or Makefile not in current directory

**Solution**:

```bash
# Check you're in project root
pwd
ls Makefile

# Check available commands
make help

# Use correct command name
make test  # not "make tests"
```

---

### Commands Run But Do Nothing

**Problem**: Package not installed

**Solution**:

```bash
# Install first
make install-dev

# Then try again
make test
```

---

## Customizing the Makefile

You can add your own commands:

```makefile
# Add to Makefile

## my-command: Description of what it does
my-command:
	@echo "Running my custom command..."
	python scripts/my_script.py
```

Then use:

```bash
make my-command
```

---

## Summary

### Most Used Commands

```bash
make install-dev   # Setup (once)
make test          # Test code (often)
make lint-fix      # Fix style (often)
make format        # Format code (often)
make check         # Before commit (always)
make docker-build  # Build container (when needed)
```

### Quick Reference Table

| Task         | Command             | Frequency               |
| ------------ | ------------------- | ----------------------- |
| Setup        | `make install-dev`  | Once                    |
| Test         | `make test`         | After every change      |
| Lint         | `make lint`         | Before commit           |
| Fix lint     | `make lint-fix`     | When lint fails         |
| Format       | `make format`       | Before commit           |
| Type check   | `make type-check`   | Before commit           |
| All checks   | `make check`        | Before every commit     |
| Docker build | `make docker-build` | When Dockerfile changes |
| Clean up     | `make clean`        | When needed             |

---

**Pro Tip**: Type `make` with no arguments to see the help message!

```bash
make
# Shows: make help
```
