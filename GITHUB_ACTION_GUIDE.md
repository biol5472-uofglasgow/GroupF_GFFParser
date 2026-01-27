# GitHub Actions CI/CD - Complete Guide

## What It Does Automatically

Every time you push code or create a PR, GitHub Actions:

### 1. **Linting** (Checks Code Style)

- Runs `ruff check src/ tests/`
- Ensures consistent code formatting
- Reports style violations

### 2. **Type Checking** (Validates Type Hints)

- Runs `mypy src/`
- Catches type errors before runtime
- Ensures type safety

### 3. **Testing** (Runs All Tests)

- Tests on Python 3.9, 3.10, and 3.11
- Generates coverage report
- Uploads to Codecov (optional)

### 4. **Integration Tests** (End-to-End)

- Tests complete workflows
- Ensures everything works together

### 5. **Docker Build** (Container Build)

- Builds Docker image
- Verifies it works
- Tests `--help` command

### 6. **Security Scan** (Finds Vulnerabilities)

- Checks dependencies for security issues
- Runs security linters

## How to Use It

### Basic Workflow

```bash
# 1. Make changes to your code
vim src/gene_summariser/qc.py

# 2. Commit and push
git add .
git commit -m "Add new QC check"
git push origin feature/new-qc-check

# 3. GitHub Actions runs automatically!
# View results at: https://github.com/your-org/your-repo/actions
```

### For Pull Requests

```bash
# 1. Create a branch
git checkout -b feature/amazing-feature

# 2. Make changes and commit
git add .
git commit -m "Add amazing feature"
git push origin feature/amazing-feature

# 3. Create PR on GitHub
# GitHub Actions runs on your PR!
# You'll see status checks at the bottom of the PR

# 4. Fix any failures
# Push fixes to the same branch
git add .
git commit -m "Fix linting errors"
git push origin feature/amazing-feature

# Actions re-run automatically!
```

### Viewing Results

#### In the Actions Tab

1. Go to **your-repo → Actions**
2. Click on a workflow run
3. See detailed results:

```
✓ Lint Code (12s)
✓ Type Check (8s)
✓ Test Python 3.10 (1m 23s)
  ✓ Install dependencies
  ✓ Run tests
  ✓ Upload coverage
✓ Docker Build (2m 15s)
```

#### In Pull Requests

At the bottom of your PR, you'll see:

```
All checks have passed ✓
- Lint Code ✓
- Type Check ✓
- Test Python 3.10 ✓
- Docker Build ✓

This branch has no conflicts with the base branch
[Merge pull request]
```

## What Each Job Does

### Lint Job

```yaml
lint:
  runs-on: ubuntu-latest
  steps:
    - Checkout code
    - Install Python
    - Install dependencies
    - Run ruff linter
    - Check formatting
```

**Purpose**: Ensures code follows style guidelines

**Fixes Issues With**:

```bash
# Locally run what CI runs
ruff check src/ tests/
ruff format src/ tests/
```

---

### Type Check Job

```yaml
type-check:
  runs-on: ubuntu-latest
  steps:
    - Checkout code
    - Install Python
    - Install dependencies
    - Run mypy
```

**Purpose**: Catches type errors

**Fixes Issues With**:

```bash
# Locally run what CI runs
mypy src/
```

---

### Test Job

```yaml
test:
  strategy:
    matrix:
      python-version: ["3.9", "3.10", "3.11"]
  steps:
    - Checkout code
    - Setup Python ${{ matrix.python-version }}
    - Install dependencies
    - Run pytest with coverage
    - Upload coverage
```

**Purpose**: Ensures all tests pass on multiple Python versions

**Fixes Issues With**:

```bash
# Locally run what CI runs
pytest --cov=gene_summariser --cov-report=term
```

---

### Docker Build Job

```yaml
docker-build:
  needs: [lint, type-check, test] # Only runs if these pass
  steps:
    - Checkout code
    - Setup Docker Buildx
    - Build Docker image
    - Test Docker image
```

**Purpose**: Ensures Docker container builds and works

**Fixes Issues With**:

```bash
# Locally run what CI runs
docker build -t gene-summariser:test .
docker run --rm gene-summariser:test --help
```

---

## Troubleshooting CI Failures

### Failure: Linting Errors

**Error Message:**

```
❌ Lint Code failed
src/gene_summariser/qc.py:45:1: E501 Line too long (120 > 100 characters)
```

**Fix:**

```bash
# Auto-fix locally
ruff check --fix src/ tests/
ruff format src/ tests/

# Commit and push
git add .
git commit -m "Fix linting errors"
git push
```

---

### Failure: Type Check Errors

**Error Message:**

```
❌ Type Check failed
src/gene_summariser/qc.py:67: error: Incompatible return value type
```

**Fix:**

```bash
# Check locally
mypy src/

# Fix the type issue in your code
# Then commit
git add .
git commit -m "Fix type annotations"
git push
```

---

### Failure: Tests Failed

**Error Message:**

```
❌ Test Python 3.10 failed
tests/test_qc.py::test_check_no_cds FAILED
```

**Fix:**

```bash
# Run failing test locally
pytest tests/test_qc.py::test_check_no_cds -vv

# Debug and fix
# Then commit
git add .
git commit -m "Fix failing test"
git push
```

---

### Failure: Docker Build

**Error Message:**

```
❌ Docker Build failed
ERROR: failed to solve: failed to compute cache key
```

**Fix:**

```bash
# Test Docker build locally
docker build -t gene-summariser:test .

# Fix Dockerfile if needed
# Then commit
git add Dockerfile
git commit -m "Fix Docker build"
git push
```

---

## Advanced: Workflow Triggers

### Current Triggers

The workflow runs on:

- Push to `main` or `develop` branches
- Pull requests to `main` or `develop`
- Release creation

```yaml
on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main, develop]
  release:
    types: [published]
```

### Custom Triggers

**Run on all branches:**

```yaml
on:
  push:
  pull_request:
```

**Run manually:**

```yaml
on:
  workflow_dispatch: # Adds "Run workflow" button in Actions tab
```

**Schedule (run daily):**

```yaml
on:
  schedule:
    - cron: "0 0 * * *" # Every day at midnight
```

---

## Status Badges

Add build status to your README:

```markdown
[![CI](https://github.com/your-org/your-repo/actions/workflows/ci.yml/badge.svg)](https://github.com/your-org/your-repo/actions/workflows/ci.yml)
```

This shows:

- ✓ Green badge when passing
- ✗ Red badge when failing

---

## Secrets Configuration

For Docker publishing (optional):

1. Go to **Settings → Secrets and variables → Actions**
2. Add secrets:
   - `DOCKER_USERNAME`: Your Docker Hub username
   - `DOCKER_PASSWORD`: Your Docker Hub password

Then the workflow can publish images on release.

---

## Best Practices

### 1. Run Checks Locally First

```bash
# Before pushing, run what CI will run
make check

# Or manually:
pytest
ruff check src/ tests/
mypy src/
```

### 2. Fix Failures Quickly

- CI failures block merging
- Fix and push to the same branch
- Actions re-run automatically

### 3. Review Workflow Logs

- Click on failed jobs
- Read error messages
- Reproduce locally
- Fix and push

### 4. Keep Workflows Fast

- Current workflow: ~5 minutes
- Optimize if it gets slower
- Use caching (already implemented)

---

## For Your Individual Report

### Document Your CI/CD Usage

**Example:**

> "I configured GitHub Actions CI/CD to automatically validate all code changes.
> The workflow runs linting, type checking, and tests on three Python versions,
> ensuring code quality before merging. I resolved 5+ CI failures during
> development, demonstrating debugging skills."
>
> Evidence:
>
> - CI configuration: `.github/workflows/ci.yml`
> - Successful workflow runs: [link to Actions tab]
> - PR with CI checks: [link to PR]

---

## Quick Reference

| Task       | Command            | CI Equivalent    |
| ---------- | ------------------ | ---------------- |
| Lint       | `ruff check src/`  | Lint job         |
| Format     | `ruff format src/` | Lint job         |
| Type check | `mypy src/`        | Type Check job   |
| Test       | `pytest`           | Test job         |
| Coverage   | `pytest --cov`     | Test job         |
| Docker     | `docker build .`   | Docker Build job |
| **All**    | `make check`       | **Full CI**      |

---

## Summary

### What CI Does For You

1. ✓ Catches errors before merging
2. ✓ Ensures tests pass on multiple Python versions
3. ✓ Maintains code quality
4. ✓ Validates Docker builds
5. ✓ Provides confidence in your code

### How to Use It

1. **Push code** → CI runs automatically
2. **View results** in Actions tab
3. **Fix failures** if any
4. **Merge** when all checks pass

### Key URLs

- **Actions**: `https://github.com/your-org/your-repo/actions`
- **Workflow file**: `.github/workflows/ci.yml`
- **Logs**: Click on any workflow run

---

**Remember**: Green checkmarks in GitHub Actions = Your code is production-ready! ✅
