# Gene Model Summariser - Complete Docker Guide

**A comprehensive guide to building, deploying, and using the Gene Model Summariser with Docker.**

---

## Table of Contents

1. [Overview](#overview)
2. [Project Structure](#project-structure)
3. [Docker Configuration Files](#docker-configuration-files)
4. [Quick Start](#quick-start)
5. [Building Images](#building-images)
6. [Running Containers](#running-containers)
7. [Docker Compose Usage](#docker-compose-usage)
8. [Volume Mapping & Data Management](#volume-mapping--data-management)
9. [Common Use Cases](#common-use-cases)
10. [Troubleshooting](#troubleshooting)
11. [Advanced Usage](#advanced-usage)
12. [Best Practices](#best-practices)
13. [CI/CD Integration](#cicd-integration)

---

## Overview

The Gene Model Summariser is a quality control tool for GFF3 gene annotations. This guide covers everything you need to know about running the tool in Docker containers, from basic usage to advanced deployment scenarios.

**Key Features:**

- Containerized Python application
- Multi-stage builds for optimization
- Support for development and production environments
- Volume mounting for flexible data handling
- Docker Compose orchestration

---

## Project Structure

```
gene-summariser/
├── src/
│   └── gene_summariser/
│       ├── __init__.py       # Package initialization
│       ├── cli.py            # Command-line interface
│       ├── parser.py         # GFF3 parsing
│       ├── models.py         # Data structures
│       ├── qc.py             # QC checks
│       ├── metrics.py        # Metrics calculation
│       ├── writer.py         # Output generation
│       └── fasta.py          # FASTA sequence handling
├── tests/
│   ├── fixtures/             # Test data files
│   │   ├── models.gff3
│   │   ├── testfasta.fasta
│   │   └── annotations.gff3
│   ├── test_models.py
│   ├── test_parser.py
│   ├── test_qc.py
│   └── test_integration.py
├── data/                     # User's input files (create this directory)
│   └── .gitkeep              # (Directory not in repo by default)
├── results/                  # Output directory (created by container)
├── Dockerfile                # Docker container definition
├── docker-compose.yml        # Docker Compose configuration
├── pyproject.toml            # Project configuration
├── README.md                 # Main documentation
└── DOCKER_USAGE.md          # This file
```

---

## Docker Configuration Files

### Dockerfile

**Purpose:** Defines how to build the Gene Model Summariser container image.

```dockerfile
# Use official Python runtime as base image
FROM python:3.11-slim

# Set metadata
LABEL maintainer="BIOL5472_GROUP_F <your.email@glasgow.ac.uk>"
LABEL description="Gene Model Summariser - QC tool for GFF3 annotations"
LABEL version="0.1.0"

# Set working directory in container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Copy dependency files first (for layer caching)
COPY pyproject.toml /app/

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -e .

# Copy source code
COPY src/ /app/src/

# Install the package
RUN pip install --no-cache-dir -e .

# Create directories for data
RUN mkdir -p /data/input /data/output

# Set environment variables
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

# Default command (shows help)
ENTRYPOINT ["gene-summariser"]
CMD ["--help"]
```

**Key Features:**

- Based on Python 3.11 slim image for minimal size
- Optimized layer caching (dependencies before source code)
- System dependencies for compilation (gcc, g++)
- Environment variables for Python optimization
- Flexible entrypoint for easy command override

### docker-compose.yml

**Purpose:** Orchestrates multiple container configurations for different use cases.

```yaml
version: "3.8"

services:
  gene-summariser:
    build:
      context: .
      dockerfile: Dockerfile
    image: gene-summariser:latest
    container_name: gene-summariser-app
    volumes:
      # Mount test data from fixtures (or change to ./data for your own files)
      - ./tests/fixtures:/data
      - ./results:/app/results
    environment:
      - PYTHONUNBUFFERED=1
    # Override default command
    # Example: Run on test data
    command: >
      --gff /data/models.gff3
      --fasta /data/testfasta.fasta
      --outdir /app/results
    # For your own data, change volumes to:
    # - ./data:/data
    # And update the command accordingly

    # For interactive use, comment out command and use:
    # stdin_open: true
    # tty: true
    # entrypoint: /bin/bash

  # Optional: Development container with tests
  gene-summariser-dev:
    build:
      context: .
      dockerfile: Dockerfile
      target: development
    image: gene-summariser:dev
    container_name: gene-summariser-dev
    volumes:
      - .:/app
      - /app/.pytest_cache
      - /app/htmlcov
    environment:
      - PYTHONUNBUFFERED=1
    command: pytest --cov=gene_summariser --cov-report=html
```

**Key Features:**

- Production service for running analyses
- Development service with test support
- Pre-configured volume mounts
- Flexible command override
- Cache directories for test artifacts

---

## Quick Start

### Prerequisites

- Docker installed (version 20.10+)
- Docker Compose installed (version 1.29+)
- Input data files (GFF3, optionally FASTA)

**Note:** The repository includes test data in `tests/fixtures/` that you can use immediately for testing. For your own analysis, you'll need to create a `data/` directory and add your files.

### 1. Build the Docker Image

```bash
docker build -t gene-summariser:latest .
```

### 2. Prepare Your Data Directory

You have two options:

**Option A: Use test data from the repository**

```bash
mkdir -p results
# Use the test fixtures that are already in the repo
docker run --rm \
  -v $(pwd)/tests/fixtures:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/models.gff3 \
  --fasta /data/testfasta.fasta \
  --outdir /app/results
```

**Option B: Create a data directory for your own files**

```bash
mkdir -p data results
cp /path/to/your/annotations.gff3 data/
cp /path/to/your/genome.fasta data/  # Optional
```

### 3. Run Analysis

**Quick test with included test data:**

```bash
docker run --rm \
  -v $(pwd)/tests/fixtures:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/models.gff3 \
  --fasta /data/testfasta.fasta \
  --outdir /app/results
```

**Or with your own data (after creating data/ directory):**

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --fasta /data/genome.fasta \
  --outdir /app/results
```

### 4. View Results

```bash
ls results/
# Expected output:
# - transcript_summary.tsv
# - run.json
# - qc_flags.gff3
# - qc_flags.bed
```

---

## Building Images

### Standard Production Build

```bash
docker build -t gene-summariser:latest .
```

This creates an optimized image for running the tool in production.

### Development Build (with Test Dependencies)

```bash
docker build --target development -t gene-summariser:dev .
```

Note: This requires a multi-stage Dockerfile with a `development` target.

### Build with Specific Python Version

```bash
docker build --build-arg PYTHON_VERSION=3.10 -t gene-summariser:py310 .
```

Requires adding `ARG PYTHON_VERSION=3.11` to your Dockerfile.

### Using BuildKit for Faster Builds

```bash
DOCKER_BUILDKIT=1 docker build -t gene-summariser:latest .
```

BuildKit provides improved caching and parallel builds.

### Verify Build Success

```bash
# Check image exists
docker images gene-summariser

# Test the image
docker run --rm gene-summariser:latest --version
docker run --rm gene-summariser:latest --help
```

---

## Running Containers

### Basic Analysis (GFF3 Only)

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

### With FASTA File (for Start/Stop Codon Checks)

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --fasta /data/genome.fasta \
  --outdir /app/results
```

### Strict Mode (Fail on Warnings)

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --strict \
  --outdir /app/results
```

### Custom Output Format (JSON)

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --format json \
  --output /app/results/qc_report.json \
  --outdir /app/results
```

### Interactive Session

```bash
docker run --rm -it \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  --entrypoint bash \
  gene-summariser:latest
```

Then run commands inside the container:

```bash
gene-summariser --gff /data/annotations.gff3 --outdir /app/results
ls /app/results
exit
```

---

## Docker Compose Usage

### Run Standard Analysis

```bash
# Edit docker-compose.yml to configure your input files
docker-compose up gene-summariser
```

### Run in Detached Mode

```bash
docker-compose up -d gene-summariser
docker-compose logs -f gene-summariser
```

### Run Tests

```bash
docker-compose up gene-summariser-dev
```

### Interactive Development Environment

```bash
docker-compose run --rm gene-summariser-dev bash
```

Inside the container:

```bash
pytest tests/ -v
ruff check src/
mypy src/
gene-summariser --help
```

### Stop and Clean Up

```bash
# Stop running containers
docker-compose down

# Stop and remove volumes
docker-compose down -v
```

### Rebuild After Code Changes

```bash
docker-compose build
docker-compose up gene-summariser
```

---

## Volume Mapping & Data Management

### Understanding Volume Mounts

Docker volumes connect directories on your host machine to directories inside the container:

```bash
-v /host/path:/container/path
```

### Standard Volume Configuration

| Host Path   | Container Path | Purpose                   |
| ----------- | -------------- | ------------------------- |
| `./data`    | `/data`        | Input files (GFF3, FASTA) |
| `./results` | `/app/results` | Output files and reports  |

### Example Directory Structure

**Repository structure (what you have):**

```
gene-summariser/
├── src/
├── tests/
│   └── fixtures/
│       ├── models.gff3        # Test GFF3 file
│       ├── testfasta.fasta    # Test FASTA file
│       └── annotations.gff3   # Sample annotations
├── Dockerfile
├── docker-compose.yml
└── README.md
```

**Working directory structure (what you create):**

```
your-project/
├── tests/
│   └── fixtures/
│       ├── models.gff3        # Test data (from repo)
│       └── testfasta.fasta    # Test data (from repo)
├── data/                      # Your input files (create this)
│   ├── your_annotations.gff3
│   └── your_genome.fasta
├── results/                   # Output directory (created by container)
│   ├── transcript_summary.tsv
│   ├── run.json
│   ├── qc_flags.gff3
│   └── qc_flags.bed
├── Dockerfile
├── docker-compose.yml
└── README.md
```

### Creating Directories

**For testing with included test data:**

```bash
# Just create results directory - test data already exists in tests/fixtures/
mkdir -p results

# Run with test data
docker run --rm \
  -v $(pwd)/tests/fixtures:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/models.gff3 \
  --fasta /data/testfasta.fasta \
  --outdir /app/results
```

**For your own data:**

```bash
# Create directories before first run
mkdir -p data results

# Copy your data files
cp /path/to/your/data/annotations.gff3 data/
cp /path/to/your/data/genome.fasta data/

# Run with your data
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --fasta /data/genome.fasta \
  --outdir /app/results
```

### Handling Permissions

If you encounter permission issues:

```bash
# Run container with your user ID
docker run --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

### Read-Only Mounts

For input data protection:

```bash
docker run --rm \
  -v $(pwd)/data:/data:ro \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

---

## Common Use Cases

### Use Case 1: Quick QC Check

**Scenario:** Rapidly check annotation quality without saving detailed reports.

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --output /dev/stdout
```

### Use Case 2: Full Analysis with Reports

**Scenario:** Complete QC with all output files.

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --fasta /data/genome.fasta \
  --format json \
  --outdir /app/results
```

### Use Case 3: Batch Processing Multiple Files

**Scenario:** Process multiple GFF3 files in one go.

Create `process_all.sh`:

```bash
#!/bin/bash

for gff in data/*.gff3; do
  filename=$(basename "$gff" .gff3)

  echo "Processing: $filename"

  docker run --rm \
    -v $(pwd)/data:/data \
    -v $(pwd)/results:/app/results \
    gene-summariser:latest \
    --gff /data/$filename.gff3 \
    --outdir /app/results/$filename
done
```

Make executable and run:

```bash
chmod +x process_all.sh
./process_all.sh
```

### Use Case 4: Pipeline Integration

**Nextflow Example:**

```groovy
process QC_CHECK {
    container 'gene-summariser:latest'

    input:
    path gff
    path fasta

    output:
    path "results/*"

    script:
    """
    gene-summariser \
      --gff ${gff} \
      --fasta ${fasta} \
      --outdir results
    """
}
```

**Snakemake Example:**

```python
rule qc_check:
    input:
        gff="data/{sample}.gff3",
        fasta="data/{sample}.fasta"
    output:
        summary="results/{sample}/transcript_summary.tsv",
        qc_flags="results/{sample}/qc_flags.gff3"
    container:
        "gene-summariser:latest"
    shell:
        """
        gene-summariser \
          --gff {input.gff} \
          --fasta {input.fasta} \
          --outdir results/{wildcards.sample}
        """
```

### Use Case 5: Continuous Monitoring

**Scenario:** Monitor a directory for new files and process automatically.

```bash
#!/bin/bash
# watch_and_process.sh

inotifywait -m data/ -e create -e moved_to |
  while read path action file; do
    if [[ "$file" =~ .*gff3$ ]]; then
      echo "New file detected: $file"
      docker run --rm \
        -v $(pwd)/data:/data \
        -v $(pwd)/results:/app/results \
        gene-summariser:latest \
        --gff /data/$file \
        --outdir /app/results/$(basename $file .gff3)
    fi
  done
```

---

## Troubleshooting

### Issue: Permission Denied Errors

**Symptoms:**

```
Permission denied: '/app/results/transcript_summary.tsv'
```

**Solution:**

```bash
# Option 1: Run with your user ID
docker run --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results

# Option 2: Fix permissions after running
sudo chown -R $(whoami):$(whoami) results/
```

### Issue: File Not Found

**Symptoms:**

```
FileNotFoundError: [Errno 2] No such file or directory: '/data/annotations.gff3'
```

**Common Mistakes:**

❌ **WRONG** - Using host path:

```bash
docker run gene-summariser:latest --gff ~/data/file.gff3
```

✅ **CORRECT** - Using container path after mounting:

```bash
docker run -v ~/data:/data gene-summariser:latest --gff /data/file.gff3
```

### Issue: Memory/Resource Limits

**Symptoms:**

```
Killed
```

**Solution:**

```bash
# Increase memory limits
docker run --rm \
  --memory="4g" \
  --memory-swap="8g" \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/large_genome.gff3
```

### Issue: Container Immediately Exits

**Debug Steps:**

```bash
# Check container logs
docker ps -a  # Find container ID
docker logs <container_id>

# Run with verbose output
docker run --rm \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --verbose

# Get interactive shell
docker run --rm -it \
  -v $(pwd)/data:/data \
  --entrypoint bash \
  gene-summariser:latest
```

### Issue: Build Failures

**Symptoms:**

```
ERROR: failed to solve: process "/bin/sh -c pip install..." did not complete successfully
```

**Solutions:**

```bash
# Clear build cache
docker builder prune

# Build without cache
docker build --no-cache -t gene-summariser:latest .

# Check Dockerfile syntax
docker build --check -t gene-summariser:latest .
```

### Issue: Slow Performance

**Optimization Tips:**

```bash
# Use tmpfs for temporary files
docker run --rm \
  --tmpfs /tmp:rw,noexec,nosuid,size=1g \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3

# Limit CPU usage
docker run --rm \
  --cpus="2.0" \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3
```

### Debugging Checklist

- [ ] Verify Docker is running: `docker info`
- [ ] Check image exists: `docker images | grep gene-summariser`
- [ ] Verify volume paths: `ls data/` and `ls results/`
- [ ] Check file permissions: `ls -la data/`
- [ ] Test with sample data first
- [ ] Review container logs: `docker logs <container_id>`
- [ ] Try interactive mode: `docker run -it --entrypoint bash ...`

---

## Advanced Usage

### Custom QC Thresholds

If you extend the CLI to accept custom thresholds:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --max-exon-count 75 \
  --min-cds-length 50 \
  --max-intron-length 100000 \
  --outdir /app/results
```

### Environment Variables for Configuration

```bash
docker run --rm \
  -e LOG_LEVEL=DEBUG \
  -e MAX_WORKERS=4 \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3
```

### Multi-Stage Processing

```bash
# Stage 1: Initial QC
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/stage1:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results

# Stage 2: Process only flagged features
docker run --rm \
  -v $(pwd)/stage1:/data \
  -v $(pwd)/stage2:/app/results \
  gene-summariser:latest \
  --gff /data/qc_flags.gff3 \
  --strict \
  --outdir /app/results
```

### Running Tests Inside Container

```bash
docker run --rm \
  -v $(pwd):/app \
  gene-summariser:dev \
  pytest tests/ -v --cov=gene_summariser
```

### Extracting Specific Files

```bash
# Run analysis
CONTAINER_ID=$(docker run -d \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results)

# Wait for completion
docker wait $CONTAINER_ID

# Copy specific file
docker cp $CONTAINER_ID:/app/results/transcript_summary.tsv ./

# Clean up
docker rm $CONTAINER_ID
```

### Network Configuration

For tools that need internet access:

```bash
docker run --rm \
  --network host \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3
```

---

## Best Practices

### 1. Version Control

✅ **DO:** Use specific version tags in production

```bash
docker build -t gene-summariser:0.1.0 .
docker build -t gene-summariser:0.1 .
docker build -t gene-summariser:latest .
```

❌ **DON'T:** Rely solely on `latest` tag

```bash
# In production scripts
docker run gene-summariser:0.1.0  # ✓ Reproducible
docker run gene-summariser:latest  # ✗ May change
```

### 2. Data Management

✅ **DO:** Use volume mounts for data

```bash
-v $(pwd)/data:/data  # ✓ Flexible, easy updates
```

❌ **DON'T:** COPY data into image

```dockerfile
# Don't do this in Dockerfile
COPY data/ /app/data/  # ✗ Bloats image, not reusable
```

### 3. Resource Limits

✅ **DO:** Set appropriate limits for production

```bash
docker run --rm \
  --memory="4g" \
  --memory-swap="8g" \
  --cpus="2.0" \
  gene-summariser:latest
```

### 4. Cleanup

✅ **DO:** Regular maintenance

```bash
# Remove stopped containers
docker container prune

# Remove unused images
docker image prune -a

# Full cleanup (careful!)
docker system prune -a --volumes
```

### 5. Security

✅ **DO:** Run as non-root user when possible

```dockerfile
# Add to Dockerfile
RUN useradd -m -u 1000 appuser
USER appuser
```

✅ **DO:** Use read-only mounts for input data

```bash
-v $(pwd)/data:/data:ro
```

### 6. Logging

✅ **DO:** Configure appropriate log drivers

```bash
docker run --rm \
  --log-driver json-file \
  --log-opt max-size=10m \
  --log-opt max-file=3 \
  gene-summariser:latest
```

### 7. Documentation

✅ **DO:** Document your container usage

```bash
# Add to your project README
docker run --rm gene-summariser:latest --help
docker run --rm gene-summariser:latest --version
```

---

## CI/CD Integration

### GitHub Actions

```yaml
name: Docker Build and Test

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]

jobs:
  docker:
    runs-on: ubuntu-latest

    steps:
      - name: Checkout code
        uses: actions/checkout@v3

      - name: Set up Docker Buildx
        uses: docker/setup-buildx-action@v2

      - name: Build Docker image
        run: |
          docker build -t gene-summariser:latest .

      - name: Test image
        run: |
          docker run --rm gene-summariser:latest --help
          docker run --rm gene-summariser:latest --version

      - name: Run tests in container
        run: |
          docker run --rm \
            -v $(pwd)/tests:/app/tests \
            gene-summariser:latest \
            pytest tests/ -v

      - name: Login to Docker Hub
        if: github.event_name == 'push' && github.ref == 'refs/heads/main'
        uses: docker/login-action@v2
        with:
          username: ${{ secrets.DOCKER_USERNAME }}
          password: ${{ secrets.DOCKER_PASSWORD }}

      - name: Push to Docker Hub
        if: github.event_name == 'push' && github.ref == 'refs/heads/main'
        run: |
          docker tag gene-summariser:latest ${{ secrets.DOCKER_USERNAME }}/gene-summariser:latest
          docker push ${{ secrets.DOCKER_USERNAME }}/gene-summariser:latest
```

### GitLab CI

```yaml
stages:
  - build
  - test
  - deploy

variables:
  IMAGE_TAG: $CI_REGISTRY_IMAGE:$CI_COMMIT_SHORT_SHA

build:
  stage: build
  script:
    - docker build -t $IMAGE_TAG .
    - docker tag $IMAGE_TAG $CI_REGISTRY_IMAGE:latest

test:
  stage: test
  script:
    - docker run --rm $IMAGE_TAG --help
    - docker run --rm $IMAGE_TAG pytest tests/ -v

deploy:
  stage: deploy
  only:
    - main
  script:
    - docker login -u $CI_REGISTRY_USER -p $CI_REGISTRY_PASSWORD $CI_REGISTRY
    - docker push $IMAGE_TAG
    - docker push $CI_REGISTRY_IMAGE:latest
```

### Jenkins Pipeline

```groovy
pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                script {
                    docker.build("gene-summariser:${env.BUILD_ID}")
                }
            }
        }

        stage('Test') {
            steps {
                script {
                    docker.image("gene-summariser:${env.BUILD_ID}").inside {
                        sh 'pytest tests/ -v'
                    }
                }
            }
        }

        stage('Push') {
            when {
                branch 'main'
            }
            steps {
                script {
                    docker.withRegistry('https://registry.hub.docker.com', 'docker-hub-credentials') {
                        docker.image("gene-summariser:${env.BUILD_ID}").push('latest')
                    }
                }
            }
        }
    }
}
```

---

## Container Registry

### Pushing to Docker Hub

```bash
# Login to Docker Hub
docker login

# Tag your image
docker tag gene-summariser:latest yourusername/gene-summariser:0.1.0
docker tag gene-summariser:latest yourusername/gene-summariser:latest

# Push to Docker Hub
docker push yourusername/gene-summariser:0.1.0
docker push yourusername/gene-summariser:latest
```

### Pulling from Docker Hub

```bash
docker pull yourusername/gene-summariser:latest
```

### Using in Projects

Once published, others can use your image:

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  yourusername/gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

### Private Registry

```bash
# Login to private registry
docker login registry.company.com

# Tag for private registry
docker tag gene-summariser:latest registry.company.com/gene-summariser:0.1.0

# Push to private registry
docker push registry.company.com/gene-summariser:0.1.0
```

---

## Performance Considerations

### Image Size Optimization

```bash
# Check image size
docker images gene-summariser:latest

# Use multi-stage builds (example)
FROM python:3.11-slim as builder
# ... build steps ...

FROM python:3.11-slim
COPY --from=builder /app /app
# ... final image ...

# Use .dockerignore
echo "*.pyc" >> .dockerignore
echo "__pycache__" >> .dockerignore
echo ".git" >> .dockerignore
echo "tests/" >> .dockerignore
```

### Build Cache Optimization

```bash
# Enable BuildKit
export DOCKER_BUILDKIT=1

# Build with cache
docker build -t gene-summariser:latest .

# Build with cache from remote
docker build --cache-from yourusername/gene-summariser:latest \
  -t gene-summariser:latest .
```

### Runtime Performance

```bash
# CPU limits
docker run --cpus="2.0" gene-summariser:latest

# Memory limits
docker run --memory="2g" --memory-swap="4g" gene-summariser:latest

# Use tmpfs for temporary files
docker run --tmpfs /tmp:rw,size=1g gene-summariser:latest
```

---

## Support and Additional Resources

### Getting Help

**Documentation in Container:**

```bash
docker run --rm gene-summariser:latest --help
docker run --rm gene-summariser:latest --version
```

**Container Inspection:**

```bash
# View image details
docker inspect gene-summariser:latest

# View container logs
docker logs <container_id>

# Check running containers
docker ps

# Check all containers
docker ps -a
```

**Interactive Debugging:**

```bash
docker run --rm -it \
  -v $(pwd)/data:/data \
  --entrypoint bash \
  gene-summariser:latest
```

### Useful Docker Commands

```bash
# List images
docker images

# Remove image
docker rmi gene-summariser:latest

# List containers
docker ps -a

# Remove container
docker rm <container_id>

# View resource usage
docker stats

# Clean up everything
docker system prune -a --volumes
```

### Additional Documentation

- Main README: See `README.md` for tool usage
- QC Documentation: See `QC_CHECKS_DOCUMENTATION.md` for quality check details
- Docker Documentation: https://docs.docker.com
- Docker Hub: https://hub.docker.com

---

## Summary

This guide has covered:

- ✅ Project structure and Docker configuration files
- ✅ Building images for different environments
- ✅ Running containers with various options
- ✅ Docker Compose for orchestration
- ✅ Volume management and data handling
- ✅ Common use cases and workflows
- ✅ Troubleshooting and debugging
- ✅ Advanced usage patterns
- ✅ Best practices and security
- ✅ CI/CD integration

**Quick Reference:**

```bash
# Build
docker build -t gene-summariser:latest .

# Run
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results

# Debug
docker run --rm -it --entrypoint bash gene-summariser:latest
```

For issues or questions, please refer to the main project documentation or open an issue on the project repository.
