# Docker Usage Guide

This guide explains how to use the Gene Model Summariser with Docker.

## Quick Start

### 1. Build the Docker Image

```bash
docker build -t gene-summariser:latest .
```

### 2. Run with Your Data

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/your_annotations.gff3 \
  --fasta /data/your_genome.fasta \
  --outdir /app/results
```

---

## Detailed Usage

### Building the Image

#### Standard Build

```bash
docker build -t gene-summariser:latest .
```

#### Development Build (with tests)

```bash
docker build --target development -t gene-summariser:dev .
```

#### Build with specific Python version

```bash
docker build --build-arg PYTHON_VERSION=3.10 -t gene-summariser:py310 .
```

---

### Running the Container

#### Basic Analysis

```bash
docker run --rm \
  -v /path/to/your/data:/data \
  -v /path/to/output:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

#### With FASTA File (for start/stop codon checks)

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --fasta /data/genome.fasta \
  --outdir /app/results
```

#### Strict Mode (fail on warnings)

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --strict \
  --outdir /app/results
```

#### Custom Output Format

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

---

### Using Docker Compose

#### Run Analysis

```bash
# Edit docker-compose.yml to set your files
docker-compose up gene-summariser
```

#### Run Tests

```bash
docker-compose up gene-summariser-dev
```

#### Interactive Development

```bash
docker-compose run --rm gene-summariser-dev bash
```

---

## Volume Mapping

The container expects data in specific locations:

| Host Path   | Container Path | Purpose                   |
| ----------- | -------------- | ------------------------- |
| `./data`    | `/data`        | Input files (GFF3, FASTA) |
| `./results` | `/app/results` | Output files              |

### Example Directory Structure

```
your-project/
├── data/
│   ├── annotations.gff3
│   └── genome.fasta
├── results/          # Created by container
│   ├── transcript_summary.tsv
│   ├── run.json
│   ├── qc_flags.gff3
│   └── qc_flags.bed
└── docker-compose.yml
```

---

## Common Use Cases

### 1. Quick QC Check

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --output /dev/stdout
```

### 2. Batch Processing

Create a script `process_all.sh`:

```bash
#!/bin/bash

for gff in data/*.gff3; do
  filename=$(basename "$gff" .gff3)

  docker run --rm \
    -v $(pwd)/data:/data \
    -v $(pwd)/results:/app/results \
    gene-summariser:latest \
    --gff /data/$filename.gff3 \
    --outdir /app/results/$filename
done
```

### 3. Integration with Pipelines

```bash
# Nextflow example
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

---

## Troubleshooting

### Permission Issues

If you encounter permission errors:

```bash
# Run with your user ID
docker run --rm \
  --user $(id -u):$(id -g) \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

### File Not Found

Ensure paths are **inside the container**:

```bash
# ✗ WRONG (host path)
docker run gene-summariser:latest --gff ~/data/file.gff3

# ✓ CORRECT (container path after mounting)
docker run -v ~/data:/data gene-summariser:latest --gff /data/file.gff3
```

### Memory Issues

For large genomes, increase Docker memory:

```bash
docker run --rm \
  --memory="4g" \
  --memory-swap="8g" \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/large_genome.gff3
```

### Debugging

Run interactively:

```bash
docker run --rm -it \
  -v $(pwd)/data:/data \
  --entrypoint bash \
  gene-summariser:latest
```

Then inside container:

```bash
gene-summariser --gff /data/annotations.gff3 --outdir /app/results
```

---

## Advanced Usage

### Custom QC Thresholds

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --max-exon-count 75 \
  --min-cds-length 50 \
  --outdir /app/results
```

Note: You'll need to modify the CLI to accept these parameters.

### Extracting Specific Output

```bash
# Get only the summary file
docker run --rm \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results \
  && docker cp <container_id>:/app/results/transcript_summary.tsv ./
```

### Running Tests Inside Container

```bash
docker run --rm \
  -v $(pwd):/app \
  gene-summariser:dev \
  pytest tests/ -v
```

---

## Container Registry

### Push to Docker Hub

```bash
# Tag the image
docker tag gene-summariser:latest yourusername/gene-summariser:0.1.0
docker tag gene-summariser:latest yourusername/gene-summariser:latest

# Push to Docker Hub
docker push yourusername/gene-summariser:0.1.0
docker push yourusername/gene-summariser:latest
```

### Pull from Docker Hub

```bash
docker pull yourusername/gene-summariser:latest
```

### Use in Your Projects

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  yourusername/gene-summariser:latest \
  --gff /data/annotations.gff3
```

---

## Performance Considerations

### Image Size

Check image size:

```bash
docker images gene-summariser:latest
```

Optimize:

```bash
# Multi-stage build (already implemented)
# Remove unnecessary files
# Use slim base image (python:3.11-slim)
```

### Build Cache

Speed up builds:

```bash
# Use BuildKit
DOCKER_BUILDKIT=1 docker build -t gene-summariser:latest .
```

### Resource Limits

```bash
docker run --rm \
  --cpus="2.0" \
  --memory="2g" \
  -v $(pwd)/data:/data \
  gene-summariser:latest \
  --gff /data/annotations.gff3
```

---

## CI/CD Integration

### GitHub Actions

```yaml
- name: Build Docker image
  run: docker build -t gene-summariser:latest .

- name: Run tests in container
  run: docker run --rm gene-summariser:latest pytest

- name: Push to registry
  run: |
    echo "${{ secrets.DOCKER_PASSWORD }}" | docker login -u "${{ secrets.DOCKER_USERNAME }}" --password-stdin
    docker push yourusername/gene-summariser:latest
```

---

## Best Practices

1. **Always use version tags** in production:

   ```bash
   docker run gene-summariser:0.1.0  # ✓ Good
   docker run gene-summariser:latest  # ⚠️ May change
   ```

2. **Use volume mounts** for data (not COPY in Dockerfile):

   ```bash
   -v $(pwd)/data:/data  # ✓ Flexible
   ```

3. **Set memory limits** for large files:

   ```bash
   --memory="4g"
   ```

4. **Clean up** unused containers:

   ```bash
   docker system prune -a
   ```

5. **Use `.dockerignore`** to keep images small

---

## Documentation in Docker Image

View help from container:

```bash
docker run --rm gene-summariser:latest --help
```

View version:

```bash
docker run --rm gene-summariser:latest --version
```

---

## Example Workflows

### Workflow 1: Simple QC

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

### Workflow 2: Full Analysis with Genome

```bash
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --fasta /data/genome.fasta \
  --format json \
  --strict \
  --outdir /app/results
```

### Workflow 3: Development

```bash
# Start development container
docker run --rm -it \
  -v $(pwd):/app \
  --entrypoint bash \
  gene-summariser:dev

# Inside container
pytest tests/
ruff check src/
mypy src/
```

---

## Support

For issues with Docker:

- Check logs: `docker logs <container_id>`
- Inspect image: `docker inspect gene-summariser:latest`
- Debug: `docker run -it --entrypoint bash gene-summariser:latest`

For issues with the tool itself, see the main README.md.
