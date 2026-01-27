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