# Gene Summariser -GFF3 QC Tool

This tool analyses a GFF/GFF3 (General Feature Format) genome annotation file and performs high quality control (QC) on gene and transcript models, and generate summary tables , plots and an HTML report .

## What the tool does :

Given:
- a **GFF / GFF3 file** (required)
- a **FASTA file** (optional)

1. Parses gene and transcript models from a GFF/GFF3 file.
2. Calculates QC metrics (exon counts, CDS presence, etc)
3. Flags suspicious gene models
4. Generates plots and an HTML report

## Folder structure

GROUPF_GFFPARSER/

├── data/ # Input files (optional convenience folder)

├── results/ # Output files (generated automatically)

├── src/ # Source code (no need to edit)

├── test/ # Automated tests

├── Dockerfile # Docker image definition

├── pyproject.toml # Project configuration

└── README.md # This file



# Option 1 : Run using Docker 🐳 (Recommended)
### Prerequisite : Docker installed
             https://www.docker.com/get-started/

## Step 1 Build docker image (One time)
    Run this **from the project root directory**
    (the folder containing `Dockerfile`):

    docker build -t gene-summariser:test .


## Step 2  Run the following commands from the project root:

    (Beginner Friendly)

    docker run --rm \
    -v "/inputpath:/data" \
    -v "/outputpath:/results" \
    gene-summariser:test \
    -g /data/test.gff -f /data/test.fasta --outdir /results

    # Replace `/inputpath` and `/outputpath` with actual folders on your computer.
    

## Using files from ANY location (Advanced)

    You are not required to use the data/ folder.

    Rule to remember
    -v HOST_PATH : CONTAINER_PATH

    Left side (HOST_PATH) → any folder on your computer (you may change this)

    Right side (CONTAINER_PATH) → path inside Docker (must match -g / -f)

    Example: GFF and FASTA in different folders (Windows)
    
      docker run --rm `
    -v "C:/Users/User/Documents/gff_files:/gff" `
    -v "C:/Users/User/Downloads/fasta_files:/fasta" `
    -v "${PWD}/results:/results" `
    gene-summariser:test `
    -g /gff/genes.gff3 `
    -f /fasta/genome.fasta `
    --outdir /results

    Example: GFF and FASTA in different folders (macOS / Linux)
    
      docker run --rm \
    -v "$HOME/gff_files:/gff" \
    -v "$HOME/fasta_files:/fasta" \
    -v "$(pwd)/results:/results" \
    gene-summariser:test \
    -g /gff/genes.gff3 \
    -f /fasta/genome.fasta \
    --outdir /results


### Having file path issues?

    This repository already includes a `data/` folder.

    If you encounter repeated errors related to file paths or file locations,
    copy your GFF/GFF3 and FASTA files into the `data/` folder and rerun
    the Quick Start Docker command.


    (Windows) Poweshell

    docker run --rm `
    -v "${PWD}/data:/data" `
    -v "${PWD}/results:/results" `
    gene-summariser:test `
    -g /data/test.gff `
    -f /data/test.fasta `
    --outdir /results

    macOS/ Linux 

    docker run --rm \
    -v "$(pwd)/data:/data" \
    -v "$(pwd)/results:/results" \
    gene-summariser:test \
    -g /data/test.gff \
    -f /data/test.fasta \
    --outdir /results

**Note:** If `${PWD}` causes issues, replace it with the full absolute path to your project directory.

          Example:   -v "C:/Users/Name/project/data:/data"



