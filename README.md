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
### Prequisite : Docker installed
    https://www.docker.com/get-started/

## Step 1 Build docker image (One time)
    Run this **from the project root directory**
    (the folder containing `Dockerfile`):

    
    docker build -t gene-summariser:latest .


## Step 2  Run the following commands from the project root:

windows:  .\run_analysis.bat <GFF_file_Path> <Fasta_File_Path> <Results_file_path> 

mac: .\run_analysis.sh <GFF_file_Path> <Fasta_File_Path> <Results_file_path>



    Example:  

    windows:    .\run_analysis.bat C:\User\Desktop\example.gff C:\User\Desktop\fasta.fasta  C:\Users\98184\remote\OneDrive\Desktop\result\
    macOS:      ./run_analysis.sh /Users/username/Desktop/example.gff /Users/username/Desktop/fasta.fasta /Users/username/OneDrive/Desktop/result\

    
# Option 2 Run without Docker (Python only)
    This option allows users to run the tool locally using Python, without Docker.

## Step 1 Create and activate a virtual environment (Optional but recommended)

    windows: python -m venv venv   venv\Scripts\activate                                                                        
    macOS:   python3 -m venv venv  source venv/bin/activate

## Step 2 

pip install .

## Step 3

gene-summariser --gff <PATH_TO_GFF> --fasta <PATH_TO_FASTA> --outdir <OUTPUT_DIRECTORY>





### Note :

1. If any file or folder path contains spaces, wrap it in double quotes.
2. While relative paths may work, absolute paths are strongly recommended, especially when running via Docker, to avoid path resolution issues.
3. FASTA file must match the reference used in the GFF (if required)
4. All output files (summary tables, QC flags, and HTML report) will be written to the specified results directory.


## Strict Mode

The tool supports an optional strict mode that enforces more stringent quality control checks during analysis.

When strict mode is enabled, the program will fail fast if critical QC issues are detected, instead of continuing with warnings.

What strict mode does

When --strict is enabled:

Example :

Windows: gene-summariser --gff "C:\Users\USERNAME\Desktop\large.gff" --fasta "C:\Users\USERNAME\Desktop\large.fasta" --outdir "C:\Users\USERNAME\Desktop\results" --strict
macOS: gene-summariser --gff "/Users/username/Desktop/large.gff" --fasta "/Users/username/Desktop/large.fasta" --outdir "/Users/username/Desktop/results" --strict


Treats selected QC issues as errors instead of warnings

Stops execution if critical problems are found in the input data

Ensures outputs are only produced for high-quality, valid inputs

Helps catch problems early in automated or production workflows

