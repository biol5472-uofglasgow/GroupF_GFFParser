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

    windows:    .\run_analysis.bat C:\User\Desktop\example.gff C:\User\Desktop\fasta.fasta  C:\Users\98184\remote\OneDrive\Desktop\TESTSOFTWARE
    macOS:      ./run_analysis.sh /Users/username/Desktop/example.gff /Users/username/Desktop/fasta.fasta /Users/username/OneDrive/Desktop/TESTSOFTWARE

    

### Note :

1. If any file or folder path contains spaces, wrap it in double quotes.
2. While relative paths may work, absolute paths are strongly recommended, especially when running via Docker, to avoid path resolution issues.
3. FASTA file must match the reference used in the GFF (if required)
4. All output files (summary tables, QC flags, and HTML report) will be written to the specified results directory..

