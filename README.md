# Gene Summariser -GFF3 QC Tool

This tool analyses a GFF/GFF3 (General Feature Format) genome annotation file and performs high quality control (QC) on gene and transcript models, and generate summary tables , plots and an HTML report .

## What the tool does :

Given:
- a **GFF / GFF3 file** (required)
- a **FASTA file** (optional)

## 🚀 Features

- Parses gene and transcript models from a GFF/GFF3 file.
- Calculates QC metrics (exon counts, CDS presence, etc)
- Detects common gene model QC issues:
  - Missing start codon
  - Missing stop codon
  - CDS phase inconsistencies
  - Single-exon transcripts
  - Transcripts without CDS
- Produces **visualisation-ready BED files**
- Outputs structured **GFF3 annotations**
- Generates a ** HTML QC report**
- Captures **full run provenance** for reproducibility

## Folder structure

GROUPF_GFFPARSER/

    ├── data/ # Input files (optional convenience folder)
    ├── results/ # Output files (generated automatically)
    ├── src/ # Source code (no need to edit)
    ├── test/ # Automated tests
    ├── Dockerfile # Docker image definition
    ├── pyproject.toml # Project configuration
    └── README.md # This file



# Option 1 : Run using Docker (Recommended)
### Prequisite : Docker installed
    https://www.docker.com/get-started/

## Step 1 Build docker image (One time)
    Run this **from the project root directory**
    (the folder containing `Dockerfile`):

    
    docker build -t gene-summariser:latest .


## Step 2  Run the following commands from the project root:

    windows:  .\run_analysis.bat <PATH_TO_GFF> <PATH_TO_FASTA> <OUTPUT_DIRECTORY> 

    mac: ./run_analysis.sh <PATH_TO_GFF> <PATH_TO_FASTA> <OUTPUT_DIRECTORY>



    Example:  

    windows:    .\run_analysis.bat .\test\fixtures\models.gff3 .\test\fixtures\testfasta.fasta  .\results\testrun
    macOS:      ./run_analysis.sh ./test/fixtures/models.gff3 ./test/fixtures/testfasta.fasta ./results/testrun

## File handeling during run
    The above commands should work out of the box assuming you are located at the root directory of the project and cloned the full repository
    Data can however be used from any location on your device, assuming you follow your operating systems notation for absolute paths

# Option 2 Run without Docker (Python only)
    This option allows users to run the tool locally using Python, without Docker.

## Step 1 Create and activate a virtual environment (Optional but recommended)

## Step 2 

    pip install .

## Step 3

    gene-summariser --gff <PATH_TO_GFF> --fasta <PATH_TO_FASTA> --outdir <OUTPUT_DIRECTORY>





### Note :

1. If any file or folder path contains spaces, wrap it in double quotes.
2. While relative paths may work, absolute paths are strongly recommended, especially when running via Docker, to avoid path resolution issues.
3. If you use relative paths when running the script, prefix them with `./` on macOS/Linux and `.\` on Windows.
4. FASTA file must match the reference used in the GFF (if required)
5. All output files (summary tables, QC flags, and HTML report) will be written to the specified results directory.
   


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




## Testing

Run the test suite locally using `pytest`. Tests include unit tests for core logic and an integration test using fixture inputs.



## Output Files

| File | Description |
|----|----|
| `qc_flags.bed` | Genomic intervals of QC-flagged transcripts (BED format) |
| `qc_flags.gff3` | Feature-rich QC annotations with biological context |
| `transcript_summary.tsv` | Tabular summary of transcript-level QC results |
| `qc_report.html` | Interactive HTML report with statistics and examples |
| `run.json` | Run metadata, parameters, and provenance |


## BED File Usage (Genome Browsers)

The `qc_flags.bed` file can be directly loaded into genome browsers for **visual inspection and validation**.

### Supported browsers
- **IGV (Integrative Genomics Viewer)**

### Viewing QC flags in IGV

QC flags are provided as a BED file and can be visualised using IGV.

1. Open IGV.
2. Load the genome:
   - `Genomes → Load Genome from Server`
   - Example: If the genome is *Plasmodium falciparum* (Pf3D7) then search for this and select **Plasmodium falciparum (Pf3D7)**.
3. Load the BED file:
   - `File → Load from File`
   - Select the generated `.bed` file.

The genome must be loaded before the BED file so that chromosome names (e.g. `Pf3D7_01_v3`) match.


> ⚠️ Note: BED files use **0-based, half-open coordinates** (BED standard).
           

## HTML QC Report

The `qc_report.html` provides:
- Overall transcript statistics
- QC flag distribution
- Examples of flagged transcripts
- Full provenance (tool version, inputs, parameters)
