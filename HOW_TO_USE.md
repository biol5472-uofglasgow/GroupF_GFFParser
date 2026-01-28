# Gene Summariser

## Quick Start

### Standard Installation & Execution

```bash
# Install
pip install -e .

# Run basic analysis
gene-summariser --gff test/fixtures/models.gff3 --outdir results/
```

### With FASTA (for sequence validation)

```bash
gene-summariser \
  --gff test/fixtures/models.gff3 \
  --fasta test/fixtures/testfasta.fasta \
  --outdir results/
```

---

## Alternative Execution Methods

### Option A: Makefile

```bash
# Run analysis
make run INPUT=data/annotations.gff3 OUTPUT=results/
```

**Makefile target:**

```makefile
## run: Run analysis on specified input
run:
	@echo "Running gene-summariser..."
	gene-summariser --gff $(INPUT) --outdir $(OUTPUT)
```

### Option B: Docker

```bash
# Build image
docker build -t gene-summariser:latest .

# Run analysis
docker run --rm \
  -v $(pwd)/data:/data \
  -v $(pwd)/results:/app/results \
  gene-summariser:latest \
  --gff /data/annotations.gff3 \
  --outdir /app/results
```

For details, see [DOCKER_USAGE.md](DOCKER_USAGE.md)

### Option C: Shell Script Wrapper

```bash
# Make executable and run
chmod +x run_analysis.sh
./run_analysis.sh --gff data.gff3 --outdir results/
```

---

## Exploring Results

### Option 1: HTML Report

Generate an interactive HTML report:

```bash
python scripts/generate_report.py results/ -o qc_report.html

# Open in browser
open report.html         # macOS
xdg-open report.html     # Linux
start report.html        # Windows
```

**Report includes:**

- Summary statistics (total transcripts, flagged, etc.)
- QC flag breakdown
- Flagged transcripts table
- Full provenance information

### Option 2: Spreadsheet Analysis

```bash
# Open TSV in Excel/LibreOffice
libreoffice results/transcript_summary.tsv
```

### Option 3: Genome Browser

Load output files in IGV, UCSC Genome Browser, or similar tools:

- `results/qc_flags.gff3`
- `results/qc_flags.bed`

### Option 4: Interactive Jupyter Notebook

```bash
# Install dependencies
pip install jupyter matplotlib seaborn

# Start notebook
jupyter notebook notebooks/explore_results.ipynb
```

**Notebook features:**

- Interactive data exploration
- Statistical visualizations
- Custom analysis capabilities

### Option 5: Command-Line Summary

```bash
python scripts/summarize_results.py results/
```

---

## Documentation

- [Docker Usage Guide](DOCKER_USAGE.md)
- Additional documentation coming soon

---

## Support

For issues or questions, please open an issue on the repository.
