# Testing Guide

Comprehensive guide for testing the Gene Model Summariser.

## Test Structure

```
tests/
├── fixtures/              # Test data files
│   ├── simple.gff3       # Basic test case
│   ├── complex.gff3      # Multiple transcripts
│   ├── edge_cases.gff3   # Problematic annotations
│   └── genome.fasta      # Test genome sequence
├── test_models.py         # Data model tests
├── test_parser.py         # Parser tests
├── test_qc.py            # QC check tests
├── test_metrics.py       # Metrics tests
├── test_writer.py        # Output writer tests
├── test_cli.py           # CLI tests
├── test_integration.py   # End-to-end tests
└── conftest.py           # Pytest configuration
```

## Running Tests

### Basic Test Execution

```bash
# Run all tests
pytest

# Verbose mode
pytest -v

# Run specific test file
pytest tests/test_qc.py

# Run specific test
pytest tests/test_qc.py::test_check_no_cds

# Stop at first failure
pytest -x

# Show local variables on failure
pytest -l
```

### Coverage

```bash
# Run with coverage
pytest --cov=gene_summariser

# HTML coverage report
pytest --cov=gene_summariser --cov-report=html

# Missing lines report
pytest --cov=gene_summariser --cov-report=term-missing

# Specific module coverage
pytest --cov=gene_summariser.qc --cov-report=term
```

### Test Selection

```bash
# Run tests matching pattern
pytest -k "test_check"

# Run only failed tests from last run
pytest --lf

# Run failed tests first
pytest --ff

# Run tests in parallel (install pytest-xdist)
pytest -n auto
```

## Test Categories

### Unit Tests

Test individual components in isolation.

**Example: QC Check Test**

```python
def test_check_no_cds():
    """Test that transcripts without CDS get flagged."""
    checker = QCChecker(fasta_file=None)
    transcript = Transcript(
        transcript_id="T1",
        gene_id="G1",
        seqid="chr1",
        start=100,
        end=500,
        strand="+",
    )
    # No CDS added

    flags = checker.check_transcript(transcript)
    assert "no_cds" in flags
```

### Integration Tests

Test complete workflows.

**Example: End-to-End Test**

```python
def test_complete_pipeline(tmp_path):
    """Test entire pipeline from GFF3 to output."""
    # Create test GFF3
    gff_path = tmp_path / "test.gff3"
    gff_path.write_text(TEST_GFF3_CONTENT)

    # Parse
    parser = ParserGFF(str(gff_path))
    transcripts = parser.parse_transcripts()

    # QC
    checker = QCChecker(fasta_file=None)
    calculator = MetricsCalculator(checker)
    summaries = calculator.calculate_summaries(transcripts)

    # Write
    writer = OutputWriter(tmp_path / "output")
    tsv_path = writer.write_transcript_summary(summaries)

    # Verify
    assert tsv_path.exists()
    df = pd.read_csv(tsv_path, sep="\t")
    assert len(df) > 0
```

### Parametrized Tests

Test multiple scenarios efficiently.

```python
@pytest.mark.parametrize("n_exons,expected_flag", [
    (0, "no_exons"),
    (1, "single_exon"),
    (101, "high_exon_count"),
])
def test_exon_count_flags(n_exons, expected_flag):
    """Test exon count QC flags."""
    checker = QCChecker(fasta_file=None, max_exon_count=100)
    transcript = create_transcript_with_n_exons(n_exons)

    flags = checker.check_transcript(transcript)
    assert expected_flag in flags
```

## Test Fixtures

### Creating Test Data

**conftest.py**

```python
import pytest
from pathlib import Path

@pytest.fixture
def test_gff3(tmp_path):
    """Create a simple test GFF3 file."""
    gff_content = """##gff-version 3
chr1	test	gene	1000	5000	.	+	.	ID=gene001
chr1	test	mRNA	1000	5000	.	+	.	ID=t001;Parent=gene001
chr1	test	exon	1000	1200	.	+	.	ID=e001;Parent=t001
chr1	test	CDS	1100	1200	.	+	0	ID=c001;Parent=t001
"""
    gff_path = tmp_path / "test.gff3"
    gff_path.write_text(gff_content)
    return gff_path

@pytest.fixture
def test_fasta(tmp_path):
    """Create a simple test FASTA file."""
    fasta_content = """>chr1
ATGGCCTACGGATAATAGACC
"""
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_content)
    return fasta_path
```

### Using Fixtures

```python
def test_with_fixtures(test_gff3, test_fasta):
    """Test using fixture data."""
    parser = ParserGFF(str(test_gff3))
    transcripts = parser.parse_transcripts()

    checker = QCChecker(fasta_file=str(test_fasta))
    flags = checker.check_transcript(transcripts[0])

    assert len(transcripts) > 0
```

## Testing Strategies

### Testing QC Checks

```python
class TestQCChecker:
    """Test suite for QC checks."""

    @pytest.fixture
    def checker(self):
        """Create QC checker instance."""
        return QCChecker(fasta_file=None)

    def test_no_cds(self, checker):
        """Test no_cds flag."""
        transcript = create_transcript(has_cds=False)
        flags = checker.check_transcript(transcript)
        assert "no_cds" in flags

    def test_with_cds(self, checker):
        """Test transcript with CDS."""
        transcript = create_transcript(has_cds=True)
        flags = checker.check_transcript(transcript)
        assert "no_cds" not in flags
```

### Testing Parser

```python
def test_parser_basic(test_gff3):
    """Test basic parsing."""
    parser = ParserGFF(str(test_gff3))
    transcripts = parser.parse_transcripts()

    assert len(transcripts) > 0
    assert transcripts[0].transcript_id is not None

def test_parser_gene_structure(test_gff3):
    """Test gene structure parsing."""
    parser = ParserGFF(str(test_gff3))
    genes = list(parser.parse_genes())

    assert len(genes) > 0
    assert genes[0].n_transcripts > 0
    assert genes[0].transcripts[0].n_exons > 0
```

### Testing CLI

```python
def test_cli_help():
    """Test CLI help command."""
    result = subprocess.run(
        ["gene-summariser", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "usage" in result.stdout.lower()

def test_cli_execution(test_gff3, tmp_path):
    """Test CLI execution."""
    outdir = tmp_path / "output"

    result = subprocess.run([
        "gene-summariser",
        "--gff", str(test_gff3),
        "--outdir", str(outdir)
    ], capture_output=True, text=True)

    assert result.returncode == 0
    assert (outdir / "transcript_summary.tsv").exists()
```

## Test Data Examples

### Simple GFF3

```gff3
##gff-version 3
chr1	test	gene	1000	5000	.	+	.	ID=gene001
chr1	test	mRNA	1000	5000	.	+	.	ID=t001;Parent=gene001
chr1	test	exon	1000	1200	.	+	.	ID=e001;Parent=t001
chr1	test	exon	2000	2300	.	+	.	ID=e002;Parent=t001
chr1	test	CDS	1100	1200	.	+	0	ID=c001;Parent=t001
chr1	test	CDS	2000	2300	.	+	2	ID=c002;Parent=t001
```

### Edge Cases GFF3

```gff3
##gff-version 3
# Transcript with no CDS
chr1	test	gene	1000	2000	.	+	.	ID=gene001
chr1	test	ncRNA	1000	2000	.	+	.	ID=t001;Parent=gene001
chr1	test	exon	1000	2000	.	+	.	ID=e001;Parent=t001

# Single exon transcript
chr1	test	gene	3000	4000	.	+	.	ID=gene002
chr1	test	mRNA	3000	4000	.	+	.	ID=t002;Parent=gene002
chr1	test	exon	3000	4000	.	+	.	ID=e002;Parent=t002
chr1	test	CDS	3100	3900	.	+	0	ID=c002;Parent=t002

# Overlapping exons (error)
chr1	test	gene	5000	6000	.	+	.	ID=gene003
chr1	test	mRNA	5000	6000	.	+	.	ID=t003;Parent=gene003
chr1	test	exon	5000	5100	.	+	.	ID=e003;Parent=t003
chr1	test	exon	5050	5150	.	+	.	ID=e004;Parent=t003
```

## Continuous Integration

### GitHub Actions

Tests run automatically on:

- Push to main/develop
- Pull requests
- Release creation

View results:

```bash
# In your repository
https://github.com/biol5472-uofglasgow/your-repo/actions
```

### Local CI Simulation

```bash
# Run all CI checks locally
make check

# Or manually
pytest --cov=gene_summariser
ruff check src/ tests/
mypy src/
```

## Debugging Tests

### Using pdb

```python
def test_with_debugging():
    """Test with debugger."""
    transcript = create_transcript()

    import pdb; pdb.set_trace()  # Debugger starts here

    flags = checker.check_transcript(transcript)
    assert "no_cds" in flags
```

### Using pytest debugging

```bash
# Drop into debugger on failure
pytest --pdb

# Drop into debugger at start of test
pytest --trace
```

### Verbose Output

```bash
# Show print statements
pytest -s

# Very verbose
pytest -vv

# Show local variables on failure
pytest -l --tb=long
```

## Performance Testing

### Timing Tests

```python
def test_parser_performance(benchmark, large_gff3):
    """Benchmark parser performance."""
    def parse():
        parser = ParserGFF(str(large_gff3))
        return parser.parse_transcripts()

    result = benchmark(parse)
    assert len(result) > 0
```

### Memory Profiling

```python
import tracemalloc

def test_memory_usage(test_gff3):
    """Test memory usage."""
    tracemalloc.start()

    parser = ParserGFF(str(test_gff3))
    transcripts = parser.parse_transcripts()

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Assert reasonable memory usage
    assert peak < 100 * 1024 * 1024  # 100 MB
```

## Coverage Goals

### Target Coverage

- **Overall**: > 90%
- **Core modules** (parser, qc, metrics): > 95%
- **CLI**: > 80%
- **Models**: 100% (simple dataclasses)

### Checking Coverage

```bash
# Generate coverage report
pytest --cov=gene_summariser --cov-report=html

# View in browser
open htmlcov/index.html

# Find untested lines
pytest --cov=gene_summariser --cov-report=term-missing
```

## Best Practices

### 1. Test Naming

```python
# ✓ Good: Descriptive
def test_check_no_cds_returns_flag_when_transcript_has_no_cds():

# ✓ Acceptable: Concise
def test_no_cds_flag():

# ✗ Bad: Unclear
def test_1():
```

### 2. Test Organization

```python
# Group related tests
class TestQCChecker:
    class TestStructuralChecks:
        def test_no_exons(self):
            pass

        def test_single_exon(self):
            pass

    class TestCDSChecks:
        def test_no_cds(self):
            pass
```

### 3. Arrange-Act-Assert Pattern

```python
def test_example():
    # Arrange: Set up test data
    transcript = create_transcript()
    checker = QCChecker(fasta_file=None)

    # Act: Execute the code under test
    flags = checker.check_transcript(transcript)

    # Assert: Verify the results
    assert "no_cds" in flags
```

### 4. Use Fixtures

```python
# ✓ Good: Reusable fixture
@pytest.fixture
def qc_checker():
    return QCChecker(fasta_file=None)

def test_something(qc_checker):
    flags = qc_checker.check_transcript(transcript)
    assert len(flags) > 0

# ✗ Bad: Repeated setup
def test_something():
    checker = QCChecker(fasta_file=None)  # Repeated
    # ...
```

### 5. Test Edge Cases

```python
# Test normal cases
def test_normal_transcript():
    pass

# Test edge cases
def test_empty_transcript():
    pass

def test_single_exon():
    pass

def test_many_exons():
    pass

def test_negative_coordinates():
    pass
```

## Troubleshooting

### Tests Won't Run

```bash
# Check pytest installation
pytest --version

# Reinstall if needed
pip install -e ".[dev]"

# Check Python path
python -c "import gene_summariser"
```

### Import Errors

```bash
# Ensure package is installed
pip install -e .

# Check PYTHONPATH
export PYTHONPATH="${PYTHONPATH}:$(pwd)/src"
```

### Fixture Not Found

```python
# Ensure conftest.py is in tests directory
# Or import fixture explicitly
from .conftest import test_gff3
```

## Additional Resources

- [pytest Documentation](https://docs.pytest.org/)
- [pytest-cov Documentation](https://pytest-cov.readthedocs.io/)
- [Testing Best Practices](https://docs.python-guide.org/writing/tests/)
- [GFF3 Specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

---

**Remember**: Good tests are documentation. Write tests that explain what the code should do.
