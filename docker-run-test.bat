@echo off
REM Docker run script for gene-summariser on Windows
REM Usage: docker-run-test.bat

echo Running Docker container with test data...

docker run --rm ^
    -v "%cd%/tests/fixtures:/data" ^
    -v "%cd%/results:/app/results" ^
    gene-summariser:latest ^
    --gff /data/models.gff3 ^
    --fasta /data/testfasta.fasta ^
    --outdir /app/results

if %ERRORLEVEL% EQU 0 (
    echo.
    echo Success! Results saved to .\results\
) else (
    echo.
    echo Error: Docker run failed. Make sure the image is built first.
    echo Run: docker build -t gene-summariser:latest .
)