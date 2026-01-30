# Docker run script for gene-summariser (PowerShell)
# Usage: .\docker-run-test.ps1

Write-Host "Running Docker container with test data..." -ForegroundColor Green

docker run --rm `
    -v "${PWD}/tests/fixtures:/data" `
    -v "${PWD}/results:/app/results" `
    gene-summariser:latest `
    --gff /data/models.gff3 `
    --fasta /data/testfasta.fasta `
    --outdir /app/results

if ($LASTEXITCODE -eq 0) {
    Write-Host "`nSuccess! Results saved to .\results\" -ForegroundColor Green
} else {
    Write-Host "`nError: Docker run failed. Make sure the image is built first." -ForegroundColor Red
    Write-Host "Run: docker build -t gene-summariser:latest ." -ForegroundColor Yellow
}