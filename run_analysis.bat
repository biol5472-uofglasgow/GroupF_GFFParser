@echo off
REM Minimal run_analysis.bat with spaces support
REM Usage: run_analysis.bat <GFF_FILE> <FASTA_FILE> <OUT_DIR>

SET GFF_FILE=%1
SET FASTA_FILE=%2
SET OUT_DIR=%3

IF "%OUT_DIR%"=="" SET OUT_DIR=results

REM Make sure output directory exists
IF NOT EXIST "%OUT_DIR%" mkdir "%OUT_DIR%"

REM Get folder of GFF file and convert backslashes to forward slashes
FOR %%F IN ("%GFF_FILE%") DO SET INPUT_DIR=%%~dpF
SET INPUT_DIR=%INPUT_DIR:~0,-1%
SET INPUT_DIR=%INPUT_DIR:\=/%

REM Print info
echo GFF file: %GFF_FILE%
echo FASTA file: %FASTA_FILE%
echo Output directory: %OUT_DIR%
echo Mounted input directory: %INPUT_DIR%

REM Run Docker
docker run --rm -v "%INPUT_DIR%:/work" -w /work gene-summariser -g "%~nx1" -f "%~nx2" --outdir "%OUT_DIR%"