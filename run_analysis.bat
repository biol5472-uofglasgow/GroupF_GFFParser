@echo off
REM ================================
REM run_analysis.bat
REM Usage: run_analysis.bat <GFF_FILE> <FASTA_FILE> [OUT_DIR] [--strict]
REM ================================

SET GFF_FILE=%1
SET FASTA_FILE=%2
SET OUT_DIR=%3
SET STRICT_FLAG=

REM Default output folder
IF "%OUT_DIR%"=="" SET OUT_DIR=results

REM Make sure output directory exists
IF NOT EXIST "%OUT_DIR%" mkdir "%OUT_DIR%"

REM Check if last argument is --strict
SET STRICT_FLAG=
IF "%~nx4"=="--strict" SET STRICT_FLAG=--strict

REM Check GFF file exists, exits if it does not
IF NOT EXIST "%GFF_FILE%" (
    echo ERROR: GFF file not found: %GFF_FILE%
    exit /b 1
)

REM Check FASTA file exists, warns if it does not
IF NOT EXIST "%FASTA_FILE%" (
    echo WARNING: FASTA file not found: %FASTA_FILE%
)

REM This gets the directory from the gff file path, this is used to mount the directory in docker
SET INPUT_DIR=%~dp1
SET INPUT_DIR=%INPUT_DIR:~0,-1%
SET INPUT_DIR=%INPUT_DIR:\=/%

REM Print logging info like imput files and output dir as well as strict flag status
echo GFF file: %GFF_FILE%
echo FASTA file: %FASTA_FILE%
echo Output directory: %OUT_DIR%
echo Mounted input directory: %INPUT_DIR%
echo Passing strict flag: %STRICT_FLAG%

REM Run Docker
docker run --rm -v "%INPUT_DIR%:/work" -w /work gene-summariser -g "%~nx1" -f "%~nx2" --outdir "%OUT_DIR%" %STRICT_FLAG%