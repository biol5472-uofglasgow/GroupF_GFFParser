@echo off
REM ================================
REM run_analysis.bat
REM Usage: run_analysis.bat <GFF_FILE> [FASTA_FILE] [OUT_DIR] [--strict]
REM ================================

SET GFF_FILE=%1
SET FASTA_FILE=
SET OUT_DIR=
SET STRICT_FLAG=

REM Determine strict flag
IF "%~nx4"=="--strict" SET STRICT_FLAG=--strict

REM Determine FASTA and OUT_DIR based on arguments
IF "%3"=="" (
    REM Only 2 arguments → GFF and OUT_DIR
    SET OUT_DIR=%2
) ELSE (
    REM 3 arguments → FASTA provided
    SET FASTA_FILE=%2
    SET OUT_DIR=%3
)


REM Create output folder if missing
IF NOT EXIST "%OUT_DIR%" mkdir "%OUT_DIR%"

REM Get directory from GFF for Docker mount
SET INPUT_DIR=%~dp1
SET INPUT_DIR=%INPUT_DIR:~0,-1%
SET INPUT_DIR=%INPUT_DIR:\=/%

REM Print info
echo GFF file: %GFF_FILE%
echo FASTA file: %FASTA_FILE%
echo Output directory: %OUT_DIR%
echo Mounted input directory: %INPUT_DIR%
echo Passing strict flag: %STRICT_FLAG%

REM Build Docker command
SET DOCKER_CMD=docker run --rm -v "%INPUT_DIR%:/work" -w /work gene-summariser -g "%~nx1" --outdir "%OUT_DIR%" %STRICT_FLAG%

REM Add FASTA if provided
IF NOT "%FASTA_FILE%"=="" SET DOCKER_CMD=%DOCKER_CMD% -f "%~nx2"

REM Run Docker
%DOCKER_CMD%