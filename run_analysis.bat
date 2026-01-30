@echo off
REM Usage: docker-run-test.bat <gff_file> <fasta_file> <outdir>
REM Must be run from the project root

SET GFF_FILE=%1
SET FASTA_FILE=%2
SET OUT_DIR=%3

IF "%OUT_DIR%"=="" SET OUT_DIR=results

IF NOT EXIST "%OUT_DIR%" mkdir "%OUT_DIR%"

docker run --rm -v %CD%:/work -w /work gene-summariser:latest --gff %GFF_FILE% --fasta %FASTA_FILE% --outdir %OUT_DIR%
