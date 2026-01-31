Commands to Run
Windows:
docker build -t gene-summariser:latest .
.\run_analysis.bat <gff_path> <fasta_path> <out_dir> <--strict>

Mac / linux:
docker build -t gene-summariser:latest .
.\run_analysis.sh <gff_path> <fasta_path> <out_dir> <--strict>

For both the path can either be an absolute path or a relative path from the root dir, in the second case all paths should start with .\
e.g. .\run_analysis.sh .\test\fixtures\models.gff3 .\test\fixtures\testfasta.fasta .\reusults\run1

fasta and --strict flags are optional, but gff and our_dir are required


Contributers:
Callum Poole
ABEL BIJU DANIEL
Nosakhare Osaro
