from Bio import SeqIO

from gene_summariser.models import Transcript


def extract_gene_sequences(fasta_file: str, transcript: Transcript) -> str:
    """
    Extracts gene sequences from a FASTA file from a Transcript object.

    Args:
        fasta_file (str): Path to the FASTA file.
        Transcript (Transcript): Transcript object containing gene information.

    Returns:
        str: The extracted gene sequence.
    """
    cds_sequence = []
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    for cds in transcript.cds_features:
        record = genome[cds.seqid]
        start = cds.start - 1
        end = cds.end
        sequence = record.seq[start:end]

        if cds.strand == "-":
            sequence = sequence.reverse_complement()

        if cds.phase is not None:
            sequence = sequence[cds.phase :]

    cds_sequence.append(str(sequence))
    cds_sequence = "".join(str(seq) for seq in cds_sequence)

    return cds_sequence
