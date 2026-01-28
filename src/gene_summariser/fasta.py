from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from gene_summariser.models import Transcript


def get_full_sequence(genome: dict[str, SeqRecord], transcript: Transcript) -> str:
    """
    Extract the full spliced CDS (no phase),
    suitable for start/stop codon QC.
    Args:
        fasta_file (str): Path to the FASTA file.
        transcript (Transcript): Transcript object containing gene information.
    Returns:
        str: The full spliced CDS sequence.
    """
    seq_chunks: list[Seq] = []

    for cds in transcript.cds_features:
        record = genome[cds.seqid]

        seq=record.seq
        if seq is None:
            raise ValueError(f"No sequence found for {cds.seqid}")

        seq_chunks.append(seq[cds.start - 1 : cds.end])

    full_sequence = sum(seq_chunks, Seq(""))

    if transcript.strand == "-":
        full_sequence = full_sequence.reverse_complement() # type: ignore[no-untyped-call]

    return str(full_sequence)
