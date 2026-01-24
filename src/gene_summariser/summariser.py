from gene_summariser.models import Transcript, TranscriptSummary

def make_summary(transcript: Transcript) -> TranscriptSummary:
    """
    Creates a TranscriptSummary object to be used in qc and output generation.
    Args:
        transcript (Transcript): A Transcript object.
    Returns:
        TranscriptSummary: A TranscriptSummary object.
    """

    return TranscriptSummary(
        gene_id=transcript.gene_id,
        transcript_id=transcript.transcript_id,
        n_exons=transcript.n_exons,
        has_cds=transcript.has_cds,
        total_exon_length=transcript.total_exon_length,
        total_cds_length=transcript.total_cds_length if transcript.has_cds else None,
    )