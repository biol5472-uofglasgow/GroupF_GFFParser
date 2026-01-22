from gene_summariser.praser import ParserGFF
from gene_summariser.summariser import make_summary
from gene_summariser.models import Transcript, TranscriptSummary

def flag_devisable_by_three(transcript: Transcript, summary: TranscriptSummary) -> TranscriptSummary:
    cds_length = transcript.total_cds_length

    if cds_length == 0:
        return summary
    if cds_length % 3 != 0:
        summary.flags.append("CDS_NOT_DEVISABLE_BY_3")

    return summary


