from gene_summariser.models import Transcript, TranscriptSummary
from gene_summariser.summariser import make_summary


def flag_devisable_by_three(transcript: Transcript, summary: TranscriptSummary) -> None:
    """
    Mutates the TranscriptSummary object to append a flag if the total CDS length is not divisible by 3.
    Args:
        transcript (Transcript): A Transcript object.
        summary (TranscriptSummary): A TranscriptSummary object to be mutated.
    Returns:
        None
    """
    cds_length = transcript.total_cds_length

    if cds_length == 0:
        return summary
    if cds_length % 3 != 0:
        summary.flags.append("CDS_NOT_DEVISABLE_BY_3")


# QC for phase incosistency
"""
This fucntion performs a quality check to ensure that all cds within a transcript maintain a 
consistent reading frame
"""


def check_cds_phase_errors(transcript: Transcript, summary: TranscriptSummary) -> None:

    cds_list = transcript.cds_features

    # if a transcript has fewer than 2 cds then it will skip it
    if len(cds_list) < 2:
        return
    # and for function with multiple cds it iterates over each consecutive pair of cds
    for i in range(len(cds_list) - 1):
        current_cds = cds_list[i]
        next_cds = cds_list[i + 1]

        cds_length = current_cds.length
        expected_next_phase = (
            current_cds.phase + cds_length
        ) % 3  # calculates the next cds phase

        if (
            next_cds.phase != expected_next_phase
        ):  # if the next phase not equal to expected phase throw an error
            summary.flags.append("CDS_Phase_Incosistent")
            break


def append_flags_to_summary(transcripts: list[Transcript]) -> list[TranscriptSummary]:
    """
    Creates a Transcript Summary for each Transcript and appends the QC flags.
    Args:
        transcripts (list[Transcript]): A list of Transcript objects.
    Returns:
        list[TranscriptSummary]: A list of TranscriptSummary objects with QC flags.
    """
    summaries = []
    for transcript in transcripts:
        summary = make_summary(transcript)
        flag_devisable_by_three(transcript, summary)
        summaries.append(summary)
    return summaries
