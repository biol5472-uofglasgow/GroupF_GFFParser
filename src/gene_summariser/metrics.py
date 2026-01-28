from gene_summariser.models import Transcript, TranscriptSummary
from gene_summariser.qc import QCChecker


class MetricsCalculator:
    """Calculate summary metrics and QC flags for transcripts.

    This class coordinates between QC checking and summary creation,
    producing TranscriptSummary objects that include both metrics
    and quality control flags.

    Attributes:
        qc_checker: QCChecker instance for running quality checks
    """

    def __init__(self, qc_checker: QCChecker | None = None) -> None:
        """Initialize metrics calculator.

        Args:
            qc_checker: QC checker instance. If None, creates a default checker.

        Example:
            >>> # Use default QC settings
            >>> calculator = MetricsCalculator()
            >>>
            >>> # Use custom QC settings
            >>> strict_checker = QCChecker(max_exon_count=50, min_cds_length=100)
            >>> calculator = MetricsCalculator(strict_checker)
        """
        self.qc_checker = qc_checker or QCChecker()

    def calculate_transcript_summary(self, transcript: Transcript) -> TranscriptSummary:
        """Calculate summary statistics for a single transcript.

        This method:
        1. Runs all QC checks on the transcript
        2. Extracts key metrics (exon count, CDS presence, lengths)
        3. Creates a TranscriptSummary object with all information

        Args:
            transcript: Transcript to summarize

        Returns:
            TranscriptSummary object with metrics and QC flags

        Example:
            >>> calculator = MetricsCalculator()
            >>> transcript = parser.parse_transcript(...)
            >>> summary = calculator.calculate_transcript_summary(transcript)
            >>> print(f"{summary.transcript_id}: {summary.n_exons} exons, "
            ...       f"flags: {summary.flags_str}")
            t001: 3 exons, flags:
            t002: 1 exon, flags: no_cds,single_exon
        """
        # Run QC checks
        flags = self.qc_checker.check_transcript(transcript)

        # Create summary with metrics and flags
        return TranscriptSummary(
            gene_id=transcript.gene_id,
            transcript_id=transcript.transcript_id,
            n_exons=transcript.n_exons,
            has_cds=transcript.has_cds,
            flags=flags,
            total_exon_length=transcript.total_exon_length,
            total_cds_length=(
                transcript.total_cds_length if transcript.has_cds else None
            ),
        )

    def calculate_summaries(
        self, transcripts: list[Transcript]
    ) -> list[TranscriptSummary]:
        """Calculate summaries for multiple transcripts.

        Convenience method for batch processing of transcripts.

        Args:
            transcripts: List of transcripts to summarize

        Returns:
            List of transcript summaries in the same order as input

        Example:
            >>> calculator = MetricsCalculator()
            >>> transcripts = [t1, t2, t3]
            >>> summaries = calculator.calculate_summaries(transcripts)
            >>> flagged = [s for s in summaries if s.flags]
            >>> print(f"Found {len(flagged)} transcripts with issues")
        """
        return [self.calculate_transcript_summary(t) for t in transcripts]
