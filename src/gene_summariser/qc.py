"""Quality control checks for transcript annotations."""

from collections.abc import Callable

from gene_summariser.fasta import get_full_sequence
from gene_summariser.models import Transcript


class QCChecker:
    """Performs quality control checks on transcripts.

    This class implements a comprehensive set of QC checks for gene annotations,
    including both structural validation and biological plausibility checks.

    Attributes:
        max_exon_count: Maximum reasonable number of exons before flagging
        min_cds_length: Minimum CDS length in base pairs
        max_exon_length: Maximum reasonable exon length in base pairs
    """

    def __init__(
        self,
        fasta_file: str,
        max_exon_count: int = 50,
        min_cds_length: int = 30,
        max_exon_length: int = 1000000,
    ) -> None:
        self.fasta_file = fasta_file
        """Initialize QC checker with configurable thresholds.

        Args:
            max_exon_count: Maximum reasonable number of exons (default: 50)
            min_cds_length: Minimum CDS length in base pairs (default: 30)
            max_exon_length: Maximum reasonable exon length (default: 1,000,000)

        Example:
            >>> # Strict checking
            >>> checker = QCChecker(max_exon_count=50, min_cds_length=100)
            >>>
            >>> # Lenient checking
            >>> checker = QCChecker(max_exon_count=200, min_cds_length=10)
        """
        self.max_exon_count = max_exon_count
        self.min_cds_length = min_cds_length
        self.max_exon_length = max_exon_length

        # Register all QC checks - automatically run by check_transcript()
        self._checks: list[Callable[[Transcript], str | None]] = [
            # Structural checks
            self.check_no_exons,
            self.check_single_exon,
            self.check_high_exon_count,
            self.check_long_exon,
            self.check_overlapping_exons,
            # CDS checks
            self.check_no_cds,
            self.check_short_cds,
            self.check_cds_not_divisible_by_3,
            self.check_cds_phase_inconsistent,
            # Gene checks
            self.check_stop_codon,
            self.check_start_codon,
        ]

    def check_transcript(self, transcript: Transcript) -> list[str]:
        """Run all registered QC checks on a transcript.

        This method automatically runs all checks registered in self._checks
        and collects any flags that are returned.

        Args:
            transcript: Transcript to check

        Returns:
            List of QC flag strings. Empty list if no issues found.

        Example:
            >>> checker = QCChecker()
            >>> transcript = parser.parse_transcript(...)
            >>> flags = checker.check_transcript(transcript)
            >>> if flags:
            ...     print(f"Found {len(flags)} issues: {', '.join(flags)}")
            Found 2 issues: no_cds, single_exon
        """
        flags: list[str] = []

        for check in self._checks:
            flag = check(transcript)
            if flag:
                flags.append(flag)

        return flags

    # ========================================================================
    # STRUCTURAL CHECKS
    # ========================================================================

    def check_no_exons(self, transcript: Transcript) -> str | None:
        """Check if transcript has no exons (critical annotation error).

        A transcript without exons indicates a severe annotation problem.
        This should never occur in properly formatted GFF3 files.

        Args:
            transcript: Transcript to check

        Returns:
            "no_exons" if transcript has no exons, None otherwise
        """
        if transcript.n_exons == 0:
            return "no_exons"
        return None

    def check_single_exon(self, transcript: Transcript) -> str | None:
        """Check if transcript has only one exon.

        Single-exon transcripts are unusual in eukaryotes (though valid for
        some genes like histones). This flag helps identify potential
        annotation issues or biologically interesting features.

        Args:
            transcript: Transcript to check

        Returns:
            "single_exon" if transcript has exactly one exon, None otherwise
        """
        if transcript.n_exons == 1:
            return "single_exon"
        return None

    def check_high_exon_count(self, transcript: Transcript) -> str | None:
        """Check if transcript has an unusually high number of exons.

        Transcripts with extremely high exon counts may indicate annotation
        errors such as fragmented exons or merged gene models.

        Args:
            transcript: Transcript to check

        Returns:
            "high_exon_count" if exon count exceeds threshold, None otherwise
        """
        if transcript.n_exons > self.max_exon_count:
            return "high_exon_count"
        return None

    def check_long_exon(self, transcript: Transcript) -> str | None:
        """Check if any exon is unusually long.

        Extremely long exons may indicate unannotated introns or merged
        exon features. Typical exons are a few hundred bp, with most
        under 10kb.

        Args:
            transcript: Transcript to check

        Returns:
            "long_exon" if any exon exceeds length threshold, None otherwise
        """
        for exon in transcript.exons:
            if exon.length > self.max_exon_length:
                return "long_exon"
        return None

    def check_overlapping_exons(self, transcript: Transcript) -> str | None:
        """Check if any exons overlap (critical annotation error).

        Overlapping exons within the same transcript indicate a serious
        annotation error. Exons should be separated by introns.

        Note: This only checks exons on the same strand. Different
        transcripts or genes can have overlapping exons.

        Args:
            transcript: Transcript to check

        Returns:
            "overlapping_exons" if any exons overlap, None otherwise
        """
        # Sort exons by start position for efficient overlap checking
        exons = sorted(transcript.exons, key=lambda e: e.start)

        # Check each adjacent pair
        for i in range(len(exons) - 1):
            if exons[i].overlaps(exons[i + 1]):
                return "overlapping_exons"

        return None

    # ========================================================================
    # CDS CHECKS
    # ========================================================================

    def check_no_cds(self, transcript: Transcript) -> str | None:
        """Check if transcript has no coding sequence.

        Transcripts without CDS are typically non-coding RNAs (ncRNAs,
        lncRNAs, etc.). This flag helps distinguish coding from non-coding
        transcripts.

        Args:
            transcript: Transcript to check

        Returns:
            "no_cds" if transcript lacks CDS features, None otherwise
        """
        if not transcript.has_cds:
            return "no_cds"
        return None

    def check_short_cds(self, transcript: Transcript) -> str | None:
        """Check if CDS is unusually short.

        Very short CDS regions may indicate incomplete annotations,
        pseudogenes, or annotation errors. Typical proteins are at least
        50-100 amino acids (150-300 bp).

        Args:
            transcript: Transcript to check

        Returns:
            "short_cds" if CDS length is below threshold, None otherwise
        """
        if transcript.has_cds and transcript.total_cds_length < self.min_cds_length:
            return "short_cds"
        return None

    def check_cds_not_divisible_by_3(self, transcript: Transcript) -> str | None:
        """Check if CDS length is not divisible by 3.

        CDS length must be divisible by 3 to encode complete codons.
        A length not divisible by 3 indicates:
        - Frameshift mutation
        - Incomplete annotation
        - Annotation error
        - Premature stop codon

        This is a critical biological validation check.

        Args:
            transcript: Transcript to check

        Returns:
            "cds_not_divisible_by_3" if CDS length not divisible by 3,
            None otherwise
        """
        if not transcript.has_cds:
            return None

        if transcript.total_cds_length % 3 != 0:
            return "cds_not_divisible_by_3"

        return None

    def check_cds_phase_inconsistent(self, transcript: Transcript) -> str | None:
        """Check if CDS phases are inconsistent across multiple CDS features.

        In a properly annotated transcript with multiple CDS features
        (split by introns), the phase of each CDS must be consistent with
        the cumulative length of previous CDS features.

        Phase calculation:
            next_phase = (current_phase + current_length) % 3

        Inconsistent phases indicate:
        - Annotation errors in phase assignment
        - Incorrectly merged/split CDS features
        - Errors in intron/exon boundaries

        This is an advanced biological validation check that ensures
        proper reading frame maintenance across introns.

        Args:
            transcript: Transcript to check

        Returns:
            "cds_phase_inconsistent" if phases are inconsistent,
            None otherwise

        Example:
            CDS 1: phase=0, length=100  → next expected phase = (0+100)%3 = 1
            CDS 2: phase=1 ✓ (correct)
            CDS 3: phase=0 ✗ (inconsistent - should be 2)
        """
        cds_list = transcript.cds_features

        # Need at least 2 CDS features to check consistency
        if len(cds_list) < 2:
            return None

        # Sort CDS features by genomic position (accounting for strand direction)
        if transcript.strand == "+":
            # Forward strand: sort by start position (5' to 3')
            sorted_cds = sorted(cds_list, key=lambda c: c.start)
        else:
            # Reverse strand: sort by start position descending (5' to 3' in reverse)
            sorted_cds = sorted(cds_list, key=lambda c: c.start, reverse=True)

        # Check phase consistency between adjacent CDS features
        for i in range(len(sorted_cds) - 1):
            current_cds = sorted_cds[i]
            next_cds = sorted_cds[i + 1]

            # Calculate expected phase based on current CDS length
            cds_length = current_cds.length
            expected_next_phase = (current_cds.phase + cds_length) % 3

            # Check if actual phase matches expected phase
            if next_cds.phase != expected_next_phase:
                return "cds_phase_inconsistent"

        return None

    def check_start_codon(self, transcript: Transcript) -> str | None:
        if not self.fasta_file:
            return None
        start_codons = {"ATG", "GTG", "TTG"}
        if transcript.has_cds:
            full_sequence = get_full_sequence(self.fasta_file, transcript)
            start = full_sequence[:3].upper()
            if start not in start_codons:
                return "missing_start_codon"
        return None

    def check_stop_codon(self, transcript: Transcript) -> str | None:
        if not self.fasta_file:
                return None
        stop_codons = {"TAA", "TAG", "TGA"}
        if transcript.has_cds:
            full_sequence = get_full_sequence(self.fasta_file, transcript)
            stop = full_sequence[-3:].upper()
            if stop not in stop_codons:
                return "missing_stop_codon"
        return None
