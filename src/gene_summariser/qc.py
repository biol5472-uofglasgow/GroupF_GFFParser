"""Quality control checks for transcript annotations."""

from typing import Callable

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
        max_exon_count: int = 100,
        min_cds_length: int = 30,
        max_exon_length: int = 1000000,
    ) -> None:
        """Initialize QC checker with configurable thresholds.
        
        Args:
            max_exon_count: Maximum reasonable number of exons (default: 100)
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