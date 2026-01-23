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