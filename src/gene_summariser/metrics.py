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