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