"""Quality control checks for transcript annotations."""

from typing import Callable

from gene_summariser.models import Transcript

class QCChecker:
    """Performs quality control checks on transcripts."""