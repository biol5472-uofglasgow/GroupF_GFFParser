"""Gene model summariser package."""
# Module-level docstring describing the purpose of this package:
# it provides functionality to summarise gene models.

__version__ = "0.1.0"
# Defines the current version of the package as a string,
# useful for tracking releases and compatibility.

from gene_summariser.models import Gene, Transcript, Exon, CDS, TranscriptSummary
# Imports core data model classes representing biological entities
# (Gene, Transcript, Exon, CDS) and a TranscriptSummary helper class.


from gene_summariser.writer import OutputWriter
# Imports the OutputWriter class, used to format and write results
# to files or other output targets.

__all__ = [
    "Gene",
    # Exposes the Gene class as part of the public package API

    "Transcript",
    # Exposes the Transcript class for external imports

    "Exon",
    # Exposes the Exon class for external imports

    "CDS",
    # Exposes the CDS (Coding Sequence) class for external imports

    "TranscriptSummary",
    # Exposes the TranscriptSummary class for summarised transcript data

    "OutputWriter",
    # Exposes the OutputWriter for generating output files
]
# The __all__ list explicitly defines what symbols are exported
# when `from gene_summariser import *` is used.
