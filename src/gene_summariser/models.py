from dataclasses import dataclass, field
from typing import Optional

@dataclass
class Feature:
    """Base class for genomic features."""
    
    seqid: str
    start: int
    end: int
    strand: str
    attributes: dict[str, str] = field(default_factory=dict)
    
    @property
    def length(self) -> int:
        """Return the length of the feature."""
        return self.end - self.start + 1
    
@dataclass
class Exon(Feature):
    """Represents an exon."""
    
    exon_id: str = ""
    
    def overlaps(self, other: "Exon") -> bool:
        """Check if this exon overlaps with another exon on the same strand."""
        if self.seqid != other.seqid or self.strand != other.strand:
            return False
        return not (self.end < other.start or self.start > other.end)
    
@dataclass
class CDS(Feature):
    """Represents a coding sequence."""
    
    phase: int = 0

@dataclass
class Transcript:
    """Represents a transcript with exons and optional CDS."""
    
    transcript_id: str
    gene_id: str
    seqid: str
    start: int
    end: int
    strand: str
    exons: list[Exon] = field(default_factory=list)
    cds_features: list[CDS] = field(default_factory=list)
    attributes: dict[str, str] = field(default_factory=dict)
    
    @property
    def n_exons(self) -> int:
        """Return the number of exons."""
        return len(self.exons)
    
    @property
    def has_cds(self) -> bool:
        """Return True if the transcript has CDS features."""
        return len(self.cds_features) > 0
    
    @property
    def total_exon_length(self) -> int:
        """Return the total length of all exons."""
        return sum(exon.length for exon in self.exons)
    
    @property
    def total_cds_length(self) -> int:
        """Return the total length of all CDS features."""
        return sum(cds.length for cds in self.cds_features)



