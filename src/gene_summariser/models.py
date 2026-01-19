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




