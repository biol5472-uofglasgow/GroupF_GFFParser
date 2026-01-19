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
