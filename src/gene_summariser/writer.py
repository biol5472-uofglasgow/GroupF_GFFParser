"""Output writers for various file formats."""

import json
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from gene_summariser.models import TranscriptSummary, Transcript

class OutputWriter:
    """Handles writing outputs to various formats."""

    def __init__(self, output_dir: Path) -> None:
        """Initialize output writer.
        
        Args:
            output_dir: Directory for output files
        """
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def write_transcript_summary(self, summaries: list[TranscriptSummary]) -> Path:
        """Write transcript summary to TSV file.
        
        Args:
            summaries: List of transcript summaries
            
        Returns:
            Path to the output file
        """
        # Convert to DataFrame for easy TSV writing
        data = []
        for summary in summaries:
            data.append({
                "gene_id": summary.gene_id,
                "transcript_id": summary.transcript_id,
                "n_exons": summary.n_exons,
                "has_cds": summary.has_cds,
                "flags": summary.flags_str,
            })
        
        df = pd.DataFrame(data)
        
        # Sort for deterministic output
        df = df.sort_values(["gene_id", "transcript_id"]).reset_index(drop=True)
        
        output_path = self.output_dir / "transcript_summary.tsv"
        df.to_csv(output_path, sep="\t", index=False)
        
        return output_path