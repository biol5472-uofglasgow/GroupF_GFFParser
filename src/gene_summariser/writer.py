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
    
    def write_provenance(
        self,
        input_file: Path,
        parameters: dict[str, Any],
        version: str = "0.1.0",
    ) -> Path:
        """Write provenance information to JSON.
        
        Args:
            input_file: Path to input GFF3 file
            parameters: Dictionary of parameters used
            version: Tool version
            
        Returns:
            Path to the output file
        """
        provenance = {
            "tool": "gene-summariser",
            "version": version,
            "timestamp": datetime.now().isoformat(),
            "inputs": {
                "gff3_file": str(input_file.absolute()),
            },
            "parameters": parameters,
        }
        
        output_path = self.output_dir / "run.json"
        with open(output_path, "w") as f:
            json.dump(provenance, f, indent=2)
        
        return output_path