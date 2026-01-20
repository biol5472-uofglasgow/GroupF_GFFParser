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
    
    def write_qc_flags_gff3(
        self,
        transcripts: list[Transcript],
        summaries: list[TranscriptSummary],
    ) -> Path:
        """Write QC flags to GFF3 format.
        
        Args:
            transcripts: Original transcript objects
            summaries: Transcript summaries with flags
            
        Returns:
            Path to the output file
        """
        # Create a mapping of transcript_id to flags
        flags_map = {s.transcript_id: s.flags_str for s in summaries if s.flags}
        
        output_path = self.output_dir / "qc_flags.gff3"
        
        with open(output_path, "w") as f:
            f.write("##gff-version 3\n")
            
            # Write only transcripts with flags
            for transcript in sorted(transcripts, key=lambda t: (t.seqid, t.start)):
                if transcript.transcript_id in flags_map:
                    flags = flags_map[transcript.transcript_id]
                    
                    # Build attributes string
                    attrs = [
                        f"ID={transcript.transcript_id}",
                        f"Parent={transcript.gene_id}",
                        f"qc_flags={flags}",
                    ]
                    attrs_str = ";".join(attrs)
                    
                    # Write GFF3 line
                    line = "\t".join([
                        transcript.seqid,
                        "gene_summariser",
                        "transcript",
                        str(transcript.start),
                        str(transcript.end),
                        ".",
                        transcript.strand,
                        ".",
                        attrs_str,
                    ])
                    f.write(line + "\n")
        
        return output_path
    
    def write_qc_flags_bed(
        self,
        transcripts: list[Transcript],
        summaries: list[TranscriptSummary],
    ) -> Path:
        """Write QC flags to BED format.
        
        Args:
            transcripts: Original transcript objects
            summaries: Transcript summaries with flags
            
        Returns:
            Path to the output file
        """
        # Create a mapping of transcript_id to flags
        flags_map = {s.transcript_id: s.flags_str for s in summaries if s.flags}
        
        output_path = self.output_dir / "qc_flags.bed"
        
        with open(output_path, "w") as f:
            # BED header
            f.write("# chrom\tstart\tend\tname\tscore\tstrand\n")
            
            # Write only transcripts with flags
            for transcript in sorted(transcripts, key=lambda t: (t.seqid, t.start)):
                if transcript.transcript_id in flags_map:
                    flags = flags_map[transcript.transcript_id]
                    
                    # BED uses 0-based coordinates
                    line = "\t".join([
                        transcript.seqid,
                        str(transcript.start - 1),  # Convert to 0-based
                        str(transcript.end),
                        f"{transcript.transcript_id}|{flags}",
                        "0",
                        transcript.strand,
                    ])
                    f.write(line + "\n")
        
        return output_path