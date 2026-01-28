# Parser module for GFF files. parser.py
import logging
import os
from collections.abc import Iterator

import gffutils

from gene_summariser.models import CDS, Exon, Gene, Transcript

logger = logging.getLogger(__name__)


class ParserGFF:
    """
    ParserGFF provides the necessary methods to parse a specified GFF file,
    extracting Gene, transcript, exon, and CDS information.

    Attributes:
        gff_path (str): Path to the GFF file to be parsed.
        db (gffutils.FeatureDB): GFF database object for efficient querying of features.

    Methods:
        parse_transcripts() -> list[Transcript]: Parses the GFF file and returns
            a list of Transcript objects, each containing associated exons and CDS
            features.
        parse_genes() -> Iterator[Gene]: Parses the GFF file and yields Gene objects.
    """

    allowed_transcript_types = [
        "mRNA",
        "transcript",
        "ncRNA",
        "lnc_RNA",
        "miRNA",
        "rRNA",
        "tRNA",
        "snoRNA",
    ]

    def __init__(self, gff_path: str) -> None:
        """
        Creating the GFF database if it does not exist, otherwise loading the existing
        database.

        Args:
            gff_path (str): Path to the GFF file to be parsed.
        Raises:
            FileNotFoundError: If the specified GFF file does not exist.
            ValueError: If the GFF file cannot be parsed.
        """
        self.gff_path = gff_path

        # Validate file exists
        if not os.path.exists(gff_path):
            raise FileNotFoundError(f"GFF3 file not found: {gff_path}")

        if not os.path.isfile(gff_path):
            raise ValueError(f"Path is not a file: {gff_path}")

        try:
            self.db = gffutils.create_db(
                self.gff_path,
                dbfn=":memory:",
                force=True,
                keep_order=True,
                merge_strategy="create_unique",
                sort_attribute_values=True,
            )

        except Exception as e:
            raise ValueError(f"Error parsing GFF3 file {gff_path}: {str(e)}") from e

    @staticmethod
    def _get_id(feature: gffutils.Feature) -> str:
        """
        Method to extract the ID from a feature, handling different possible attribute
        names.

        Args:
            feature (gffutils.Feature): The feature from which to extract the ID.
        Returns:
            str: The extracted ID.
        """
        if "gene_id" in feature.attributes:
            gene_ids = feature.attributes["gene_id"]
            gene_id = gene_ids[0] if isinstance(gene_ids, list) else gene_ids

        elif "Parent" in feature.attributes:
            parents = feature.attributes["Parent"]
            gene_id = parents[0] if isinstance(parents, list) else parents

        else:
            gene_id = feature.id.split(".")[0]

        return gene_id

    def _get_exons(self, transcript_feature: gffutils.Feature) -> list[Exon]:
        """
        Parses the exon features associated with a given transcript feature.

        Args:
            transcript_feature (gffutils.Feature): The transcript feature for which to ext
            ract exons.
        Returns:
            list[Exon]: A list of Exon objects associated with the transcript.
        """
        exons = []
        exon_children = list(self.db.children(transcript_feature, featuretype="exon"))

        if not exon_children:
            logger.warning(
                f"No exon features found for this transcript {transcript_feature.id}"
            )

        # Sort by position, accounting for strand
        exon_children.sort(
            key=lambda x: x.start,
            reverse=(transcript_feature.strand == "-"),
        )

        for exon in exon_children:
            exons.append(
                Exon(
                    seqid=exon.seqid,
                    start=exon.start,
                    end=exon.end,
                    strand=exon.strand,
                    exon_id=exon.id,
                    attributes=dict(exon.attributes),
                )
            )

        return exons

    def _get_cdss(self, transcript_feature: gffutils.Feature) -> list[CDS]:
        """
        Parses the CDS features associated with a given transcript feature.

        Args:
            transcript_feature (gffutils.Feature): The transcript feature for which to
            extract CDS features.
        Returns:
            list[CDS]: A list of CDS objects associated with the transcript.
        """
        cds_features = []
        cds_children = list(self.db.children(transcript_feature, featuretype="CDS"))

        if not cds_children:
            logger.debug(
                f"No CDS features found for this transcript {transcript_feature.id}"
            )

        # Sort by position, accounting for strand
        cds_children.sort(
            key=lambda x: x.start,
            reverse=(transcript_feature.strand == "-"),
        )

        for cds in cds_children:
            # Safely parse phase
            try:
                if cds.frame in (".", None, ""):
                    phase = 0
                else:
                    phase = int(cds.frame)
            except (ValueError, TypeError, AttributeError):
                phase = 0
                logger.warning(f"Invalid phase for CDS {cds.id}, using default 0")

            cds_features.append(
                CDS(
                    seqid=cds.seqid,
                    start=cds.start,
                    end=cds.end,
                    strand=cds.strand,
                    phase=phase,
                    attributes=dict(cds.attributes),
                )
            )
        return cds_features

    def parse_transcript(self, transcript_feature: gffutils.Feature) -> Transcript:
        """
        Parses the GFF file and returns a single Transcript object, format of which was
        defined in models.py.

        Args:
            transcript_feature: gffutils Feature object for a transcript

        Returns:
            Transcript: A single transcript object
        """
        transcript = Transcript(
            transcript_id=transcript_feature.id,
            gene_id=self._get_id(transcript_feature),
            seqid=transcript_feature.seqid,
            start=transcript_feature.start,
            end=transcript_feature.end,
            strand=transcript_feature.strand,
            exons=self._get_exons(transcript_feature),
            cds_features=self._get_cdss(transcript_feature),
            attributes=dict(transcript_feature.attributes),
        )

        return transcript

    def parse_transcripts(self) -> list[Transcript]:
        """
        Creates a list of Transcript objects by parsing the GFF file.

        Returns:
            list[Transcript]: A list of Transcript objects.
        """
        transcripts = []

        for transcript_feature in self.db.features_of_type(
            self.allowed_transcript_types
        ):
            transcript = self.parse_transcript(transcript_feature)
            transcripts.append(transcript)

        return transcripts

    def parse_genes(self) -> Iterator[Gene]:
        """
        Creates an iterator of Gene objects by parsing the GFF file.

        Yields genes one at a time for memory efficiency.

        Yields:
            Gene: Gene objects with their associated transcripts

        Example:
            >>> parser = ParserGFF("annotations.gff3")
            >>> for gene in parser.parse_genes():
            ...     print(f"Gene {gene.gene_id}: {gene.n_transcripts} transcripts")
        """
        for gene_feature in self.db.features_of_type("gene"):
            transcripts = []

            for transcript_feature in self.db.children(
                gene_feature, featuretype=self.allowed_transcript_types
            ):
                transcript = self.parse_transcript(transcript_feature)
                transcripts.append(transcript)

            gene = Gene(
                gene_id=gene_feature.id,
                seqid=gene_feature.seqid,
                start=gene_feature.start,
                end=gene_feature.end,
                strand=gene_feature.strand,
                transcripts=transcripts,
                attributes=dict(gene_feature.attributes),
            )

            yield gene  # FIX: Changed from "yield genes" to "yield gene"


# Remove the print statement for production
# print("ParserGFF module loaded successfully.")
