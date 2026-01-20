import os

import gffutils

from gene_summariser.models import Transcript, Exon, CDS, Feature


class ParserGFF:
    """
    PraserGFF provides the nessesary methods to parse a specified GFF file, extracting Gene, transcript, exon, and CDS information.

    Attributes:
        gff_path (str): Path to the GFF file to be parsed.
        db (gffutils.FeatureDB): GFF database object for efficient querying of features.

    Methods:
        parse_transcripts() -> list[Transcript]: Parses the GFF file and returns a list of Transcript objects, each containing associated exons and CDS features.
        parse_genes() -> list[Gene]: Parses the GFF file and returns a list of Gene features. (To be implemented)
    """

    def __init__(self, gff_path: str) -> None:
        """
        Creating the GFF database if it does not exist, otherwise loading the existing database.

        Args:
            gff_path (str): Path to the GFF file to be parsed.
        Raises:
            FileNotFoundError: If the specified GFF file does not exist.
        """
        self.gff_path = gff_path
        # Creating the database path from the GFF file path, this is used to both find or create the database.
        db_path = os.path.splitext(self.gff_path)[0] + ".db"

        if os.path.exists(db_path):
            self.db = gffutils.FeatureDB(db_path, keep_order=True)
        else:
            try:
                self.db = gffutils.create_db(
                    self.gff_path,
                    dbfn=db_path,
                    force=True,
                    keep_order=True,
                    sort_attribute_values=True,
                )
            except FileNotFoundError as e:
                raise FileNotFoundError(
                    f"Error creating GFF database from {self.gff_path}: {e}"
                )

    def parse_transcripts(self) -> list[Transcript]:
        """
        Parses the GFF file and returns a list of Transcript objects, format of which was defined in models.py.
        Returns:
            list[Transcript]: A list of Transcript objects with associated exons and CDS features.
        Raises:
            (To be implemented)
        """
        transcripts: list[Transcript] = []
        # In GFF3 files, mRNA features represent transcripts.
        for feature in self.db.features_of_type("mRNA"):

            exons = []
            cds_features = []

            # Creating the list of Exon and CDS objects associated to the current transcript.
            for exon in self.db.children(feature, featuretype="exon"):
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

            for cds in self.db.children(feature, featuretype="CDS"):
                cds_features.append(
                    CDS(
                        seqid=cds.seqid,
                        start=cds.start,
                        end=cds.end,
                        strand=cds.strand,
                        phase=int(cds.frame) if cds.frame is not None else 0,
                        attributes=dict(cds.attributes),
                    )
                )

            transcript = Transcript(
                transcript_id=feature.id,
                gene_id=feature.attributes.get("gene_id", [""])[0],
                seqid=feature.seqid,
                start=feature.start,
                end=feature.end,
                strand=feature.strand,
                exons=exons,
                cds_features=cds_features,
                attributes=dict(feature.attributes),
            )

            transcripts.append(transcript)

        return transcripts
