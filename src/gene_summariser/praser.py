import os

import gffutils

from gene_summariser.models import Gene, Transcript, Exon, CDS, Feature


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

    def parse_transcript(self, transcript_feature) -> Transcript:
        """
        Parses the GFF file and returns a single Transcript objects, format of which was defined in models.py.
        Returns:
            Transcript: A single transcript object
        Raises:
            (To be implemented)
        """

        exons = []
        cds_features = []

        # Creating the list of Exon and CDS objects associated to the current transcript.
        for exon in self.db.children(transcript_feature, featuretype="exon"):
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

        for cds in self.db.children(transcript_feature, featuretype="CDS"):
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
            transcript_id=transcript_feature.id,
            gene_id=transcript_feature.attributes.get("gene_id", [""])[0],
            seqid=transcript_feature.seqid,
            start=transcript_feature.start,
            end=transcript_feature.end,
            strand=transcript_feature.strand,
            exons=exons,
            cds_features=cds_features,
            attributes=dict(transcript_feature.attributes),
        )

        return transcript

    def parse_transcripts(self) -> list[Transcript]:
        """
        Creates a list of Transcript objects by parsing the GFF file.
        Returns:
            list[Transcript]: A list of Transcript objects.
        Raises:
            (To be implemented)
        """
        transcripts = []

        for transcript_feature in self.db.features_of_type("mRNA"):
            transcript = self.parse_transcript(transcript_feature)
            transcripts.append(transcript)

        return transcripts
    
    def parse_genes(self) -> list[Gene]:
        
        genes = []
        
        for gene_feature in self.db.features_of_type("gene"):
            transcripts = []

            for transcript_feature in self.db.children(gene_feature, featuretype="mRNA"):
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
            genes.append(gene)

        return genes
            

