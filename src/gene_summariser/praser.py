import os

import gffutils

import logging

from gene_summariser.models import Gene, Transcript, Exon, CDS, Feature

logger = logging.getLogger(__name__)


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

    allowed_transcript_types = ["mRNA", "transcript", "ncRNA", "lnc_RNA", 
                       "miRNA", "rRNA", "tRNA", "snoRNA"]
    

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
        if not os.path.exists(gff_path):
            raise FileNotFoundError(f"GFF3 file not found: {gff_path}")

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
            raise ValueError(
                f"Error parsing GFF3 file {gff_path}: {str(e)}"
            ) from e
    
    def _get_id(self, feature) -> str:
        """
        Method to extract the ID from a feature, handling different possible attribute names.
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
        # Sorting them based on their start position, taking into account the strand of the transcript. This wille make gene flags using FASTA data easier later on.
        exon_children = list(self.db.children(transcript_feature, featuretype="exon"))
        
        # This can be changed when doing flags, however for now its usful for testing that its caught, same goes for no CDS
        if not exon_children:
            logger.warning(f"No exon features found for this transcript {transcript_feature.id}")
        
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


        cds_children = list(self.db.children(transcript_feature, featuretype="CDS"))
        
        if not cds_children:
            logger.warning(f"No CDS features found for this transcript {transcript_feature.id}")

        cds_children.sort(
            key=lambda x: x.start,
            reverse=(transcript_feature.strand == "-"),
        )

        for cds in cds_children:
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
            gene_id=self._get_id(transcript_feature),
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

        for transcript_feature in self.db.features_of_type(self.allowed_transcript_types):
            transcript = self.parse_transcript(transcript_feature)
            transcripts.append(transcript)

        return transcripts
    
    def parse_genes(self) -> list[Gene]:
        """
        Creates a list of Gene objects by parsing the GFF file.
        Returns:
            list[Gene]: A list of Gene objects.
        Raises:
            (To be implemented)
        """
        genes = []
        
        for gene_feature in self.db.features_of_type("gene"):
            transcripts = []

            for transcript_feature in self.db.children(gene_feature, featuretype=self.allowed_transcript_types):
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
            

