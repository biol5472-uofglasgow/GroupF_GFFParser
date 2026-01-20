import os

import gffutils

from gene_summariser.models import Transcript, Exon, CDS, Feature


class ParserGFF:
    def __init__(self, gff_file: str) -> None:
        self.gff_path = gff_file
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
                    f"Error creating GFF database from {gff_file}: {e}"
                )
            
    def parse_transcripts(self) -> list[Transcript]:

        

        transcripts: list[Transcript] = []
        for feature in self.db.features_of_type("mRNA"):

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