import os

import gffutils

from gene_summariser.models import Transcript, Exon, CDS, Feature


class ParserGFF:
    def __init__(self, gff_path: str) -> None:
        self.gff_path = gff_path
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

        transcripts: list[Transcript] = []
        for feature in self.db.features_of_type("mRNA"):

            exons = []
            cds_features = []

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
