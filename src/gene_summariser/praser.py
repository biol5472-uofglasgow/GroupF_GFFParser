import os

import gffutils


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
