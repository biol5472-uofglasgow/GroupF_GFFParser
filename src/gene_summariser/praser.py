import os

import gffutils


class ParserGFF:
    def __init__(self, gff_file: str) -> None:
        self.gff_path = os.path.join(os.getcwd(), gff_file)

        if os.path.exists(os.path.splitext(self.gff_path)[0] + ".db"):
            self.db = gffutils.FeatureDB(
                os.path.splitext(self.gff_path)[0] + ".db", keep_order=True
            )
        else:
            try:
                self.db = gffutils.create_db(
                    self.gff_path,
                    dbfn=":memory:",
                    force=True,
                    keep_order=True,
                    sort_attribute_values=True,
                )
            except FileNotFoundError as e:
                raise FileNotFoundError(
                    f"Error creating GFF database from {gff_file}: {e}"
                )
