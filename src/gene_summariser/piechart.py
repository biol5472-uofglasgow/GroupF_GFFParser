from pathlib import Path

import matplotlib.pyplot as plt

from gene_summariser.models import TranscriptSummary


class PieChart:
    """
    Generate a pie chart showing the distribution of flags in transcript summaries.
    These flags are stored within TranscriptSummary.flags as a list of strings.
    """

    def __init__(
        self,
        transcripts: list[TranscriptSummary],
        title="Flag Distribution",
        output_file="pie_chart.png",
        output_dir: Path = Path("."),
    ):
        self.data = self._process_transcripts(transcripts)
        self.title = title
        self.output_file = output_file
        self.output_dir = output_dir

    def _process_transcripts(
        self, transcripts: list[TranscriptSummary]
    ) -> dict[str, int]:
        flag_counts: dict[str, int] = {}

        for transcript in transcripts:
            for flag in transcript.flags:
                if flag in flag_counts:
                    flag_counts[flag] += 1
                else:
                    flag_counts[flag] = 1

        return flag_counts

    def generate_pie_chart(self):
        labels = list(self.data.keys())
        values = list(self.data.values())

        plt.figure(figsize=(8, 8))
        plt.pie(values, labels=labels)
        plt.title(self.title)
        plt.axis("equal")

        plt.savefig(self.output_dir / self.output_file)
        plt.close()
