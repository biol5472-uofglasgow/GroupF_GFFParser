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
    ):
        self.data = self._process_transcripts(transcripts)
        self.title = title
        self.output_file = output_file

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
