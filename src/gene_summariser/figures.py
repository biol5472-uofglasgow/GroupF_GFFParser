from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # MUST come before pyplot

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
        self.output_dir = Path(output_dir) / "figures"

    def _process_transcripts(
        self, transcripts: list[TranscriptSummary]
    ) -> dict[str, int]:
        flag_counts: dict[str, int] = {"no_flag": 0}

        for transcript in transcripts:
            if not transcript.flags:
                flag_counts["no_flag"] += 1
                continue
            for flag in transcript.flags:
                if flag in flag_counts:
                    flag_counts[flag] += 1
                else:
                    flag_counts[flag] = 1

        return flag_counts

    def generate_pie_chart(self) -> None:
        labels = list(self.data.keys())
        values = list(self.data.values())

        plt.figure(figsize=(8, 8))
        plt.pie(values, startangle=90)
        plt.legend(
            labels,
            loc="center left",
            bbox_to_anchor=(1, 0.5),
        )
        plt.title(self.title)
        plt.axis("equal")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(
            self.output_dir / self.output_file,
            bbox_inches="tight",
        )
        plt.close()


class FlaggedBarChart:
    """
    Generate a bar chart showing the number of flagged and unflagged transcripts.
    """

    def __init__(
        self,
        transcripts: list[TranscriptSummary],
        title="Flagged vs Unflagged Transcripts",
        output_file="flagged_vs_unflagged_bar.png",
        output_dir: Path = Path("."),
    ):
        self.counts = self._count_flags(transcripts)
        self.title = title
        self.output_file = output_file
        self.output_dir = Path(output_dir) / "figures"

    def _count_flags(self, transcripts) -> dict[str, int]:
        counts = {"flagged": 0, "unflagged": 0}
        for transcript in transcripts:
            if transcript.flags:
                counts["flagged"] += 1
            else:
                counts["unflagged"] += 1
        return counts

    def generate_bar_plot(self) -> None:
        group = ["flagged", "unflagged"]
        counts = [self.counts["flagged"], self.counts["unflagged"]]
        labels = ["Flagged", "Unflagged"]
        bar_colours = ["#409ab6", "#ce6868"]

        plt.figure(figsize=(4, 6))
        plt.bar(group, counts, label=labels, color=bar_colours)

        plt.xlabel("Transcript Status")
        plt.ylabel("Number of transcripts")
        plt.title(self.title)
        plt.legend()

        plt.tight_layout()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(self.output_dir / self.output_file)
        plt.close()


class ExonCountHistogram:
    """
    Generate a histogram showing the distribution of exon counts across transcripts.
    """

    def __init__(
        self,
        transcripts: list[TranscriptSummary],
        title="Exon Count Distribution",
        output_file="exon_count_distribution.png",
        output_dir: Path = Path("."),
    ):
        self.exon_counts = [transcript.n_exons for transcript in transcripts]
        self.title = title
        self.output_file = output_file
        self.output_dir = Path(output_dir) / "figures"

    def generate_histogram(self) -> None:
        plt.figure(figsize=(8, 6))
        plt.hist(
            self.exon_counts,
            bins=range(1, max(self.exon_counts) + 2),
            edgecolor="black",
            alpha=0.7,
        )
        plt.xlabel("Number of Exons")
        plt.ylabel("Number of Transcripts")
        plt.title(self.title)
        # Forcing x-ticks to be every integer exon count, instead of matplotlibs default
        plt.xticks(range(1, max(self.exon_counts) + 1))
        # Adding a grid and setting its transparency
        plt.grid(axis="y", alpha=0.75)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(self.output_dir / self.output_file)
        plt.close()


class CDSLengthHistogram:
    """
    Generate a histogram showing the distribution of CDS lengths across transcripts.
    """

    def __init__(
        self,
        transcripts: list[TranscriptSummary],
        title="CDS Length Distribution",
        output_file="cds_length_distribution.png",
        output_dir: Path = Path("."),
    ):
        self.cds_lengths = [
            transcript.total_cds_length
            for transcript in transcripts
            if transcript.total_cds_length is not None
            and transcript.total_cds_length > 0
        ]
        self.title = title
        self.output_file = output_file
        self.output_dir = Path(output_dir) / "figures"

    def generate_histogram(self) -> None:
        plt.figure(figsize=(8, 6))
        plt.hist(
            self.cds_lengths,
            edgecolor="black",
            alpha=0.7,
        )
        plt.xlabel("CDS Length")
        plt.ylabel("Number of Transcripts")
        plt.title(self.title)
        # Adding a grid and setting its transparency
        plt.grid(axis="y", alpha=0.75)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(self.output_dir / self.output_file)
        plt.close()
