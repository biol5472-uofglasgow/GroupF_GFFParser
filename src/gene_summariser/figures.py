from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # MUST come before pyplot

import matplotlib.pyplot as plt

from gene_summariser.models import TranscriptSummary


class PieChart:
    """
    Generate a pie chart showing the distribution of flags in transcript summaries.
    This will allow visualisation of the most common flags / issues across transcripts.
    """

    def __init__(
        self,
        transcripts: list[TranscriptSummary],
        title: str="Flag Distribution",
        output_file: str="pie_chart.png",
        output_dir: Path = Path("."),
    ):
        # Prossess the flag counts in a Dict using the helper method
        self.data = self._process_transcripts(transcripts)
        self.title = title
        self.output_file = output_file
        # Set output directory for figures
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

        # Setting a consistent size for all charts
        plt.figure(figsize=(8, 8))
        plt.pie(values, startangle=90)
        # Ensures the legend is not overlapping with the piechart
        plt.legend(
            labels,

            bbox_to_anchor=(1, 0.5),
        )
        plt.title(self.title)
        # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.axis("equal")
        # Creating output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        # Saving the figure, setting bbox_inches to tight to ensure nothing is cut off
        plt.savefig(
            self.output_dir / self.output_file,
            bbox_inches="tight",
        )
        plt.close()


class FlaggedBarChart:
    """
    Generate a bar chart showing the number of flagged and unflagged transcripts.
    Gives a more general overview of the QC status of the transcripts.
    Allowing identification of "good" vs "problematic" transcripts.
    """

    def __init__(
        self,
        transcripts: list[TranscriptSummary],
        title:str="Flagged vs Unflagged Transcripts",
        output_file:str="flagged_vs_unflagged_bar.png",
        output_dir: Path = Path("."),
    ):
        # Count flagged vs unflagged transcripts using helper method
        self.counts = self._count_flags(transcripts)
        self.title = title
        self.output_file = output_file
        self.output_dir = Path(output_dir) / "figures"

    def _count_flags(
            self, transcripts: list[TranscriptSummary]) -> dict[str, int]:
        counts = {"flagged": 0, "unflagged": 0}
        for transcript in transcripts:
            if transcript.flags:
                counts["flagged"] += 1
            else:
                counts["unflagged"] += 1
        return counts

    # Setting up the lists used for plotting the bar plot
    def generate_bar_plot(self) -> None:
        group = ["flagged", "unflagged"]
        counts = [self.counts["flagged"], self.counts["unflagged"]]
        labels = ["Flagged", "Unflagged"]
        bar_colours = ["#409ab6", "#ce6868"]

        # Setting a consistent size for all charts
        plt.figure(figsize=(8, 8))
        plt.bar(group, counts, label=labels, color=bar_colours)

        plt.xlabel("Transcript Status")
        plt.ylabel("Number of transcripts")
        plt.title(self.title)
        plt.legend()

        plt.tight_layout()
        # Creating output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(self.output_dir / self.output_file)
        plt.close()


class ExonCountHistogram:
    """
    Generate a histogram showing the distribution of exon counts across transcripts.
    Allows for visualisation of how exons are distibuted within the dataset
    """

    def __init__(
        self,
        transcripts: list[TranscriptSummary],
        title:str="Exon Count Distribution",
        output_file:str="exon_count_distribution.png",
        output_dir: Path = Path("."),
    ):
        # Extract exon counts from transcripts
        self.exon_counts = [transcript.n_exons for transcript in transcripts]
        self.title = title
        self.output_file = output_file
        self.output_dir = Path(output_dir) / "figures"

    def generate_histogram(self) -> None:
        # Setting a consistent size for all charts
        plt.figure(figsize=(8, 8))
        # Bins set to give one bin per exon count value, eg 1,2,3,4...
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
        # Creating output directory if it doesn't exist
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
        title:str="CDS Length Distribution",
        output_file:str="cds_length_distribution.png",
        output_dir: Path = Path("."),
    ):
        # Extract CDS lengths from transcripts, ensures the len is not None or 0
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
        # Setting a consistent size for all charts
        plt.figure(figsize=(8, 8))
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
        # Creating output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(self.output_dir / self.output_file)
        plt.close()
