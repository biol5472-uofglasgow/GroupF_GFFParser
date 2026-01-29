import sys

from gene_summariser.cli import main

"""
Unit tests for the CLI work slightly different due the argparse usage.
"""


def test_cli_correct_gff(monkeypatch) -> None:
    """
    Test that cli runs with correct gff file
    """
    monkeypatch.setattr(
        sys,
        "argv",
        ["gene_summariser", "--gff", "test/fixtures/models.gff3"],
    )

    try:
        main()
    except SystemExit as e:
        assert e.code == 0  # Expecting normal exit


def test_cli_incorrect_gff(monkeypatch) -> None:
    """
    Test that cli does not run with incorrect gff file
    """
    monkeypatch.setattr(
        sys,
        "argv",
        ["gene_summariser", "--gff", "test/fixtures/incorrect.gff3"],
    )

    try:
        main()
    except SystemExit as e:
        assert e.code == 1  # Expecting error exit
