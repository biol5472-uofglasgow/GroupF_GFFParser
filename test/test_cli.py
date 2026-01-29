import sys

from gene_summariser.cli import main

"""
Unit tests for the CLI work slightly different due the argparse usage.
As alot of the actual output from the CLI is tested in the integration test
this mainly tests if the CLI runs or exits correctly without worrying about output.
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
        assert e.code == 0


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
        assert e.code == 1


def test_cli_help(monkeypatch) -> None:
    """
    Test that cli help works
    """
    monkeypatch.setattr(
        sys,
        "argv",
        ["gene_summariser", "--help"],
    )

    try:
        main()
    except SystemExit as e:
        assert e.code == 0


def test_cli_strict_mode(monkeypatch) -> None:
    """
    Test that cli fails in strict mode when flags are present
    2 transcripts in the test gff have qc flags so should exit with code 2
    """
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gene_summariser",
            "--gff",
            "test/fixtures/models.gff3",
            "--strict",
        ],
    )

    try:
        main()
    except SystemExit as e:
        assert e.code == 2
