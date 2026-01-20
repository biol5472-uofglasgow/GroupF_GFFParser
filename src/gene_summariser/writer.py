"""Output writers for various file formats."""

import json
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd

from gene_summariser.models import TranscriptSummary, Transcript

class OutputWriter:
    """Handles writing outputs to various formats."""