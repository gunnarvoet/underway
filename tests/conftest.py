from pathlib import Path

import pytest


@pytest.fixture
def data():
    """Directory with truncated real raw files."""
    return Path(__file__).parent / "data"
