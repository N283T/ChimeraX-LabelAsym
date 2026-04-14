"""Make the bundle's ``src`` directory importable as ``hook``."""

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent.parent.joinpath("src")
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
