"""Make ``semiactive`` and ``activation`` importable from any working directory.

Anchored to this file, never to the cwd.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "layer2"))
