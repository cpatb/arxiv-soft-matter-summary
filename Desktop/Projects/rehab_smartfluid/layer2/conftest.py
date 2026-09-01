"""Make ``activation`` importable no matter which directory pytest was started from.

Anchored to this file rather than the cwd: a cwd-relative path literal makes
"the tests pass" true only from one directory, which has bitten this repo's
owner before (seven frisbee P5 tests were green only from P5_Synthesis).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
