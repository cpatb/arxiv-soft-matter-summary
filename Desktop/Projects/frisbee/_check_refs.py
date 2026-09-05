import re
import sys
from pathlib import Path

src = Path(sys.argv[1]).read_text(encoding="utf-8")
labels = set(re.findall(r"\\label\{([^}]+)\}", src))
refs = set(re.findall(r"\\(?:ref|eqref)\{([^}]+)\}", src))
print("dangling:", sorted(refs - labels))
for m in re.finditer(r"Qbyheading", src):
    ln = src[: m.start()].count("\n") + 1
    seg = src[max(0, m.start() - 100): m.start() + 40].replace("\n", " | ")
    print("  line %d: ...%s" % (ln, seg))
