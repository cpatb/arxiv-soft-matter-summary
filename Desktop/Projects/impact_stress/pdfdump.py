import sys, fitz
p = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else None
d = fitz.open(p)
txt = []
for i, pg in enumerate(d):
    txt.append(f"\n===== PAGE {i+1} =====\n" + pg.get_text())
s = "".join(txt)
if out:
    open(out, "w", encoding="utf-8").write(s)
    print(f"{p}: {len(d)} pages -> {out} ({len(s)} chars)")
else:
    print(s)
