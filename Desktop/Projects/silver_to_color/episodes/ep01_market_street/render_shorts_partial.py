"""Render ep01's Shorts from the partial colour master (film 0-470 s).

All four Shorts take footage from before 300 s, so the frames are identical to what the
finished master holds there. Delete this script once color_master.mp4 exists and use
`python -m studio shorts episodes/ep01_market_street/episode.yaml`.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from studio import episode as E, render as R  # noqa: E402

PARTIAL = Path("E:/silver_to_color/ep01/color_partial_0-470s.mp4")
ep = E.load(Path(__file__).with_name("episode.yaml"))
ep.color = PARTIAL
for s in ep.shorts:
    for c in s.cues:
        assert c.at is None or c.at < 440, (s.id, c.id, c.at)
    final, length = R.render_short(ep, s)
    print(f"{s.id}: {final} ({length:.1f}s)", flush=True)
