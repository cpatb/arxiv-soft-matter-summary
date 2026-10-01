"""Render alignment: real ffmpeg on a synthetic numbered clip.

Trick: feed the SAME clip as both the B&W and the colour master. A wipe that keeps both
sides on the same film frame is then pixel-identical to a plain single-stream segment;
any cadence misalignment shows up as a whole-frame difference.
"""

import subprocess
from pathlib import Path

import cv2
import numpy as np
import pytest

from studio import colorize as C
from studio import render as R
from studio import timeline as T
from studio.episode import Cue, Episode

pytestmark = pytest.mark.skipif(
    subprocess.run(["ffmpeg", "-version"], capture_output=True).returncode,
    reason="ffmpeg not available")


def _clip(path: Path, n=320, fps=16):
    proc = C.writer(path, 320, 240, fps, "gray", crf=0)
    for i in range(n):
        img = np.full((240, 320), 30, np.uint8)
        cv2.putText(img, str(i), (10, 190), cv2.FONT_HERSHEY_SIMPLEX, 5, 230, 14)
        proc.stdin.write(img.tobytes())
    proc.stdin.close()
    assert proc.wait() == 0


def _frames(path: Path):
    return list(C.read_gray(path, 320, 180))


@pytest.mark.parametrize("film_at,speed,dur", [(10.0, 1.0, 6.0), (10.03, 1.0, 6.0),
                                               (7.0, 0.5, 6.0), (10.0, 1.0, 3.0)])
def test_wipe_is_frame_aligned(tmp_path, film_at, speed, dur):
    _clip(tmp_path / "m.mp4")
    ep = Episode(path=tmp_path / "e.yaml", id="t", title="t", bw=tmp_path / "m.mp4",
                 color=tmp_path / "m.mp4", credit="", output=tmp_path, voice={},
                 recordings=tmp_path, facts_path=None, cues=[], shorts=[])
    out = {}
    for mode in ("wipe", "color"):
        cue = Cue(id=mode, at=film_at, mode=mode, speed=speed, hold=dur)
        placed = T.layout([cue], {})[0]
        cmd = R.segment_cmd(ep, placed, R.FULL, f"{mode}.mp4")
        R.run(cmd, cwd=tmp_path)
        out[mode] = _frames(tmp_path / f"{mode}.mp4")
    a, b = out["wipe"], out["color"]
    assert len(a) == len(b) == round(dur * R.FPS)
    bad = sum(np.abs(x.astype(int) - y.astype(int)).max() > 40 for x, y in zip(a, b))
    assert bad == 0, f"{bad} of {len(a)} frames differ"
    # Discriminating negative: shifting by one output frame must be caught by the same bar
    shifted = sum(np.abs(x.astype(int) - y.astype(int)).max() > 40 for x, y in zip(a[1:], b))
    assert shifted > len(a) // 4
