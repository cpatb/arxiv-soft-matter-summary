"""Tests for studio.colorize that need no model: estimators, decimation, chunk alignment.

The chunk-alignment test is the important one. Chunks are read by seeking, and a seek
that lands one frame off would shift colour against luminance by a frame without any
error being raised. A fake colorizer whose chroma encodes the frame's own content makes
any misalignment visible as a pixel difference against a single-chunk run.
"""

import subprocess

import cv2
import numpy as np
import pytest

from studio import colorize as C
import temporal_chroma as tc


def test_masked_median_matches_nanmedian():
    rng = np.random.default_rng(0)
    stack = rng.normal(size=(5, 16, 16, 2)).astype(np.float32)
    mask = rng.random((5, 16, 16)) > 0.4
    mask[0] = True                                  # centre always valid
    ref = np.nanmedian(np.where(mask[..., None], stack, np.nan), axis=0)
    np.testing.assert_allclose(C.masked_median(stack, mask), ref, rtol=0, atol=1e-6)


def test_masked_median_all_valid_is_plain_median():
    rng = np.random.default_rng(1)
    stack = rng.normal(size=(4, 8, 8, 2)).astype(np.float32)
    np.testing.assert_allclose(C.masked_median(stack, np.ones((4, 8, 8), bool)),
                               np.median(stack, axis=0), atol=1e-6)


def test_duplicate_mask_finds_one_in_three():
    # 16->24 recipe: every third diff is a duplicate (small), the others real motion.
    rng = np.random.default_rng(2)
    diffs = rng.uniform(5, 11, size=90)
    diffs[::3] = rng.uniform(0.8, 1.5, size=30)
    mask = C.duplicate_mask(diffs)
    assert mask[::3].all()
    assert not mask[1::3].any() and not mask[2::3].any()


def test_duplicate_mask_keeps_everything_when_nothing_repeats():
    # Negative case: a genuinely 24 fps film must lose no frames.
    rng = np.random.default_rng(3)
    assert not C.duplicate_mask(rng.uniform(5, 11, size=200)).any()


def test_temporal_median_tracks_motion_like_reference():
    """On a translating textured scene, our 256-grid flow estimator should agree with the
    reference 512-grid tc.flow_warped_median much better than a naive median does."""
    rng = np.random.default_rng(4)
    base = cv2.GaussianBlur((rng.random((560, 560)) * 255).astype(np.uint8), (0, 0), 3)
    base = cv2.normalize(base, None, 0, 255, cv2.NORM_MINMAX)
    grays, abs_ = [], []
    for k in range(5):
        g = base[8 * k:8 * k + 512, 8 * k:8 * k + 512].copy()      # 8 px/frame drift
        grays.append(g)
        a = (g.astype(np.float32) - 128) / 4
        ab = np.stack([a, -a], axis=-1)
        ab += rng.normal(scale=6, size=ab.shape).astype(np.float32)  # per-frame chroma noise
        abs_.append(ab.astype(np.float32))
    flows = C.pair_flows(grays, 2, threads=2)
    ours = C.temporal_median(abs_, flows, 2, 2)
    ref = tc.flow_warped_median(np.stack(abs_), np.stack(grays), 2)
    naive = np.median(np.stack(abs_), axis=0)
    inner = (slice(40, -40), slice(40, -40))
    err_ours = np.abs(ours - ref)[inner].mean()
    err_naive = np.abs(naive - ref)[inner].mean()
    assert err_ours < 0.5 * err_naive, (err_ours, err_naive)


class FakeColorizer:
    """Chroma is a deterministic function of the frame, so misalignment shows."""
    backend = "fake"

    def predict_ab(self, frames):
        out = []
        for f in frames:
            s = cv2.resize(f, (C.INPUT_SIZE, C.INPUT_SIZE)).astype(np.float32)
            out.append(np.stack([(s - 128) / 3, (128 - s) / 5], axis=-1))
        return out


def _make_clip(path, n=37, fps=16):
    """A grey clip whose frames are all very different (frame number as big digits)."""
    proc = C.writer(path, 320, 240, fps, "gray", crf=0)
    for i in range(n):
        img = np.full((240, 320), 40, np.uint8)
        cv2.putText(img, str(i), (20, 200), cv2.FONT_HERSHEY_SIMPLEX, 6, 230, 12)
        proc.stdin.write(img.tobytes())
    proc.stdin.close()
    assert proc.wait() == 0


@pytest.mark.skipif(subprocess.run(["ffmpeg", "-version"], capture_output=True).returncode,
                    reason="ffmpeg not available")
def test_chunked_colorize_equals_single_chunk(tmp_path):
    src = tmp_path / "bw.mp4"
    _make_clip(src)
    one = C.colorize(src, tmp_path / "one.mp4", FakeColorizer(), chunk_frames=1000,
                     workdir=tmp_path / "w1", crf=0)
    many = C.colorize(src, tmp_path / "many.mp4", FakeColorizer(), chunk_frames=7,
                      workdir=tmp_path / "w2", crf=0)
    assert one["source_frames"] == many["source_frames"] == 37
    a = list(C.read_gray(tmp_path / "one.mp4", 320, 240))
    b = list(C.read_gray(tmp_path / "many.mp4", 320, 240))
    assert len(a) == len(b) == 37
    worst = max(np.abs(x.astype(int) - y.astype(int)).max() for x, y in zip(a, b))
    assert worst <= 2, worst          # lossless encode; 2 levels allows yuv420p rounding
    # Discriminating negative: a one-frame shift must fail the same bar.
    shifted = max(np.abs(x.astype(int) - y.astype(int)).max() for x, y in zip(a[1:], b))
    assert shifted > 50, shifted
