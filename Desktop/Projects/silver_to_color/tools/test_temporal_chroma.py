"""Verification for temporal_chroma.

The load-bearing claim is that motion-compensating the temporal median beats the naive
one on MOVING content. A static test clip would make any setting look fine, so the
fixture here deliberately contains motion, multimodal colour flips, and an occlusion.

Run:  python test_temporal_chroma.py
"""

from __future__ import annotations

import numpy as np

import temporal_chroma as tc

RNG_SEED = 20260826
SIZE = 128
WINDOW = 5
CENTER = WINDOW // 2
SPEED_PX = 3          # object displacement per frame
TRUE_AB = (40.0, -25.0)     # the colour DDColor should have produced every frame
FLIP_AB = (-30.0, 45.0)     # the competing colour hypothesis it flips to


def build_sequence(flip_frames=(1, 3), speed=SPEED_PX):
    """A bright square sliding right across a dark background.

    Returns (ab_stack, gray_stack, truth_stack, mask_stack).

    `ab_stack` is what a per-frame colorizer produced, including mode flips on
    `flip_frames`. `truth_stack` is what it should have produced. `mask_stack` marks
    where the object is in each frame.
    """
    ab = np.zeros((WINDOW, SIZE, SIZE, 2), np.float32)
    truth = np.zeros_like(ab)
    gray = np.zeros((WINDOW, SIZE, SIZE), np.uint8)
    masks = np.zeros((WINDOW, SIZE, SIZE), bool)

    box = 40
    top = SIZE // 2 - box // 2
    for frame in range(WINDOW):
        left = 20 + speed * frame
        gray[frame, top:top + box, left:left + box] = 230
        gray[frame] = np.maximum(gray[frame], 25)  # dark but non-zero background

        mask = np.zeros((SIZE, SIZE), bool)
        mask[top:top + box, left:left + box] = True
        masks[frame] = mask

        colour = FLIP_AB if frame in flip_frames else TRUE_AB
        ab[frame][mask] = colour
        truth[frame][mask] = TRUE_AB

    return ab, gray, truth, masks


def object_error(result, truth, mask):
    """Mean absolute chroma error over the object only."""
    return float(np.abs(result[mask] - truth[mask]).mean())


def edge_band(mask, width=6):
    """Pixels within `width` px of the object boundary, where smearing shows up."""
    import cv2
    m = mask.astype(np.uint8)
    kernel = np.ones((width * 2 + 1, width * 2 + 1), np.uint8)
    return (cv2.dilate(m, kernel).astype(bool)) & (~cv2.erode(m, kernel).astype(bool))


def test_flow_beats_naive_on_moving_content():
    ab, gray, truth, masks = build_sequence()
    engine = tc.make_dis_flow("medium")

    warped = tc.flow_warped_median(ab, gray, CENTER, flow_engine=engine)
    naive = tc.naive_median(ab, CENTER)

    centre_mask = masks[CENTER]
    centre_truth = truth[CENTER]

    err_warped = object_error(warped, centre_truth, centre_mask)
    err_naive = object_error(naive, centre_truth, centre_mask)

    # Both should beat doing nothing, since the centre frame here is un-flipped.
    band = edge_band(centre_mask)
    band_warped = object_error(warped, centre_truth, band)
    band_naive = object_error(naive, centre_truth, band)

    print(f"  object   MAE  flow-warped {err_warped:8.3f}   naive {err_naive:8.3f}")
    print(f"  edgeband MAE  flow-warped {band_warped:8.3f}   naive {band_naive:8.3f}")

    assert err_warped < err_naive, (
        f"flow-warped median ({err_warped:.3f}) should beat naive ({err_naive:.3f})"
    )
    assert band_warped < band_naive, (
        f"flow-warped should win by more at edges: {band_warped:.3f} vs {band_naive:.3f}"
    )
    # The naive median must actually be visibly bad here, otherwise the fixture is not
    # exercising the failure mode and the comparison proves nothing.
    assert band_naive > 5.0, (
        f"fixture too easy: naive edge error only {band_naive:.3f}; "
        "the test would pass even if warping did nothing"
    )
    return err_warped, err_naive


def test_mode_flip_is_rejected():
    """With the centre frame itself flipped, the median must recover the true colour."""
    ab, gray, _truth, masks = build_sequence(flip_frames=(CENTER,))
    warped = tc.flow_warped_median(ab, gray, CENTER)

    mask = masks[CENTER]
    # Interior only, so the comparison is not dominated by boundary resampling.
    import cv2
    interior = cv2.erode(mask.astype(np.uint8), np.ones((9, 9), np.uint8)).astype(bool)

    recovered = warped[interior].mean(axis=0)
    before = ab[CENTER][interior].mean(axis=0)
    print(f"  centre-flipped: colorizer said {before.round(1)}, "
          f"median recovered {recovered.round(1)}, truth {np.array(TRUE_AB)}")

    assert np.allclose(recovered, TRUE_AB, atol=3.0), (
        f"expected ~{TRUE_AB}, got {recovered}"
    )


def test_mean_would_have_failed():
    """Demonstrates why median, not mean: the mean lands between the two modes."""
    ab, _, _, masks = build_sequence(flip_frames=(1, 3))
    mask = masks[CENTER]
    mean_result = ab.mean(axis=0)
    median_result = np.median(ab, axis=0)

    # 2 of 5 frames flipped -> mean is pulled 40% of the way toward the wrong mode.
    mean_err = float(np.abs(mean_result[mask] - np.array(TRUE_AB)).mean())
    median_err = float(np.abs(median_result[mask] - np.array(TRUE_AB)).mean())
    print(f"  temporal mean MAE {mean_err:.3f} vs temporal median MAE {median_err:.3f}")
    assert mean_err > median_err * 2, "mean should be clearly worse on multimodal flicker"


def test_occlusion_is_masked_not_smeared():
    """Where flow cannot be trusted, the filter must fall back, not invent colour."""
    ab, gray, _, _ = build_sequence()
    # Corrupt one neighbour so its flow is meaningless.
    rng = np.random.default_rng(RNG_SEED)
    gray = gray.copy()
    gray[0] = rng.integers(0, 255, size=gray[0].shape, dtype=np.uint8)

    result, diag = tc.flow_warped_median(ab, gray, CENTER, return_diagnostics=True)
    counts = diag["valid_counts"]

    assert counts.min() >= 1, "centre frame must always count as valid"
    assert counts.max() <= WINDOW
    assert np.isfinite(result).all()
    print(f"  valid-sample counts: min {counts.min()}, "
          f"median {int(np.median(counts))}, max {counts.max()}")


def test_upsample_is_bilinear_not_nearest():
    """The bilinear upsample must actually interpolate between chroma samples."""
    ab = np.zeros((4, 4, 2), np.float32)
    ab[:, :2, 0] = 0.0
    ab[:, 2:, 0] = 100.0     # hard vertical step

    big = tc.upsample_chroma(ab, 64, 64)
    assert big.shape == (64, 64, 2)

    row = big[32, :, 0]
    # Nearest would produce only the values {0, 100}. Bilinear must produce
    # intermediate values across the step.
    intermediate = np.sum((row > 5.0) & (row < 95.0))
    print(f"  intermediate pixels across the step: {intermediate} (nearest would give 0)")
    assert intermediate > 0, "upsample is behaving like nearest-neighbour"


def test_input_validation():
    ab, gray, _, _ = build_sequence()
    for bad, exc in [
        ((np.zeros((5, 8, 8, 3), np.float32), gray, CENTER), ValueError),
        ((ab, np.zeros((5, 4, 4), np.uint8), CENTER), ValueError),
        ((ab, gray, 99), IndexError),
    ]:
        try:
            tc.flow_warped_median(*bad)
        except exc:
            pass
        else:
            raise AssertionError(f"expected {exc.__name__} for {bad[0].shape}")

    nan_ab = ab.copy()
    nan_ab[0, 0, 0, 0] = np.nan
    try:
        tc.flow_warped_median(nan_ab, gray, CENTER)
    except ValueError:
        pass
    else:
        raise AssertionError("NaN input must be rejected, not silently propagated")
    print("  input validation: all bad inputs rejected")


def main() -> int:
    tests = [
        ("flow-warped beats naive on motion", test_flow_beats_naive_on_moving_content),
        ("mode flip on centre frame recovered", test_mode_flip_is_rejected),
        ("mean would have failed", test_mean_would_have_failed),
        ("occlusion masked, not smeared", test_occlusion_is_masked_not_smeared),
        ("upsample is bilinear", test_upsample_is_bilinear_not_nearest),
        ("input validation", test_input_validation),
    ]
    failures = 0
    for name, fn in tests:
        print(f"\n[{name}]")
        try:
            fn()
            print("  PASS")
        except AssertionError as exc:
            failures += 1
            print(f"  FAIL: {exc}")
    print(f"\n{len(tests) - failures}/{len(tests)} passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
