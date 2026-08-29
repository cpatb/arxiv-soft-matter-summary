"""Analysis of catfish fry-ball schooling from handheld overhead video.

No metric calibration exists for these recordings (no ruler, no calibration
target, unknown and drifting camera height). Every length in this package is
therefore either in pixels (internal) or in body lengths, BL (reported). The
body length is measured from the fish themselves, per frame, which has the
useful side effect of absorbing camera-distance drift: if the camera moves
closer, both the school and the fish grow in pixels by the same factor, so the
ratio is preserved.
"""

__all__ = ["frames", "segment", "calibrate", "topology", "flow", "track"]
