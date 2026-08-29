"""Frame access and per-clip metadata.

Frames were extracted with ffmpeg, which bakes the container's display rotation
into the pixels. Working from extracted frames rather than decoding the .MOV
directly avoids the well-known ambiguity in how different OpenCV versions handle
QuickTime rotation metadata -- a silent orientation bug that would corrupt every
direction-dependent result downstream.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from skimage import color, io

ROOT = Path(__file__).resolve().parents[2]
FRAME_ROOT = ROOT / "data" / "frames"

CLIPS = ("IMG_0544", "IMG_0545", "IMG_0546", "IMG_0547")

# All four clips share this rate (30000/1001 from the container).
FPS = 30000.0 / 1001.0


@dataclass(frozen=True)
class Clip:
    name: str
    paths: tuple[Path, ...]
    fps: float = FPS

    def __len__(self) -> int:
        return len(self.paths)

    @property
    def duration_s(self) -> float:
        return len(self) / self.fps

    def time_s(self, index: int) -> float:
        return index / self.fps

    def rgb(self, index: int) -> np.ndarray:
        """Float RGB in [0, 1], shape (H, W, 3)."""
        img = io.imread(self.paths[index])
        if img.ndim != 3 or img.shape[2] != 3:
            raise ValueError(f"{self.paths[index]}: expected RGB, got shape {img.shape}")
        if img.dtype != np.uint8:
            raise ValueError(f"{self.paths[index]}: expected uint8, got {img.dtype}")
        return img.astype(np.float64) / 255.0

    def gray_u8(self, index: int) -> np.ndarray:
        """Luminance as uint8 -- the form OpenCV routines expect."""
        img = io.imread(self.paths[index])
        return (color.rgb2gray(img.astype(np.float64) / 255.0) * 255.0).astype(np.uint8)


@functools.lru_cache(maxsize=8)
def load_clip(name: str) -> Clip:
    d = FRAME_ROOT / name
    if not d.is_dir():
        raise FileNotFoundError(f"no extracted frames at {d}")
    paths = tuple(sorted(d.glob("*.jpg")))
    if not paths:
        raise FileNotFoundError(f"no frames in {d}")
    return Clip(name=name, paths=paths)


def all_clips() -> list[Clip]:
    return [load_clip(c) for c in CLIPS]
