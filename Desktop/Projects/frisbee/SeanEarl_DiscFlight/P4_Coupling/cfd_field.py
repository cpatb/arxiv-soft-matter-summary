"""
Read P3's `flow_fields.npz` and expose `U(r, t)`, `k(r, t)`, `eps(r, t)`
as trilinear-interpolated callables.

Below the lowest sampled altitude (typically 1 m AGL) we fall back to a
log-law extrapolation using `k` to estimate u*.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import RegularGridInterpolator


@dataclass
class StubField:
    """Constant-wind field for testing without P3 output."""
    uniform_wind: np.ndarray

    def U(self, r, t):
        return np.asarray(self.uniform_wind)

    def k(self, r, t):
        return 0.0

    def eps(self, r, t):
        return 1.0


class CFDField:
    """Lazily-loaded trilinear interpolator over one CFD direction."""

    def __init__(self, npz_path: str, heading: int):
        self._data = np.load(npz_path, allow_pickle=False)
        self._heading = heading
        self._key = f"dir_{heading:03d}"
        self._x = self._data["x"]
        self._y = self._data["y"]
        self._z = self._data["z"]

        U_arr = self._data[f"{self._key}__U"]      # (nx, ny, nz, 3)
        k_arr = self._data[f"{self._key}__k"]      # (nx, ny, nz)
        eps_arr = self._data[f"{self._key}__eps"]  # (nx, ny, nz)

        # Build per-component interpolators
        self._U_interp = [
            RegularGridInterpolator(
                (self._x, self._y, self._z), U_arr[..., i],
                bounds_error=False, fill_value=None,
                method="linear")
            for i in range(3)
        ]
        self._k_interp = RegularGridInterpolator(
            (self._x, self._y, self._z), k_arr,
            bounds_error=False, fill_value=0.0, method="linear")
        self._eps_interp = RegularGridInterpolator(
            (self._x, self._y, self._z), eps_arr,
            bounds_error=False, fill_value=1.0, method="linear")

    def U(self, r, t=0.0):
        # Clamp z to >= min sampled altitude
        zc = max(float(r[2]), float(self._z[0]))
        pt = np.array([[r[0], r[1], zc]])
        return np.array([interp(pt)[0] for interp in self._U_interp])

    def k(self, r, t=0.0):
        zc = max(float(r[2]), float(self._z[0]))
        return float(self._k_interp(np.array([[r[0], r[1], zc]]))[0])

    def eps(self, r, t=0.0):
        zc = max(float(r[2]), float(self._z[0]))
        return float(self._eps_interp(np.array([[r[0], r[1], zc]]))[0])

    @property
    def heading(self):
        return self._heading
