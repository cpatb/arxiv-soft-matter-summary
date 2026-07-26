"""
Read P3's `flow_fields.npz` and expose `U(r, t)`, `k(r, t)`, `eps(r, t)`
as trilinearly-interpolated callables.

PERFORMANCE NOTE
----------------
These are called once per right-hand-side evaluation of the trajectory ODE,
i.e. O(10^3) times per throw and O(10^9) times over a full Monte Carlo sweep.
The obvious implementation -- one `scipy.interpolate.RegularGridInterpolator`
per velocity component plus one each for k and eps -- costs five scipy calls
per evaluation, each allocating arrays and running generic N-D machinery. That
was measured to dominate the entire sweep.

We therefore locate the containing cell once with `np.searchsorted` and combine
the eight corners by hand, sharing the cell location and weights across all
five fields. The arithmetic is identical to trilinear
`RegularGridInterpolator`, which is asserted in `test_coupling.py`.

Outside the sampled box, coordinates are CLAMPED rather than extrapolated.
Linear extrapolation of a wind field beyond the CFD domain produces
unphysical, unbounded velocities; clamping degrades gracefully to the nearest
sampled condition. Trajectories that leave the box are flagged via
`fraction_clamped`.
"""

from __future__ import annotations

from dataclasses import dataclass, field as _dc_field

import numpy as np


@dataclass
class StubField:
    """Constant-wind field for testing without P3 output."""
    uniform_wind: np.ndarray

    def U(self, r, t=0.0):
        return np.asarray(self.uniform_wind, dtype=float)

    def k(self, r, t=0.0):
        return 0.0

    def eps(self, r, t=0.0):
        return 1.0


class CFDField:
    """Fast trilinear interpolator over one CFD inflow direction."""

    def __init__(self, npz_path: str, heading: int):
        data = np.load(npz_path, allow_pickle=False)
        self._heading = int(heading)
        self._key = f"dir_{self._heading:03d}"

        if f"{self._key}__U" not in data.files:
            raise KeyError(
                f"{npz_path} has no entry '{self._key}__U'; available "
                f"directions: "
                f"{sorted(n[4:7] for n in data.files if n.endswith('__U'))}"
            )

        self._x = np.ascontiguousarray(data["x"], dtype=float)
        self._y = np.ascontiguousarray(data["y"], dtype=float)
        self._z = np.ascontiguousarray(data["z"], dtype=float)

        self._U = np.ascontiguousarray(data[f"{self._key}__U"], dtype=float)
        self._k = np.ascontiguousarray(data[f"{self._key}__k"], dtype=float)
        self._eps = np.ascontiguousarray(data[f"{self._key}__eps"], dtype=float)

        for name, arr, shape in (
            ("U", self._U, (len(self._x), len(self._y), len(self._z), 3)),
            ("k", self._k, (len(self._x), len(self._y), len(self._z))),
            ("eps", self._eps, (len(self._x), len(self._y), len(self._z))),
        ):
            if arr.shape != shape:
                raise ValueError(
                    f"{self._key}__{name} has shape {arr.shape}, expected "
                    f"{shape} from the x/y/z axes"
                )

        self.is_stub = bool(data["is_stub"]) if "is_stub" in data.files else False

        self._n_eval = 0
        self._n_clamped = 0

    # -- cell location -----------------------------------------------------

    @staticmethod
    def _locate(axis: np.ndarray, v: float):
        """Return (i, w) with axis[i-1] <= v <= axis[i] and w the fraction.

        Values outside the axis are clamped to the end cells (w in [0, 1]),
        never extrapolated. Returns a third element flagging that clamping
        occurred.
        """
        if v <= axis[0]:
            return 1, 0.0, v < axis[0]
        if v >= axis[-1]:
            return len(axis) - 1, 1.0, v > axis[-1]
        i = int(np.searchsorted(axis, v))
        i = min(max(i, 1), len(axis) - 1)
        w = (v - axis[i - 1]) / (axis[i] - axis[i - 1])
        return i, w, False

    def _weights(self, r):
        ix, wx, cx = self._locate(self._x, float(r[0]))
        iy, wy, cy = self._locate(self._y, float(r[1]))
        iz, wz, cz = self._locate(self._z, float(r[2]))
        self._n_eval += 1
        if cx or cy or cz:
            self._n_clamped += 1
        return (ix, wx), (iy, wy), (iz, wz)

    @staticmethod
    def _blend(c000, c001, c010, c011, c100, c101, c110, c111, wx, wy, wz):
        c00 = c000 * (1 - wz) + c001 * wz
        c01 = c010 * (1 - wz) + c011 * wz
        c10 = c100 * (1 - wz) + c101 * wz
        c11 = c110 * (1 - wz) + c111 * wz
        c0 = c00 * (1 - wy) + c01 * wy
        c1 = c10 * (1 - wy) + c11 * wy
        return c0 * (1 - wx) + c1 * wx

    def _interp(self, arr, r):
        (ix, wx), (iy, wy), (iz, wz) = self._weights(r)
        i0, j0, k0 = ix - 1, iy - 1, iz - 1
        return self._blend(
            arr[i0, j0, k0], arr[i0, j0, iz], arr[i0, iy, k0], arr[i0, iy, iz],
            arr[ix, j0, k0], arr[ix, j0, iz], arr[ix, iy, k0], arr[ix, iy, iz],
            wx, wy, wz,
        )

    # -- public interface --------------------------------------------------

    def U(self, r, t=0.0):
        return self._interp(self._U, r)

    def U_and_k(self, r, t=0.0):
        """Mean velocity and TKE from a SINGLE cell location.

        The composed wind callable needs both at the same point on every
        right-hand-side evaluation (U for the mean, k to scale the gust). Doing
        them separately repeats the searchsorted and the eight-corner weights;
        sharing them halves the interpolation cost of the whole sweep.
        """
        (ix, wx), (iy, wy), (iz, wz) = self._weights(r)
        i0, j0, k0 = ix - 1, iy - 1, iz - 1
        U = self._blend(
            self._U[i0, j0, k0], self._U[i0, j0, iz],
            self._U[i0, iy, k0], self._U[i0, iy, iz],
            self._U[ix, j0, k0], self._U[ix, j0, iz],
            self._U[ix, iy, k0], self._U[ix, iy, iz],
            wx, wy, wz)
        kk = self._blend(
            self._k[i0, j0, k0], self._k[i0, j0, iz],
            self._k[i0, iy, k0], self._k[i0, iy, iz],
            self._k[ix, j0, k0], self._k[ix, j0, iz],
            self._k[ix, iy, k0], self._k[ix, iy, iz],
            wx, wy, wz)
        return U, float(kk)

    def k(self, r, t=0.0):
        return float(self._interp(self._k, r))

    def eps(self, r, t=0.0):
        return float(self._interp(self._eps, r))

    @property
    def heading(self):
        return self._heading

    @property
    def fraction_clamped(self) -> float:
        """Fraction of queries that fell outside the sampled box."""
        return self._n_clamped / self._n_eval if self._n_eval else 0.0

    @property
    def bounds(self):
        return ((self._x[0], self._x[-1]),
                (self._y[0], self._y[-1]),
                (self._z[0], self._z[-1]))
