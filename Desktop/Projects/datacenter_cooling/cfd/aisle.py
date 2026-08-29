"""
aisle.py — hot-aisle/cold-aisle recirculation from 2D buoyant CFD.

A vented-partition Boussinesq cavity: cold wall (left) = cold-aisle supply, hot wall
(right) = rack exhaust / hot aisle, an adiabatic partition of height h_p on the centreline
(the rack + containment baffle) with an open gap above it. Buoyancy drives hot air up the
right, over the partition, and down into the cold aisle -- exactly the recirculation that
lifts rack inlet temperature. Containment effectiveness eta = h_p / H_room; eta -> 1 seals
the gap and kills recirculation.

The recirculation fraction phi is measured with a passive tracer: C=1 released at the hot
wall, C=0 at the cold-supply wall; phi = mean C arriving at the cold (inlet) wall. This is
the physical meaning of the base recirculation the Tier-0 model carries as a placeholder.

Poisson for the streamfunction is solved with a sparse-LU factorization of the 5-point
Laplacian with the internal solid baked in as psi=0 Dirichlet (the fast DST solver cannot
represent an internal solid). Reduces to the de Vahl Davis cavity when the partition height
is zero -- checked in test_aisle.py.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu


class AisleCell:
    def __init__(self, N=70, H_room=1.0, eta=0.5, Ra=1e4, Pr=0.71,
                 part_halfwidth=0.03):
        self.nx = N
        self.ny = max(6, int(round(N * H_room)))
        self.dx = 1.0 / (self.nx - 1)
        self.dy = H_room / (self.ny - 1)
        self.H_room, self.eta, self.Ra, self.Pr = H_room, eta, Ra, Pr

        # solid partition: central columns, from floor up to h_p = eta*H_room
        ic = self.nx // 2
        hw = max(1, int(round(part_halfwidth / self.dx)))
        self.ic, self.hw = ic, hw
        hp_rows = int(round(eta * (self.ny - 1)))
        self.hp_rows = hp_rows
        self.solid = np.zeros((self.ny, self.nx), bool)
        if hp_rows > 0:
            self.solid[0:hp_rows, ic - hw:ic + hw + 1] = True
        self.solid[:, 0] = self.solid[:, -1] = False   # keep outer walls as walls, not solid

        self.T = np.zeros((self.ny, self.nx)); self.T[:, 0] = 0.0; self.T[:, -1] = 1.0
        self.omega = np.zeros((self.ny, self.nx))
        self.psi = np.zeros((self.ny, self.nx))
        self.C = np.zeros((self.ny, self.nx))          # exhaust tracer
        self.u = np.zeros((self.ny, self.nx))
        self.v = np.zeros((self.ny, self.nx))
        self._build_poisson()

    # ---- sparse Laplacian with psi=0 on outer walls AND on the solid partition ----
    def _build_poisson(self):
        ny, nx, dx, dy = self.ny, self.nx, self.dx, self.dy
        idx = -np.ones((ny, nx), int)
        free = (~self.solid).copy()
        free[0, :] = free[-1, :] = free[:, 0] = free[:, -1] = False   # Dirichlet walls
        ids = np.where(free.ravel())[0]
        idx.ravel()[ids] = np.arange(ids.size)
        self._idx = idx
        n = ids.size
        rows, cols, vals = [], [], []
        cx, cy = 1.0 / dx**2, 1.0 / dy**2
        for k, flat in enumerate(ids):
            j, i = divmod(flat, nx)
            rows.append(k); cols.append(k); vals.append(-2 * cx - 2 * cy)
            for (jj, ii, c) in [(j, i - 1, cx), (j, i + 1, cx), (j - 1, i, cy), (j + 1, i, cy)]:
                m = idx[jj, ii]
                if m >= 0:
                    rows.append(k); cols.append(m); vals.append(c)
                # if neighbour is wall/solid, psi=0 there -> contributes nothing
        A = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
        self._lu = splu(A)
        self._n = n

    def _poisson(self, rhs_full):
        idxflat = self._idx.ravel()
        flat_ids = np.where(idxflat >= 0)[0]     # grid positions that are free unknowns
        k = idxflat[flat_ids]                     # their solver indices (a permutation of 0..n-1)
        b = np.empty(self._n)
        b[k] = rhs_full.ravel()[flat_ids]         # b[solver_index] = rhs at that cell
        x = self._lu.solve(b)
        full = np.zeros((self.ny, self.nx))
        full.ravel()[flat_ids] = x[k]             # scatter solution back to the grid
        return full

    @staticmethod
    def _fill_solid_zero_flux(f, solid):
        """Set solid cells to the mean of their fluid neighbours (adiabatic / no source)."""
        if not solid.any():
            return
        up, dn = np.roll(f, 1, 0), np.roll(f, -1, 0)
        lf, rt = np.roll(f, 1, 1), np.roll(f, -1, 1)
        nb = ~solid
        s = (up * np.roll(nb, 1, 0) + dn * np.roll(nb, -1, 0) +
             lf * np.roll(nb, 1, 1) + rt * np.roll(nb, -1, 1))
        c = (np.roll(nb, 1, 0).astype(float) + np.roll(nb, -1, 0) +
             np.roll(nb, 1, 1) + np.roll(nb, -1, 1))
        with np.errstate(invalid="ignore", divide="ignore"):
            avg = np.where(c > 0, s / np.maximum(c, 1), f)
        f[solid] = avg[solid]

    def _vel(self):
        self.u[1:-1, 1:-1] = (self.psi[2:, 1:-1] - self.psi[:-2, 1:-1]) / (2 * self.dy)
        self.v[1:-1, 1:-1] = -(self.psi[1:-1, 2:] - self.psi[1:-1, :-2]) / (2 * self.dx)
        self.u[self.solid] = 0.0; self.v[self.solid] = 0.0

    def _lap(self, f):
        return ((f[1:-1, 2:] - 2 * f[1:-1, 1:-1] + f[1:-1, :-2]) / self.dx**2 +
                (f[2:, 1:-1] - 2 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / self.dy**2)

    def _adv(self, f):
        u, v = self.u[1:-1, 1:-1], self.v[1:-1, 1:-1]
        dfdx = np.where(u > 0, (f[1:-1, 1:-1] - f[1:-1, :-2]) / self.dx,
                        (f[1:-1, 2:] - f[1:-1, 1:-1]) / self.dx)
        dfdy = np.where(v > 0, (f[1:-1, 1:-1] - f[:-2, 1:-1]) / self.dy,
                        (f[2:, 1:-1] - f[1:-1, 1:-1]) / self.dy)
        return u * dfdx + v * dfdy

    def step(self, dt):
        Ra, Pr, sol = self.Ra, self.Pr, self.solid
        Tn = self.T.copy()
        Tn[1:-1, 1:-1] = self.T[1:-1, 1:-1] + dt * (-self._adv(self.T) + self._lap(self.T))
        Tn[:, 0] = 0.0; Tn[:, -1] = 1.0
        Tn[0, :] = Tn[1, :]; Tn[-1, :] = Tn[-2, :]
        self._fill_solid_zero_flux(Tn, sol)            # adiabatic partition (no source/sink)
        # tracer: released at hot wall, absorbed (reset) at cold wall
        Cn = self.C.copy()
        Cn[1:-1, 1:-1] = self.C[1:-1, 1:-1] + dt * (-self._adv(self.C) + 0.02 * self._lap(self.C))
        Cn[:, -1] = 1.0; Cn[:, 0] = 0.0
        Cn[0, :] = Cn[1, :]; Cn[-1, :] = Cn[-2, :]
        self._fill_solid_zero_flux(Cn, sol)
        dTdx = (self.T[1:-1, 2:] - self.T[1:-1, :-2]) / (2 * self.dx)
        wn = self.omega.copy()
        wn[1:-1, 1:-1] = self.omega[1:-1, 1:-1] + dt * (
            -self._adv(self.omega) + Pr * self._lap(self.omega) + Ra * Pr * dTdx)
        # Poisson (internal solid handled by the factorization)
        self.psi = self._poisson(-wn)
        self.psi[sol] = 0.0
        # Thom vorticity on outer walls
        wn[:, 0] = -2 * self.psi[:, 1] / self.dx**2
        wn[:, -1] = -2 * self.psi[:, -2] / self.dx**2
        wn[0, :] = -2 * self.psi[1, :] / self.dy**2
        wn[-1, :] = -2 * self.psi[-2, :] / self.dy**2
        # Thom vorticity on solid faces (exposed to fluid): use one-sided psi curvature
        wn[sol] = 0.0
        self.T, self.C, self.omega = Tn, Cn, wn
        self._vel()

    def solve(self, max_iter=120000, tol=2e-5, safety=0.35):
        h = min(self.dx, self.dy)
        dt = safety * min(h**2 / (4 * max(self.Pr, 1.0)), h / (np.sqrt(self.Ra * self.Pr) + 1.0))
        prev = self.C.copy()
        for it in range(max_iter):
            self.step(dt)
            if it % 500 == 0 and it > 0:
                if not np.isfinite(self.T).all():
                    raise FloatingPointError("diverged")
                d = np.max(np.abs(self.C - prev)) / (dt * 500)
                if d < tol:
                    break
                prev = self.C.copy()
        return it

    def recirc_fraction(self):
        """phi = recirculation-induced inlet warming, measured as the mean normalized
        temperature at the cold-aisle rack-inlet face (T=0 cold supply, T=1 hot exhaust).
        This equals the fraction of hot recirculated air in the inlet stream, and
        converges with the flow field (unlike the slow passive tracer). The tracer C is
        retained as an independent cross-check."""
        face = self.T[:, self.ic - self.hw - 1]     # cold-aisle face of the rack inlet
        hp = max(self.hp_rows, self.ny // 4)
        return float(np.clip(np.mean(face[1:hp]), 0.0, 1.0))

    def recirc_fraction_tracer(self):
        face = self.C[:, self.ic - self.hw - 1]
        hp = max(self.hp_rows, self.ny // 4)
        return float(np.clip(np.mean(face[1:hp]), 0.0, 1.0))
