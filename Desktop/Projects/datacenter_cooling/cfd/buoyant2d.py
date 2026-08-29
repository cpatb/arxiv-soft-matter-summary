"""
buoyant2d.py — compact 2D Boussinesq buoyant-convection solver
(streamfunction-vorticity), used to calibrate the data-hall recirculation matrix.

Non-dimensional governing equations (Pr = Prandtl, Ra = Rayleigh):
    vorticity:   d(omega)/dt + u.grad(omega) = Pr * lap(omega) + Ra*Pr * dT/dx
    energy:      d(T)/dt     + u.grad(T)      = lap(T)
    streamfn:    lap(psi) = -omega ,   u = d(psi)/dy ,  v = -d(psi)/dx

The Poisson solve for psi (psi=0 on all walls) uses a fast 2D discrete sine
transform, so each step is O(N^2 log N). Time-marched to steady state.

Validated in test_buoyant2d.py against the de Vahl Davis (1983) differentially
heated square cavity Nusselt-number benchmarks.
"""
from __future__ import annotations
import numpy as np
from scipy.fft import dstn, idstn


def _poisson_dst(rhs, dx, dy):
    """Solve lap(psi) = rhs on interior with psi=0 Dirichlet, via 2D DST-I."""
    ny, nx = rhs.shape                      # interior points
    rhat = dstn(rhs, type=1)
    i = np.arange(1, nx + 1)
    j = np.arange(1, ny + 1)
    lam_x = (2 * np.cos(np.pi * i / (nx + 1)) - 2) / dx**2
    lam_y = (2 * np.cos(np.pi * j / (ny + 1)) - 2) / dy**2
    denom = lam_x[None, :] + lam_y[:, None]
    psihat = rhat / denom
    return idstn(psihat, type=1)


class Cavity:
    """Differentially heated cavity on [0,1]x[0,H]; hot left wall, cold right wall,
    adiabatic top/bottom. Used for solver validation."""

    def __init__(self, N=64, H=1.0, Ra=1e4, Pr=0.71):
        self.N, self.H, self.Ra, self.Pr = N, H, Ra, Pr
        self.nx = N
        self.ny = max(4, int(round(N * H)))
        self.dx = 1.0 / (self.nx - 1)
        self.dy = H / (self.ny - 1)
        self.T = np.zeros((self.ny, self.nx))
        self.T[:, 0] = 1.0                    # hot left
        self.omega = np.zeros((self.ny, self.nx))
        self.psi = np.zeros((self.ny, self.nx))
        self.u = np.zeros((self.ny, self.nx))
        self.v = np.zeros((self.ny, self.nx))

    def _velocities(self):
        self.u[1:-1, 1:-1] = (self.psi[2:, 1:-1] - self.psi[:-2, 1:-1]) / (2 * self.dy)
        self.v[1:-1, 1:-1] = -(self.psi[1:-1, 2:] - self.psi[1:-1, :-2]) / (2 * self.dx)

    def _lap(self, f):
        return ((f[1:-1, 2:] - 2 * f[1:-1, 1:-1] + f[1:-1, :-2]) / self.dx**2 +
                (f[2:, 1:-1] - 2 * f[1:-1, 1:-1] + f[:-2, 1:-1]) / self.dy**2)

    def _upwind(self, f):
        u, v = self.u[1:-1, 1:-1], self.v[1:-1, 1:-1]
        dfdx = np.where(u > 0, (f[1:-1, 1:-1] - f[1:-1, :-2]) / self.dx,
                        (f[1:-1, 2:] - f[1:-1, 1:-1]) / self.dx)
        dfdy = np.where(v > 0, (f[1:-1, 1:-1] - f[:-2, 1:-1]) / self.dy,
                        (f[2:, 1:-1] - f[1:-1, 1:-1]) / self.dy)
        return u * dfdx + v * dfdy

    def step(self, dt):
        Ra, Pr = self.Ra, self.Pr
        # --- energy ---
        Tn = self.T.copy()
        Tn[1:-1, 1:-1] = self.T[1:-1, 1:-1] + dt * (-self._upwind(self.T) + self._lap(self.T))
        Tn[:, 0] = 1.0; Tn[:, -1] = 0.0             # hot/cold walls
        Tn[0, :] = Tn[1, :]; Tn[-1, :] = Tn[-2, :]  # adiabatic top/bottom
        # --- vorticity transport (buoyancy source ~ dT/dx) ---
        dTdx = (self.T[1:-1, 2:] - self.T[1:-1, :-2]) / (2 * self.dx)
        wn = self.omega.copy()
        wn[1:-1, 1:-1] = self.omega[1:-1, 1:-1] + dt * (
            -self._upwind(self.omega) + Pr * self._lap(self.omega) + Ra * Pr * dTdx)
        # --- streamfunction ---
        self.psi[1:-1, 1:-1] = _poisson_dst(-wn[1:-1, 1:-1], self.dx, self.dy)
        # --- Thom's vorticity wall BC ---
        wn[:, 0] = -2 * self.psi[:, 1] / self.dx**2
        wn[:, -1] = -2 * self.psi[:, -2] / self.dx**2
        wn[0, :] = -2 * self.psi[1, :] / self.dy**2
        wn[-1, :] = -2 * self.psi[-2, :] / self.dy**2
        self.T, self.omega = Tn, wn
        self._velocities()

    def solve(self, max_iter=200000, tol=1e-5, safety=0.4):
        h = min(self.dx, self.dy)
        # explicit-stability CFL: min of diffusive and advective/buoyancy limits.
        # nondim velocity scale ~ sqrt(Ra*Pr); source term ~ Ra*Pr.
        dt_diff = h**2 / (4.0 * max(self.Pr, 1.0))
        dt_adv = h / (np.sqrt(self.Ra * self.Pr) + 1.0)
        dt = safety * min(dt_diff, dt_adv)
        prevT = self.T.copy()
        for it in range(max_iter):
            self.step(dt)
            if it % 500 == 0 and it > 0:
                if not np.isfinite(self.T).all():
                    raise FloatingPointError("solver diverged")
                d = np.max(np.abs(self.T - prevT)) / (dt * 500)
                if d < tol:
                    break
                prevT = self.T.copy()
        return it

    def nusselt_hot(self):
        """Average Nusselt number on the hot wall: Nu = -integral dT/dx |_{x=0} dy."""
        dTdx = (self.T[:, 1] - self.T[:, 0]) / self.dx
        return -np.trapezoid(dTdx, dx=self.dy) / self.H
