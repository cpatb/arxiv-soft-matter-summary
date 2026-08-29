"""
Differentiable calibration of the Tier-0 Wyart-Cates closure (PyTorch autodiff).

This realises the top simulation-frontier recommendation: make the closure
differentiable so its parameters can be *calibrated by gradient descent* against
data, instead of hand-tuned. Here we (1) reimplement the steady flow curve in
torch, (2) generate synthetic stress-controlled "measurements" from known TRUE
parameters + noise, and (3) recover those parameters by backprop, validating
recovered-vs-true. This is the ground-truth test that must pass before trusting
the method on real rheometry.

Framework note: the numpy version lives in `wyart_cates.py`; this is the
differentiable twin. Same physics, gradients for free. (JAX/diffrax would add
autodiff *through the transient ODE* too; that's the next step.)

CST regime (phi < phi_m) is used so eta stays finite everywhere during the fit.
"""
from __future__ import annotations

import numpy as np
import torch
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

torch.set_default_dtype(torch.float64)  # rheology wants float64 (cf. JAX x64)
PURPLE, GOLD, INK = "#4E2A84", "#B6862C", "#2B2B2B"

PHI0 = 0.65          # frictionless jamming (known geometry; fixed)
ETA_S = 1.0e-3       # solvent viscosity (known; fixed)

# ground-truth parameters to recover
TRUE = dict(phi=0.520, phi_m=0.560, sigma_star=8.0)


def eta_app(sigma_h, phi, phi_m, sigma_star, phi0=PHI0, eta_s=ETA_S):
    """Differentiable steady apparent viscosity at applied (hydrodynamic) stress.
    eta_app = eta_r * eta_s, with the Wyart-Cates frictional-contact closure.
    """
    f = torch.exp(-sigma_star / sigma_h)              # fraction frictional
    phi_J = (1.0 - f) * phi0 + f * phi_m
    gap = torch.clamp(1.0 - phi / phi_J, min=1e-3)    # keep finite during fit
    eta_r = gap ** (-2.0)
    return eta_r * eta_s


def make_data(seed=0, n=40):
    """Synthetic stress-controlled measurement: eta_app vs applied stress + noise."""
    rng = np.random.default_rng(seed)
    sigma = torch.logspace(-1, 4, n)
    with torch.no_grad():
        eta = eta_app(sigma, TRUE["phi"], TRUE["phi_m"], TRUE["sigma_star"])
    # 5% multiplicative lognormal noise (viscometry is ~log-normal)
    noise = torch.tensor(rng.lognormal(mean=0.0, sigma=0.05, size=n))
    return sigma, eta * noise


class Closure(torch.nn.Module):
    """Learnable closure parameters, reparametrised to stay physical."""

    def __init__(self):
        super().__init__()
        # start deliberately OFF the true values
        self._phi = torch.nn.Parameter(torch.tensor(0.0))     # -> phi via sigmoid
        self._phi_m = torch.nn.Parameter(torch.tensor(0.0))
        self._log_sig = torch.nn.Parameter(torch.tensor(np.log(20.0)))  # true is 8

    def _bounded(self, raw, lo, hi):
        return lo + (hi - lo) * torch.sigmoid(raw)

    def params(self):
        return (self._bounded(self._phi, 0.40, 0.62),
                self._bounded(self._phi_m, 0.45, 0.63),
                torch.exp(self._log_sig))

    def forward(self, sigma):
        phi, phi_m, sig_star = self.params()
        return eta_app(sigma, phi, phi_m, sig_star)


def main():
    print("=" * 62)
    print("DIFFERENTIABLE CALIBRATION of the Tier-0 closure (torch)")
    print("=" * 62)
    sigma, eta_data = make_data()

    model = Closure()
    opt = torch.optim.Adam(model.parameters(), lr=0.05)
    log_eta_data = torch.log(eta_data)

    history = []
    for step in range(4000):
        opt.zero_grad()
        pred = model(sigma)
        loss = torch.mean((torch.log(pred) - log_eta_data) ** 2)  # log-space MSE
        loss.backward()
        opt.step()
        if step % 500 == 0 or step == 3999:
            history.append((step, loss.item()))

    with torch.no_grad():
        phi, phi_m, sig_star = (p.item() for p in model.params())

    print(f"\n  {'param':<12}{'true':>10}{'recovered':>12}{'abs err':>10}")
    rows = [("phi", TRUE["phi"], phi), ("phi_m", TRUE["phi_m"], phi_m),
            ("sigma_star", TRUE["sigma_star"], sig_star)]
    for name, t, r in rows:
        print(f"  {name:<12}{t:>10.3f}{r:>12.3f}{abs(t - r):>10.3f}")

    # gradient sanity: finite-difference check on sigma_star loss
    def loss_at_sigstar(val):
        with torch.no_grad():
            e = eta_app(sigma, torch.tensor(phi), torch.tensor(phi_m), torch.tensor(val))
            return torch.mean((torch.log(e) - log_eta_data) ** 2).item()
    fd = (loss_at_sigstar(sig_star + 1e-3) - loss_at_sigstar(sig_star - 1e-3)) / 2e-3
    print(f"\n  d(loss)/d(sigma_star) at optimum ~ {fd:+.2e}  (expect ~0)")

    # -------- validation --------
    ok_phi = abs(phi - TRUE["phi"]) < 0.02
    ok_phim = abs(phi_m - TRUE["phi_m"]) < 0.03
    ok_sig = abs(sig_star - TRUE["sigma_star"]) / TRUE["sigma_star"] < 0.20
    ok_grad = abs(fd) < 1e-2
    allpass = ok_phi and ok_phim and ok_sig and ok_grad

    # figure
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4), dpi=140)
    with torch.no_grad():
        fit = model(sigma)
    axA.loglog(sigma.numpy(), eta_data.numpy() * 1e3, "o", ms=5, color=GOLD,
               alpha=0.8, label="noisy 'measured' data")
    axA.loglog(sigma.numpy(), fit.numpy() * 1e3, "-", color=PURPLE, lw=2.4,
               label="calibrated closure")
    axA.set_xlabel(r"applied stress  $\sigma$  (Pa)", color=INK)
    axA.set_ylabel(r"apparent viscosity  (mPa$\cdot$s)", color=INK)
    axA.set_title("Gradient-descent calibration fit", fontsize=11)
    axA.legend(fontsize=9, frameon=False)
    axA.grid(True, which="both", color="#E7E7E7", lw=0.7)

    steps, losses = zip(*history)
    axB.semilogy(steps, losses, "-o", color=PURPLE, ms=4)
    axB.set_xlabel("Adam step", color=INK)
    axB.set_ylabel("log-space MSE loss", color=INK)
    axB.set_title("Convergence", fontsize=11)
    axB.grid(True, color="#E7E7E7", lw=0.7)
    for ax in (axA, axB):
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig("diff_calibration.png", facecolor="white")
    plt.close(fig)

    print("\n  wrote: diff_calibration.png")
    print("-" * 62)
    print(f"  CALIBRATION VALIDATION: {'PASS' if allpass else 'FAIL'}")
    print("-" * 62)


if __name__ == "__main__":
    main()
