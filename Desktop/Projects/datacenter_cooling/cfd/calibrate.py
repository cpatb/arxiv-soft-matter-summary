"""Fit the CFD recirculation data phi(eta), calibrate the Tier-0 model with it,
re-run the optimization, and compare to the placeholder-B result."""
import json, sys, os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(".."))
from dc_thermal import Hall, generic_hall          # Tier-0 model
from optimize_cooling import optimize

# ---- load CFD points, average over the (weak) Ra dependence ----
data = json.load(open("cfd_recirc.json"))
etas = sorted({d["eta"] for d in data})
phi = np.array([np.mean([d["phi"] for d in data if d["eta"] == e]) for e in etas])
etas = np.array(etas)

# ---- fit phi(eta) = floor + amp * exp(-eta/tau) ----
def model(e, floor, amp, tau):
    return floor + amp * np.exp(-e / tau)
(floor, amp, tau), _ = curve_fit(model, etas, phi, p0=[0.07, 0.17, 0.4],
                                 bounds=([0, 0, 0.05], [0.3, 0.5, 5]))
print(f"CFD fit: phi(eta) = {floor:.3f} + {amp:.3f} exp(-eta/{tau:.3f})")
print(f"  phi(0)={model(0,floor,amp,tau):.3f}  phi(1)={model(1,floor,amp,tau):.3f}  "
      f"(placeholder linear went to 0 at eta=1)")

def phi_cfd(eta):
    return float(model(np.clip(eta, 0, 1), floor, amp, tau))

# ---- CFD-calibrated Hall: override the recirculation base with phi_cfd ----
class CalibratedHall(Hall):
    def influence_matrix(self, provisioning):
        p = float(provisioning)
        under = max(0.0, 1.0 - p)
        phi = phi_cfd(self.containment) * (1.0 + under)      # CFD base, keep under-provisioning amp
        phi = np.clip(phi, 0.0, 0.95)
        return phi * self._W

def make(cls, **kw):
    import numpy as np
    rng = np.random.default_rng(1)
    base = np.full(20, 200.0 / 20)
    load = base * rng.uniform(0.6, 1.0, 20)
    hot = rng.choice(20, 3, replace=False); load[hot] *= 1.8
    load *= 200.0 / load.sum()
    h = cls.__new__(cls); Hall.__init__(h, load, containment=kw["containment"])
    return h

# ---- compare optimized PUE vs containment: placeholder vs CFD-calibrated ----
etar = np.linspace(0, 0.95, 8)
pue_old, pue_new = [], []
for e in etar:
    ho = make(Hall, containment=e); pue_old.append(optimize(ho).fun)
    hn = make(CalibratedHall, containment=e); pue_new.append(optimize(hn).fun)
pue_old, pue_new = np.array(pue_old), np.array(pue_new)

print("\n containment | PUE placeholder | PUE CFD-calibrated")
for e, a, b in zip(etar, pue_old, pue_new):
    print(f"    {e:4.2f}     |     {a:.3f}      |     {b:.3f}")

# ---- figures ----
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ee = np.linspace(0, 1, 100)
ax[0].plot(etas, phi, "o", ms=7, label="CFD (de Vahl Davis-validated solver)")
ax[0].plot(ee, model(ee, floor, amp, tau), "-", label="fit: floor + amp·e$^{-\\eta/\\tau}$")
ax[0].plot(ee, 0.28 * (1 - ee), "--", color="grey", label="Tier-0 placeholder  0.28(1-η)")
phi1 = model(1.0, floor, amp, tau)
ax[0].axhline(phi1, ls=":", color="C1", alpha=0.7)
ax[0].annotate(f"~{phi1:.02f} residual recirc at η=1\n(placeholder wrongly →0)",
               (0.95, phi1), xytext=(0.45, 0.155), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=0.7))
ax[0].set_xlabel("containment effectiveness  η"); ax[0].set_ylabel("recirculation fraction  φ")
ax[0].set_title("CFD-calibrated recirculation vs the placeholder"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.25)

ax[1].plot(etar, pue_old, "s--", color="grey", label="placeholder 0.28(1-η)")
ax[1].plot(etar, pue_new, "o-", color="C3", label="CFD-calibrated φ(η)")
ax[1].set_xlabel("containment effectiveness  η"); ax[1].set_ylabel("optimized PUE")
ax[1].set_title("Optimized PUE: the calibration changes the containment ROI")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.25)
fig.tight_layout(); fig.savefig("fig_cfd_calibration.png", dpi=140)

json.dump(dict(floor=floor, amp=amp, tau=tau,
               phi0=model(0, floor, amp, tau), phi1=model(1, floor, amp, tau),
               pue_placeholder=list(map(float, pue_old)),
               pue_calibrated=list(map(float, pue_new)), eta=list(map(float, etar))),
          open("cfd_calibration.json", "w"), indent=2)
print("\nwrote fig_cfd_calibration.png, cfd_calibration.json")
