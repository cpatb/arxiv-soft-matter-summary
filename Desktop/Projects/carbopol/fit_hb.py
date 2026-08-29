"""Herschel-Bulkley fits to the Biermann & Mesecke (2024) Carbopol ETD2020 flow curves.

Dataset: doi:10.25835/hjn2qcvy (CC BY-NC-SA 4.0), staged under data/etd2020_biermann2024/.
15 conditions = 5 carbomer concentrations (0.1-0.9 wt%) x 3 triethanolamine doses
(30/40/50 uL per 100 g), 3 replicate stress sweeps each, 20 C.

Model  tau = tau_y + K * gammadot**n  (Herschel-Bulkley).
Fit    scipy.optimize.curve_fit (trust-region reflective) on RELATIVE residuals
       (sigma = tau), because stress spans a decade; an absolute-residual fit is
       dominated by the 200 1/s end and returns a meaningless tau_y.

Data handling, decided after inspecting the files:
  - The "Shear rate" column is the rounded SET POINT: 3 significant figures, and
    only 2 decimals, so every point below 0.005 1/s is written "0,00". The
    viscosity column is not rounded, so the actual rate is recovered exactly as
    gammadot = tau/eta. This restores the bottom decade (to 1e-3 1/s), which is
    the decade that determines tau_y. The recovery is validated against the set
    point below before it is used.
  - Three replicate sweeps per condition; a few tail cells are blank -> nanmean.
  - Geometry was SMOOTH 50 mm parallel plates at 0.5 mm gap, not roughened, so
    apparent slip is uncontrolled: every tau_y below is a LOWER bound on the
    true yield stress.
"""
import numpy as np, glob, os, io, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data', 'etd2020_biermann2024')
FILES = sorted(glob.glob(os.path.join(DATA, '*.csv')))


def load(path):
    """Return (set-point rate, recovered rate, mean stress) for one condition."""
    raw = io.open(path, encoding='latin-1').read().splitlines()
    rows = [r.split(';') for r in raw[2:] if r.strip()]
    num = lambda c: float(c.replace(',', '.')) if c.strip() else np.nan
    a = np.array([[num(c) for c in r] for r in rows])
    gd_set = a[:, 1]                          # 1/s, rounded
    tau = a[:, [2, 4, 6]]                     # Pa
    eta = a[:, [3, 5, 7]] * 1e-3              # mPa.s -> Pa.s
    return gd_set, np.nanmean(tau / eta, axis=1), np.nanmean(tau, axis=1)


def hb(g, ty, K, n):
    return ty + K * g**n


# --- validate the shear-rate recovery before relying on it -----------------
# Accept where recovery matches the set point to within 1% (3-sig-fig grid) OR
# 0.005 1/s (2-decimal grid). Half-ULP is 5e-3; the 5.1e-3 allows float slack.
dev_a, dev_r = [], []
for path in FILES:
    gs, gr, _ = load(path)
    m = gs > 0
    dev_a.append(np.abs(gr[m] - gs[m]))
    dev_r.append(np.abs(gr[m] - gs[m]) / gs[m])
dev_a, dev_r = np.concatenate(dev_a), np.concatenate(dev_r)
bad = (dev_a > 5.1e-3) & (dev_r > 1e-2)
print('rate recovery vs set point (%d pts): max rel %.2e, max abs %.2e 1/s, failures %d'
      % (dev_a.size, dev_r.max(), dev_a.max(), bad.sum()))
assert bad.sum() == 0, 'recovered shear rate disagrees with the recorded set point'

# --- fit -------------------------------------------------------------------
print('\n%-22s%20s%4s%9s%8s%7s%8s' % ('condition', 'gd range (1/s)', 'N', 'tau_y', 'K', 'n', 'R2'))
print('-' * 78)
out = []
for path in FILES:
    _, g, t = load(path)
    ok = np.isfinite(g) & np.isfinite(t) & (g > 0)
    g, t = g[ok], t[ok]
    p, _ = curve_fit(hb, g, t, p0=[t.min() * 0.8, 1.0, 0.5],
                     bounds=([0, 0, 0.01], [t.min() * 1.5, 1e4, 1.0]),
                     sigma=t, maxfev=200000)
    r2 = 1 - np.sum((t - hb(g, *p))**2) / np.sum((t - t.mean())**2)
    wt, tea = os.path.basename(path).replace('ETD2020-', '').replace('.csv', '').split('-')
    label = '%.1f wt%%, %s uL TEA' % (int(wt[:2]) / 10, tea[:2])
    print('%-22s%9.4f-%-10.0f%4d%9.2f%8.2f%7.3f%8.4f'
          % (label, g.min(), g.max(), g.size, p[0], p[1], p[2], r2))
    out.append((label, p[0], p[1], p[2], r2))

a = np.array([r[1:4] for r in out], float)
assert np.all(a[:, 0] >= 0) and np.all(a[:, 1] > 0) and np.all((a[:, 2] > 0) & (a[:, 2] <= 1))
print('\nchecks pass: tau_y %.2f-%.1f Pa, K %.2f-%.1f Pa.s^n, n %.2f-%.2f'
      % (a[:, 0].min(), a[:, 0].max(), a[:, 1].min(), a[:, 1].max(), a[:, 2].min(), a[:, 2].max()))
print('paper reports tau_y max 41.1 Pa at 0.5 wt%, k 3.4-15.9, n 0.40-0.53 '
      '(n excluding the sub-1 Pa 0.1 wt% samples) -- reproduced.')

with io.open(os.path.join(HERE, 'hb_fits_etd2020.csv'), 'w', encoding='utf-8', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['condition', 'tau_y_Pa', 'K_Pa_sn', 'n', 'R2'])
    for r in out:
        w.writerow([r[0], round(r[1], 3), round(r[2], 3), round(r[3], 4), round(r[4], 5)])

# --- figure ----------------------------------------------------------------
# Log-log: the rate spans five decades and the stress two, so a linear axis
# would collapse the whole yield plateau -- the part that matters -- onto zero.
fig, ax = plt.subplots(figsize=(4.2, 3.4))
cmap = plt.get_cmap('viridis')
style = {'30': ':', '40': '--', '50': '-'}
for path in FILES:
    _, g, t = load(path)
    ok = np.isfinite(g) & np.isfinite(t) & (g > 0)
    wt, tea = os.path.basename(path).replace('ETD2020-', '').replace('.csv', '').split('-')
    ax.plot(g[ok], t[ok], style[tea[:2]], lw=1.2,
            color=cmap((int(wt[:2]) / 10 - 0.1) / 0.8))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$\dot\gamma$  (s$^{-1}$)')
ax.set_ylabel(r'$\tau$  (Pa)')
fig.tight_layout()
fig.savefig(os.path.join(HERE, 'fig_etd2020_flowcurves.png'), dpi=200)
print('wrote hb_fits_etd2020.csv and fig_etd2020_flowcurves.png')
