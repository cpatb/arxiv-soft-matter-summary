"""Due diligence on a CANDIDATE dataset. This creates no record and ships no verdict.

`data/dmem_sints2024/` holds the first openly-licensed (CC BY 4.0) rheometer dataset found
for this project (NEXT_STEPS item 20). Before anything is built on it, the question is
whether its curves can support a record at all -- and the source paper makes a claim this
data can be asked about directly.

**The paper's claim** (Sints et al., arXiv:2409.13650, abstract, verbatim): DMEM
supplemented with FBS *"requires a more complex treatment (i.e., is a non-Newtonian
fluid)"*, with *"observed shear stress values significantly exceeding literature
reports"*.

**The observation that motivated this scan.** In every interval the SHEAR STRESS column is
close to affine in shear rate with a non-zero intercept: `tau ~ tau_0 + eta_p gd`. A
viscosity computed as `tau/gd` from such a series falls like `tau_0/gd + eta_p`, which
looks exactly like shear thinning on a log-log plot whatever `tau_0` is. So the question
this scan asks is narrow and answerable:

    does a two-parameter BINGHAM form (tau_0, eta_p) describe these curves as well as, or
    better than, a two-parameter shear-thinning POWER LAW -- and what plastic viscosity
    is left once the intercept is taken out?

That is not the same question as "is the fluid non-Newtonian". A Bingham fluid IS
non-Newtonian, so a Bingham fit **agrees** with the paper's headline. What it would change
is the MECHANISM and everything downstream of it: a stress intercept and a shear-thinning
exponent extrapolate to a microfluidic channel completely differently.

**THE ANSWER IS NO, AND IT REFUTES THE HYPOTHESIS THAT MOTIVATED THIS SCAN.** The power
law fits better than Bingham on **7 of the 8 intervals**, on the same error measure and at
the same parameter count. The affine-looking stress column is a two-parameter power law
with `n` around 0.3-0.7, not a yield stress plus a plastic viscosity. The paper's
shear-thinning reading is the one the data supports and mine was not; this docstring keeps
the refuted hypothesis rather than deleting it, because the scan is only worth reading if
what it was looking for is visible.

Since `tau_0` belongs to the losing model it is not a yield stress and must not be quoted
as one. The question of what it would have been -- real serum structure, a free-surface
effect at the cone edge, or an instrument offset -- does not arise, and the scan reports
the numbers and stops.

The paper was read on 2026-09-08 and `../data/dmem_sints2024/SOURCE.md` records what it
settles: room temperature is 20 degrees C, the authors fit this same power law and publish
their own indices (which this scan reproduces without having seen them), the startup
transient is handled by a declared 2000 s pre-shear, and the high low-rate viscosity is
the paper's headline result rather than a defect.

Run:  python scripts/dmem_candidate_scan.py
Writes `data/dmem_sints2024/candidate_scan.md`. No record is created.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "dmem_sints2024"
ZIP = DATA / "DMEM_FBS_rheology.zip"
OUT = DATA / "candidate_scan.md"

# The instrument writes a header block, then a `Meas. Pts.` table, and a file may hold
# more than one interval. Every interval is parsed; none is dropped for looking wrong.
_HDR = re.compile(r"^Meas\. Pts\.\s")


def intervals(text: str) -> list[dict]:
    lines = text.splitlines()
    out, i = [], 0
    profile = None
    while i < len(lines):
        if "d(gamma)/dt" in lines[i]:
            profile = " ".join(lines[i].split("\t")[-1].split()) or profile
        if _HDR.match(lines[i]):
            rows = []
            j = i + 2                      # skip the units line
            while j < len(lines):
                c = lines[j].split("\t")
                if len(c) < 4 or not c[0].strip().isdigit():
                    break
                try:
                    rows.append((float(c[1]), float(c[2]), float(c[3]), float(c[5])))
                except (ValueError, IndexError):
                    break
                j += 1
            if rows:
                a = np.array(rows, dtype=float)
                out.append({"gd": a[:, 0], "tau": a[:, 1], "eta": a[:, 2],
                            "torque_uNm": a[:, 3], "profile": profile})
            i = j
        else:
            i += 1
    return out


def _rmse_log(eta_obs, eta_mod):
    return float(np.sqrt(np.mean((np.log(eta_obs) - np.log(eta_mod)) ** 2)))


def fits(gd, tau) -> dict:
    """Three two-or-fewer-parameter descriptions, compared on the SAME quantity
    (`rmse` in log viscosity), so the comparison is not a units artefact."""
    eta = tau / gd
    out = {}
    # Newtonian: one parameter, the null model
    mu = float(np.mean(tau / gd))
    out["newtonian"] = {"params": {"mu": mu}, "rmse_log_eta": _rmse_log(eta, mu + 0 * gd)}
    # Bingham: tau = tau_0 + eta_p gd
    (t0, ep), _ = curve_fit(lambda x, a, b: a + b * x, gd, tau,
                            p0=[float(tau.min()), float(np.median(tau / gd))], maxfev=20000)
    out["bingham"] = {"params": {"tau_0_Pa": float(t0), "eta_p_Pa_s": float(ep)},
                      "rmse_log_eta": _rmse_log(eta, t0 / gd + ep)}
    # power law: tau = k gd^n
    (k, n), _ = curve_fit(lambda x, a, b: a * x ** b, gd, tau,
                          p0=[float(np.median(tau)), 0.5], maxfev=20000)
    out["power_law"] = {"params": {"k": float(k), "n": float(n)},
                        "rmse_log_eta": _rmse_log(eta, k * gd ** (n - 1.0))}
    return out


def main() -> int:
    z = zipfile.ZipFile(ZIP)
    names = sorted(n for n in z.namelist() if n.endswith(".txt") and "readme" not in n)
    blocks, rows = [], []
    for n in names:
        for idx, iv in enumerate(intervals(z.read(n).decode("utf-8", "replace"))):
            # the instrument writes a first point at essentially zero rate in one file
            # (gd = 2.24e-7 1/s, eta = 8710 Pa s). It is EXCLUDED, by a stated rule --
            # gd below 1e-3 1/s -- and the exclusion is counted, never silent.
            keep = iv["gd"] > 1e-3
            n_drop = int((~keep).sum())
            gd, tau = iv["gd"][keep], iv["tau"][keep]
            if gd.size < 4:
                continue
            f = fits(gd, tau)
            label = f"{Path(n).stem} #{idx + 1}"
            rows.append({
                "label": label, "n_points": int(gd.size), "n_dropped_below_1e-3": n_drop,
                "gd_min": float(gd.min()), "gd_max": float(gd.max()),
                "eta_min": float((tau / gd).min()), "eta_max": float((tau / gd).max()),
                "torque_min_uNm": float(iv["torque_uNm"][keep].min()),
                "profile": iv["profile"], **f})
            blocks.append(label)

    best = [r for r in rows if r["bingham"]["rmse_log_eta"] < r["power_law"]["rmse_log_eta"]]
    L = [
        "# DMEM + FBS (Sints et al. 2024) — candidate scan, NOT a registry record",
        "",
        "Written by `scripts/dmem_candidate_scan.py`. **No record is created by this file "
        "and no verdict is shipped.** It exists so that a later session decides on evidence "
        "rather than on a first impression. Read the script's docstring before the table.",
        "",
        "## Why it was scanned",
        "",
        "This is the first openly-licensed (CC BY 4.0) rheometer dataset found for this "
        "project. Its source paper states that DMEM + FBS *\"requires a more complex "
        "treatment (i.e., is a non-Newtonian fluid)\"*. In every interval the **shear stress "
        "column is close to affine in shear rate with a non-zero intercept**, and a "
        "viscosity formed as `tau/gd` from such a series falls like `tau_0/gd + eta_p` — "
        "which looks like shear thinning on a log-log plot whatever `tau_0` is. So the "
        "narrow, answerable question is whether a two-parameter **Bingham** form describes "
        "these curves as well as a two-parameter shear-thinning **power law**.",
        "",
        "**A Bingham fluid is also non-Newtonian, so this does not contradict the paper's "
        "headline.** What it would change is the mechanism, and everything downstream: a "
        "stress intercept and a thinning exponent extrapolate to a channel completely "
        "differently.",
        "",
        f"## Result: **the hypothesis is REFUTED.** Bingham beats the power law on only "
        f"{len(best)} of {len(rows)} intervals",
        "",
        f"The scan was run to test whether an affine stress column was masquerading as "
        f"shear thinning. It is not. **The two-parameter power law fits better on "
        f"{len(rows) - len(best)} of {len(rows)} intervals**, at the same parameter count "
        f"and on the same error measure, with `n` between "
        f"{min(r['power_law']['params']['n'] for r in rows):.2f} and "
        f"{max(r['power_law']['params']['n'] for r in rows):.2f}. "
        "**The paper's shear-thinning reading is the one this data supports.** The "
        "hypothesis that motivated the scan is kept in the script's docstring rather than "
        "deleted, because a scan is only worth reading if what it was looking for is visible.",
        "",
        "All three descriptions are compared on the same quantity — rms error in "
        "log viscosity — so the ranking is not a units artefact. `newtonian` is the "
        "one-parameter null; the other two have two parameters each. Both two-parameter "
        f"forms beat the Newtonian null on {sum(1 for r in rows if min(r['bingham']['rmse_log_eta'], r['power_law']['rmse_log_eta']) < r['newtonian']['rmse_log_eta'])} "
        f"of {len(rows)}, which is the least this data could have shown and does not by "
        "itself establish a form.",
        "",
        "| interval | pts | rate range (1/s) | eta range (Pa s) | Newtonian rmse | "
        "Bingham rmse | power-law rmse | `tau_0` (Pa) | `eta_p` (Pa s) | `n` |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        L.append(
            f"| `{r['label']}` | {r['n_points']} | {r['gd_min']:.3g}–{r['gd_max']:.3g} | "
            f"{r['eta_min']:.3g}–{r['eta_max']:.3g} | "
            f"{r['newtonian']['rmse_log_eta']:.3f} | "
            f"**{r['bingham']['rmse_log_eta']:.3f}** | "
            f"{r['power_law']['rmse_log_eta']:.3f} | "
            f"{r['bingham']['params']['tau_0_Pa']:.4g} | "
            f"{r['bingham']['params']['eta_p_Pa_s']:.4g} | "
            f"{r['power_law']['params']['n']:.3f} |")

    eps = [r["bingham"]["params"]["eta_p_Pa_s"] for r in rows]
    t0s = [r["bingham"]["params"]["tau_0_Pa"] for r in rows]
    dropped = sum(r["n_dropped_below_1e-3"] for r in rows)
    L += [
        "",
        f"`{dropped}` point(s) were excluded across all intervals by a stated rule — shear "
        "rate below 1e-3 1/s — and the count is printed rather than the exclusion being "
        "silent. One file records a first point at 2.24e-7 1/s with an apparent viscosity of "
        "8710 Pa s.",
        "",
        "## What the numbers say, and what they do not",
        "",
        f"- The fitted intercept `tau_0` runs **{min(t0s):.4g} to {max(t0s):.4g} Pa** — "
        "millipascals, which is small in absolute terms and is most of the stress at the "
        "bottom of a rate range that starts at 0.01 1/s.",
        f"- The plastic viscosity left after removing it runs **{min(eps):.4g} to "
        f"{max(eps):.4g} Pa s**. Water is about 1e-3 Pa s, so the coldest, most "
        f"serum-rich interval still sits about {max(eps) / 1e-3:.0f}x water while the "
        f"warmest sits about {min(eps) / 1e-3:.0f}x.",
        "- `tau_0` is **not constant across the four samples**: it tracks FBS content and "
        "falls with temperature. That would have argued against a fixed instrument offset — "
        "but since the Bingham form is the WORSE description here, `tau_0` is a fitted "
        "parameter of a rejected model and should not be quoted as a yield stress at all.",
        "",
        "**The magnitude is the PAPER'S OWN HEADLINE, not an oddity this scan found.** "
        "The room-temperature intervals sit 14–32x water and the 37 °C ones 2–5x, for a "
        "fluid that is ~99 % water. Reading the source (2026-09-08) settles it: the paper "
        "reports that *\"while the measured values do converge to the literature value as the "
        "shear rate increases, the viscosity at lower shear rates is much grater [sic] than "
        "the literature value\"*, against a baseline of 9.3e-4 Pa s. That gap is the result "
        "the paper exists to report.",
        "",
        "**The scan's fitted indices reproduce the paper's own, without having seen them.** "
        "The authors fit the Ostwald-de Waele power law and quote n = 0.32 (10 % FBS at "
        "20 °C), 0.31 (10 % and 5 % at 20 °C) and 0.5–0.54 (10 % at 37 °C, and 1 % at "
        "20 °C). This scan gets 0.301/0.302 for 10 % at 20 °C, 0.345/0.310 for 5 % at "
        "20 °C, and 0.492 and 0.507 for the two the paper groups at 0.5–0.54. That "
        "agreement is the best available evidence that the parse, the interval splitting and "
        "the zero-rate exclusion are all correct. **The one that agrees least is named "
        "rather than omitted**: the 37 °C 20-point interval gives 0.700, and it is the only "
        "interval running to 10 1/s rather than 1 1/s, so it is fitted over a decade the "
        "others do not cover and a single power law is not expected to hold across both.",
        "",
        "**What a later session should do:** the form is settled (`power_law`, the family the "
        "authors themselves fit). Fit it through `biorheo.identify` for intervals and an "
        "identifiability verdict, pre-register a gate against the paper's published indices, "
        "and write the record. See `SOURCE.md` for the ordered list and for the pre-shear "
        "protocol any record must inherit.",
        "",
        "## The other reason this dataset cannot carry an injectability claim",
        "",
        "The measured rates stop at 1 1/s (10 1/s on one interval). A fine needle imposes "
        "1e4–1e5 1/s. Any force computed from a record built here would be **four to five "
        "decades outside the measured range**, against 1.42 decades for the blood records. "
        "That makes it an excellent demonstration of the extrapolation flag and a poor basis "
        "for a number anyone should act on.",
        "",
    ]
    OUT.write_text("\n".join(L), encoding="utf-8")
    print(f"scanned {len(rows)} intervals; Bingham beats power law on {len(best)}; "
          f"wrote {OUT.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
