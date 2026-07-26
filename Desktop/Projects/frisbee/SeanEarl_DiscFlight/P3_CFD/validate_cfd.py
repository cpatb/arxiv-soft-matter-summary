"""
CFD validation for P3. Two tests, both required before any flow field is
trusted downstream.

(a) ABL HORIZONTAL HOMOGENEITY -- Hargreaves & Wright (2007) protocol
    An equilibrium atmospheric boundary layer imposed at the inlet should
    advect across an EMPTY domain unchanged. Any drift between the inlet
    profile and the profile far downstream is spurious, caused by the
    inlet/wall-function/turbulence-model triple being mutually inconsistent,
    and it contaminates every result. We run a buildings-free case and compare
    U(z), k(z) at several fetches against the analytic Richards-Hoxey profile

        U(z) = (u*/kappa) ln((z + z0)/z0),   k = u*^2/sqrt(Cmu)

    and report the maximum deviation over the disc-altitude band.

(b) GRID INDEPENDENCE
    Two meshes differing ONLY in the fieldBox refinement level (a factor of 2
    in the region the paper reports) are compared at disc altitude. The
    manuscript's acceptance criterion is agreement in mean velocity at disc
    height within ~10%.

Both tests write a figure and print a PASS/FAIL table. Nothing here tunes
anything: it measures and reports.
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

KAPPA = 0.41
CMU = 0.09
FOAM_BASHRC = "/usr/lib/openfoam/openfoam2606/etc/bashrc"

# Disc-altitude band the paper cares about.
Z_DISC = (1.0, 6.0)

PROFILE_DICT = """/*--------------------------------*- C++ -*----------------------------------*/
FoamFile {{ version 2.0; format ascii; class dictionary; object sampleProfiles; }}

type            sets;
libs            (sampling);
setFormat       raw;
interpolationScheme cellPoint;
fields          (U k omega);

sets
{{
{sets}
}}
"""


def sample_profiles(case: Path, stations, z_max=120.0, n_z=120):
    """Sample vertical profiles at the given (x, y) stations."""
    blocks = []
    for i, (sx, sy) in enumerate(stations):
        blocks.append(f"""    prof{i}
    {{
        type    uniform;
        axis    z;
        start   ({sx} {sy} 0.5);
        end     ({sx} {sy} {z_max});
        nPoints {n_z};
    }}""")
    (case / "system" / "sampleProfiles").write_text(
        PROFILE_DICT.format(sets="\n".join(blocks)))

    cmd = (f"source {FOAM_BASHRC} && cd {case} && "
           "postProcess -func sampleProfiles -latestTime")
    r = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"profile sampling failed in {case}:\n"
                         f"{r.stdout[-2500:]}\n{r.stderr[-1500:]}")

    base = case / "postProcessing" / "sampleProfiles"
    tdir = sorted(base.glob("*"), key=lambda p: float(p.name))[-1]
    out = []
    for i in range(len(stations)):
        aU = np.loadtxt(sorted(tdir.glob(f"prof{i}_U.xy"))[0])
        ak = np.loadtxt(sorted(tdir.glob(f"prof{i}_k.xy"))[0])
        out.append(dict(z=aU[:, 2], U=aU[:, 3:6], k=ak[:, 3]))
    return out


def richards_hoxey(z, u_star, z0):
    return (u_star / KAPPA) * np.log((np.maximum(z, 0.0) + z0) / z0)


def test_abl_homogeneity(case: Path, outdir: Path):
    info = json.loads((case / "case_info.json").read_text())
    u_star, z0 = info["u_star"], info["z0_inlet"]
    dom = info["domain"]

    # Fetches measured from the inlet.
    x0 = dom["xmin"]
    fetches = [50.0, 300.0, 600.0, 900.0, 1400.0]
    stations = [(x0 + f, 0.0) for f in fetches]
    profs = sample_profiles(case, stations)

    print("\n" + "=" * 74)
    print("(a) ABL HORIZONTAL HOMOGENEITY  (Hargreaves & Wright 2007)")
    print("=" * 74)
    print(f"  empty domain, z0 = {z0} m, u* = {u_star:.4f} m/s")
    print(f"  {'fetch':>8s} {'max |dU|/U':>12s} {'max |dk|/k':>12s}   "
          f"(over z = {Z_DISC[0]}-{Z_DISC[1]} m)")

    k_ref = u_star ** 2 / np.sqrt(CMU)
    rows = []
    worst = 0.0
    for f, p in zip(fetches, profs):
        band = (p["z"] >= Z_DISC[0]) & (p["z"] <= Z_DISC[1])
        Umag = np.linalg.norm(p["U"], axis=1)
        Uref = richards_hoxey(p["z"], u_star, z0)
        with np.errstate(divide="ignore", invalid="ignore"):
            dU = np.abs(Umag - Uref) / np.maximum(Uref, 1e-9)
            dk = np.abs(p["k"] - k_ref) / k_ref
        mU, mk = float(np.nanmax(dU[band])), float(np.nanmax(dk[band]))
        worst = max(worst, mU)
        rows.append((f, p, mU, mk))
        print(f"  {f:7.0f}m {100*mU:11.1f}% {100*mk:11.1f}%")

    verdict = "PASS" if worst < 0.10 else ("MARGINAL" if worst < 0.20
                                           else "FAIL")
    print(f"  -> worst velocity deviation {100*worst:.1f}%  [{verdict}]")
    print("     (Hargreaves & Wright report a few % for a well-posed setup;")
    print("      >10% means the inlet/wall-function pair is inconsistent.)")

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    zz = np.linspace(0.5, 120, 300)
    ax[0].plot(richards_hoxey(zz, u_star, z0), zz, "k--", lw=2,
               label="Richards-Hoxey (analytic)")
    for f, p, _, _ in rows:
        ax[0].plot(np.linalg.norm(p["U"], axis=1), p["z"], lw=1,
                   label=f"fetch {f:.0f} m")
        ax[1].plot(p["k"], p["z"], lw=1, label=f"fetch {f:.0f} m")
    ax[1].axvline(k_ref, color="k", ls="--", lw=2, label="$u_*^2/\\sqrt{C_\\mu}$")
    for a, xl in zip(ax, ["|U| (m/s)", "k (m$^2$/s$^2$)"]):
        a.axhspan(Z_DISC[0], Z_DISC[1], color="tab:orange", alpha=0.15)
        a.set_xlabel(xl); a.set_ylabel("z (m)"); a.set_ylim(0, 120)
        a.grid(alpha=0.3); a.legend(fontsize=7)
    fig.suptitle("ABL horizontal homogeneity over an empty fetch\n"
                 f"worst deviation in the disc band: {100*worst:.1f}% "
                 f"[{verdict}]")
    fig.tight_layout()
    fig.savefig(outdir / "validation_abl_homogeneity.png", dpi=150)
    plt.close(fig)
    print(f"  [figure: {outdir/'validation_abl_homogeneity.png'}]")
    return dict(test="abl_homogeneity", worst_dU=worst, verdict=verdict,
                fetches=fetches,
                per_fetch=[dict(fetch=f, dU=a, dk=b) for f, _, a, b in rows])


def test_grid_independence(fine: Path, coarse: Path, outdir: Path):
    print("\n" + "=" * 74)
    print("(b) GRID INDEPENDENCE  (fieldBox refinement level 3 vs 2)")
    print("=" * 74)

    info_f = json.loads((fine / "case_info.json").read_text())
    info_c = json.loads((coarse / "case_info.json").read_text())
    print(f"  fine  : {info_f['resolution']:11s} field cell "
          f"{info_f['effective_field_cell_m'][0]:.2f} x "
          f"{info_f['effective_field_cell_m'][1]:.2f} m")
    print(f"  coarse: {info_c['resolution']:11s} field cell "
          f"{info_c['effective_field_cell_m'][0]:.2f} x "
          f"{info_c['effective_field_cell_m'][1]:.2f} m")

    # Profiles over the field itself, in SOLVER coordinates (both cases use
    # the same rotation for the same heading, so no rotation is needed here).
    stations = [(0.0, 0.0), (-30.0, 0.0), (30.0, 0.0), (0.0, -20.0),
                (0.0, 20.0)]
    pf = sample_profiles(fine, stations, z_max=20.0, n_z=80)
    pc = sample_profiles(coarse, stations, z_max=20.0, n_z=80)

    print(f"\n  {'station':>16s} {'|U| fine':>9s} {'|U| coarse':>11s} "
          f"{'diff':>8s}")
    diffs = []
    for (sx, sy), a, b in zip(stations, pf, pc):
        band = (a["z"] >= Z_DISC[0]) & (a["z"] <= Z_DISC[1])
        ua = float(np.mean(np.linalg.norm(a["U"], axis=1)[band]))
        ub = float(np.mean(np.linalg.norm(b["U"], axis=1)[band]))
        d = abs(ua - ub) / max(ua, 1e-9)
        diffs.append(d)
        print(f"  ({sx:6.1f},{sy:6.1f}) {ua:8.3f} {ub:10.3f} "
              f"{100*d:7.1f}%")

    worst = max(diffs)
    verdict = "PASS" if worst < 0.10 else ("MARGINAL" if worst < 0.20
                                           else "FAIL")
    print(f"  -> worst disc-band mean-velocity difference {100*worst:.1f}%  "
          f"[{verdict}]  (criterion: 10%)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 5))
    for i, ((sx, sy), a, b) in enumerate(zip(stations, pf, pc)):
        c = f"C{i}"
        ax.plot(np.linalg.norm(a["U"], axis=1), a["z"], c, lw=2,
                label=f"({sx:.0f},{sy:.0f}) fine")
        ax.plot(np.linalg.norm(b["U"], axis=1), b["z"], c, lw=1, ls="--",
                label=f"({sx:.0f},{sy:.0f}) coarse")
    ax.axhspan(Z_DISC[0], Z_DISC[1], color="tab:orange", alpha=0.15)
    ax.set_xlabel("|U| (m/s)"); ax.set_ylabel("z (m)"); ax.set_ylim(0, 20)
    ax.grid(alpha=0.3); ax.legend(fontsize=6, ncol=2)
    ax.set_title(f"Grid independence at disc altitude\n"
                 f"worst difference {100*worst:.1f}% [{verdict}]")
    fig.tight_layout()
    fig.savefig(outdir / "validation_grid_independence.png", dpi=150)
    plt.close(fig)
    print(f"  [figure: {outdir/'validation_grid_independence.png'}]")
    return dict(test="grid_independence", worst_dU=worst, verdict=verdict,
                stations=[list(s) for s in stations], diffs=diffs)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--empty", default="cases/dir_090_empty")
    ap.add_argument("--fine", default="cases/dir_090")
    ap.add_argument("--coarse", default="cases/dir_090_coarse")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--out-json", default="cfd_validation.json")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results = []

    if Path(args.empty).exists():
        results.append(test_abl_homogeneity(Path(args.empty), outdir))
    else:
        print(f"skipping ABL test: {args.empty} not found")

    if Path(args.fine).exists() and Path(args.coarse).exists():
        results.append(test_grid_independence(Path(args.fine),
                                              Path(args.coarse), outdir))
    else:
        print("skipping grid-independence test: need both meshes")

    Path(args.out_json).write_text(json.dumps(results, indent=2))
    print(f"\nwrote {args.out_json}")
    print("\nSUMMARY")
    for r in results:
        print(f"  {r['test']:22s} worst {100*r['worst_dU']:5.1f}%  "
              f"{r['verdict']}")


if __name__ == "__main__":
    main()
