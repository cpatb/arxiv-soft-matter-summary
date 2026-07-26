"""
Test the manuscript's stated hypotheses against the computed results.

The draft's "Expected results / placeholder findings" section makes three
falsifiable predictions. This script evaluates each against the actual CFD +
Monte Carlo output and reports SUPPORTED / CONTRADICTED / NOT TESTABLE, with
the numbers behind the verdict. A contradiction is a finding, not a failure,
and is reported as such.

H1. Easterly cases (90 +/- 30 deg) show a sharp Q gradient across the field's
    eastern half corresponding to a developing internal boundary layer;
    backhand throws from the windward sideline are most consistent; hammer
    throws systematically degraded.

H2. Westerly cases (270 +/- 30 deg) show patchier Q with sheltered "throw
    lanes" downstream of the Norville Center; forehand throws benefit most.

H3. The best-condition wind direction for canonical 25 m throws is a light
    SSW (~210 deg), aligning with prevailing climatology.
"""

from __future__ import annotations

import argparse
import json

import numpy as np


def verdict(ok, msg):
    tag = {True: "SUPPORTED", False: "CONTRADICTED", None: "NOT TESTABLE"}[ok]
    return f"[{tag:13s}] {msg}"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Q", default="Q.npz")
    ap.add_argument("--fields", default="../P3_CFD/flow_fields.npz")
    ap.add_argument("--landings", default="../P4_Coupling/landing_stats.npz")
    ap.add_argument("--out-json", default="hypothesis_check.json")
    args = ap.parse_args()

    Qd = np.load(args.Q, allow_pickle=False)
    F = np.load(args.fields, allow_pickle=False)
    L = np.load(args.landings, allow_pickle=False)

    Q = Qd["Q"]
    x, y, headings = Qd["x"], Qd["y"], Qd["headings"]
    throws = [str(t) for t in L["throws"]]
    results = []

    print("=" * 78)
    print("HYPOTHESIS CHECK against the draft's predicted findings")
    print("=" * 78)
    print(f"computed headings: {list(int(h) for h in headings)}")
    print(f"throw types: {throws}")

    # ---------------- H1: easterly Q gradient across the field ------------
    print("\n" + "-" * 78)
    print("H1  easterly: sharp Q gradient across the field's eastern half")
    print("-" * 78)
    h_east = [i for i, h in enumerate(headings) if 60 <= int(h) <= 120]
    if not h_east:
        results.append(dict(hypothesis="H1", verdict="NOT TESTABLE",
                            reason="no easterly direction computed"))
        print(verdict(None, "no easterly case computed"))
    else:
        ih = h_east[0]
        Qh = Q[..., ih]
        # East half vs west half of the launch grid.
        mid = len(x) // 2
        west, east = np.nanmean(Qh[:mid]), np.nanmean(Qh[mid:])
        # Gradient magnitude along x, normalised by the field mean.
        gx = np.abs(np.gradient(np.nanmean(Qh, axis=1), x))
        rel_grad = float(np.nanmax(gx) * (x[-1] - x[0]) /
                         max(np.nanmean(Qh), 1e-12))
        ok = rel_grad > 0.5
        msg = (f"theta={int(headings[ih])}: Q west half {west:.4f}, "
               f"east half {east:.4f}; max normalised gradient "
               f"{rel_grad:.2f} across the field")
        print(verdict(ok, msg))
        results.append(dict(hypothesis="H1a_gradient",
                            verdict="SUPPORTED" if ok else "CONTRADICTED",
                            Q_west=float(west), Q_east=float(east),
                            normalised_gradient=rel_grad))

        # Mean speed at 2 m across the field: is there an IBL-like ramp?
        key = f"dir_{int(headings[ih]):03d}"
        if f"{key}__U" in F.files:
            zk = int(np.argmin(np.abs(F["z"] - 2.0)))
            mag = np.linalg.norm(F[f"{key}__U"][:, :, zk, :], axis=-1)
            prof = np.nanmean(mag, axis=1)
            print(f"    |U| at 2 m along x: {prof.min():.2f} .. "
                  f"{prof.max():.2f} m/s "
                  f"(west end {prof[0]:.2f}, east end {prof[-1]:.2f})")
            results.append(dict(hypothesis="H1_speed_profile",
                                U_west_end=float(prof[0]),
                                U_east_end=float(prof[-1]),
                                U_min=float(prof.min()),
                                U_max=float(prof.max())))

        # Hammer degraded?
        if "hammer" in throws and "right_backhand" in throws:
            sr = L["success_rate"]
            ih_t = throws.index("hammer")
            ib_t = throws.index("right_backhand")
            s_h = float(np.nanmean(sr[:, :, ih, ih_t, :]))
            s_b = float(np.nanmean(sr[:, :, ih, ib_t, :]))
            # The hammer's still-air range is 12.1 m; the scoring target is
            # 25 m downrange. Its success rate is therefore structurally zero
            # whatever the wind does, so this comparison cannot test the
            # hypothesis -- it measures the target choice.
            print(verdict(None, f"hammer success {s_h:.4f} vs backhand "
                                f"{s_b:.4f} -- NOT INTERPRETABLE: the "
                                f"hammer's 12.1 m still-air range cannot "
                                f"reach the 25 m scoring target, so 0% "
                                f"measures the target, not the wind"))
            results.append(dict(hypothesis="H1b_hammer_degraded",
                                verdict="NOT TESTABLE",
                                reason="hammer still-air range 12.1 m vs 25 m "
                                       "scoring target; success rate is "
                                       "structurally zero",
                                hammer=s_h, backhand=s_b))

    # ---------------- H2: westerly sheltered lanes, forehand gains --------
    print("\n" + "-" * 78)
    print("H2  westerly: patchier Q, sheltered lanes, forehand benefits most")
    print("-" * 78)
    h_west = [i for i, h in enumerate(headings) if 240 <= int(h) <= 300]
    if not h_west:
        print(verdict(None, "no westerly direction computed"))
        results.append(dict(hypothesis="H2", verdict="NOT TESTABLE"))
    else:
        iw = h_west[0]
        Qw = Q[..., iw]
        cv_w = float(np.nanstd(Qw) / max(np.nanmean(Qw), 1e-12))
        cvs = {int(headings[i]): float(np.nanstd(Q[..., i]) /
                                       max(np.nanmean(Q[..., i]), 1e-12))
               for i in range(len(headings))}
        ok = cv_w >= max(cvs.values()) - 1e-9
        print(verdict(ok, f"spatial CV of Q by heading: "
                          + ", ".join(f"{k}:{v:.3f}" for k, v in cvs.items())))
        results.append(dict(hypothesis="H2a_patchier",
                            verdict="SUPPORTED" if ok else "CONTRADICTED",
                            cv_by_heading=cvs))

        if "right_forehand" in throws and "right_backhand" in throws:
            sr = L["success_rate"]
            i_f = throws.index("right_forehand")
            i_b = throws.index("right_backhand")
            gain_w = (float(np.nanmean(sr[:, :, iw, i_f, :])) -
                      float(np.nanmean(sr[:, :, iw, i_b, :])))
            gains = {int(headings[i]):
                     float(np.nanmean(sr[:, :, i, i_f, :])) -
                     float(np.nanmean(sr[:, :, i, i_b, :]))
                     for i in range(len(headings))}
            ok = gain_w >= max(gains.values()) - 1e-9
            print(verdict(ok, "forehand-minus-backhand success by heading: "
                              + ", ".join(f"{k}:{v:+.4f}"
                                          for k, v in gains.items())))
            results.append(dict(hypothesis="H2b_forehand_benefits",
                                verdict="SUPPORTED" if ok else "CONTRADICTED",
                                forehand_minus_backhand=gains))

    # ---------------- H3: best direction is the SSW prevailing ------------
    print("\n" + "-" * 78)
    print("H3  best-condition heading is a light SSW (~210 deg)")
    print("-" * 78)
    mean_by_h = {int(h): float(np.nanmean(Q[..., i]))
                 for i, h in enumerate(headings)}
    best = max(mean_by_h, key=mean_by_h.get)
    ok = abs(best - 210) <= 30
    print(verdict(ok, "mean Q by heading: "
                      + ", ".join(f"{k}:{v:.4f}" for k, v in mean_by_h.items())
                      + f"  -> best = {best} deg"))
    print("    NOTE: 'best' means best OF THE COMPUTED HEADINGS only.")
    results.append(dict(hypothesis="H3_best_is_SSW",
                        verdict="SUPPORTED" if ok else "CONTRADICTED",
                        mean_Q_by_heading=mean_by_h, best_heading=best,
                        caveat="best of the computed subset only"))

    with open(args.out_json, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
