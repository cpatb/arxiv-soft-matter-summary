"""Print every headline number from every results file, in one compact block.

Exists so the written report can be assembled from the actual outputs rather
than from numbers copied by hand out of scattered JSON.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from catfish import frames  # noqa: E402

RESULTS = frames.ROOT / "results"


def load(name: str):
    p = RESULTS / name
    if not p.exists():
        return None
    with open(p) as fh:
        return json.load(fh)


def section(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def main() -> None:
    section("BODY LENGTH  (the unit everything is divided by)")
    bl = load("bodylength_summary.json")
    if bl:
        for k, v in bl.items():
            print(f"  {k}: BL = {v['body_px']:.2f} px   n = {v['n_silhouettes']}   "
                  f"IQR {v['iqr_px'][0]:.1f}-{v['iqr_px'][1]:.1f}   "
                  f"aspect {v['aspect_median']:.2f}   "
                  f"threshold systematic {v['systematic_pct']:.1f}%")

    section("PER-CLIP SUMMARY")
    p = RESULTS / "summary_by_clip.csv"
    if p.exists():
        df = pd.read_csv(p)
        pd.set_option("display.width", 250, "display.max_columns", 60)
        print(df.to_string(index=False))

    section("COLLECTIVE STATE  (polarisation / rotation / rind vs core)")
    cs = load("collective_state.json")
    if cs:
        for k, v in cs.items():
            print(f"  {k}: Phi = {v['polarisation_median']:.3f} "
                  f"[{v['polarisation_iqr'][0]:.2f}-{v['polarisation_iqr'][1]:.2f}]   "
                  f"M = {v['rotation_median']:.3f} "
                  f"[{v['rotation_iqr'][0]:.2f}-{v['rotation_iqr'][1]:.2f}]   "
                  f"-> {v['state'].upper()}")
            print(f"      Phi rind {v['pol_rind_median']:.3f} vs core "
                  f"{v['pol_core_median']:.3f};  speed core/rind "
                  f"{v['speed_core_over_rind']:.3f}")

    section("CORRELATION  (scale-free test)")
    an = load("analysis_summary.json")
    if an and an.get("correlation_scaling"):
        cc = an["correlation_scaling"]
        for k, v in cc.get("per_clip", {}).items():
            flag = "  <-- SATURATED" if v.get("saturated") else ""
            # xi_decay_bl, not xi_bl: C(r) never crosses zero in this dataset,
            # so what is reported everywhere is the 1/e decay length.
            print(f"  {k}: xi = {v['xi_decay_bl']:.2f} BL (1/e decay)   "
                  f"r_max = {v['r_max_bl']:.1f} BL   "
                  f"xi/r_max = {v['xi_over_rmax']:.2f}   D = {v['median_size_bl']:.1f} BL"
                  f"{flag}")
        if "loglog_slope" in cc:
            print(f"\n  size bins D (BL): {[round(x,1) for x in cc['bin_size_bl']]}")
            print(f"  xi per bin  (BL): {[round(x,2) for x in cc['bin_xi_bl']]}")
            print(f"  curves per bin  : {cc['n_curves_per_bin']}")
            print(f"  log-log slope   : {cc['loglog_slope']:.3f}   r = {cc['loglog_r']:.3f}")
            print(f"  xi/D median     : {cc['xi_over_L_median']:.3f}")
            print(f"  size range      : {cc['size_range_bl'][0]:.1f}-{cc['size_range_bl'][1]:.1f} BL "
                  f"(dynamic range {cc['size_dynamic_range']:.2f}x)")
            print("  NOTE: slope ~1 = xi proportional to D (scale-free); "
                  "slope ~0 = fixed interaction range.")

    section("CORRELATION NULL  (is xi real, or an artefact of mean subtraction?)")
    cn = load("correlation_null.json")
    if cn:
        for k, v in cn.items():
            print(f"  {k}: observed xi = {v['xi_bl_observed']:.2f} BL   "
                  f"shuffled xi = {v['xi_bl_shuffled']:.2f} BL   "
                  f"ratio = {v['ratio_observed_over_shuffled']:.1f}x   "
                  f"p = {v['permutation_p']:.3f}")

    section("PACKING PROFILE  (rind -> core)")
    if an and an.get("packing_profile"):
        for k, v in an["packing_profile"].items():
            print(f"  {k}: edge {v['edge_packing']:.3f} -> core {v['core_packing']:.3f}   "
                  f"(core - edge = {v['core_minus_edge']:+.3f}, n = {v['n_profiles']})")

    section("FISSION / FUSION")
    if an and an.get("events"):
        for k, v in an["events"].items():
            print(f"  {k}: {v['n_fission_frames']} fission, {v['n_fusion_frames']} fusion")

    section("FRINGE INDIVIDUALS  (independent speed check)")
    fa = load("fringe_analysis.json")
    if fa:
        for k, v in fa.items():
            print(f"  {k}: {v['n_tracks']} tracks   "
                  f"track speed {v['track_speed_bl_s_median']:.2f} BL/s   "
                  f"flow speed {v['flow_speed_bl_s_median']:.2f} BL/s   "
                  f"ratio {v['speed_ratio_track_over_flow']:.2f}")
            print(f"      fringe NND {v['nnd_bl_median']:.2f} BL "
                  f"(upper bound; detection incomplete)")

    section("CAMERA MOTION QA  (how much apparent order was the photographer?)")
    cq = load("camera_qa.json")
    if cq:
        print(f"  n frames: {cq['n_frames']}")
        print(f"  polarisation  uncorrected {cq['pol_raw_median']:.3f}  ->  "
              f"corrected {cq['pol_corrected_median']:.3f}")
        print(f"  median inflation {cq['median_inflation']:+.3f}, "
              f"max {cq['max_inflation']:+.3f}")
        print(f"  frames inflated by >0.2: {100*cq['frac_frames_inflated_by_over_0.2']:.1f}%")
        w = cq["worst_case"]
        print(f"  worst case: {w['clip']} frame {w['frame']}, camera "
              f"{w['cam_shift_px']:.1f} px/frame, Phi {w['pol_raw']:.3f} -> {w['pol_cor']:.3f}")

    section("RESOLUTION FLOORS  (what cannot be measured)")
    if bl:
        for k, v in bl.items():
            print(f"  {k}: optical-flow window 25 px = {25/v['body_px']:.2f} BL  "
                  f"-> xi below ~{25/v['body_px']:.1f} BL is the estimator, not the fish")
    print()


if __name__ == "__main__":
    main()
