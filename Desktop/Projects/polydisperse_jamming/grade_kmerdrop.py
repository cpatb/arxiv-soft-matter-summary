"""grade_kmerdrop.py -- PREREG_KMERDROP_2026-09-22.md. Writes out/grade_kmerdrop.txt.
Input: out/kmerdrop/<label>/ from `scripts/pull_jamwin.sh kmerdrop`. Written before the data.
"""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "kmerdrop")
LADDER = (("w0.00", 0.00, 0.54932, 1.0), ("w0.25", 0.25, 0.53258, 1.0),
          ("w0.50", 0.50, 0.51708, 0.9036), ("w0.75", 0.75, 0.50282, 0.7937),
          ("w1.00", 1.00, 0.48856, 0.7075))


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("grade_kmerdrop -- PREREG_KMERDROP_2026-09-22.md (x = 0.85, sigma* 300, delta 1.3169)")
    rows = []
    for lab, w, pm, jf in LADDER:
        d = os.path.join(ROOT, lab)
        if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            say(f"  {lab}: not landed")
            continue
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        m = re.search(r"phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+) sigma\*=([\d.e+]+) Pa "
                      r"delta=([\d.]+)", ban)
        jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
        prov = bool(m and jb) and all(math.isclose(float(v), t, rel_tol=1e-6) for v, t in
                                      zip(m.groups(), (0.527, 0.62, pm, 300.0, 1.3169))) \
            and math.isclose(float(jb.group(1)), jf, rel_tol=1e-6)
        a = np.loadtxt(os.path.join(d, "spread.dat"), comments="#", ndmin=2) \
            if os.path.exists(os.path.join(d, "spread.dat")) else np.zeros((0, 21))
        healthy = st["rc"] == 0 and a.size and a[-1, 1] >= 5.9
        r = {"lab": lab, "w": w, "prov": prov, "healthy": bool(healthy), "rc": st["rc"],
             "t_end": st["last_t_over_tau"]}
        if a.size:
            r.update(beta=float(a[:, 2].max()), J=float(a[:, 7].max()),
                     fmean=float(a[:, 6].max()),
                     urms_late=float(a[a[:, 1] >= 5.0, 19].mean()) if (a[:, 1] >= 5.0).any()
                     else math.nan)
        rows.append(r)
        say(f"  {lab}: prov {prov} healthy {bool(healthy)} (rc {st['rc']}, t/tau {st['last_t_over_tau']})"
            + (f"  beta_max {r['beta']:.3f}  J(f_c={jf}) {r['J']:.4f}  <fhat>max {r['fmean']:.3f}"
               f"  u_rms(t/tau>=5) {r['urms_late']:.4f}" if a.size else ""))
    ok = [r for r in rows if r["prov"] and r["healthy"]]
    say(f"D3 health: {len(ok)}/{len(rows)} (of 5)")
    if len(ok) == 5:
        b = [r["beta"] for r in ok]
        d1 = all(b[i] > b[i + 1] for i in range(4))
        say(f"D1 beta_max monotone falling in w: {[round(v, 3) for v in b]}  {'PASS' if d1 else 'FAIL'}")
        j = [r["J"] for r in ok[2:]]
        d2 = ok[0]["J"] == 0 and ok[1]["J"] == 0 and j[0] > 0 and j[0] <= j[1] <= j[2]
        say(f"D2 J at w 0.5/0.75/1 = {[round(v, 4) for v in j]} (w 0, 0.25: "
            f"{ok[0]['J']}, {ok[1]['J']})  {'PASS' if d2 else 'FAIL'}")
    else:
        say("D1/D2: not graded -- fewer than 5 healthy cells; rows above are descriptive")
    # ---- Amendment 1: all channels (chi + sheared sigma* + sheared delta) -----------------
    say("")
    say("=== Amendment 1: all doublet channels vs window-only ===")
    say("   registered: J_all / J_window in [1/1.5, 1.5]; beta_max within one grid quantum (1/64)")
    a1 = []
    for lab, ref in (("all_w0.50", "w0.50"), ("all_w1.00", "w1.00")):
        pa = os.path.join(HERE, "out", "kmerdrop2", lab, "spread.dat")
        pb = os.path.join(ROOT, ref, "spread.dat")
        sa = os.path.join(HERE, "out", "kmerdrop2", lab, "RUN_STATUS.json")
        if not (os.path.exists(pa) and os.path.exists(pb) and os.path.exists(sa)):
            say(f"  {lab}: not landed")
            continue
        st = json.load(open(sa, encoding="utf-8"))
        a = np.loadtxt(pa, comments="#", ndmin=2)
        b = np.loadtxt(pb, comments="#", ndmin=2)
        healthy = st["rc"] == 0 and a[-1, 1] >= 5.9
        ratio = a[:, 7].max() / b[:, 7].max()
        dbeta = abs(a[:, 2].max() - b[:, 2].max())
        ok_ = (1 / 1.5 <= ratio <= 1.5) and dbeta <= 1 / 64 + 1e-12
        # A death is EXCLUDED, never a FAIL (the D3 health rule of this prereg; corrected
        # 2026-09-23 after all_w1.00 died -- the first version appended False for a dead cell).
        if healthy:
            a1.append(ok_)
        say(f"  {lab}: healthy {healthy}  J {a[:, 7].max():.4f} vs {b[:, 7].max():.4f} (ratio "
            f"{ratio:.2f})  beta {a[:, 2].max():.4f} vs {b[:, 2].max():.4f}  Rg "
            f"{a[:, 11].max():.4f} vs {b[:, 11].max():.4f}  "
            + (('PASS' if ok_ else 'FAIL') if healthy else f"EXCLUDED (rc {st['rc']}, t/tau {a[-1, 1]})"))
        if not healthy:
            # POST-HOC, labelled: the same comparison truncated before the first discontinuity in
            # this cell's J series -- a jump > 0.2 between consecutive samples AFTER the spread has
            # peaked (argmax of Rg_r), when the drop is at rest and cannot produce one (the impact
            # rise itself is excluded by construction) -- applied to BOTH records at the same t.
            # Not a registered verdict.
            dj = np.diff(a[:, 7])
            i_rest = int(np.argmax(a[:, 11]))
            jump = np.where((dj > 0.2) & (np.arange(dj.size) >= i_rest))[0]
            if jump.size:
                tcut = a[jump[0] + 1, 1]
                ma, mb = a[:, 1] < tcut, b[:, 1] < tcut
                r2 = a[ma, 7].max() / b[mb, 7].max()
                db2 = abs(a[ma, 2].max() - b[mb, 2].max())
                say(f"    POST-HOC truncated at t/tau {tcut} (first J jump {dj[jump[0]]:+.3f}): J "
                    f"{a[ma, 7].max():.4f} vs {b[mb, 7].max():.4f} (ratio {r2:.3f}), beta "
                    f"{a[ma, 2].max():.4f} vs {b[mb, 2].max():.4f}: "
                    f"{'within the registered bars' if (1 / 1.5 <= r2 <= 1.5 and db2 <= 1 / 64 + 1e-12) else 'outside'}")
    say(f"Amendment 1: {sum(a1)}/{len(a1)} healthy cells graded" + ("  PASS" if len(a1) == 2 and all(a1) else
                                              "  FAIL" if a1 and not all(a1) else "  INCOMPLETE"))
    with open(os.path.join(HERE, "out", "grade_kmerdrop.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
