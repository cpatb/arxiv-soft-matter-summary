"""diag_carve_release.py — is the low-J spreading the IMPACT, or the carving?

Phase 5 varies the impact number J = rho u0^2 / p_carve, and it reaches low J by
raising the carve pressure. But a pack compressed to a high p_carve and then cut out of
its periodic cell carries stored elastic energy, and can fly apart on release with no
impact at all. Phase 4 hit exactly this and lowered CARVE_P_HAT from 1e-3 to 1e-4
because the free aggregate lost every contact during a settle.

If that is happening again, then `corr(log J, spread) = -0.946` for caged CV = 0.5 is
measuring the carving, not the collision, and the J axis is confounded with it.

The discriminator is already in the records: `z_initial` is read AS CARVED and
`z_launch` after a 300-step settle under the boundary condition, both BEFORE the
aggregate is given any velocity. A pack that is disintegrating on release loses
coordination between those two, with no impact involved.

Prints the drop for every cell, grouped by carve pressure.
"""

import collections
import glob
import json

import numpy as np


def main() -> None:
    recs = []
    for f in glob.glob("out/phase5*.jsonl") + glob.glob("out/phase4*.jsonl"):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if not r.get("failed") and r.get("z_initial") is not None:
                    recs.append(r)
    if not recs:
        print("no records")
        return

    g = collections.defaultdict(list)
    for r in recs:
        cp = r.get("carve_p_hat", float("nan"))
        g[(r["bc"], round(r["cv"], 3), cp)].append(r)

    print("Coordination BEFORE any impact: as carved -> after a 300-step settle.")
    print("A large drop means the aggregate is coming apart on RELEASE, not on impact.")
    print()
    print(f"{'bc':>10} {'CV':>5} {'p_carve':>9} {'n':>3} {'Z_carved':>9} "
          f"{'Z_launch':>9} {'drop':>8} {'drop %':>8} {'spread':>8}")
    for k in sorted(g, key=lambda x: (x[0], x[1], x[2])):
        bc, cv, cp = k
        rows = g[k]
        zi = np.array([r["z_initial"] for r in rows], float)
        zl = np.array([r.get("z_launch", np.nan) for r in rows], float)
        sp = np.array([r.get("spread", np.nan) for r in rows], float)
        ok = np.isfinite(zi) & np.isfinite(zl)
        if not ok.any():
            continue
        d = float(np.mean(zi[ok] - zl[ok]))
        pct = 100.0 * d / float(np.mean(zi[ok])) if np.mean(zi[ok]) else float("nan")
        print(f"{bc:>10} {cv:5.2f} {cp:9.1e} {ok.sum():3d} {np.mean(zi[ok]):9.3f} "
              f"{np.mean(zl[ok]):9.3f} {d:8.3f} {pct:7.1f}% "
              f"{np.nanmean(sp):8.3f}")

    print()
    print("Reading: if the drop grows with p_carve, the high-p_carve aggregates are")
    print("disintegrating on release and the J axis is confounded with the carving,")
    print("because low J is reached BY raising p_carve.")


if __name__ == "__main__":
    main()
