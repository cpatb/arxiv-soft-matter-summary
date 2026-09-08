#!/usr/bin/env python3
"""Rebuild ``cfd_fit_summary.json``'s ``grid_independence`` block from the CONVERGED
grid family, replacing values the 2026-09-03 grid study retracted.

Why this exists. The block shipped ``refined_x2.C_P = 1.5178`` and a 12.28 % change,
both from a solve that stopped at its 3000-iteration cap at 3e-5 against a 1e-5 target.
Restarted, that grid converges at iteration 5105 and ``C_P`` moves to 1.39279. The stale
value reached ``model/run_bidet.py``'s printed output and the manuscript's grid-independence
sentence. Nothing else in the summary is touched: the standoff and Reynolds sweeps, the
free-jet validation and ``all_cases`` are unchanged and are re-checked for identity.

The numbers come from ``convergence/gci_Re10000_H6_fine_conv.csv``, which the GCI tool
writes, NOT from this file. If the family is re-run, re-run this.

Run:  python refresh_grid_block.py [--check]
"""
from __future__ import annotations

import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUMMARY = os.path.join(HERE, "cfd_fit_summary.json")
GCI = os.path.join(HERE, "convergence", "gci_Re10000_H6_fine_conv.csv")

# phi1 / phi2 / phi3 in Celik's ordering are the FINEST first: refine 4 / 2 / 1.
REFINE_OF = {"phi1": 4.0, "phi2": 2.0, "phi3": 1.0}


def read_gci() -> dict[str, dict]:
    out = {}
    with open(GCI, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            out[row["quantity"]] = row
    return out


def build_block(g: dict[str, dict]) -> dict:
    def v(q: str, key: str) -> float:
        return float(g[q][key])

    base = {"n_t": v("n_t", "phi3"), "C_P": v("C_P", "phi3"),
            "tau15": v("tau_at_15mm_Pa", "phi3"), "tau_peak": v("tau_peak_Pa", "phi3")}
    x2 = {"n_t": v("n_t", "phi2"), "C_P": v("C_P", "phi2"),
          "tau15": v("tau_at_15mm_Pa", "phi2"), "tau_peak": v("tau_peak_Pa", "phi2")}
    x4 = {"n_t": v("n_t", "phi1"), "C_P": v("C_P", "phi1"),
          "tau15": v("tau_at_15mm_Pa", "phi1"), "tau_peak": v("tau_peak_Pa", "phi1")}

    def rel(a: dict, b: dict) -> dict:
        return {k: (a[k] - b[k]) / b[k] for k in b}

    return {
        "source": "convergence/gci_Re10000_H6_fine_conv.csv, family fine_conv "
                  "(refine 4/2/1, 491520/122880/30720 cells)",
        "iterations_to_residualControl": {"refine_4": 11552, "refine_2": 5105,
                                          "refine_1": 2261},
        "base_refine_1": base,
        "refined_x2": x2,
        "refined_x4": x4,
        "rel_change": rel(x2, base),
        "rel_change_x4": rel(x4, base),
        "verdict": {q: g[q]["verdict"] for q in g},
        "RETRACTED_2026-09-03": {
            "what": "refined_x2 previously read C_P = 1.5178 with rel_change 0.1228 "
                    "(12.28 %). That solve stopped at its 3000-iteration cap at 3e-5 "
                    "against a 1e-5 residualControl target; restarted it converges at "
                    "iteration 5105 and C_P is 1.39279, i.e. 3.03 %.",
            "old_C_P": 1.5178,
            "old_rel_change_C_P": 0.1228265858821802,
            "where_it_leaked": ["model/run_bidet.py printed line",
                                "bidet_optimization_1.JOURNAL.tex grid-independence "
                                "sentence (b)"],
        },
        "WARNING": "These per-grid percentages are NOT a convergence study. C_P DIVERGES "
                   "under refinement (|eps21/eps32| = 4.589), n_t and tau_at_15mm are "
                   "non-monotone, and the only quantity with an observed order "
                   "(tau_peak, p = 0.195) is outside the asymptotic band. Quote them as "
                   "sensitivities, never as an error estimate. See GRID_STUDY_2026-09-03.md.",
    }


def main() -> int:
    check = "--check" in sys.argv
    with open(SUMMARY, encoding="utf-8") as f:
        s = json.load(f)
    before = json.dumps({k: v for k, v in s.items() if k != "grid_independence"},
                        sort_keys=True)

    block = build_block(read_gci())
    if check:
        cur = s.get("grid_independence", {})
        same = json.dumps(cur, sort_keys=True) == json.dumps(block, sort_keys=True)
        print("grid_independence is up to date" if same else
              "grid_independence is STALE; run without --check")
        return 0 if same else 1

    s["grid_independence"] = block
    after = json.dumps({k: v for k, v in s.items() if k != "grid_independence"},
                       sort_keys=True)
    assert before == after, "this script must touch grid_independence and nothing else"
    with open(SUMMARY, "w", encoding="utf-8") as f:
        json.dump(s, f, indent=1)
        f.write("\n")

    r = block["rel_change"]
    print("rebuilt grid_independence from", os.path.basename(GCI))
    print("  x2 refinement (converged): n_t %+.2f %%, C_P %+.2f %%, tau15 %+.2f %%"
          % (r["n_t"] * 100, r["C_P"] * 100, r["tau15"] * 100))
    r4 = block["rel_change_x4"]
    print("  x4 refinement (converged): n_t %+.2f %%, C_P %+.2f %%, tau15 %+.2f %%"
          % (r4["n_t"] * 100, r4["C_P"] * 100, r4["tau15"] * 100))
    print("  RETRACTED: C_P = 1.5178 / 12.28 %% (unconverged refine = 2)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
