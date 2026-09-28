"""so_watch.py -- PREREG_so_triplet Amendment 2, health rule H1, applied operationally.

simpleFoam's residualControl cannot fire on grids where the round-off Uz residual stalls
(see the amendment). Every POLL seconds this reads each running case under ~/bidet_so and,
at the first iteration where every initial residual EXCEPT Uz is below TARGET, sets
`stopAt writeNow;` in that case's system/controlDict (runTimeModifiable is on), and logs
the iteration and residuals to ~/bidet_so/lane.log. H2 (Uz is round-off) is checked by
the grader on the final field, not here.

Stops when ~/bidet_so/watch.stop exists or after MAX_HOURS. Run inside WSL.
"""
import glob
import os
import re
import time

TARGET = 1e-5
POLL = 60
MAX_HOURS = 9.0
ROOT = os.path.expanduser("~/bidet_so")
pat = re.compile(r"Solving for (\w+),\s+Initial residual = ([0-9.eE+-]+)")


def last_complete_iteration(log):
    """(iteration, {field: first initial residual}) for the last COMPLETE iteration."""
    it, cur, done = None, {}, None
    with open(log, errors="replace") as fh:
        for line in fh:
            if line.startswith("Time = "):
                if it is not None and cur:
                    done = (it, cur)
                it, cur = int(float(line.split("=")[1])), {}
            else:
                m = pat.search(line)
                if m and m.group(1) not in cur:
                    cur[m.group(1)] = float(m.group(2))
    return done


def main():
    t0 = time.time()
    stopped = set()
    while time.time() - t0 < MAX_HOURS * 3600 and not os.path.exists(os.path.join(ROOT, "watch.stop")):
        for case in sorted(glob.glob(os.path.join(ROOT, "case*"))):
            log = os.path.join(case, "log.simpleFoam")
            cd = os.path.join(case, "system", "controlDict")
            if case in stopped or not (os.path.exists(log) and os.path.exists(cd)):
                continue
            with open(cd) as fh:
                txt = fh.read()
            if "stopAt          writeNow;" in txt:
                stopped.add(case)
                continue
            r = last_complete_iteration(log)
            if r is None:
                continue
            it, res = r
            others = {k: v for k, v in res.items() if k != "Uz"}
            if others and max(others.values()) < TARGET:
                txt = re.sub(r"^stopAt\s+.*$", "stopAt          writeNow;", txt, flags=re.M)
                with open(cd + ".tmp", "w", newline="\n") as fh:
                    fh.write(txt)
                os.replace(cd + ".tmp", cd)
                stopped.add(case)
                msg = (f"[{os.path.basename(case)}] H1 met at iteration {it}: worst non-Uz "
                       f"{max(others.values()):.3e} ({max(others, key=others.get)}), "
                       f"Uz {res.get('Uz', float('nan')):.3e} -> stopAt writeNow "
                       f"{time.strftime('%F %T')}")
                print(msg, flush=True)
                with open(os.path.join(ROOT, "lane.log"), "a") as fh:
                    fh.write(msg + "\n")
        time.sleep(POLL)


if __name__ == "__main__":
    main()
