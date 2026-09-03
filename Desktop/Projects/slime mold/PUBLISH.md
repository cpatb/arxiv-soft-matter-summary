# Publishing the repository — the last posting blocker, reduced to one decision

**2026-09-02.** `check_ready.py` had two failures. One is now closed; this document exists
so the other costs the author a decision and about five minutes, rather than an evening.

```
BEFORE (2026-08-31)          AFTER (2026-09-02)
FAIL p.1 affiliation         ✅ closed — see §1
FAIL p.8 repository URL      ⬜ open — needs the author, see §2
```

## 1. The affiliation was not a decision, it was a lookup

The flag in `doc/chicago_slime_mold.REVIEWED.tex` said the department was "left off
deliberately: this tree gives three different ones for the same author", and the Status
card carried it as an author decision. **The author settled it one day after the flag was
written and it was never propagated here.** `bidet_optimization_1.JOURNAL.tex` and
`.REVIEWED.tex` both carry, at line 72:

> AUTHOR-CONFIRMED 2026-07-23: Department of Physics and Astronomy is the canonical
> affiliation (chosen over the Mechanical Eng. / Chemical & Biomolecular Eng. variants that
> appear in other files). Email cpb@u.northwestern.edu confirmed for publication.

Corroborated independently by the drop-impact collaboration manuscript
(`Github/Drop Impact/lit_rheology_impact_2026-08-25/pdfs/src_2601/manuscript.tex`), which
gives *Department of Physics and Astronomy, Northwestern University, Evanston, IL 60208,
USA* for this author. Applied verbatim; the manuscript rebuilds at 10 pages and
`check_ready.py` no longer flags page 1.

⚠ If that 07-23 confirmation is wrong, this is one line to change and it is marked in the
`.tex` with its provenance.

## 2. What is left: the repository has to exist

`gh repo list` on this account has no slime-mold repository. The data-availability
statement on page 8 promises one, so the placeholder cannot be removed by editing prose —
**something has to be published.** That is an outward-facing, hard-to-undo action and it is
the author's call, not an agent's.

**The decision is only: public GitHub repo, Zenodo DOI, or both.** Everything else is done.

### The payload is ready and it is small

```
code/       116K   common_terminal_md.py  data_chicago.py  physarum.py  run_experiments.py
data/        37K   _headline_D.npy  _metrics.json  results.csv  results.json  (+2 stale, see below)
figures/    1.2M   fig1 … fig6 .png
doc/               chicago_slime_mold.REVIEWED.tex/.pdf
requirements.txt   pinned, and honest about what it cannot promise
README.md
bike/       1.2M   the bike-gap side study — DECIDE whether this ships (see below)
```

### It reproduces byte-identically, verified in this session

Not asserted — run. The local environment matches `requirements.txt` exactly (networkx
3.6.1, numpy 2.4.2, scipy 1.17.1, matplotlib 3.10.8, pandas 2.3.3, shapely 2.1.2), and both
stages were re-run from a clean tree:

```
python run_experiments.py 1      # rewrote _headline_D.npy _metrics.json results.csv fig1–fig4
python run_experiments.py 2      # rewrote results.json fig5 fig6
git status -- data/ figures/     # (no output)
```

Every regenerated file is byte-identical to the committed one — confirmed by file mtimes
moving while `git status` stays silent, so this is a real re-run and not a no-op. The
Steiner row comes back at **TL = 147.0 km, MD = 26.9 km**, exactly the networkx-version
artifact `requirements.txt` documents and the manuscript footnotes, rather than the
tabulated 145.2. **That is the file working as intended, not a discrepancy.**

Multi-seed headline reproduces at TL = 247.5 ± 11.1, TC = 175.7 ± 6.4, MD = 18.69 ± 0.34,
FT = 0.884 ± 0.019.

### Two small things to fix on the way out

1. **`data/gamma_sweep.json` and `data/multiseed.json` are stale.** The current
   `run_experiments.py` does not write them — they date from an earlier version and survive
   only because they are committed. Everything they held now lives in `results.json`.
   Delete them, or a reader will diff them against a run that never touches them.
2. **Decide whether `bike/` ships.** It is a separate study (`bike/REVIEW_bike.md`) with its
   own code, data and figures, and the manuscript does not cite it. Shipping it is fine;
   shipping it *unannounced* invites a reader to think it is part of the paper.

### Then, in order

```bash
# 1. author decision: name and visibility
gh repo create cpatb/chicago-slime-mold --public --description "Slime-mold physics model for Chicago 'L' expansion"

# 2. the payload (from a clean copy, not this working tree)
#    include: code/ data/ figures/ doc/ requirements.txt README.md
#    exclude: __pycache__/ run.log *.aux *.fls *.fdb_latexmk *.out *.log

# 3. add a LICENSE — no licence means "all rights reserved", which contradicts
#    "all code, data, and figure scripts ... are available"
#    MIT for the code, CC BY 4.0 for the data and figures, is the usual pairing

# 4. tag it, so the paper cites a fixed state rather than a moving branch
git tag -a v1.0-preprint -m "State of the code at preprint submission"
git push --tags

# 5. OPTIONAL but better: link the repo to Zenodo and cut a release, which mints a DOI.
#    The paper says "URL/DOI"; a DOI is the half that survives a repo being renamed.

# 6. fill the placeholder in doc/chicago_slime_mold.REVIEWED.tex (one line, ~line 445),
#    rebuild, and confirm
python check_ready.py        # must print PASS
```

## 3. What check_ready.py will still not tell you

It reads the compiled PDF for placeholders and for the affiliation. It does **not** check
that the URL you paste resolves, that the repository is public, or that the tag exists.
After step 6, open the link from a logged-out browser.
