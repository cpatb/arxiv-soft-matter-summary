"""One-shot prose pass on README, 01, 02 (humanizer skill). Exact-match replacements;
asserts each old string occurs once. Kept for provenance."""
from pathlib import Path

ROOT = Path(r"C:\Users\compl\Desktop\Projects\biofluids_venture")

EDITS = {
"README.md": [
("# biofluids_venture — foundation for a bodily-fluid simulation and testing company",
 "# biofluids_venture: foundation for a bodily-fluid simulation and testing company"),
("The error this introduces is not small and, more importantly, it is\nalmost never quantified.",
 "The error this introduces is not small, and it is almost never quantified."),
("phases, gates and budgets under the 5–10 hr/week constraint",
 "phases, gates and budgets under the 5 to 10 hr/week constraint"),
],
"01_LANDSCAPE.md": [
("# 01 — Landscape: who does what, and where the gaps are",
 "# 01. Landscape: who does what, and where the gaps are"),
("Herschel–Bulkley triple", "Herschel-Bulkley triple"),
("the formulations at 100–300 mg/mL", "the formulations at 100 to 300 mg/mL"),
("This is a small, unglamorous and completely open modelling problem with a\nmedical and a consumer market attached.",
 "The modelling problem is open, nobody prominent is working on it, and it has a\nmedical and a consumer market attached."),
("## 4. The empty slots, stated plainly", "## 4. The empty slots"),
("Slots 1–4 are software.", "Slots 1 to 4 are software."),
],
"02_NICHE_AND_POSITIONING.md": [
("# 02 — Niche and positioning", "# 02. Niche and positioning"),
("The word carrying the weight is *credible*.", "The important word is *credible*."),
("Requires having been burned by it. The `yfrac-is-a-clamp-readout` episode (a month of a yield-stress observable measuring the regularisation ceiling) is the product requirement, written in blood.",
 "Requires having been burned by it. The `yfrac-is-a-clamp-readout` episode (a month of a yield-stress observable measuring the regularisation ceiling) is where this requirement came from."),
("A process, not a feature. Most consultancies cannot afford to tell a client the model cannot answer the question. Making that the brand is a moat precisely because it costs revenue in the short run.",
 "A process, not a feature. Most consultancies cannot afford to tell a client the model cannot answer the question. It is hard to copy because it costs revenue in the short run."),
("**W1 — Injectability and drug-product fluid handling (primary).**",
 "**W1. Injectability and drug-product fluid handling (primary).**"),
("Herschel–Bulkley fluids (Buckingham–Reiner)", "Herschel-Bulkley fluids (Buckingham-Reiner)"),
("**W2 — Airway, nasal and mucus (second).**", "**W2. Airway, nasal and mucus (second).**"),
("**W3 — Device-induced thrombosis and hemolysis (flagship research asset, not first revenue).**",
 "**W3. Device-induced thrombosis and hemolysis (flagship research asset, not first revenue).**"),
("The `thrombosis_jamming` closure is genuinely novel: occlusion as",
 "The `thrombosis_jamming` closure is novel as far as the prior-art check has gone (an arXiv\nprecedent is still flagged for verification in `04_ASSET_MAP.md`): occlusion as"),
("**W4 — Urology, continence, surgical splash and aerosol (public case studies).**",
 "**W4. Urology, continence, surgical splash and aerosol (public case studies).**"),
("Uroflow splashback, bidet jets, drop impact on wet surfaces, catheter flow. Directly finished\nor nearly so.",
 "Uroflow splashback, bidet jets, drop impact on wet surfaces, catheter flow. Finished or\nnearly finished."),
("**W5 — Ophthalmic, joint, GI, dysphagia (later, same stack).**",
 "**W5. Ophthalmic, joint, GI, dysphagia (later, same stack).**"),
("highest value per study and the strongest moat, but requires",
 "highest value per study and the hardest to copy, but requires"),
("""- Not a general CFD shop. Newtonian device flows are served, and competing there is a price war.
- Not a preclinical image-analysis service. Declined on record (`monetization-direction`).
- Not a "rigour for AI-generated science" brand. The rigour is a feature of the named work, not
  the identity.
- Not a solver vendor. Solvers are open (OpenFOAM, Basilisk) or the client's; the product is
  the material layer, the bridge, the validation, and the credibility trail.""",
"""- It does not compete as a general CFD shop. Newtonian device flows are already served, and
  competing there is a price war.
- It does not do preclinical image analysis for hire. That was declined on record
  (`monetization-direction`).
- It is not branded as "rigour for AI-generated science". The rigour is a feature of the named
  work rather than the identity.
- It does not sell a solver. Solvers are open (OpenFOAM, Basilisk) or the client's; the product
  is the material layer, the bridge, the validation, and the credibility trail."""),
],
}

for fname, edits in EDITS.items():
    p = ROOT / fname
    t = p.read_text(encoding="utf-8")
    for old, new in edits:
        n = t.count(old)
        assert n == 1, f"{fname}: {n} matches for {old[:50]!r}"
        t = t.replace(old, new)
    p.write_text(t, encoding="utf-8", newline="\n")
    rem = t.count("\u2014") + t.count("\u2013")
    print(f"{fname}: {len(edits)} edits, {rem} dashes remaining")
