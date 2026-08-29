# Candidate papers and datasets for gate G4

> **Superseded in part, 2026-08-27.** Datasets 1 to 3 below are now downloaded and
> surveyed; see `DATA_INVENTORY.md` for what is actually in them. One correction from that
> survey: **Guy, Hermes & Poon is PMMA, not silica.** The system is PHSA-stabilised PMMA,
> d = 268 to 3770 nm, in a solvent of 2.4 mPa.s. Silica appears only as literature
> comparison points. The genuine silica in the download is Lin et al., Fig 2a, 2 um silica,
> and that is a single volume fraction rather than a series.

G4 wants two or more published silica DST papers, with `phi_m` recovered by `dst-fit`
landing within 0.02 of the paper's own stated value. That needs three things at once, and
most papers give only the first:

1. flow curves for a series of volume fractions of one material,
2. a sweep that **brackets the frictional transition**, not just its onset,
3. the paper's own `phi_m`, stated numerically, to compare against.

Requirement 2 comes out of gate G5. An unconstrained `sigma*` does not run to a bound, it
migrates into whatever stress window was measured, so a published curve that stops partway
up the thickening branch will return a confident-looking `sigma*` that means nothing.
Check the frictional fraction runs roughly 0 to 1 across the plotted range before using a
paper as a reference.

## Verified: real numerical data deposited

The Poon group at Edinburgh deposits the numbers behind their figures in Edinburgh
DataShare. That is better than digitizing: it is the plotted data, not pixels.

The whole collection is worth browsing: https://datashare.ed.ac.uk/handle/10283/833

| # | Paper | Data | Contents |
|---|---|---|---|
| 1 | **Guy, Hermes & Poon**, *Towards a Unified Description of the Rheology of Hard-Particle Suspensions*, PRL **115**, 088304 (2015). arXiv:1506.09041 | [10.7488/ds/286](https://doi.org/10.7488/ds/286) → [handle 10283/834](https://datashare.ed.ac.uk/handle/10283/834) | `Fig 1.zip` 52.4 KB, `Fig 2.zip` 2.25 KB. Hard spheres from about 1 to 50 um, spanning colloidal to granular |
| 2 | **Lin, Guy, Hermes, Ness, Sun, Poon & Cohen**, *Hydrodynamic and Contact Contributions to Continuous Shear Thickening*, PRL **115**, 228304 (2015). arXiv:1509.02750 | [10.7488/ds/324](https://doi.org/10.7488/ds/324) → [handle 10283/886](https://datashare.ed.ac.uk/handle/10283/886) | Shear-reversal experiments on micron silica and latex, separating hydrodynamic from contact stress |
| 3 | **Richards, Royer, Liebchen, Guy & Poon**, *Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions*, PRL **123**, 038004 (2019) | [handle 10283/3357](https://datashare.ed.ac.uk/handle/10283/3357) | `Figure1.zip` 6.29 MB, `Figure2.zip` 1.52 MB, `Figure3.zip` 27.74 MB, `Figure4.zip` 199 KB |
| 4 | **Sun & Ness**, *Shear thickening regimes of dense non-Brownian suspensions* | Edinburgh DataShare, 2015 | Simulation rather than experiment. Useful as a second opinion on `phi_m`, not as an experimental reference |
| 5 | **Blanco et al.**, conching chocolate, 2019 | [10.7488/ds/2515](https://doi.org/10.7488/ds/2515) | Frictionally jammed solid to flowable suspension at maximal solid content. Not silica, but a hard-particle system at high `phi` |

**Start with #1.** It is the canonical experimental test of the frictional-contact picture
and it is the paper whose parameters the model was written against.

**#3 is already in your memory** as the `WC_CLOCK` source (Richards PRL 123 038004, the
measured rest clock). Worth pulling for a second reason: the oscillations give access to a
shear-jamming portion of the flow curve that conventional rheometry forbids, so those
curves may extend further up the thickening branch than an ordinary sweep, which is
exactly the coverage problem G5 flagged.

## Worth checking, data availability not yet confirmed

- **Royer, Blair & Hudson**, *Rheological Signature of Frictional Interactions in Shear
  Thickening Suspensions*, PRL **116**, 188301 (2016). arXiv:1603.00422. Silica spheres over
  a wide range of volume fractions, and it reports **first normal stress difference N1**
  alongside viscosity. Your corn-starch export has a normal-stress column too, so this is
  the one that could support a second observable rather than just a second material. NIST
  is a co-affiliation, so check the NIST public data repository.
- **Comtet, Chatté, Niguès, Bocquet, Siria & Colin**, *Pairwise frictional profile between
  particles determines discontinuous shear thickening transition*, Nat. Commun. **8**, 15633
  (2017). Nature journals require a data availability statement.
- **Hsu, Ramakrishna, Zanini, Spencer & Isa**, roughness-dependent shear thickening, PNAS
  (2018). Silica, and roughness is the variable that sets `phi_m`.

## Directly relevant to why N4 failed

**Clavaud, Singh et al.**, *Shear-thickening in presence of adhesive contact forces: the
singularity of cornstarch*, arXiv:2303.15915. The title is the finding. If corn starch
thickens through adhesive contacts rather than purely frictional ones, that is a
mechanism for `phi_m` not behaving as a constant, and it is a fifth candidate to add to
the four already in the spec. Read this before concluding anything about the 2023 series.

Also already on disk and relevant: `Desktop\Papers\LITERATURE_REVIEW\paper_synopses\`
contains synopses of `wyart-cates_2014_frictional-DST-model`,
`brown-jaeger_2014_rpp_shear-thickening-review`,
`reversible-shear-thickening_monodisperse_bender-wagner-1996` and
`quick-starch-guide_clavaud-singh_2025_ST-perspective`. The Bender and Wagner paper is
monodisperse silica and predates the frictional picture, so its flow curves are a clean
test set uncontaminated by having been fitted with the model already.

## If a paper has no deposit

Digitize with WebPlotDigitizer. Log-log axes need care: calibrate on two decades, not two
adjacent gridlines, and check a known point round-trips. Record the digitized points as a
CSV in `data/literature/` with the figure number and the axis calibration used, so the
comparison is reproducible rather than a one-off.
