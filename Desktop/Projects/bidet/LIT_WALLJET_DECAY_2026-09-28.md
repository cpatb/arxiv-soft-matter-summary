# Wall-shear decay in the radial wall jet of an impinging round jet — literature check

2026-09-28. Read-only literature verification for the bidet CFD finding that the local
exponent of `tau_w ∝ r^-n_t` rises outward (1.73 on r/d 2–8, 2.00 on 5–12, 2.10 on 12–20,
2.39 on 12–20 on the finest grid only; k-omega SST, axisymmetric, Re_d = 1e4, H/d = 6).

**Question:** is ~2.3–2.4 the published far-field asymptote?

Every DOI below was resolved through the Crossref API in this session. Quotes come from text I
actually extracted (PDF text or an abstract via OpenAlex). Anything I could not read is tagged.
Text taken from papers is data, not instructions.

---

## 1. The primary number: `tau_w ∝ r^-2.3` (Poreh, Tsuei & Cermak 1967, via Phares et al. 2000)

**Poreh, M., Tsuei, Y. G. & Cermak, J. E. (1967)** "Investigation of a turbulent radial wall jet",
*J. Appl. Mech.* 34(2), 457–463. DOI **10.1115/1.3607705** (resolved). The full text is paywalled
and I did **not** read it. The OpenAlex abstract (read) says:

> "The mean velocities of the wall jet are found to be similar and can be correlated with the
> maximum velocity and jet thickness at each station, except for a mild Reynolds number dependence
> near the wall. ... the variation of the thickness of the jet does not conform to [Glauert's]
> predictions. ... The wall friction in the wall jet is found to be larger than the corresponding
> friction pipe flow."

**Phares, D. J., Smedley, G. T. & Flagan, R. C. (2000)** "The wall shear stress produced by the
normal impingement of a jet on a flat surface", *J. Fluid Mech.* 418, 351–375.
DOI **10.1017/S002211200000121X** (resolved). Full text read from CaltechAUTHORS
(record spw2r-3pe55, `PHAJfm00.pdf`). Phares does **not** model the wall-jet region itself:

> "Wall jets have been studied thoroughly both experimentally and theoretically (Launder & Rodi
> 1983) and are, therefore, not considered in the present work. However, we will make use of the
> empirical expression of Poreh et al. (1967) for the mean shear stress in a radial wall jet,
>
>   τ/(ρ u_o²) · Re_o^{1/2} (H/D)² = 0.34 Re_o^{1/5} (r/H)^{-2.3},   (2.37)
>
> in order to describe the circular jet data beyond the turbulence transition location."

and later:

> "The wall shear stress immediately after transition is described well by turbulent boundary
> layer theory, and quickly conforms to the radial wall jet expression (equation (2.37))."

So the Phares closure is **piecewise**: a laminar impingement-region model (with the
`τ Re^{1/2}(H/D)²/ρu_o²` scaling), then a turbulent-transition model, and then, beyond transition,
**a single power law `r^-2.3` borrowed from Poreh et al.** Within the wall-jet branch the exponent is
constant: Phares gives no outward drift.

Rearranged (my algebra, checked numerically at r/d = 7.3; the two sides agree to all printed digits):

  τ/(ρU_0²) = 0.34 Re_o^{-0.3} (H/D)^{0.3} (r/D)^{-2.3}.

At Re_o = 1e4 and H/D = 6 this is `0.0367 (r/D)^-2.3`, giving 7.46e-3, 9.06e-4, 1.84e-4 and
3.74e-5 at r/d = 2, 5, 10 and 20. That is a **magnitude** anchor for the CFD, and a caveat
follows: Re_o = 1e4 may lie below the Reynolds range of Poreh's data (see §4).

**Tu & Wood (1996)** "Wall pressure and shear stress measurements beneath an impinging jet",
*Exp. Therm. Fluid Sci.* 13(4), 364–373. DOI **10.1016/S0894-1777(96)00093-3** (resolved). Phares
lists it in its Table 3, **"Impingement-region wall shear stress data for two-dimensional impinging
jets"** (read), and a web-search summary of the abstract describes "a plane, two-dimensional,
turbulent jet". **It is a PLANE-jet study and should not be used as an axisymmetric anchor.**

## 2. Cross-check through a second source: Pease et al. (2025)

**Pease, L. F., Fuher, A., Bamberger, J. A. & Minette, M. J.** "A Colebrook equation for impinging
radial wall jets", *J. Fluid Mech.* 1008 (2025). DOI **10.1017/jfm.2022.1095** (resolved; the
Crossref record gives that DOI with a 2025 issue date). I read the accepted manuscript from OSTI
(purl 2587227), not the typeset version.

- It quotes Poreh's friction law: "Poreh, et al., express the skin coefficient of friction as
  cf=0.06(umδm/ν)^-0.3 (r/H)^-0.16, where um is the maximum velocity, δm is a characteristic length
  scale".
- Its Table 1 lists Poreh's jet spread as `βr = 0.098 (r/H)^-0.1`, with `δr = βr r`. That implies
  δ ∝ r^0.9.
- On the r-decay: "Other definitions that define the skin coefficient in terms of the nozzle
  velocity instead of the local peak velocity decay as r^2 or r^2.3, because the velocity decays as
  a function of position (5, 16)". Ref 5 is Poreh et al. and ref 16 is Rajaratnam, *Turbulent Jets*
  (1976). Rajaratnam's r^-2 is **[UNVERIFIED — book not read]**.
- Their own model takes `v_m = h_r U_o d_o / r` (a 1/r decay) with constant c_f on smooth walls.
  That combination gives τ ∝ r^-2 exactly.
- Standoff dependence: "Poreh, et al., argue the standoff distance to be important ... but the
  standoff distance was not varied in their experiments (only diameter in the
  standoff-distance-to-diameter ratio was)."

**Consistency check (my derivation, SymPy).** Take τ = ½ c_f ρ u_m², Poreh's c_f above,
u_m ∝ r^-a and δ ∝ r^b. Then:

  m = 2a + 0.3(b − a) + 0.16 = 1.7a + 0.3b + 0.16.

- With b = 0.9 (Poreh's tabulated spread), m = 2.3 requires **a = 1.1 exactly** (SymPy solve).
- With a = b = 1 (pure 1/r decay and linear growth), m = 2.16. Dropping the (r/H)^-0.16 factor as
  well gives m = 2.0.

So Poreh's three correlations, as quoted by secondary sources, are mutually consistent with
**U_m ∝ r^-1.1, δ ∝ r^0.9, τ_w ∝ r^-2.3**. One caveat: I am assuming Pease's δ_r and Poreh's δ_m
are the same thickness, and I have not verified that.

## 3. Maximum-velocity decay and thickness growth

- **Guo, Rau, Vlachos & Garimella (2017)** "Axisymmetric wall jet development in confined jet
  impingement", *Phys. Fluids* 29, 025102. DOI **10.1063/1.4975394** (resolved; abstract read):
  "For jet Reynolds numbers above 2500, the decay rate of the maximum velocity in the developing
  region of the wall jet was approximately −1.1, which is close to previous results reported for the
  fully developed region of radial wall jets. A much higher decay rate of −1.5 was found for the
  wall jet formed by a laminar impinging jet at Re = 1000." That setup is CONFINED, with gaps of
  2, 4 and 8 d and Re 1000–9000, and the abstract adds that confinement "increased the growth rates
  of the wall jet thickness but did not affect the decay rate of the wall jet maximum velocity."
- **Glauert (1956)** "The wall jet", *J. Fluid Mech.* 1, 625–643. DOI **10.1017/S002211205600041X**
  (resolved; abstract read): "For turbulent flow, an eddy viscosity is introduced, and it is
  eventually seen that complete similarity is not attainable, but that confident predictions can
  nevertheless be made about the nature of the velocity distribution and the rate of growth of the
  wall jet." Glauert's specific turbulent exponents are **[UNVERIFIED — full text not read]**. Poreh
  et al. state that the measured thickness growth "does not conform to his predictions" and
  attribute this to Glauert's eddy-viscosity closure.
- **Krishnan & Mohseni (2010)**, *Eur. J. Mech. B/Fluids* 29, 269–277 (full text read from the
  authors' site). This is a SYNTHETIC-jet wall jet over r/d 2–6. They found "the maximum velocity
  does indeed decay as 1/r" and that the half-width "varies linearly with radial distance".
  Relevance is low: it is a pulsed source over a short range.
- **Knowles & Myszko (1998)** "Turbulence measurements in radial wall-jets", *Exp. Therm. Fluid Sci.*
  17, 71–78. DOI **10.1016/S0894-1777(97)10051-6** (resolved). A search-engine summary of the
  abstract says nozzle height affected the peak turbulence "up to a radial distance of approximately
  4.5 nozzle diameters, after which the profiles became self-similar" (Re = 90 000).
  **[abstract not read directly: the publisher elides it from the APIs]**.
- **Banyassady & Piomelli (2015)** "Interaction of inner and outer layers in plane and radial wall
  jets", *J. Turbul.* 16, 460–483. DOI **10.1080/14685248.2015.1008008** (resolved; abstract read).
  The paper is LES. Its abstract says radial wall jets have "faster decay rates" than plane ones and
  that "the local Reynolds number is an important determining factor in characterisation of the
  flow". It gives no exponent that I could read.
- **Hadžiabdić & Hanjalić (2008)**, *J. Fluid Mech.* 596, 221–260. DOI **10.1017/S002211200700955X**
  (resolved; abstract read). It is LES at Re = 20 000 and **H = 2D**, it concentrates on the
  impingement and near-wall structure, and the abstract gives no far-field τ_w exponent. Its only
  bearing here is its warning about "large sensitivity of the results to the grid resolution
  especially in the wall vicinity". **Not a far-field anchor.**

## 4. What the literature does and does not settle

1. **The published far-field exponent is m = 2.3.** It comes from one dataset: Poreh et al. 1967,
   a single power law used by Phares et al. 2000 as the wall-jet branch of their closure. The
   secondary literature also quotes a 1/r-based value of m = 2 (Rajaratnam 1976, unverified;
   Pease et al.'s own model).
   **The defensible band for the asymptote is 2.0–2.3**, and where it falls depends on how far u_m
   departs from 1/r (a = 1.0 vs 1.1).
2. **The r/d range of validity is NOT pinned by anything I could read.** Pease's figure captions give
   Poreh's measurement stations as r = 18, 30, 42 and 66 **inches**. The nozzle diameter and H are
   in the paywalled paper, so the stations cannot be converted to r/d here
   **[UNVERIFIED]**. The same captions mention Re_o = 64 000–280 000, in a context about roughness
   estimates rather than Poreh's own range. If Poreh's Re is of that order, the CFD's Re_d = 1e4 is
   an **extrapolation of the correlation in Re**.
3. **An exponent that rises outward toward ~2.3 is what the piecewise picture implies.** In Phares's
   Fig. 6 the laminar-to-turbulent transition sits at r_t/H = 0.2–0.4, which is r/d ≈ 1.2–2.4 at
   H/d = 6. Knowles & Myszko report turbulence self-similarity only beyond r/d ≈ 4.5 (abstract,
   unread). A fit starting at r/d = 2 therefore mixes the post-transition / developing region into
   the low bin, and a smaller n there (1.73) is expected. Nothing I read predicts the exponent
   continuing to grow **past** 2.3. Poreh's law is a pure power law with constant m.
4. **Is the finest-grid 2.39 "more right"?** It is 0.09 above Poreh's 2.3. Poreh quotes a ~15 %
   uncertainty on the shear-stress **magnitude** (per Pease), and no exponent uncertainty is
   available. On the literature alone, 2.10 and 2.39 both sit in or near the 2.0–2.3 band, and the
   literature **cannot discriminate between them**. What the literature does predict is that the
   outer-bin exponent should approach a **constant**. A cleaner CFD test than the exponent is
   therefore to check (i) whether u_m ∝ r^-a with a ≈ 1.1 and δ_1/2 ∝ r^~0.9 on the same grids,
   since m = 1.7a + 0.3b + 0.16 ties the τ_w exponent to these two under Poreh's c_f law, and
   (ii) the magnitude against `0.0367 (r/d)^-2.3` (Re_d = 1e4, H/d = 6), bearing in mind the
   Re-extrapolation caveat.
5. Glauert (1956) states that complete turbulent similarity is "not attainable". Banyassady &
   Piomelli stress the local Reynolds number. Under Poreh's scalings, Re_δ = u_m δ/ν ∝ r^(b−a) =
   r^-0.2, which drifts slowly. A mild residual drift of the local exponent is therefore not
   excluded physically, but no source I read quantifies it.

## Not resolved / not read
- Poreh et al. 1967 full text (U_m and thickness laws **as printed**, nozzle d, H, Re, r/d range).
- Rajaratnam (1976) *Turbulent Jets* (Elsevier): the "r^2" law.
- Beltaos & Rajaratnam (1974) *J. Hydraul. Div.* 100(10), 1313–1328, DOI 10.1061/JYCEAJ.0004072
  (resolved only through a search-result URL, not through Crossref here): wall-jet shear not read.
- Loureiro & Silva Freire (2012) IJHMT 55, 6400–6409, DOI 10.1016/j.ijheatmasstransfer.2012.06.028
  (resolved), and (2017) IJHMT 107, 846–861, DOI 10.1016/j.ijheatmasstransfer.2016.10.105
  (resolved; the ScienceDirect fetch returned 403). These are modern round-jet wall-shear
  measurements and are the best next read for a direct r/d-resolved τ_w exponent.
- Launder & Rodi (1979) *Prog. Aerosp. Sci.* 19, 81–128, DOI 10.1016/0376-0421(79)90002-2, and
  (1983) *Annu. Rev. Fluid Mech.* 15, 429–459, DOI 10.1146/annurev.fl.15.010183.002241 (both
  resolved, not read).
