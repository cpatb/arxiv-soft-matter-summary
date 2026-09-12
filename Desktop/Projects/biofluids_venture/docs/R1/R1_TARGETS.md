# R1 target sheet — the ten calls, named

Date: 2026-09-12. Purpose: remove the activation energy from `00_SYNTHESIS.md` §7, which has
been the binding gate since 2026-09-07 with no target list, no draft and no log on disk.

**Sourcing rule.** Every name below is traceable to a URL already in
`docs/market_research_2026-09-06/` or to a lookup run on 2026-09-12 and recorded here.
**No contact detail is invented.** Each row's `contact` column says exactly what is known:
*published* (printed in the paper), *snippet* (from a search result, needs one confirming
click), or *unresolved* (a form or a LinkedIn search the founder must do).

The segment quotas come from §7: 4 formulation, 3 device/combination-product, 2
regulatory/CM&S, 1 platform.

---

## Tier A — the two groups who published the failure (send these first)

These are not cold sales calls. **Both groups' own papers conclude that the extrapolation
this venture sells fails.** §7 question 1 ("when your injection-force model disagreed with
the bench, what did you do, and who did it?") is literally a question about their own
result, asked by a physicist who works on the same fluids. F9's documented failure was a
year of cold *sales* calls (CFD Support); a research question to a corresponding author is
a different act with a different base rate.

| # | who | why them | source | contact |
|---|---|---|---|---|
| A1 | **Georgina B. Armstrong** — Drug Substance Development, GlaxoSmithKline, Stevenage; also Strathclyde Institute of Pharmacy & Biomedical Sciences | First author and co-corresponding of the paper whose conclusion is *"Our findings caution against extrapolation from low concentration measurements, as they lack predictive ability for ultra-high concentration regimes"* (% differences **−186 % to +177 %**). Industry side of the wedge. | *MAbs* 16(1):2438172 (2024), **10.1080/19420862.2024.2438172**; PMC11790245 | **published**: Georgina.x.armstrong@gsk.com |
| A2 | **Prof. Zahra Rattray** — Strathclyde Institute of Pharmacy & Biomedical Sciences, University of Glasgow/Strathclyde | Senior/co-corresponding author of the same paper; academic side, so free to discuss method rather than product. Also the natural route to a co-authored short course or a joint measurement. | same | **published**: Zahra.rattray@strath.ac.uk |
| A3 | **Manasi Nawathe** — Johns Hopkins University (MS thesis, 20 Jul 2023) | *"Factors Impacting Injection Forces of High Concentration Monoclonal Antibody Formulations: The Significance of Needle Dimensions and Fluid Viscosity"* — the source of F4's central number: all three models (modified Hagen–Poiseuille + two non-Newtonian) hold to <25 % error but **deviate above 230 mg/mL and in 27G/30G needles**, attributed to "insufficient characterization of high shear viscosity behavior". That sentence is the venture's thesis, written by someone else. | jscholarship.library.jhu.edu item `32605470-29f4-42c5-a309-a530919d5d64` | **snippet** — author name from a 2026-09-12 search result; the repository page returns HTTP 403 to an automated fetcher. **Confirm the name and find current affiliation with one browser click before sending.** |

⚠ A3's name has not been read off the primary page by this session. Do not put it in an
email until the page is opened by hand.

## Tier B — formulation / CDMO (fills the §7 quota of 4 with A1)

| # | who | why them | source | contact |
|---|---|---|---|---|
| B1 | **Flora Felsovalyi — Lonza** | Moderated the PDA UPFS Vienna (21–22 Oct 2025) session on "the complex interplay between formulation science, delivery mechanics, and patient experience in the context of high-viscosity, large-volume subcutaneous administration". Someone who convenes that panel knows who has paid for what — §7 question 8 (referral) is the real ask here. | report 01 line 73; pda.org event agenda | **unresolved** (LinkedIn / Lonza) |
| B2 | **Simon Moss — Pfizer** | Spoke on "Exploring Opportunities for Large Volume Subcutaneous (LVSC) Injections for Acute Illness Treatment" at the same session. Buyer-side LVSC programme. | report 01 line 73 | **unresolved** |
| B3 | **Coriolis Pharma (CRO, Germany) — injectability & syringeability group** | The **only** seller found that states the problem in the venture's own words: *"biopharmaceutical products rarely show Newtonian flow behavior, their rheological properties (governing the injection process) depend on the force applied"*, with "In Silico Modeling" listed as a service and **no detail and no price**. Either the closest competitor or the most natural subcontract client; the call decides which. | report 01 line 37; coriolis-pharma.com | **unresolved** (web form) |

## Tier C — device / combination product (§7 quota of 3)

| # | who | why them | source | contact |
|---|---|---|---|---|
| C1 | **Barichello — Stevanato Group** and **Fetz — SHL Medical** | Co-authors of the Jan-2026 joint high-viscosity study: 1–50 cP, 2.0–2.25 mL, injection time ">25 s to 8–12 s", glide force "<15 N", 27G 8 mm — using **"Hydrodynamic pressure drop modeling using Hagen–Poiseuille equation"**, i.e. a Newtonian model, on the vendors' own page. Ask question 4 (would a reviewer ask for an error bar on that?). Count as **one** call if they answer jointly, two if separately. | report 01 line 53; shl-medical.com blog 2026-01-05 | **unresolved** |
| C2 | **Veryst Engineering** (~18 staff, Needham MA) | Sells a COMSOL subcutaneous-injection model — injection force, delivery time, bleb formation, bolus dispersion, poroelastic tissue — with **non-Newtonian drug rheology, validation and UQ all absent from the page**. The exact hole this venture fills, at a firm small enough to subcontract. This is the rank-2 (V&V 40) conversation. | report 01 line 47; report 02 lines 11, 57 | **unresolved** |
| C3 | **Springboard** (UK) or **Team Consulting** (UK) | Autoinjector device consultancies with no rheology-specific claim found. Same subcontract logic as C2, UK time zone. | report 01 lines 48–49 | **unresolved** |

## Tier D — regulatory / CM&S (§7 quota of 2)

| # | who | why them | source | contact |
|---|---|---|---|---|
| D1 | **Numalogics** (Canada, founded 2010 by three spinal surgeons; ASME V&V 40 contributor) | Holds the **first MDDT qualification of a mechanical test as an FE model** (numaScrew, ASTM F543 screw pullout, qualified 10 Apr 2025) — and needed **Sawbones as a physical-comparator partner** to get it. That is the precedent that says a small shop can reach regulatory-grade work, and the partner structure is the template. Ask question 6. | report 02 lines 11, 29, 113 | **unresolved** (has an open FE engineer role, so the site has a live contact) |
| D2 | **Jeff Bodner — Medtronic Distinguished Scientist, ASME V&V 40 subcommittee** (or a boutique from the 2024 FDA/MDIC credibility panels: CardioMed, NC MedTech Solutions, Lighthouse Regulatory) | Report 02's own open question 147: ask what fraction of blood-contacting CFD submissions use a non-Newtonian viscosity model, and **whether viscosity choice has ever appeared in a deficiency letter**. A yes converts F7 from a gap into a requirement; a no is close to fatal for rank 2 and should be recorded as such. | report 02 lines 58, 36, 147 | **unresolved** |

## Tier E — platform company (§7 quota of 1)

| # | who | why them | source | contact |
|---|---|---|---|---|
| E1 | **Halozyme** (acquired Surf Bio for up to **$400M**, Dec 2025, for 500 mg/mL hyperconcentration) or **Elektrofi** | At 400–500 mg/mL the working fluid is a dense suspension — the PhD — and F5 records **zero services** for it. Rank 3's kill test is exactly this: "a platform company takes a call and names an unsolved problem." | report 01 line 13; synthesis F5 | **unresolved** |

---

## Order of sending, and why

1. **A1, A2, A3** first, in one afternoon. Verified or near-verified contacts, the warmest
   possible pretext, and the two answers that most change the plan (did anyone pay to fix
   this? is it still open?).
2. **C1, C2** next. C2 is the one conversation that tests §5.2 of `08_AI_EXPOSURE.md` — the
   recommendation that rank 2 be promoted above rank 1.
3. **D1, D2** third. D2's deficiency-letter question is the single highest-information
   question in the whole script and is cheap to ask.
4. **B1, B2, E1** last, and treat question 8 (referral) as their main purpose.

**Do not send more than three at once.** The first replies will change the script, and a
batch of fifteen identical notes is the cold-outreach failure mode F9 documents.

## What this sheet deliberately does not contain

- **No invented email addresses, phone numbers or LinkedIn URLs.** Two addresses are here
  because the authors printed them in a journal; every other row says *unresolved* and
  costs the founder one search each.
- **No claim that any of these people will reply.** The base rate is unknown and
  `00_SYNTHESIS.md` §8 already flags that no buyer was ever found.
- **No pricing.** Prices are discussed only if the other side asks; `OFFER.md` §5 is marked
  INFERRED and must be corrected before any call quotes a number.
