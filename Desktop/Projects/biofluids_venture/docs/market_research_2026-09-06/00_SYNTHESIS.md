# 00 — Synthesis: what to build, ranked

Date: 2026-09-07. Inputs: reports 01–06 in this folder (1,642 lines, every claim URL-sourced
and dated; vendor market sizes discarded). This file carries only conclusions; the evidence and
the "not found" flags are in the numbered reports. Everything the founder must decide is in §5.

## 1. The one-paragraph answer

The stack you have built (cited registry → identification with identifiability verdict and
corrected confidence intervals → solver export with round trip → closed-form device prediction
with UQ) is the artefact nobody sells, in any industry, at any price. Six contract-rheology
labs, three rheometer vendors, two physics-AI unicorns and the whole V&V 40 consultancy tier
were checked; none states an uncertainty on a fitted constitutive law, none checks that the
solver reproduces it, and none tells the buyer which parameters the answer does not depend on.
That is a real gap. What the research does **not** find is a single buyer who has paid for it.
Every documented small physics shop reached revenue by selling fixed-price studies and teaching,
1.5–5 years in, and moved to product afterward if at all. So the plan is: the stack is the
differentiator, not the product; the product for the first 24 months is fixed-price studies that
each leave a registry record behind, in one wedge, chosen by whoever answers ten discovery
calls first.

## 2. Settled findings

| # | finding | source |
|---|---|---|
| F1 | No seller anywhere delivers a constitutive fit with uncertainty, identifiability, or a solver round trip. Closest: Centre for Industrial Rheology's "modelling-ready data for CFD" article — a claim, not a deliverable. | 03, 04 §Q3/Q5, 06 §2 |
| F2 | AI commoditises case-file writing (Foam-Agent 88 %), curve fitting (TRIOS ships nine models), generic UQ, hourly meshing. It raises the value of verdicts, envelope guards, provenance, credibility trails. | 04 table |
| F3 | Physics-AI surrogates (PhysicsX/Neural Concept/Luminary, $470M between them) are Newtonian single-phase RANS; their own 2026 posts ask for failure diagnostics. Not a threat for 5+ years; a possible customer for a validation harness. | 04 §Q1 |
| F4 | Injectability prediction is a documented failure at the high-shear end: JHU and GSK papers report Hagen–Poiseuille models breaking above ~230 mg/mL and in 27–30G needles because rheometer data stops at 1e4 1/s and the needle runs at 1e5. No injectability service (Coriolis, QATCH, Eurofins, Smithers, DDL) offers a predictive layer or an error bar. | 01 |
| F5 | The working fluid is becoming a dense suspension (Surf Bio → Halozyme $400M Dec 2025; Elektrofi; 400–500 mg/mL spray-dried particles). That is jamming/shear-thickening physics — the PhD — and no service exists for it. | 01 §Unmet 7 |
| F6 | FDA's Nov-2023 guidance adopts V&V 40; MDIC 2021: 72 % of device sponsors do CM&S in-house, 23 % use consultants, 28 % have no internal resource; top barrier is "what does the regulator expect". Small firms do win (Numalogics MDDT, April 2025, with a physical-comparator partner). No public prices. No standalone V&V 40 authoring product exists. | 02 |
| F7 | Blood-contacting CFD's credibility gap is rheology: the FDA nozzle benchmarks used one viscosity at 200 1/s, models spread ~6× at low Re, the May-2026 hemolysis RST dataset says nothing about viscosity. | 02 §win 2 |
| F8 | Outside biomedicine, three markets buy the skillset: battery electrode slurry/coating (only dated non-vendor cost-of-problem: Fraunhofer FFB 2024 ~€30k/day per scrap point), agrochemical spray retention (CropLife Europe pays SSAU, a small UK lab), inkjet/dispensing of particle-loaded inks (most small-shop precedent). | 03 |
| F9 | Nobody is hired as a "rheologist" (176 Glassdoor hits, all US industries). Buyers hire process engineers; contracting is not their default. Cold outreach has the only documented failure (CFD Support: a year of calls, no reply, a 400-hour free pilot). Teaching is the only channel with repeated evidence (CoolSim webinars, CFD Direct 81 training days/yr, Cunningham ~3,000 trained). | 03 §2, 05 |
| F10 | Price floor and ceiling: university cores rent rheometers to industry at $43–64/h; per-sample contract rheology $170–340; PhD consultant median $200/h (IEEE-USA); specialist day rates $1,200–5,000; US CFD studies $1.5k–20k fixed. Offshore CFD from $300–500 per case. | 05, 06 |
| F11 | Non-dilutive money is large but slow: NSF SBIR Phase I $305k at ~9–12 % end-to-end (Project Pitch gate), NIBIB $256k–$323k; Indiana match $50k, Illinois $75k; PI must be >50 % employed at award. Consultant day caps and subaward fractions apply. | 05 |
| F12 | A credible four-rig bench (needle force, capillary, spray/drop imaging, syringe pump) is $30–70k used, $150k+ new with a Phantom. Regulated nasal-spray and ISO 11608 testing is owned by GMP/ISO 17025 incumbents and not enterable by a 1–3 person shop. Simulants only → no OSHA 1910.1030, no exposure-control plan, no IRB. | 06 |
| F13 | Open-source maintenance does not pay (OpenFOAM Foundation ~€250k/yr against €500k need; OpenFOAM Ltd 5 people, £874k). No solo founder with disclosed revenue from an agent-run simulation business was found. | 05 |

## 3. Candidate products, ranked

Scoring on: dated non-vendor demand evidence; reuse of the existing stack; time to first invoice
at 5–10 h/week; durability against AI over 5–10 years; whether the founder can do it alone.

| rank | product | buyer | demand evidence | stack reuse | first invoice | durability | kill test |
|---|---|---|---|---|---|---|---|
| **1** | **Injectability calibration study, fixed price**: client rheometry (or a simulant measured here) → identified law with CI → needle force/time vs concentration, temperature, gauge with an error bar and an identifiability verdict → worst-case corner selection for ISO 11608 / EDDO test matrices | formulation groups at biologics companies and CDMOs; device makers (Stevanato, Nemera, BD platform partners) | F4, F5: documented model failure, no seller with an error bar, consolidating $400M market | 95 % (L0–L2, L6 built; item 12/14/15 already the case study) | fastest — the deliverable exists today for blood; needs one mAb-like dataset | high: the verdict and CI are what AI cannot fake | 10 discovery calls; ≥3 say "we paid someone / we would" for the high-shear extrapolation. If 0/10, drop W1. |
| **2** | **Rheology-aware V&V 40 credibility package for blood-contacting or drug-delivery CFD**: subcontracted OpenFOAM + full viscosity curve propagated with UQ + the QOI/COU/credibility-factor plan authored | the 28 % of sponsors without CM&S, start-ups, and two established consultancies as subcontractor | F6, F7: FDA-recognised standard, named barrier, small-firm precedent; no prices | 70 % (L3 export + round-trip harness; L5 manifest is the missing layer) | medium — needs one retroactive manifest as a portfolio piece, then a subcontract | high: knowledge product, regulator-facing | one consultancy agrees to subcontract or one sponsor buys a plan. Not realistic alone: PMA, MDDT, in silico trials. |
| **3** | **Dense-suspension injectability** (spray-dried particle formulations, 400–500 mg/mL): jamming/shear-thickening closure, needle-clogging risk, bench simulant | Halozyme/Elektrofi-class platforms and their licensees | F5: acquisition signal; zero services found | 60 % (dst_fit Wyart–Cates, drop-impact jamming code) plus new: particle friction in a needle | slow — research-led sale, a paper first | very high: it is the PhD, and it is unopened | a platform company takes a call and names an unsolved problem |
| **4** | **Battery slurry / coating-window study** with fitted rheology and Basilisk slot-die/bead physics | cell start-ups, coater OEMs, binder suppliers | F8: only market with a dated cost-of-problem | 40 % (fitting yes; coating-bead solver new) | medium | high | one paid pilot from a dry-electrode or Si-anode start-up. Compete: process engineers in-house (F9). |
| **5** | **Cited registry + solver export as a product** (per record / per seat / MCP-API) | anyone running non-Newtonian CFD; physics-AI vendors as training data | F1: absence-confirmed; no positive buyer; F13: nobody paid for open data | 90 % | slow — database cold start, and every precedent says product after services | highest — provenance + data compounds | 3 of the first 10 study clients ask to keep the record and the export. Build it as the by-product of ranks 1–4, never as the first product. |
| 6 | Validation harness for physics-AI in complex fluids | PhysicsNeMo/SimAI vendors | F3: they ask for it; no budget line | 70 % | 2–4 yr | high | a vendor funds one certified run |
| 7 | Agrochemical spray retention with audited imaging | formulation/adjuvant groups | F8 contract-lab precedent | 30 % | needs rig ($15–25k imaging) | high | SSAU-class lab is the incumbent; only if a client brings the fluid |
| 8 | Bench line as a stand-alone service | — | F10/F12: commodity, rate capped at $64/h by university cores | — | — | low | Do not build. A rig exists only to validate a prediction. |
| — | Hourly CFD, generic rheology testing, open-core maintenance, regulated GMP spray testing, in silico trials | — | F2, F9, F12, F13 | — | — | commoditised or closed | Excluded. |

**Ranks 1–3 are one wedge (W1 biologics delivery) at three depths; rank 2 is the same
deliverable wrapped for a regulated buyer.** Rank 4 is the hedge if W1 returns 0/10. Rank 5 is
the endgame every other rank feeds.

## 4. Contradictions between reports, resolved

| contradiction | resolution |
|---|---|
| 04: services are being commoditised. 05: every viable shop started with services. | Both true. What is commoditised is *hourly* CFD and *raw* rheometry. Fixed-price studies whose deliverable is a verdict + CI + record are not on any commoditisation list. Rule: never sell hours; every study leaves a registry record; the record catalogue is the product that services fund. |
| 03: battery has the strongest cost-of-problem. 01/02: biopharma has a compliance budget line and a documented model failure. | Resolved on *who has already paid an outside shop*: biopharma pays Coriolis/Eurofins/consultants today and has a regulator asking for credibility; battery buys via in-house process engineers (F9). Biopharma first; battery is the hedge. |
| 02: "start-ups without CM&S are underserved." 05: cold outreach to start-ups failed for a year. | Channel, not market. Enter through teaching (a webinar/short course on high-shear extrapolation and identifiability) and through subcontracting to two consultancies, never cold calls. |
| 06: a bench costs $30–70k. 05: the founder's ceiling is ~$200/h and no revenue exists. | Do not buy a bench before an invoice. Rent university core time ($43–64/h) for the first validation; buy a used force stand (~$1.5–9k) only when a paid study needs it. Blood work uses simulants only (F12). |
| 04: registry-as-product is the durable asset. 05: product-first has no precedent. | Sequence, not either/or. The registry is built as the by-product of paid studies and sold separately only after ≥3 clients ask to keep the export. |

## 5. Decisions only the founder can make

| id | decision | recommendation | why |
|---|---|---|---|
| D-A | Wedge: W1 injectability (ranks 1–3) vs battery (rank 4) vs keep both open | **W1, with battery as the pre-registered fallback if R1 returns < 3/10** | F4–F7 give W1 a documented failure, a regulator, a consolidating market and a PhD-shaped extension; battery has a cost-of-problem but no contracting habit |
| D-B | First revenue form | **Fixed-price studies ($5–20k) + one paid short course; refuse hourly and free pilots** | F9, F10; the only documented failures are hourly/cold/free |
| D-C | Whether to spend money on a bench before revenue | **No.** Core-facility rental for the first validation | F10 caps bench rates below what a shop can charge; F12 |
| D-D | Company formation timing vs PhD end | **Incorporate as soon as the first study is signed or the SBIR Project Pitch is invited, whichever first; the >50 % employment rule binds at award, not proposal** | F11 |
| D-E | Human blood or simulants | **Simulants only, as a design rule** | F12: avoids OSHA 1910.1030 and IRB entirely |
| D-F | dst_fit item 17 upstream fix (D5) | unchanged; blocks nothing commercial | — |
| D-G | Name (D1) | needed before the course/webinar, not before the calls | — |

## 6. Twelve-month plan, 5–10 h/week + overnight agents

Founder hours go to calls, teaching, and verdicts; agents do everything else. Each month has a
gate; a failed gate is a recorded result and triggers the named fallback.

| month | founder (≤10 h/wk) | overnight agents | gate |
|---|---|---|---|
| 1 | Write the ten-call script (§7). Send 15 warm asks via alumni/advisor/LinkedIn to formulation scientists, device engineers, one regulatory consultant. Name the company (D1). | Build the mAb-like demo: fit a published high-concentration viscosity curve (JHU/GSK data from 01), needle force vs concentration/temperature/gauge with CI; one-page "what this force does not depend on". | 5 calls booked; demo reproduces the paper's reported failure point with an error bar |
| 2 | Hold calls 1–5. Log verbatim: who paid, how much, for what. | Retroactive V&V 40 manifest (L5 skeleton) on the needle study; short-course outline: "high-shear extrapolation and identifiability for injectability". | Notes on 5 calls; manifest lists every credibility factor as missing-or-met, none by default |
| 3 | Calls 6–10. **R1 verdict**: ≥3/10 have paid or would pay for the high-shear/identifiability layer. | Registry to 25 records (mAb curves, simulants, blood); ISO 11608/EDDO corner-selection tool. | **R1 pass → W1. Fail → switch month 4–6 to battery: same plan, coating-bead physics, Fraunhofer FFB contacts.** |
| 4 | Pitch a fixed-price study ($5–10k) to the two warmest calls. Free 45-min webinar on the course topic (the acquisition channel, F9). | Deliverable template: 90-second summary, CI table, identifiability verdict, exported OpenFOAM dictionary + round-trip report. | 1 proposal out; ≥15 webinar registrants |
| 5–6 | Deliver study 1 using core-facility rheometry if the client has none. Ask: "was the verdict worth paying for?" Draft SBIR Project Pitch (NSF or NIBIB) on dense-suspension injectability (rank 3). | Dense-suspension closure: Wyart–Cates → needle; clogging criterion; paper draft. | **Invoice 1 paid.** Pitch submitted. |
| 7–8 | Approach two device consultancies with the manifest as the pitch (rank 2, subcontract). Paid short course #1 ($300–800/seat, 10–20 seats). | Blood-contacting OpenFOAM case (FDA nozzle) with full viscosity-curve UQ — the F7 gap, as a portfolio piece. | 1 subcontract conversation past intro; course revenue ≥ $3k |
| 9–10 | Study 2 (repeat client or webinar lead). Incorporate if not yet (D-D). E&O insurance before any device-tier work. | Registry export licence terms drafted; MCP/API prototype for records. | **Invoice 2.** ≥1 client asks to keep the export (rank 5 signal) |
| 11–12 | Decide product vs services mix from evidence: count who asked for the record. SBIR Phase I outcome. Re-anchor phases to the PhD end date. | Rank-3 paper submitted; benchmark suite v0 public. | Two invoices, one course, one pitch, one paper: the 05 "1.5–5 yr" clock starts from here with revenue evidence, not hope |

Money: $0 capex through month 6; core rental ≤ $1k; used force stand ≤ $2k only inside a paid
study. Target month-12 revenue $15–30k, which is below living income and is the honest number:
F10/05 say the viable-consultancy point is 18 months at the earliest.

## 7. What the web cannot answer — the ten calls

Ask each of ten people (target: 4 formulation scientists at biologics companies or CDMOs,
3 device/combination-product engineers, 2 regulatory or CM&S consultants, 1 platform company
such as Halozyme/Elektrofi), in this order:

1. When your injection-force model disagreed with the bench, what did you do, and who did it?
2. Has your group ever paid an outside party for rheology-to-force prediction? How much? What
   did the deliverable look like?
3. Where does your viscosity data stop, and how do you get from there to the needle?
4. Does anyone put an error bar on the predicted force? Would a reviewer or FDA ask for one?
5. For suspensions (spray-dried, >400 mg/mL): what breaks first, and is it modelled at all?
6. Who owns the credibility argument for your CM&S — internal, consultant, nobody?
7. If a one-page verdict said "your force depends on n and eta_inf and not on tau_y, and here is
   the 95 % band", what would that be worth to you in a design-verification matrix?
8. Who else should I talk to?

Score each call: paid before (Y/N), would pay (Y/N), price named, deliverable named, referral.
R1 = number with paid-before OR would-pay ≥ 3.

## 8. What this research did not find (flags)

- No public price for any injectability, V&V 40, or credibility-plan service; the $5–20k
  study price is inferred from generic CFD-study quotes (F10), not from this segment.
- No buyer quote anywhere saying "we wanted an uncertainty on the constitutive law".
- No solo founder with disclosed revenue from an agent-run simulation business.
- Upwork per-profile rates (403); the ISO 11608-1:2022 text itself (paywalled).
- Battery cost-of-problem is one source (Fraunhofer FFB 2024); not corroborated.
