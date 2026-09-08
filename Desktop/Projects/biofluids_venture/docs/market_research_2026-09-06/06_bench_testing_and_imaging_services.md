# 06 — Bench testing and imaging as a service: market, rate cards, set-up cost

Research date: 2026-09-06/07. Method: ~40 web searches and page fetches. Every number below carries the URL it came from and the date it was read. Where a page returned 403/404 to the fetcher, the number is taken from the search-engine snippet and marked **[snippet]** — treat as lower-confidence until opened in a browser. "Not found" means it was searched for and no public figure exists; it is not a guess.

Scope caveat: this is a desk study of PUBLIC prices. Most contract labs quote per project, so the rate-card table is thin by construction. That thinness is itself a finding (Q1).

---

## Summary

1. Contract rheology is a quote-per-project commodity with almost no public rate cards; the one found is ASM Services (UK): £90/h instrument time + £200–300 fixed for calibration/report (asm-serv.co.uk, read 2026-09-06). US university cores rent the same class of rheometer to external users at $43–63.90/h unassisted and $147–192/h assisted (Notre Dame, UF), which caps what a small shop can charge for raw instrument time.
2. No contract rheology provider found (Intertek, Measurlabs, Rheology Testing Services, Centre for Industrial Rheology, Rheology Solutions, JK Analytical) advertises model fitting with uncertainty, identifiability, or solver-ready parameters; the closest is a Centre for Industrial Rheology article titled "Modelling-Ready Rheology Data for CFD", which is a claim, not a documented deliverable. That gap is the founder's wedge, not bench time.
3. High-speed spray/drop imaging services exist (Oxford Lasers 20+ years contract imaging; Proveris GMP/ISO 9001 nasal-spray characterisation; Spraying Systems Co.; Silsoe SSAU for agrochemical drift) and none publish prices. The regulated pharma segment (spray pattern/plume geometry under the FDA nasal-spray guidance) is owned by GMP-registered incumbents and is not enterable by a 1–3 person shop.
4. Injection-force / ISO 11608 / ISO 11040-4 testing is sold by Eurofins MDS, Smithers, DDL, Nelson-class labs — accredited (ISO 17025, UKAS, cGMP), no public prices, and none mention viscosity- or rheology-based prediction of the force. The prediction layer is missing from the whole segment.
5. Equipment for a minimal lab: used MCR 302 €8,000–17,000 [snippet, LabX]; used TA DHR ~$40,000 with install/warranty [snippet, LabX]; used Instron 5944 ~$9,200 [snippet, eBay]; Mark-10 ESM303 stand ~$1,420 + M5-50 gauge $600–1,670; Chronos 2.1-HD new $4,995–6,000; used Phantom Miro M320S ~$13,250 [snippet]; NE-1000 syringe pump $932 new; calibrated glass capillary viscometers $221–521 each. A credible bench for the four rigs named is roughly $30k–70k used, $150k+ new with a Phantom.
6. Space: Illinois Research Park EnterpriseWorks offers 570 sq ft wet labs with hoods (price on request); Purdue Technology Center lab space listed at $38.09/RSF/yr [snippet, Merrillville]; 16 Tech Indianapolis build-to-suit "at least $40/sq ft/yr" [snippet, IBJ]; mHUB Chicago $335–1,350+/mo but it is a prototyping shop, not a wet lab.
7. Compliance: OSHA 29 CFR 1910.1030 defines "blood" as HUMAN blood/components/products, and OSHA's own guidance lists "synthetic blood" among alternatives to OPIM. Using simulants only avoids the bloodborne-pathogen standard, an exposure-control plan, and IRB entirely. That is a real cost avoidance and should be a design rule.
8. Evidence that buyers pay a premium for sim + bench from one vendor: Veryst (COMSOL blog) is the canonical example but publishes no client outcomes, revenue, or lab size. No case of a SMALL shop winning specifically because it had both was found. This question is open.
9. Academic-adjacent revenue: NSF SBIR caps consultants at $1,000/day (8 h) and requires ≥2/3 (Phase I) / ≥50 % (Phase II) of the budget to stay with the small business, so subaward work to a founder-lab is structurally capped at ≤1/3 of a Phase I. No survey of small testing labs' revenue mix by customer type was found.
10. Verdict: a stand-alone bench-testing line at core-facility rates ($45–190/h) cannot carry a 1–3 person company; a bench that exists to PRODUCE solver-ready, uncertainty-quantified parameters and validation data for the simulation product is defensible, because that deliverable is absent from every rate card found.

---

## Evidence by question

### Q1. Contract rheology testing — who, what they charge, what they don't deliver

**Providers found and what their pages say (read 2026-09-06):**

| Provider | Type | Pricing statement | Turnaround | Deliverables named | Source |
|---|---|---|---|---|---|
| ASM Services (UK) | independent | "£90 per hour" instrument time + "£200-£300" fixed for calibration, report, analysis | not stated | free feasibility study; report | https://www.asm-serv.co.uk/contract-testing |
| Rheology Testing Services (Chapel Hill NC) | independent, one-person by appearance (gmail contact) | quote on request | "approximately one week for draft report" | overlay plots, summary tables, Excel data; "does not explicitly mention model fitting, uncertainty quantification, or simulation parameters" | https://www.rheologytestingservices.com/applications-and-support-overview |
| Rheology Solutions (AU) | independent | quote after 48 h triage; "estimate of the projected time scale and associated fee" | not stated | contract testing kits by post | http://www.rheologysolutions.com/company/applications-laboratory/ |
| Centre for Industrial Rheology (UK, rheologylab.com) | independent | none public | none | "full report, including graphics, data tables and explanations"; ongoing Q&A; no model fitting or uncertainty mentioned; one article "Rheology and Viscosity Measurements for CFD" | https://www.rheologylab.com/services/viscosity-testing-lab/ ; https://www.rheologylab.com/articles/modelling-ready-rheology-for-computational-fluid-dynamics/ |
| Measurlabs (FI, broker) | marketplace | "Prices excluding VAT" but none listed; quote in one business day | not stated | "a clear report on the results" | https://measurlabs.com/methods/rheology-testing/ |
| Intertek (Allentown PA and EU) | multinational | none public | not stated | ASTM-standard rotational rheology to 400 °C | https://www.intertek.com/analytical-laboratories/rotational-rheology/ |
| JK Analytical | independent | page title says "affordable prices"; page not fetched (limit) | — | — | https://www.jkanalyticallab.com/rheology-viscosity-testing-services-6018325.html |
| Rheology Labs (IEH company) | independent | not fetched (limit) | — | — | https://rheologylabs.com/ |
| Anton Paar | vendor | no contract-testing rate found; vendor pages sell throughput to contract labs (sample changers) | — | — | https://www.anton-paar.com/us-en/rheology-instruments/ |

**University core facilities (the price floor for instrument time):**

| Facility | Instrument | Rate | Source |
|---|---|---|---|
| Notre Dame MCF | TA HR-2 Discovery Hybrid | internal unassisted $26–29/h; internal assisted $77–133/h; external unassisted $43–45/h; external assisted $147–192/h; $258/day external assisted | https://mcf.nd.edu/internal-users/rate-sheet/ and https://mcf.nd.edu/capabilities-instruments/rheometry-dmta-tribology/hybrid-rheometer/ (two pages, two numbers each — both quoted) |
| U. Florida NRF | ARES LS1 | academic/gov $16.20/h; affiliated commercial $47.90/h; commercial $63.90/h; international academic $21.60/h [snippet; page 404 on fetch] | https://nrf.aux.eng.ufl.edu/ccb/resource.asp?id=105 |
| U. Illinois MRL | TA DHR-3 | listed, rates on a separate "usage rates" page not fetched | https://mrl.illinois.edu/facilities/equipment/ta-dhr-3-rheometer |
| Purdue Birck / Chem RIC | none rheometer-specific found | — | https://www.chem.purdue.edu/ric/rates.html |
| UChicago MRSEC | "Rheometry Facility" and "High Speed Imaging Facility" exist; rates not on the MPML page | — | https://mrsec.uchicago.edu/facilities/materials-preparation-and-measurement-laboratory/rates/ |

**What nobody delivers.** Across nine provider pages, deliverables stop at plots, tables and a written interpretation. Not one page names a constitutive model fit, a parameter confidence interval, an identifiability statement, or a file in a solver's input format. The only phrase pointing that way is the Centre for Industrial Rheology's "Modelling-Ready Rheology Data for CFD" article, and the services page it belongs to lists no such deliverable.

**Is it a commodity?** Yes on price signals: no differentiation beyond temperature range and turnaround; university cores rent the same instruments at $43–64/h external unassisted; the one public independent rate is £90/h. A shop charging materially more must sell something the rate card does not contain.

### Q2. High-speed imaging / spray and drop characterisation services

| Provider | What they sell | Price | Buyers evidence | Source |
|---|---|---|---|---|
| Oxford Lasers Imaging Division (UK) | contract high-speed imaging, particle/droplet size and velocity, inhaler (pMDI/DPI/SMI) characterisation, spray pattern and plume geometry; in-house, on-site, rental | none published | 20+ years; cites five papers, no named clients or case studies | https://oxfordlasers.com/services/laser-imaging-services/ |
| Proveris Laboratories (Hudson MA) | spray pattern, plume geometry, droplet size, APSD, IVBE for OINDP; GMP + ISO 9001, FDA-registered | none published | first SprayVIEW installed at FDA St. Louis; FDA draft guidance took technical input from founder (marketing claim) | https://proverislabs.com/ ; https://www.proveris.com/about-us/ |
| Spraying Systems Co. Spray Analysis & Research Services | drop size, velocity, impact, coverage, evaporation; wind tunnel; high-speed imaging | none published | vendor-affiliated lab (nozzle manufacturer) | https://www.spray.com/services/advanced-spray-testing-services ; https://www.sprayanalysis.com/testing-services/instruments-and-facilities/high-speed-imaging |
| Silsoe Spray Applications Unit Ltd (UK) | agrochemical spray/drift wind tunnel, droplet sizing, field trials | none published | a CropLife Europe-commissioned contract report is public — direct evidence an industry body buys | https://croplifeeurope.eu/wp-content/uploads/2021/01/S0252-contract-report-final-complete.pdf ; https://ukagritechcentre.com/community-directory/silsoe-spray-applications-unit-ltd/ |
| ImageXpert (NH) | inkjet jetting tests, waveform optimisation, drop volume/velocity/satellites; also sells JetXpert dropwatchers | none published | Quantica partnership (3D printing) | https://imagexpert.com/process-development/ ; https://imagexpert.com/drop-analysis/dropwatcher/ |
| Melbourn Scientific (UK) | first UK CRO to offer FDA spray pattern/plume geometry | none | — | https://www.cambridgenetwork.co.uk/news/melbourn-first-offer-fda-approved-spray-pattern-and-plume-geometry-analysis |
| Dantec, LaVision, TSI, Malvern Spraytec, Sympatec, InnovaPrep | **not found** as service labs with public rates in the searches run; they sell instruments | — | — | — |

**Camera rental as the shadow price of imaging time:** Phantom rentals ~$2,400–3,000+/day (cinema market) — https://beverlyboy.com/how-much-does-it-cost/how-much-does-a-phantom-camera-cost/ ; ATEC rents Phantom/Photron/Chronos for test and measurement — https://www.atecorp.com/ndt-inspection/visual-inspection/high-speed-cameras. No scientific imaging-service day rate was found.

**Room for a small lab?** The regulated pharma niche (FDA nasal-spray CMC guidance) needs GMP/ISO 17025 and is held by Proveris/Melbourn-class labs. Ag-spray drift is held by wind-tunnel facilities (Silsoe, Spraying Systems). Inkjet has ImageXpert selling both the instrument and the service. The unoccupied slot is drop impact / jet / needle-spray imaging with audited calibration delivered as validation data for a model — no provider found sells that.

### Q3. Injection-force / syringeability test houses

| Provider | Tests | Accreditation | Price / turnaround | Rheology-linked prediction? | Source |
|---|---|---|---|---|---|
| Eurofins Medical Device Services | break-loose, glide, Luer, hub retention, accuracy, siliconisation, ISO 11608 | ISO/ASTM/AAMI | none public; "plan testing timelines" | not mentioned | https://www.eurofins.com/medical-device/testing/combination-products/syringe-testing/ |
| Smithers | ISO 11608 parts 1–6, cartridge testing | UKAS ISO 17025, GMP | none public (page 403 on fetch) | not mentioned | https://www.smithers.com/services/testing/standard/iso/iso-11608-1 |
| DDL (Minnetonka MN) | ISO 11608 NIS testing | ISO 17025, FDA-registered cGMP | none public | not mentioned | https://www.ddltesting.com/medical-device-testing/needle-based-injection-system-testing/ |
| Nelson Labs, West/BD/Stevanato service arms | **not found** in searches run (limit) | — | — | — | — |
| Instrument vendors (Instron, ZwickRoell, Mecmesin, Stable Micro Systems) | publish ISO 11040/11608 method guides to sell testers, not services | — | — | — | https://www.instron.com/en/testing-solutions/iso-standards/iso-11040/ ; https://www.zwickroell.com/industries/medicalpharmaceutical/therapy-systems/syringe-testing-iso-11040/ ; https://www.mecmesin.com/publications/break-loose-and-sustaining-force-tests-pre-filled-drug-cartridges-iso-11608-3 |

Rates: **not found** for any ISO 11608/11040 test house. General ISO 17025 mechanical-test hourly rates were also not public in the searches run.

What is missing is uniform across the segment: every page frames the force as a device pass/fail measurement; none links it to the fluid's flow curve or offers a prediction before the device is built. That is exactly the biorheo needle-force calculation the founder already has (see `biorheo` item 12: family choice moves the force up to 30 %, parameter uncertainty a median 29 %).

### Q4. Equipment cost for a minimal lab (2026)

See the equipment table below. Notes on provenance:
- LabX, eBay and Exapro listing pages returned 403 to the fetcher; those numbers are from search snippets and are marked. Open them in a browser before budgeting.
- New rheometer prices (MCR 72/92/302e, DHR, Kinexus) are quote-only; Anton Paar's page shows an "online price" widget that did not render. **Not found.**
- A Harvard PHD Ultra "new" figure of $21,646.99 appeared on a medical reseller (graylinemedical.com) and is inconsistent with the used range ($1,300–1,500); treat as unreliable.
- Capillary rheometer (pressure-driven, e.g. Dynisco LCR): **not found**; glass capillary viscometers are cheap and listed.
- Lighting/optics: **not found** in this pass (backlight LED panels and macro lenses are commodity; budget separately).

### Q5. Space and compliance

**Space (Indiana/Illinois):**

| Site | Offer | Price | Source |
|---|---|---|---|
| Illinois Research Park, EnterpriseWorks (Champaign) | 23 labs, each 570 sq ft, fume hoods, bench gas/vacuum, 3-phase, eyewash/shower, shared equipment lab; campus core access at faculty-staff rates | on request (cfaullin@illinois.edu) | https://researchpark.illinois.edu/locate-here/lab-spaces/ |
| Purdue Technology Center (Merrillville IN) | lab space | $38.09/RSF/yr [snippet] | https://www.commercialsearch.com/commercial-property/us/in/merrillville/purdue-technology-center/ |
| Purdue Research Park West Lafayette | Class-A office and wet lab, gross leases | rate not found | https://purdueresearchpark.com/locations/west-lafayette/ |
| 16 Tech / Waterway Labs (Indianapolis) | build-to-suit lab | "at least $40 per square foot per year" [snippet, IBJ] | https://www.ibj.com/articles/16-tech-pencils-in-plans-for-next-two-buildings-after-early-success ; https://www.waterway-labs.com/ (403) |
| IU LAB Bio Start-up Center (Indianapolis) | incubator programme, up to 40 companies/yr | price not found | https://iulab.iu.edu/bio-start-up-center/index.html |
| mHUB (Chicago) | co-working $335/mo; reserved desk $450/mo; office $1,350+/mo; shop punch card $450/10 uses | prototyping shop, no wet-lab restrictions stated but no wet lab offered | https://mhub.org/become-a-member |
| Boston reference point | LabCentral bench $4,600/mo (2019), ~$100/sq ft | for contrast only | https://www.excedr.com/blog/how-much-does-lab-space-cost |

At $38–40/sq ft/yr, a 570 sq ft lab is ~$1,800–1,900/mo before services. That is the right order for a 1–3 person bench.

**Compliance — simulants vs human fluids:**
- OSHA 29 CFR 1910.1030 covers employees "reasonably expected to have exposure to human blood or other potentially infectious materials"; "blood" is defined as human blood, components, and products made from human blood — https://www.osha.gov/bloodborne-pathogens/general ; https://www.osha.gov/laws-regs/standardinterpretations/1993-02-01-0
- OSHA's own compliance guidance lists "non-infectious animal blood, synthetic blood or computer simulations" as alternatives to OPIM — https://www.osha.gov/laws-regs/standardinterpretations/2009-02-27
- Plasma-derived products ARE "blood" under the standard — https://www.osha.gov/laws-regs/standardinterpretations/2002-12-27 — so buying pooled human plasma for rheology brings the standard back in.
- Consequence: a simulants-only rule (glycerol-water, xanthan, Carbopol, PEG/PEO blood analogues, commercial synthetic blood) means no exposure-control plan, no HBV vaccination offer, no biosafety cabinet, no IRB, no BSL-2 lease clause. Animal blood (non-infectious) is also outside 1910.1030 but brings USDA/vendor and disposal rules — **not researched here**.

### Q6. Do clients pay for sim + bench from one vendor?

- Veryst Engineering (Needham MA) is the reference model: "combines material testing and simulation", tests with SHPB, calibrates with MCalibration, delivers the calibrated model; software later spun out as PolymerFEM (an Ansys technology partner) — https://www.comsol.com/blogs/veryst-combines-material-testing-and-simulation-for-reliable-results ; https://www.veryst.com/services/testing/material-model-calibration ; https://www.veryst.com/news-events/veryst-software-spun-new-company-polymerfem
- The COMSOL piece contains **no client names, no outcomes, no lab size, no revenue**. It is vendor marketing.
- Veryst's revenue, headcount and testing share: **not found**.
- A small (1–3 person) shop that won specifically because it had both: **not found**. The nearest evidence is the Centre for Industrial Rheology positioning "modelling-ready rheology for CFD", i.e. a rheology lab reaching toward the simulation buyer from the other side.
- Structural point that IS supported: Veryst's bench exists to calibrate a model it sells the client (the deliverable is the calibrated model, not the test report). That is the same shape as the proposed venture.

### Q7. Academic-adjacent revenue

- NSF SBIR/STTR budget rules (read 2026-09-07): consultants "may not exceed $1,000 per day" (8 h); Phase I SBIR "at least two-thirds of the budget be allocated to the small business"; Phase II SBIR at least 50 %; STTR at least 40 % small business + 30 % research institution; fee up to 7 % (Phase I) / 10 % (Phase II); indirect without a negotiated rate up to 50 % of salaries or 15 % de minimis — https://seedfund.nsf.gov/solicitation-budget/
- Implication: as a SUBCONTRACTOR to someone else's NSF Phase I, a founder-lab can take at most ~1/3 of the award; as a CONSULTANT, $125/h equivalent. Being the PRIME is the only way SBIR is a major line.
- Core-facility-style hourly access: the market rate is set by the universities themselves ($43–64/h unassisted external, $147–192/h assisted; Notre Dame, UF above). A private lab cannot undercut a subsidised core and cannot exceed it by much for raw hours.
- Revenue-mix fraction for small testing labs (industry vs academic vs government): **not found**. IBISWorld sizes US "Laboratory Testing Services" at $33.6 bn (2025, +2.4 %) across ~12,038 businesses (2026), i.e. a mean of ~$2.8 M/business, skewed by large firms — https://www.ibisworld.com/united-states/market-size/laboratory-testing-services/1408/ ; https://www.ibisworld.com/united-states/number-of-businesses/laboratory-testing-services/1408/

---

## Rate-card table

| Service | Provider | Price | Source | Date read |
|---|---|---|---|---|
| Rheology/thermal instrument time | ASM Services (UK) | £90/h + £200–300 fixed (calibration, report, analysis) | https://www.asm-serv.co.uk/contract-testing | 2026-09-06 |
| Rheometer, external unassisted | Notre Dame MCF (HR-2) | $43/h (rate sheet) / $45/h (instrument page) | https://mcf.nd.edu/internal-users/rate-sheet/ | 2026-09-06 |
| Rheometer, external assisted | Notre Dame MCF (HR-2) | $147/h or $258/day (rate sheet); $192/h (instrument page) | https://mcf.nd.edu/capabilities-instruments/rheometry-dmta-tribology/hybrid-rheometer/ | 2026-09-06 |
| Rheometer, internal | Notre Dame MCF | $26–29/h unassisted; $77–133/h assisted; $156/day | same | 2026-09-06 |
| Rheometer ARES LS1 | U. Florida NRF | $16.20/h academic; $47.90/h affiliated commercial; $63.90/h commercial [snippet] | https://nrf.aux.eng.ufl.edu/ccb/resource.asp?id=105 | 2026-09-06 |
| Rheology report, turnaround | Rheology Testing Services (NC) | quote only; ~1 week to draft | https://www.rheologytestingservices.com/applications-and-support-overview | 2026-09-06 |
| Rheology | Intertek, Measurlabs, Centre for Industrial Rheology, Rheology Solutions | quote only, none public | see Q1 table | 2026-09-06 |
| Spray/plume/droplet imaging | Oxford Lasers, Proveris, Spraying Systems, Silsoe SSAU, ImageXpert | none public | see Q2 table | 2026-09-06 |
| ISO 11608 / 11040 injection force | Eurofins MDS, Smithers, DDL | none public | see Q3 table | 2026-09-06 |
| High-speed camera rental (shadow price) | Phantom rental houses | ~$2,400–3,000+/day | https://beverlyboy.com/how-much-does-it-cost/how-much-does-a-phantom-camera-cost/ | 2026-09-06 |
| SBIR consultant cap | NSF | $1,000/day (8 h) | https://seedfund.nsf.gov/solicitation-budget/ | 2026-09-07 |
| Lab space | Purdue Technology Center, Merrillville | $38.09/RSF/yr [snippet] | https://www.commercialsearch.com/commercial-property/us/in/merrillville/purdue-technology-center/ | 2026-09-06 |
| Lab space | 16 Tech Indianapolis, build-to-suit | ≥$40/sq ft/yr [snippet] | https://www.ibj.com/articles/16-tech-pencils-in-plans-for-next-two-buildings-after-early-success | 2026-09-06 |
| Maker/prototype membership | mHUB Chicago | $335/mo co-working; $450/mo desk; $1,350+/mo office | https://mhub.org/become-a-member | 2026-09-06 |

## Equipment cost table

| Item | New range | Used range | Source | Date read |
|---|---|---|---|---|
| Rotational rheometer, Anton Paar MCR 302 | quote only (not found) | €8,000 (Altamira, LabX) – €17,000 (Zeetree, LabX) [snippet] | https://www.labx.com/item/anton-paar-mcr-302-rheometer/scp-236757-5394916 ; https://www.labx.com/item/anton-paar-mcr-302/scp-233109-5250273 | 2026-09-06 |
| Rotational rheometer, TA DHR (model unstated) | quote only | $39,999.99 incl. shipping, install, warranty, training (Liberty Instruments, LabX) [snippet] | https://www.labx.com/item/ta-instruments-dhr-rheometer-free-shipping-install-warranty-training/scp-220751-262ed154-1789-4fd8-afac-5beef71b500f | 2026-09-06 |
| Rotational rheometer, TA HR-2 | — | Exapro listing "price on request", sold, as-is | https://www.exapro.com/ta-instruments-discovery-hr-2-p250214033/ | 2026-09-06 |
| Entry rheometer, Anton Paar MCR 72/92; Kinexus; DHR-1 | not found (quote only) | not found | https://www.anton-paar.com/us-en/products/details/rheometer-mcr-72-and-mcr-92/ | 2026-09-06 |
| Brookfield DV2T viscometer | listed by Fisher/CMTP, price not rendered — not found | not found | https://certifiedmtp.com/brookfield-dv2t-viscometer-lv-torque-threaded-coupling-w-point-jewel/ | 2026-09-07 |
| Glass capillary viscometer, Cannon-Fenske routine (calibrated) | $221.15–223.45 (Cole-Parmer); $521 (YSI) | — | https://www.coleparmer.com/p/cannon-fenske-calibrated-glass-capillary-routine-viscometers/7145 ; https://www.ysi.com/product/id-285403753/cannon-fenske-routine-viscometer | 2026-09-07 |
| Glass capillary viscometer, Cannon-Ubbelohde (calibrated) | $314.60–469.70 (Cole-Parmer/VWR) | — | https://www.coleparmer.com/p/cannon-ubbelohde-calibrated-glass-capillary-viscometers/7143 | 2026-09-07 |
| Pressure capillary rheometer (Dynisco-class) | not found | not found | — | — |
| High-speed camera, Kron Chronos 2.1-HD (1080p @ 1000 fps, 24 kfps reduced) | $4,995 (8 GB); $5,400 (16 GB); $6,000 (32 GB) | — | https://www.newsshooter.com/2021/01/09/krontech-chronos-2-1-hd-1000fps-for-4995-usd/ ; https://www.krontech.ca/product/chronos-2-1-hd-high-speed-camera/ | 2026-09-07 |
| High-speed camera, Kron Chronos 4K12 | $14,495 (64 GB); $15,955 (128 GB) | — | https://petapixel.com/2023/11/15/the-new-affordable-chronos-4k12-and-q12-shoot-at-nearly-30000-fps/ | 2026-09-06 |
| High-speed camera, Phantom T1340 | ~$80,000+ | — | https://www.redsharknews.com/phantom-t1340-a-potentially-lower-cost-way-of-shooting-high-frame-rates | 2026-09-06 |
| High-speed camera, Phantom Miro M320S | — | $13,250 complete kit (eBay) [snippet] | https://kitplus.com/for-sale/accessories-35/phantom/high-speed-miro-m320s-manaudio891-at-gmailcom/32989 ; https://www.ebay.com/shop/phantom-vision-research?_nkw=phantom+vision+research | 2026-09-07 |
| High-speed camera, Phantom v711 | — | eBay listing, price not shown in snippet (403 on fetch) | https://www.ebay.com/itm/406290245513 | 2026-09-06 |
| High-speed camera, Photron Mini AX / Nova | not found (quote only) | not found | https://photron.com/fastcam-mini-ax/ | 2026-09-06 |
| Universal tester, Instron 5944 (2 kN) | quote only | $9,200 (eBay) [snippet] | https://www.ebay.com/itm/205608131725 | 2026-09-07 |
| Motorised test stand, Mark-10 ESM303 (1.5 kN) | $1,420.25 (QC Solutions) | $999 bundle with M5-50 gauge (eBay) [snippet] | https://www.qc-solutions.com/product/mark-10-esm303-motorized-test-stand-force-gauge ; https://www.ebay.com/itm/256929030572 | 2026-09-06 |
| Force gauge, Mark-10 M5-50 (250 N) | $1,670 (Radwell new) | $599.67 refurbished (Radwell) | https://www.radwell.com/Buy/MARK%2010/MARK%2010/M5-50 | 2026-09-06 |
| Mecmesin MultiTest-i / dV | not found (quote only) | — | https://www.mecmesin.com/software-controlled-force-systems/single-column-force-tester | 2026-09-06 |
| Syringe pump, New Era NE-1000 | $932.25 (Fisher) | eBay listing exists, price not captured | https://www.fishersci.com/shop/products/syringe-pump-ne-1000/NC9355767 | 2026-09-06 |
| Syringe pump, Chemyx Fusion 200-X | $2,594.82 (Fisher) | $450–519 (Fusion 100/200, LabX/eBay) [snippet] | https://www.fishersci.com/shop/products/fusion-200-x-syringe-pump/NC2235109 ; https://www.labx.com/item/chemyx-fusion-100-syringe-pump/DIS-87999-3373978 | 2026-09-06 |
| Syringe pump, Harvard PHD Ultra | $21,646.99 on a medical reseller — inconsistent, treat as unreliable | $1,300–1,500 (LabX/Lab World Group) [snippet] | https://www.graylinemedical.com/products/harvard-apparatus-standard-infusion-phd-ultra-syringe-pumps-pump-syringe-infusion-phd-ultra-703005 ; https://www.labx.com/item/harvard-apparatus-phd-ultra-syringe-pump/DIS-50372-412808 | 2026-09-06 |
| Lighting / optics / lenses | not researched | — | — | — |

**Indicative totals (arithmetic on the rows above, not a quote):** used MCR 302 (~$9–19k at current EUR/USD, unknown accessories) + Chronos 2.1-HD ($5–6k) + Mark-10 stand and gauge ($2–3k) + NE-1000 pump ($0.9k) + two capillary viscometers ($0.5k) ≈ **$18k–30k** for a bare bench with a consumer-grade high-speed camera; substitute a used Phantom Miro (~$13k) and a used DHR with warranty (~$40k) and the total is ≈ **$55k–70k**; a new Phantom T1340 alone is ~$80k. Accessories (Peltier plate, cone/plate geometries, solvent trap) are not priced here and on a used rheometer are often the deciding line.

---

## Is a bench line viable for a 1–3 person shop, and in which form?

**As a stand-alone service: no.** The price ceiling for raw rheometer or camera time is set by subsidised university cores at $43–64/h unassisted and $147–192/h assisted (Notre Dame, UF), and the one public independent rate is £90/h plus a £200–300 report fee. At 1,000 billable instrument-hours a year (an optimistic utilisation for one person also doing sales and reporting) that is $50k–190k gross, before the rheometer's service contract, calibration fluids, space at ~$1,900/mo, and insurance. Regulated segments that command real money (FDA nasal-spray CMC, ISO 11608 release testing) require GMP / ISO 17025 accreditation and are held by Proveris, Eurofins, Smithers, DDL.

**As the validation arm of a simulation product: yes, in three specific forms.**
1. **Parameter-and-uncertainty packages.** Every rate card found ends at "plots, tables, report". None ends at "Herschel-Bulkley (k, n, tau_y) with profile-likelihood intervals, identifiability flags, exported as an OpenFOAM `transportProperties` block". The founder's `biorheo`/`dst_fit` code already produces that; the bench exists to generate the flow curves it needs and to charge for the package, not the hours.
2. **Device-level validation data with audited calibration.** Needle force vs flow rate, drop/jet imaging with px/mm and fps traceability — sold as the evidence that a simulation is right, priced per study, delivered alongside the simulation. This is the Veryst shape (test → calibrated model → simulation) scaled down; Veryst's existence is the only market evidence for it, and Veryst publishes no numbers.
3. **Pre-accreditation R&D testing for device developers who are not yet at release testing.** The ISO 11608 labs found do batch release under GMP; a startup at prototype stage wanting "will this formulation be injectable through 27G at 10 s?" has no cheap answer today. Simulants-only keeps it out of 1910.1030 and IRB.

**Form to avoid:** hourly instrument access, GMP/ISO 17025 release testing, and anything that requires human blood or plasma.

## Evidence AGAINST

- No public evidence that any buyer pays more for combined sim + bench than for the two separately. Veryst is one firm, and its testing revenue share is unknown.
- No small shop found that won on having both. The search "small consultancy CFD experimental validation in-house lab" returned only academic papers and one rheology lab's article.
- University cores make raw instrument time nearly free to academics ($16–29/h) and cheap to industry ($43–64/h); academic customers are therefore not a bench revenue source.
- Contract rheology's incumbents include very small operations (Rheology Testing Services in NC has a gmail contact and one-week turnaround) — meaning the segment is already served at the one-person scale, and they survive on quote-only pricing, which suggests thin margins rather than a gap.
- SBIR subaward math caps a subcontractor at ≤1/3 of a Phase I; consultant cap $1,000/day.
- Used rheometer prices in the snippets (€8,000) are almost certainly bare heads; a Peltier plate and geometries can double the number, and no accessory prices were captured.
- The best-priced imaging camera (Chronos 2.1-HD, $5k) is 1000 fps at 1080p; splash and jet break-up work in this lab has used far higher frame rates, and a research-grade camera is $13k used to $80k new.

## Open questions

1. Actual quotes: request three (ASM, Rheology Testing Services, Centre for Industrial Rheology) for a 5-sample flow-curve + amplitude-sweep job to convert "quote only" into numbers.
2. Open the LabX/eBay listings in a browser: what accessories come with the €8,000 and €17,000 MCR 302s and the $40k DHR?
3. Illinois Research Park EnterpriseWorks lab rent (email in Q5) and Purdue West Lafayette wet-lab rate — both on request.
4. Nelson Labs, West, BD, Stevanato service arms: not reached in this pass.
5. Does any ISO 11608 lab sell pre-release R&D syringeability studies at a lower tier? Not found.
6. Veryst revenue split and headcount (private; try Dun & Bradstreet or state filings).
7. Animal-blood (porcine/bovine) handling rules if a simulant is ever insufficient — not researched.
8. Lighting, optics, Peltier, geometries, calibration oils, rheometer service contract — none priced.
9. Fraction of small testing-lab revenue from academia/government — no survey found; IBISWorld report (paywalled) may have a buyer segmentation.

## Sources

Rheology providers and cores
- https://www.asm-serv.co.uk/contract-testing
- https://www.rheologytestingservices.com/applications-and-support-overview
- https://www.rheologytestingservices.com/
- http://www.rheologysolutions.com/company/applications-laboratory/
- http://www.rheologysolutions.com/resources/industry-contract-testing-kits/
- https://www.rheologylab.com/services/viscosity-testing-lab/
- https://www.rheologylab.com/articles/modelling-ready-rheology-for-computational-fluid-dynamics/
- https://measurlabs.com/methods/rheology-testing/
- https://www.intertek.com/analytical-laboratories/rotational-rheology/
- https://www.intertek.com/analytical-laboratories/allentown/rheology-testing/
- https://www.jkanalyticallab.com/rheology-viscosity-testing-services-6018325.html
- https://rheologylabs.com/
- https://www.anton-paar.com/us-en/rheology-instruments/
- https://mcf.nd.edu/internal-users/rate-sheet/
- https://mcf.nd.edu/capabilities-instruments/rheometry-dmta-tribology/hybrid-rheometer/
- https://nrf.aux.eng.ufl.edu/ccb/resource.asp?id=105
- https://mrl.illinois.edu/facilities/equipment/ta-dhr-3-rheometer
- https://mrsec.uchicago.edu/facilities/materials-preparation-and-measurement-laboratory/rates/
- https://www.chem.purdue.edu/ric/rates.html

Imaging / spray
- https://oxfordlasers.com/services/laser-imaging-services/
- https://proverislabs.com/
- https://www.proveris.com/about-us/
- https://www.cambridgenetwork.co.uk/news/melbourn-first-offer-fda-approved-spray-pattern-and-plume-geometry-analysis
- https://www.spray.com/services/advanced-spray-testing-services
- https://www.sprayanalysis.com/testing-services/instruments-and-facilities/high-speed-imaging
- https://croplifeeurope.eu/wp-content/uploads/2021/01/S0252-contract-report-final-complete.pdf
- https://ukagritechcentre.com/community-directory/silsoe-spray-applications-unit-ltd/
- https://imagexpert.com/process-development/
- https://imagexpert.com/drop-analysis/dropwatcher/
- https://beverlyboy.com/how-much-does-it-cost/how-much-does-a-phantom-camera-cost/
- https://www.atecorp.com/ndt-inspection/visual-inspection/high-speed-cameras

Injection force
- https://www.eurofins.com/medical-device/testing/combination-products/syringe-testing/
- https://www.smithers.com/services/testing/standard/iso/iso-11608-1
- https://www.ddltesting.com/medical-device-testing/needle-based-injection-system-testing/
- https://www.instron.com/en/testing-solutions/iso-standards/iso-11040/
- https://www.zwickroell.com/industries/medicalpharmaceutical/therapy-systems/syringe-testing-iso-11040/
- https://www.mecmesin.com/publications/break-loose-and-sustaining-force-tests-pre-filled-drug-cartridges-iso-11608-3

Equipment
- https://www.labx.com/item/anton-paar-mcr-302-rheometer/scp-236757-5394916
- https://www.labx.com/item/anton-paar-mcr-302/scp-233109-5250273
- https://www.labx.com/item/ta-instruments-dhr-rheometer-free-shipping-install-warranty-training/scp-220751-262ed154-1789-4fd8-afac-5beef71b500f
- https://www.exapro.com/ta-instruments-discovery-hr-2-p250214033/
- https://www.anton-paar.com/us-en/products/details/rheometer-mcr-72-and-mcr-92/
- https://www.newsshooter.com/2021/01/09/krontech-chronos-2-1-hd-1000fps-for-4995-usd/
- https://www.krontech.ca/product/chronos-2-1-hd-high-speed-camera/
- https://petapixel.com/2023/11/15/the-new-affordable-chronos-4k12-and-q12-shoot-at-nearly-30000-fps/
- https://www.redsharknews.com/phantom-t1340-a-potentially-lower-cost-way-of-shooting-high-frame-rates
- https://kitplus.com/for-sale/accessories-35/phantom/high-speed-miro-m320s-manaudio891-at-gmailcom/32989
- https://www.ebay.com/itm/406290245513
- https://www.ebay.com/itm/205608131725
- https://www.qc-solutions.com/product/mark-10-esm303-motorized-test-stand-force-gauge
- https://www.ebay.com/itm/256929030572
- https://www.radwell.com/Buy/MARK%2010/MARK%2010/M5-50
- https://www.fishersci.com/shop/products/syringe-pump-ne-1000/NC9355767
- https://www.fishersci.com/shop/products/fusion-200-x-syringe-pump/NC2235109
- https://www.labx.com/item/harvard-apparatus-phd-ultra-syringe-pump/DIS-50372-412808
- https://www.coleparmer.com/p/cannon-fenske-calibrated-glass-capillary-routine-viscometers/7145
- https://www.coleparmer.com/p/cannon-ubbelohde-calibrated-glass-capillary-viscometers/7143
- https://www.ysi.com/product/id-285403753/cannon-fenske-routine-viscometer

Space and compliance
- https://researchpark.illinois.edu/locate-here/lab-spaces/
- https://www.commercialsearch.com/commercial-property/us/in/merrillville/purdue-technology-center/
- https://purdueresearchpark.com/locations/west-lafayette/
- https://www.ibj.com/articles/16-tech-pencils-in-plans-for-next-two-buildings-after-early-success
- https://iulab.iu.edu/bio-start-up-center/index.html
- https://mhub.org/become-a-member
- https://www.excedr.com/blog/how-much-does-lab-space-cost
- https://www.osha.gov/bloodborne-pathogens/general
- https://www.osha.gov/laws-regs/standardinterpretations/1993-02-01-0
- https://www.osha.gov/laws-regs/standardinterpretations/2009-02-27
- https://www.osha.gov/laws-regs/standardinterpretations/2002-12-27

Sim + bench, grants, market size
- https://www.comsol.com/blogs/veryst-combines-material-testing-and-simulation-for-reliable-results
- https://www.veryst.com/services/testing/material-model-calibration
- https://www.veryst.com/news-events/veryst-software-spun-new-company-polymerfem
- https://seedfund.nsf.gov/solicitation-budget/
- https://www.ibisworld.com/united-states/market-size/laboratory-testing-services/1408/
- https://www.ibisworld.com/united-states/number-of-businesses/laboratory-testing-services/1408/
