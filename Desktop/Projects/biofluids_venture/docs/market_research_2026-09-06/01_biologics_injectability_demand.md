# Demand for injectability / drug-product fluid-handling expertise in biopharma

Market-research dimension 01. Compiled 2026-09-06 from ~40 web searches and page fetches.
Every factual claim carries a URL. "Vendor" marks vendor marketing. "Not found" means the
search was run and returned nothing usable, not that the thing does not exist. Several
primary pages (ScienceDirect, ONdrugDelivery, Kymanox, Regeneron careers, Indeed) returned
HTTP 403/401 to the fetcher; where only a search-engine snippet was available, that is said.

---

## Summary

1. The pull is real and growing: SC volumes are moving from 1-2 mL to "up to 10 and even 25 mL" (Allmendinger/Mahler/Bittner, ADDR, July 2026), Halozyme's SC-enabling royalties hit $868M in 2025 (+52 %) and it paid up to $400M for Surf Bio's 500 mg/mL hyperconcentration tech in Dec 2025.
2. Money flows to **enablers** (excipients, devices, viscometers), not to injectability *analysis* as a line item; no public price for an injectability study, ISO 11608 campaign or CFD study was found.
3. Sellers of *measurement* are plentiful (West ALS, Smithers, Cormica, Coriolis, Lonza, Stevanato, WuXi); sellers of *prediction with stated error and non-Newtonian rheology* are few (Veryst: COMSOL tissue model, no non-Newtonian, no validation on page; QATCH visQ.AI: ML, "<=20 % error", pilot stage).
4. Practitioner-stated gaps match the founder's skills: viscosity models fail on extrapolation from low to ultra-high concentration (-186 % to +177 %, GSK/Strathclyde 2024); injection-force models drift above 230 mg/mL and in 27G/30G needles (JHU thesis); cold-chain viscosity (43 cP at 5 C vs 24 cP at 15 C, vendor); needle shear ~1e5 1/s; tissue back-pressure 24 kPa at 1 mL/10 s in minipig.
5. Regulatory pull exists but is indirect: FDA's June 2024 EDDO draft guidance and the 2013 pen/jet guidance make injection force, time and dose accuracy verified design outputs; ISO 11608-1:2022 and ISO 11040-8 define the tests. This funds test labs, not modellers, unless a model shortens verification.
6. Job-posting evidence is weak: formulation roles cite "viscosity, self-association ... injectability" as qualifications, but no posting titled or centred on injectability modelling / CFD for injection was found; career pages blocked the fetcher.
7. Conference signal is strong and consistent 2025-26 (PDA UPFS Vienna session on high-viscosity large-volume SC moderated by Lonza; PODD 2025 with 28 named pharma; SHL/Stevanato 50 cP, 2.25 mL data, Jan 2026).
8. Evidence against: the buyer for analysis is a pharma device/formulation group that already owns Instron/Zwick rigs and a Hagen-Poiseuille spreadsheet; consultancies (Veryst, Springboard, Team) already sell modelling; the physics is 80 % solved by a 2018 Roche model with <12 % error.
9. Credible wedge for 1-3 people: rheology-to-needle extrapolation with uncertainty and identifiability, cold/warm injection-time envelopes, and shear/aggregation risk screening, sold as a fixed-price study to formulation groups without in-house rheologists (mid-size biotech, CDMOs).
10. Must be settled by conversation, not search: what a device group actually pays for an external injectability assessment, and whether a model output is admissible in a design-verification file.

---

## Evidence by question

### Q1. Who currently sells injectability testing / modelling?

**Measurement services (someone sells this).**

| Seller | What they say they offer | Non-Newtonian / prediction / uncertainty? | Source |
|---|---|---|---|
| West Pharmaceutical, Analytical Laboratory Services | Automated autoinjector testing to ISO 11608-5, "independent analytical laboratory", tests other manufacturers' devices | Not mentioned | https://www.westpharma.com/blog/2023/january/auto-injector-iso-11608-5-regulatory-requirements-analytical-services (2023, vendor) |
| Smithers | Cartridge (ISO 11608), hypodermic needle, pen/auto-injector, PFS testing | Not mentioned | https://www.smithers.com/industries/medical-device/products-testing/injectables/cartridge-testing (vendor) |
| Cormica | "GMP, GLP, and ISO 17025 accredited laboratory with specific accreditation for ISO 11608-1:2022" | Not mentioned | https://www.cormica.com/services/physical-and-packaging-testing/physical-mechanical-testing-for-medical-devices/iso-11608-needle-based-injection-system-testing/ (vendor) |
| Coriolis Pharma (CRO, Germany) | Syringeability: break-loose and gliding force in representative PFS; capillary viscometry | **Yes, explicitly**: "biopharmaceutical products rarely show Newtonian flow behavior, their rheological properties (governing the injection process) depend on the force applied." "In Silico Modeling" listed as a separate service; no detail, no price | https://www.coriolis-pharma.com/our-services/injectability-amp-syringeability/ (vendor, fetched 2026-09-06) |
| Lonza, Biologics Special Services | "drug and device functionality testing (such as injectability)" with "access to the required equipment" | Not mentioned | https://www.lonza.com/integrated-biologics/parenteral-drug-product/biologics-special-services (vendor) |
| WuXi Biologics, WuXiHigh 2.0 | High-concentration formulation platform "up to 230mg/mL", "over 24 proprietary excipient combinations to reduce viscosity by up to 90%" | Formulation platform, not analysis | https://www.wuxibiologics.com/drug-product-development-services/ and https://pharmasource.global/content/news/cdmo-news/wuxi-biologics-launches-high-concentration-formulation-platform-wuxihigh2-0/ (vendor) |
| Stevanato Group analytical services | Container selection, drug compatibility, labs in Italy and Boston | Not mentioned | https://www.stevanatogroup.com/en/analytical-testing-services/analytical-services/ (vendor) |
| ZwickRoell, Instron | Test *instruments* and fixtures for ISO 11608 / ISO 11040 (cap removal, activation force, delivery time, delivered mass, needle guard) | Instrument vendors, not services | https://www.zwickroell.com/industries/medicalpharmaceutical/therapy-systems/autoinjectors-iso-11608-5/ ; https://www.instron.com/en/testing-solutions/iso-standards/iso-11608/ ; https://www.ondrugdelivery.com/navigating-iso-11608-with-automated-equipment-for-customised-pen-and-autoinjector-testing/ |

**Modelling / simulation services.**

| Seller | Offer | Gap visible on the page | Source |
|---|---|---|---|
| Veryst Engineering | COMSOL multiphysics model of SC injection: "injection force, delivery time at specified flow rates, surface rise ('bleb' formation), and three-dimensional drug bolus dispersion"; poroelastic tissue | Non-Newtonian drug rheology **not mentioned**; experimental validation and uncertainty quantification **not mentioned**; no clients, no price | https://www.veryst.com/case-studies/optimizing-subcutaneous-drug-delivery-multiphysics-simulation (fetched 2026-09-06, vendor) ; intravitreal sibling https://www.veryst.com/case-studies/optimizing-intravitreal-drug-delivery-multiphysics-simulation |
| Springboard (UK device consultancy) | "concept design, feasibility and modeling, Usability and Human Factors Engineering, and performance testing" for autoinjectors; Gx Inbeneo case study | No rheology-specific claim found | https://www.springboard.pro/case-studies/gx-inbeneo-autoinjector-design/ ; https://www.springboard.pro/autoinjectors-vs-on-body-systems-the-future-of-large-volume-delivery/ (vendor) |
| Team Consulting | Attends PDA UPFS 2026; device consultancy | Page fetched only as event listing | https://www.team-consulting.com/us/insights/pda-universe-of-pre-filled-syringes-and-injection-devices-conference/ |
| QATCH Technologies, visQ.AI | ML prediction of viscosity/injectability trained on nanovisQ data; claims "100 - 15,000,000 s-1 shear-rate range", "R2 > 0.85", "<= 20 % error" on untested formulations, 4 uL; "pilot collaborations" only | Instrument-tied; accuracy claim is vendor; no physics of the device | https://qatchtech.com/visq-ai (fetched 2026-09-06, vendor) |
| RheoSense (instrument vendor) | Blog series on injection-force prediction from VROC viscosity; "understanding the viscosity at the shear rate inside your needle is crucial" | Newtonian/power-law formulas; no uncertainty; instrument marketing | https://blog.rheosense.com/predicting-injection-force-theoretical-foundation (2021-03-09, vendor) ; https://www.rheosense.com/m-vrocii/application-brief-auto-injector-high-viscosity-testing |
| Roche (in-house, published) | Physics-based autoinjector injection-time model, "modeling error for all investigated device configurations was smaller than 12% and the error range was less than 6%"; used for statistical simulation over component variability | Shows the incumbent buyer can build this in-house | https://www.tandfonline.com/doi/full/10.2147/MDER.S151727 (2018) |
| SHL Medical + Stevanato (device makers) | Jan 2026 joint study: 1-50 cP, 2.0-2.25 mL, injection time cut from ">25 seconds to 8-12 seconds", glide force "<15N", 27G 8 mm; "Hydrodynamic pressure drop modeling using Hagen-Poiseuille equation" | Newtonian model only, on the vendors' own page | https://www.shl-medical.com/news-insights/blog/high-viscosity-delivery-molly-autoinjector (2026-01-05, vendor) |

**What is missing across all sellers:** no page found that offers (a) rheology measured at or extrapolated to needle shear rates with a stated model and confidence band, (b) identifiability of the fitted constitutive parameters, (c) a validated non-Newtonian injection-force/time prediction with reported error against a hold-out, or (d) temperature-dependent (2-25 C) injection envelopes as a deliverable. Exponent, Thornton Tomasetti, Nemera, BD, Ypsomed service arms: **not found** offering injectability *analysis* as a service (BD/Ypsomed appear as device collaborators; see Q7).

### Q2. Job postings 2025-2026

Weak evidence; career sites blocked the fetcher.

- Indeed aggregate: "646 pharmaceutical biologics formulation scientist jobs" (search-snippet count for a broad query, not injectability-specific). Snippet text from one posting on that results page: "deep understanding of challenges at elevated protein concentrations such as viscosity, self-association, aggregation, and familiarity with subcutaneous drug delivery fundamentals including injection volume/concentration constraints and injectability" and "proficiency with rheology". https://www.indeed.com/q-Pharmaceutical-Biologics-Formulation-Scientist-jobs.html and https://www.indeed.com/q-injectable-formulation-scientist-jobs.html (snippets only; the specific employer was not resolvable; direct fetch returned HTTP 401).
- AbbVie, Senior Scientist I, Biologics Drug Product Development (South San Francisco): fetched; **no** mention of injectability, viscosity, rheology, device or modelling. https://careers.abbvie.com/en/job/senior-scientist-i-biologics-drug-product-development-in-south-san-francisco-ca-jid-10300
- AbbVie, Senior Device Development Engineer (North Chicago): "Late-Stage Injection Device Development team with end-to-end accountability for commercial development, submission, industrialization and marketed product support for injection devices" (snippet). https://careers.abbvie.com/en/job/senior-device-development-engineer-in-north-chicago-il-jid-26583
- Regeneron, Sr Scientist Formulations (Tarrytown) and two sibling postings exist; page returned 403, content not verified. https://careers.regeneron.com/en/jobs/r43812/sr-scientist-formulations-drug-product-development/
- Eurofins, Biologics Formulation Scientist (Temporary), LinkedIn listing. https://www.linkedin.com/jobs/view/biologics-formulation-scientist-temporary-at-eurofins-4448630800 (not fetched).
- CFD-for-injection roles: a search returned only inhalation CFD ("predict aerosol deposition ... ANSYS Fluent") and generic device roles. https://www.indeed.com/q-Biomedical-Engineer-Drug-Delivery-System-jobs.html
- `site:linkedin.com/jobs "injectability"`: returned no relevant listings.

**Count:** not obtainable from public search. No posting was found whose title contains injectability, injection modelling, or CFD for injection devices.

### Q3. Conference and industry-body signals

- **PDA Universe of Pre-filled Syringes and Injection Devices 2025** (Vienna, 21-22 Oct 2025): a session with "three experts exploring the complex interplay between formulation science, delivery mechanics, and patient experience in the context of high-viscosity, large-volume subcutaneous administration", moderated by Flora Felsovalyi (Lonza); talks "Exploring Opportunities for Large Volume Subcutaneous (LVSC) Injections for Acute Illness Treatment" (Simon Moss, Pfizer) and "Enabling High-Viscosity Biologic Delivery with SHL Medical's Molly 2.25 mL Autoinjector and Stevanato Group's Syringes" (Barichello, Stevanato; Fetz, SHL). Call for abstracts lists "High Concentration, High Viscosity, Large Volume, Shear Sensitive". https://www.pda.org/rest-api/event-agenda/pda-universe-of-pre-filled-syringes-and-injection-devices-conference-2025 (agenda API; 403 on direct fetch, content from search snippet) ; https://www.pda.org/global-event-calendar/abstract-submissions/call-for-abstracts/pda-universe-of-pre-filled-syringes-and-injection-devices-conference-2025
- **PDA UPFS 2026** (Palm Springs, 6-7 Oct 2026): topic list is human factors, manufacturing, lifecycle; high viscosity not in the ten headline topics on the third-party summary. https://intuitionlabs.ai/conferences/pda-universe-of-pre-filled-syringes-injection-devices ; https://www.pda.org/global-event-calendar/event-detail/pda-universe-of-pre-filled-syringes-and-injection-devices-conference-2026
- **PODD 2025** (Boston, 27-28 Oct 2025): 28 named pharma presenting (AbbVie, Amgen, ..., UCB); agenda covers high viscosity and large-volume SC. https://www.pharmtech.com/view/podd-2025-unveils-agenda-highlighting-innovation-in-drug-delivery-technologies ; https://www.poddconference.com/
- **AAPS PharmSci 360 2025** (San Antonio, 9-12 Nov 2025): poster on a "low-volume, broad-shear viscometry platform" and a poster on HP-beta-CD as viscosity-reducing excipient. https://aaps2025.eventscribe.net/fsPopup.asp?PosterID=764032&mode=posterInfo
- **PDA 2024 poster** (Wu, Li, Wang, Zhuang, Chen, Guo): injection-force model with Newtonian and shear-thinning viscosity and injection temperature; published EJPB 2024, art. 114221. https://www.pda.org/docs/default-source/posters/universe-of-prefilled-syringes-and-injection-devices/chen-wei.pdf ; https://www.sciencedirect.com/science/article/abs/pii/S093964112400047X
- **Stevanato webinar** "Overcoming Volume and Viscosity Barriers in Biologic Formulations". https://www.stevanatogroup.com/en/news-events/webinars/overcoming-volume-and-viscosity-barriers-in-biologic-formulations/ (vendor)
- **ONdrugDelivery "Trends for 2026"** exists but returned 403. https://www.ondrugdelivery.com/trends-for-2026/
- IPAC-RS: not searched (inhalation body; out of scope for this dimension). USP <1>: not found relevant.

Hot topics 2025-26, consistently: large-volume SC (2-10+ mL), >15-50 cP formulations, injection time >20-30 s as the usability barrier, on-body injectors vs 2.25 mL autoinjectors, ultra-high concentration (>400 mg/mL) suspensions.

### Q4. Price points

- Injectability study price, ISO 11608 campaign price, CFD study price: **not found**. No rate card from any lab above. A generic bioanalytical-CRO article notes rate cards hide "$85,000-180,000 in method work" before study samples (different service; context only). https://biobridgex.com/costs/bioanalytical-testing-costs
- Standard document price: ISO 11608-1:2022, CHF 225. https://www.iso.org/standard/70733.html
- Public grant sizes for injectability *measurement* tech (QATCH Technologies, Duke spin-out, RTP NC):
  - NIH SBIR Phase II, "Developing prototype injectability and developability testing system based on microfluidic quartz resonators", **$1,500,001**, 2022. https://www.sbir.gov/awards/200235
  - NIH SBIR Direct Phase II, high-throughput nanovisQ, **$1,883,069**, 2023-09-15 to 2025-09-14. https://ars-pharma.inknowvation.com/sbir/awards/nih-2023-high-throughput-injectability-screening-high-concentration-protein-formulations
  - Six SBIR/STTR awards totalling "approximately $3M"; NCBiotech $250,000 (2023). https://sbtdc.org/news/nc-sbir-success-story-qatch-technologies ; https://www.ncbiotech.org/news/ncbiotech-awards-22-million-grants-loans-latest-quarter
- NIIMBL: ~$10M across 8 projects (Apr 2024) and $15.8M across 14 projects; a high-concentration objective ">=150 grams/liter" appears, but no award specifically for injectability modelling was found. https://www.prnewswire.com/news-releases/niimbl-announces-10m-for-8-new-biopharmaceutical-manufacturing-projects-302112132.html ; https://www.niimbl.org/news/niimbl-announces-15-8m-to-fund-14-new-biopharmaceutical-manufacturing-projects/

### Q5. Regulatory pull

- **FDA draft guidance, "Essential Drug Delivery Outputs for Devices Intended to Deliver Drugs and Biological Products"** (June 2024; Federal Register 2024-07-01; comments to 2024-09-30). EDDO = "design outputs necessary to ensure delivery of the intended drug dose to the intended delivery site", including "initiation, progression, completion of the dose". Appendix gives example EDDOs for PFS, autoinjectors, on-body injectors, pen injectors. West: "all verification and validation activities for EDDOs should be conducted after stress conditioning the product". https://www.federalregister.gov/documents/2024/07/01/2024-14409/essential-drug-delivery-outputs-for-devices-intended-to-deliver-drugs-and-biological-products-draft ; https://www.thefdalawblog.com/2024/07/fda-releases-draft-guidance-on-essential-drug-delivery-outputs/ ; https://www.westpharma.com/blog/2024/september/fda-guidance-essential-drug-delivery-output-requirements ; https://www.kymanox.com/thought-leadership-series/blog/fda-guidance-essential-drug-delivery-outputs-kymanox-what-you-need-to-know/
- **FDA guidance (final, June 2013), "Technical Considerations for Pen, Jet, and Related Injectors"**: dose accuracy per ISO 11608-1 at 23 +/- 5 C; breakaway/glide force; notes viscous formulations demand more force and stronger drive. https://www.fda.gov/files/about%20fda/published/Technical-Considerations-for-Pen--Jet--and-Related-Injectors-Intended-for-Use-with-Drugs-and-Biological-Products.pdf ; https://www.federalregister.gov/documents/2013/06/07/2013-13484/
- **FDA "Glass Syringes ... Technical Information to Supplement ISO 11040-4"**. https://www.fda.gov/media/85748/download (PDF fetched but not text-extractable here; contents unverified)
- **ISO 11608-1:2022** (dose accuracy, now with dose-profile/time evaluation) and **ISO 11608-5:2022** (autoinjectors). https://www.iso.org/standard/70733.html ; https://cdn.standards.iteh.ai/samples/76627/440310c2283e45dea7ca5908c4ae7975/ISO-11608-5-2022.pdf
- **ISO 11040-8** (finished PFS as combination product: break-loose/extrusion force, penetration force, leakage, burst, administration time). https://www.instron.com/en/testing-solutions/iso-standards/iso-11040/ ; https://www.instron.com/en/resources/blog/2026/june/iso-11040-8-2026-pre-filled-syringe-testing-guide/
- BioProcess International on PFS EPRs: "Performance testing should take into consideration the viscosity of the drug or biological product, resistance created by the fluid contact path, and the back pressure/shear force generated during the injection." https://www.bioprocessintl.com/fill-finish/essential-performance-requirement-assessment-of-prefilled-syringes-break-loose-and-gliding-force-measurement-and-analysis

Reading: the compliance budget line is for *physical verification* of force/time/dose under stressed conditions. A model enters only if it reduces the test matrix (worst-case selection, temperature envelope) or explains a failure. USP <1>: not found relevant.

### Q6. Technical unmet needs stated by practitioners

- **Concentration extrapolation fails.** GSK/Strathclyde (mAbs 2024): "Our findings caution against extrapolation from low concentration measurements, as they lack predictive ability for ultra-high concentration regimes"; "% differences -186% to +177%"; "Predictive regression models from small datasets are generally overfitted and lack generalizability." https://pmc.ncbi.nlm.nih.gov/articles/PMC11790245/
- **Force models drift at high concentration and fine needles.** JHU thesis: "Deviations from theoretical models were identified in drug product concentrations above 230 mg/mL and in thin needles (27G, 30G) potentially due to ... insufficient characterization of high shear viscosity behavior"; models "<25% error"; rule of thumb "27-gauge thin walled needle, the liquid viscosity should be maintained at or below 20 cP in order to not exceed the glide force of 25 Newton". https://jscholarship.library.jhu.edu/items/32605470-29f4-42c5-a309-a530919d5d64
- **Needle shear rate.** "a 10-second injection of 1 mL through a 27 G needle presents shear rates up to 10^5 s-1". https://pmc.ncbi.nlm.nih.gov/articles/PMC13119207/
- **Shear-thinning matters in the model.** "using a Newtonian model will overestimate injection time should shear thinning occur"; narrower needle "increases dynamic pressure for Newtonian fluids but decreases for shear-thinning". https://www.tandfonline.com/doi/full/10.2147/MDER.S151727 ; EJPB 2024 model with Newtonian and shear-thinning branches and injection temperature https://www.sciencedirect.com/science/article/abs/pii/S093964112400047X
- **Component attributes.** WuXi Biologics (Biopharmaceutics, Oct 2025): 29G to 27G "reduced glide force by 59%"; viscosity "nonlinear increase ... can impair injectability". https://pubmed.ncbi.nlm.nih.gov/41093078/
- **Cold-chain viscosity.** Vendor app brief: "viscosity at 5 C is 43 cP whereas that at 15 C is 24 cP"; 30-min equilibration "often not feasible or adhered to"; must "gauge the viscosity at temperatures ranging from 2-25 C". https://www.rheosense.com/m-vrocii/application-brief-auto-injector-high-viscosity-testing (vendor) ; https://www.labmanager.com/rheological-properties-of-injectable-drugs-34794
- **Tissue back-pressure.** Minipig: "1 mL injection in 10 seconds (360 mL/h) generated 24.0 +/- 3.4 kPa, whereas 10 mL in 10 minutes (60 mL/h) generated 7.4 +/- 7.8 kPa". https://link.springer.com/article/10.1007/s11095-020-02906-9 ; review of 31 LVSC studies https://www.tandfonline.com/doi/full/10.1080/17425247.2023.2219891 ; https://www.sciencedirect.com/science/article/abs/pii/S0022354916413791
- **Volume ceiling moving.** ADDR (July 2026): "SC administration volumes of well above 2 mL, and up to 10 and even 25 mL are increasingly pursued". https://pubmed.ncbi.nlm.nih.gov/41935636/ ; BMS review (2026): "assumptions regarding SC volume limits (commonly cited up to 2 ml) are increasingly challenged"; AI/PFS "5-10 s", OBDS "5-60 min"; **no** shear-rate or non-Newtonian discussion in that review. https://pmc.ncbi.nlm.nih.gov/articles/PMC13248910/
- **Shear / aggregation.** 2025 J Pharm Sci: "elevated sliding speeds promoted proteinaceous subvisible and/or visible particle formation ... presumed to be caused by fine air bubbles generated due to rapid pressure changes arising from shear stress during needle passage". https://pubmed.ncbi.nlm.nih.gov/39862973/ ; interfacial mechanism https://www.researchgate.net/publication/261568698_Protein_Aggregation_and_Particle_Formation_in_Prefilled_Glass_Syringes
- **TFF/UF of viscous mAbs.** "high viscosity may lead to high feed pressure, exceeding upper limit for membrane integrity". https://pubmed.ncbi.nlm.nih.gov/41707938/ ; https://www.sciencedirect.com/science/article/abs/pii/S0376738816300886
- **Ophthalmic.** Intravitreal: ophthalmic-specific syringes "required significantly lower injection forces than standard syringes" for high-viscosity pegcetacoplan. https://pmc.ncbi.nlm.nih.gov/articles/PMC13147785/
- **Survey signal.** "A recent survey from formulation experts reported that challenges in developing high-concentration formulations ... could result in delays in clinical trials or product launches" (AAPS Open 2026). https://link.springer.com/article/10.1186/s41120-026-00160-8
- Nasal injection: not searched. Depot gels: literature uses texture-analyser force tests with 23G needles; no service market found. https://www.researchgate.net/publication/51108443_Injectability_Evaluation_An_Open_Issue

### Q7. Startups / companies 2020-2026 and funding

| Company | Founded | What | Money | Source |
|---|---|---|---|---|
| Surf Bio (acquired by Halozyme) | preclinical | Spray-dry + protective excipient "hyperconcentration" to 500 mg/mL | Acquired Dec 2025: $305M consideration, up to $100M milestones ("up to $400M") | https://www.fiercepharma.com/pharma/halozyme-catches-ma-wave-snap-surf-bio-delivery-tech-deal-worth-400m ; https://www.contractpharma.com/breaking-news/halozyme-acquires-surf-bio-for-up-to-400m/ |
| Halozyme (ENHANZE) | incumbent | rHuPH20 for large-volume SC | FY2025 revenue $1.397B, royalties $868M (+52 %); guides >$1B royalties 2026 | https://www.prnewswire.com/news-releases/halozyme-reports-full-year-2025-record-revenue-of-1-4-billion-and-reiterates-strong-2026-financial-guidance-302689822.html (2026-02-17) |
| Elektrofi | 2015 (Boston) | Non-aqueous microparticle suspensions >400 mg/mL, "<20 N" force | $40M Series B, June 2022 (Marshall Wace, BVF) | https://www.businesswire.com/news/home/20220613005094/en/ |
| Comera Life Sciences | public then private | Caffeine-class viscosity-reducing excipients (SQore); deals with Astellas, Bayer, KBI; MilliporeSigma licence 2020 | Taken private Jan 2024 by bridge investors | https://ir.comeralifesciences.com/news-releases/news-release-details/comera-life-sciences-announces-expansion-patent-portfolio ; https://stockanalysis.com/quote/otc/CMRA/ |
| QATCH Technologies | 2016 (Duke spin-out) | nanovisQ microliter viscometer, visQ.AI | ~$3M in 6 SBIR/STTR; $1.5M + $1.88M Phase IIs | https://sbtdc.org/news/nc-sbir-success-story-qatch-technologies ; https://www.sbir.gov/awards/200235 |
| Persist AI | 2022 | Robotic cloud lab for formulation; Lilly investor | $12M Series A (Spero), ~$16.5M total; Cloud Lab May 2025 | https://intuitionlabs.ai/articles/top-biotech-startups-2026 (secondary; primary release not fetched) |
| Altruist Biologics (Innovent subsidiary) | 2026 launch | Altru-CON platform to 200 mg/mL, ">two-fold viscosity reduction" | n/a | https://tools.prnewswire.com/en-us/live/20813/release/20260304EN99994 (2026-03-04, vendor) |
| Excelse Bio | -- | Amino-acid blend high-concentration formulation | "has not raised any funding rounds yet" (Tracxn) | https://tracxn.com/d/companies/excelse-bio/__yqh99WrqmBjSLsg4plP_l1Ecfa7hwE5X0FMhnFUVn4g |
| Arecor | UK, listed | High-concentration formulation tech | no 2024-25 figure found | https://arecor.com/technology/high-concentration-antibodies/ |
| BD + Ypsomed | 2025 collaboration | Neopak XtraFlow PFS + YpsoMate 2.25 for ">15cP" biologics | n/a | https://drug-dev.com/bd-ypsomed-collaborate-to-advance-self-injection-systems-for-high-viscosity-biologics/ |
| Nemera Symbioze | -- | Reusable on-body injector, motorised for viscous fluids | n/a | https://drug-dev.com/special-feature-injection-devices-faster-integrated-and-flexible-for-the-iv-to-sc-transition/ (vendor) |

Injectability-*modelling* software startups (2020-2026): **not found** other than QATCH's visQ.AI pilot. RheoSense: no 2025 funding or acquisition announcement found.

---

## Named buyers and named sellers

| Role | Names | Evidence |
|---|---|---|
| Buyers (pharma with stated high-viscosity/LVSC programmes) | Pfizer (LVSC talk), Roche (in-house injection-time model), GSK (viscosity prediction paper), BMS (SC high-concentration review), WuXi Biologics (PFS injectability paper), AbbVie, Amgen, Lilly, Regeneron, Novartis, Sanofi, Takeda, UCB, argenx (PODD presenters / ENHANZE licensees) | Q3, Q6, Q7 sources |
| Buyers (device makers who need drug-side data) | SHL, Ypsomed, BD, Stevanato, Gerresheimer, Nemera, West | Q1, Q7 |
| Sellers, measurement | West ALS, Smithers, Cormica, Coriolis, Lonza, Stevanato, Nelson Labs (exhibits), Eurofins (hiring) | Q1 |
| Sellers, instruments | RheoSense, QATCH, ZwickRoell, Instron, Stable Micro Systems | Q1, Q4 |
| Sellers, modelling | Veryst, Springboard, Team Consulting, Cambridge Consultants (not verified), Kymanox (regulatory/engineering) | Q1 |
| Sellers, enablers (where the money actually goes) | Halozyme, Surf Bio, Elektrofi, Comera, WuXi (WuXiHigh), Altruist | Q7 |

"Someone sells this" is established for measurement and generic modelling. "Someone buys *analysis* from outside" is established only indirectly (CRO service pages exist and persist; Coriolis lists non-Newtonian awareness). No buyer quote of the form "we paid X for an injectability assessment" was found.

---

## Unmet needs a 1-3 person company could credibly address

1. **Needle-shear extrapolation with a confidence band.** Fit Cross/Carreau/Casson-type laws to 1e0-1e4 1/s rheometer data, propagate profile-likelihood uncertainty to injection force/time at 1e4-1e5 1/s in 27G/29G/30G and thin-wall geometries, and state identifiability (which parameters the force actually depends on). This is exactly where the JHU and GSK papers report failure and where no seller states an error bar. (The founder's biorheo codebase already does the fitting and profiling.)
2. **Temperature envelope.** Deliver force/time vs temperature 2-25 C with the 30-min equilibration assumption removed; the vendor app brief shows a ~1.8x viscosity swing over 10 C.
3. **Worst-case selection for design verification.** Use the model to pick which (viscosity, temperature, needle lot) corners must be physically tested under ISO 11608-1/11040-8 and the EDDO draft, shrinking a test matrix rather than replacing it.
4. **Family-choice sensitivity.** Show how much the predicted force moves between constitutive families fitted to the same data (biorheo item 12 found up to 30 % for blood); pharma groups using a single Hagen-Poiseuille spreadsheet have no such number.
5. **Shear-exposure and bubble/cavitation screening.** Compute wall shear stress, residence time and pressure transients in the fluid path for aggregation-risk ranking; pair with a benchtop particle count.
6. **Large-volume SC coupling.** Couple a non-Newtonian fluid-path model to a poroelastic tissue back-pressure model (Veryst does tissue; nobody found does both with drug rheology), calibrated on the published minipig kPa data.
7. **Suspensions (400-500 mg/mL spray-dried particles).** Surf Bio/Elektrofi/Halozyme now make the working fluid a dense suspension, i.e. the founder's PhD physics (jamming, shear thickening, particle friction). No injectability service found for suspension formulations.

---

## Evidence AGAINST this being a good market

- **No price signal.** Not one public rate for an injectability study, ISO 11608 campaign or CFD study; the only dollar figures are grants to an instrument company and M&A for enablers.
- **The core physics is cheap and old.** Roche's 2018 model hits <12 % error; the SHL/Stevanato 2026 study used plain Hagen-Poiseuille and got a saleable result. Buyers may see uncertainty and non-Newtonian corrections as refinement, not need.
- **Buyers are self-sufficient.** Big pharma device groups own Instron/Zwick rigs, VROC viscometers and in-house modellers (Roche, GSK, BMS, WuXi all publish their own).
- **Established modelling consultancies exist** (Veryst, Springboard, Team) with COMSOL/ANSYS credentials and device pedigrees; a newcomer competes on rigour, not presence.
- **Regulatory pull funds tests, not models.** EDDO/ISO require physical verification after stress conditioning; a model is not evidence unless the sponsor argues it is.
- **Job postings do not name the skill.** Injectability appears as one qualification among many in formulation roles; no dedicated role was found.
- **The value pool went to excipients and devices**, which change the fluid rather than analyse it (Surf Bio up to $400M; Elektrofi $40M; WuXiHigh; Halozyme $868M royalties).
- **Access:** many primary pages (careers, ScienceDirect, ONdrugDelivery, Kymanox) are gated; several claims above rest on search snippets.

---

## Open questions requiring a human conversation

1. What does a mid-size biotech or CDMO actually pay a CRO for a syringeability/injectability package (per formulation, per device)? Ask Coriolis, West ALS, Cormica for a quote on a 3-formulation, 2-needle study.
2. Has any sponsor used a rheology-based injection model as leverage in a design-verification file or a regulatory response, and did FDA accept it? (Ask a Kymanox / former FDA CDRH reviewer.)
3. Do formulation groups measure above 1e4 1/s routinely (VROC/nanovisQ), or do they extrapolate? If they measure, the extrapolation wedge shrinks to uncertainty and identifiability.
4. Is the 20-cP / 25-N / 10-s rule of thumb still what device teams use for go/no-go, and who owns the decision (formulation vs device)?
5. For 400-500 mg/mL suspensions (Halozyme/Surf Bio, Elektrofi), who is characterising suspension rheology and needle flow, and is that work outsourced?
6. Would a device maker (SHL, Ypsomed, Stevanato) resell an independent drug-side model to its pharma customers as a service?
7. Cold-chain: do sponsors verify injection time at 5 C or only at 23 +/- 5 C per the 2013 guidance, and has FDA ever asked for the cold case?

---

## Sources

Fetched or cited pages (date accessed 2026-09-06 unless a publication date is given):

- https://www.westpharma.com/blog/2023/january/auto-injector-iso-11608-5-regulatory-requirements-analytical-services
- https://www.smithers.com/industries/medical-device/products-testing/injectables/cartridge-testing
- https://www.cormica.com/services/physical-and-packaging-testing/physical-mechanical-testing-for-medical-devices/iso-11608-needle-based-injection-system-testing/
- https://www.zwickroell.com/industries/medicalpharmaceutical/therapy-systems/autoinjectors-iso-11608-5/
- https://www.instron.com/en/testing-solutions/iso-standards/iso-11608/
- https://www.instron.com/en/testing-solutions/iso-standards/iso-11040/
- https://www.instron.com/en/resources/blog/2026/june/iso-11040-8-2026-pre-filled-syringe-testing-guide/
- https://www.ondrugdelivery.com/navigating-iso-11608-with-automated-equipment-for-customised-pen-and-autoinjector-testing/
- https://www.coriolis-pharma.com/our-services/injectability-amp-syringeability/
- https://www.coriolis-pharma.com/our-services/capillary-viscometry/
- https://www.lonza.com/integrated-biologics/parenteral-drug-product/biologics-special-services
- https://www.wuxibiologics.com/drug-product-development-services/
- https://pharmasource.global/content/news/cdmo-news/wuxi-biologics-launches-high-concentration-formulation-platform-wuxihigh2-0/
- https://www.stevanatogroup.com/en/analytical-testing-services/analytical-services/
- https://www.stevanatogroup.com/en/news-events/webinars/overcoming-volume-and-viscosity-barriers-in-biologic-formulations/
- https://www.veryst.com/case-studies/optimizing-subcutaneous-drug-delivery-multiphysics-simulation
- https://www.veryst.com/case-studies/optimizing-intravitreal-drug-delivery-multiphysics-simulation
- https://www.springboard.pro/case-studies/gx-inbeneo-autoinjector-design/
- https://www.springboard.pro/autoinjectors-vs-on-body-systems-the-future-of-large-volume-delivery/
- https://www.team-consulting.com/us/insights/pda-universe-of-pre-filled-syringes-and-injection-devices-conference/
- https://qatchtech.com/visq-ai
- https://qatchtech.com/
- https://blog.rheosense.com/predicting-injection-force-theoretical-foundation
- https://blog.rheosense.com/newtonian-and-non-newtonian-injection-force
- https://www.rheosense.com/m-vrocii/application-brief-auto-injector-high-viscosity-testing
- https://www.tandfonline.com/doi/full/10.2147/MDER.S151727
- https://www.shl-medical.com/news-insights/blog/high-viscosity-delivery-molly-autoinjector
- https://www.indeed.com/q-Pharmaceutical-Biologics-Formulation-Scientist-jobs.html
- https://www.indeed.com/q-injectable-formulation-scientist-jobs.html
- https://www.indeed.com/q-Biomedical-Engineer-Drug-Delivery-System-jobs.html
- https://careers.abbvie.com/en/job/senior-scientist-i-biologics-drug-product-development-in-south-san-francisco-ca-jid-10300
- https://careers.abbvie.com/en/job/senior-device-development-engineer-in-north-chicago-il-jid-26583
- https://careers.regeneron.com/en/jobs/r43812/sr-scientist-formulations-drug-product-development/
- https://www.linkedin.com/jobs/view/biologics-formulation-scientist-temporary-at-eurofins-4448630800
- https://www.pda.org/rest-api/event-agenda/pda-universe-of-pre-filled-syringes-and-injection-devices-conference-2025
- https://www.pda.org/global-event-calendar/abstract-submissions/call-for-abstracts/pda-universe-of-pre-filled-syringes-and-injection-devices-conference-2025
- https://www.pda.org/global-event-calendar/event-detail/pda-universe-of-pre-filled-syringes-and-injection-devices-conference-2026
- https://intuitionlabs.ai/conferences/pda-universe-of-pre-filled-syringes-injection-devices
- https://www.pharmtech.com/view/podd-2025-unveils-agenda-highlighting-innovation-in-drug-delivery-technologies
- https://www.poddconference.com/
- https://aaps2025.eventscribe.net/fsPopup.asp?PosterID=764032&mode=posterInfo
- https://www.pda.org/docs/default-source/posters/universe-of-prefilled-syringes-and-injection-devices/chen-wei.pdf
- https://www.sciencedirect.com/science/article/abs/pii/S093964112400047X
- https://www.ondrugdelivery.com/trends-for-2026/
- https://biobridgex.com/costs/bioanalytical-testing-costs
- https://www.iso.org/standard/70733.html
- https://cdn.standards.iteh.ai/samples/76627/440310c2283e45dea7ca5908c4ae7975/ISO-11608-5-2022.pdf
- https://www.sbir.gov/awards/200235
- https://ars-pharma.inknowvation.com/sbir/awards/nih-2023-high-throughput-injectability-screening-high-concentration-protein-formulations
- https://sbtdc.org/news/nc-sbir-success-story-qatch-technologies
- https://www.ncbiotech.org/news/ncbiotech-awards-22-million-grants-loans-latest-quarter
- https://www.prnewswire.com/news-releases/niimbl-announces-10m-for-8-new-biopharmaceutical-manufacturing-projects-302112132.html
- https://www.niimbl.org/news/niimbl-announces-15-8m-to-fund-14-new-biopharmaceutical-manufacturing-projects/
- https://www.federalregister.gov/documents/2024/07/01/2024-14409/essential-drug-delivery-outputs-for-devices-intended-to-deliver-drugs-and-biological-products-draft
- https://www.thefdalawblog.com/2024/07/fda-releases-draft-guidance-on-essential-drug-delivery-outputs/
- https://www.westpharma.com/blog/2024/september/fda-guidance-essential-drug-delivery-output-requirements
- https://www.kymanox.com/thought-leadership-series/blog/fda-guidance-essential-drug-delivery-outputs-kymanox-what-you-need-to-know/
- https://www.kymanox.com/thought-leadership-series/blog/large-volume-subcutaneous-injection/
- https://www.fda.gov/files/about%20fda/published/Technical-Considerations-for-Pen--Jet--and-Related-Injectors-Intended-for-Use-with-Drugs-and-Biological-Products.pdf
- https://www.federalregister.gov/documents/2013/06/07/2013-13484/
- https://www.fda.gov/media/85748/download
- https://www.bioprocessintl.com/fill-finish/essential-performance-requirement-assessment-of-prefilled-syringes-break-loose-and-gliding-force-measurement-and-analysis
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11790245/
- https://jscholarship.library.jhu.edu/items/32605470-29f4-42c5-a309-a530919d5d64
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13119207/
- https://pubmed.ncbi.nlm.nih.gov/41093078/
- https://www.labmanager.com/rheological-properties-of-injectable-drugs-34794
- https://link.springer.com/article/10.1007/s11095-020-02906-9
- https://www.tandfonline.com/doi/full/10.1080/17425247.2023.2219891
- https://www.sciencedirect.com/science/article/abs/pii/S0022354916413791
- https://pubmed.ncbi.nlm.nih.gov/41935636/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13248910/
- https://pubmed.ncbi.nlm.nih.gov/39862973/
- https://www.researchgate.net/publication/261568698_Protein_Aggregation_and_Particle_Formation_in_Prefilled_Glass_Syringes
- https://pubmed.ncbi.nlm.nih.gov/41707938/
- https://www.sciencedirect.com/science/article/abs/pii/S0376738816300886
- https://pmc.ncbi.nlm.nih.gov/articles/PMC13147785/
- https://link.springer.com/article/10.1186/s41120-026-00160-8
- https://www.researchgate.net/publication/51108443_Injectability_Evaluation_An_Open_Issue
- https://www.tandfonline.com/doi/full/10.1080/19420862.2025.2550757
- https://www.fiercepharma.com/pharma/halozyme-catches-ma-wave-snap-surf-bio-delivery-tech-deal-worth-400m
- https://www.contractpharma.com/breaking-news/halozyme-acquires-surf-bio-for-up-to-400m/
- https://www.prnewswire.com/news-releases/halozyme-reports-full-year-2025-record-revenue-of-1-4-billion-and-reiterates-strong-2026-financial-guidance-302689822.html
- https://www.businesswire.com/news/home/20220613005094/en/
- https://ir.comeralifesciences.com/news-releases/news-release-details/comera-life-sciences-announces-expansion-patent-portfolio
- https://stockanalysis.com/quote/otc/CMRA/
- https://intuitionlabs.ai/articles/top-biotech-startups-2026
- https://tools.prnewswire.com/en-us/live/20813/release/20260304EN99994
- https://tracxn.com/d/companies/excelse-bio/__yqh99WrqmBjSLsg4plP_l1Ecfa7hwE5X0FMhnFUVn4g
- https://arecor.com/technology/high-concentration-antibodies/
- https://drug-dev.com/bd-ypsomed-collaborate-to-advance-self-injection-systems-for-high-viscosity-biologics/
- https://drug-dev.com/special-feature-injection-devices-faster-integrated-and-flexible-for-the-iv-to-sc-transition/
