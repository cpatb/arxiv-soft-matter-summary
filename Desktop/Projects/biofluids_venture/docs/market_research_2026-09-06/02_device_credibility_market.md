# 02 — The market for computational-model credibility work in medical devices and pharma

Research date: 2026-09-06/07. Method: ~35 web searches and ~15 page/PDF fetches. Every factual claim below carries a URL; "not found" means I searched and could not locate public evidence, not that it does not exist. Vendor pages are marked **[vendor marketing]**. Market-size press releases are marked **[analyst PR, unaudited]**. No dollar figure below is estimated; each is quoted from its source.

## Summary

1. There is a real, regulator-endorsed demand pull: FDA's final guidance (16 Nov 2023) adopts ASME V&V 40-2018 as the credibility framework; the EU MDR mentions CM&S "only superficially" and notified bodies still prefer standard-based bench testing (46-respondent survey, 2025).
2. Volume is documented only indirectly: 220 OSEL consults involved CM&S in 2017 (44 % of surveyed consults, 9 % of all OSEL consults); FDA says V&V 40 has been used in "scores" of device submissions (2022); 65.7 % of 229 symposium attendees said adoption grew over five years (2024). No FDA year-by-year count for 2020-2026 exists in public: **not found**.
3. The seller side is crowded at the top (Ansys/Synopsys, Exponent, Thornton Tomasetti, Stress Engineering, Veryst, SimuTech) and at the platform end (Virtonomy, InSilicoTrials, ELEM, Adsilico, Numalogics); all use V&V 40 language, none publishes prices.
4. Buyers mostly do it in-house: in MDIC's 2021 survey 72 % had internal resources and only 23 % used consultants when required; the top barriers were "uncertainty in regulator expectations" (61 %) and "lack of expertise" (49 %), which is the consultant's opening.
5. Small firms demonstrably win regulatory-grade work: Numalogics (Canada, founded by three surgeons) obtained the first MDDT-qualified virtual mechanical test (qualified 10 Apr 2025); Veryst (~18 staff) cites 510(k)/PMA support; a boutique (CardioMed, founded 2003) sat on the FDA/MDIC credibility panel.
6. Non-Newtonian rheology is a stated, quantified problem in the FDA nozzle benchmark (Newtonian model underestimates NIH by 8 %; up to ~6x spread among non-Newtonian models at low Re), but hemolysis validation is limited by experimental scatter (47-94 % intra-lab CV) more than by the CFD.
7. No public price for a regulatory CFD study or credibility package was found; the only hard money is grant ceilings (NIBIB SBIR Phase I up to $256K, Phase II up to $1.7M) and start-up rounds (Adsilico $4.4M, ELEM $2.27M).
8. Software that assembles the evidence package exists inside solver vendors (Ansys Minerva V&V 40 templates; Virtonomy auto-generated "FDA 1807 & ASME V&V40" reports); a standalone, solver-independent "V&V 40 manifest" product was **not found**.
9. Supplier-side barriers are ordinary: ISO 13485 requires risk-proportionate supplier controls, several simulation vendors advertise ISO 13485 certification as a differentiator, and E&O insurance is a standard contract clause; no evidence that a solo shop is formally excluded (and none that it is routinely admitted either: **not found**).
10. Verdict: a 1-3 person shop can plausibly sell subcontracted CFD + credibility-plan authoring to small/mid sponsors and to the platform companies, but must expect price opacity, in-house competition, a relationship-driven sales cycle, and a bench-validation requirement it does not yet own.

## Evidence by question

### Q1. Who sells this today?

Named sellers with public V&V 40 / regulatory-simulation language (size and lab ownership in the table below):

- **Stress Engineering Services** publishes its full internal QMS procedure "V&V 40: Computational Modeling for Medical Devices, L2-CIN-PR-22 Rev 0, 17-Aug-2022" (26 pp), with 5x5 model-influence x decision-consequence risk tables, credibility-factor tables (SQA, numerical code verification, discretization error, solver error, use error; model form, inputs, test samples, test conditions, equivalency of inputs, output comparison; relevance of QOIs and validation activities) and the project manager as responsible party. It is the only public example of what a subcontractor's V&V 40 deliverable structure looks like: https://sesforce4.com/wp-content/uploads/2024/04/L2-CIN-PR-22_Rev0-VV40.pdf . Its services page claims "multiple successful regulatory submissions to the FDA and EU MDR using CM&S based approaches" for in silico design verification **[vendor marketing]**: https://www.stress.com/in-silico-design-verification-testing-of-medical-devices/ ; also https://www.stress.com/fda-cranks-up-pressure-for-modeling-and-simulation/
- **Thornton Tomasetti (Applied Science)** "supports clients from verification and validation (ASME V&V 40) over risk management to their regulatory submissions" and describes benchtop tests used to validate simplified models **[vendor marketing]**: https://www.thorntontomasetti.com/news/fda-guidance-computational-modeling-medical-device-submissions ; https://www.thorntontomasetti.com/solution/biomechanics-life-sciences ; CFD capability page: https://www.thorntontomasetti.com/capability/computational-fluid-dynamics . Firm-wide 144 open positions (ZipRecruiter, 2026): https://www.ziprecruiter.com/Jobs/Thornton-Tomasetti
- **Veryst Engineering** (Needham MA): FEA support "for Class I, II, and III device regulatory documents", "design validation experiments for fatigue and strength for FDA submissions", 510(k)/PMA/CE; a case study simulating the FDA benchmark nozzle with several hemolysis models **[vendor marketing]**: https://www.veryst.com/services/specialized-expertise/medical-devices ; https://www.veryst.com/services/simulation-analysis/computational-fluid-dynamics-cfd . Headcount ~18 (third-party estimate, July 2026): https://leadiq.com/c/veryst-engineering-llc/5a1d83d424000024005e5717
- **Exponent**: FEA for medical devices/biomaterials and "biofluid dynamics" including "assessment of dead zones, hemolysis levels, and pressures on valve leaves in artificial heart valves" **[vendor marketing]**: https://www.exponent.com/capabilities/medical-devices ; https://www.exponent.com/capabilities/finite-element-analysis-medical-devices-biomaterials ; https://www.exponent.com/capabilities/biofluid-dynamics-biomedical-devices-products
- **SimuTech Group** (Ansys channel partner): principal CFD engineer Mark Goodin (25 years cardiovascular devices) sits on the ASME V&V 40 committee; webinar claims a V&V 40-framed model "accepted by the FDA as part of a successful premarket submission" and "reduced the number of required physical tests by up to 75 %" **[vendor marketing]**: https://go.simutechgroup.com/medical-device-1 ; https://simutechgroup.com/asme-vv-40-standard-medical-devices/ ; https://simutechgroup.com/industry-solutions/healthcare/
- **Numalogics** (Canada, founded 2010 by three spinal surgeons; ASME V&V 40 contributor): with Sawbones obtained the first MDDT qualification of an orthopaedic mechanical test as an FE model (ASTM F543 screw pullout), qualified 10 Apr 2025, announced 20 May 2025: https://numalogics.com/wp-content/uploads/2025/05/PR-Screw-Pullout-MDDT-2025-05-20.pdf ; https://sawbones.com/blog/firstever-fda-qualification-of-an-orthopedic-device-mechanical-test-under-mddt-program/ ; Ansys profile: https://www.ansys.com/blog/democratizing-in-silico-medical-device-testing
- **Virtonomy** (Munich, founded 2019): v-Patients platform, "FDA-compliant reporting (FDA 1807 & ASME V&V40)", "Validated ISO benchtop test simulations" for ISO 5840 / 7198 / 25539; two named customers (Medira GmbH, Scandinavian RealHeart); no pricing **[vendor marketing]**: https://www.virtonomy.io/product/ . Investor note (403 on fetch): https://medium.com/pace-ventures/why-we-invested-in-virtonomy-53c0f14aa273
- **ELEM Biotech** (Barcelona, founded 2018, Alya-based virtual humans for cardiac and respiratory/inhaler): $2.27M raised across 3 rounds, latest an EIC grant in 2023 (third-party database; may be incomplete): https://tracxn.com/d/companies/elem-biotech/__hHHepAcMcrWoxjD_nFg2gsdy3Qr8LMc1TPCapR2vsfg ; https://www.elem.bio/
- **Adsilico** (UK): raised $4.4M (GBP 2M Northern Gritstone + GBP 1.5M Parkwalk), Feb 2024: https://www.medicaldevice-network.com/news/uks-adsilico-raises-4-4m-to-develop-virtual-trial-solution/
- **InSilicoTrials** (cloud platform for device and drug in silico trials): https://mtec-sc.org/life-sciences/insilicotrials ; comparison with Virtonomy: https://www.cbinsights.com/compare/insilico-trials-vs-virtonomy
- **Ansys (Synopsys)**: Ansys Minerva has a "template-based, end-to-end workflow" for V&V 40 credibility documentation (Q5); Synopsys Simpleware Medical holds 510(k) clearance and CE mark for image-to-model, developed under ISO 13485: https://www.synopsys.com/simpleware/resources/regulatory.html
- **PADT** (Ansys partner) advertises ISO 13485 certification with FEA/CFD services: https://padtinc.com/blog/the-focus/fda-opening-to-simulation-supported-verification-and-validation-for-medical-devices ; **Simq GmbH** advertises ISO 13485-certified software development and sat on the 2024 FDA/MDIC panel: https://simq.de/en/software-development/
- Regulatory boutiques on the 2024 FDA/MDIC credibility panels (they broker or buy this work rather than run solvers): NC MedTech Solutions LLC, MCRA, Lighthouse Regulatory, CardioMed Device Consultants (founded Jan 2003, Tentamus Group), Voisin Consulting; notified bodies TÜV SÜD and BSI also attended: https://iopscience.iop.org/article/10.1088/2516-1091/ae1c05 ; https://www.cardiomedllc.com/ ; https://www.tentamus.com/locations/cardiomed-device-consultants-llc/
- Content-marketing pages without a stated V&V 40 service: StarFish Medical https://starfishmedical.com/resource/fda-guidance-on-cms-in-medical-device-submissions/
- **Not found** in a regulatory-CFD capacity: CFD Direct, Sim&Cure, Materialise (V&V 40 services), Confluent Medical.

Bench labs: no simulation consultancy page states that it owns a blood flow loop or hemolysis lab. Hemocompatibility testing is sold by contract labs, e.g. Nelson Labs (ISO 10993-4 mechanically induced hemolysis): https://www.nelsonlabs.com/mechanically-induced-hemolysis-iso-10993-4/ . Thornton Tomasetti and Veryst describe designing or running benchtop validation tests; ownership of a fluids lab is not confirmed on their pages.

### Q2. How often do FDA submissions include CM&S?

- 2017 OSEL survey (Morrison et al., Frontiers in Medicine 2018): "Of the 500 consults performed by the scientists surveyed, 220 (44 %) included computational modeling and digital evidence in the submission", equivalent to "9 % of all expert consults performed by OSEL"; of the 220: 36 % 510(k), 25 % clinical-trial applications, 24 % pre-submissions, 13 % PMA; largest areas neurological, then cardiovascular and orthopaedic. Primary use: "to identify the appropriate bench testing configurations"; second: MRI RF safety. https://pmc.ncbi.nlm.nih.gov/articles/PMC6167449/
- Morrison (J. Med. Devices 2017, quoted by StarFish): CM&S "is currently 15 % of gathered evidence for medical device submissions, but is projected to increase to nearly 50 %" — a projection, not a count: https://starfishmedical.com/resource/computational-modelling-simulation-medical-device-development/ ; original: https://www.feops.com/src/Frontend/Files/MediaLibrary/08/morrison-j-med-devices-11-2017-role-of-computational-modeling.pdf
- FDA Senior Science Council M&S Working Group report (2022), p. 13: ASME V&V 40 "has been used in scores of medical device regulatory submissions and a handful of CDER Investigational New Drug Applications (IND) and fit-for-purpose applications": https://www.fda.gov/media/163156/download
- 2024 FDA/MDIC symposium (229 attendees): 65.7 % said CM&S "has become more widely adopted for use in medical device regulatory submissions over the past five years", 4.9 % said no; top impediments "Uncertainty in CM&S Expectations by Regulatory Bodies", "Challenges Demonstrating Credibility", "Lack of Model Credibility Assessment Best Practices or Examples": https://iopscience.iop.org/article/10.1088/2516-1091/ae1c05 ; https://pubmed.ncbi.nlm.nih.gov/41192059/ ; presentations: https://mdic.org/resources/2024-fda-mdic-symposium-on-computational-modeling-and-simulation-presentations/
- Final guidance issued 16 Nov 2023 (docket FDA-2021-D-0980): https://www.federalregister.gov/documents/2023/11/17/2023-25470/assessing-the-credibility-of-computational-modeling-and-simulation-in-medical-device-submissions ; PDF: https://www.fda.gov/media/154985/download ; FDA webinar: https://www.fda.gov/medical-devices/medical-devices-news-and-events/webinar-assessing-credibility-computational-modeling-and-simulation-medical-device-submissions-final
- FDA's own framing (Pathmanathan, ICCVAM CoP, Jan 2024): the guidance "extends approach of ASME V&V40 2018"; "ASME V&V40 implicitly assumes validation against prospective well-controlled bench tests"; the deck contains no submission statistics: https://ntp.niehs.nih.gov/sites/default/files/2024-02/2_Pathmanathan_ICCVAM-CoP-2024_508.pdf
- EU MDR (survey of 46 respondents: 26 manufacturers, 11 notified bodies, 8 consultants; TIRS, July 2025): CM&S "were not widely observed, but seem to get more relevance in higher risk products such as class IIb implants and class III products (more than 30 %)" versus under 30 % for class I; "The MDR has addressed CM&S methods only superficially and without specifying criteria for acceptance or credibility"; notified bodies favour "standard-based bench testing"; CM&S considered appropriate by >60 % of consultants, >40 % of notified bodies, ~20 % of manufacturers: https://pmc.ncbi.nlm.nih.gov/articles/PMC12579633/
- Year-by-year FDA counts of CM&S-containing submissions 2020-2026: **not found** (searched OSEL pages, the Jan 2024 FDA deck, the 2024 FDA/MDIC proceedings; none gives a count).
- Drug analogue: QSP submissions to FDA "more than doubled since 2021" (2024 analysis): https://www.drugdiscoverynews.com/how-modeling-and-simulation-are-reshaping-regulatory-science-16617

### Q3. Do small firms/subcontractors get this work?

- Buyer-side: MDIC industry survey (fielded 2021, published 2024): 72 % "had resources available internally", 64 % "personally apply CM&S regularly", 36 % "rely on software company expertise", 23 % "use consultants or apply technology when specifically required"; 69 % agree CM&S reduces time-to-market; barriers: 61 % regulator-expectation uncertainty, 49 % "lack of expertise in CM&S", 39 % business case unclear: https://medtechintelligence.com/news_article/medical-device-computational-modeling-and-simulation-reports-highlights-best-uses-and-barriers-to-adoption/ ; report landing page (PDF is form-gated): https://mdic.org/resources/landscape-report-industry-survey-on-the-use-of-computational-modeling-simulation-in-medical-device-development/
- Evidence small firms win: Numalogics + Sawbones MDDT qualification (Q1); Veryst at ~18 staff cites 510(k)/PMA support (Q1); FDA/MDIC 2024 panels seated boutiques (CardioMed, NC MedTech, Lighthouse) beside Medtronic, J&J, BD, Stryker: https://iopscience.iop.org/article/10.1088/2516-1091/ae1c05 . The same proceedings note "challenges that startups and innovators face in securing funding and support for CM&S" — start-ups are underserved buyers.
- In-house teams at large sponsors: Medtronic Distinguished Scientist Jeff Bodner served on the V&V 40 subcommittee (DMD 2025 panel with Gore, Boston Scientific, Numalogics, Ansys): https://dmd.umn.edu/2025/vv40 ; Medtronic "Computational Modeling Engineer" posting (May 2024, expired): https://www.jobzmall.com/medtronic/job/computational-modeling-engineer ; Abbott Heart Failure "Senior Product Development Engineer – Computational Modeling" (FEA/CFD/thermal; removed 9 Jan 2025): https://builtin.com/job/senior-product-development-engineer-computational-modeling/3757574
- Consultancy hiring: Numalogics open FE simulation engineer role: https://numalogics.com/careers/finite-element-simulation-engineer/ ; ZipRecruiter average "computational modeling" pay $54.93/h (range $46.88-$73.56/h, Aug 2026) — a cost-of-labour proxy, not a billing rate: https://www.ziprecruiter.com/Jobs/Computational-Modeling
- Job postings containing the literal phrase "V&V 40": **not found** by web search (the phrase appears in company bios and standards pages, not in advertisements).

### Q4. Where is non-Newtonian rheology a stated problem in regulatory CFD?

- FDA nozzle, non-Newtonian effects (PLOS ONE, 25 Mar 2014): Newtonian vs Carreau-Yasuda differences up to 10 % in axial velocity, 3-4 % in pressure drop, and an "8 % underestimation of normalized index of hemolysis" by the Newtonian model — "almost half of the improvement" obtained by doing CFD at all rather than analytical estimates; "non-Newtonian effects are actually amplified in blood damage indices"; "no significant additional computational cost": https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0092638
- Same benchmark, later work: "Significant prediction variations in blood damage (by nearly a factor of six) among the Non-Newtonian models have been observed at low Reynolds numbers"; "models with higher viscosity at low shear rates (e.g., Casson, Carreau-Yasuda) predicted elevated hemolysis index values" (University of North Dakota thesis; quoted from search abstract, full text not read — **[unverified wording]**): https://commons.und.edu/theses/7509/
- FDA multi-laboratory hemolysis experiment on the nozzle: blood viscosity measured "each test day at a shear rate of 200 s-1 using a cone-and-plate rheometer" at 36 % haematocrit, i.e. a single-rate Newtonian proxy; interlaboratory variability 54-85 %, intralaboratory "47 to 94 %"; "CFD predictions of blood damage are heavily dependent on which empirical shear-stress exposure time model is used and how it is implemented": https://pmc.ncbi.nlm.nih.gov/articles/PMC4381842/
- FDA blood-pump interlaboratory CFD study (Annals of Biomedical Engineering, 2022): 24 submissions (14 steady, 10 transient); "not a strong correlation between mesh resolution or turbulence model with CFD accuracy"; need to validate "across the entire range of operating conditions and for all quantities of interest"; only 37 % of submissions included hemolysis predictions: https://link.springer.com/article/10.1007/s10439-022-03105-w ; ASAIO abstract: https://asaio.org/Conference/program/2022/BIO10.cgi
- FDA's own worked V&V 40 example is hemolysis in a centrifugal blood pump (Morrison et al., ASAIO J 2019): https://pubmed.ncbi.nlm.nih.gov/30973403/
- FDA Regulatory Science Tool "Hemolysis Dataset for Validating CFD Simulation through Generalized Cardiovascular Medical Device Geometries", RST26CV02.01, dated 5 May 2026 (for hemodialysis catheters, ECMO cannulae, VADs, bypass pumps); the catalogue entry makes no mention of viscosity or non-Newtonian properties: https://cdrh-rst.fda.gov/hemolysis-dataset-validating-computational-fluid-dynamic-cfd-simulation-through-generalized ; raw data: https://github.com/OSEL-DAM/CFD-and-Blood-Damage-Benchmarks
- Hemolysis-model uncertainty: UQ of power-law hemolysis models (arXiv 2024): https://arxiv.org/pdf/2407.18757 ; dependence on the empirical coefficient set: https://pubmed.ncbi.nlm.nih.gov/39177925/
- Standards: ISO 14708-5:2020 Annex C covers in vitro and in silico evaluation of circulatory support devices: https://www.iso.org/standard/69898.html ; ISO/AWI 26302 (in progress): https://www.iso.org/standard/90299.html
- Inhaled/nasal (pharma): FDA GDUFA regulatory-science reports state CFD has been used "to predict regional drug deposition for inhalers and nasal sprays" and that "for nasal suspension drug products, in silico models can predict the influence of device and formulation differences on regional deposition and absorption": https://www.fda.gov/media/166594/download ; FY2018 report: https://www.fda.gov/media/130621/download ; FY2025 workshop: https://www.fda.gov/media/188015/download and https://www.fda.gov/media/188040/download ; CFD-PK for nasal suspensions: https://www.sciencedirect.com/science/article/abs/pii/S0378517324008949 . Mucus rheology as a *stated* regulatory-CFD problem in FDA text: **not found**.
- Dissolution hydrodynamics: CFD of the PhEur/USP disintegration tester with Newtonian and non-Newtonian fed-state media found a 1 Pa s medium gives a "10-fold increase in shear stress on a tablet" and that the pharmacopoeial apparatus does not reproduce in vivo hydrodynamics: https://pubmed.ncbi.nlm.nih.gov/26017815/ ; USP basket apparatus CFD: https://www.sciencedirect.com/science/article/abs/pii/S0022354919307506 . A regulatory requirement to model dissolution hydrodynamics: **not found**. Viscous long-acting injectables as a regulatory CFD topic: **not found**.
- Practitioner statements on which blood-viscosity model is used in actual submissions: **not found** (submission content is confidential; SSED PDFs found do not describe viscosity models).

### Q5. Standards, benchmarks, tooling

- ASME V&V 40-2018 (FDA-recognised consensus standard): https://www.asme.org/codes-standards/find-codes-standards/assessing-credibility-of-computational-modeling-through-verification-and-validation-application-to-medical-devices ; committee page (VVUQ 40): https://cstools.asme.org/csconnect/CommitteePages.cfm?Committee=100108782&Action=42086 ; VVUQ overview and annual symposium: https://www.asme.org/codes-standards/publications-information/verification-validation-uncertainty
- Pipeline (DMD 2025 session): "ASME VVUQ 40.1 Technical Report", an end-to-end tibial-tray example, was expected early 2025; a patient-specific-modelling working group is active: https://dmd.umn.edu/2025/vv40 . FDA's Jan 2024 deck lists an in silico clinical trials playbook and a "CM&S mock submission" as ongoing, a "credibility assessment toolbox" and "detailed credibility assessment examples" as possible future efforts, and states "There is a need for future field/device-specific prescriptive standards": https://ntp.niehs.nih.gov/sites/default/files/2024-02/2_Pathmanathan_ICCVAM-CoP-2024_508.pdf
- FDA Regulatory Science Tools catalogue: https://cdrh-rst.fda.gov/ ; OSEL credibility research programme: https://www.fda.gov/medical-devices/medical-device-regulatory-science-research-programs-conducted-osel/credibility-computational-models-program-research-computational-models-and-simulation-associated
- MDDT programme (voluntary, no fee): https://www.fda.gov/medical-devices/medical-device-development-tools-mddt . Recent qualifications listed on that page: WOUND-Q (8 Nov 2024), numaScrew virtual pullout (10 Apr 2025), MolecuLightDX (26 Jan 2026), IVIES (3 Feb 2026); an MRI RF-heating FEM tool is also qualified. Overview: https://meddeviceguide.com/blog/fda-medical-device-development-tools-mddt-program-guide
- EU: Team-NB position paper on MDR Article 61(10) (Apr 2026): https://www.team-nb.org/wp-content/uploads/2026/04/Team-NB-PositionPaper-Article61-10-V1-20260421.pdf ; the 2025 Commission proposal to simplify MDR/IVDR references new approach methodologies including in silico: https://healthcarelifesciences.bakermckenzie.com/2025/12/19/the-eus-2025-proposal-to-simplify-the-medical-and-in-vitro-diagnostic-devices-regulations-mdr-ivdr/
- Academic quality frameworks: "Ensuring the quality of in silico evidence" (Briefings in Bioinformatics, 2026): https://academic.oup.com/bib/article/27/4/bbag442/8762911 ; systematic review of regulatory adoption of in silico trials (TIRS 2025): https://link.springer.com/article/10.1007/s43441-025-00871-2
- Evidence-package software: Ansys Minerva "template-based, end-to-end workflow" for V&V 40 credibility documentation **[vendor marketing; page returned 403 on fetch, text from search snippet]**: https://ansys.synopsys.com/blog/ansys-minerva-streamlines-credibility-assessment-for-healthcare-in-silico-testing ; Virtonomy automated "FDA 1807 & ASME V&V40" reports: https://www.virtonomy.io/product/ . A standalone, solver-independent "V&V 40 manifest" or evidence-assembly SaaS: **not found**.
- MDIC CM&S Reference Hub and pre-symposium training deck: https://mdic.org/resources/mdic-cm_s-reference-hub/ ; https://mdic.org/wp-content/uploads/2024/06/MDIC-Pre-Symposium-2.pdf

### Q6. Price points

- Public price for a regulatory CFD study, credibility assessment or in silico trial: **not found** (Stress Engineering, Veryst, Virtonomy, Thornton Tomasetti, SimuTech pages carry no prices).
- Grant ceilings: NIBIB SBIR/STTR Phase I "up to $256K for 6-12 months", Phase II "up to $1.7M for up to 2 years": https://www.nibib.nih.gov/programs/small-business-programs/project-structure
- Only CFD-hemolysis SBIR award surfaced: Research South Inc., NIH Phase I, 2002, $98,322, "Modern Design Software for Biomedical Devices" (CFD shape optimisation to minimise hemolysis): https://www.sbir.gov/awards/99339 . Awards 2022-2026 on this topic: **not found** by web search; sbir.gov advanced search is needed.
- Start-up funding: Adsilico $4.4M (Feb 2024); ELEM $2.27M total to 2023 (Q1 URLs).
- Analyst market sizing for "in silico clinical trials" **[analyst PR, unaudited]**: $6.68B by 2032 (SNS Insider, Nov 2025): https://www.globenewswire.com/news-release/2025/11/28/3195999/0/en/In-Silico-Clinical-Trials-Market-Size-to-Reach-6-68-Billion-by-2032-SNS-Insider.html ; $3.97B (2025) to $8.51B (2035) (July 2026 release, GlobeNewswire). These are unrelated to device-CFD subcontracting revenue and should not be used for sizing.
- Labour-cost proxy: $46.88-$73.56/h employee pay for "computational modeling" (ZipRecruiter, Aug 2026): https://www.ziprecruiter.com/Jobs/Computational-Modeling

### Q7. Insurance / liability / QMS barriers

- ISO 13485 supplier control is risk-proportionate: organisations must apply controls to "suppliers, contractors, and consultants based on the effect their products or services can have on the quality of the medical device" (Johner Institute): https://blog.johner-institute.com/quality-management-iso-13485/supplier-evaluation-supplier-selection-supplier-audits/ ; practical guide (6 May 2026): https://avendium.com/2026/05/06/supplier-qualification-under-iso-13485-a-practical-guide/
- Simulation vendors that advertise ISO 13485 certification as a selling point: PADT, Simq, Synopsys Simpleware (Q1). Stress Engineering's V&V 40 procedure sits inside a QMS document hierarchy (L1-QMS-PR-09 Deliverables, L2-CIN-PR-08 Test and Measurement, L3-CIN-GD-13 Good Documentation Practices), which is what a sponsor's supplier audit would look at: https://sesforce4.com/wp-content/uploads/2024/04/L2-CIN-PR-22_Rev0-VV40.pdf
- Professional liability: E&O coverage is routinely required by contract "with specified amounts per occurrence ... with annual aggregate limits" for engineering services (Oregon DAS clause page): https://www.oregon.gov/das/risk/pages/insclausesproferroromis.aspx ; general engineer PL coverage descriptions: https://www.thehartford.com/professional-liability-insurance/engineers ; https://www.insureon.com/consulting-business-insurance/professional-liability
- Device-sponsor-specific requirements for a simulation subcontractor (required limits, audit expectations, a case where a solo shop was refused): **not found** in public text or forums.
- Boston Engineering's page on what sponsors expect from outsourced development partners (design-control familiarity): https://blog.boston-engineering.com/what-medical-device-teams-should-know-before-hiring

## Named sellers table

| Name | Size (evidence) | Offering | Has bench lab? |
|---|---|---|---|
| Stress Engineering Services | Large multi-office US firm (headcount not quoted) | V&V 40-procedure-based CM&S, in silico DVT, FDA/EU MDR submissions [vendor] | Has a "Test and Measurement" procedure; fluids/blood lab not confirmed |
| Thornton Tomasetti Applied Science | 144 open jobs firm-wide (2026) | V&V 40 through submission; CFD; benchtop validation design [vendor] | Describes benchtop validation; ownership not confirmed |
| Veryst Engineering | ~18 staff (LeadIQ, Jul 2026) | FEA/CFD for 510(k)/PMA; FDA-nozzle hemolysis case study [vendor] | Not confirmed |
| Exponent | Public company | FEA, biofluid dynamics incl. hemolysis; regulatory consulting [vendor] | Corporate physical labs exist; blood flow loop not confirmed |
| SimuTech Group | Ansys channel partner | V&V 40-framed CFD/FEA; committee member on staff [vendor] | Not found |
| Numalogics | Small (Canada, founded 2010) | FEA apps with embedded ASTM/ISO tests; MDDT-qualified pullout model | Partner Sawbones supplies the physical comparator |
| Virtonomy | Start-up (Munich, 2019) | v-Patients SaaS, auto V&V 40 / FDA-1807 reports, simulated ISO bench tests [vendor] | No (simulated bench tests) |
| ELEM Biotech | Seed, $2.27M (Tracxn) | Alya virtual humans, cardiac + inhaler | No |
| Adsilico | $4.4M raised (Feb 2024) | In silico trials for devices | No |
| InSilicoTrials | Not found | Cloud in silico trial platform | No |
| Ansys / Synopsys | Large | Minerva V&V 40 templates; Simpleware Medical 510(k) | No |
| PADT | SME, ISO 13485 | FEA/CFD services [vendor] | Not found |
| Simq GmbH | SME, ISO 13485 | Simulation software/services | Not found |
| CardioMed Device Consultants | Boutique, founded 2003 (Tentamus) | Regulatory/pre-clinical consulting; panelist on FDA credibility case studies | Not found |
| Nelson Labs (contract lab, not a modeller) | Large | ISO 10993-4 mechanically induced hemolysis testing | Yes |

"Sells" vs "buys": every row above sells. Documented buyers are the two Virtonomy customers (Medira, Scandinavian RealHeart), Sawbones (Numalogics), and the unnamed sponsors behind SimuTech's and Stress Engineering's claims. The MDIC survey is the only quantitative buyer-side source, and it says most buyers self-serve.

## What a 1-3 person shop can realistically win here

- **Subcontracted CFD plus credibility-plan authoring for sponsors without in-house CM&S** (the 28 % of MDIC respondents without internal resources, and the start-ups FDA itself says are underserved). The deliverable structure is public (Stress Engineering's procedure), the guidance is non-prescriptive, and the top buyer barrier is "what does the regulator expect", so authoring the QOI/COU/risk/credibility-factor plan is a knowledge product rather than a headcount product.
- **Blood-contacting CFD where rheology is the credibility gap.** The FDA nozzle literature quantifies the Newtonian error (8 % NIH, ~6x spread among models at low Re), the multi-lab experiment characterised blood by one viscosity at 200 1/s, and FDA's May 2026 hemolysis RST dataset says nothing about viscosity. A shop that measures the full viscosity curve of the actual test blood or analogue and propagates it through the model, with uncertainty, does something the benchmarks themselves did not.
- **Physical comparators for platform companies** (Virtonomy, Adsilico, ELEM) that sell simulated bench tests and own no lab; the "validation" column of V&V 40 needs measured data, and FDA states the standard "implicitly assumes" prospective bench tests.
- **Inhaled/nasal deposition and dissolution hydrodynamics** for generics (GDUFA-funded science where FDA already uses CFD). No evidence that non-Newtonian mucus or fed-state media is yet a regulatory ask, so this is a research-led sale, not a compliance-led one.
- **Not realistic alone**: full in silico clinical trials, an MDDT qualification (multi-year; Numalogics needed a physical-comparator partner), or competing with Exponent/Thornton Tomasetti on Class III PMA programmes.

## Evidence AGAINST

- Buyers mostly self-serve: 72 % have internal resources; only 23 % use consultants when required (2021 MDIC survey).
- FDA states the standard "implicitly assumes validation against prospective well-controlled bench tests" (Pathmanathan, Jan 2024): a simulation-only subcontractor without a lab must partner or buy tests, and a bench-only one is competing with Nelson Labs.
- Hemolysis validation is limited by 47-94 % intra-lab experimental scatter; a better viscosity model cannot buy credibility that the comparator cannot resolve.
- EU notified bodies still favour standard-based bench testing; the MDR has no credibility criteria; no EU growth number exists.
- No public price anywhere; the sales cycle is relationship-driven (panels, committees), and every incumbent already uses V&V 40 vocabulary.
- FDA's 2024 symposium lists "Uncertainty in CM&S Expectations by Regulatory Bodies" as the top impediment: the regulator itself has not converged, so a "credibility package" product may keep being re-scoped; FDA also says prescriptive device-specific standards are still needed.
- No annual FDA count of CM&S submissions was found: the growth story rests on a 2017 consult survey, a 2022 "scores" phrase, and a 2024 opinion poll.
- The only sizing numbers are analyst press releases for "in silico clinical trials", which are not this market.

## Open questions requiring a human conversation

1. Ask an FDA OSEL reviewer or an ASME VVUQ 40 member (SimuTech's Goodin, Numalogics, Medtronic's Bodner) what fraction of blood-contacting CFD submissions they see using a non-Newtonian viscosity model, and whether viscosity choice has ever appeared in a deficiency letter.
2. Ask two small Class II cardiovascular sponsors what they paid for outsourced CFD plus a credibility plan, how long the engagement ran, and what supplier-qualification paperwork and insurance limits they required.
3. Ask Virtonomy, Adsilico or ELEM whether they buy physical comparator testing, from whom, and at what price.
4. Obtain the gated MDIC Landscape Report PDF (form-gated) for the breakdown of outsourcing by company size.
5. Query sbir.gov advanced search for 2022-2026 awards mentioning hemolysis, CFD, V&V 40 or in silico device testing (web search surfaced only a 2002 award).
6. Confirm whether ASME VVUQ 40.1 was published in 2025 and whether a V&V 40 revision is scheduled.
7. Ask a contract hemocompatibility lab (Nelson Labs or equivalent) whether they are ever asked to supply rheometry of the test blood alongside hemolysis results.

## Sources

- https://www.federalregister.gov/documents/2023/11/17/2023-25470/assessing-the-credibility-of-computational-modeling-and-simulation-in-medical-device-submissions
- https://www.fda.gov/media/154985/download
- https://www.fda.gov/medical-devices/medical-devices-news-and-events/webinar-assessing-credibility-computational-modeling-and-simulation-medical-device-submissions-final
- https://www.fda.gov/media/163156/download
- https://ntp.niehs.nih.gov/sites/default/files/2024-02/2_Pathmanathan_ICCVAM-CoP-2024_508.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/PMC6167449/
- https://iopscience.iop.org/article/10.1088/2516-1091/ae1c05
- https://pubmed.ncbi.nlm.nih.gov/41192059/
- https://mdic.org/resources/2024-fda-mdic-symposium-on-computational-modeling-and-simulation-presentations/
- https://mdic.org/wp-content/uploads/2024/06/MDIC-Pre-Symposium-2.pdf
- https://medtechintelligence.com/news_article/medical-device-computational-modeling-and-simulation-reports-highlights-best-uses-and-barriers-to-adoption/
- https://mdic.org/resources/landscape-report-industry-survey-on-the-use-of-computational-modeling-simulation-in-medical-device-development/
- https://mdic.org/resources/mdic-cm_s-reference-hub/
- https://starfishmedical.com/resource/computational-modelling-simulation-medical-device-development/
- https://starfishmedical.com/resource/fda-guidance-on-cms-in-medical-device-submissions/
- https://www.feops.com/src/Frontend/Files/MediaLibrary/08/morrison-j-med-devices-11-2017-role-of-computational-modeling.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12579633/
- https://link.springer.com/article/10.1007/s43441-025-00871-2
- https://www.drugdiscoverynews.com/how-modeling-and-simulation-are-reshaping-regulatory-science-16617
- https://sesforce4.com/wp-content/uploads/2024/04/L2-CIN-PR-22_Rev0-VV40.pdf
- https://www.stress.com/in-silico-design-verification-testing-of-medical-devices/
- https://www.stress.com/fda-cranks-up-pressure-for-modeling-and-simulation/
- https://www.stress.com/services/medical-pharmaceutical/analysis-simulation/
- https://www.thorntontomasetti.com/news/fda-guidance-computational-modeling-medical-device-submissions
- https://www.thorntontomasetti.com/solution/biomechanics-life-sciences
- https://www.thorntontomasetti.com/capability/computational-fluid-dynamics
- https://www.ziprecruiter.com/Jobs/Thornton-Tomasetti
- https://www.veryst.com/services/specialized-expertise/medical-devices
- https://www.veryst.com/services/simulation-analysis/computational-fluid-dynamics-cfd
- https://leadiq.com/c/veryst-engineering-llc/5a1d83d424000024005e5717
- https://www.exponent.com/capabilities/medical-devices
- https://www.exponent.com/capabilities/finite-element-analysis-medical-devices-biomaterials
- https://www.exponent.com/capabilities/biofluid-dynamics-biomedical-devices-products
- https://go.simutechgroup.com/medical-device-1
- https://simutechgroup.com/asme-vv-40-standard-medical-devices/
- https://simutechgroup.com/industry-solutions/healthcare/
- https://numalogics.com/wp-content/uploads/2025/05/PR-Screw-Pullout-MDDT-2025-05-20.pdf
- https://sawbones.com/blog/firstever-fda-qualification-of-an-orthopedic-device-mechanical-test-under-mddt-program/
- https://www.ansys.com/blog/democratizing-in-silico-medical-device-testing
- https://numalogics.com/careers/finite-element-simulation-engineer/
- https://www.virtonomy.io/product/
- https://medium.com/pace-ventures/why-we-invested-in-virtonomy-53c0f14aa273
- https://tracxn.com/d/companies/elem-biotech/__hHHepAcMcrWoxjD_nFg2gsdy3Qr8LMc1TPCapR2vsfg
- https://www.elem.bio/
- https://www.medicaldevice-network.com/news/uks-adsilico-raises-4-4m-to-develop-virtual-trial-solution/
- https://mtec-sc.org/life-sciences/insilicotrials
- https://www.cbinsights.com/compare/insilico-trials-vs-virtonomy
- https://ansys.synopsys.com/blog/ansys-minerva-streamlines-credibility-assessment-for-healthcare-in-silico-testing
- https://www.synopsys.com/simpleware/resources/regulatory.html
- https://padtinc.com/blog/the-focus/fda-opening-to-simulation-supported-verification-and-validation-for-medical-devices
- https://simq.de/en/software-development/
- https://www.cardiomedllc.com/
- https://www.tentamus.com/locations/cardiomed-device-consultants-llc/
- https://www.nelsonlabs.com/mechanically-induced-hemolysis-iso-10993-4/
- https://www.jobzmall.com/medtronic/job/computational-modeling-engineer
- https://builtin.com/job/senior-product-development-engineer-computational-modeling/3757574
- https://dmd.umn.edu/2025/vv40
- https://www.ziprecruiter.com/Jobs/Computational-Modeling
- https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0092638
- https://commons.und.edu/theses/7509/
- https://pmc.ncbi.nlm.nih.gov/articles/PMC4381842/
- https://link.springer.com/article/10.1007/s10439-022-03105-w
- https://asaio.org/Conference/program/2022/BIO10.cgi
- https://pubmed.ncbi.nlm.nih.gov/30973403/
- https://cdrh-rst.fda.gov/hemolysis-dataset-validating-computational-fluid-dynamic-cfd-simulation-through-generalized
- https://github.com/OSEL-DAM/CFD-and-Blood-Damage-Benchmarks
- https://arxiv.org/pdf/2407.18757
- https://pubmed.ncbi.nlm.nih.gov/39177925/
- https://www.iso.org/standard/69898.html
- https://www.iso.org/standard/90299.html
- https://www.fda.gov/media/166594/download
- https://www.fda.gov/media/130621/download
- https://www.fda.gov/media/188015/download
- https://www.fda.gov/media/188040/download
- https://www.sciencedirect.com/science/article/abs/pii/S0378517324008949
- https://pubmed.ncbi.nlm.nih.gov/26017815/
- https://www.sciencedirect.com/science/article/abs/pii/S0022354919307506
- https://www.asme.org/codes-standards/find-codes-standards/assessing-credibility-of-computational-modeling-through-verification-and-validation-application-to-medical-devices
- https://cstools.asme.org/csconnect/CommitteePages.cfm?Committee=100108782&Action=42086
- https://www.asme.org/codes-standards/publications-information/verification-validation-uncertainty
- https://cdrh-rst.fda.gov/
- https://www.fda.gov/medical-devices/medical-device-regulatory-science-research-programs-conducted-osel/credibility-computational-models-program-research-computational-models-and-simulation-associated
- https://www.fda.gov/medical-devices/medical-device-development-tools-mddt
- https://meddeviceguide.com/blog/fda-medical-device-development-tools-mddt-program-guide
- https://www.team-nb.org/wp-content/uploads/2026/04/Team-NB-PositionPaper-Article61-10-V1-20260421.pdf
- https://healthcarelifesciences.bakermckenzie.com/2025/12/19/the-eus-2025-proposal-to-simplify-the-medical-and-in-vitro-diagnostic-devices-regulations-mdr-ivdr/
- https://academic.oup.com/bib/article/27/4/bbag442/8762911
- https://www.nibib.nih.gov/programs/small-business-programs/project-structure
- https://www.sbir.gov/awards/99339
- https://www.globenewswire.com/news-release/2025/11/28/3195999/0/en/In-Silico-Clinical-Trials-Market-Size-to-Reach-6-68-Billion-by-2032-SNS-Insider.html
- https://blog.johner-institute.com/quality-management-iso-13485/supplier-evaluation-supplier-selection-supplier-audits/
- https://avendium.com/2026/05/06/supplier-qualification-under-iso-13485-a-practical-guide/
- https://www.oregon.gov/das/risk/pages/insclausesproferroromis.aspx
- https://www.thehartford.com/professional-liability-insurance/engineers
- https://www.insureon.com/consulting-business-insurance/professional-liability
- https://blog.boston-engineering.com/what-medical-device-teams-should-know-before-hiring
