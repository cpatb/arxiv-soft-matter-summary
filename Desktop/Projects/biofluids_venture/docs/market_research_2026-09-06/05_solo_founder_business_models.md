# 05 — How small physics / CFD / rheology companies actually become viable

Research date: 2026-09-06/07. Method: ~35 web searches and ~20 page fetches. Every
factual claim carries a URL; where a page was fetched the number is quoted from that page.
Where I could not find evidence I say "not found". Search budget was rate-limited twice,
so several items (offshore CFD rate cards, the CFD-online pricing thread, Upwork profile
rates, Robin Knowles's about page) are marked as blocked (HTTP 403/404) rather than absent.

Context for the reader: the founder is a finishing PhD (dense suspensions, drop impact,
rheology, Basilisk/OpenFOAM DNS, agentic software), 5–10 h/week now, Indiana/Illinois,
target: full-time founder income ($150–300k/yr) in a 1–5 person niche technical company.

---

## Summary

1. Every small CFD/rheology company with a documented history started by **selling hours** (consulting, contract testing, training) and only later, if ever, moved to product. Flow Science (1980→1984), Convergent Science (1997→2008), CFD Support (2009→~2014), CoolSim (2005→2008) all followed that path; none went product-first.
2. **Time to a viable consultancy was 1.5–5 years** where it is documented: CFD Support's first paid project came 18 months in, after a year of cold calls that "no one ever called back" and a free pilot that ate 400+ hours; Flow Science left pure contract work after 4 years because "a small company cannot easily exist on contract work".
3. The rheology side has a **stable, low-ceiling, solo-viable pattern**: one PhD + one rheometer + contract testing + training (Cunningham/Centre for Industrial Rheology since 1998/2012, Gammadot since 2006, Rheology Testing Services, solo, Chapel Hill). None publishes revenue; all are still trading.
4. **Pricing floor and ceiling are both visible.** University core rheometers rent to industry at $50–64/h instrument time plus ~$50–55/h staff; per-sample contract rheology is ~$170–340; US CFD consulting quotes are $80–250/h and $1.5k–20k per study; IEEE-USA's survey put PhD consultants' median at $200/h; specialist day rates are quoted at $1,200–3,000+.
5. **Offshore pressure is real but not quantified for CFD specifically.** India-based CFD shops publish mesh-count pricing with $300 (2D) / $500 (3D) minimums; general Indian engineering outsourcing runs $25–90/h. Upwork rate pages were blocked (403) so per-profile rates are not found.
6. **Non-dilutive money is large relative to the target income but slow and gated.** NSF SBIR Phase I is now $305k (6–18 months), requires a Project Pitch invitation, and the end-to-end pitch-to-award rate is ~9–12%. NIH Phase I caps are $323k–$700k by institute. Indiana matches $0.50/federal dollar up to $50k; Illinois matches up to $75k on Phase I. The PI must be >50% employed by the company at award, not at proposal.
7. **The one channel with repeated evidence is teaching**: webinars were CoolSim's strongest acquisition channel; CFD Direct delivered 81 training days and sold 3,500 books in year 10; Cunningham has trained "an estimated three thousand" people. Cold outreach has the only documented failure (CFD Support, year 1).
8. **Nobody in the sample got rich on open-source maintenance alone.** OpenFOAM Foundation maintenance funding was ~€250k/yr against a stated need of €500k; OpenFOAM Ltd (ESI) was 5 people and £874k revenue in 2020/21.
9. **AI-era leverage evidence is anecdotal, not financial.** There are 2026 case studies of physicists building simulation software with Claude Code in days (arXiv 2605.30353; Medium/Towards AI posts), but no solo founder with disclosed revenue from an agent-run simulation business was found. The funded "AI-native simulation" players (Quanscient, €10M Series A, 40–50 staff) are not solo.
10. **Net for this founder:** the evidence supports contract rheology/simulation studies sold at fixed project prices, a training/webinar channel from day one, and an SBIR Project Pitch as soon as a company exists — and it argues against product-first, hourly billing, and open-core as the first revenue.

---

## Evidence by question

### Q1. Case studies: what did they sell first, how long to full-time, what killed the ones that died?

**CFD Direct (UK, 2015).** Launched 25 March 2015 by Henry Weller, Chris Greenshields and Jenya Collings; revenue streams are training, cloud (CFDDFC on AWS), support/consultancy and OpenFOAM development. Year-10 (2024–25) numbers from their own page: 81 training days, 1,350 AWS cloud subscribers, 3,500 copies of their $29.95 book sold, 961 development hours and 1,415 maintenance hours. Maintenance "required designated funding" so they launched a funding campaign in 2017 and an Open Day in 2018. Team size not stated; two named hires in years 3 and 9. (https://cfd.direct/about/year-10-2025/ ; https://cfd.direct/about/)

**Wikki Ltd (UK, 2004).** Founded by Hrvoje Jasak when FOAM went open source; predecessor Nabla Ltd (2000–04) was funded mainly by Fluent Inc. and wound up in 2006. Henrik Rusche joined as partner in 2008. Revenue and headcount not found. (https://wikki.co.uk/index.php/about-us/)

**Engys (UK, 2009).** Privately held; HELYX-Core developed since 2010 on OpenFOAM; "hundreds of consultancy projects"; offices UK/Germany/Italy/USA/Australia. Revenue and headcount not found. (https://old.engys.com/about-us ; https://engys.com/)

**Flow Science (USA, 1980).** Hirt left Los Alamos in 1980 and "initially performed contract work"; in 1984 switched to selling FLOW-3D licenses because "a small company cannot easily exist on contract work because it either has too much work and not enough workers, or too many workers and not enough work." Chose a niche (free-surface flows) deliberately. (https://www.flow3d.com/a-life-in-cfd-by-dr-c-w-tony-hirt/ ; https://www.flow3d.com/flow-science/)

**Convergent Science (USA, 1997).** Founded by UW-Madison graduate students; "initially ... a CFD consulting company"; founders noticed most consulting time went into meshing and built the autonomous-meshing product; **first CONVERGE licenses sold in 2008, eleven years after founding**, during the financial crisis. Later 60 staff in Madison plus 20 elsewhere. (https://en.wikipedia.org/wiki/Convergent_Science ; https://news.wisc.edu/engine-software-from-uw-spinoff-being-used-around-the-world/)

**CFD Support (Czechia, 2009).** Two fresh PhDs from the von Karman Institute. Founder's own account (LinkedIn, Oct 2020): year 1 "complete failure" — cold calling produced dozens of meetings and "no one ever called back"; free pilots (one 400+ hours) did not convert; first paid project (a steam-turbine stage) came ~18 months in and made subsequent sales "easier and easier"; consultancy viable years 2–5; then pivoted to product (TCFD→TCAE) because raw-OpenFOAM consulting "didn't scale"; 30+ man-years invested over seven years. No revenue figures. (https://www.linkedin.com/pulse/cfd-support-another-dream-come-true-lubos-pirkl?articleId=6724652159088001024 ; https://www.cfdsupport.com/about-us/)

**Applied CCM (Australia 2011 / Canada 2015).** Founded by Darrin Stephens (Melbourne) and Chris Sideroff (Ottawa); OpenFOAM support, consulting, and the Caelus fork. Podcast episode page (Dec 2016) carried no financial detail. Viability timeline not found. (https://www.appliedccm.com/about/ ; https://cfdengine.com/podcast/darrin-stephens-chris-sideroff-caelus/)

**CFD Engine (UK, solo).** Robin Knowles, one-person OpenFOAM/AWS consultancy for motorsport aero, sells **fixed-price "done-for-you" aero programmes** (not hourly), a paid course ("Beginner's Guide to OpenFOAM on AWS"), a newsletter and the Talking CFD podcast. "7+ years" of solo delivery claimed on LinkedIn. Prices, subscriber counts and revenue not published; about page returned 404. (https://cfdengine.com/ ; https://uk.linkedin.com/in/robinknowles ; https://learn.cfdengine.com/)

**CoolSim / Applied Math Modeling (USA, 2008).** Team left Ansys/Fluent in 2008 to sell a niche SaaS (data-centre airflow) on top of Fluent; founder's stated lessons: niche domain expertise, client-server SaaS for support, and **webinars as one of the strongest acquisition channels**. (https://www.coolsimsoftware.com/about/ ; https://cfdengine.com/podcast/paul-bemis-coolsim/)

**Cascade Technologies (USA).** Stanford CTR spin-off (CharLES LES code), SBIR-funded, acquired by Cadence in Nov 2022. Founding date and early revenue not found. (https://community.cadence.com/cadence_blogs_8/b/cfd/posts/cadence-welcomes-cascade-technologies-1221740980 ; https://www.sbir.gov/portfolio/118337)

**Veryst Engineering (USA, 2005).** Founded by Stuart Brown; product design/manufacturing/failure analysis consulting; heavy COMSOL partner. The "Exponent spinout" premise in the brief is **not confirmed** by anything found. (https://www.veryst.com/about-us ; https://www.veryst.com/about/our-team/stuart-brown)

**Fluid Codes (UAE, 1997).** Ansys channel partner (reseller + training + consulting), four offices. It is a reseller model, not an independent consultancy. (https://fluidcodes.com/about/)

**Rheology Solutions (Australia, 1998).** Instrument sales + service + training + contract testing; owner Pat Griffin; now trading as Online Rheometer Pty Ltd. (https://www.labonline.com.au/directories/rheology-solutions-pty-ltd-422767475)

**Centre for Industrial Rheology (UK, 2012) / Neil Cunningham.** Sales engineer at TA/Rheometric 1995–98, independent consultant since 1998, lab founded 2012; testing, consulting, training; claims "many hundreds of companies" and ~3,000 people trained; named clients GSK/Haleon, Pfizer, Estée Lauder, Chanel. (https://www.rheologylab.com/about/ ; https://sites.rutgers.edu/centerfordermalresearch/wp-content/uploads/sites/271/2022/07/Neil-Cunning-bio_photo_abstract.pdf)

**Gammadot Rheology (UK, 2006).** Mark Edwards, 18 years at Rapra before going independent; testing + failure analysis consultancy. No prices published. (https://www.gammadot.com/About_Us/gammadot.htm)

**Rheology Testing Services (USA, Chapel Hill).** Mark Patrick PhD, 24 years at GSK; solo; assay menu (flow, yield, thixotropy, sweeps, tack, tribology); pricing not disclosed; offers a free intro presentation as lead magnet. (https://www.rheologytestingservices.com/rheology-testing-overview)

**RheoSense (USA, 2001).** Instrument company (VROC chip viscometer), founder Seonggi Baek; only $780k total disclosed funding per CB Insights — an instrument company grown on revenue, not VC. (https://www.rheosense.com/about-rheosense ; https://www.cbinsights.com/company/rheosense)

**Formulaction (France, 1994).** Turbiscan inventor; 36 employees; acquired by Verder/Microtrac in May 2023. (https://www.verder-scientific.com/news/news/details/news_detail/verder-group-acquires-formulaction ; https://www.oaklins.com/nl/en/deals/109018/)

**What killed the ones that died:** the only documented death in the sample is Nabla Ltd (2000–06), which was funded by Fluent and wound up after FOAM went open source (https://wikki.co.uk/index.php/about-us/). The rest of the sample survives, which is a survivorship bias: dead solo consultancies do not leave web pages. **Not found:** any post-mortem of a failed solo CFD or rheology consultancy.

### Q2. Pricing benchmarks

| Item | Number | Source / date |
|---|---|---|
| University rheometer, external commercial | $63.90/h instrument + $55/h staff (UF ARES LS1); $50/h (UW Anton Paar MCR-301) | https://nrf.aux.eng.ufl.edu/resources/resource.asp?id=105 ; https://www.cheme.washington.edu/research/facilities/bindra/rheometer |
| Contract rheology per sample | $170 basic curve, $340 advanced (UW Clean Energy Testbeds) | https://www.wcet.washington.edu/per-sample-analysis/ |
| US CFD consulting hourly | $80–120/h (Cad Crowd, Feb 2023); $100–200/h co-sourced, $100–250/h outsourced (Resolved Analytics, Durham NC) | https://www.cadcrowd.com/blog/what-are-cfd-engineering-rates-consulting-services-costs-company-pricing/ ; https://www.resolvedanalytics.com/ebook/cfd-costs |
| US CFD per-study | $1,500–4,500 per analysis; $2,000 simple 2D to ~$20,000 complex transient (Cad Crowd); $2,500–9,500 (CFD Source) | same ; https://www.cfdsource.com/ShowBlog/whats-the-real-cost-of-cfd-consulting/ |
| UK fixed-price CFD packages | £1,990 / £3,990 / £5,990 by complexity (8020 Engineering) | https://8020engineering.com/consulting-pricing-and-packages/ |
| Loaded cost of an in-house CFD analyst | ~$95k salary, ~$190k loaded (2.0x) | https://www.resolvedanalytics.com/ebook/cfd-costs |
| PhD consultant median hourly | $200/h (IEEE-USA fee survey, as reported by Physics Today, Feb 2023) | https://physicstoday.aip.org/news/consulting-can-be-stimulating-and-lucrative-for-physicists |
| Specialist technical day rate | $1,200–3,000+/day; senior SME $2,000–5,000+/day (secondary summary, not a survey) | https://www.quora.com/What-is-the-consulting-daily-rate-for-a-PhD-holder-with-more-than-20-years-experience |
| Freelance scientific consultants (Kolabtree, Feb 2024) | $50–150/h, biology/pharma heavy | https://www.kolabtree.com/blog/scientific-consultants-hire-freelance-guide/ |
| Generic independent consultant (Payscale, 2026) | $76.30/h mean, $25–194 range | https://www.payscale.com/research/US/Job=Independent_Consultant/Hourly_Rate |
| India CFD shop | mesh-tiered; minimum $300 (2D) / $500 (3D); tier prices only in an image | https://cfdconsultantsindia.com/cfd-consulting-rates/ |
| India engineering outsourcing, general | junior $25–45/h, senior $60–90/h; low tier $15–30/h | https://qubit-labs.com/average-hourly-rates-offshore-development-services-software-development-costs-guide/ |

Arithmetic the founder should do himself: $150k/yr at $150/h is 1,000 billable hours, i.e. ~20 h/week billed at a 50% utilisation — that is a full-time consultancy, not a side business. At $200/h and fixed-price studies of $5–10k, it is 15–30 studies a year.

**Blocked / not found:** the CFD-online "CFD consulting work pricing" thread (403), Upwork OpenFOAM profile rates (403), any Glassdoor rate card for CFD consultants, any published rheology consultancy day rate.

### Q3. Productized service vs SaaS vs open-core vs material library vs training vs grants

- **Consulting-first, product-later is the universal documented sequence** (Flow Science, Convergent, CFD Support, CoolSim, Q1). The product came out of a bottleneck seen repeatedly in consulting (meshing, workflow), not from a plan.
- **Training is a real revenue line for small OpenFOAM firms**: CFD Direct 81 training days in one year and a 3,500-copy book (https://cfd.direct/about/year-10-2025/); Engys, CFD Support and Applied CCM all list training; Cunningham's ~3,000 trainees (Q1).
- **Open-core / maintenance funding is thin**: OpenFOAM Foundation maintenance plans at €5k/€25k/€100k tiers totalled ~€250k/yr, "of the order of 0.1% of the revenue of big commercial CFD", against a stated need of €500k (https://openfoam.org/ as summarised in search; see also https://cfd.direct/openfoam/management/). OpenFOAM Ltd under ESI: 5 employees, £874k revenue (Dec 2020) (https://tracxn.com/d/companies/openfoam/__8Dfq1x_cN9UYJ8i86QHjdpm19Ts7i9Q9psm4lTlUg84).
- **Niche SaaS on top of a commercial solver** worked for CoolSim (2008–present), with webinars as the channel (https://www.coolsimsoftware.com/about/).
- **Licensing a material/rheology library:** **not found** — no small company selling a rheological material database as its primary product turned up. Closest are instrument vendors' application databases (Anton Paar, https://www.anton-paar.com/us-en/services-support/anton-paar-technical-centers/).
- **Value pricing vs hourly:** Jonathan Stark's claim of doubling income in year one by dropping hourly billing and an effective $2,000+/h is self-reported and generic, not deep-tech-specific (https://jonathanstark.com/hbin). CFD Engine's fixed-price programmes are the one CFD instance of the same idea (https://cfdengine.com/).
- **Solo indie revenue benchmarks (non-scientific):** TinyPilot $997k revenue / $236k profit in 2023 as a bootstrapped solo hardware+software business (https://mtlynch.io/solo-developer-year-6/). No comparable disclosed figures for a scientific-software solo founder were found.
- **Which reaches $150–300k soonest:** no dataset answers this for deep-tech. The rates in Q2 make consulting the only route whose unit economics reach $150k inside 12–24 months on published numbers; SaaS/product histories in Q1 took 4–11 years to first license.

### Q4. Non-dilutive funding fit (see funding table below)

- **NSF SBIR/STTR (NSF 26-510):** Phase I up to $305,000 for 6–18 months; Phase II up to $1.25M; Fast-Track up to $1.555M; Project Pitch mandatory (opened 2 June 2026, rolling); full-proposal deadlines 27 July 2026, 4 Nov 2026, 4 Mar 2027 (https://grantedai.com/news/nsf-sbir-sttr-250m-2026-05-31 ; https://www.bwcoconsulting.com/fod/nsfsbir). **NSF 26-511** is a new $40M scientific-instrumentation pilot (microfluidic platforms named) — relevant if the founder ever builds a measurement device.
- **Eligibility for a not-yet-founder:** PI must be >50% employed by the small business **at award and during the project, not at proposal**; no PhD or ownership requirement (https://seedfund.nsf.gov/solicitation-eligibility/ ; https://sbir.org/getting-started/pi-key-staff-eligibility-sbir-sttr/). So a graduating PhD can pitch before the company is fully staffed, but must be able to commit >20 h/week at award — incompatible with a 5–10 h/week side status.
- **Success rates:** Phase I full-proposal rate ~15–20% (historically 10–20%); 50–65% of pitches invited; end-to-end ~9–12% (https://www.getcada.com/insights/nsf-sbir-success-rates ; https://fundunlocked.com/learn/nsf-sbir-success-rate/ ; https://www.ebhoward.com/success-rates-of-the-national-science-foundations-sbir-sttr/).
- **NIH SBIR (PA-27-100):** Phase I caps of the SBA guideline ($323,090), $400k or $700k depending on institute; first receipt 5 Sept 2026 (https://grantedai.com/news/pa-27-100-1-4b-2026-05-31 ; https://www.ebhoward.com/nih-sbir-sttr-budget-limits-vary-more-than-most-applicants-realize/). Precedent in the founder's exact niche: QATCH's nanovisQ SBIR Phase I on viscosity for injectability/manufacturability of protein therapeutics (https://www.sbir.gov/portfolio/1168615).
- **DOE:** a $10M FY2026 "Genesis Mission" Phase I opportunity, ~40 awards, mandatory pitch by **10 Sept 2026** (https://minnesotasbir.org/federal-funding/doe-opens-fy2026-genesis-mission-sbir-sttr-phase-i-opportunity-under-new-pitch-first-process/). Whether battery-slurry rheology is one of the four areas: **not found**.
- **USDA NIFA SBIR:** Phase I up to $175,000, 8 months, TABA $6,500; last deadline 15 July 2026 has passed; agrochemical spray topic wording for FY27 **not found** (https://grantedai.com/grants/usda-small-business-innovation-research-sbir-phase-i-u-s-department-of-agriculture-usda-5918acdc ; https://www.nifa.usda.gov/application-deadline-small-business-innovation-research-program-phase-i-sbir-funding-opportunity).
- **NSF I-Corps Teams:** $50,000; entrepreneurial lead is typically a PhD student/postdoc; without NSF-award lineage you must go through a regional hub first — Indiana and Illinois are both in the Great Lakes Hub (UIUC and Purdue co-lead) (https://www.nsf.gov/funding/opportunities/nsf-national-innovation-corps-teams-nsf-national-i-corps-tm ; https://greatlakesicorps.org/programs/national-nsf-i-corps/ ; https://icorps.illinois.edu/).
- **Indiana:** SBIR/STTR match $0.50 per federal dollar up to $50k per award, Phase I and II, now run by the Applied Research Institute; Elevate Ventures $3M pre-seed fund, $20–100k cheques, no match required (https://www.insideindianabusiness.com/articles/elevate-ventures-matching-funds-for-small-tech-startups ; https://iedc.in.gov/events/news/details/2023/12/19/indiana-commits-increased-funding-resources-to-catalyze-entrepreneurship-innovation).
- **Illinois:** DCEO match up to $75k on Phase I and $250k on Phase II; must be IL-headquartered, have a Phase I within 12 months, ≥51% of work in IL; ~100 awards/yr; rolling to 30 June 2026 (https://dceo.illinois.gov/whyillinois/sbir-sttr-phase1.html ; https://grantedai.com/news/illinois-5-million-sbir-sttr-state-matching-program-2026).
- **University IP:** default is institution-owned under its IP policy (Bayh-Dole lets it retain title on federally funded inventions); typical spin-out terms ~5% or less equity plus low-single-digit royalties (https://www.ropesgray.com/en/insights/alerts/2025/02/considerations-when-contracting-with-academic-spinouts ; https://medium.com/included-vc/startups-science-ip-intellectual-property-licensing-101-for-deep-tech-founders-0685a13851fd). Purdue/Illinois/Notre Dame specific terms: **not found**. Practical point: code and know-how written on personal time and equipment, outside a sponsored project, is the cleanest asset a graduating PhD can carry out; anything in the thesis repo should be assumed university-owned until the TTO says otherwise.

### Q5. Channel realities

- **Teaching/webinars** — the only channel named as "strongest" by a founder (CoolSim, https://cfdengine.com/podcast/paul-bemis-coolsim/); CFD Direct's 81 training days and book (https://cfd.direct/about/year-10-2025/); Cunningham's trainees-to-clients funnel (https://www.rheologylab.com/about/).
- **Cold outreach** — documented failure in year 1 at CFD Support; first paid job came via a specific technical demonstration, after which "acquiring subsequent projects became easier and easier" (https://www.linkedin.com/pulse/cfd-support-another-dream-come-true-lubos-pirkl?articleId=6724652159088001024).
- **Newsletter + podcast + LinkedIn** — CFD Engine's whole front end (https://cfdengine.com/podcast/ ; https://uk.linkedin.com/in/robinknowles); general freelance surveys put LinkedIn as #1 source and 77% of freelancers getting >51% of work from repeat clients (https://themightymarketer.com/marketing-survey/) — generic, not CFD-specific.
- **Instrument-vendor referral** — Anton Paar runs technical centres, ASPIRE (academic instrument placement) and the Rheo Academy; **no formal consultant-referral programme was found** (https://www.anton-paar.com/us-en/services-support/anton-paar-technical-centers/ ; https://www.anton-paar.com/us-en/anton-paar-aspire/). Cunningham's path (vendor sales engineer → independent) is the informal version of that channel.
- **Subcontracting to larger consultancies** — only generic consulting-practice evidence found (https://victorvalentineromo.com/articles/subcontracting-consulting-practice/); no CFD-specific case.
- **Open-source reputation** — Wikki, CFD Direct, Applied CCM and CFD Support all trade on it; the revenue it produces directly is thin (Q3).
- **Free presentation as lead magnet** — Rheology Testing Services offers a "free, no obligation presentation on Rheology Principles" (https://www.rheologytestingservices.com/rheology-testing-overview).

### Q6. Failure modes

- **Feast/famine is the stated reason Flow Science left contract work** in 1984 (https://www.flow3d.com/a-life-in-cfd-by-dr-c-w-tony-hirt/). MBO Partners 2025: 76% of independents satisfied, but nearly half cite stress and inconsistent workload (https://www.mbopartners.com/state-of-independence).
- **Year-one dead zone**: CFD Support's 18 months to first paid project with no revenue (Q1). Physics Today's profiles list variable early income, bookkeeping, rate-setting, health insurance and 80-hour weeks (https://physicstoday.aip.org/news/consulting-can-be-stimulating-and-lucrative-for-physicists).
- **"Consulting doesn't scale"** — CFD Support's stated reason for the product pivot; clients "struggled implementing raw OpenFOAM without clear workflows" (Q1).
- **Clients wanting Ansys, liability**: no consultant testimony found. What exists is vendor/comparison prose that commercial codes bring "vendor training and enterprise support [that] reduce project risk" and that the buying question is which workflow the team can "verify, validate, automate, support, and defend" (https://www.resolvedanalytics.com/theflux/comparing-cfd-software ; https://www.bohrium.com/en/blog/tutorials/openfoam-vs-fluent/). The Fluid Codes model (Ansys reseller + services) shows where the money in "clients want Ansys" actually sits.
- **Offshore undercutting**: Indian CFD shops publish $300/$500 minimums and mesh-tiered pricing (https://cfdconsultantsindia.com/cfd-consulting-rates/); Indian senior engineers at $60–90/h (https://qubit-labs.com/average-hourly-rates-offshore-development-services-software-development-costs-guide/). Upwork profile rates were blocked. No case of a US/UK CFD shop closing because of offshore competition was found.
- **Selling hours as a PhD**: the only quantitative anchor is that a PhD lifts the median consulting rate by $20/h to $200/h (IEEE-USA via Physics Today) — a credential premium, not a business model.
- **Regularly documented failure of the free-pilot strategy**: 400+ unpaid hours, no conversion (CFD Support).

### Q7. AI-era angle

- **Physicist-supervised agent development is documented** but as engineering, not as a business: CLAX-PT built with Claude Code over 12 work days / 57 sessions (arXiv 2605.30353, May 2026, https://arxiv.org/abs/2605.30353); a solo developer's fluid-thermal simulation platform on Claude Code (https://pub.towardsai.net/how-i-built-a-production-grade-fluid-thermal-simulation-platform-using-claude-code-0f9176eabe72 — page blocked, only the search abstract read); a simulator plus validation agent in an afternoon (https://medium.com/@vivek-karmarkar/i-built-a-physics-simulator-with-claude-then-built-an-ai-agent-to-prove-it-was-correct-30c2dfcf3a58); Anthropic's own "vibe physics" and long-running scientific-computing write-ups (https://www.anthropic.com/research/vibe-physics ; https://www.anthropic.com/research/long-running-Claude).
- **Funded "agentic simulation" is a venture play, not a solo one**: Quanscient, founded 2021, 40–50 staff, €10M Series A in May 2026 on ~€9M prior (https://www.hpcwire.com/off-the-wire/quanscient-raises-e10m-series-a-to-advance-ai-quantum-and-multiphysics-simulation-platform/ ; https://quanscient.com/blog/agentic-ai-for-multiphysics-simulation-is-here).
- **Generic solo-founder AI SaaS claims** ($2M ARR in 8 months etc.) exist in marketing blogs (https://www.justdoers.com/blog/how-a-solo-founder-used-ai-agents-to-build-a-2m-arr-saas-in-8-months) and are unverifiable; none is scientific software.
- **Not found:** any solo founder with disclosed revenue from an agent-operated simulation, rheology or scientific-software business in 2025–26. The honest reading is that the leverage is real on the build side (days instead of months) and undemonstrated on the sell side.

---

## Case-study table

| Company | Founded | First product/service | Years to viability / to product | Current state | Source |
|---|---|---|---|---|---|
| Flow Science | 1980 | Contract CFD work | 4 yrs to licence sales (1984) | Ongoing, Santa Fe, global | https://www.flow3d.com/a-life-in-cfd-by-dr-c-w-tony-hirt/ |
| Convergent Science | 1997 | CFD consulting | 11 yrs to first licence (2008) | ~80 staff, market leader in IC engines | https://en.wikipedia.org/wiki/Convergent_Science |
| Nabla Ltd | 2000 | Commercial FOAM, Fluent-funded | — | Wound up 2006 after open-sourcing | https://wikki.co.uk/index.php/about-us/ |
| RheoSense | 2001 | VROC chip viscometer | not found | Ongoing; $780k total disclosed funding | https://www.cbinsights.com/company/rheosense |
| Wikki Ltd | 2004 | OpenFOAM consulting/dev | not found | Ongoing, UK/DE | https://wikki.co.uk/ |
| Veryst Engineering | 2005 | Engineering consulting (design, failure) | not found | Ongoing, Boston | https://www.veryst.com/about-us |
| Gammadot Rheology | 2006 | Rheology testing + failure analysis | founder had 18 yrs at Rapra | Ongoing, UK | https://www.gammadot.com/About_Us/gammadot.htm |
| CoolSim / Applied Math Modeling | 2008 | Niche SaaS on Fluent | product day 1 (built in-house at Fluent 2005) | Ongoing | https://www.coolsimsoftware.com/about/ |
| CFD Support | 2009 | OpenFOAM consulting | 1.5 yrs to first paid job; product TCFD after ~5 yrs | "Stable", agents on 3 continents | https://www.linkedin.com/pulse/cfd-support-another-dream-come-true-lubos-pirkl?articleId=6724652159088001024 |
| Engys | 2009 | OpenFOAM consulting + HELYX | not found | 5-country offices | https://engys.com/ |
| Applied CCM | 2011 / 2015 | OpenFOAM support, Caelus | not found | Ongoing, AU/CA/US | https://www.appliedccm.com/about/ |
| Centre for Industrial Rheology | 2012 (consulting since 1998) | Testing, consulting, training | consulting from 1998; lab 14 yrs later | Ongoing, UK | https://www.rheologylab.com/about/ |
| CFD Direct | 2015 | OpenFOAM dev, training, cloud, support | funded maintenance from 2017–18 | 10 yrs; 1,350 cloud subs, 81 training days/yr | https://cfd.direct/about/year-10-2025/ |
| CFD Engine (solo) | ~2012 (Co. no. 08271608) | Fixed-price aero CFD programmes | not found | Solo, UK | https://cfdengine.com/ |
| Rheology Testing Services (solo) | not found | Contract rheology assays | founder ex-GSK 24 yrs | Solo, Chapel Hill NC | https://www.rheologytestingservices.com/rheology-testing-overview |
| Rheology Solutions | 1998 | Instrument sales + testing + training | not found | Ongoing (Online Rheometer Pty) | https://www.labonline.com.au/directories/rheology-solutions-pty-ltd-422767475 |
| Formulaction | 1994 | Turbiscan instrument | not found | 36 staff; acquired by Verder 2023 | https://www.oaklins.com/nl/en/deals/109018/ |
| Cascade Technologies | not found | CharLES LES code, SBIR | not found | Acquired by Cadence Nov 2022 | https://community.cadence.com/cadence_blogs_8/b/cfd/posts/cadence-welcomes-cascade-technologies-1221740980 |
| Fluid Codes | 1997 | Ansys reselling + services | not found | 4 offices, MEA | https://fluidcodes.com/about/ |
| Quanscient | 2021 | Cloud multiphysics SaaS (VC) | — | 40–50 staff, €10M Series A May 2026 | https://tech.eu/2026/05/26/quanscient-lands-eur10m-to-advance-ai-and-quantum-native-hardware-engineering/ |

---

## Funding routes table

| Route | Amount | Duration | Next deadline | Eligibility notes | Success rate | Source |
|---|---|---|---|---|---|---|
| NSF SBIR/STTR Phase I (NSF 26-510) | up to $305,000 | 6–18 mo | Pitch rolling; full proposals 4 Nov 2026, 4 Mar 2027 | ≤500 staff; PI >50% employed at award; Pitch invitation required; STTR needs 60/40 with a research institution | ~15–20% of full proposals; ~9–12% pitch-to-award | https://grantedai.com/news/nsf-sbir-sttr-250m-2026-05-31 ; https://www.getcada.com/insights/nsf-sbir-success-rates |
| NSF SBIR Phase II | up to $1,250,000 | 24 mo | after Phase I | — | ~60% (all agencies) | https://www.ebhoward.com/success-rates-of-the-national-science-foundations-sbir-sttr/ |
| NSF 26-511 instrumentation pilot | $40M pool | — | see solicitation | "next-gen scientific instrumentation" incl. microfluidic platforms | new | https://grantedai.com/news/nsf-sbir-sttr-250m-2026-05-31 |
| NIH SBIR Phase I (PA-27-100, R43) | $323,090 / $400k / $700k by institute | 6–12 mo typ. | 5 Sept 2026 first receipt, then standard cycles | same SBA rules | not found for FY26 | https://grantedai.com/news/pa-27-100-1-4b-2026-05-31 |
| DOE SBIR Phase I (Genesis Mission) | $10M pool, ~40 awards | — | Pitch 10 Sept 2026 | pitch mandatory | — | https://minnesotasbir.org/federal-funding/doe-opens-fy2026-genesis-mission-sbir-sttr-phase-i-opportunity-under-new-pitch-first-process/ |
| USDA NIFA SBIR Phase I | up to $175,000 (+$6,500 TABA) | 8 mo | FY27 not yet announced (FY26 was 15 July 2026) | — | not found | https://grantedai.com/grants/usda-small-business-innovation-research-sbir-phase-i-u-s-department-of-agriculture-usda-5918acdc |
| NSF I-Corps Teams | $50,000 | ~7 weeks | rolling | needs NSF lineage or regional-hub letter; EL is typically PhD student/postdoc | — | https://www.nsf.gov/funding/opportunities/nsf-national-innovation-corps-teams-nsf-national-i-corps-tm |
| Great Lakes I-Corps regional (IN/IL) | training, no cash | ~4 weeks | cohorts | UIUC/Purdue co-lead | — | https://greatlakesicorps.org/programs/national-nsf-i-corps/ |
| Indiana SBIR/STTR match (ARI) | $0.50/federal $ up to $50k per award | — | rolling | IN company, Phase I or II | — | https://www.insideindianabusiness.com/articles/elevate-ventures-matching-funds-for-small-tech-startups |
| Elevate Ventures pre-seed | $20k–100k | — | — | IN, no match required | — | https://www.insideindianabusiness.com/articles/state-boosting-funding-to-support-entrepreneurship-innovation |
| Illinois DCEO SBIR match | up to $75k (Ph I), $250k (Ph II) | — | rolling to 30 June 2026 (FY) | IL HQ; Phase I within 12 mo; ≥51% work in IL | ~100 awards/yr | https://dceo.illinois.gov/whyillinois/sbir-sttr-phase1.html |

---

## Revenue path most consistent with the evidence, for this founder

The evidence is thin, biased toward survivors, and mostly non-quantitative; the path below is what it supports, not what it proves.

**Months 0–12 (still 5–10 h/week, PhD finishing)**

1. **Sell fixed-price studies, not hours, from the first client.** The one solo CFD shop with a public model (CFD Engine) sells fixed-price programmes; 8020 Engineering's £1,990/£3,990/£5,990 tiers and Cad Crowd's $2k–20k range bracket the price. A rheology-plus-simulation study for a formulation, device or slurry client at $4–10k is inside the published band. At 5–10 h/week the capacity is 2–4 such studies a year, i.e. $10–40k — proof of demand, not income.
2. **Start the teaching channel now, because it is the only channel with a founder calling it "strongest".** A monthly webinar or a written series on the founder's real niche (dense-suspension rheology, injectability, slurry flow) is what CoolSim, CFD Direct and Cunningham all did; CFD Support's cold-calling year is the documented control.
3. **Do not run free pilots.** 400+ hours, zero conversions (CFD Support).
4. **Incorporate early, but keep the SBIR clock in view.** A pitch to NSF 26-510 can go in before the PI is >50% employed; the award requires it. The 4 Nov 2026 and 4 Mar 2027 full-proposal dates mean a Project Pitch in autumn 2026 lands an award decision around graduation, which is the one moment the >50% rule is easy to meet. Precedent in-niche: QATCH's viscosity-for-injectability Phase I. Expect ~10% odds end-to-end; treat it as an option, not a plan.
5. **Register with the Great Lakes I-Corps hub** as a PhD-student entrepreneurial lead (no NSF lineage needed at the regional stage) — the customer-discovery interviews double as the outreach the cold-call approach failed at.
6. **Settle IP with the TTO in writing** before the thesis is deposited: which code is personal, which is university, what a non-exclusive licence costs. Default terms elsewhere are ≤5% equity plus low-single-digit royalties, which is survivable; ambiguity is not.

**Months 12–24 (full-time)**

7. **Target 15–30 fixed-price studies a year at $5–10k**, which is $150–300k gross on published rates; at 50% utilisation and $150–200/h the same number falls out of hourly arithmetic. This is the Flow Science / CFD Support "viable consultancy" stage and every documented company passed through it.
8. **Use a state match if a Phase I lands**: +$50k (Indiana) or +$75k (Illinois) — the location decision is worth $25k of grant money and the Illinois match also requires ≥51% of the work in-state.
9. **Watch for the bottleneck that recurs across clients** — that is where Convergent (meshing) and CFD Support (workflow) found their product. For this founder the candidates already exist in the biorheo repo (record→solver export, profile intervals, needle-force study). Do not build it before three clients have paid for the manual version.
10. **Keep the training line as ~20–30% of revenue** the way CFD Direct does (81 days/yr), because it is counter-cyclical to project famine and it is the lead source.
11. **Expect year 2 to be the first year that looks like a business.** No documented case in the sample reached stable revenue in under 18 months; Flow Science and CFD Support both said contract work alone did not hold, and both fixed it with a product only after years 4–5.

What the AI leverage buys, on the evidence: build time (days per module in the arXiv case study), which lets a one-person shop deliver the study *and* the reusable tool from the same engagement. It does not, on any evidence found, shorten the sales cycle.

---

## Evidence AGAINST

- **Survivorship**: every case study is a company that still exists. Dead solo consultancies leave no trace; the failure rate is unknown.
- **The one solo CFD shop found (CFD Engine) publishes no revenue**, and its founder also sells courses and runs a podcast, which may be the income rather than a channel.
- **Time to product was 4–11 years** everywhere it is documented, and those were 2–5 person teams working full-time, not a 5–10 h/week founder.
- **$150–300k from studies needs 15–30 paid engagements a year**, and CFD Support needed 18 months to land the first one, with a full-time two-person team.
- **Offshore pricing exists at $300–500 per study minimum** (India). The founder's defence is niche depth (rheology + simulation + validation), which none of the offshore pages claim — but no evidence shows that defence holding for a US solo shop.
- **Rheology contract testing has a visible price ceiling** ($170–340/sample; $50–64/h instrument time at universities) that a solo shop with one instrument cannot scale past.
- **Open-core / maintenance money is ~0.1% of commercial CFD revenue** and the Foundation itself calls it inadequate.
- **SBIR odds are ~10% end-to-end** and the PI employment rule is incompatible with keeping a day job; a Phase I is a lottery ticket that also forces a full-time decision.
- **No AI-era solo scientific-software business with disclosed revenue was found**; the funded "agentic simulation" company has 40–50 staff and €19M raised.
- **"Clients want Ansys"**: no consultant testimony found either way, so the risk is unmeasured, not absent.

---

## Open questions

1. What do the Indian CFD shops actually charge per tier (image-only price chart), and what do PhD-level OpenFOAM freelancers list on Upwork? (Both pages returned 403.)
2. What did the CFD-online "CFD consulting work pricing" thread say — it is the only primary-source peer discussion found and it was blocked.
3. Does any instrument vendor (Anton Paar, TA, Malvern) run a formal consultant/partner referral list? Nothing public found; a direct ask to a US application lab would settle it.
4. Which four technology areas does the DOE FY26 Genesis Mission SBIR cover, and does slurry/electrode rheology fit? Pitch deadline 10 Sept 2026.
5. What are Purdue / Illinois / Notre Dame's actual spin-out equity and royalty terms for student-founded companies?
6. What does the Talking CFD archive (cfdengine.com/newsletter) say episode by episode about client acquisition? Only two episodes were readable.
7. Is there any disclosed P&L for a 1–3 person rheology testing lab? None found; the UK trio (Cunningham, Edwards, Patrick) are the natural people to ask.
8. Does the NIH parent SBIR list a specific injectability/formulation-modelling topic for FY27, or is it investigator-initiated only?

---

## Sources

- https://cfd.direct/about/year-10-2025/ (fetched; year-10 numbers)
- https://cfd.direct/about/
- https://cfd.direct/openfoam/management/
- https://wikki.co.uk/index.php/about-us/
- https://old.engys.com/about-us ; https://engys.com/
- https://en.wikipedia.org/wiki/Convergent_Science
- https://news.wisc.edu/engine-software-from-uw-spinoff-being-used-around-the-world/ (fetched)
- https://www.flow3d.com/a-life-in-cfd-by-dr-c-w-tony-hirt/ (fetched)
- https://www.flow3d.com/flow-science/
- https://www.linkedin.com/pulse/cfd-support-another-dream-come-true-lubos-pirkl?articleId=6724652159088001024 (fetched)
- https://www.cfdsupport.com/about-us/
- https://www.appliedccm.com/about/ ; https://cfdengine.com/podcast/darrin-stephens-chris-sideroff-caelus/ (fetched)
- https://cfdengine.com/ (fetched) ; https://cfdengine.com/podcast/ ; https://uk.linkedin.com/in/robinknowles ; https://learn.cfdengine.com/
- https://cfdengine.com/podcast/paul-bemis-coolsim/ (fetched) ; https://www.coolsimsoftware.com/about/
- https://community.cadence.com/cadence_blogs_8/b/cfd/posts/cadence-welcomes-cascade-technologies-1221740980 ; https://www.sbir.gov/portfolio/118337
- https://www.veryst.com/about-us ; https://www.veryst.com/about/our-team/stuart-brown
- https://fluidcodes.com/about/
- https://www.labonline.com.au/directories/rheology-solutions-pty-ltd-422767475
- https://www.rheologylab.com/about/ ; https://sites.rutgers.edu/centerfordermalresearch/wp-content/uploads/sites/271/2022/07/Neil-Cunning-bio_photo_abstract.pdf
- https://www.gammadot.com/About_Us/gammadot.htm
- https://www.rheologytestingservices.com/rheology-testing-overview (fetched)
- https://www.rheosense.com/about-rheosense ; https://www.cbinsights.com/company/rheosense
- https://www.verder-scientific.com/news/news/details/news_detail/verder-group-acquires-formulaction ; https://www.oaklins.com/nl/en/deals/109018/
- https://nrf.aux.eng.ufl.edu/resources/resource.asp?id=105
- https://www.cheme.washington.edu/research/facilities/bindra/rheometer
- https://www.wcet.washington.edu/per-sample-analysis/
- https://www.cadcrowd.com/blog/what-are-cfd-engineering-rates-consulting-services-costs-company-pricing/ (fetched)
- https://www.resolvedanalytics.com/ebook/cfd-costs (fetched)
- https://www.cfdsource.com/ShowBlog/whats-the-real-cost-of-cfd-consulting/
- https://8020engineering.com/consulting-pricing-and-packages/ (fetched)
- https://cfdconsultantsindia.com/cfd-consulting-rates/ (fetched)
- https://qubit-labs.com/average-hourly-rates-offshore-development-services-software-development-costs-guide/
- https://physicstoday.aip.org/news/consulting-can-be-stimulating-and-lucrative-for-physicists (fetched)
- https://www.kolabtree.com/blog/scientific-consultants-hire-freelance-guide/ (fetched)
- https://www.payscale.com/research/US/Job=Independent_Consultant/Hourly_Rate
- https://www.quora.com/What-is-the-consulting-daily-rate-for-a-PhD-holder-with-more-than-20-years-experience
- https://www.cfd-online.com/Forums/main/91105-cfd-consulting-work-pricing.html (403)
- https://www.upwork.com/hire/openfoam-freelancers/ ; https://www.upwork.com/hire/openfoam-freelancers/in/ (403)
- https://tracxn.com/d/companies/openfoam/__8Dfq1x_cN9UYJ8i86QHjdpm19Ts7i9Q9psm4lTlUg84
- https://openfoam.org/
- https://jonathanstark.com/hbin
- https://mtlynch.io/solo-developer-year-6/
- https://grantedai.com/news/nsf-sbir-sttr-250m-2026-05-31 (fetched)
- https://www.bwcoconsulting.com/fod/nsfsbir
- https://seedfund.nsf.gov/solicitation-eligibility/ ; https://sbir.org/getting-started/pi-key-staff-eligibility-sbir-sttr/
- https://www.ebhoward.com/success-rates-of-the-national-science-foundations-sbir-sttr/ (fetched)
- https://www.getcada.com/insights/nsf-sbir-success-rates ; https://fundunlocked.com/learn/nsf-sbir-success-rate/
- https://grantedai.com/news/pa-27-100-1-4b-2026-05-31 ; https://www.ebhoward.com/nih-sbir-sttr-budget-limits-vary-more-than-most-applicants-realize/
- https://www.sbir.gov/portfolio/1168615 (QATCH)
- https://minnesotasbir.org/federal-funding/doe-opens-fy2026-genesis-mission-sbir-sttr-phase-i-opportunity-under-new-pitch-first-process/
- https://grantedai.com/grants/usda-small-business-innovation-research-sbir-phase-i-u-s-department-of-agriculture-usda-5918acdc
- https://www.nifa.usda.gov/application-deadline-small-business-innovation-research-program-phase-i-sbir-funding-opportunity
- https://www.nsf.gov/funding/opportunities/nsf-national-innovation-corps-teams-nsf-national-i-corps-tm
- https://greatlakesicorps.org/programs/national-nsf-i-corps/ ; https://icorps.illinois.edu/
- https://www.insideindianabusiness.com/articles/elevate-ventures-matching-funds-for-small-tech-startups
- https://www.insideindianabusiness.com/articles/state-boosting-funding-to-support-entrepreneurship-innovation
- https://iedc.in.gov/events/news/details/2023/12/19/indiana-commits-increased-funding-resources-to-catalyze-entrepreneurship-innovation
- https://dceo.illinois.gov/whyillinois/sbir-sttr-phase1.html ; https://grantedai.com/news/illinois-5-million-sbir-sttr-state-matching-program-2026
- https://www.ropesgray.com/en/insights/alerts/2025/02/considerations-when-contracting-with-academic-spinouts
- https://medium.com/included-vc/startups-science-ip-intellectual-property-licensing-101-for-deep-tech-founders-0685a13851fd
- https://www.anton-paar.com/us-en/services-support/anton-paar-technical-centers/ ; https://www.anton-paar.com/us-en/anton-paar-aspire/
- https://themightymarketer.com/marketing-survey/
- https://victorvalentineromo.com/articles/subcontracting-consulting-practice/
- https://www.resolvedanalytics.com/theflux/comparing-cfd-software ; https://www.bohrium.com/en/blog/tutorials/openfoam-vs-fluent/
- https://www.mbopartners.com/state-of-independence
- https://arxiv.org/abs/2605.30353
- https://pub.towardsai.net/how-i-built-a-production-grade-fluid-thermal-simulation-platform-using-claude-code-0f9176eabe72 (403)
- https://medium.com/@vivek-karmarkar/i-built-a-physics-simulator-with-claude-then-built-an-ai-agent-to-prove-it-was-correct-30c2dfcf3a58
- https://www.anthropic.com/research/vibe-physics ; https://www.anthropic.com/research/long-running-Claude
- https://www.hpcwire.com/off-the-wire/quanscient-raises-e10m-series-a-to-advance-ai-quantum-and-multiphysics-simulation-platform/
- https://tech.eu/2026/05/26/quanscient-lands-eur10m-to-advance-ai-and-quantum-native-hardware-engineering/
- https://quanscient.com/blog/agentic-ai-for-multiphysics-simulation-is-here
- https://www.justdoers.com/blog/how-a-solo-founder-used-ai-agents-to-build-a-2m-arr-saas-in-8-months (unverified marketing claim)
