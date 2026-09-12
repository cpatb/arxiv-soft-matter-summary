# R1 outreach drafts — send-ready

Date: 2026-09-12. Companion to `R1_TARGETS.md`. Four drafts, one per target type.

## Design rules these drafts obey

1. **It is a question, not a pitch.** F9's only documented failure is a year of cold sales
   calls plus a 400-hour free pilot (CFD Support). None of these asks for a meeting about a
   service; each asks one thing the sender genuinely does not know.
2. **No company name.** D-G/D1 is undecided, and signing with an invented company is worse
   than signing as a person. Sign as a physicist. Add the company later.
3. **No price, no deliverable, no attachment** on first contact. `OFFER.md` §5 prices are
   marked INFERRED and must not reach a prospect.
4. **Under 150 words.** Every draft below is counted.
5. **One concrete number from their own work**, so it is unmistakably not a mail-merge.
6. **The ask is 15 minutes or a reply-by-email answer** — the lower-friction option is
   offered explicitly, because a written answer scores the R1 sheet just as well as a call.

⚠ Fill the bracketed slots by hand. Do not send with a bracket left in.

---

## Draft 1 — paper author who published the failure (A1, A2, A3)

> Subject: your [MAbs 2024 / JHU 2023] result on extrapolating to needle shear rates
>
> Dear [Dr Armstrong / Prof. Rattray / Ms Nawathe],
>
> I'm a physicist finishing a PhD on dense suspensions and high-shear rheology, and I'm
> trying to work out whether one specific gap is worth building a small practice around.
> Your [2024 MAbs paper / thesis] is the clearest statement of it I've found — [the −186 %
> to +177 % spread when extrapolating from low-concentration measurements / that all three
> force models hold to <25 % error but break above 230 mg/mL and in 27G and 30G needles].
>
> One question, and an emailed sentence is a complete answer:
>
> **When that disagreement showed up, what happened next — was it modelled again in-house,
> measured differently, or left as a known limitation?**
>
> If you have 15 minutes some time I'd also be glad to hear how you'd want the high-shear
> end handled. Happy to send you what I've built on public blood and microgel data if it's
> of any use either way.
>
> Best regards,
> [Name] · [PhD candidate, institution] · [phone]

*(141 words.)* Scores §7 questions 1 and 3. If the reply is "left as a known limitation",
that is a **direct hit on the venture's premise** and should be logged verbatim.

---

## Draft 2 — device maker or device consultancy (C1, C2, C3)

> Subject: error bars on Hagen–Poiseuille glide-force predictions — one question
>
> Dear [name],
>
> I'm a physicist finishing a PhD on suspension rheology. [Your Molly 2.25 mL study /
> Your subcutaneous-injection case study] models the pressure drop with [Hagen–Poiseuille /
> a COMSOL poroelastic tissue model], which is also what most of the published work does.
>
> My question is narrow: **has a reviewer, a notified body or an internal design-control
> process ever asked you for an uncertainty on a predicted injection force — or for a
> justification of the viscosity model used?**
>
> I ask because I build identified constitutive models with confidence intervals and a
> statement of which parameters the answer does not depend on, and I cannot tell from the
> outside whether anyone is required to care. A one-line answer would settle it.
>
> Best regards,
> [Name] · [PhD candidate, institution]

*(126 words.)* Scores questions 4 and 6. For C2 (Veryst) add one sentence: *"If it is ever
useful to subcontract the non-Newtonian rheology layer of a study, I'd like to be on your
list."* That is the rank-2 test.

---

## Draft 3 — regulatory / CM&S (D1, D2)

> Subject: has viscosity-model choice ever appeared in a deficiency letter?
>
> Dear [name],
>
> I'm a physicist finishing a PhD on suspension rheology, and I'm trying to answer one
> question that I don't think is answerable from the public record.
>
> **In blood-contacting or drug-delivery CFD submissions you've seen: roughly what fraction
> use a non-Newtonian viscosity model — and has the choice of viscosity model ever been
> raised in a deficiency letter or an information request?**
>
> The FDA nozzle benchmark work used a single viscosity at 200 1/s and the reported spread
> between models is about 6× at low Reynolds number, so I can see why it *should* matter. I
> can't tell whether it *does*.
>
> [For Numalogics only: I also wanted to say that the numaScrew MDDT qualification with
> Sawbones as the physical comparator is the clearest precedent I've found that a small
> group can reach regulatory-grade work — if you'd ever share how the comparator
> partnership was structured, I'd be grateful.]
>
> Best regards,
> [Name] · [PhD candidate, institution]

*(148 words without the bracketed paragraph.)* Scores question 6 and settles whether rank 2
has a requirement behind it or only a gap.

---

## Draft 4 — platform / hyperconcentration company (E1)

> Subject: what breaks first above 400 mg/mL?
>
> Dear [name],
>
> I'm a physicist finishing a PhD on dense suspensions — jamming and shear thickening — and
> at 400–500 mg/mL a protein formulation stops being a viscous liquid and starts being a
> suspension, which is a different physics with different failure modes (clogging, local
> jamming at the needle entrance, a viscosity that depends on how fast you push).
>
> **My question: at those concentrations, what fails first in practice, and is any of it
> modelled at all?**
>
> I'm not selling anything today; I'm trying to find out whether the problem I can work on
> is a problem anyone has. If it is, I'd like to hear what the honest version of it is.
>
> Best regards,
> [Name] · [PhD candidate, institution]

*(126 words.)* Scores question 5 and is rank 3's kill test verbatim.

---

## If someone asks "what would this cost?"

Do not answer with a number on the first exchange. Say: *"It depends on whether you have a
flow curve already — let me send you the one-page deliverable and the scope, and a fixed
price with it."* Then correct `OFFER.md` §5 (currently INFERRED, no public comparable exists
for this segment) before sending anything.

## If someone asks "what have you actually done?"

The three honest artefacts, in order of strength:

1. `biorheo/data/foam_koponen2023/foam_pipe_anchor.md` — a model identified on a 19 mm pipe
   predicts an 8.15 mm pipe to **18.7 %**, against a positive control at 2.4 %, on 660
   measured rows with independently imaged wall slip. **This is the only external
   measured-flow anchor in the repo.** Say it is a foam, not a biofluid.
2. `biorheo/case_studies/blood_injectability/` — 20 real donors, five constitutive
   families, profile-likelihood intervals, and the finding that **a donor's own parameter
   uncertainty (median 29 % of the force) is larger than the effect of which model family
   you pick (12.5 %)**. That is the sentence the service exists to say.
3. `biorheo/case_studies/deliverable_example/` — the page a client receives, including the
   versions that **refuse** to answer and say why.

And state the limits without being asked, because they will be found anyway: the blood
records extrapolate **1.42 decades** past their measured envelope to reach a 27G needle, and
the published force numbers carry **zero needle/plunger friction** (`friction: float = 0.0`),
so they are the hydrodynamic term only.
