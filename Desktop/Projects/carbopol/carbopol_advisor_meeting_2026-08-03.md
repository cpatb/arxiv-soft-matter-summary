# Carbopol bath experiment — advisor meeting, 3 Aug 2026

## The one sentence

*"Thrust 3 of the proposal names Carbopol substrates and there's no bench data behind
it. I can get you that data in two weeks on the rig we already have."*

Lead with her proposal. The physics comes second, your own curiosity third. An advisor
says yes fastest to work that reduces her own risk.

---

## 1. Why it's worth doing, in the order that works

**It de-risks the proposal.** Thrust 3 proposes Carbopol at 0.1–1 % w/w, τ_y ≈ 3–125 Pa.
A panel reviewing a substrate-mechanics proposal will want to know the group can make
and control that substrate. One figure answers it: crater depth against impact speed
across a τ_y ladder, with a plain-water bath as the control. Preliminary data is cheap
to produce now and expensive to be missing in January.

**The gap is a clean empty cell.** Draw this on her whiteboard:

|                       | rigid impactor                        | deformable drop        |
|-----------------------|---------------------------------------|------------------------|
| **Newtonian pool**    | sphere entry, classical               | drop-into-pool, vast   |
| **yield-stress bath** | sphere penetration in Carbopol: characterized (de Bruyn, Coussot) | **empty** |

Both neighbours are well studied. Yield-stress *drops* onto rigid solids are done too
(Luu & Forterre). The crossed cell, a Newtonian drop into a yield-stress bath, is where
deformability meets a substrate that keeps a permanent record. That one is thin.

**The yield stress hands you a measurement water can't give you.** The crater freezes.
Every drop-into-pool experiment ever run loses its cavity to gravity inside ~100 ms.
Here it stays put and you photograph or section it at leisure. That turns a hard
high-speed measurement into an easy static one, and it creates an observable with no
Newtonian analogue:

> **M = h_frozen / h_cavity,max**, how faithfully the frozen crater remembers the
> deepest transient cavity.

For a rigid sphere M ≈ 1 by construction: the sphere digs, the crater records the dig.
A drop spreads and splashes, so energy that would have gone into depth goes sideways.
M should fall as We rises, and crater depth against We should turn over rather than
climb. Neither behaviour is available to a rigid impactor. Say this before the data
exists, so it's on record as a prediction.

**It's cheap and already in the building.** Carbopol 940/980, NaOH, the existing backlit
high-speed rig at 36 kfps, the lab rheometer. No drop tower, no new capability,
consumables under a hundred dollars.

**It carries broader impacts.** The protocol was written to be runnable by a high-school
student: one batch, one ladder, a ruler and a camera. That is mentoring language for the
same proposal, out of the same work.

---

## 2. The plan you're proposing (say the boundary out loud)

**Week 1, is the effect there at all.** One Carbopol batch at mid τ_y, characterized on
the rheometer. Water drops, 5 heights, 3 reps. Side-view high-speed plus frozen-crater
photos. Plot h_frozen and h_cavity,max against We, compute M.

**Kill criterion, stated up front:** if M doesn't move with We, the deformability story
is dead and I stop at the end of week 2.

**Weeks 2–3, only if week 1 says go.** The ladder: 4 Carbopol concentrations × 4 speeds,
plus a water bath at the same 4 speeds as control, 20 runs, all read off the same
backlit frames.

Volunteering the stopping rule is your strongest scope-control move. It reads as rigour,
and it puts an end date in the room before she can supply one.

---

## 3. What to actually ask her for

1. **~4 h of rheometer time**, or permission to bank aliquots and batch the flow curves.
2. **Two weeks of protected bench time**, in the gaps between imaging sessions.
3. **One design decision, and this is the one to hand her.**

> *"Do I index the ladder on wt %, or on measured yield stress?"*

You already know the answer is measured τ_y: Carbopol's yield stress depends on
neutralization and mixing history, so wt % is a recipe label, not a material property.
Ask anyway. She will say τ_y, the rheometer time then follows from her own answer
instead of from your request, and the design becomes something she chose. You aren't
staging it. It's a real question, and letting her make that call is what makes the rest
of the meeting easy.

**One design point to volunteer, unprompted:** the bath needs ≥ 110 mm inner diameter
filled ≥ 55 mm deep, because the yielded region around the impact scales as τ_y^(−1/3)
and the weakest gel at the fastest drop reaches ~50 mm. Bring the number. Having already
sized the container is what separates this from an idea.

---

## 4. How to run the room

- **Order:** her proposal, then the empty cell, then the frozen crater, then the
  two-week plan, then the three asks. Roughly ten minutes to the first ask.
- **Bring a physical prop if you can.** One Carbopol crater in a dish, made tonight,
  beats any slide. She can see the substrate keeping the record.
- **Come with a proposal, not a question.** An open "what do you think I should do?"
  invites her to design it, and a designed-by-advisor experiment grows. A bounded plan
  invites yes or no.
- **When she adds scope** (glycerol drops, surfactant, bottom-view, a second Carbopol
  grade) write it down visibly and say *"that's the version after this one works."*
  Recording it honours the idea without accepting it this month.
- **Close by naming the next checkpoint yourself:** *"I'll bring you M vs. We in two
  weeks."* Otherwise she sets the date.

---

## 5. Questions she will ask

**"Hasn't someone done this?"** The 2×2 above. Then be straight: *"My search covers
recent arXiv, so it's a recency check, not a full novelty check. If you know something
older I'd rather hear it now than after the runs."* Inviting the counterexample costs
nothing and protects you from the one thing that can sink the meeting.

**"How does this fit Paper 1?"** It doesn't compete. Different substrate class, two
weeks, in the gaps. Say plainly that Paper 1 is still the thesis, you aren't asking to
trade it, and this is the piece with a grant behind it.

**"What sets the yield stress?"** Neutralization to pH 6–7 with NaOH, degassed overnight,
flow curve per batch, Herschel–Bulkley fit, index on τ_y.

**"How do you know the crater is arrested and not just slow?"** Photograph it again at
5 min and at 1 h. A yield stress means it stops. Anything still creeping is a
shear-thinning fluid and belongs in the Newtonian control.

**"Have you modelled it?"** *"A scaling argument, not a model. That's where the tank size
comes from."* True, and a complete answer.

If she asks directly whether you're running simulations, don't deny it flatly. Repos,
clusters and collaborators all leak, and a denial costs far more later than a deferral
does now. *"Nothing I'd want to show yet"* is true, closes the topic, and costs you
nothing.

---

## 6. Do not say

These numbers only exist because of the simulations. Each one invites *"where's that
from?"*, and there's no sim-free answer:

- crater moves ~80×, λ moves 33 % of its own mean
- the water control rises 91 % from U₀ = 1 → 3 m/s
- fluidized fraction, yielded fraction, J^(−0.88)
- anything about λ collapsing (or not) on J
- "which observable reports on the substrate" as a framing

The deformability framing in §1 needs none of them. Keep the observable-selection
question in your pocket. It becomes the *second* paper, and by then you'll have bench
data that motivates it on its own.
