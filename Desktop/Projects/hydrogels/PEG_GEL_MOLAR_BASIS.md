# PEG gels: stating the recipe on a molar basis

## Short version

The spreadsheet already works in moles. The written protocol describes it in the wrong
units, which is why the recipe looks like a mass or volume recipe when it isn't.

Line 7 of the sheet, "80:20", is a **mole ratio** of PEGMEA to PEGDA. The numbers
0.7328 and 0.2672 are **mass fractions** (dimensionless) that the sheet computes from
that mole ratio, and they are not "% gram/mL".

Check, with PEGMEA Mn = 480 g/mol and PEGDA Mn = 700 g/mol:

    w_PEGMEA = (0.80 x 480) / (0.80 x 480 + 0.20 x 700) = 384 / 524 = 0.7328
    w_PEGDA  = (0.20 x 700) / 524                       = 140 / 524 = 0.2672

Those reproduce the sheet to four decimals, so the conversion is settled. Going the
other way from the worked example confirms it: 21.984 g / 480 = 0.04580 mol PEGMEA,
8.016 g / 700 = 0.011451 mol PEGDA, a ratio of exactly 4.000 : 1.

Two consequences for the protocol text:

- The basis is 30 **grams**, not 30 mL. Multiplying a mass fraction by a volume is a
  units error that happens to be invisible here only because the numbers sum to 1. At a
  mixture density near 1.11 g/mL, 30 g of this blend is about 27 mL. The DMPA step
  already says "30 g x 0.5%", so the rest of the protocol should say grams too.
- Both Mn values are load-bearing. The recipe is only reproducible if the bottle's Mn
  is written down beside it.

## Why moles and not mass

A gel's stiffness and mesh size are set by network topology: how many crosslink
junctions there are per unit volume, and how long the strands between them are. Both of
those are counts of molecules. Chemistry counts molecules; a balance counts grams. The
two are proportional only within a single molecular weight.

The gap is large in practice. One gram of PEGDA at Mn 700 is 1.43 mmol of crosslinker.
One gram of PEGDA at Mn 6000 is 0.167 mmol. Same mass on the balance, roughly nine times
fewer junctions in the gel. "20% PEGDA by mass" therefore does not specify a gel; it
specifies a gel only once you also name the Mn. A recipe passed to another lab, or into
a paper, has to carry the molar ratio or it cannot be reproduced.

Mass tells you how much stuff went in. Moles tell you how many reactive ends went in, and it is the reactive ends that build the network.

## What each ingredient does

**PEGDA** (PEG diacrylate) is the **crosslinker**. It carries an acrylate group at *both*
ends of the PEG chain, so once each end has been incorporated into a growing chain, that
molecule is a bridge tying two chains together. Those bridges are the junctions that turn
a solution into a solid. Raising the PEGDA mole fraction adds junctions, shortens the
average strand between them, and stiffens the gel.

**PEGMEA** (PEG methyl ether acrylate) is **monofunctional**, not a crosslinker. One end
is a capped methyl ether, which is unreactive; only the acrylate end can react. It can be
built into a chain but cannot bridge to a second one, so it terminates as a dangling
side arm. It is what dilutes the junctions and softens the gel. Increasing PEGMEA is how
you move down in modulus without changing the chemistry.

**DMPA** (2,2-dimethoxy-2-phenylacetophenone, sold as Irgacure 651) is the
**photoinitiator**. It is not part of the network in any load-bearing sense; it is the
starter. It absorbs near 365 nm and undergoes Norrish Type I alpha-cleavage, splitting
into a benzoyl radical and a dimethoxybenzyl radical. Each radical adds across an
acrylate C=C, and that opened double bond becomes the next radical, which adds to the
next acrylate, and so on. This is ordinary free-radical chain-growth polymerisation; the
only job of the UV lamp and the DMPA is to decide *when* it starts.

So DMPA sets the kinetics, not the topology. More DMPA means more chains started at once,
so a faster cure but shorter kinetic chains; less DMPA means a slower cure, and if it is
too low, oxygen inhibition at the free surface can leave a tacky uncured skin. It does not
appear in the crosslink count, which is why 0.5 wt% is a reasonable way to dose it.

## Rewriting the protocol

Keep the spreadsheet; it is doing the right arithmetic. Change the labels and add the
molecular weights, so the molar step is visible instead of buried.

Specify a batch by three numbers:

1. **Total mass** m_total (g), not volume.
2. **Mole fraction of PEGDA**, x_DA = n_DA / (n_DA + n_MEA). The example above is
   x_DA = 0.20.
3. **DMPA loading**, 0.5 wt% of m_total.

Then

    w_DA  = x_DA*Mn_DA / (x_DA*Mn_DA + (1 - x_DA)*Mn_MEA)
    m_DA  = w_DA * m_total
    m_MEA = m_total - m_DA
    m_DMPA = 0.005 * m_total

For x_DA = 0.20 and a 30 g batch this returns 8.016 g PEGDA, 21.984 g PEGMEA, 0.15 g
DMPA, which is the existing recipe. Nothing at the bench changes, but the number you quote, x_DA = 0.20, now means the same thing to anyone who reads it.

## Optional: the functional-group basis

If you ever want to compare against the acrylate literature, the more informative
variable is the share of acrylate groups contributed by the crosslinker, because PEGDA
brings two reactive ends per molecule and PEGMEA brings one:

    acrylate equivalents = 2*n_DA + 1*n_MEA
    f_crosslink = 2*n_DA / (2*n_DA + n_MEA)

The example recipe is x_DA = 0.20 on a molar basis but f_crosslink = 0.333 on an
acrylate basis. Both are correct; they answer different questions. Quote x_DA for "what
did I weigh out", f_crosslink for "how crosslinked is the network". Whichever you use,
say which one in the caption, because a bare "20% PEGDA" is ambiguous across all three
bases and the three differ by a lot: 26.7% by mass, 20% by mole, 33.3% by acrylate.

## Things to nail down before this goes in a methods section

- The actual Mn on each bottle, from the certificate of analysis, not the catalogue name.
  Mn 480 and Mn 700 are inferred here by matching the spreadsheet, and they fit exactly,
  but they should be read off the label to confirm.
- Whether the sheet was ever used with the "x 30 mL" reading literally, i.e. whether any
  past batch was made by measuring 30 mL of total volume instead of 30 g. Those batches
  would be off in total amount, though the PEGMEA:PEGDA ratio would still be right.
- Measured densities, if the volume of the finished gel matters for mould filling.
