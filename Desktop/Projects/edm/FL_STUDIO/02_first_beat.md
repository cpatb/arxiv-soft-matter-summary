# Your first house beat, in twenty minutes

No theory required. Follow the steps, and read `THEORY/03_rhythm.md` afterwards to find out
what you just did.

Everything here uses samples that are already on your machine. Nothing to download.

## Setup

1. **File > New**. Remember that Ctrl+N is "save new version" in FL, not "new file".
2. Right-click the tempo readout in the top toolbar, choose **Type in value**, enter **124**.
3. Press **F6** to open the Channel Rack. If it already has channels loaded, right-click each
   name and choose **Delete** to clear it out.
4. Make sure you are in **Pattern mode**. The toolbar shows PAT or SONG, and **L** toggles it.

## Load four sounds

Open the Browser on the left (**Alt+F8** if it is hidden) and navigate to **Packs > Drums**.

Drag each of these onto an empty slot in the Channel Rack:

| Sound | Browser location |
|---|---|
| Kick | `Packs > Drums > Kicks > 909 Kick` |
| Clap | `Packs > Drums > Percussion > 909 Clap` |
| Open hat | `Packs > Drums > Hats > 909 OH` |
| Closed hat | `Packs > Drums > Hats > 909 CH 1` |

These are samples of the Roland TR-909, the drum machine that Chicago house and techno were
built on. Single-clicking a sample in the browser previews it, so you can hear them before
you commit.

## Program the pattern

Each channel has sixteen step buttons. They are one bar, divided into sixteenth notes.
Left-click to turn a step on, right-click to turn it off.

Turn on these steps:

```
Step        1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
Count       1   e   &   a   2   e   &   a   3   e   &   a   4   e   &   a

909 Kick    X   .   .   .   X   .   .   .   X   .   .   .   X   .   .   .
909 Clap    .   .   .   .   X   .   .   .   .   .   .   .   X   .   .   .
909 OH      .   .   X   .   .   .   X   .   .   .   X   .   .   .   X   .
909 CH 1    X   X   .   X   X   X   .   X   X   X   .   X   X   X   .   X
```

Press **Space**.

That is a house beat. Kick on every beat, clap on 2 and 4, open hi-hat on every offbeat,
closed hats filling the rest.

Build it one channel at a time and listen after each. In particular, program the kick alone,
then add the open hat, and notice how much the open hat does. It is the element that makes
the loop move rather than just plod.

## Make it feel human

Two adjustments, both of which matter more than they sound like they should.

**Swing.** Find the swing knob on the Channel Rack toolbar and turn it up to somewhere
between 15% and 25%. This delays every second sixteenth slightly. The pattern stops sounding
like a machine and starts sounding like a groove. Turn it back to zero and compare.

**Velocity.** Select the closed hat channel, then click the **graph editor** button at the
top-left of the Channel Rack (the small bar-chart icon). You now see a bar for each step.
Drag them so the hats alternate loud and quiet, or so that every fourth one is loudest.

Real drummers do not hit every note identically. Programmed hats at uniform velocity are the
clearest giveaway of an amateur track, and this takes thirty seconds to fix.

## Route to the mixer

Right now every channel goes straight to the master, which means you cannot process them
separately.

Select all the channels in the Channel Rack (Ctrl+A works, or click the first channel button
and Ctrl+click the others), then press **Ctrl+L**. Each channel is assigned to its own free
mixer track.

Press **F9** to open the Mixer and confirm you now have four labelled tracks.

This is the most useful shortcut in FL Studio and it is very hard to discover on your own.

## Balance it

In the mixer, set rough levels. The kick should be the loudest thing. A starting point:

- Kick: leave at 0 dB, and check the master is peaking around −6 dB
- Clap: 3 to 6 dB below the kick
- Closed hat: quiet, 10 to 15 dB below the kick
- Open hat: slightly louder than the closed hat

Then add **Fruity Parametric EQ 2** to both hat channels and highpass them at around 400 Hz.
Hats have no useful information down low, and removing it cleans up the whole mix. Do the
same to the clap at around 200 Hz.

## Now make it yours

The four-line pattern above is correct and generic. Everything interesting comes from what
you add to it. Try these one at a time:

**A shaker or tambourine on sixteenths.** `Packs > Drums > Percussion > 707 Tamb`. Set it
quiet, and use the graph editor so it pulses.

**A rimshot on the "a" of beat 4** (step 16). `Packs > Drums > Snares > 909 Rim`. One extra
hit at the end of the bar makes the loop turn over instead of just repeating.

**An 808 cowbell** somewhere unexpected. `Packs > Drums > Percussion > 808 Cowbell`. Almost
comically effective in house.

**A second kick sample layered underneath** for weight, or an 808 kick with its decay
shortened.

**Drop the closed hats for two bars.** You will need to extend the pattern length first
(the pattern length control is at the top of the Channel Rack, set it to 2 bars and you get
32 steps).

**Take the kick out entirely for four bars** and hear what happens when it comes back.

## Save

**Ctrl+S**, and give it a real name. Not "untitled".

## What you have

> **What you have.** A loop, not a track. Twelve to sixteen hits, four samples, one bar.
> Turning it into a track is `03_first_track.md`, and it is a different and harder skill.

Before you go there, spend some time on this loop. Try to rebuild the drum pattern of a
Chicago house record by ear, using only these four sounds. Put on Marshall Jefferson's "Move
Your Body" or anything from `LISTENING/chicago_canon.md`, and copy what you hear. That single
exercise teaches more than any tutorial.
