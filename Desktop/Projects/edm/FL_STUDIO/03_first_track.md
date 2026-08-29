# From loop to finished track

Continues from `02_first_beat.md`. By the end you will have an exported MP3. It will not be
good. That is not the point of the exercise.

Budget two to three hours. Do it in one sitting if you can, because the failure mode is
stopping halfway.

We are in **A minor** at **124 BPM**. Everything below uses FL Studio's note naming, in which
middle C is labelled **C5**.

## 1. Organise your patterns first

You will end up with several patterns, and naming them now saves real pain later.

Press **F4** to create a new named pattern. Make these:

- `DRUMS` (the one you already made)
- `DRUMS MINIMAL` (kick and hats only, for the intro)
- `BASS`
- `CHORDS`

**F2** renames the currently selected pattern if you need to fix one.

## 2. The chords

Add a channel: click the **+** at the bottom of the Channel Rack and choose **FL Keys**. It
is FL's electric piano, and a Rhodes-type sound is the correct instrument for deep house.

Select the CHORDS pattern, select the FL Keys channel, press **F7** for the piano roll.

Set the pattern length to 2 bars, then draw two chords, one per bar:

| Bar | Chord | Notes (FL naming) |
|---|---|---|
| 1 | A minor 7 | A4, C5, E5, G5 |
| 2 | D minor 7 | A4, C5, D5, F5 |

Draw each chord as four notes stacked vertically, each lasting a full bar.

> **Why this works.** A4 and C5 stay exactly where they are. E5 moves down to D5, and G5 moves
> down to F5. Two notes move, by one step each. That is voice leading, and it is why this
> sounds smooth instead of blocky. The lazy version (A-C-E-G then D-F-A-C, everything jumping)
> sounds noticeably worse. Play both.

Both chords are minor sevenths rather than plain triads, which is what gives house its
warm, unresolved, loops-forever quality. See `THEORY/04_harmony.md`.

### Making it sound like house rather than like a piano

Right now it is a sustained chord, which is fine for a pad but dull as the main hook. Two
options:

**The stab.** Shorten every note to a sixteenth and move them onto the offbeats (steps 3, 7,
11, 15). Now the chord punches between the kicks. This is the Chicago sound.

**The pad.** Leave them long, add **Fruity Reeverb 2** to the channel, and turn the mix up
to about 30%. Now it is atmosphere.

Do both, on two separate channels, and use the stab in the busy sections and the pad in the
breakdown.

## 3. The bass

Add another channel. Two good options on your install:

**FLEX** with a bass preset, which is the easy route.

**Transistor Bass**, which is FL's Roland TB-303 emulation, and is the historically correct
choice. If you want acid, this is the instrument. Turn the resonance up high and sweep the
cutoff while it plays.

Select the BASS pattern and open the piano roll. Set the length to 2 bars.

Play the root of whatever chord is above it, on the offbeats, so the bass never collides with
the kick:

```
Count      1   e   &   a   2   e   &   a   3   e   &   a   4   e   &   a
Kick       X   .   .   .   X   .   .   .   X   .   .   .   X   .   .   .
Bass       .   .   X   .   .   .   X   .   .   .   X   .   .   .   X   .
```

Bar 1 plays **A2** (55 Hz). Bar 2 plays **D3** (73 Hz).

That is the entire disco and house bass tradition in one pattern, and it solves a mixing
problem for free: because the bass never sounds at the same instant as the kick, they cannot
fight each other.

Once that works, try the octave bounce. Alternate A2 and A3 instead of repeating A2. It
immediately has more movement.

Make sure the bass is **monophonic**, one note at a time. Chords in the bass register turn to
mud.

## 4. Route everything and set levels

Select all channels, press **Ctrl+L**, open the mixer with **F9**.

Highpass the chords at around 200 Hz with **Fruity Parametric EQ 2**. They do not need low
end, and removing it gives the bass room.

Set the kick as the loudest element and balance everything against it. Aim for the master
peaking around **−6 dB** with everything playing.

## 5. Sidechain the bass to the kick

This is the pumping effect, and it is what makes the low end work.

1. In the Mixer, select the **kick's** track.
2. Look at the small send arrows at the bottom of the mixer tracks. Right-click the arrow
   underneath the **bass's** track and choose **Sidechain to this track**.
3. Select the bass's track. Add **Fruity Limiter** to an effect slot.
4. In Fruity Limiter, switch to the **COMP** tab.
5. Set the **SIDECHAIN** knob to the sidechain number you just created (it will be 1 if this
   is your first).
6. Set the ratio to about 4:1, threshold down until you see 3 to 6 dB of gain reduction on
   each kick, attack fast, and release tuned so the bass recovers just before the next kick.

If you hear the kick doubling through the bass channel, turn that send's level knob to zero.
You want the kick's signal used as a trigger, not mixed in.

Do the same for the pad. Not for the kick itself.

If this feels like too much routing for now, the simpler alternative: right-click the bass
channel's volume knob, choose **Create automation clip**, and draw a shape that dips on every
beat. Less elegant, entirely functional, and you get exact control of the curve.

## 6. Arrange it

Switch to **Song mode** (press **L**) and open the Playlist (**F5**).

Place pattern clips using the paint tool. Here is a 128-bar arrangement, which runs about
4 minutes 8 seconds at 124 BPM. Copy it directly:

| Bars | DRUMS MINIMAL | DRUMS | BASS | CHORDS (stab) | CHORDS (pad) |
|---|---|---|---|---|---|
| 1 to 16 | yes | | | | |
| 17 to 32 | yes | | yes | | |
| 33 to 64 | | yes | yes | yes | |
| 65 to 80 | | | | | yes |
| 81 to 96 | yes | | | | yes |
| 97 to 128 | | yes | yes | yes | yes |

Bars 65 to 80 are the breakdown: the drums stop, the pad carries it. Bars 81 to 96 are the
build. Bar 97 is the drop, and it is the loudest moment in the track.

**Ctrl+B** duplicates a selection to the right, which is how you fill 32 bars without a
hundred clicks.

Notice every section boundary is on a multiple of 16 bars. This is not optional. See
`THEORY/06_arrangement.md`.

## 7. Make the build build

At the moment bars 81 to 96 are just quieter. Add tension:

**A riser.** Drag something from `Packs > Risers` onto the playlist so it ends exactly on
bar 97.

**A filter opening.** Right-click the cutoff knob on the pad's filter, create an automation
clip, and draw it rising across those 16 bars.

**A hat roll.** Make a new pattern with closed hats every 1/8 for four bars, then 1/16 for
four, then 1/32 for the last two.

**Silence.** Delete everything from the last half bar before the drop. This is the single
most effective trick available and it is free.

**A crash on bar 97.** `Packs > Drums > Cymbals`. Marks the boundary.

## 8. Master, lightly

Select the Master track in the mixer. Add **Fruity Limiter**. Set the ceiling to **−1 dB**
and push the gain until it is comfortably loud.

That is enough. Do not spend three hours here. A well-arranged track at moderate loudness
beats a badly arranged one that is very loud.

## 9. Export

First make sure the playlist knows how long the song is. If you have a time selection active
FL will export only that, so click in an empty part of the timeline to clear it.

**Ctrl+R** exports a WAV. **Ctrl+Shift+R** exports an MP3.

For sharing, MP3 at 320 kbps. For archiving or sending to anyone who will work on it further,
WAV at 24-bit.

In the export dialog, set the **Tail** option to **Leave remainder** so your reverb tails do
not get chopped off at the end.

## 10. Then listen to it somewhere else

Put it on your phone. Play it in the car. Play it on a bluetooth speaker.

It will sound different and worse, and the specific ways it sounds worse are your next
lesson. Write them down.

## What to expect

Your first track will have most of these problems: the mix is muddy, the arrangement is
repetitive, the drums are stiff, and it does not sound like the records you like. All of that
is normal and none of it is a reason to stop.

You have now done the whole loop once, from empty project to exported file. That is the thing
most people never get to. Everything from here is refinement.

`PRACTICE/8_week_plan.md` for what to do next.
