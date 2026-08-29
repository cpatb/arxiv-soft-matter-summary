# Mixing, and enough about mastering

Mixing is making the parts audible together. It is the least glamorous part of production
and the one that most obviously separates a demo from a record.

The good news for you specifically: Chicago house was mixed badly, on purpose and by
accident, and it sounds fantastic. You do not need a professional mix to make good house.
You need to not make the four or five mistakes that render a track unlistenable on a club
system.

## Gain staging

Before anything else. Every channel should be at a sensible level, with the master well
below clipping.

Set your individual channels so the **master output peaks around −6 dBFS** with everything
playing. Not −0.1. Digital clipping above 0 dBFS is ugly, and leaving headroom means the
mastering stage has room to work.

If your master is clipping, turn things **down** rather than turning the master fader down.
The relationship between the channels is what you are trying to control.

Set the kick first, at around −8 to −6 dBFS on its own channel, and balance everything else
against it. The kick is the loudest thing in a house track and everything is measured
relative to it.

## The frequency map

There are only about 10 octaves of audible spectrum and every instrument is competing for
them. Knowing roughly where things live is most of EQ.

| Range | Name | What lives there |
|---|---|---|
| 20 to 60 Hz | Sub | Felt more than heard. Kick fundamental, sub bass. Club systems reproduce this; laptop speakers do not |
| 60 to 120 Hz | Bass | Kick body, bass fundamentals. The weight of the track |
| 120 to 250 Hz | Low mid | Warmth. Also where things get boomy |
| 250 to 500 Hz | Mud | Almost nothing needs boosting here. Cutting here cleans up a mix faster than anything else |
| 500 Hz to 2 kHz | Mid | Body of most instruments. Vocals, chords, snare fundamentals |
| 2 to 6 kHz | Presence | Attack, clarity, intelligibility. The ear is most sensitive here, so boosts get harsh fast |
| 6 to 12 kHz | Brightness | Hats, cymbals, sibilance, air on vocals |
| 12 to 20 kHz | Air | Sparkle. Subtle |

## EQ, in three rules

**Cut before you boost.** If two things clash, remove frequencies from the less important
one rather than boosting the more important one. Boosting adds energy and eats headroom.

**Highpass almost everything.** Every element except the kick and the bass should have a
highpass filter removing energy below where that element actually lives. Hats at 300 to 500
Hz. Chords and pads at 150 to 300 Hz. Vocals at 80 to 120 Hz. Those rumbling low frequencies
you cannot hear on your monitors are still consuming headroom and turning the club system to
soup. This single habit is the biggest improvement available to a beginner mixer.

**Carve, do not sculpt.** A narrow cut at a specific problem frequency is surgical and
useful. A broad boost across a whole region is usually a sign you picked the wrong sound.

To find a problem frequency: boost narrowly and sweep until the annoying thing gets much
worse, then cut at that frequency instead.

FL Studio's tool for this is **Fruity Parametric EQ 2**, and it is genuinely good. There is
nothing you need to buy.

## The kick and bass problem

These two live in the same octave and this is the central mixing challenge of the genre.

Four solutions, and you will usually use two of them at once:

**Separate them in frequency.** Let the kick own 50 to 80 Hz and the bass own 80 to 200 Hz,
or vice versa. Highpass the bass above the kick's fundamental.

**Separate them in time.** Put the bass on the offbeats so it never plays at the same
instant as the kick. See `03_rhythm.md`. This is the classic disco and house answer and it
is the easiest.

**Sidechain compression.** Covered below.

**Tune the kick to the key.** If your track is in A minor, tune the kick's fundamental to
an A (55 Hz, or FL note A2). A kick whose pitch clashes with the bass note produces a
dissonant beating you will hear as "muddy" without knowing why. Pitch the sample up or down
in the channel settings until it sits.

## Sidechain compression, the house sound

You know the effect already: everything ducks down when the kick hits and swells back in
between. That pumping is not a side effect. It is one of the most recognisable production
choices in dance music, and it does two jobs at once.

Functionally, it clears space in the low end at the exact moments the kick needs it, so the
kick punches through without the bass having to be quiet.

Musically, it creates rhythmic movement in sustained sounds. A pad that would otherwise sit
still starts breathing in time with the track.

Apply it to the bass, the pads, and often the chords. Never to the kick itself.

Three ways to do it in FL Studio, easiest first:

**Volume automation.** Draw a repeating shape on the channel's volume that dips on each
beat. Completely predictable, no routing, and you get exact control of the curve. Start
here.

**Fruity Limiter in COMP mode** with the kick routed to its sidechain input. The standard
method. Ratio around 4:1, fast attack, release tuned by ear so the sound recovers just
before the next kick.

**Fruity Peak Controller** on the kick channel, linked to the volume of the target channel.
Flexible, slightly more setup.

The amount matters. Deep house wants a subtle 2 to 3 dB of duck. Big-room and French house
want it obvious and extreme. Chicago house from the 1980s has essentially none of it,
because the technique came later.

## Compression, generally

A compressor reduces dynamic range: it turns down the loud parts so you can turn the whole
thing up.

**Threshold** is the level above which it starts working. **Ratio** is how hard it clamps
(4:1 means 4 dB in produces 1 dB out above the threshold). **Attack** is how fast it reacts,
and a slower attack lets the transient through, which keeps drums punchy. **Release** is how
fast it recovers.

Uses in house: glue a drum bus together, control an inconsistent vocal, add sustain and
punch to a bass, and sidechain.

The beginner error is compressing everything by default. If a sound is already consistent, a
compressor does nothing except make it less lively. Use it when you can name the problem.

## Reverb and delay, on sends

Do not put a separate reverb on every channel. It wastes CPU and, more importantly, it puts
each element in a different acoustic space, which sounds incoherent.

Instead, set up one or two reverbs on **send tracks**, and route a little of each channel to
them. Everything is then in the same room, and you control depth by how much you send.

A useful default: one short reverb (around 1 second) for drums and percussion, one long
reverb (3 seconds or more) for pads and atmosphere.

Highpass the reverb return above 300 Hz. Low frequencies in a reverb are pure mud.

Keep the kick and bass dry. Almost always.

Delay set to a dotted 1/8 note (363 ms at 124 BPM) is the classic dance music echo. FL's
**Fruity Delay 3** has tempo-synced presets.

## Stereo

The rule that matters: **keep everything below about 150 Hz in mono.** Stereo bass causes
phase problems, sounds unstable, and can partially cancel on a mono club system. FL's
**Fruity Stereo Enhancer** has a "separation" control that can do this, and Fruity
Parametric EQ 2 can be set to affect only the sides.

Above that, spread things out. Pan percussion elements slightly left and right. Widen pads
with chorus or a stereo delay. Keep the kick, bass, snare and lead vocal centred, because
those are the anchors.

Check your mix in mono occasionally. If something disappears when you sum to mono, it has a
phase problem.

## Referencing

The most valuable habit in mixing, and it costs nothing.

Load a commercially released track you admire into your project on a separate channel.
Level-match it to your track (turn the reference down, since it is mastered and much louder).
Then A/B constantly while you mix.

> **Why this works.** You will immediately hear that your kick is too quiet, your hats are too
> loud, your bass has no weight, and your mix is muddier than theirs. Referencing converts a
> vague feeling into a specific comparison, which is why it beats any amount of theory.

Pick references from the same subgenre and the same era you are aiming at.

## Loudness and what mastering is

Mastering is the final stage: making the track loud, consistent and translatable across
systems. It is a real specialisation and people do it for a living.

What you need to know now:

**LUFS** is the standard loudness measurement, roughly matched to human perception.
Streaming services normalise playback, so making your track brutally loud gains you nothing
there. Spotify and YouTube normalise to about −14 LUFS integrated, Apple Music to about −16.

**Club and DJ-oriented masters** are usually louder, often −9 to −6 LUFS integrated, because
they are competing against other records in a DJ's playlist rather than being normalised.

**True peak should stay at or below −1 dBTP** to avoid distortion when the file is converted
to MP3 or AAC.

For your first year: put **Fruity Limiter** or **Maximus** on the master, pull the ceiling to
−1 dB, push the gain until it is loud enough, and stop. Do not spend three hours on this. A
well-arranged track at a modest volume beats a badly arranged track that is very loud.

## Monitoring

You will mix on whatever you have, and that is fine, but know what your system lies about.

Laptop speakers reproduce almost nothing below 150 Hz, so you are mixing the low end blind.
Consumer headphones usually exaggerate bass, so you will under-mix it. Neither is a
disaster if you reference against commercial tracks on the same system.

Check on multiple systems: headphones, laptop, car, phone speaker. A mix that works
everywhere is the goal, and phone speakers are a brutal but honest test of whether your
groove survives without any bass at all.

Room acoustics matter more than speaker price, and a bedroom corner is the worst possible
place to put a desk. If you get monitors later, that is the point to read about acoustic
treatment.
