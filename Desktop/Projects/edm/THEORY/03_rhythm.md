# Rhythm, the grid, and the drum machine

This is the file to read first if you only read one. Harmony is optional in house. Rhythm
is not.

## The grid

Everything is built on a grid of sixteen equal slots per bar. Producers call them steps.
Four beats to a bar, four steps to a beat.

Counted aloud, drummers say it like this:

```
Step      1   .   .   .   5   .   .   .   9   .   .   .  13   .   .   .
Count     1   e   &   a   2   e   &   a   3   e   &   a   4   e   &   a
```

The numbers are the beats. The "&" is halfway between two beats, which musicians call the
offbeat or the upbeat. The "e" and "a" are the sixteenth-note subdivisions.

When you open the Channel Rack in FL Studio you get exactly this: sixteen buttons in a row.
That is not a simplification for beginners. It is how the original drum machines worked and
it is how most house is still programmed.

## Tempo, and what it means in seconds

Tempo is beats per minute. House lives between 118 and 128 BPM, and 122 to 126 is the
middle of the road. Pick 124 and stop thinking about it for now.

You will need the conversion from BPM to milliseconds constantly, because delay times,
reverb decays and LFO rates all get set in musical units that you sometimes have to
translate:

    milliseconds per beat = 60000 / BPM

| At 124 BPM | Duration |
|---|---|
| One bar (4 beats) | 1935 ms |
| One beat (1/4 note) | 484 ms |
| One eighth (1/8) | 242 ms |
| One sixteenth (1/16) | 121 ms |
| Dotted eighth | 363 ms |
| Eighth triplet | 161 ms |
| 8 bars | 15.5 s |
| 16 bars | 31.0 s |
| 32 bars | 62.0 s |

That last block is why house tracks are six or seven minutes long. A track with an intro, a
build, two main sections, a breakdown and an outro, all in 32-bar blocks, is going to run
past six minutes whether you planned it or not.

The dotted eighth is worth remembering separately. A delay set to a dotted eighth produces
the classic offbeat echo you hear on everything from dub to The Edge's guitar, because it
cycles against the beat instead of on it.

## Four on the floor

The defining pattern of house, techno, disco and trance is a kick drum on every single
beat. Four kicks per bar, evenly spaced. Hence the name.

```
Count      1   e   &   a   2   e   &   a   3   e   &   a   4   e   &   a
Kick       X   .   .   .   X   .   .   .   X   .   .   .   X   .   .   .
```

It is metronomic and it is supposed to be. It gives the dancer an unambiguous pulse that
survives a loud room, and it gives the DJ something to beatmatch. It also does something
psychoacoustically useful: it masks the low end four times a bar, which is why sidechain
compression became standard (see `07_mixing.md`).

## The full basic house pattern

Add a clap on beats 2 and 4, and an open hi-hat on every offbeat.

```
Count      1   e   &   a   2   e   &   a   3   e   &   a   4   e   &   a
Kick       X   .   .   .   X   .   .   .   X   .   .   .   X   .   .   .
Clap       .   .   .   .   X   .   .   .   .   .   .   .   X   .   .   .
Open hat   .   .   X   .   .   .   X   .   .   .   X   .   .   .   X   .
```

Those three lines are a house record. Genuinely. Load a 909 kick, a 909 clap and a 909 open
hat, program that, and you have something you could play out. Everything else in this
folder is decoration on top of those twelve hits.

Two things to notice.

The clap on 2 and 4 is the **backbeat**, and it is inherited from soul, funk and rock. It is
the part of the pattern that makes house feel like Black American dance music rather than
like European electronic music. Techno often drops it entirely, and that is a large part of
why techno feels colder.

> **Why this works.** The open hat on the offbeat is the part that makes it *move*. The kick
> pulls you down four times a bar; the hat pulls you up in between. That tension is the whole
> engine. Mute the offbeat hat and listen to the loop go flat.

## The rest of the kit

Closed hi-hat, on sixteenths or eighths, fills in the texture and controls how busy the
track feels. Sixteenths make it drive. Eighths make it breathe.

Ride and crash mark transitions rather than keeping time. A crash on the downbeat of a new
section is one of the oldest tricks there is and it still works.

Percussion, meaning congas, bongos, shakers, rimshots, cowbell, tambourine, is where the
personality lives. Latin percussion has been part of house since the beginning because
disco was full of it. A shaker on sixteenths and a conga pattern that resolves every two
bars will do more for a boring loop than any amount of new synth work.

Snare rolls and hat rolls are build-up tools, not groove tools. See `06_arrangement.md`.

## The classic drum machines

You will see these names constantly and they are worth being able to identify by ear.

**Roland TR-909 (1983).** The house and techno machine. Punchy analogue kick with a very
audible click at the front, sharp snappy snare, and sampled metal for the hats and cymbals,
which is why 909 hats sound thin and crisp rather than round. The kick has a pitch envelope
you can tune, which is why 909 kicks can be made to sit in key with the bass. If you learn
one drum machine, learn this one.

**Roland TR-808 (1980).** Longer, deeper, boomier kick that sustains rather than punches.
Cowbell, congas, that specific clap. It went to hip-hop and electro more than to house. An
808 kick in a house track will muddy your low end unless you shorten it, because the tail
overlaps the next kick.

**Roland TR-707 (1985).** Sampled acoustic-ish sounds, cheaper and thinner than either.
Used all over early Chicago records because it was affordable. Worth trying when you want
something drier.

You have all three sampled in FL Studio already. `FL_STUDIO/02_first_beat.md` gives the
exact file paths.

**Roland TB-303 (1981).** Not a drum machine, a bassline synth, but it belongs in this list
because it is the other founding instrument of Chicago music. See `05_sound_design.md` and
the acid section of `02_chicago_house.md`.

## Swing

A straight grid divides each beat into two exactly equal halves. Swing delays the second
half.

The amount is usually given as a percentage, where 50% is straight and 66.7% is a full
triplet feel. House usually sits somewhere between 52% and 58%, which is a small enough
shift that you might not consciously hear it but will absolutely feel it. Garage and UK
house push much further, to 60% and beyond.

Swing is the single biggest lever between "this sounds like a machine" and "this sounds like
a groove". If your loop feels stiff, try swing before you try anything else. In FL Studio
there is a swing knob right on the Channel Rack.

> **Watch out.** Apply swing to hats and percussion, not usually to the kick. A swung kick in
> a four-on-the-floor track stops being four on the floor.

## Velocity and ghost notes

Every hit has a velocity, a number from 0 to 127 that usually maps to loudness and often
also to brightness. A pattern where every hit is at 127 sounds like a fax machine.

Real drummers hit the accented notes hard and fill the gaps with quiet ones. Copy this.
Sixteenth hats with alternating loud and quiet hits, or a repeating loud-quiet-medium-quiet
cycle, immediately sound more like playing and less like programming. Those quiet in-between
hits are called ghost notes.

The Channel Rack in FL Studio has a graph editor underneath the steps for exactly this. Use
it on every hat pattern you ever write.

## Syncopation

Syncopation means placing accents where the meter does not expect them: on offbeats, on the
"e" and "a", across the barline.

House uses it in a very specific way. The kick is completely unsyncopated, and everything
else is syncopated against it. The bass plays on the offbeats. The chord stabs land on the
"&". The percussion falls between the kick hits. This is why house feels simultaneously
rigid and loose.

The classic offbeat bass pattern, where the bass note sits exactly between the kicks:

```
Count      1   e   &   a   2   e   &   a   3   e   &   a   4   e   &   a
Kick       X   .   .   .   X   .   .   .   X   .   .   .   X   .   .   .
Bass       .   .   X   .   .   .   X   .   .   .   X   .   .   .   X   .
```

That is most of disco and a huge fraction of house. It also solves a mixing problem for
free, because the bass never plays at the same instant as the kick and so they never fight.

## Polyrhythm and the three-against-four trick

Put a percussion element on every third sixteenth instead of every fourth. It will cycle
against the bar and take three bars to line up again.

```
Count      1   e   &   a   2   e   &   a   3   e   &   a   4   e   &   a
Perc       X   .   .   X   .   .   X   .   .   X   .   .   X   .   .   X
```

This is a cheap and extremely effective way to make a loop feel like it is evolving when
nothing is actually changing. Use it sparingly. Two competing polyrhythms turn into mud.

## How to practise this

Put on any house record and clap the offbeat, the "&" between the kicks, for a full three
minutes. Not the kick. The gap. It is harder than it sounds and it will rewire how you hear
the genre.

Then open FL Studio and rebuild the drum pattern of a track you like by ear, one element at
a time. `PRACTICE/exercises.md` has a structured version of this.
