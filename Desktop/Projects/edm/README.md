# EDM and Chicago house: a starting kit

Written for someone who has never made a track, does not read music, and lives in the city
where this whole thing was invented.

You played guitar and trumpet as a kid, which is more useful than you think. You already
know what a beat feels like and you have held an instrument. What you are missing is a
vocabulary, a mental model of how a dance record is built, and about six hours of clicking
around FL Studio until the layout stops feeling hostile. That is all this is.

## What is here

**THEORY/** is the reading. Nine files, roughly an hour each if you read slowly and
actually stop to listen to the examples.

| File | What it covers |
|---|---|
| `01_what_edm_is.md` | What dance music is *for*, and why that explains almost every formal choice in it |
| `02_chicago_house.md` | The Warehouse, the Music Box, Trax, acid, ghetto house, footwork. Your local history |
| `03_rhythm.md` | The grid, four on the floor, swing, the drum machine vocabulary |
| `04_harmony.md` | Notes, keys, scales, chords, bass lines. Music theory, but only the parts house uses |
| `05_sound_design.md` | Oscillators, filters, envelopes, LFOs. How a synth actually makes a sound |
| `06_arrangement.md` | Why everything happens in multiples of 8 bars, and how a track is shaped over six minutes |
| `07_mixing.md` | Gain staging, EQ, compression, sidechain, loudness |
| `08_genre_map.md` | Deep house through footwork through techno through drum and bass, with tempos and tells |
| `09_glossary.md` | Every term you will hear, defined in one line |

**FL_STUDIO/** is the software. Written against the version you actually have installed
(FL Studio 2025, v25.1.4, All Plugins bundle), with the real sample paths and the real
keyboard shortcuts checked against Image-Line's manual.

Start with `00_setup.md`, then `02_first_beat.md`. You can make a recognisable Chicago
house loop in about twenty minutes without understanding any of the theory. Do that first.
The reading makes far more sense once your hands know where things are.
`07_sampling_from_youtube.md` is separate from the sequence — read it when you want to pull
source material off YouTube, which in this genre is sooner rather than later.

**PRACTICE/** is an eight week plan and a list of drills. The plan is the single most
important file in here, because the failure mode for beginners is not lack of knowledge.
It is making forty eight-bar loops and never finishing anything.

**LISTENING/** is a canon and a scene guide. Chicago house is one of the few genres where
the founding documents are all on streaming and most of the founders are, or were until
recently, still playing in the city.

## How to use it

Suggested first week, in order:

1. `FL_STUDIO/00_setup.md`, then `FL_STUDIO/02_first_beat.md`. Get a kick, clap and hats
   looping. One hour.
2. `THEORY/01_what_edm_is.md` and `THEORY/03_rhythm.md`. Now the loop you just made will
   have names attached to it.
3. `LISTENING/chicago_canon.md`. Put on the first five tracks while you read.
4. `THEORY/02_chicago_house.md`.
5. `FL_STUDIO/03_first_track.md`. Add chords and bass, arrange it, export an MP3.

After that, `PRACTICE/8_week_plan.md` takes over.

## Two things worth saying up front

**You do not need to learn theory before you can make music.** House was built almost
entirely by people who could not read music and were working out chords by ear on
second-hand equipment. Marshall Jefferson bought a piano and taught himself the riff to
"Move Your Body" by trial and error. The theory in here is a shortcut for understanding why
the things you like sound good, not a licence you have to earn.

**The gap between your ears and your hands is the whole difficulty.** You have excellent
taste already, because you have been listening to music your entire life. Your first
tracks will be bad, and you will hear exactly how bad, and that gap is discouraging in a
way that has nothing to do with talent. Everyone goes through it. The only known treatment
is finishing things.

## The LaTeX edition

`tex/` holds a typeset version of everything above: 21 standalone chapter PDFs and one
86-page bound volume, `tex/edm_guide.pdf`.

The markdown is the source of truth. `tex/build.py` converts it with pandoc and wraps it in
`tex/preamble.tex`, so editing a markdown file and re-running the script updates the
standalone chapter and the bound volume together:

    cd tex
    python build.py --pdf

## Where the facts come from

Historical dates and details were checked against Wikipedia, Image-Line's manual, and the
City of Chicago's landmark designation report for the Warehouse. FL Studio keyboard
shortcuts come from Image-Line's official shortcut reference, not memory. Venue and record
shop details change, so anything in `LISTENING/chicago_scene_today.md` is worth a search
before you get on the train.
