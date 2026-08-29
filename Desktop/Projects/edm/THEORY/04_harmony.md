# Music theory, but only the parts house uses

You do not need this to start. You need it to stop guessing.

I have written this assuming you are comfortable with logarithms and ratios, because that
makes the first section take five minutes instead of an hour.

## Notes and why there are twelve of them

Pitch is frequency. Double a frequency and you get a note that sounds like the same note,
higher. That 2:1 ratio is the octave, and it is the one interval that is essentially
universal across musical cultures.

Western music divides the octave into twelve equal steps on a logarithmic scale. One step
is a semitone, and the frequency ratio is

    2^(1/12) = 1.059463...

so twelve of them multiply to exactly 2. This system is called twelve-tone equal
temperament. Every semitone is the same size, which means you can transpose anything to any
key and it still works. The price is that no interval except the octave is exactly a
small-integer ratio: a perfect fifth in equal temperament is 2^(7/12) = 1.4983 instead of
1.5. Nobody can hear the difference in a club.

The reference point is A4 = 440 Hz. From any MIDI note number n:

    f = 440 × 2^((n − 69) / 12)

The twelve notes are named A, A#/Bb, B, C, C#/Db, D, D#/Eb, E, F, F#/Gb, G, G#/Ab. The
sharps and flats are the black keys. The reason the naming is irregular (no E# or B#) is
historical and not worth your time.

### An FL Studio trap, worth knowing now

> **The trap.** FL Studio labels middle C as **C5**. Most other software, and most theory
> writing, calls it C4. FL's octave numbers are one higher than the scientific convention
> throughout, so every note name you read in a tutorial written for another program is off by
> an octave.

This matters when you follow a tutorial written for Ableton, or when you read a synth manual.
It is not a bug, just an inherited quirk. Practical version, in FL's own labels:

| Layer | Rough range in FL note names |
|---|---|
| Sub bass and kick fundamentals | C2 to C3 |
| Bassline | C2 to C4 |
| Chords | C4 to C6 |
| Lead / melody | C5 to C7 |

And for tuning a kick drum or a sub, the frequencies in FL's labelling:

| FL note | Hz |
|---|---|
| C2 | 32.7 |
| D2 | 36.7 |
| E2 | 41.2 |
| F2 | 43.7 |
| G2 | 49.0 |
| A2 | 55.0 |
| B2 | 61.7 |
| C3 | 65.4 |

A typical house kick has its fundamental between 45 and 65 Hz, so somewhere around G2 to
C3. See `07_mixing.md` for why tuning it to your track's key is worth doing.

## Intervals

The distance between two notes, counted in semitones. These are the building blocks of
everything else, and the names are worth memorising because chord names are built from them.

| Semitones | Name | Sounds like |
|---|---|---|
| 1 | minor 2nd | tense, adjacent |
| 2 | major 2nd | a step |
| 3 | minor 3rd | sad |
| 4 | major 3rd | happy |
| 5 | perfect 4th | open, suspended |
| 6 | tritone | unstable |
| 7 | perfect 5th | hollow, strong, stable |
| 8 | minor 6th | melancholy |
| 9 | major 6th | warm |
| 10 | minor 7th | jazzy, unresolved |
| 11 | major 7th | dreamy, floating |
| 12 | octave | the same note |

The minor third versus major third distinction is the single most important one. It is the
difference between sad and happy, and it is one semitone.

## Scales

A scale is a subset of the twelve notes that you agree to mostly stay inside. Picking one
is what makes a piece sound coherent rather than random.

Written as semitone steps between consecutive notes:

| Scale | Steps | Character |
|---|---|---|
| Major | 2 2 1 2 2 2 1 | bright, resolved |
| Natural minor | 2 1 2 2 1 2 2 | sad, the default for house |
| Dorian | 2 1 2 2 2 1 2 | minor but with a lifted 6th; funky, soulful |
| Phrygian | 1 2 2 2 1 2 2 | dark, Spanish, used in hard techno |
| Harmonic minor | 2 1 2 2 1 3 1 | dramatic, exotic; use with care |

Starting on A, the natural minor scale is A B C D E F G, all white keys. That is why A
minor is where beginners should start: no accidentals to keep track of.

Dorian is worth special attention. It is natural minor with the sixth note raised one
semitone (A B C D E **F#** G). That single note is most of the difference between "sad" and
"cool". A huge amount of deep house and disco sits in Dorian.

## Keys

A key is a scale plus a home note. "A minor" means the natural minor scale starting from A,
with A as the note everything wants to return to.

Two facts that save time:

Every major key has a **relative minor** using exactly the same notes, starting three
semitones lower. C major and A minor contain identical notes. So do F major and D minor, G
major and E minor.

House is overwhelmingly in **minor keys**. A minor, F minor, G minor, C minor, D minor.
Partly because minor sounds better against dark club sound systems, partly because the whole
genre inherited its emotional register from disco and soul, which are full of yearning
rather than triumph.

Practical advice: write everything in A minor for your first month. When it starts feeling
constraining, move to F minor, which is the other workhorse key of the genre.

## Chords

A chord is three or more notes at once. Named by their root and their pattern of intervals
above it, in semitones:

| Chord | Semitones above root | Example on A |
|---|---|---|
| minor triad | 0, 3, 7 | A C E |
| major triad | 0, 4, 7 | A C# E |
| diminished | 0, 3, 6 | A C Eb |
| sus2 | 0, 2, 7 | A B E |
| sus4 | 0, 5, 7 | A D E |
| minor 7th (m7) | 0, 3, 7, 10 | A C E G |
| major 7th (maj7) | 0, 4, 7, 11 | A C# E G# |
| dominant 7th (7) | 0, 4, 7, 10 | A C# E G |
| minor 9th (m9) | 0, 3, 7, 10, 14 | A C E G B |
| major 9th (maj9) | 0, 4, 7, 11, 14 | A C# E G# B |
| add9 | 0, 4, 7, 14 | A C# E B |
| minor 6th (m6) | 0, 3, 7, 9 | A C E F# |

### Why house lives on 7ths and 9ths

A plain triad is stark and finished. It states a fact and sits down.

Add the 7th and the chord becomes slightly unresolved: pleasant, hanging, in no hurry. Add
the 9th and it becomes lush. This is the sound of Rhodes pianos, of Larry Heard, of
basically all deep house.

> **Why this works.** House harmony usually does not *move* anywhere. A rock song sets up
> tension and resolves it. A house track loops two chords for six minutes. A plain triad loop
> gets boring because each chord keeps announcing that it has arrived. A 7th chord never quite
> lands, so it can repeat indefinitely without demanding a resolution that never comes.

If your chords sound amateurish, this is almost always the reason. Add the 7th.

## Building chords from a scale

Take your scale, and build a chord on each note using only notes from that scale. In A
natural minor:

| Degree | Chord | Notes | Roman numeral |
|---|---|---|---|
| 1 | A minor 7 | A C E G | i7 |
| 2 | B minor 7 flat 5 | B D F A | ii°7 |
| 3 | C major 7 | C E G B | IIImaj7 |
| 4 | D minor 7 | D F A C | iv7 |
| 5 | E minor 7 | E G B D | v7 |
| 6 | F major 7 | F A C E | VImaj7 |
| 7 | G dominant 7 | G B D F | VII7 |

Roman numerals are the standard shorthand. Lowercase means minor, uppercase means major.
The point of the notation is that it is key-independent: "i, VI, III, VII" describes a
progression you can transpose anywhere.

## Progressions that actually appear in house

Written in A minor. Transpose freely.

**The two-chord vamp: Am7 to Dm7 (i to iv).** Two bars each. This is the deep house default
and it is not a cop-out. "Can You Feel It" territory. Because both chords are minor 7ths and
they share notes, the loop feels like breathing rather than like a chord change.

**i, VI, III, VII: Am7, Fmaj7, Cmaj7, G.** The most common minor progression in all popular
music. Melancholy but forward-moving. One bar each, or two if you want it to feel patient.

**i, VII, VI, VII: Am7, G, Fmaj7, G.** Descends then turns around. Slightly more anthemic.

**ii, V, i: Bm7b5, E7, Am7.** Borrowed from jazz. Use it as a four-bar turnaround at the end
of an eight-bar phrase to make things feel less mechanical.

**ii, V, I in the relative major: Dm7, G7, Cmaj7.** The single most common progression in
jazz, and it shows up all over deep house and nu-disco. Warmer than the minor options.

**i to v: Am7 to Em7.** Very static, very hypnotic. Good for techno-leaning house where you
want almost no harmonic information.

**In F minor, the other house key: Fm7 to Bbm7,** or **Fm7, Dbmaj7, Abmaj7, Eb.** Same
shapes, darker register.

## Voicing and inversions

Which notes you play matters less than where you put them.

An **inversion** is the same chord with a different note at the bottom. Am7 as A-C-E-G is
root position. C-E-G-A is the first inversion. Same chord, different colour, and crucially,
smoother movement between chords.

The rule that fixes most beginner chord work: **when moving between chords, move each note
as little as possible.** Instead of playing Am (A C E) then F (F A C) with both hands
jumping down, play Am as A-C-E and F as A-C-F. Only one note moved. This is called voice
leading and it is why professionally written chords sound smooth and yours sound blocky.

Two more practical rules:

Do not stack notes below about C4 in FL's labelling. Low intervals turn to mud fast. If you
want low weight, that is the bass part's job, and it should be a single note.

Leave out notes. A four-note chord played by four fingers is a piano habit. On a synth pad
you often want just the root, the 3rd and the 7th, with the 5th dropped entirely, because
the 5th adds almost no information and takes up space.

## The bass

In house the bassline is the most important melodic element and it usually does one of
three jobs.

**Root notes on the offbeat.** The bass plays the root of the current chord, on the "&" of
every beat, out of the way of the kick. Simple, disco, works every time. See the pattern in
`03_rhythm.md`.

**Octave bounce.** Alternate the root note between two octaves. Root low, root high, root
low, root high. Instantly adds movement without adding harmony.

**A proper line.** Root, 5th, octave, minor 7th, back down. Now the bass is a melody. This
is where funk and disco basslines come from, and it is the hardest of the three to do
without it becoming busy.

Whatever you do, the bass should be **monophonic**, meaning one note at a time. Chords in
the bass register are muddy. Almost every "why does my low end sound bad" question comes
down to either two-note bass chords or a bass that overlaps the kick in time.

## Melody, briefly

House often has no melody at all, and this is fine. When it does, the melody is usually
short, repetitive, and made of notes already in the chords underneath.

The most reliable way to write one that works: pick notes from the current chord for the
strong beats, and pass through other scale notes on the way between them. That is the entire
technique.

The **pentatonic scale** is a cheat code worth knowing. It is the minor scale with the two
awkward notes removed: A C D E G. Every note in it sounds acceptable over every chord in the
key. If you improvise using only those five notes in A minor you will not play a wrong note.
Set your FL piano roll scale highlight to it and click around.

## Picking a key, and harmonic mixing

For your own tracks, key choice barely matters musically. A minor and F minor are fine
forever.

It matters for DJs, because two records in clashing keys sound bad played together. The
standard tool is the **Camelot wheel**, a relabelling of the circle of fifths where each key
gets a code like 8A or 5B. Keys with adjacent numbers, or the same number with the other
letter, mix well. Software like Mixed In Key does the detection automatically, and Beatport
lists keys on every track.

You do not need this yet. Know that it exists so the numbers on Beatport stop being
mysterious.

## The honest summary

Almost every house track ever made can be described as: one minor key, two to four seventh
chords looping every two or four bars, a monophonic bassline on the roots, and a melody made
of chord tones. That is the whole harmonic language.

The interest is elsewhere. It is in the groove, the sound design and the arrangement. Do not
spend six months on theory when the actual differentiator is whether your hi-hats swing.
