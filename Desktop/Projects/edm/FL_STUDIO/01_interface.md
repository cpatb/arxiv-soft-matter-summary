# The interface

FL Studio has a reputation for being confusing, and the reason is that it works differently
from every other DAW in one specific way. Once you understand that difference the rest is
just windows.

## The one idea that makes FL make sense

In most DAWs, a track is a lane that contains both an instrument and its notes. Ableton,
Logic and Pro Tools all work this way.

FL Studio splits these apart:

- A **Channel** is an instrument or a sample. It lives in the Channel Rack.
- A **Pattern** is a container of note data. It can hold notes for *any number of channels
  at once*.
- The **Playlist** is where you place pattern clips on a timeline.
- A **Mixer track** is a destination for audio, where effects live. It is a separate
  numbering system from channels, and channels have to be routed to it manually.

So one pattern might contain the kick, clap and hats together. You then drop that single
pattern clip onto the playlist wherever you want drums.

This is genuinely more flexible than the alternative and genuinely more confusing at first.

> **The catch.** Patterns are global. Which playlist lane you drop a pattern into makes no
> difference to what it sounds like. The lanes are purely visual organisation, and treating
> them as tracks is the first mistake everyone coming from Ableton makes.

## The five windows

| Window | Key | What it is for |
|---|---|---|
| Channel Rack | F6 | Your list of instruments, plus the step sequencer for drums |
| Piano Roll | F7 | Drawing notes for one channel |
| Playlist | F5 | Arranging the whole track on a timeline |
| Mixer | F9 | Levels, effects, routing |
| Browser | Alt+F8 | Samples, presets, plugins, project contents |

F12 closes all of them at once, which is the fastest way out of a mess.

### Channel Rack (F6)

A vertical list of channels. Each row has, from left to right: a green LED for mute, a pan
and volume knob, the channel's name button, and then sixteen **step buttons**.

Left-click a step to turn it on, right-click to turn it off. Those sixteen buttons are one
bar divided into sixteenth notes, exactly the grid described in `THEORY/03_rhythm.md`.

Clicking a channel's **name button** opens its settings and the plugin's interface.

At the bottom is a **+** button to add a new channel. At the top-left is a **graph editor**
toggle (a small bar-chart icon) which reveals a velocity lane for the selected channel. Use
it constantly on hi-hats.

There is a **swing** knob on the Channel Rack toolbar. It applies to the step sequencer.
Somewhere between 15% and 25% is a good house shuffle.

The step sequencer only shows one pattern at a time. Which one is set by the pattern selector
in the main toolbar.

### Piano Roll (F7)

For anything that is not a simple drum pattern. Select a channel and press F7, or right-click
the channel name and choose **Piano roll**.

Notes are drawn left to right over time, bottom to top in pitch. Left-click with the pencil
to draw, right-click to delete, drag the right edge of a note to change its length.

Remember from `THEORY/04_harmony.md` that **FL labels middle C as C5**, one octave number
higher than most other software.

Two features worth finding early. The **scale highlighting** option (in the piano roll's
options menu, under "Helpers") greys out the notes that are not in your chosen key, which
makes it very hard to play a wrong note. And the **stamp** tool in the piano roll toolbar
drops in ready-made chord shapes, including the seventh and ninth chords house is built on.

### Playlist (F5)

The timeline. Horizontal lanes hold **pattern clips**, **audio clips** and **automation
clips**.

Select a pattern in the toolbar, pick the paint tool, and click in a lane to place it. Drag
its edge to extend it. Ctrl+B duplicates whatever is selected to the right, which is how you
build 32 bars quickly.

Set up your lanes with a consistent layout: drums on lane 1, bass on 2, chords on 3, and so
on. It is only cosmetic but it stops you losing things at bar 90.

### Mixer (F9)

The Master track on the left, then numbered **Insert** tracks. Each has a fader, a
pan control, and eight **effect slots** on the right-hand panel.

> **Watch out.** Routing is manual. A new channel goes to the master by default, which means
> you cannot EQ it separately. Select the channels in the Channel Rack and press **Ctrl+L**,
> which routes each one to its own free mixer track. This is the single most useful shortcut
> in FL Studio and almost nobody finds it on their own.

At the bottom of each mixer track are small **send arrows** that route audio to other tracks.
That is how you build reverb sends and how you set up sidechain compression.

### Browser (Alt+F8)

The docked panel on the left. Its most useful folder for you is **Packs**, which is the
bundled sample library:

    Packs
      Drums          Kicks, Hats, Snares, Percussion, Cymbals, Toms, Kits
      Drums (ModeAudio)
      Instruments    Bass, Guitar, Keyboard, Orchestral
      Legacy         the older FL sample library
      Loops
      Risers         build-up sweeps
      SFX            impacts and noise
      Shapes         automation curve presets
      Vocals
      FLEX           presets for the FLEX synth

You can drag anything from the browser straight onto the Channel Rack or the playlist.
Single-clicking a sample previews it.

## Pattern mode versus Song mode

The toolbar has a **PAT / SONG** switch, and **L** toggles it.

**Pattern mode** loops the currently selected pattern. This is where you work on a groove.

**Song mode** plays the playlist from the timeline. This is where you arrange.

If you press play and hear only drums when your project has chords, you are almost certainly
in Pattern mode with the drum pattern selected. This is the most common confusion in the
entire program.

## Transport and patterns

**Space** starts and stops. **Ctrl+Space** pauses instead of returning to the start.

**F4** creates the next empty pattern and asks you to name it. Name your patterns. A project
with "Pattern 1" through "Pattern 23" is unusable by the time you are arranging.

**F2** renames the selected pattern.

**Ctrl+H** is panic, and stops all sound. Useful when something gets stuck.

## Tools

In the Playlist and Piano Roll, the toolbar has a set of tools with single-key shortcuts:

**P** draw (pencil), **B** paint (brush, for repeated placement), **D** delete, **E** select,
**Z** zoom, **C** slice, **T** mute, **Y** playback, and **S** for slip edit in the playlist.

Hold **Alt** while dragging to bypass snapping, which you will need constantly for
fine-tuning timing.

## What to ignore for now

FL has a great deal you do not need yet: Patcher, the piano roll's arpeggiator and Riff
Machine, Gross Beat, Edison, the video player, ReWire, performance mode, and about thirty
of the synths.

None of it is wasted, but exploring it is not learning. Ignore it until you have finished
three tracks.
