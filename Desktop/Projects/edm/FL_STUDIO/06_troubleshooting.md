# Troubleshooting

The problems that stop beginners, and their fixes. Most of them are settings, not skill.

## No sound at all

**Check the audio device.** F10, Audio tab. Device should be **FL Studio ASIO**. If it says
something else, or "no device", that is your problem.

**Check you are in the right mode.** The toolbar shows PAT or SONG, and **L** toggles.
Pattern mode plays only the selected pattern. If your notes are in the playlist and you are
in Pattern mode, you hear nothing.

**Check the master fader** in the Mixer (F9) is not at zero.

**Check the channel is not muted.** The green LED to the left of each channel name in the
Channel Rack.

**Check Windows is not sending audio elsewhere.** Bluetooth headphones connecting mid-session
will silently steal the output.

If all else fails, **Ctrl+H** is panic, and clears stuck notes.

## Crackling, popping, or stuttering

Your buffer is too small for what your CPU is doing.

F10, Audio tab, **Show ASIO panel**, and raise the buffer to 1024 samples. Latency goes up,
which only matters when you are playing something in live.

If it still crackles, look at the CPU meter in the top-right of FL. If it is pinned:

- Right-click a channel and choose **Freeze** to render it to audio temporarily
- Reduce reverb instances by using sends instead of per-channel reverbs
- Increase the "Underruns" tolerance in audio settings
- Close Chrome, which is almost certainly using more CPU than FL is

## Everything I play is late

Latency. Lower the buffer to 256 or 128 samples in the ASIO panel, and make sure you are on
FL Studio ASIO rather than a generic Windows driver.

There is a tradeoff: lower buffer means lower latency and higher CPU load. Use a small buffer
while recording or playing in, and a large one while mixing.

## I hear only drums when I press play

Pattern mode. Press **L** to switch to Song mode.

This is the single most common confusion in FL Studio and it catches everyone at least twice.

## My pattern is in the playlist but nothing plays

Check the playlist track is not muted (the small LED on each lane).

Check the clip is actually on the timeline where the playhead is, not sitting before bar 1.

Check the pattern actually contains notes. An empty named pattern looks identical to a full
one from the outside.

## I cannot EQ one sound without affecting everything

Your channels are all going to the Master. Select them in the Channel Rack and press
**Ctrl+L** to give each its own mixer track.

## My kick and bass sound muddy together

Three fixes, in order of effort:

Move the bass notes to the offbeats so they never sound at the same instant as the kick.

Highpass the bass above the kick's fundamental with Fruity Parametric EQ 2.

Sidechain the bass to the kick. See `03_first_track.md`.

Also check that the kick is tuned near your track's key. A kick whose pitch clashes with the
bass produces a beating you will hear as mud. See `THEORY/07_mixing.md`.

## The exported file sounds different from the project

Almost always one of three things.

**The tail was cut.** In the export dialog, set **Tail** to **Leave remainder** so reverb and
delay tails finish.

**Only part exported.** If you have a time selection in the playlist, FL exports only that
range. Click an empty area of the timeline to clear the selection first.

**It is quieter.** Your project has no limiter on the master. Add Fruity Limiter with the
ceiling at −1 dB.

## The exported file is distorted

Your master is clipping. Turn individual channels down, not the master fader. Aim for the
master peaking at around −6 dB before you add a limiter.

## I lost my project

FL autosaves aggressively into `Documents\Image-Line\FL Studio\Projects` with names like
`untitled (autosaved at 21h49).flp`. There is also a `Backup` folder in there.

You already have several of these from May. Long term, save with real names and use
**Ctrl+N** deliberately (it saves an incremented version) at milestones.

## A plugin window disappeared

**F12** closes all windows, and it is easy to hit by accident. To get a plugin back, click
the channel's name button in the Channel Rack, or find it in the Browser under **Current
project**.

## The project opens with missing samples

FL stores sample paths, and moving or renaming files breaks them. When prompted, point FL at
the new location.

To avoid it: **Tools > Macros > Prepare for sharing** collects every sample the project uses
into the project folder.

## Everything is too quiet compared to real records

Commercial tracks are mastered. Yours is not.

Put Fruity Limiter on the master with the ceiling at −1 dB and raise the gain. If you want a
shortcut while learning, Soundgoodizer on the master, knob at about 30%.

Do not chase commercial loudness. Streaming services normalise it away, and a well-arranged
quiet track beats a loud badly-arranged one. See `THEORY/07_mixing.md`.

## It does not sound like the records I like

This one is not a settings problem, and it is worth saying plainly: this gap does not close
by finding the right plugin. It closes by finishing tracks, referencing against commercial
releases while you mix, and accumulating hours.

The specific thing that helps fastest is **referencing**. Load a track you admire into your
project, level-match it, and A/B constantly. You will hear exactly which of your elements is
wrong within about two minutes, which is much faster than guessing.

## Where to get help

Image-Line's own manual is genuinely good and it is installed locally (F1 opens it, though
in this version most pages redirect to the web).

The Image-Line forums are active and the moderators are the actual developers.

For anything conceptual rather than mechanical, the files in `THEORY/` in this folder will
usually cover it.
