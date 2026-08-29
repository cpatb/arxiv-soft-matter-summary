# Setup: the first twenty minutes

Written against what is actually installed on this machine: **FL Studio 2025, version
25.1.4**, with the All Plugins bundle. Every path and shortcut below has been checked
against your install and against Image-Line's official shortcut reference rather than
recalled from memory.

You already opened FL once, on 8 May, and worked for about two hours before stopping. Worth
knowing that those projects are still there under
`Documents\Image-Line\FL Studio\Projects` if you want to look at what past-you was doing.

## 1. Audio settings, before anything else

If the audio is misconfigured, everything you play will be late, crackly, or silent, and
you will conclude that you are bad at this rather than that a dropdown is wrong.

Press **F10** to open Settings, then click the **Audio** tab.

Set **Device** to **FL Studio ASIO**. It is already installed on this machine. ASIO is the
low-latency Windows audio driver standard, and the generic Windows driver will give you
audible delay between pressing a key and hearing a note. FL Studio ASIO has the useful
property of letting other applications keep using your sound card at the same time, so you
can have a YouTube tutorial open while FL is running.

If you later buy an audio interface, use its own ASIO driver instead.

Click **Show ASIO panel** and set the buffer to **512 samples**, which is about 11 ms of
latency. If you get crackling, go up to 1024. If you are playing an instrument in and it
feels sluggish, go down to 256.

Set the sample rate to **44100 Hz**. There are arguments for 48 kHz. None of them will
matter to you for a year.

## 2. Project settings

**Tempo.** The number in the top toolbar. Right-click it and choose **Type in value**, then
enter **124**. That is the middle of the house range and a good default.

**Save your project immediately** with Ctrl+S, and give it a real name. FL autosaves, but
it autosaves into files called "untitled (autosaved at 21h49)", which is how you end up with
the folder you currently have.

A warning about a shortcut that will surprise you: in FL Studio **Ctrl+N is not "new
file"**. It is "save new version", which increments your filename. Use **File > New** from
the menu instead.

## 3. Make a starter template

You will do the same setup at the start of every project. Do it once.

Build a project containing:

- Tempo at 124
- Four drum channels loaded (kick, clap, closed hat, open hat) and routed to their own
  mixer tracks
- A Fruity Limiter on the master with the ceiling at −1 dB
- Two mixer tracks set up as reverb sends, one short and one long

Then save it into:

    C:\Users\compl\Documents\Image-Line\FL Studio\Projects\Templates\

That folder exists and is currently empty. Anything you save there appears under **File >
New from template**. Call it something like `house_124.flp`.

This is a bigger time-saver than it looks. Removing twenty minutes of setup friction is the
difference between opening FL when you have a spare hour and not bothering.

## 4. Where your files live

| What | Where |
|---|---|
| The program | `C:\Program Files\Image-Line\FL Studio 2025` |
| Bundled samples | `...\FL Studio 2025\Data\Patches\Packs` (shown as **Packs** in the browser) |
| Your projects | `C:\Users\compl\Documents\Image-Line\FL Studio\Projects` |
| Your templates | `...\FL Studio\Projects\Templates` |
| Rendered audio | `...\FL Studio\Audio` by default |

Worth adding your own sample folder to the browser at some point: **Settings (F10) > File**,
then point one of the extra search folders at wherever you keep downloaded packs. They then
show up in the browser alongside the bundled ones.

## 5. Turn on two things that will save you

**Options > Undo history** is not enabled by default in a useful depth. Set the undo levels
higher in Settings > General.

**The metronome** is Ctrl+M. You will not need it for step-sequenced drums, but you will the
moment you try to play anything in.

## 6. A note on editions

FL Studio comes in Fruity, Producer, Signature and All Plugins editions. Your install has
the All Plugins generator set: Sytrus, Harmor, Morphine, Ogun, Toxic Biohazard, Transistor
Bass, Kepler, Drumaxx, Slicex, DirectWave and the rest are all present.

That means every tutorial you find online will work for you, and you never need to buy a
plugin to follow along. It also means you have about forty synthesisers and will be tempted
to try all of them. Do not. `05_plugins.md` lists the six that matter.

## Next

`01_interface.md` to learn the windows, or skip straight to `02_first_beat.md` and make
noise. Making noise first is the better order.
