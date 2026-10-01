"""Pure timing logic: cue durations -> absolute timeline, captions, chapters, ASS text.

No ffmpeg, no network, so all of it is unit-tested directly.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from .episode import Cue

FRAME_RATE = 24              # delivery frame rate; render.FPS reads this

# Caption chunking. Full video: subtitle-style lines. Shorts: big 1-3 word bursts.
FULL_MAX_CHARS = 42          # Netflix/BBC single-line guidance
FULL_MAX_SECONDS = 6.0
SHORT_MAX_WORDS = 3
SHORT_MAX_CHARS = 18


@dataclass
class Word:
    text: str
    start: float      # seconds, relative to the start of its own audio file
    end: float


@dataclass
class Placed:
    cue: Cue
    start: float               # timeline seconds
    film_at: float             # film seconds where this cue's footage starts
    duration: float
    voice_start: float         # timeline seconds where the line begins
    voice_duration: float
    words: list[Word] = field(default_factory=list)   # absolute timeline times

    @property
    def end(self) -> float:
        return self.start + self.duration


def layout(cues: list[Cue], voice: dict[str, tuple[float, list[Word]]]) -> list[Placed]:
    """Place cues back to back. A cue lasts lead + line + pad, or `hold` if longer."""
    t = 0.0
    cursor = 0.0          # where the footage has got to, in film seconds
    placed = []
    for c in cues:
        dur_v, words = voice.get(c.id, (0.0, []))
        dur = max(c.hold, (c.lead + dur_v + c.pad) if dur_v > 0 else 0.0)
        # Whole output frames: each segment is cut with -t at FRAME_RATE, which rounds,
        # and the rounding errors would random-walk picture away from the audio and
        # captions placed at exact times (~0.1 s rms after 60 cues, per review).
        dur = math.ceil(dur * FRAME_RATE - 1e-6) / FRAME_RATE
        if dur <= 0:
            raise ValueError(f"cue {c.id} has zero duration")
        film_at = cursor if c.at is None else c.at
        v0 = t + c.lead
        placed.append(Placed(c, t, film_at, dur, v0, dur_v,
                             [Word(w.text, v0 + w.start, v0 + w.end) for w in words]))
        if c.mode != "card":           # a card's background does not advance the ride
            cursor = film_at + dur * c.speed
        t += dur
    return placed


def total(placed: list[Placed]) -> float:
    return placed[-1].end if placed else 0.0


def estimate_words(text: str, duration: float) -> list[Word]:
    """Word timings for a recording with no aligner: proportional to character count.

    An estimate, and labelled as one wherever it is used. Good enough for Shorts burst
    captions, where a chunk is on screen ~1 s; upload the SRT to YouTube and let its own
    aligner re-time it for the full video.
    """
    tokens = text.split()
    if not tokens or duration <= 0:
        return []
    weights = [len(tok) + 1 for tok in tokens]
    scale = duration / sum(weights)
    out, t = [], 0.0
    for tok, w in zip(tokens, weights):
        out.append(Word(tok, t, t + w * scale))
        t += w * scale
    return out


# ----------------------------------------------------------------------------- captions

def _chunks(words: list[Word], max_chars: int, max_words: int | None,
            max_seconds: float) -> list[list[Word]]:
    chunks, cur = [], []
    for w in words:
        trial = " ".join(x.text for x in cur + [w])
        too_long = len(trial) > max_chars or (max_words and len(cur) >= max_words)
        too_slow = cur and (w.end - cur[0].start) > max_seconds
        if cur and (too_long or too_slow):
            chunks.append(cur)
            cur = []
        cur.append(w)
        filled = len(" ".join(x.text for x in cur))
        if (re.search(r"[.!?;:]$", w.text) and filled > max_chars // 3) or \
                (w.text.endswith(",") and filled > 0.6 * max_chars):
            chunks.append(cur)
            cur = []
    if cur:
        chunks.append(cur)
    return chunks


def caption_chunks(placed: list[Placed], short: bool) -> list[list[Word]]:
    out = []
    for p in placed:
        if short:
            out += _chunks(p.words, SHORT_MAX_CHARS, SHORT_MAX_WORDS, 2.0)
        else:
            out += _chunks(p.words, FULL_MAX_CHARS, None, FULL_MAX_SECONDS)
    return out


def srt_time(t: float) -> str:
    ms = int(round(max(0.0, t) * 1000))
    h, ms = divmod(ms, 3_600_000)
    m, ms = divmod(ms, 60_000)
    s, ms = divmod(ms, 1000)
    return f"{h:02d}:{m:02d}:{s:02d},{ms:03d}"


def to_srt(placed: list[Placed]) -> str:
    lines = []
    for i, chunk in enumerate(caption_chunks(placed, short=False), 1):
        lines += [str(i), f"{srt_time(chunk[0].start)} --> {srt_time(chunk[-1].end)}",
                  " ".join(w.text for w in chunk), ""]
    return "\n".join(lines)


def chapters(placed: list[Placed]) -> list[tuple[float, str]]:
    """YouTube chapter list. The first chapter is forced to 0:00, as YouTube requires."""
    out = [(p.start, p.cue.chapter) for p in placed if p.cue.chapter]
    if out:
        out[0] = (0.0, out[0][1])
    return out


def yt_time(t: float) -> str:
    t = int(t)
    h, rem = divmod(t, 3600)
    m, s = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


def chapter_problems(chs: list[tuple[float, str]], length: float) -> list[str]:
    """YouTube shows chapters only if: first at 0:00, >= 3 of them, each >= 10 s."""
    probs = []
    if len(chs) < 3:
        probs.append(f"{len(chs)} chapters; YouTube needs >= 3")
    bounds = [t for t, _ in chs] + [length]
    for (t, name), nxt in zip(chs, bounds[1:]):
        if nxt - t < 10:
            probs.append(f"chapter {name!r} is {nxt - t:.1f}s; YouTube needs >= 10 s")
    return probs


# ----------------------------------------------------------------------------- ASS

def ass_time(t: float) -> str:
    cs = int(round(max(0.0, t) * 100))
    h, cs = divmod(cs, 360000)
    m, cs = divmod(cs, 6000)
    s, cs = divmod(cs, 100)
    return f"{h}:{m:02d}:{s:02d}.{cs:02d}"


def ass_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("{", "(").replace("}", ")").replace("\n", "\\N")


def ass_colour(hex_rgb: str, alpha: int = 0) -> str:
    """'#RRGGBB' -> ASS '&HAABBGGRR'."""
    h = hex_rgb.lstrip("#")
    r, g, b = h[0:2], h[2:4], h[4:6]
    return f"&H{alpha:02X}{b}{g}{r}".upper()


@dataclass
class Theme:
    serif: str = "Georgia"
    sans: str = "Segoe UI"
    accent: str = "#E8C170"      # warm silver-gelatin gold
    ink: str = "#F4EFE6"
    shadow: str = "#000000"


def ass_header(width: int, height: int, theme: Theme, short: bool) -> str:
    a, ink, sh = ass_colour(theme.accent), ass_colour(theme.ink), ass_colour(theme.shadow)
    s = height / 1080.0 if not short else width / 1080.0
    cap_size = int(round((84 if short else 44) * s))
    styles = [
        # name, font, size, primary, secondary, outline, back, bold, italic, border, outline, shadow, align, mL, mR, mV
        ("Caption", theme.sans, cap_size, ink, a, sh, ass_colour("#000000", 0x60), -1, 0, 1,
         int(6 * s) if short else int(3 * s), int(2 * s), 2, 60, 60,
         int(840 * s) if short else int(60 * s)),
        ("Title", theme.serif, int((96 if short else 88) * s), ink, a, sh, sh, -1, 0, 1, int(3 * s),
         int(3 * s), 5, 80, 80, 0),
        ("Subtitle", theme.serif, int((52 if short else 44) * s), a, a, sh, sh, 0, -1, 1, int(2 * s),
         int(2 * s), 5, 80, 80, 0),
        ("Lower", theme.serif, int(46 * s), ink, a, sh, ass_colour("#000000", 0x40), 0, -1, 3,
         int(14 * s), 0, 1, int(70 * s), 60, int(90 * s)),
        ("Label", theme.sans, int((40 if short else 30) * s), a, a, sh, sh, -1, 0, 1, int(2 * s),
         int(2 * s), 7, 0, 0, 0),
        ("Hook", theme.serif, int(70 * s), ink, a, sh, ass_colour("#000000", 0x30), -1, 0, 3,
         int(18 * s), 0, 8, 50, 50, int(70 * s)),
        ("Credit", theme.sans, int(24 * s), ink, a, sh, sh, 0, 0, 1, int(1.5 * s), 0, 3, 30, 30,
         int(24 * s)),
    ]
    lines = ["[Script Info]", "ScriptType: v4.00+", f"PlayResX: {width}", f"PlayResY: {height}",
             "WrapStyle: 0", "ScaledBorderAndShadow: yes", "", "[V4+ Styles]",
             "Format: Name, Fontname, Fontsize, PrimaryColour, SecondaryColour, OutlineColour, "
             "BackColour, Bold, Italic, Underline, StrikeOut, ScaleX, ScaleY, Spacing, Angle, "
             "BorderStyle, Outline, Shadow, Alignment, MarginL, MarginR, MarginV, Encoding"]
    for (name, font, size, pri, sec, out, back, bold, ital, border, outline, shadow,
         align, ml, mr, mv) in styles:
        lines.append(f"Style: {name},{font},{size},{pri},{sec},{out},{back},{bold},{ital},0,0,"
                     f"100,100,0,0,{border},{outline},{shadow},{align},{ml},{mr},{mv},1")
    lines += ["", "[Events]",
              "Format: Layer, Start, End, Style, Name, MarginL, MarginR, MarginV, Effect, Text"]
    return "\n".join(lines) + "\n"


def _event(start, end, style, text, layer=0, margin_v=0):
    return (f"Dialogue: {layer},{ass_time(start)},{ass_time(end)},{style},,0,0,"
            f"{margin_v},,{text}\n")


def ass_events(placed: list[Placed], width: int, height: int, theme: Theme, short: bool,
               hook: str = "", credit: str = "", split_y: int | None = None,
               stack_rects: tuple | None = None) -> str:
    """All burned-in text for one output: captions, titles, lower thirds, labels."""
    ev = []
    acc = ass_colour(theme.accent)
    ink = ass_colour(theme.ink)
    end_all = total(placed)

    # Cards.
    for p in placed:
        c = p.cue
        if c.mode == "card":
            fade = r"{\fad(500,400)}"
            if c.card_subtitle:
                ev.append(_event(p.start, p.end, "Title",
                                 fade + rf"{{\pos({width // 2},{int(height * 0.44)})}}"
                                 + ass_escape(c.card_title or ""), 2))
                ev.append(_event(p.start + 0.4, p.end, "Subtitle",
                                 fade + rf"{{\pos({width // 2},{int(height * 0.56)})}}"
                                 + ass_escape(c.card_subtitle), 2))
            else:
                ev.append(_event(p.start, p.end, "Title", fade + ass_escape(c.card_title or ""), 2))
        if c.lower_third:
            ev.append(_event(p.start + 0.3, min(p.end, p.start + 6.0), "Lower",
                             r"{\fad(400,400)}" + ass_escape(c.lower_third), 2))

    # Before/after labels, wherever both pictures are on screen.
    for p in placed:
        if p.cue.mode == "split" and not short:
            y = int(height * 0.16)
            ev.append(_event(p.start, p.end, "Label",
                             rf"{{\fad(300,300)\pos({int(width * 0.03)},{y})}}AS SCANNED", 1))
            ev.append(_event(p.start, p.end, "Label",
                             rf"{{\fad(300,300)\pos({int(width * 0.51)},{y})}}RESTORED + COLORIZED", 1))
        if short and stack_rects and p.cue.mode != "card":
            (bx, by), (cx, cy) = stack_rects
            ev.append(_event(p.start, p.end, "Label", rf"{{\pos({bx},{by})}}ORIGINAL", 1))
            ev.append(_event(p.start, p.end, "Label", rf"{{\pos({cx},{cy})}}RESTORED", 1))

    # Captions. Shorts: burst chunks with the spoken word in the accent colour.
    if short:
        chunks = caption_chunks(placed, short=True)
        for ci, chunk in enumerate(chunks):
            # A chunk may linger 0.15 s after its last word, but never into the next
            # chunk: two bursts on screen at once render as one garbled line.
            nxt = chunks[ci + 1][0].start if ci + 1 < len(chunks) else float("inf")
            for k, w in enumerate(chunk):
                t1 = chunk[k + 1].start if k + 1 < len(chunk) else min(w.end + 0.15, nxt)
                parts = []
                for m, x in enumerate(chunk):
                    tok = ass_escape(x.text.upper())
                    parts.append(rf"{{\c{acc}}}{tok}{{\c{ink}}}" if m == k else tok)
                pos = rf"{{\an5\pos({width // 2},{split_y})}}" if split_y else ""
                ev.append(_event(w.start, t1, "Caption", pos + " ".join(parts), 3))
        if hook:
            ev.append(_event(0, end_all, "Hook", ass_escape(hook), 4))
    else:
        for chunk in caption_chunks(placed, short=False):
            ev.append(_event(chunk[0].start, chunk[-1].end, "Caption",
                             ass_escape(" ".join(w.text for w in chunk)), 3))
    if credit:
        ev.append(_event(0, end_all, "Credit", r"{\alpha&H40&}" + ass_escape(credit), 1))
    return "".join(ev)
