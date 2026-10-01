"""Render an episode: full 16:9 video, vertical Shorts, captions, chapters, thumbnail.

Every cue becomes one short video segment (no audio), rendered from the frame-aligned
B&W and colour masters. Segments are concatenated, then one final pass burns the text
layer (ASS: titles, lower thirds, labels, Shorts captions) and muxes the narration mix.

Paths given to ffmpeg filters are always RELATIVE to the work directory and ffmpeg runs
with that directory as cwd. Filter arguments use ':' as a separator, so an absolute
Windows path ('C:/...') inside -vf breaks the parse -- the same bug tools/pipeline.py
hit and documented.
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from . import timeline as tl
from . import voice as vo
from .episode import Episode, Short

FPS = tl.FRAME_RATE           # delivery rate; the ~16 fps film is cadenced up by fps=24
SR = 48000
BG = "0x0b0b0b"
LUFS = -14.0                  # YouTube's playback reference loudness
TRUE_PEAK = -1.5
# Deliver BT.709, tagged. The scale filter reads each input's own tag (the colour master
# says BT.709, or SMPTE 170M for a master made before 2026-10-01) and converts; the grey
# B&W master has neutral chroma, so its matrix cannot change a pixel.
TO_709 = "scale=out_color_matrix=bt709:out_range=tv,format=yuv420p"
TAGS_709 = ["-colorspace", "bt709", "-color_primaries", "bt709", "-color_trc", "bt709"]
FONT_FILES = ["georgia.ttf", "georgiab.ttf", "georgiai.ttf", "georgiaz.ttf",
              "segoeui.ttf", "segoeuib.ttf"]


@dataclass
class Canvas:
    w: int
    h: int
    short: bool
    layout: str = "single"


FULL = Canvas(1920, 1080, short=False)

# Shorts stack layout: two 3:2 panels (4:3 film cropped 11 % top and bottom) between a
# hook band at the top and the bottom ~15 % that YouTube's own UI covers.
STACK_PANEL = (1080, 720)
STACK_TOP = 200
STACK_SEAM = STACK_TOP + STACK_PANEL[1]


def run(cmd, cwd=None):
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"ffmpeg failed ({proc.returncode}):\n{' '.join(map(str, cmd))}\n"
                           f"{proc.stderr[-3000:]}")
    return proc


def film_duration(path: Path) -> float:
    return vo.duration_of(path)


# ----------------------------------------------------------------------------- video

def _src(label: str, speed: float, dur: float) -> str:
    return (f"[{label}]setpts=(PTS-STARTPTS)/{speed:.6f},fps={FPS},"
            f"tpad=stop_mode=clone:stop_duration={dur + 1:.3f}")


def _fit(canvas: Canvas, push: float, dur: float, pan: float) -> str:
    """Filter chain taking a 4:3 film frame to the canvas (no label)."""
    if canvas.short:            # crop layout: fill 1080x1920 from the centre (or pan)
        base = (f"scale=-2:{canvas.h}:flags=lanczos,"
                f"crop={canvas.w}:{canvas.h}:(iw-{canvas.w})*{pan:.3f}:0")
        fw, fh = canvas.w, canvas.h
    else:
        fh = canvas.h
        fw = int(round(fh * 4 / 3 / 2) * 2)
        base = f"scale={fw}:{fh}:flags=lanczos"
    if push > 0:
        frames = max(1, int(round(dur * FPS)))
        base += (f",scale={fw * 2}:{fh * 2}:flags=bicubic,"
                 f"zoompan=z='1+{push:.4f}*on/{frames}':x='iw/2-(iw/zoom/2)':"
                 f"y='ih/2-(ih/zoom/2)':d=1:s={fw}x{fh}:fps={FPS}")
    if not canvas.short:
        base += f",pad={canvas.w}:{canvas.h}:(ow-iw)/2:(oh-ih)/2:color={BG}"
    return base + ",setsar=1," + TO_709


def segment_cmd(ep: Episode, p: tl.Placed, canvas: Canvas, out: str) -> list[str]:
    c, dur = p.cue, p.duration
    span = dur * c.speed + 1.0
    head = ["ffmpeg", "-v", "error", "-y"]

    def inp(path, start, length):
        return ["-ss", f"{max(0.0, start):.4f}", "-t", f"{length:.4f}", "-i", str(path)]

    tail = ["-t", f"{dur:.4f}", "-an", "-c:v", "libx264", "-crf", "14", "-preset", "veryfast",
            "-pix_fmt", "yuv420p", "-r", str(FPS), "-g", str(2 * FPS)] + TAGS_709 + [out]

    if c.mode == "card":
        # Background: the colour footage, slowed, blurred and darkened under the text.
        chain = (f"{_src('0:v', 0.5, dur)},{_fit(canvas, 0, dur, 0.5)},"
                 f"gblur=sigma=24,eq=brightness=-0.22:saturation=0.55[v]")
        return head + inp(ep.color, p.film_at, dur * 0.5 + 1) + ["-filter_complex", chain,
                                                            "-map", "[v]"] + tail

    if canvas.layout == "stack":
        pw, ph = STACK_PANEL
        crop = f"scale={pw}:-2:flags=lanczos,crop={pw}:{ph}"
        chain = (f"color=c={BG}:s={canvas.w}x{canvas.h}:r={FPS}:d={dur + 1:.3f}[bg];"
                 f"{_src('0:v', c.speed, dur)},{crop}[b];"
                 f"{_src('1:v', c.speed, dur)},{crop}[c];"
                 f"[bg][b]overlay=0:{STACK_TOP}[t];[t][c]overlay=0:{STACK_SEAM},"
                 f"setsar=1,{TO_709}[v]")
        return head + inp(ep.bw, p.film_at, span) + inp(ep.color, p.film_at, span) + [
            "-filter_complex", chain, "-map", "[v]"] + tail

    if c.mode in ("color", "bw"):
        src = ep.color if c.mode == "color" else ep.bw
        chain = f"{_src('0:v', c.speed, dur)},{_fit(canvas, c.push, dur, c.pan)}[v]"
        return head + inp(src, p.film_at, span) + ["-filter_complex", chain, "-map", "[v]"] + tail

    if c.mode == "split":
        pw, ph = 936, 702
        y = (canvas.h - ph) // 2
        chain = (f"color=c={BG}:s={canvas.w}x{canvas.h}:r={FPS}:d={dur + 1:.3f}[bg];"
                 f"{_src('0:v', c.speed, dur)},scale={pw}:{ph}:flags=lanczos[b];"
                 f"{_src('1:v', c.speed, dur)},scale={pw}:{ph}:flags=lanczos[c];"
                 f"[bg][b]overlay=16:{y}[t];[t][c]overlay={canvas.w - pw - 16}:{y},"
                 f"setsar=1,{TO_709}[v]")
        return head + inp(ep.bw, p.film_at, span) + inp(ep.color, p.film_at, span) + [
            "-filter_complex", chain, "-map", "[v]"] + tail

    if c.mode == "wipe":
        # B&W plays, then colour sweeps across left to right while the film keeps
        # running. xfade consumes its second input from that input's own t=0, so the
        # colour stream is trimmed by `offset` AFTER both streams went through the same
        # seek and the same fps cadence. Seeking the colour input separately to
        # film_at + offset put the two cadences out of phase: 18-38 of 56 frames one
        # source frame off in review (2026-10-01). Offsets are whole output frames.
        offset = round(min(1.2, dur * 0.2) * FPS) / FPS
        width = round(min(2.5, max(0.6, dur * 0.4)) * FPS) / FPS
        fit = _fit(canvas, 0, dur, c.pan)
        chain = (f"{_src('0:v', c.speed, dur)},{fit},trim=0:{offset + width:.6f}[a];"
                 f"{_src('1:v', c.speed, dur)},trim=start={offset:.6f},"
                 f"setpts=PTS-STARTPTS,{fit}[b];"
                 f"[a][b]xfade=transition=wiperight:duration={width:.6f}:offset={offset:.6f}[v]")
        return head + inp(ep.bw, p.film_at, span) + inp(ep.color, p.film_at, span) + [
            "-filter_complex", chain, "-map", "[v]"] + tail

    raise ValueError(f"unhandled mode {c.mode}")


# ----------------------------------------------------------------------------- audio

def _decode(path: Path) -> np.ndarray:
    raw = subprocess.run(["ffmpeg", "-v", "error", "-i", str(path), "-ac", "1", "-ar", str(SR),
                          "-f", "f32le", "-"], capture_output=True, check=True).stdout
    return np.frombuffer(raw, np.float32)


def mix(placed: list[tl.Placed], voices: dict, length: float, out: Path,
        music: Path | None = None, music_db: float = -20.0, duck_db: float = -10.0) -> dict:
    """Narration (+ optional self-made music bed, ducked under speech) -> loudnormed WAV."""
    n = int(round(length * SR)) + SR
    track = np.zeros(n, np.float32)
    speech = np.zeros(n, np.float32)
    for p in placed:
        if p.cue.id not in voices:
            continue
        x = _decode(voices[p.cue.id][0])
        i = int(round(p.voice_start * SR))
        x = x[: max(0, n - i)]
        track[i:i + len(x)] += x
        speech[i:i + len(x)] = 1.0
    if music and Path(music).exists():
        m = _decode(Path(music))
        reps = int(np.ceil(n / max(1, len(m))))
        m = np.tile(m, reps)[:n]
        # Duck envelope: speech mask smoothed over ~300 ms so the bed breathes, not pumps.
        k = int(0.3 * SR)
        env = np.convolve(speech, np.ones(k, np.float32) / k, mode="same")
        gain = 10 ** (music_db / 20) * 10 ** (duck_db * np.clip(env, 0, 1) / 20)
        track += m * gain.astype(np.float32)
    raw = out.with_suffix(".raw.wav")
    pcm = (np.clip(track, -1, 1) * 32767).astype("<i2")
    _write_wav(raw, pcm)
    return loudnorm(raw, out)


def _write_wav(path: Path, pcm: np.ndarray):
    import wave
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes(pcm.tobytes())


def loudnorm(src: Path, dst: Path) -> dict:
    """EBU R128 two-pass (measure, then linear correction) to LUFS / TRUE_PEAK."""
    target = f"loudnorm=I={LUFS}:TP={TRUE_PEAK}:LRA=11"
    proc = run(["ffmpeg", "-hide_banner", "-i", str(src), "-af", target + ":print_format=json",
                "-f", "null", "-"])
    m = json.loads(re.search(r"\{[^{}]*\"input_i\"[^{}]*\}", proc.stderr, re.S).group(0))
    second = (f"{target}:measured_I={m['input_i']}:measured_TP={m['input_tp']}:"
              f"measured_LRA={m['input_lra']}:measured_thresh={m['input_thresh']}:"
              f"offset={m['target_offset']}:linear=true:print_format=json")
    proc = run(["ffmpeg", "-hide_banner", "-y", "-i", str(src), "-af", second, "-ar", str(SR),
                "-ac", "2", str(dst)])
    # loudnorm silently falls back to dynamic (pumping) mode when a linear gain would
    # break the true-peak ceiling; say so instead of shipping it unannounced.
    m2 = re.search(r"\{[^{}]*\"normalization_type\"[^{}]*\}", proc.stderr, re.S)
    m["normalization_type"] = json.loads(m2.group(0))["normalization_type"] if m2 else "unknown"
    if m["normalization_type"] != "linear":
        print(f"WARN loudnorm used {m['normalization_type']} mode for {dst.name} "
              f"(input {m['input_i']} LUFS, true peak {m['input_tp']} dBTP)")
    return m


# ----------------------------------------------------------------------------- outputs

def _fonts(workdir: Path):
    dst = workdir / "fonts"
    dst.mkdir(exist_ok=True)
    for name in FONT_FILES:
        src = Path("C:/Windows/Fonts") / name
        if src.exists() and not (dst / name).exists():
            shutil.copy2(src, dst / name)


def _concat(parts: list[str], workdir: Path, out: str):
    (workdir / "segments.txt").write_text("".join(f"file '{p}'\n" for p in parts), "utf-8")
    run(["ffmpeg", "-v", "error", "-y", "-f", "concat", "-safe", "0", "-i", "segments.txt",
         "-c", "copy", out], cwd=workdir)


def render_placed(ep: Episode, placed: list[tl.Placed], voices: dict, canvas: Canvas,
                  workdir: Path, final: Path, ass_text: str, jobs: int = 3) -> Path:
    workdir.mkdir(parents=True, exist_ok=True)
    _fonts(workdir)
    parts = []
    cmds = []
    for i, p in enumerate(placed):
        name = f"seg_{i:03d}_{p.cue.id}.mp4"
        cmd = segment_cmd(ep, p, canvas, name)
        stamp = workdir / (name + ".cmd")
        key = json.dumps(cmd)
        if not ((workdir / name).exists() and stamp.exists() and stamp.read_text("utf-8") == key):
            cmds.append((cmd, stamp, key))
        parts.append(name)
    from concurrent.futures import ThreadPoolExecutor

    def go(item):
        cmd, stamp, key = item
        run(cmd, cwd=workdir)
        stamp.write_text(key, "utf-8")
    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(go, cmds))
    length = tl.total(placed)
    _concat(parts, workdir, "video.mp4")
    (workdir / "text.ass").write_text(ass_text, "utf-8")
    loud = mix(placed, voices, length, workdir / "audio.wav", ep.music)
    fade_out = max(0.0, length - 0.8)
    vf = (f"fade=t=in:st=0:d=0.4,fade=t=out:st={fade_out:.3f}:d=0.8,"
          f"ass=text.ass:fontsdir=fonts")
    final.parent.mkdir(parents=True, exist_ok=True)
    run(["ffmpeg", "-v", "error", "-y", "-i", "video.mp4", "-i", "audio.wav", "-vf", vf,
         "-map", "0:v", "-map", "1:a", "-t", f"{length:.3f}", "-c:v", "libx264", "-crf", "18",
         "-preset", "medium", "-pix_fmt", "yuv420p"] + TAGS_709 + ["-c:a", "aac", "-b:a", "192k",
         "-movflags", "+faststart", str(final.resolve())], cwd=workdir)
    (workdir / "render.json").write_text(json.dumps({
        "final": str(final), "length_s": round(length, 3), "canvas": canvas.__dict__,
        "loudness_in": loud, "cues": [
            {"id": p.cue.id, "start": round(p.start, 3), "dur": round(p.duration, 3),
             "mode": p.cue.mode, "film_at": round(p.film_at, 3), "speed": p.cue.speed,
             "voice": voices.get(p.cue.id, (None, 0, [], "none"))[3]} for p in placed]},
        indent=2), "utf-8")
    return final


def _voice_map(ep: Episode, cues) -> dict:
    return vo.synthesize(cues, ep.recordings, ep.output / "voice_cache", ep.voice)


def _placed(cues, voices):
    return tl.layout(cues, {k: (v[1], v[2]) for k, v in voices.items()})


def render_full(ep: Episode, theme: tl.Theme | None = None, burn_captions: bool = False):
    theme = theme or tl.Theme()
    voices = _voice_map(ep, ep.cues)
    placed = _placed(ep.cues, voices)
    length = tl.total(placed)
    ass = tl.ass_header(FULL.w, FULL.h, theme, short=False)
    events = tl.ass_events(placed, FULL.w, FULL.h, theme, short=False)
    if not burn_captions:      # full video: captions go up as a separate SRT track
        events = "".join(e for e in events.splitlines(True) if ",Caption," not in e)
    credit = ep.credit
    if credit:                 # credit on screen for the first and last 8 seconds
        for a, b in ((0, min(8, length)), (max(0, length - 8), length)):
            events += tl._event(a, b, "Credit", r"{\fad(500,500)}" + tl.ass_escape(credit), 1)
    final = ep.output / f"{ep.id}_full.mp4"
    render_placed(ep, placed, voices, FULL, ep.output / "work_full", final, ass + events)
    stem = ep.output / f"{ep.id}_full"
    Path(str(stem) + ".srt").write_text(tl.to_srt(placed), "utf-8")
    chs = tl.chapters(placed)
    Path(str(stem) + ".description.txt").write_text(description(ep, chs), "utf-8")
    return final, placed, tl.chapter_problems(chs, length)


def render_short(ep: Episode, short: Short, theme: tl.Theme | None = None):
    theme = theme or tl.Theme()
    voices = _voice_map(ep, short.cues)
    placed = _placed(short.cues, voices)
    canvas = Canvas(1080, 1920, short=True, layout=short.layout)
    ass = tl.ass_header(canvas.w, canvas.h, theme, short=True)
    if short.layout == "stack":
        rects = ((24, STACK_TOP + 16), (24, STACK_SEAM + 16))
        events = tl.ass_events(placed, canvas.w, canvas.h, theme, short=True, hook=short.hook,
                               split_y=STACK_SEAM, stack_rects=rects)
    else:
        events = tl.ass_events(placed, canvas.w, canvas.h, theme, short=True, hook=short.hook,
                               split_y=int(canvas.h * 0.62))
    final = ep.output / "shorts" / f"{ep.id}_{short.id}.mp4"
    render_placed(ep, placed, voices, canvas, ep.output / f"work_{short.id}", final, ass + events)
    length = tl.total(placed)
    return final, length


# ----------------------------------------------------------------------------- text

def fact_sources(facts_path: Path | None) -> dict[str, list[str]]:
    """F-id -> URLs found in that fact's section of FACTS.md."""
    """F-id -> source URLs. A fact's section may give URLs directly or cite S-ids that a
    bibliography table (`| S1 | title | url |`) resolves."""
    if not facts_path or not facts_path.exists():
        return {}
    text = facts_path.read_text(encoding="utf-8")
    biblio = {m.group(1): m.group(2).rstrip(".,;") for m in
              re.finditer(r"^\|\s*(S\d+)\s*\|[^|]*\|\s*(https?://\S+?)\s*\|", text, re.M)}
    out, cur = {}, None
    for line in text.splitlines():
        m = re.match(r"^\s*(?:#+\s*|[-*]\s*)?\**\[?(F\d{2,3})\b", line)
        if m:
            cur = m.group(1)
            out.setdefault(cur, [])
        elif re.match(r"^(#+\s|---|\*\*[HM]\d)", line):
            cur = None
        if not cur:
            continue
        found = [u.rstrip(".,;") for u in re.findall(r"https?://[^\s)\]>\"'|]+", line)]
        found += [biblio[s] for s in re.findall(r"\b(S\d+)\b", line) if s in biblio]
        for url in found:
            if url not in out[cur]:
                out[cur].append(url)
    return out


def description(ep: Episode, chs) -> str:
    cited = []
    for c in ep.cues:
        for f in c.facts:
            if f not in cited:
                cited.append(f)
    sources = fact_sources(ep.facts_path)
    urls = []
    for f in cited:
        for u in sources.get(f, []):
            if u not in urls:
                urls.append(u)
    lines = [ep.description.strip(), ""] if ep.description else []
    if chs:
        lines += ["Chapters"] + [f"{tl.yt_time(t)} {name}" for t, name in chs] + [""]
    lines += ["About this restoration",
              ep.credit,
              "Colour is an AI estimate (DDColor), not a historical record: the original "
              "film is black and white, and the colours of specific objects are unknown "
              "unless a source says otherwise. Frame rate corrected to the film's true "
              "speed; flicker reduced. Narration researched from the sources below.", ""]
    if urls:
        lines += ["Sources"] + [f"- {u}" for u in urls] + [""]
    return "\n".join(lines).strip() + "\n"


def thumbnail(ep: Episode, out: Path, at: float | None = None, text: str | None = None) -> Path:
    """1280x720: B&W left, colour right, slanted seam, big serif text."""
    from PIL import Image, ImageDraw, ImageFont
    at = float(at if at is not None else ep.thumbnail.get("at", 60))
    text = text if text is not None else ep.thumbnail.get("text", ep.title)
    tmp = out.parent
    tmp.mkdir(parents=True, exist_ok=True)
    frames = {}
    for kind, src in (("bw", ep.bw), ("color", ep.color)):
        f = tmp / f"_thumb_{kind}.png"
        run(["ffmpeg", "-v", "error", "-y", "-ss", f"{at:.3f}", "-i", str(src), "-frames:v", "1",
             "-vf", "scale=1280:-2,crop=1280:720", str(f)])
        frames[kind] = Image.open(f).convert("RGB")
    W, H = 1280, 720
    mask = Image.new("L", (W, H), 0)
    d = ImageDraw.Draw(mask)
    d.polygon([(int(W * 0.55), 0), (W, 0), (W, H), (int(W * 0.40), H)], fill=255)
    img = Image.composite(frames["color"], frames["bw"], mask)
    draw = ImageDraw.Draw(img)
    draw.line([(int(W * 0.55), 0), (int(W * 0.40), H)], fill=(232, 193, 112), width=8)
    font = ImageFont.truetype("C:/Windows/Fonts/georgiab.ttf", 92)
    words = text.split()
    lines, cur = [], ""
    for w in words:
        if len(cur + " " + w) > 14 and cur:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    lines.append(cur)
    y = 40
    for line in lines[:3]:
        draw.text((48, y), line, font=font, fill=(244, 239, 230), stroke_width=7,
                  stroke_fill=(0, 0, 0))
        y += 104
    img.save(out, quality=92)
    for f in tmp.glob("_thumb_*.png"):
        f.unlink()
    return out
