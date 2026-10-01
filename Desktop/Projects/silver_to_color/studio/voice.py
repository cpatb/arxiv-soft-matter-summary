"""Narration audio for each cue, plus word timings for captions.

Priority per cue:
  1. Your recording: <recordings>/<cue id>.wav|.flac|.m4a|.mp3  (always wins)
  2. Scratch TTS via edge-tts (Microsoft neural voices, needs network)

The TTS is a SCRATCH track for timing the cut, not the channel's voice. The monetization
research (research/monetization.md §1.4, §3.2) says narration is the single biggest lever
and that an AI *persona* presenting as a historian is the one thing to avoid. Record the
lines yourself (`python -m studio script` prints a teleprompter), drop the files into
recordings/, and re-render: the timeline re-flows to your timing automatically.

edge-tts output is cached by (text, voice, rate), so re-rendering after editing one line
only re-synthesises that line.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import subprocess
from pathlib import Path

from .timeline import Word, estimate_words

AUDIO_EXTS = (".wav", ".flac", ".m4a", ".mp3", ".ogg")


def duration_of(path: Path) -> float:
    out = subprocess.run(["ffprobe", "-v", "error", "-show_entries", "format=duration",
                          "-of", "csv=p=0", str(path)], capture_output=True, text=True,
                         check=True).stdout.strip()
    return float(out)


def find_recording(recordings: Path, cue_id: str) -> Path | None:
    for ext in AUDIO_EXTS:
        p = recordings / f"{cue_id}{ext}"
        if p.exists():
            return p
    return None


async def _edge(text: str, voice: str, rate: str, mp3: Path) -> list[Word]:
    import edge_tts
    comm = edge_tts.Communicate(text, voice, rate=rate, boundary="WordBoundary")
    words = []
    with open(mp3, "wb") as f:
        async for chunk in comm.stream():
            if chunk["type"] == "audio":
                f.write(chunk["data"])
            elif chunk["type"] == "WordBoundary":
                start = chunk["offset"] / 1e7          # 100 ns units
                words.append(Word(chunk["text"], start, start + chunk["duration"] / 1e7))
    return words


def _attach_punctuation(text: str, words: list[Word]) -> list[Word]:
    """edge-tts reports bare words; put the script's punctuation back for captions."""
    tokens = text.split()
    out, i = [], 0
    for w in words:
        while i < len(tokens) and w.text.lower() not in tokens[i].lower():
            i += 1
        if i < len(tokens):
            out.append(Word(tokens[i], w.start, w.end))
            i += 1
        else:
            out.append(w)
    return out


def synthesize(cues, recordings: Path, cache: Path, voice: dict) -> dict:
    """-> {cue_id: (audio_path, duration_s, [Word], source)} for every cue with narration."""
    cache.mkdir(parents=True, exist_ok=True)
    result = {}
    todo = []
    for c in cues:
        if not c.say:
            continue
        rec = find_recording(recordings, c.id)
        if rec:
            d = duration_of(rec)
            result[c.id] = (rec, d, estimate_words(c.say, d), "recording")
            continue
        key = hashlib.sha1(f"{c.say}|{voice['name']}|{voice['rate']}".encode()).hexdigest()[:12]
        mp3, meta = cache / f"{c.id}-{key}.mp3", cache / f"{c.id}-{key}.json"
        if mp3.exists() and meta.exists():
            words = [Word(**w) for w in json.loads(meta.read_text("utf-8"))]
            result[c.id] = (mp3, duration_of(mp3), words, "tts-cached")
        else:
            todo.append((c, mp3, meta))

    async def run_all():
        sem = asyncio.Semaphore(4)

        async def one(c, mp3, meta):
            async with sem:
                for attempt in range(3):
                    try:
                        words = await _edge(c.say, voice["name"], voice["rate"], mp3)
                        break
                    except Exception:
                        if attempt == 2:
                            raise
                        await asyncio.sleep(2 + 3 * attempt)
                words = _attach_punctuation(c.say, words)
                meta.write_text(json.dumps([w.__dict__ for w in words]), "utf-8")
                return c.id, mp3, words
        return await asyncio.gather(*(one(*t) for t in todo))

    if todo:
        for cid, mp3, words in asyncio.run(run_all()):
            result[cid] = (mp3, duration_of(mp3), words, "tts")
    return result
