"""Episode spec: one YAML file describes a full video and its Shorts.

Minimal example (see episodes/ep01_market_street/episode.yaml for a real one):

    id: ep01_market_street
    title: A Trip Down Market Street, 1906
    film:
      bw: E:/silver_to_color/ep01/bw_master.mp4
      color: E:/silver_to_color/ep01/color_master.mp4
      credit: "Film: Miles Brothers, 1906. Scan: Prelinger Archives / Internet Archive"
    voice: {engine: edge, name: en-US-AndrewNeural, rate: "-4%"}
    output: E:/silver_to_color/ep01/renders
    facts: FACTS.md
    full:
      cues:
        - id: open01
          say: "Four days before the city burned, a camera rolled down Market Street."
          at: 120          # film seconds where this cue's footage starts
          mode: wipe       # color | bw | split | wipe | card
          chapter: Cold open
          facts: [F02]
    shorts:
      - id: s01
        title: Filmed four days before it burned
        layout: stack      # stack | crop
        cues: [...]

Every field a cue can carry is in `Cue` below. Paths are relative to the YAML file.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import yaml

MODES_FULL = {"color", "bw", "split", "wipe", "card"}
MODES_SHORT = {"color", "bw", "wipe", "card"}
LAYOUTS = {"stack", "crop"}
CUE_ID = re.compile(r"^[A-Za-z0-9_\-]+$")


class SpecError(ValueError):
    pass


@dataclass
class Cue:
    id: str
    say: str = ""                 # narration; empty = no voice (music/ambience beat)
    at: float | None = None       # film time (s) where footage starts; None = carry on
                                  # from where the previous cue's footage ended
    mode: str = "color"
    speed: float = 1.0            # footage playback speed (0.5 = slow motion)
    hold: float = 0.0             # minimum cue duration (s)
    lead: float = 0.25            # silence before the line (s)
    pad: float = 0.5              # silence after the line (s)
    push: float = 0.0             # slow zoom-in over the cue, e.g. 0.06 = 6 %
    pan: float = 0.5              # Shorts crop layout: horizontal position 0..1
    chapter: str | None = None
    lower_third: str | None = None
    card_title: str | None = None
    card_subtitle: str | None = None
    facts: list[str] = field(default_factory=list)


@dataclass
class Short:
    id: str
    title: str
    hook: str = ""
    layout: str = "stack"
    cues: list[Cue] = field(default_factory=list)


@dataclass
class Episode:
    path: Path
    id: str
    title: str
    bw: Path
    color: Path
    credit: str
    output: Path
    voice: dict
    recordings: Path
    facts_path: Path | None
    cues: list[Cue]
    shorts: list[Short]
    description: str = ""
    thumbnail: dict = field(default_factory=dict)
    music: Path | None = None

    @property
    def root(self) -> Path:
        return self.path.parent


def _resolve(base: Path, value) -> Path:
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else (base / p)


def _cue(raw: dict, where: str, modes: set[str]) -> Cue:
    if not isinstance(raw, dict):
        raise SpecError(f"{where}: cue must be a mapping, got {type(raw).__name__}")
    known = set(Cue.__dataclass_fields__) | {"card"}
    unknown = set(raw) - known
    if unknown:
        raise SpecError(f"{where}: unknown cue keys {sorted(unknown)}")
    data = dict(raw)
    card = data.pop("card", None) or {}
    if card:
        data["card_title"] = card.get("title")
        data["card_subtitle"] = card.get("subtitle")
        data.setdefault("mode", "card")
    data["facts"] = [str(f) for f in (data.get("facts") or [])]
    cue = Cue(**data)
    cue.id = str(cue.id)
    if not CUE_ID.match(cue.id):
        raise SpecError(f"{where}: cue id {cue.id!r} must be letters, digits, _ or -")
    if cue.mode not in modes:
        raise SpecError(f"{where} ({cue.id}): mode {cue.mode!r} not in {sorted(modes)}")
    if not 0.1 <= cue.speed <= 4:
        raise SpecError(f"{where} ({cue.id}): speed {cue.speed} outside [0.1, 4]")
    if cue.at is not None and cue.at < 0:
        raise SpecError(f"{where} ({cue.id}): at must be >= 0")
    if not 0 <= cue.pan <= 1:
        raise SpecError(f"{where} ({cue.id}): pan must be in [0, 1]")
    if cue.mode == "card" and not cue.card_title:
        raise SpecError(f"{where} ({cue.id}): a card needs card: {{title: ...}}")
    if not cue.say.strip() and cue.hold <= 0:
        raise SpecError(f"{where} ({cue.id}): a cue with no narration needs hold > 0")
    cue.say = " ".join(cue.say.split())
    return cue


def load(path) -> Episode:
    path = Path(path).resolve()
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise SpecError(f"{path}: top level must be a mapping")
    base = path.parent
    for key in ("id", "title", "film", "output"):
        if key not in raw:
            raise SpecError(f"{path}: missing required key {key!r}")
    film = raw["film"]
    voice = {"engine": "edge", "name": "en-US-AndrewNeural", "rate": "+0%",
             **(raw.get("voice") or {})}
    full = raw.get("full") or {}
    cues = [_cue(c, f"full.cues[{i}]", MODES_FULL) for i, c in enumerate(full.get("cues") or [])]
    shorts = []
    for j, s in enumerate(raw.get("shorts") or []):
        layout = s.get("layout", "stack")
        if layout not in LAYOUTS:
            raise SpecError(f"shorts[{j}]: layout {layout!r} not in {sorted(LAYOUTS)}")
        shorts.append(Short(
            id=str(s["id"]), title=str(s.get("title", s["id"])), hook=str(s.get("hook", "")),
            layout=layout,
            cues=[_cue(c, f"shorts[{j}].cues[{i}]", MODES_SHORT)
                  for i, c in enumerate(s.get("cues") or [])]))
    ep = Episode(
        path=path, id=str(raw["id"]), title=str(raw["title"]),
        bw=_resolve(base, film["bw"]), color=_resolve(base, film["color"]),
        credit=str(film.get("credit", "")), output=_resolve(base, raw["output"]),
        voice=voice, recordings=_resolve(base, voice.get("recordings", "recordings")),
        facts_path=_resolve(base, raw["facts"]) if raw.get("facts") else None,
        cues=cues, shorts=shorts, description=str(raw.get("description", "")),
        thumbnail=raw.get("thumbnail") or {},
        music=_resolve(base, raw["music"]) if raw.get("music") else None)
    ids = [c.id for c in ep.cues] + [c.id for s in ep.shorts for c in s.cues]
    dupes = sorted({i for i in ids if ids.count(i) > 1})
    if dupes:
        raise SpecError(f"duplicate cue ids {dupes}: each id names its own recording file")
    return ep


FACT_ID = re.compile(r"^\s*(?:#+\s*|\*\*|[-*]\s*\**)?\**\[?(F\d{2,3})\]?", re.M)
NEEDS_SOURCE = re.compile(r"\d|\b(?:first|only|never|every|always|largest|oldest)\b", re.I)


def lint(ep: Episode, film_duration: float | None = None) -> list[str]:
    """Problems that would not crash a render but would mislead a viewer or YouTube.

    - a cue citing a fact ID that is not in the fact file
    - a cue stating a number or a superlative with no fact cited
    - footage requested past the end of the film
    - YouTube chapter rules: first at 0:00 (automatic), >= 3 chapters
    """
    out = []
    known: set[str] = set()
    if ep.facts_path and ep.facts_path.exists():
        known = set(re.findall(r"\b(F\d{2,3})\b", ep.facts_path.read_text(encoding="utf-8")))
    elif ep.facts_path:
        out.append(f"fact file {ep.facts_path} does not exist")
    groups = [("full", ep.cues)] + [(f"short {s.id}", s.cues) for s in ep.shorts]
    for name, cues in groups:
        for c in cues:
            for f in c.facts:
                if known and f not in known:
                    out.append(f"{name}/{c.id}: cites {f}, which is not in the fact file")
            if c.say and NEEDS_SOURCE.search(c.say) and not c.facts:
                out.append(f"{name}/{c.id}: states a number or superlative with no fact cited")
            if film_duration and c.mode != "card" and c.at is not None and c.at >= film_duration:
                out.append(f"{name}/{c.id}: at={c.at} is past the film's end ({film_duration:.1f}s)")
    chapters = [c for c in ep.cues if c.chapter]
    if ep.cues and len(chapters) < 3:
        out.append(f"full: {len(chapters)} chapter(s); YouTube needs at least 3 to show them")
    return out
