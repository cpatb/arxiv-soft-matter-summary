"""Timeline, captions, chapters and spec validation. No ffmpeg, no network."""

import textwrap

import pytest

from studio import episode as E
from studio import render as R
from studio import timeline as T
from studio.episode import Cue


def words(text, start=0.0, step=0.4):
    return [T.Word(w, start + i * step, start + i * step + 0.3) for i, w in enumerate(text.split())]


def test_layout_durations_and_footage_continuation():
    cues = [Cue(id="a", say="x", at=100), Cue(id="b", say="y"),
            Cue(id="card", mode="card", card_title="T", at=5, hold=3), Cue(id="c", say="z")]
    voice = {"a": (2.0, []), "b": (1.0, []), "c": (4.0, [])}
    p = T.layout(cues, voice)
    # duration = lead 0.25 + line + pad 0.5
    assert [round(x.duration, 3) for x in p] == [2.75, 1.75, 3.0, 4.75]
    assert [round(x.start, 3) for x in p] == [0.0, 2.75, 4.5, 7.5]
    # b carries on from a's footage; the card does not advance the ride, so c picks up
    # where b ended, not at the card's background time.
    assert p[1].film_at == pytest.approx(102.75)
    assert p[3].film_at == pytest.approx(102.75 + 1.75)
    assert T.total(p) == pytest.approx(12.25)


def test_layout_speed_scales_footage_consumed():
    cues = [Cue(id="a", say="x", at=10, speed=0.5), Cue(id="b", say="y")]
    p = T.layout(cues, {"a": (3.25, []), "b": (1.0, [])})
    assert p[1].film_at == pytest.approx(10 + 4.0 * 0.5)


def test_hold_wins_over_short_line_and_silent_cue_needs_hold():
    p = T.layout([Cue(id="a", say="x", hold=6)], {"a": (1.0, [])})
    assert p[0].duration == 6
    with pytest.raises(ValueError):
        T.layout([Cue(id="a", say="")], {})


def test_words_are_shifted_to_timeline_time():
    p = T.layout([Cue(id="a", say="one two", at=0), Cue(id="b", say="three four")],
                 {"a": (1.0, words("one two")), "b": (1.0, words("three four"))})
    assert p[1].words[0].start == pytest.approx(p[1].start + 0.25)


def test_estimate_words_partitions_duration():
    w = T.estimate_words("a bb ccc", 3.0)
    assert w[0].start == 0 and w[-1].end == pytest.approx(3.0)
    assert all(a.end == pytest.approx(b.start) for a, b in zip(w, w[1:]))
    assert (w[2].end - w[2].start) > (w[0].end - w[0].start)   # longer word, longer time
    assert T.estimate_words("", 3.0) == []


def test_srt_format_and_chunk_limits():
    text = ("For most of the last century this film had the wrong date and nobody "
            "noticed the puddles between the tracks until much later.")
    p = T.layout([Cue(id="a", say=text, at=0)], {"a": (10.0, words(text))})
    srt = T.to_srt(p)
    blocks = [b for b in srt.strip().split("\n\n")]
    assert blocks[0].startswith("1\n00:00:00,250 --> ")
    for b in blocks:
        line = b.split("\n")[2]
        assert len(line) <= T.FULL_MAX_CHARS
    assert T.srt_time(3723.4567) == "01:02:03,457"


def test_short_caption_chunks_respect_word_limit():
    p = T.layout([Cue(id="a", say="x", at=0)], {"a": (3.0, words("one two three four five six seven"))})
    for chunk in T.caption_chunks(p, short=True):
        assert len(chunk) <= T.SHORT_MAX_WORDS


def test_short_caption_events_never_overlap():
    """The 'RESTORATIONF,' bug: a chunk's last-word event must end before the next chunk."""
    text = "Top is the original print. Bottom is the restoration, colorized on a laptop."
    w = [T.Word(t, i * 0.3, i * 0.3 + 0.5) for i, t in enumerate(text.split())]   # overlapping words
    p = T.layout([Cue(id="a", say=text, at=0)], {"a": (5.0, w)})
    ev = T.ass_events(p, 1080, 1920, T.Theme(), short=True)
    spans = []
    for line in ev.splitlines():
        if ",Caption," in line:
            f = line.split(",")
            spans.append((f[1], f[2]))

    def sec(s):
        h, m, rest = s.split(":")
        return int(h) * 3600 + int(m) * 60 + float(rest)
    spans = sorted((sec(a), sec(b)) for a, b in spans)
    assert all(b <= c + 1e-9 for (_, b), (c, _) in zip(spans, spans[1:]))


def test_chapters_first_forced_to_zero_and_rules():
    cues = [Cue(id="a", say="x", at=0, chapter="One"), Cue(id="b", say="y", chapter="Two"),
            Cue(id="c", say="z", chapter="Three")]
    p = T.layout(cues, {"a": (20, []), "b": (3, []), "c": (20, [])})
    chs = T.chapters(p)
    assert chs[0][0] == 0.0
    probs = T.chapter_problems(chs, T.total(p))
    assert any("Two" in x for x in probs)          # 3.75 s chapter flagged
    assert T.chapter_problems(chs[:2], 100)[0].startswith("2 chapters")
    assert T.yt_time(3725) == "1:02:05" and T.yt_time(65) == "1:05"


def test_ass_colour_and_escape():
    assert T.ass_colour("#E8C170") == "&H0070C1E8"
    assert T.ass_escape("a{b}\\c\nd") == "a(b)\\\\c\\Nd"


def _write(tmp_path, body):
    p = tmp_path / "ep.yaml"
    p.write_text(textwrap.dedent(body), encoding="utf-8")
    return p


BASE = """
id: t
title: T
film: {bw: bw.mp4, color: c.mp4}
output: out
facts: FACTS.md
"""


def test_spec_rejects_unknown_keys_and_bad_modes(tmp_path):
    with pytest.raises(E.SpecError, match="unknown cue keys"):
        E.load(_write(tmp_path, BASE + "full:\n  cues:\n    - {id: a, say: x, colour: red}\n"))
    with pytest.raises(E.SpecError, match="mode"):
        E.load(_write(tmp_path, BASE + "full:\n  cues:\n    - {id: a, say: x, mode: sepia}\n"))
    with pytest.raises(E.SpecError, match="duplicate cue ids"):
        E.load(_write(tmp_path, BASE + "full:\n  cues:\n    - {id: a, say: x}\n    - {id: a, say: y}\n"))
    with pytest.raises(E.SpecError, match="hold"):
        E.load(_write(tmp_path, BASE + "full:\n  cues:\n    - {id: a, say: ''}\n"))


def test_lint_flags_unsourced_numbers_and_unknown_facts(tmp_path):
    (tmp_path / "FACTS.md").write_text("**F01: thing.** VERIFIED\n", encoding="utf-8")
    ep = E.load(_write(tmp_path, BASE + """
full:
  cues:
    - {id: a, say: "It was 1906.", chapter: A}
    - {id: b, say: "It was 1906.", facts: [F01], chapter: B}
    - {id: c, say: "Fine words.", facts: [F99], chapter: C}
    - {id: d, say: "The largest city.", at: 999}
"""))
    probs = E.lint(ep, film_duration=700)
    assert any("full/a" in x and "no fact cited" in x for x in probs)
    assert not any("full/b" in x for x in probs)
    assert any("F99" in x for x in probs)
    assert any("full/d" in x and "superlative" in x for x in probs)
    assert any("full/d" in x and "past the film" in x for x in probs)


def test_fact_sources_resolves_bibliography(tmp_path):
    f = tmp_path / "FACTS.md"
    f.write_text(textwrap.dedent("""
        **F01: Makers.** VERIFIED
        - S1 Kiehn [exact]: "..."
        - direct https://example.org/a.
        **F02: Other.** VERIFIED
        - S2 [exact]
        ## Bibliography
        | ID | Source | URL |
        |---|---|---|
        | S1 | Kiehn essay | https://loc.gov/k.pdf |
        | S2 | Catalog | https://loc.gov/item/1 |
        """), encoding="utf-8")
    src = R.fact_sources(f)
    assert src["F01"] == ["https://loc.gov/k.pdf", "https://example.org/a"]   # order of appearance
    assert src["F02"] == ["https://loc.gov/item/1"]
