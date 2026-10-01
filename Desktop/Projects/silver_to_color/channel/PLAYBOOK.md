# Silver to Color: production playbook

How an episode goes from an archive scan to a published long video plus four Shorts,
using the tools in this folder. The research behind the format choices is in
`research/monetization.md` and `research/workflow.md`; this file is the how-to.

## The format

**Premise.** Old film, restored properly and explained. Each episode is one film: we ride
through it in order, at its true speed, and the narration tells you what you're looking at,
how the film was made, and what happened next. The restoration is the hook, but the story
is what people stay for.

**Full video (8-12 min, 16:9).** Eight minutes is the mid-roll threshold.

| Block | Length | Screen |
|---|---|---|
| Cold open: the striking moment + the stakes | 20-40 s | `wipe` from B&W into colour, then `split` |
| Title card | 4-5 s | `card` |
| The making of the film | 1-2 min | the ride, in order (`color`) |
| 2-4 story chapters | 1-2 min each | the ride continues; jump with `at:` only for callouts |
| What happened next | ~1 min | `bw` for gravity, `wipe` back to colour |
| About this restoration (always) | 40-60 s | `split`: the print beside the restoration |
| Sign-off + end card | 10 s | `wipe`, `card` |

Rules that come straight from the research:
- **Narrate everything.** Narration is the reused-content defence, the RPM multiplier and the
  reason anyone watches to the end. A silent before/after is the format to avoid.
- **Show the transformation in the first 10 seconds** (cold-open wipe). A reviewer deciding
  "reused content" sees seconds, not your effort.
- **Every factual line cites a fact ID** in the episode's `FACTS.md`; `lint` enforces it.
  Disputed facts are narrated with their source named ("according to Kiehn").
- **Say the colour is an estimate**, every episode, in the restoration chapter and in the
  description. This is also the honest thing to do.
- **No library music.** It halves Shorts revenue and is an originality red flag. Silence
  under narration is fine; a self-made bed can be added with `music:` in the spec (it is
  ducked under the voice automatically).

**Shorts (60-75 s, 9:16).** One surprising, sourced fact per Short.
- The hook is on screen from frame one (`hook:`), and the first spoken line restates it.
- `stack` layout (original over restored) for "look at this" Shorts; `crop` layout with
  `wipe` cues for story Shorts.
- At least 60 s, because TikTok's Creator Rewards needs it (`shorts` prints a warning below that).
- End by pointing to the full episode. Shorts are the funnel; the long video is where the
  money is (Shorts RPM is 40-300x lower).

## The pipeline, per episode

```bash
# 0. Pick the film: channel/SLATE.md, or scout a new one
python tools/scout_films.py --era 1900 1930 --deep
python tools/make_metadata.py <archive-id> --outdir metadata      # rights/provenance record

# 1. Research: episodes/<ep>/FACTS.md, every fact with a fetched source (see ep01 for the shape)

# 2. Restore (local, Intel iGPU; ~2-3 h for a 12-minute film, resumable if interrupted)
python -m studio.colorize prepare <source.mov> E:/silver_to_color/<ep>/bw_master.mp4
python -m studio.colorize color   E:/silver_to_color/<ep>/bw_master.mp4 E:/silver_to_color/<ep>/color_master.mp4

# 3. Write episodes/<ep>/episode.yaml (copy ep01's), then
python -m studio lint   episodes/<ep>/episode.yaml   # facts cited, footage exists, chapters valid
python -m studio script episodes/<ep>/episode.yaml   # teleprompter for recording

# 4. Render (scratch voice first to check the cut, then again after recording)
python -m studio all episodes/<ep>/episode.yaml
```

Outputs land in the spec's `output:` folder: `<id>_full.mp4`, `.srt`, `.description.txt`,
`<id>_thumb.jpg`, and `shorts/*.mp4`.

**Recording the narration.** The edge-tts voice is a scratch track for timing the cut, not
the channel's voice. Record each line as `episodes/<ep>/recordings/<cue id>.wav` (the
`script` command lists the IDs), any decent USB mic in a soft-furnished room, ~0.3 s of
silence at each end. Re-render and the timeline re-flows to your timing. Lines you have not
recorded yet still use the scratch voice, so you can record in batches. Final loudness is
normalised to -14 LUFS for you.

## Upload checklist (per video)

1. Title, description (`.description.txt`: chapters, credit, colour disclaimer, sources).
2. Upload the `.srt` as English subtitles. Captions are not burned into the long video.
3. Thumbnail (`_thumb.jpg`; re-pick the frame with `thumbnail: {at: ..., text: ...}`).
4. **Altered or synthetic content: Yes** (the colour is AI-generated). Disclosure costs
   nothing; non-disclosure risks the channel (monetization §1.5).
5. Keep the provenance record from `tools/make_metadata.py` with the episode. It is what you
   cite if a Content ID claim arrives (monetization §1.6).
6. Shorts: post the same day as the long video, then one every 2-3 days, each linking the
   full episode (YouTube's "related video" field).

## Cadence

YPP thresholds double for new applicants on 2027-02-01, and review takes about a month, so
the practical target is 1,000 subscribers and 4,000 watch hours by the end of December 2026.
One episode plus four Shorts every two weeks gives six episodes and 24 Shorts by then. The
restore step runs unattended overnight, so the binding cost is research and recording
(about a day per episode).

## Episode 1 status

`episodes/ep01_market_street/`: FACTS.md (51 sourced facts, 11 myths), episode.yaml (8.4-min
full cut in 9 chapters + 4 Shorts), restore.sh. Renders on E:\silver_to_color\ep01\renders.
Next films: `channel/SLATE.md` (top pick: the 1915 Panama-Pacific Exposition opening parade,
which is the same city nine years after the fire).
