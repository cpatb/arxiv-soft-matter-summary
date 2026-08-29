# Silver to Color: the ten highest-leverage actions

Ordered by impact x ease. Each links to the reasoning.

| # | Action | Where |
|---|---|---|
| 1 | **Narrate every video.** One decision satisfies four separate pressures at once (YouTube's reused-content test, its generic-content test, TikTok's 60-second minimum, and Meta's originality carve-out) and is worth **3-5x on RPM** ($7-9 narrated vs $1.50-4 unnarrated, because silent footage is language-neutral and pulls low-RPM global traffic). | [monetization §3.2](research/monetization.md) |
| 2 | **Get *accepted* into YPP before 2027-02-01**, not merely applied. Entry thresholds double for new applicants that day and existing members are grandfathered, so with a ~1-month review queue the real deadline is around end of December 2026. | [monetization §2.2](research/monetization.md) |
| 3 | **Stand up the commissioned-restoration service now**, before the channel exists: it is the only stream that pays at zero subscribers, the only one the AI-copyright problem cannot touch, and it is priced into a vacuum rather than against an incumbent. | [monetization §5.1](research/monetization.md) |
| 4 | **Decimate duplicate frames before the GPU stages.** The 4K master pads 16 fps to 24 fps by doubling every second frame, so one frame in three is an exact duplicate; **measured at 33.3% of all GPU work saved**, for one ffmpeg pass. | [pipeline §5.6](research/pipeline-optimization.md) |
| 4b | **Rank candidate scans by bytes-per-pixel-second, not resolution.** Anyone can upscale a file to 4K; the bitrate behind the pixels is what separates a genuine scan (~10 bpps) from a modern AI upscale (~0.14). Adding it moved a 3840x2160 upscale from 2nd place to 6th and promoted a real 1080p scan from 8th to 3rd. | [workflow §2](research/workflow.md) |
| 5 | **Use the 2 GB cropped 2K source for the first film**, not the 85.4 GB master: it is already de-sprocketed and contrast-adjusted, the uploader recommends it, and it makes the 85 GB-file-on-a-20 GB-disk problem disappear entirely. | [pipeline §5.6](research/pipeline-optimization.md) |
| 6 | **Stop upscaling and deliver 2160p as a cleanup pass.** A 4K scan of 1906 stock holds no sub-grain detail to recover, and the reason to deliver 4K is YouTube's bitrate ladder (35-45 Mbps vs 8 Mbps at 1080p), not resolution. | [pipeline §0.2](research/pipeline-optimization.md) |
| 7 | **Change the chroma upsample from nearest to bilinear.** It is a one-word fix for a whole class of chroma jitter that DDColor's shipped `F.interpolate` call creates and that would otherwise be measured as, and blamed on, the model. | [pipeline §2.1](research/pipeline-optimization.md) |
| 8 | **Deflicker with `mode=am`, not `mode=median`.** Measured across three flicker profiles and eight seeds, `am`/`gm` cut flicker 78-88% while `median` managed −0.5% to 6.5%, so the plausible outlier-rejection argument is simply wrong here. | [pipeline §2.3](research/pipeline-optimization.md) |
| 9 | **Motion-compensate the temporal chroma median.** Warping neighbours by optical flow before taking the median gives **3.4x lower chroma error at moving edges (measured)**, and it is nearly free because the flow only needs computing at DDColor's native 512x512. | [pipeline §2.4](research/pipeline-optimization.md) |
| 10 | **Disclose AI colorization in Studio and keep a provenance record per film.** Disclosure is stated to cost nothing in reach or eligibility while non-disclosure risks YPP suspension, and "public domain" is *not* an enumerated Content ID dispute reason, so the record is what you actually file under. | [monetization §1.5, §1.6](research/monetization.md) |

---

## What is in this folder

```
research/monetization.md            YPP policy, RPM evidence, six revenue streams, staged roadmap
research/pipeline-optimization.md   Ranked changes with gain / cost / risk, and what to avoid
research/workflow.md                End-to-end map, bottlenecks, tool reference
market_street_pipeline.ipynb        Reference notebook: chunked, resumable, 25 cells
tools/                              Scouting, metadata, pipeline modules + their tests
candidates_demo.csv                 Live scout run, 12 real candidates, 20 columns
metadata/                           Generated description + provenance + checklist for film 1
```

**On the notebook:** the file named in the brief was not on this machine. A forced recursive
scan of the whole user profile found no match, so `market_street_pipeline.ipynb` here is a new
reference implementation rather than an edit of yours. Diff it against your Kaggle copy and take
what you want. Every code cell is syntax-checked; the non-GPU logic it orchestrates lives in
`tools/` and is covered by tests.

## Verification

Nothing here is asserted from a blog post. Model VRAM and parameter counts come from repository
source, releases and issue threads; the empirical claims were measured on this machine.

```
tools/test_pipeline.py           17/17   real ffmpeg, real files, simulated session kill
tools/test_temporal_chroma.py     6/6    motion + mode flips + occlusion fixture
tools/scout_films.py --selftest   PASS   plus a live run against 983 archive.org matches
tools/make_metadata.py --selftest PASS   plus a live run on the real item
market_street_pipeline.ipynb     15/15   code cells parse
```

Three bugs surfaced in testing rather than mid-render: ffmpeg's filter parser rejecting Windows
absolute paths (its `:` option separator), a missing container-level duration fallback for
Matroska, and a runt final chunk that decoded as garbage.

## Four things that are genuinely undecided

Flagged rather than smoothed over, because each would change what you do:

1. **Whether an AI restoration pipeline counts as "substantive modifications" under YouTube's
   reused-content policy.** No first-hand account from anyone applying to YPP with AI-colorized
   archival film could be found, and that absence *is* the finding. The strongest evidence either
   way is NASS: a pure AI restoration channel, started 2020, 787K subscribers, running ads, not
   grandfathered. Against that, a musician was ejected from YPP partly over Prelinger
   public-domain footage, which is the same archive as your first film.
2. **Which direction Real-ESRGAN's `-dn` knob runs.** Two independent analyses disagreed; the
   README's wording and an empirical probe of both checkpoints point opposite ways. Both end at
   the same instruction (A/B at 1.0 / 0.5 / 0.2 / 0.0), but label the renders by flag value, not
   by "more grain".
3. **Whether AI colorization requires YouTube's synthetic-content disclosure.** Upscaling is
   exempt by name; colorization is unaddressed, sitting between an exempt bullet and a triggering
   one. The sources are silent, not contradictory. Disclose anyway; the cost is zero.
4. **RIFE frame interpolation.** A taste decision, not a technical one. Silent film was
   hand-cranked at 12-40 fps, so there is no single "true" rate to restore.

## One thing that is not undecided, and is worth stating plainly

**A Public Domain Mark on an archive.org item is an uploader's assertion, not a determination.**
Archive.org says so itself: it "does not make guarantees as to the copyright status of items"
and has no copyright experts. Its own staff maintain a NOT-OKAY-to-post list with renewal
registration numbers at `archive.org/details/movie_research`, and films on that list are live on
the site today under Public Domain Marks, two uploaded in 2025.

Three corrections follow, all in `research/workflow.md`: the renewal band is **1931**-1963 in
calendar 2026 and moves forward every January; renewals appear only in the **annual** Catalog of
Copyright Entries, Parts 12-13, never the cumulative volumes; and films had the **highest renewal
rate of any class of work**, roughly 74% against about 7% for books, so "probably never renewed"
is a heuristic that works for books and inverts for film.

## The through-line

The same answer keeps arriving from unrelated directions. The fix for the monetization-policy
risk, the fix for the short-form originality rules, the fix for the RPM problem, and the thing
that makes the videos worth watching are **all the same fix**: a narrated restoration with
historical context, cut down into short-form teasers, rather than a silent before/after with a music
bed.

That is convenient, and it is worth acting on before publishing anything, because enforcement is
channel-level and retroactive. A back catalogue of silent loops is a liability, not an asset.
