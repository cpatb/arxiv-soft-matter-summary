# Workflow: pick a film → published video

The end-to-end map, where the time actually goes, and what has been automated.

---

## 1. The current shape, and where it hurts

| # | Stage | Mostly | Bottleneck | Automated? |
|---|---|---|---|---|
| 1 | Pick a film | browsing archive.org by hand | **No way to see scan resolution without opening each item** | ✅ `tools/scout_films.py` |
| 2 | Verify rights | ad-hoc | No record kept, so a Content ID claim means reconstructing months later | ✅ `tools/make_metadata.py` |
| 3 | Fetch + preprocess | ffmpeg by hand | Source is 85 GB; assumptions about fps unverified | ✅ notebook Stage 0-2 |
| 4 | GPU render | notebook | **Session limits, weekly quota, no resume** | ✅ notebook manifest |
| 5 | Grade + titles | DaVinci Resolve | Genuinely manual, and should stay that way | ❌ by design |
| 6 | Upload + metadata | typing into Studio | Retyping the same disclosure and credits every time | ✅ `tools/make_metadata.py` |
| 7 | Shorts cuts | manual | Manual, but see §5 | ⚠ partial |
| 8 | Tracking | nothing | Nothing recorded, so nothing learned | ⚠ §6 |

**The three bottlenecks that actually cost time**, in order:

1. **Finding a film worth restoring.** Item size is a bad proxy for scan quality: a 40 GB item
   may be a long, low-resolution film. Resolution lives in the per-item file manifest, so
   comparing candidates by hand means opening dozens of pages.
2. **GPU work not surviving a session.** A 12-hour cap, a 60-minute idle timeout and ~30 GPU-hr
   per week mean any film is a multi-session job. Without chunk-level resume, an interrupted
   session is lost work.
3. **Metadata retyping.** Same disclosure, same credit structure, same checklist every video,
   and the provenance record, the one thing that matters in a dispute, is the piece most likely
   to be skipped because it has no immediate payoff.

Stages 5 and the taste parts of 7 should stay manual. Grading is where the creative judgement
lives, and automating it would be automating the reason to do this at all.

---

## 2. `tools/scout_films.py`: find the films

Queries archive.org, reads each candidate's file manifest, ranks on **real pixel dimensions**,
and writes a CSV.

```bash
python scout_films.py --era 1900 1930 --min-size-gb 5 --sort "item_size desc" --deep
python scout_films.py --collection silent_films --deep --out silent.csv
python scout_films.py --topic "san francisco" --era 1900 1930 --deep
python scout_films.py --selftest        # offline
```

Verified live on 2026-08-26 with a real run against 983 matches:

```
 score   resolution      size  year  title
  74.5    4096x3112    98.9GB  1929  "The Desert Song" 1929 Restored Complete From 35mm N
        flags: too-big-for-free-tier;duration-mismatch(check-for-padded-frames)
  69.1    3840x2160     0.6GB  1930  Mickey Mouse Season 2 (1930-1932)
  61.4    2880x2160    42.3GB  1929  The Great Gabbo ( 1929)
  48.1    1920x1080   124.6GB  1906  A Trip Down Market Street Before the Fire (uncompres
```

Nineteen CSV columns, so you can re-sort on whatever you actually care about rather than
trusting the score.

### Ranking on bitrate, not just pixels

Resolution alone promotes junk, because anyone can upscale a low-resolution file to 4K. The
discriminating metric is **bytes per pixel-second**, `size / (width * height * length)`,
computed from the master's size but the *derivative's* length (a master's own length can be an
fps artifact). Observed values:

| bpps | What it means |
|---|---|
| ~80 | uncompressed 1080p |
| ~10 | genuine 4K film scan (Market Street 10.16, Code_Of_The_Sea 10.62) |
| 0.6-1.5 | compressed re-encode |
| ~0.14 or below | modern AI upscale of a low-resolution source |

The upscale case is the one worth catching: restoring someone else's upscale means inheriting
their artifacts and adding your own on top.

Adding this changed the live ranking substantially. "Mickey Mouse Season 2", 3840x2160, fell
from **2nd place to 6th** and now carries `LIKELY-UPSCALE-NOT-SCAN(bpps=0.01)`, while the
uncompressed 1080p Market Street scan at 80.01 bpps rose from 8th to 3rd. **Five of the top
twelve results were resolution-only artifacts.**

### Three API traps it defends against

**It does not use the Scrape API, and that is deliberate.** `services/search/v1/scrape` has a
live, reproducible **server-side cache whose key omits the query string** for `count` ≤ 250.
Distinct queries return byte-identical results, and during testing a response belonging to
*another archive.org user's* query leaked in: a 5.2-million-hit result set arriving for a
`collection:prelinger` query. `advancedsearch.php` never showed this and is used throughout. If
you ever need Scrape (collections over 10,000 items), pass `count ≥ 251` and verify the payload
against the query.

**`source: "original"` does not mean "not derived".** The Market Street HEVC file is labelled
`source: "original"` while carrying an `original` pointer to the `.mov` it was made from. The
correct test needs both conditions:

```python
is_derived = f.get("source") == "derivative" or "original" in f
```

Testing `source` alone misclassifies the flagship item's own derivative as a master.

**`year` is uploader-supplied and unvalidated.** A Prelinger industrial short was observed dated
**2026**, and an item dated **1903** turned out to be 540 GB of *21 Jump Street* episodes. Use
`year` to shortlist; never treat it as evidence of publication date. The script prints this
reminder on every run.

### Collection identifiers, verified by query

Five commonly-cited collection names do not exist. Querying them returns zero results, which
looks like "no matching films" rather than "wrong name":

| Does not exist | Use instead | Items |
|---|---|---|
| `newsreels` | `universal_newsreels` | 611 |
| `sabucat` | `sabucat_trailers` | 103 |
| `usnationalarchives` | `us-national-archives-youtube` | 3,808 |
| `AVGeeks` | `av_geeks_video_2017` (exists but is **empty**) | 0 |
| `film_noir` | `Film_Noir` (capital F and N) | 1,673 |

Verified and working: `prelinger` (10,466), `feature_films` (28,470), `silent_films` (3,527),
`academic_films` (494), `home_movies` (1,898), `nasa` (13,733 movies of 208,851 items).
`opensource_movies` holds 2.6 million items with uploader-asserted rights, which is the lowest
trust tier on the site.

The correct names are in `--collection`'s help text so they are visible at the point of use.

### How it ranks

Resolution dominates, deliberately: a 12.6 MP master outranks a 3.1 MP scan even when it is too
big to stage locally, because it genuinely is the better source and URL segmentation handles the
size. A too-big-for-free-tier penalty exists but is a tiebreak, not a veto. Popularity is a weak
term, because the goal is finding good scans, not re-restoring the ten most-uploaded films on the site.
An obscure 4K scan outranks a famous 640x480 one, and there is a test asserting exactly that.

`duration-mismatch(check-for-padded-frames)` in the flags column is worth acting on: it means
master and derivative disagree on runtime, which usually means padded/duplicated frames and
therefore free GPU savings (see the pipeline doc §5.6).

---

## 3. `tools/make_metadata.py`: rights record and description

```bash
python make_metadata.py MarketStreet19064KScan20181016 --outdir ../metadata \
    --title "A Trip Down Market Street, 1906 - Restored and Colorized in 4K"
```

Writes three files per film:

- `<slug>.description.txt`: paste-ready description plus a tag block that fits YouTube's
  500-character budget
- `<slug>.provenance.json`: the record you cite in a dispute
- `<slug>.checklist.md`: the per-video checklist from §4, pre-filled

**The metadata endpoint has three states, not two**, and all of them return HTTP 200: `{}` means
the item does not exist, **`is_dark: true` means it was withdrawn**, typically after a copyright
claim, and a `metadata` block means it is live. Collapsing "dark" into "missing" would discard
the single most useful rights signal archive.org emits: a dark item is one somebody successfully
claimed against. `ia_meta.fetch_item()` raises on it, and `ia_meta.item_status()` reports it
without raising.

**Why the provenance file earns its keep.** "It's in the public domain" is **not** one of
YouTube's enumerated Content ID dispute reasons. You file under *"misidentified or an error was
made"*, which means naming the specific archive item, its rights metadata, and the film's year.
Writing that down at publish time costs seconds; reconstructing it under a 30-day dispute clock
does not.

The record deliberately stores the **absence** of rights metadata as an explicit null plus a
note, because a missing key reads as "not checked", which is exactly the ambiguity it exists to remove.
For this item, archive.org carries **no `licenseurl`, no `rights`, no `possible-copyright-status`
at all**. Absence is not evidence either way; the public-domain position rests on the 1906
publication date alone.

The generated description also carries the uploader's own reuse grant, which is worth quoting
because it is unusually explicit:

> "Anyone may reproduce or reuse this scan. Please attribute it to its source: Prelinger
> Archives."

**A film that is not automatically public domain produces a warning instead of a rights claim.**
Passing a 1999 film emits `!! VERIFY BEFORE PUBLISHING` rather than asserting PD status, and so
does a film with no readable year. There is a test for both, because getting this backwards is
the one failure mode with legal consequences.

### Rights checking: three corrections worth having

**The renewal band starts at the PD cutoff and moves every year.** In calendar 2026 that is
**1931** to 1963, not 1929 to 1963. The cutoff advances each January.

**Renewals are in the ANNUAL Catalog of Copyright Entries, Parts 12-13.** The *cumulative*
volumes contain no renewal entries at all. Searching a cumulative volume and finding nothing is
the usual reason a film looks unrenewed when it is not.

**"Probably never renewed" is a bad prior for film.** Films had the **highest renewal rate of any
class of work**, roughly 74% and 65% in the studied samples, against about 7% for books. The
heuristic that works for books inverts for movies.

Also worth knowing: archive.org staff maintain **OKAY / NOT-OKAY / Unknown lists with renewal
registration numbers** at `archive.org/details/movie_research`. Films on the NOT-OKAY list are
live on the site today under Creative Commons Public Domain Marks, two of them uploaded in 2025.
Archive.org states it "does not make guarantees as to the copyright status of items" and that it
has no copyright experts. **A Public Domain Mark on an archive.org item is an uploader's
assertion, not a determination.**

---

## 4. The upload checklist

Generated per film; `tools/upload_checklist.md` is the master copy.

The two items under "the two that protect the channel" are the ones to never skip. Everything
else costs you a video; those two cost you the channel.

---

## 5. Shorts cuts

The pipeline renders `compare_sidebyside.mp4` and `compare_wipe.mp4` automatically (notebook
Stage 5), which is the raw material. What it deliberately does not do is auto-cut them to
length, because **the 60-second minimum on TikTok changes what the clip is**, not just how long
it is.

A 15-second silent before/after reveal fails on three counts at once: it is under TikTok's
minimum, it is the "filter applied to someone else's footage" shape that every platform's
originality rule targets, and it earns roughly nothing. The same restoration with 60 seconds of
narration about what the footage shows clears all three. That is a scripting decision, not a
cropping one.

---

## 6. Tracking

Not built, and worth being honest about why: with nothing published there is no data, and a
tracking system designed against imagined data usually measures the wrong thing.

What is worth recording from video one, in a plain spreadsheet:

| Column | Why |
|---|---|
| Archive identifier | Ties the video back to its provenance record |
| bytes-per-pixel-second | Whether you restored a scan or someone else's upscale |
| Narrated? Y/N | The one variable with a **3-5x RPM effect** (monetization doc §3.2) |
| Length, mid-roll count | The second-largest RPM lever |
| RPM, and audience-geography split | RPM is meaningless without the geography that produced it |
| Content ID claims | Which archives and which scans attract them |
| GPU-hours spent | Against the 30 h/week budget, to see which films are affordable |

The narrated/unnarrated column is the important one. It is the only variable here worth a
deliberate A/B, and the effect size is large enough to be visible after a handful of videos.

---

## 7. What a film costs, end to end

Rough, and dominated by choices rather than by the tools:

| Stage | Time | Notes |
|---|---|---|
| Scout + verify rights | 15-30 min | Mostly the manual renewal check, which is the part that cannot be automated |
| Preprocess (CPU) | 20-60 min | Unattended |
| Colorize + cleanup (GPU) | **4-10 GPU-hr** | The binding constraint. Decimation removes a third of it |
| Grade + titles | 1-3 h | Manual, and the part worth doing well |
| **Script + record narration** | **2-4 h** | The largest manual block, and the highest-value one |
| Metadata + upload | 10 min | Automated |
| Shorts cuts | 30-60 min | Per clip |

**Narration is the biggest time cost and the biggest return.** It is simultaneously the
reused-content mitigation, the generic-content mitigation, TikTok's 60-second requirement,
Meta's originality carve-out, and a 3-5x RPM multiplier. Five pressures, one piece of work.

A realistic cadence on a free tier is **one film every two to three weeks**, which against the
2027-02-01 YPP deadline (monetization doc §2.2) is the number that should drive planning.

---

## 8. Tools reference

| File | What it does | Verification |
|---|---|---|
| `ia_meta.py` | archive.org metadata access, file ranking | offline unit tests + live check |
| `scout_films.py` | Find and rank candidate films → CSV | `--selftest`, plus a live 983-match run |
| `make_metadata.py` | Description, tags, provenance, checklist | `--selftest`, plus a live run |
| `pipeline.py` | Chunking, manifest, deflicker, decimate, assemble | `test_pipeline.py`, **17/17** |
| `temporal_chroma.py` | Motion-compensated chroma median | `test_temporal_chroma.py`, **6/6** |
| `build_notebook.py` | Generates the notebook; syntax-checks every cell | 15/15 code cells parse |

Run everything:

```bash
cd tools
python test_pipeline.py && python test_temporal_chroma.py
python scout_films.py --selftest && python make_metadata.py --selftest
```

Standard library plus `numpy`, `opencv-python` and `ffmpeg` for the pipeline modules. The
archive.org tools are **stdlib only**, no `internetarchive` install needed, since the JSON
endpoints are public and unauthenticated.

Three bugs were found by these tests rather than in production, which is the argument for having
them: ffmpeg's filter parser rejecting Windows absolute paths (its `:` option separator), a
missing container-level duration fallback for Matroska, and a runt final chunk that decoded as
garbage. All three would have surfaced mid-render otherwise.
