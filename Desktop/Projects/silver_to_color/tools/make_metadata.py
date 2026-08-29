"""Generate YouTube metadata and a provenance record for one archive.org film.

Produces three files per film:

  <slug>.description.txt  - paste-ready YouTube description
  <slug>.provenance.json  - the record you cite if a Content ID claim lands
  <slug>.checklist.md     - per-video upload checklist, pre-filled

Why the provenance file matters: "it's in the public domain" is NOT one of YouTube's
enumerated Content ID dispute reasons. You dispute under "misidentified or an error was
made", which means naming the specific archive item, its stated rights status, and the
film's release year. Writing that down at publish time is far easier than reconstructing
it months later under a 30-day dispute clock.

Standard library only.

Usage:
    python make_metadata.py MarketStreet19064KScan20181016 --outdir ../metadata
    python make_metadata.py <id> --title "A Trip Down Market Street, 1906" --year 1906
    python make_metadata.py --selftest
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import textwrap
from datetime import date, timezone, datetime
from pathlib import Path

import ia_meta

# ----------------------------------------------------------------------------
# Editable defaults. Change these once; every generated description follows.
# ----------------------------------------------------------------------------

CHANNEL_NAME = "Silver to Color"

# The AI disclosure block. Two things are deliberate here:
#  - it names colorization as a reconstruction, not a record. That is honest, and it is
#    also the "educational value" argument that YouTube's reused-content policy rewards.
#  - it separates upscaling (exempt from YouTube's synthetic-content disclosure) from
#    colorization (unaddressed by YouTube, hence disclosed voluntarily).
AI_DISCLOSURE = """\
ABOUT THE COLOR
The original film is black and white. Every color you see here is a reconstruction
inferred by an AI model from the brightness of the original frames -- it is an informed
guess, not a record of what these colors actually were. Treat the color as illustration
and the shapes, motion and composition as the historical document. Where a color is
known from other sources it is noted in the video.

WHAT WAS DONE
  * Speed corrected from the scan rate to the film's intended projection rate
  * Brightness flicker reduced (the source is {source_year}-era nitrate stock)
  * Colorized per frame, then stabilized across frames so colors do not flicker
  * Upscaled and lightly denoised
  * Scored (the original is silent)

No frames were generated, removed, or reordered. Nothing was added to the scene."""

SOURCE_BLOCK = """\
SOURCE
{title}
Archive.org item: {details_url}
Collection: {collection}
File used: {filename} ({resolution}, {size})
{uploader_line}
Downloaded and processed {process_date}."""

PD_BLOCK = """\
RIGHTS
This film was published in {year} and is in the public domain in the United States
(published before {pd_cutoff}; see the Cornell copyright term chart at
https://copyright.cornell.edu/publicdomain). The restoration was made from the archive
scan linked above. The colorization is machine-generated and, per U.S. Copyright Office
guidance on automated restoration, is not claimed as a copyrightable contribution --
this video is offered without a copyright claim over the underlying footage or its
colors. Please credit the archive, not me, if you reuse the source."""

# Tags: YouTube allows 500 characters total across all tags.
BASE_TAGS = [
    "restored film", "colorized", "public domain film", "silent film",
    "film restoration", "history", "archival footage", "AI colorization",
    "4k restoration", "historical footage",
]

MAX_TAG_CHARS = 500
MAX_TITLE_CHARS = 100
MAX_DESCRIPTION_CHARS = 5000

# U.S. copyright: works published 95+ years ago are in the public domain.
# In calendar year Y, that means published before Y - 95.
PD_TERM_YEARS = 95


def pd_cutoff_year(today: date | None = None) -> int:
    """The first year that is NOT automatically public domain in the US.

    In 2026 this returns 1931: works published 1930 and earlier are PD.
    """
    today = today or date.today()
    return today.year - PD_TERM_YEARS


def slugify(text: str, max_length: int = 60) -> str:
    """Filesystem-safe lowercase slug."""
    if not text:
        return "untitled"
    slug = re.sub(r"[^\w\s-]", "", text.lower())
    slug = re.sub(r"[\s_-]+", "-", slug).strip("-")
    return (slug[:max_length].rstrip("-")) or "untitled"


def extract_year(metadata: dict) -> int | None:
    """Pull a four-digit year out of archive.org's date/year fields.

    archive.org is inconsistent here: 'year', 'date' and 'publicdate' may hold
    '1906', '1906-01-01', or a full timestamp. Returns None rather than guessing.
    """
    for field in ("year", "date"):
        raw = metadata.get(field)
        if raw is None:
            continue
        match = re.search(r"\b(1[6-9]\d{2}|20\d{2})\b", str(raw))
        if match:
            return int(match.group(1))
    return None


def build_tags(extra: list[str] | None = None, year: int | None = None) -> list[str]:
    """Assemble a tag list that fits YouTube's 500-character budget."""
    tags = list(BASE_TAGS)
    if year:
        tags.insert(0, str(year))
        tags.insert(1, f"{(year // 10) * 10}s")
    for tag in extra or []:
        if tag and tag not in tags:
            tags.append(tag)

    kept: list[str] = []
    used = 0
    for tag in tags:
        # +1 for the separating comma between tags.
        cost = len(tag) + (1 if kept else 0)
        if used + cost > MAX_TAG_CHARS:
            continue
        kept.append(tag)
        used += cost
    return kept


def build_description(
    item: dict,
    identifier: str,
    chosen_file: dict | None,
    title: str,
    year: int | None,
    extra_notes: str = "",
) -> str:
    """Assemble the full YouTube description."""
    metadata = item.get("metadata", {})
    collections = ia_meta.as_list(metadata.get("collection")) or ["(none listed)"]
    uploader = metadata.get("uploader")

    if chosen_file:
        resolution = (
            f"{chosen_file['width']}x{chosen_file['height']}"
            if chosen_file.get("width") and chosen_file.get("height")
            else "resolution not stated in metadata"
        )
        filename = chosen_file["name"]
        size = ia_meta.human_size(chosen_file.get("size"))
    else:
        resolution, filename, size = "unknown", "(no video file found)", "?"

    source = SOURCE_BLOCK.format(
        title=metadata.get("title", identifier),
        details_url=ia_meta.details_url(identifier),
        collection=", ".join(collections),
        filename=filename,
        resolution=resolution,
        size=size,
        uploader_line=f"Scanned/uploaded by: {uploader}" if uploader else
                      "Uploader not stated in item metadata.",
        process_date=date.today().isoformat(),
    )

    disclosure = AI_DISCLOSURE.format(source_year=year or "early-cinema")

    if year:
        cutoff = pd_cutoff_year()
        if year < cutoff:
            rights = PD_BLOCK.format(year=year, pd_cutoff=cutoff)
        else:
            rights = (
                f"RIGHTS\n"
                f"!! VERIFY BEFORE PUBLISHING. This film is dated {year}, which is NOT\n"
                f"automatically public domain in the US (the {PD_TERM_YEARS}-year cutoff for\n"
                f"{date.today().year} is {cutoff}). Check renewal records before you upload."
            )
    else:
        rights = (
            "RIGHTS\n"
            "!! VERIFY BEFORE PUBLISHING. No publication year could be read from the\n"
            "archive.org metadata, so public-domain status has not been established."
        )

    parts = [
        title,
        "",
        disclosure,
        "",
        source,
        "",
        rights,
    ]
    if extra_notes.strip():
        parts += ["", "NOTES", extra_notes.strip()]
    parts += [
        "",
        f"-- {CHANNEL_NAME}",
    ]

    description = "\n".join(parts)
    if len(description) > MAX_DESCRIPTION_CHARS:
        description = (
            description[: MAX_DESCRIPTION_CHARS - 40].rstrip() + "\n[description truncated]"
        )
    return description


def build_provenance(
    item: dict, identifier: str, chosen_file: dict | None, year: int | None
) -> dict:
    """The record to cite in a Content ID dispute.

    Deliberately records the ABSENCE of rights metadata as an explicit null plus a
    note, rather than omitting the key -- a missing key reads as "not checked", which
    is exactly the ambiguity this file exists to remove.
    """
    metadata = item.get("metadata", {})
    rights_fields = {
        field: metadata.get(field)
        for field in ("licenseurl", "rights", "possible-copyright-status")
    }
    has_rights_metadata = any(value is not None for value in rights_fields.values())

    return {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "archive_identifier": identifier,
        "archive_details_url": ia_meta.details_url(identifier),
        "archive_title": metadata.get("title"),
        "archive_collections": ia_meta.as_list(metadata.get("collection")),
        "archive_uploader": metadata.get("uploader"),
        "archive_publicdate": metadata.get("publicdate"),
        "archive_addeddate": metadata.get("addeddate"),
        "film_year": year,
        "us_pd_cutoff_year_at_generation": pd_cutoff_year(),
        "pd_basis": (
            f"Published {year}, before the {pd_cutoff_year()} cutoff "
            f"({PD_TERM_YEARS}-year term)."
            if year and year < pd_cutoff_year()
            else "NOT ESTABLISHED - verify renewal records before publishing."
        ),
        "archive_rights_metadata": rights_fields,
        "archive_rights_metadata_present": has_rights_metadata,
        "rights_metadata_note": (
            "archive.org populated at least one rights field; it is still not "
            "authoritative and should be corroborated."
            if has_rights_metadata
            else "archive.org carries NO rights metadata for this item. Absence of a "
                 "licenseurl is not evidence either way; PD status rests on the "
                 "publication date alone."
        ),
        "source_file": chosen_file,
        "source_file_url": (
            ia_meta.download_url(identifier, chosen_file["name"]) if chosen_file else None
        ),
        "audio_note": (
            "Original film is silent. Any audio in the published video is added by the "
            "channel and must be self-composed or provably clear -- an inherited modern "
            "score is a valid Content ID claim and also risks YouTube membership "
            "eligibility."
        ),
        "dispute_note": (
            "'Public domain' is NOT an enumerated YouTube dispute reason. File under "
            "'misidentified or an error was made', and cite: this item URL, the film "
            "year above, and the fact that PD footage is ineligible as a Content ID "
            "reference (support.google.com/youtube/answer/2605065)."
        ),
        "verification_still_required": [
            "Confirm the film's publication year from a second source, not archive.org alone.",
            "Confirm no separately-copyrighted score is embedded in the file you used.",
            "Confirm this scan is not itself a third party's claimed restoration.",
            "For works published from the cutoff year through 1963, check renewal in the "
            "ANNUAL Catalog of Copyright Entries, Parts 12-13. The CUMULATIVE volumes "
            "contain no renewal entries at all.",
            "Films had the HIGHEST renewal rate of any class (roughly 74%/65%, against "
            "about 7% for books), so 'probably never renewed' is a bad prior for film.",
            "Check archive.org/details/movie_research: staff maintain a NOT-OKAY-to-post "
            "list with renewal numbers, and films on it are live on the site anyway "
            "under Public Domain Marks.",
        ],
    }


CHECKLIST_TEMPLATE = """\
# Upload checklist - {title}

Archive item: {details_url}
Generated: {today}

## Before render
- [ ] Publication year confirmed from a second source (not archive.org alone)
- [ ] Item is not dark (withdrawn after a claim) - the tools raise on this
- [ ] Not on archive.org's own NOT-OKAY list at archive.org/details/movie_research
- [ ] Source file is the rawest available scan, not someone else's restoration
- [ ] Any embedded audio muted; score is self-composed or provably clear
- [ ] Shot boundaries detected (deflicker and temporal median must run per shot)

## Render
- [ ] Speed correction verified against the film's intended rate, not assumed
- [ ] Colorization reviewed for the known-wrong colors (skin, sky, signage)
- [ ] Before/after comparison rendered for the first 20 seconds

## Metadata
- [ ] Title under {max_title} characters
- [ ] Description pasted from {slug}.description.txt
- [ ] Tags pasted (under 500 characters total)
- [ ] Chapters added if the video is over 8 minutes
- [ ] End screen and cards set

## The two that protect the channel
- [ ] "Altered or synthetic content" disclosed in Studio
      (colorization is not on YouTube's exemption list; upscaling is. Disclosing costs
      nothing and non-disclosure risks YPP suspension.)
- [ ] Narration present, and it is the spine of the video rather than a caption over it
      (this is the reused-content policy's strongest allowed-side bullet)

## After publish
- [ ] {slug}.provenance.json filed somewhere you can find it in a hurry
- [ ] Checked for Content ID claims within 24h
- [ ] Short-form cut is 60s+ (TikTok minimum) and carries the narration, not just the reveal
"""


def write_outputs(
    outdir: Path, slug: str, description: str, tags: list[str],
    provenance: dict, title: str, identifier: str,
) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    description_path = outdir / f"{slug}.description.txt"
    description_path.write_text(
        description + "\n\n" + "TAGS (paste into the tags field)\n" + ", ".join(tags) + "\n",
        encoding="utf-8",
    )
    written.append(description_path)

    provenance_path = outdir / f"{slug}.provenance.json"
    provenance_path.write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    written.append(provenance_path)

    checklist_path = outdir / f"{slug}.checklist.md"
    checklist_path.write_text(
        CHECKLIST_TEMPLATE.format(
            title=title,
            details_url=ia_meta.details_url(identifier),
            today=date.today().isoformat(),
            slug=slug,
            max_title=MAX_TITLE_CHARS,
        ),
        encoding="utf-8",
    )
    written.append(checklist_path)

    return written


def selftest() -> int:
    """Offline checks. No network."""
    fake = {
        "metadata": {
            "identifier": "TestFilm1906",
            "title": "A Test Film",
            "year": "1906",
            "collection": ["prelinger", "silent_films"],
            "uploader": "someone@example.com",
        },
        "files": [
            {"name": "master.mov", "format": "QuickTime", "size": "900",
             "width": "4096", "height": "3072", "source": "original", "length": "404.5"},
            {"name": "small.mp4", "format": "h.264", "size": "50",
             "width": "640", "height": "480", "source": "derivative"},
        ],
    }

    year = extract_year(fake["metadata"])
    assert year == 1906, year

    best = ia_meta.best_video(fake)
    assert best is not None, "fixture has a video file, so this cannot be None"
    assert best["name"] == "master.mov", best

    cutoff = pd_cutoff_year(date(2026, 8, 26))
    assert cutoff == 1931, cutoff

    description = build_description(fake, "TestFilm1906", best, "A Test Film", year)
    assert "public domain in the United States" in description
    assert "archive.org/details/TestFilm1906" in description
    assert "4096x3072" in description
    assert len(description) <= MAX_DESCRIPTION_CHARS

    # A film too recent to be automatically PD must produce a warning, not a claim.
    recent = dict(fake)
    recent["metadata"] = dict(fake["metadata"], year="1999")
    warned = build_description(recent, "TestFilm1906", best, "T", 1999)
    assert "VERIFY BEFORE PUBLISHING" in warned
    assert "public domain in the United States" not in warned

    # Missing year must also warn rather than assert PD.
    noyear = build_description(fake, "TestFilm1906", best, "T", None)
    assert "VERIFY BEFORE PUBLISHING" in noyear

    tags = build_tags(year=1906)
    assert "1906" in tags and "1900s" in tags
    assert len(", ".join(tags)) <= MAX_TAG_CHARS, len(", ".join(tags))

    # Absent rights metadata must be recorded explicitly, not omitted.
    provenance = build_provenance(fake, "TestFilm1906", best, year)
    assert provenance["archive_rights_metadata_present"] is False
    assert "NO rights metadata" in provenance["rights_metadata_note"]
    assert provenance["archive_rights_metadata"]["licenseurl"] is None
    assert json.loads(json.dumps(provenance))  # must be serialisable

    withrights = dict(fake)
    withrights["metadata"] = dict(fake["metadata"], licenseurl="http://example/pd")
    p2 = build_provenance(withrights, "TestFilm1906", best, year)
    assert p2["archive_rights_metadata_present"] is True

    assert slugify("A Trip Down Market Street! (1906)") == "a-trip-down-market-street-1906"
    assert slugify("") == "untitled"

    # An item with no video at all must not crash.
    empty = build_description({"metadata": {}}, "X", None, "T", 1906)
    assert "no video file found" in empty

    print("make_metadata selftest: PASS")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate YouTube metadata + provenance for an archive.org film.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(__doc__ or ""),
    )
    parser.add_argument("identifier", nargs="?", help="archive.org item identifier")
    parser.add_argument("--outdir", default="metadata", help="output directory")
    parser.add_argument("--title", help="video title (defaults to the archive title)")
    parser.add_argument("--year", type=int, help="override the film's publication year")
    parser.add_argument("--file", help="use this filename instead of the largest video")
    parser.add_argument("--tags", default="", help="extra tags, comma separated")
    parser.add_argument("--notes", default="", help="extra description notes")
    parser.add_argument("--selftest", action="store_true", help="run offline checks and exit")
    args = parser.parse_args(argv)

    if args.selftest:
        return selftest()
    if not args.identifier:
        parser.error("identifier is required (or pass --selftest)")

    try:
        item = ia_meta.fetch_item(args.identifier)
    except (ia_meta.ArchiveError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    metadata = item.get("metadata", {})
    year = args.year if args.year is not None else extract_year(metadata)
    title = args.title or metadata.get("title") or args.identifier

    if len(title) > MAX_TITLE_CHARS:
        print(
            f"warning: title is {len(title)} chars, over YouTube's {MAX_TITLE_CHARS} limit",
            file=sys.stderr,
        )

    videos = ia_meta.video_files(item)
    if args.file:
        chosen = next((f for f in videos if f["name"] == args.file), None)
        if chosen is None:
            print(f"error: {args.file!r} is not a video file in this item.", file=sys.stderr)
            print("available:", file=sys.stderr)
            for f in videos:
                print(f"  {f['name']}", file=sys.stderr)
            return 1
    else:
        chosen = videos[0] if videos else None

    extra_tags = [t.strip() for t in args.tags.split(",") if t.strip()]
    description = build_description(
        item, args.identifier, chosen, title, year, args.notes
    )
    tags = build_tags(extra_tags, year)
    provenance = build_provenance(item, args.identifier, chosen, year)

    slug = slugify(title)
    written = write_outputs(
        Path(args.outdir), slug, description, tags, provenance, title, args.identifier
    )

    print(f"{title}")
    if chosen:
        print(
            f"  source file : {chosen['name']} "
            f"({chosen['width']}x{chosen['height']}, {ia_meta.human_size(chosen['size'])})"
        )
    print(f"  film year   : {year if year else 'UNKNOWN - verify before publishing'}")
    if year and year >= pd_cutoff_year():
        print(f"  !! {year} is not automatically PD (cutoff {pd_cutoff_year()}). Verify renewal.")
    for path in written:
        print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
