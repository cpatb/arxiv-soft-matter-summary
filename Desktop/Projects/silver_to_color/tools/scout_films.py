"""Find high-resolution public-domain film scans on archive.org and rank them.

Outputs a CSV you can sort in a spreadsheet, plus a short summary on stdout.

    python scout_films.py --era 1900 1930 --limit 60 --deep --out candidates.csv
    python scout_films.py --collection silent_films --min-size-gb 5 --deep
    python scout_films.py --topic "san francisco" --era 1900 1930
    python scout_films.py --selftest          # offline, no network

WHY IT USES advancedsearch AND NOT THE SCRAPE API
-------------------------------------------------
archive.org's Scrape API (`services/search/v1/scrape`) has a live, reproducible
server-side response cache whose key OMITS the query string, for `count` <= 250. Distinct
queries return byte-identical result sets, and responses belonging to *other users'*
queries have been observed leaking in. `advancedsearch.php` never exhibited this in
testing and is used here for everything. If you ever do need Scrape (collections over
10,000 items), pass `count` >= 251 and verify the payload against the query.

WHAT "HIGH RESOLUTION" ACTUALLY MEANS HERE
------------------------------------------
Item size alone is a weak proxy: a big item may just be a long, low-resolution film.
With `--deep` the script reads each candidate's file manifest and ranks on the real pixel
dimensions of its best *master* file, which is the number that matters.

Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
import urllib.parse
from dataclasses import dataclass, asdict, field

import ia_meta

ADVANCED_SEARCH = "https://archive.org/advancedsearch.php"

# advancedsearch refuses to page beyond this offset.
MAX_DEEP_PAGING = 10000

# Fields worth pulling back for every hit.
SEARCH_FIELDS = [
    "identifier", "title", "year", "date", "collection",
    "item_size", "downloads", "licenseurl", "subject",
]

# archive.org's own transcode ladder. Seeing only these means no master is present.
IA_DERIVATIVE_FORMATS = {
    "h.264 ia", "h.264", "ogg video", "512kb mpeg4", "256kb mpeg4",
    "thumbnail", "jpeg thumb", "item tile", "animated gif",
}

# Formats an uploader supplies, i.e. plausible masters.
MASTER_FORMATS = {
    "quicktime", "matroska", "mpeg4", "mpeg2", "mpeg2-ts", "mpeg1",
    "divx", "cinepack", "webm", "windows media", "flash video", "3gp", "hevc",
}

# In calendar 2026, US works published before this year are public domain (95-year term).
PD_CUTOFF_YEAR = 1931

# Bytes per pixel-second = size / (width * height * length). It separates a genuine film
# scan from a low-bitrate file that has been upscaled to look like one, which raw
# resolution cannot do. Observed values:
#     ~80    uncompressed 1080p
#     ~10    genuine 4K scan (Market Street 10.16, Code_Of_The_Sea 10.62)
#     0.6-1.5  compressed re-encode
#     ~0.14  modern AI upscale of a low-resolution source
# The upscale case is the one that matters: restoring someone else's upscale means
# inheriting their artifacts and adding your own.
BPPS_GENUINE_SCAN = 4.0     # at or above this, the pixels are probably real
BPPS_UPSCALE_SUSPECT = 0.5  # below this, the resolution is probably synthetic


@dataclass
class Candidate:
    identifier: str
    title: str = ""
    year: int | None = None
    collections: str = ""
    item_size_gb: float = 0.0
    downloads: int = 0
    licenseurl: str = ""
    # Filled by --deep
    best_file: str = ""
    best_format: str = ""
    width: int | None = None
    height: int | None = None
    megapixels: float | None = None
    master_size_gb: float | None = None
    duration_s: float | None = None
    bytes_per_pixel_second: float | None = None
    has_master: bool | None = None
    fits_kaggle_temp: bool | None = None
    score: float = 0.0
    flags: str = ""
    url: str = ""
    notes: list[str] = field(default_factory=list)


def build_query(
    era: tuple[int, int] | None = None,
    collection: str | None = None,
    topic: str | None = None,
    min_size_gb: float | None = None,
    exclude_restricted: bool = True,
) -> str:
    """Compose a Lucene query for archive.org's search endpoint."""
    parts = ["mediatype:movies"]

    if era:
        start, end = era
        if start > end:
            raise ValueError(f"era start {start} is after end {end}")
        # `year:` is uploader-supplied and slightly more inclusive than `date:`.
        parts.append(f"year:[{start} TO {end}]")
    if collection:
        parts.append(f"collection:{collection}")
    if topic:
        escaped = topic.replace('"', '\\"')
        parts.append(f'(title:"{escaped}" OR description:"{escaped}" OR subject:"{escaped}")')
    if min_size_gb:
        if min_size_gb < 0:
            raise ValueError("min_size_gb must be non-negative")
        parts.append(f"item_size:[{int(min_size_gb * 1e9)} TO 999999999999]")
    if exclude_restricted:
        parts.append("NOT collection:stream_only")
        parts.append("NOT access-restricted-item:true")

    return " AND ".join(parts)


def search(query: str, limit: int = 100, sort: str = "downloads desc",
           timeout: float = 60.0) -> list[dict]:
    """Run a search and return raw documents, paging as needed."""
    if limit < 1:
        raise ValueError("limit must be >= 1")

    documents: list[dict] = []
    page = 1
    rows = min(limit, 500)

    while len(documents) < limit:
        if (page - 1) * rows >= MAX_DEEP_PAGING:
            print(f"  note: stopping at archive.org's {MAX_DEEP_PAGING}-result paging limit",
                  file=sys.stderr)
            break

        params = [("q", query), ("rows", str(rows)), ("page", str(page)),
                  ("output", "json"), ("sort[]", sort)]
        params += [("fl[]", f) for f in SEARCH_FIELDS]
        url = f"{ADVANCED_SEARCH}?{urllib.parse.urlencode(params)}"

        payload = ia_meta.fetch_json(url, timeout=timeout)
        # advancedsearch signals deep-paging refusal as HTTP 200 with NO `response`
        # key, so a plain .get() would silently return zero results and look like a
        # finished sweep rather than a truncated one.
        if "response" not in payload:
            print(f"  note: archive.org refused page {page} "
                  f"(deep-paging limit ~{MAX_DEEP_PAGING}); returning what was collected",
                  file=sys.stderr)
            break
        response = payload["response"]
        batch = response.get("docs", [])
        if not batch:
            break
        documents.extend(batch)

        found = response.get("numFound", 0)
        if len(documents) >= found:
            break
        page += 1
        time.sleep(ia_meta.DEFAULT_DELAY_S)

    return documents[:limit]


def count_matches(query: str, timeout: float = 60.0) -> int:
    """Total matches, via the cheapest possible query (`rows=0`)."""
    url = f"{ADVANCED_SEARCH}?{urllib.parse.urlencode([('q', query), ('rows', '0'), ('output', 'json')])}"
    payload = ia_meta.fetch_json(url, timeout=timeout)
    return int(payload.get("response", {}).get("numFound", 0))


def is_derived(file_entry: dict) -> bool:
    """Whether a file was generated from another file in the same item.

    Two conditions, not one: archive.org labels some derived files
    ``source: "original"`` while still carrying an ``original`` pointer to their parent.
    The Market Street HEVC file does exactly that, so testing `source` alone
    misclassifies the flagship item's own derivative as a master.
    """
    return file_entry.get("source") == "derivative" or "original" in file_entry


def pick_master(document: dict) -> dict | None:
    """The best true master video file in an item, or None.

    Prefers a genuine uploader master; falls back to the largest video of any kind so a
    candidate is never silently dropped for having unusual metadata.
    """
    files = document.get("files", [])
    videos = ia_meta.video_files(document)
    if not videos:
        return None

    by_name = {f.get("name"): f for f in files}
    masters = []
    for video in videos:
        raw = by_name.get(video["name"], {})
        fmt = str(video.get("format", "")).lower()
        if is_derived(raw):
            continue
        if fmt in IA_DERIVATIVE_FORMATS and fmt not in MASTER_FORMATS:
            continue
        masters.append(video)

    return masters[0] if masters else videos[0]


def true_duration(document: dict) -> float | None:
    """Runtime in seconds, preferring a derivative's value.

    A master's `length` can be an fps-dependent artifact: the Market Street 4K master
    reports 404.5 s against its own derivative's 716.79 s, because the master pads
    16 fps content to 24 fps by doubling every second frame.
    """
    files = document.get("files", [])
    by_name = {f.get("name"): f for f in files}

    derivative_lengths, other_lengths = [], []
    for video in ia_meta.video_files(document):
        if video.get("length") is None:
            continue
        if is_derived(by_name.get(video["name"], {})):
            derivative_lengths.append(video["length"])
        else:
            other_lengths.append(video["length"])

    if derivative_lengths:
        return max(derivative_lengths)
    return max(other_lengths) if other_lengths else None


def score_candidate(candidate: Candidate) -> float:
    """Rank candidates. Resolution dominates; popularity is a weak tiebreak.

    Deliberately NOT a probability of anything -- it is an ordering heuristic, and the
    CSV keeps every input column so you can re-sort on whatever you actually care about.
    """
    score = 0.0

    if candidate.megapixels:
        # 2 MP (HD) scores 20, 12.6 MP (the 4096x3072 master) scores ~59.
        score += min(60.0, 20.0 * (candidate.megapixels / 2.0) ** 0.5)
    elif candidate.item_size_gb:
        # No dimensions known: fall back to size, capped well below a real measurement
        # so a measured item always outranks a guessed one.
        score += min(15.0, candidate.item_size_gb)

    if candidate.has_master:
        score += 15.0

    # Resolution is only worth what the bitrate behind it can support. An AI upscale
    # has the pixel count of a 4K scan and none of the information, so it must not
    # outrank a genuine lower-resolution master.
    bpps = candidate.bytes_per_pixel_second
    if bpps is not None:
        if bpps >= BPPS_GENUINE_SCAN:
            score += 15.0
        elif bpps >= BPPS_UPSCALE_SUSPECT:
            score += 3.0
        else:
            score -= 25.0

    if candidate.year and candidate.year < PD_CUTOFF_YEAR:
        score += 10.0

    if candidate.downloads > 0:
        # Diminishing returns: popularity signals a watchable print, but the goal is
        # finding good scans, not re-restoring the ten most-uploaded films on the site.
        score += min(10.0, 2.5 * (candidate.downloads / 1000.0) ** 0.3)

    if candidate.licenseurl and "publicdomain" in candidate.licenseurl.lower():
        score += 5.0

    # Practical penalty: an item that cannot be staged on a free-tier disk costs real
    # extra work (URL segmentation), so rank it below an equally good smaller scan.
    if candidate.fits_kaggle_temp is False:
        score -= 8.0

    return round(score, 2)


def assess(document: dict, base: Candidate) -> Candidate:
    """Fill in the deep fields from a full metadata document."""
    master = pick_master(document)
    files = document.get("files", [])
    by_name = {f.get("name"): f for f in files}

    if master:
        base.best_file = master["name"]
        base.best_format = str(master.get("format", ""))
        base.width = master.get("width")
        base.height = master.get("height")
        if base.width and base.height:
            base.megapixels = round(base.width * base.height / 1e6, 2)
        if master.get("size"):
            base.master_size_gb = round(master["size"] / 1e9, 2)
            # /kaggle/temp is ~60 GB and non-persistent; leave headroom for output.
            base.fits_kaggle_temp = base.master_size_gb < 50.0
        base.has_master = not is_derived(by_name.get(master["name"], {}))

    base.duration_s = true_duration(document)

    # Compute bpps from the master's SIZE but the derivative's LENGTH: a master's own
    # length can be an fps-dependent artifact (see true_duration).
    if (master and master.get("size") and base.width and base.height
            and base.duration_s and base.duration_s > 0):
        pixels_seconds = base.width * base.height * base.duration_s
        if pixels_seconds > 0:
            base.bytes_per_pixel_second = round(master["size"] / pixels_seconds, 2)

    flags = []
    if base.has_master is False:
        flags.append("no-master")
    if base.width and base.width < 1280:
        flags.append("low-res")
    if base.fits_kaggle_temp is False:
        flags.append("too-big-for-free-tier")
    bpps = base.bytes_per_pixel_second
    if bpps is not None:
        if bpps < BPPS_UPSCALE_SUSPECT:
            flags.append(f"LIKELY-UPSCALE-NOT-SCAN(bpps={bpps})")
        elif bpps < BPPS_GENUINE_SCAN:
            flags.append(f"re-encode-not-master(bpps={bpps})")
    if base.year and base.year >= PD_CUTOFF_YEAR:
        flags.append("VERIFY-RIGHTS")
    if not base.year:
        flags.append("no-year")

    # A duration disagreement between master and derivative is the padded-frames
    # signature -- worth surfacing, because it means free GPU savings via mpdecimate.
    lengths = {
        round(v["length"], 1)
        for v in ia_meta.video_files(document) if v.get("length")
    }
    if len(lengths) > 1:
        flags.append("duration-mismatch(check-for-padded-frames)")

    base.flags = ";".join(flags)
    base.url = ia_meta.details_url(base.identifier)
    base.score = score_candidate(base)
    return base


def to_candidate(document: dict) -> Candidate:
    """Convert one search hit into a Candidate."""
    identifier = document.get("identifier", "")
    year = None
    for key in ("year", "date"):
        value = document.get(key)
        if value is None:
            continue
        text = str(value if not isinstance(value, list) else value[0])
        digits = "".join(c for c in text[:4] if c.isdigit())
        if len(digits) == 4:
            year = int(digits)
            break

    size = document.get("item_size") or 0
    try:
        size_gb = round(float(size) / 1e9, 2)
    except (TypeError, ValueError):
        size_gb = 0.0

    try:
        downloads = int(document.get("downloads") or 0)
    except (TypeError, ValueError):
        downloads = 0

    title = document.get("title", "")
    if isinstance(title, list):
        title = title[0] if title else ""

    return Candidate(
        identifier=identifier,
        title=str(title)[:160],
        year=year,
        collections=",".join(ia_meta.as_list(document.get("collection")))[:120],
        item_size_gb=size_gb,
        downloads=downloads,
        licenseurl=str(document.get("licenseurl") or ""),
        url=ia_meta.details_url(identifier) if identifier else "",
    )


CSV_COLUMNS = [
    "score", "identifier", "title", "year", "width", "height", "megapixels",
    "duration_s", "bytes_per_pixel_second", "master_size_gb", "item_size_gb",
    "has_master", "best_format",
    "best_file", "fits_kaggle_temp", "downloads", "licenseurl", "collections",
    "flags", "url",
]


def console_safe(text: str) -> str:
    """Make a string printable on the active console encoding.

    archive.org titles carry accents, CJK and typographic quotes. A Windows console
    defaults to cp1252, where printing those raises UnicodeEncodeError and kills the run
    after the work is already done. The CSV is always written as UTF-8 and keeps the
    original characters; only this on-screen preview is degraded.
    """
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    try:
        text.encode(encoding)
        return text
    except (UnicodeEncodeError, LookupError):
        return text.encode(encoding, errors="replace").decode(encoding, errors="replace")


def write_csv(candidates: list[Candidate], path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(asdict(candidate))


def selftest() -> int:
    """Offline checks against synthetic documents. No network."""
    query = build_query(era=(1900, 1930), collection="silent_films", min_size_gb=5)
    assert "mediatype:movies" in query
    assert "year:[1900 TO 1930]" in query
    assert "collection:silent_films" in query
    assert "item_size:[5000000000 TO 999999999999]" in query
    assert "NOT access-restricted-item:true" in query

    try:
        build_query(era=(1930, 1900))
    except ValueError:
        pass
    else:
        raise AssertionError("reversed era accepted")

    # Trap 2: source == "original" plus an `original` pointer still means derived.
    assert is_derived({"source": "derivative"})
    assert is_derived({"source": "original", "original": "master.mov"})
    assert not is_derived({"source": "original"})
    assert not is_derived({})

    document = {
        "metadata": {"identifier": "X", "title": "T", "year": "1906"},
        "files": [
            {"name": "master.mov", "format": "QuickTime", "source": "original",
             "size": "91674081757", "width": "4096", "height": "3072", "length": "404.5"},
            {"name": "web.mp4", "format": "h.264", "source": "derivative",
             "original": "master.mov", "size": "62548751",
             "width": "640", "height": "480", "length": "716.79"},
            {"name": "x.hevc.mp4", "format": "HEVC", "source": "original",
             "original": "master.mov", "size": "15428269828"},
        ],
    }

    master = pick_master(document)
    assert master is not None and master["name"] == "master.mov", master
    # The HEVC must not be chosen despite source=="original".
    assert master["name"] != "x.hevc.mp4"

    # Trap 3: prefer the derivative's duration over the master's padded one.
    assert true_duration(document) == 716.79, true_duration(document)

    base = to_candidate({"identifier": "X", "title": "T", "year": "1906",
                         "item_size": "102400000000", "downloads": "5000",
                         "collection": ["prelinger"]})
    assert base.year == 1906
    assert base.item_size_gb == 102.4

    assessed = assess(document, base)
    # bpps uses the master's size with the DERIVATIVE's length (716.79 s, not 404.5).
    expected_bpps = round(91674081757 / (4096 * 3072 * 716.79), 2)
    assert assessed.bytes_per_pixel_second == expected_bpps, assessed.bytes_per_pixel_second
    assert 9.0 < expected_bpps < 12.0, expected_bpps   # genuine-scan signature
    assert "UPSCALE" not in assessed.flags
    assert assessed.width == 4096 and assessed.height == 3072
    assert assessed.megapixels == 12.58, assessed.megapixels
    assert assessed.has_master is True
    assert assessed.fits_kaggle_temp is False, "91 GB must not read as free-tier friendly"
    assert "too-big-for-free-tier" in assessed.flags
    assert "duration-mismatch" in assessed.flags
    assert assessed.score > 0

    # A modern item must be flagged for rights, not silently ranked.
    modern = to_candidate({"identifier": "Y", "year": "1999", "item_size": "1000000000"})
    assert assess({"files": []}, modern).flags.count("VERIFY-RIGHTS") == 1

    # Resolution is meant to dominate: the 12.6 MP master outranks a 3.1 MP scan even
    # though it is too big to stage locally, because it genuinely is the better source
    # and URL segmentation handles the size. The size penalty is a tiebreak, not a veto.
    smaller = assess(
        {"files": [{"name": "a.mov", "format": "QuickTime", "source": "original",
                    "size": "2500000000", "width": "2048", "height": "1536"}]},
        to_candidate({"identifier": "S", "year": "1906", "downloads": "5000"}),
    )
    assert smaller.fits_kaggle_temp is True
    assert assessed.score > smaller.score, (assessed.score, smaller.score)

    # But the penalty must actually bite: the same item, small enough to stage, scores
    # higher than its oversized twin.
    twin = assess(
        {"files": [{"name": "a.mov", "format": "QuickTime", "source": "original",
                    "size": "9000000000", "width": "4096", "height": "3072"}]},
        to_candidate({"identifier": "T", "year": "1906", "downloads": "5000"}),
    )
    oversized = assess(
        {"files": [{"name": "a.mov", "format": "QuickTime", "source": "original",
                    "size": "91000000000", "width": "4096", "height": "3072"}]},
        to_candidate({"identifier": "U", "year": "1906", "downloads": "5000"}),
    )
    assert twin.fits_kaggle_temp is True and oversized.fits_kaggle_temp is False
    assert twin.score == oversized.score + 8.0, (twin.score, oversized.score)

    # Resolution must outrank popularity: an obscure 4K scan beats a famous 640x480 one.
    obscure = assess(
        {"files": [{"name": "m.mov", "format": "QuickTime", "source": "original",
                    "size": "9000000000", "width": "4096", "height": "3072"}]},
        to_candidate({"identifier": "O", "year": "1906", "downloads": "12"}),
    )
    famous = assess(
        {"files": [{"name": "m.mp4", "format": "MPEG4", "source": "original",
                    "size": "300000000", "width": "640", "height": "480"}]},
        to_candidate({"identifier": "F", "year": "1906", "downloads": "800000"}),
    )
    assert obscure.score > famous.score, (obscure.score, famous.score)

    # An AI upscale has 4K pixels and no information behind them. It must be flagged
    # and must NOT outrank a genuine lower-resolution master.
    upscale = assess(
        {"files": [
            {"name": "up.mp4", "format": "MPEG4", "source": "original",
             "size": "600000000", "width": "3840", "height": "2160"},
            {"name": "web.mp4", "format": "h.264", "source": "derivative",
             "original": "up.mp4", "size": "50000000", "length": "5000"},
        ]},
        to_candidate({"identifier": "UP", "year": "1930", "downloads": "9000"}),
    )
    assert upscale.bytes_per_pixel_second is not None
    assert upscale.bytes_per_pixel_second < BPPS_UPSCALE_SUSPECT, upscale.bytes_per_pixel_second
    assert "LIKELY-UPSCALE-NOT-SCAN" in upscale.flags, upscale.flags

    genuine_hd = assess(
        {"files": [
            {"name": "m.mov", "format": "QuickTime", "source": "original",
             "size": "40000000000", "width": "1920", "height": "1080"},
            {"name": "w.mp4", "format": "h.264", "source": "derivative",
             "original": "m.mov", "size": "5000000", "length": "1200"},
        ]},
        to_candidate({"identifier": "HD", "year": "1930", "downloads": "9000"}),
    )
    assert genuine_hd.bytes_per_pixel_second >= BPPS_GENUINE_SCAN
    assert genuine_hd.score > upscale.score, (genuine_hd.score, upscale.score)

    # No video at all must not crash.
    empty = assess({"files": []}, to_candidate({"identifier": "Z", "year": "1910"}))
    assert empty.width is None and empty.score >= 0

    # Non-ASCII titles must not be able to kill a run after the work is done.
    assert console_safe("Café — 東京") is not None
    assert console_safe("plain ascii") == "plain ascii"

    assert len(CSV_COLUMNS) == len(set(CSV_COLUMNS))
    for column in CSV_COLUMNS:
        assert column in asdict(assessed), f"CSV column {column} missing from Candidate"

    print("scout_films selftest: PASS")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Find and rank high-resolution public-domain film scans on archive.org."
    )
    parser.add_argument("--era", nargs=2, type=int, metavar=("START", "END"),
                        default=[1900, 1930], help="publication year range")
    parser.add_argument(
        "--collection",
        help="archive.org collection. VERIFIED to exist: prelinger, feature_films, "
             "silent_films, academic_films, home_movies, universal_newsreels, nasa, "
             "sabucat_trailers, us-national-archives-youtube, Film_Noir (capital F and N). "
             "These do NOT exist: newsreels, sabucat, AVGeeks, usnationalarchives, film_noir.")
    parser.add_argument("--topic", help="keyword in title/description/subject")
    parser.add_argument("--min-size-gb", type=float,
                        help="only items at least this large (a crude scan-quality proxy)")
    parser.add_argument("--limit", type=int, default=50, help="candidates to retrieve")
    parser.add_argument("--sort", default="downloads desc",
                        help="archive.org sort, e.g. 'item_size desc'")
    parser.add_argument("--deep", action="store_true",
                        help="fetch each item's file manifest for real resolution "
                             "(one request per candidate; much better ranking)")
    parser.add_argument("--out", default="candidates.csv", help="output CSV path")
    parser.add_argument("--selftest", action="store_true", help="offline checks, then exit")
    args = parser.parse_args(argv)

    if args.selftest:
        return selftest()

    try:
        query = build_query(
            era=(args.era[0], args.era[1]),
            collection=args.collection,
            topic=args.topic,
            min_size_gb=args.min_size_gb,
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    print(console_safe(f"query: {query}"))
    try:
        total = count_matches(query)
        print(f"matches: {total:,}")
        if total == 0:
            print("No matches. Try widening the era or dropping --min-size-gb.")
            return 0
        documents = search(query, limit=args.limit, sort=args.sort)
    except ia_meta.ArchiveError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"retrieved: {len(documents)}")
    candidates = [to_candidate(d) for d in documents]

    if args.deep:
        print(f"fetching file manifests for {len(candidates)} items "
              f"(~{len(candidates) * ia_meta.DEFAULT_DELAY_S:.0f}s)...")
        for index, candidate in enumerate(candidates, 1):
            try:
                document = ia_meta.fetch_item(candidate.identifier)
                assess(document, candidate)
            except (ia_meta.ArchiveError, ValueError) as exc:
                candidate.flags = (candidate.flags + ";metadata-failed").strip(";")
                candidate.notes.append(str(exc)[:200])
            if index % 10 == 0:
                print(f"  {index}/{len(candidates)}")
            time.sleep(ia_meta.DEFAULT_DELAY_S)
    else:
        for candidate in candidates:
            candidate.score = score_candidate(candidate)

    candidates.sort(key=lambda c: c.score, reverse=True)
    write_csv(candidates, args.out)

    print(f"\nwrote {args.out}\n")
    header = (f"{'score':>6}  {'resolution':>11}  {'size':>8}  {'bpps':>6}  "
              f"{'year':>4}  title")
    print(header)
    print("-" * min(110, len(header) + 50))
    for candidate in candidates[:15]:
        resolution = (f"{candidate.width}x{candidate.height}"
                      if candidate.width and candidate.height else "-")
        size = f"{candidate.master_size_gb or candidate.item_size_gb:.1f}GB"
        bpps = (f"{candidate.bytes_per_pixel_second:.2f}"
                if candidate.bytes_per_pixel_second is not None else "-")
        print(console_safe(
            f"{candidate.score:>6.1f}  {resolution:>11}  {size:>8}  {bpps:>6}  "
            f"{candidate.year or '-':>4}  {candidate.title[:52]}"))
        if candidate.flags:
            print(f"{'':>6}  flags: {candidate.flags}")

    print("\nRights reminder: `year` on archive.org is uploader-supplied and unvalidated "
          "(a 2026-dated\nPrelinger short was observed). Treat it as a shortlist filter, "
          "never as evidence of\npublication date. Verify against the Catalog of Copyright "
          "Entries or AFI before publishing.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
