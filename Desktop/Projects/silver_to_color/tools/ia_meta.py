"""Shared archive.org metadata access for the Silver-to-Color tools.

Standard library only. The one network dependency is archive.org's public metadata
endpoint, which is documented, unauthenticated, and returns JSON:

    https://archive.org/metadata/<identifier>

Every function here degrades to a clear exception rather than a silent wrong answer.
Nothing in this module writes files or mutates global state.
"""

from __future__ import annotations

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request

# archive.org asks that automated clients identify themselves.
USER_AGENT = "silver-to-color/1.0 (public-domain film restoration; contact via archive.org)"

METADATA_URL = "https://archive.org/metadata/{identifier}"
DETAILS_URL = "https://archive.org/details/{identifier}"
DOWNLOAD_URL = "https://archive.org/download/{identifier}/{filename}"

# Politeness delay between successive archive.org requests, in seconds.
# archive.org publishes no hard rate limit for the metadata endpoint; this is a
# self-imposed floor so a scouting sweep cannot hammer the service.
DEFAULT_DELAY_S = 0.5

# archive.org identifiers are constrained to these characters.
_IDENTIFIER_RE = re.compile(r"^[A-Za-z0-9._-]{1,100}$")


class ArchiveError(RuntimeError):
    """Raised when archive.org cannot be reached or returns something unusable."""


def validate_identifier(identifier: object) -> str:
    """Return the identifier unchanged, or raise if it cannot be a real one.

    Catching this here means a typo fails immediately with a clear message rather
    than becoming a confusing 404 several steps later.
    """
    if not isinstance(identifier, str):
        raise TypeError(f"identifier must be str, got {type(identifier).__name__}")
    identifier = identifier.strip()
    if not identifier:
        raise ValueError("identifier is empty")
    if not _IDENTIFIER_RE.match(identifier):
        raise ValueError(
            f"{identifier!r} is not a valid archive.org identifier "
            "(expected 1-100 chars of A-Z a-z 0-9 . _ -)"
        )
    return identifier


def fetch_json(url: str, timeout: float = 30.0, retries: int = 3) -> dict:
    """GET a URL and parse JSON, retrying on transient failures.

    Retries on 5xx and on connection errors with exponential backoff. Does NOT retry
    on 4xx, because a 404 will still be a 404 the second time.
    """
    if retries < 1:
        raise ValueError("retries must be >= 1")

    last_error: Exception | None = None
    for attempt in range(retries):
        request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                raw = response.read()
            break
        except urllib.error.HTTPError as exc:
            # Client errors are permanent; server errors are worth retrying.
            if exc.code < 500:
                raise ArchiveError(f"HTTP {exc.code} for {url}") from exc
            last_error = exc
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last_error = exc

        if attempt < retries - 1:
            time.sleep(2.0**attempt)
    else:
        raise ArchiveError(f"could not fetch {url} after {retries} attempts: {last_error}")

    try:
        return json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ArchiveError(f"response from {url} was not valid JSON: {exc}") from exc


def fetch_item(identifier: str, timeout: float = 30.0) -> dict:
    """Fetch one archive.org item's full metadata document.

    An identifier that does not exist returns ``{}`` from archive.org rather than a
    404, so that case is turned into an explicit error here.
    """
    identifier = validate_identifier(identifier)
    document = fetch_json(METADATA_URL.format(identifier=identifier), timeout=timeout)

    # The metadata endpoint has THREE states, not two, and all return HTTP 200:
    #   {}                  -> the item does not exist
    #   {"is_dark": true}   -> the item was withdrawn, typically after a DMCA claim
    #   metadata present    -> live
    # Collapsing the middle case into "missing" would lose the single most useful
    # rights signal archive.org emits: a dark item is one somebody successfully
    # claimed, which is the opposite of a safe public-domain candidate.
    if not document:
        raise ArchiveError(
            f"archive.org has no item {identifier!r} (empty response). "
            "Check the identifier spelling."
        )
    if document.get("is_dark"):
        raise ArchiveError(
            f"item {identifier!r} is DARK: archive.org has withdrawn it from public "
            "access, which usually follows a copyright claim. Do not restore it."
        )
    if "metadata" not in document:
        raise ArchiveError(
            f"archive.org returned no metadata block for {identifier!r}"
        )
    return document


def item_status(identifier: str, timeout: float = 30.0) -> str:
    """'live', 'dark' or 'missing', without raising.

    Useful for sweeping a candidate list: a dark item is not merely unavailable, it is
    evidence that a claim succeeded against it.
    """
    try:
        document = fetch_json(
            METADATA_URL.format(identifier=validate_identifier(identifier)),
            timeout=timeout,
        )
    except ArchiveError:
        return "missing"
    if not document:
        return "missing"
    if document.get("is_dark"):
        return "dark"
    return "live" if "metadata" in document else "missing"


def _to_int(value) -> int | None:
    """Best-effort int conversion; archive.org returns numbers as strings."""
    if value is None:
        return None
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None


def _to_float(value) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    # Durations sometimes arrive as "HH:MM:SS" or "MM:SS" instead of seconds.
    if ":" in text:
        parts = text.split(":")
        try:
            seconds = 0.0
            for part in parts:
                seconds = seconds * 60.0 + float(part)
            return seconds
        except ValueError:
            return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def as_list(value) -> list[str]:
    """archive.org returns single-valued fields as str and multi-valued as list."""
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def video_files(document: object) -> list[dict]:
    """Return the item's video files, normalised and richest-first.

    Each returned dict has: name, format, size (bytes, int or None), width, height
    (int or None), length (seconds, float or None), source ('original'/'derivative'/
    unknown), and pixels (width*height, or None).

    Sorting is by pixel count then file size, so ``video_files(doc)[0]`` is the best
    available scan. Files with no resolution metadata sort last but are still
    returned -- archive.org does not always populate width/height, and discarding
    those would silently hide real masters.
    """
    if not isinstance(document, dict):
        raise TypeError("document must be the dict returned by fetch_item()")
    item: dict = document

    video_extensions = (
        ".mp4", ".mkv", ".avi", ".mov", ".m4v", ".mpg", ".mpeg",
        ".ogv", ".webm", ".ts", ".wmv", ".flv",
    )

    results: list[dict] = []
    for entry in item.get("files", []):
        name = entry.get("name", "")
        if not name.lower().endswith(video_extensions):
            continue

        width = _to_int(entry.get("width"))
        height = _to_int(entry.get("height"))
        pixels = width * height if (width and height) else None

        results.append(
            {
                "name": name,
                "format": entry.get("format", "unknown"),
                "size": _to_int(entry.get("size")),
                "width": width,
                "height": height,
                "pixels": pixels,
                "length": _to_float(entry.get("length")),
                "source": entry.get("source", "unknown"),
            }
        )

    # -1 sorts unknown resolution/size last without dropping the entry.
    results.sort(key=lambda f: (f["pixels"] or -1, f["size"] or -1), reverse=True)
    return results


def best_video(document: dict) -> dict | None:
    """The highest-resolution video file in an item, or None if it has none."""
    files = video_files(document)
    return files[0] if files else None


def download_url(identifier: str, filename: str) -> str:
    """Direct download URL for one file inside an item."""
    identifier = validate_identifier(identifier)
    return DOWNLOAD_URL.format(
        identifier=identifier, filename=urllib.parse.quote(filename)
    )


def details_url(identifier: str) -> str:
    """Human-facing archive.org page for an item."""
    return DETAILS_URL.format(identifier=validate_identifier(identifier))


def human_size(num_bytes: int | None) -> str:
    """Format a byte count for display. Returns '?' for unknown."""
    if num_bytes is None:
        return "?"
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024.0 or unit == "TB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{int(size)} B"
        size /= 1024.0
    return f"{size:.1f} TB"


def human_duration(seconds: float | None) -> str:
    """Format a duration as H:MM:SS or M:SS. Returns '?' for unknown."""
    if seconds is None or seconds < 0:
        return "?"
    total = int(round(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    if hours:
        return f"{hours}:{minutes:02d}:{secs:02d}"
    return f"{minutes}:{secs:02d}"
