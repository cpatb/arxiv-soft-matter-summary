"""DOI syntax and resolution, with a shipped cache so the registry loads offline.

A record's source exists only once its DOI has resolved through the DOI handle
system (``doi.org/api/handles``, which covers Crossref and DataCite alike; the
Carbopol dataset is DataCite and 404s on Crossref). Metadata is then fetched from
whichever registrar owns the prefix, purely for the human-readable label. Every
resolution is written to ``doi_cache.json`` with the date it happened, and the
loader's default ``allow_network=False`` reads only that cache, so a test run
cannot silently depend on the network and a record cannot cite a DOI nobody
looked up.
"""

from __future__ import annotations

import json
import re
import urllib.request
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path

__all__ = ["DOI_RE", "ResolvedDOI", "is_doi_syntax", "resolve_doi", "CACHE_PATH"]

# Crossref's own recommended pattern is looser; this is the structural minimum: the
# directory indicator 10, a 4-9 digit registrant, a slash, a non-empty suffix. It
# deliberately rejects "doi:" and "https://doi.org/" prefixes so records store the
# bare identifier.
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")

CACHE_PATH = Path(__file__).with_name("doi_cache.json")
_UA = "biorheo/0.0.1 (https://github.com/flippinlizard; mailto:flippinlizard@gmail.com)"


@dataclass(frozen=True)
class ResolvedDOI:
    doi: str
    resolved_on: str  # ISO date
    url: str  # what the handle system pointed at when resolved
    registrar: str  # "crossref" | "datacite" | "unknown"
    label: str  # human-readable citation line, from the registrar's metadata


def is_doi_syntax(s: str) -> bool:
    return isinstance(s, str) and DOI_RE.match(s) is not None


def _load_cache() -> dict[str, dict]:
    if not CACHE_PATH.exists():
        return {}
    return json.loads(CACHE_PATH.read_text(encoding="utf-8"))


def _save_cache(cache: dict[str, dict]) -> None:
    CACHE_PATH.write_text(json.dumps(cache, indent=1, sort_keys=True) + "\n",
                          encoding="utf-8", newline="\n")


def _get_json(url: str, timeout: float = 20.0) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": _UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.load(r)


def _handle_url(doi: str) -> str | None:
    """The URL the DOI handle system holds for ``doi``, or None if unregistered."""
    try:
        d = _get_json(f"https://doi.org/api/handles/{doi}")
    except Exception:
        return None
    if d.get("responseCode") != 1:
        return None
    for v in d.get("values", []):
        if v.get("type") == "URL":
            return str(v["data"]["value"])
    return None


def _label(doi: str) -> tuple[str, str]:
    """(registrar, one-line citation) from Crossref, then DataCite; never raises."""
    try:
        m = _get_json(f"https://api.crossref.org/works/{doi}")["message"]
        authors = ", ".join(
            f"{a.get('family', '')}" for a in m.get("author", [])[:3]) or "(no author)"
        year = (m.get("issued", {}).get("date-parts") or [[None]])[0][0]
        title = (m.get("title") or [""])[0]
        venue = (m.get("container-title") or [""])[0]
        return "crossref", f"{authors} ({year}). {title}. {venue}".strip()
    except Exception:
        pass
    try:
        a = _get_json(f"https://api.datacite.org/dois/{doi}")["data"]["attributes"]
        authors = ", ".join(c.get("familyName") or c.get("name", "") for c in a.get("creators", [])[:3])
        title = (a.get("titles") or [{}])[0].get("title", "")
        return "datacite", f"{authors} ({a.get('publicationYear')}). {title}. {a.get('publisher', '')}".strip()
    except Exception:
        return "unknown", ""


def resolve_doi(doi: str, allow_network: bool = False) -> ResolvedDOI | None:
    """Return the cached resolution of ``doi``; look it up and cache it if allowed.

    ``None`` means "not resolved": bad syntax, absent from the cache with the network
    disallowed, or the handle system has no such DOI.
    """
    if not is_doi_syntax(doi):
        return None
    cache = _load_cache()
    if doi in cache:
        return ResolvedDOI(**cache[doi])
    if not allow_network:
        return None
    url = _handle_url(doi)
    if url is None:
        return None
    registrar, label = _label(doi)
    r = ResolvedDOI(doi=doi, resolved_on=date.today().isoformat(), url=url,
                    registrar=registrar, label=label)
    cache[doi] = asdict(r)
    _save_cache(cache)
    return r
