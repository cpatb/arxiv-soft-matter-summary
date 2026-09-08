"""Search DataCite for licensable rheology datasets that could close capability gap 3.

WHY DataCite AND NOT ZENODO. `NEXT_STEPS.md` item 20 records the Zenodo REST API as the
search route ("returns HTTP 504 to a default python-urllib User-Agent and works with an
explicit one"). On 2026-09-08 it returns 504 to an explicit User-Agent as well, from two
independent clients. DataCite is the DOI registration agency BEHIND Zenodo, so the same
records are indexed there with the same licence metadata, plus figshare, Dryad, 4TU and
every other DataCite member repository. It answers.

WHAT THIS DOES AND DOES NOT DECIDE. This is a SEARCH, not a gate. It writes a ranked
candidate list with each record's licence, repository and file manifest where the API
exposes one. Nothing here admits a dataset: that requires reading the deposit, reading
the source paper, and a pre-registration, exactly as items 20/22/24/25 did.

Run:  python scripts/datacite_candidate_search.py
Out:  data/candidate_search_2026-09-08.md  (and .json beside it)
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

UA = "biorheo/0.1 (research; mailto:flippinlizard@gmail.com)"
API = "https://api.datacite.org/dois"
OUT_DIR = Path(__file__).resolve().parent.parent / "data"

# Licences this project may build a public demo on. Anything else is not licensable
# here even if the data is public: item 18's finding was that "posted" is not "permitted".
OK_LICENCES = {"cc-by-4.0", "cc-by-3.0", "cc0-1.0", "cc-by-sa-4.0", "mit", "apache-2.0"}

# Each query names the capability gap it would close if it hit.
QUERIES: list[tuple[str, str]] = [
    (
        "protein-flowcurve",
        '(titles.title:(rheolog* OR viscosit* OR "flow curve") AND '
        'titles.title:(protein OR antibody OR "monoclonal" OR mAb OR biologic*))',
    ),
    (
        "protein-any",
        '(descriptions.description:("high concentration" AND (antibody OR protein)) AND '
        "descriptions.description:(viscosity OR rheolog*))",
    ),
    (
        "injection-force",
        '(titles.title:(injectab* OR "injection force" OR syringe OR needle) AND '
        "descriptions.description:(viscosity OR rheolog*))",
    ),
    (
        "blood-flowcurve",
        '(titles.title:(blood OR plasma OR "whole blood" OR haemorheolog* OR hemorheolog*) '
        'AND descriptions.description:(viscosity OR "shear rate"))',
    ),
    (
        "mucus-saliva-synovial",
        "(titles.title:(mucus OR saliva OR synovial OR sputum OR "
        "\"vitreous\" OR semen) AND descriptions.description:(rheolog* OR viscosit*))",
    ),
    (
        "yield-stress-highrate",
        '(titles.title:("shear rate" OR "flow curve" OR rheometry) AND '
        "descriptions.description:(\"yield stress\" AND (suspension OR gel OR microgel)))",
    ),
]

FILE_HINT = re.compile(r"\.(csv|txt|dat|tsv|xlsx|json)\b", re.I)


def fetch(url: str, tries: int = 4) -> dict:
    last: Exception | None = None
    for k in range(tries):
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:  # noqa: PERF203
            last = exc
            time.sleep(2.0 * (k + 1))
    raise RuntimeError(f"DataCite unreachable after {tries} tries: {last}")


def licence_of(attrs: dict) -> str:
    rights = attrs.get("rightsList") or []
    ids = [str(r.get("rightsIdentifier") or "").lower() for r in rights]
    ids = [i for i in ids if i]
    if not ids:
        # A deposit with no stated licence is NOT permission -- item 18's whole finding.
        return "NONE-STATED"
    return ids[0]


def search(name: str, query: str, size: int = 50) -> list[dict]:
    url = f"{API}?query={urllib.parse.quote(query)}&page[size]={size}&resource-type-id=dataset"
    data = fetch(url)
    rows = []
    for hit in data.get("data", []):
        a = hit["attributes"]
        title = (a.get("titles") or [{}])[0].get("title", "")
        desc = " ".join(
            (d.get("description") or "") for d in (a.get("descriptions") or [])
        )
        rows.append(
            {
                "query": name,
                "doi": a.get("doi", ""),
                "title": title.strip(),
                "year": a.get("publicationYear"),
                "client": (hit.get("relationships", {}).get("client", {}).get("data") or {}).get("id"),
                "licence": licence_of(a),
                "url": a.get("url", ""),
                "desc": re.sub(r"<[^>]+>", " ", desc)[:900].strip(),
                "has_file_hint": bool(FILE_HINT.search(desc)),
            }
        )
    return rows


def main() -> int:
    all_rows: list[dict] = []
    per_query_total: dict[str, int] = {}
    for name, q in QUERIES:
        rows = search(name, q)
        per_query_total[name] = len(rows)
        all_rows.extend(rows)
        print(f"  {name:24s} {len(rows):3d} hits", flush=True)
        time.sleep(0.5)

    # De-duplicate on DOI, keeping the first query that found it.
    seen: set[str] = set()
    uniq: list[dict] = []
    for r in all_rows:
        if r["doi"] in seen:
            continue
        seen.add(r["doi"])
        uniq.append(r)

    licensable = [r for r in uniq if r["licence"] in OK_LICENCES]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "candidate_search_2026-09-08.json").write_text(
        json.dumps(
            {
                "generated_by": "scripts/datacite_candidate_search.py",
                "api": API,
                "queries": {n: q for n, q in QUERIES},
                "hits_per_query": per_query_total,
                "n_raw": len(all_rows),
                "n_unique_doi": len(uniq),
                "n_licensable": len(licensable),
                "rows": uniq,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    lines = [
        "# Candidate search for a licensable biofluid flow curve (2026-09-08)",
        "",
        "Generated by `scripts/datacite_candidate_search.py`. Every number below is printed",
        "by that script.",
        "",
        "**Route.** The Zenodo REST API returned HTTP 504 to an explicit User-Agent from two",
        "independent clients on this date, so the route recorded in item 20 is dead. DataCite",
        "is the registration agency behind Zenodo and indexes the same records with the same",
        "licence metadata, plus figshare, Dryad and 4TU. Queries are restricted to",
        "`resource-type-id=dataset`.",
        "",
        "**This search decides nothing.** A hit here is a candidate to READ. Admission still",
        "needs the deposit, the source paper, and a pre-registration written before the",
        "generator, as items 20/22/24/25 all did.",
        "",
        f"- queries run: **{len(QUERIES)}**",
        f"- rows returned: **{len(all_rows)}**",
        f"- unique DOIs: **{len(uniq)}**",
        f"- of those, licence in {sorted(OK_LICENCES)}: **{len(licensable)}**",
        "",
        "## Hits per query",
        "",
        "| query | what it would close | hits |",
        "|---|---|---|",
    ]
    gap = {
        "protein-flowcurve": "capability gap 3: the wedge fluid with a CURVE",
        "protein-any": "same, via description rather than title",
        "injection-force": "an external force anchor for the needle chain",
        "blood-flowcurve": "a licensable replacement for the Horner records",
        "mucus-saliva-synovial": "a second biofluid family",
        "yield-stress-highrate": "a licensable microgel at needle rates",
    }
    for n, _q in QUERIES:
        lines.append(f"| `{n}` | {gap[n]} | {per_query_total[n]} |")

    lines += ["", "## Licensable candidates", ""]
    if not licensable:
        lines.append("**None.** Every unique hit carried a licence outside the allowed set.")
    else:
        lines.append("| DOI | year | licence | repo | title |")
        lines.append("|---|---|---|---|---|")
        for r in sorted(licensable, key=lambda x: -(x["year"] or 0)):
            t = r["title"].replace("|", "/")[:110]
            lines.append(
                f"| `{r['doi']}` | {r['year']} | {r['licence']} | {r['client']} | {t} |"
            )

    nonstated = [r for r in uniq if r["licence"] == "NONE-STATED"]
    lines += [
        "",
        f"## Deposits with NO stated licence: {len(nonstated)}",
        "",
        "Recorded because item 18's finding was that a posted file with no licence is not",
        "permission. They are not candidates.",
        "",
    ]

    (OUT_DIR / "candidate_search_2026-09-08.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(f"\nunique DOIs {len(uniq)}, licensable {len(licensable)}")
    print(f"wrote {OUT_DIR / 'candidate_search_2026-09-08.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
