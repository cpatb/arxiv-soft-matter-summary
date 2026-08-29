# Which PDF is current

Until 2026-08-29 the newest builds of both documents sat in the **Projects root**
while their `.tex` sources sat here, so `paper1/*.pdf` was the older file and the
root copy was the one anyone had actually been reading.

Both are kept. The plain name is the current build; the date-stamped name is the
older one that used to occupy it.

| document | current (was in Projects root) | superseded (was here) |
|---|---|---|
| referee review | `paper1_referee_review.pdf` (2026-06-30) | `paper1_referee_review.2026-04-29.pdf` |
| research approach | `paper1_research_approach.pdf` (2026-06-30) | `paper1_research_approach.2026-05-13.pdf` |

Neither pair is byte-identical. Rebuild either from its `.tex` with
`latexmk -pdf` in this directory.
