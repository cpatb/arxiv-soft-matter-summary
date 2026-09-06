# Copy-paste session prompt (biorheo / biofluids venture)

Paste everything below the line into a fresh Claude Code session started in `C:\Users\compl`.

---

Autonomous work session on `biorheo` (Desktop\Projects\biorheo), the product code of the biofluids venture. Work until you exhaust the session limit; do not stop early, do not ask me questions, do not wait for me. Be efficient with tokens but never trade correctness for speed: a wrong answer that runs is worse than no answer.

Resume protocol, in this order, before touching anything:
1. Read `Documents\Obsidian Vault\Status\Biofluids Venture.md` §1–3 (state, findings, next actions).
2. Read `Desktop\Projects\biorheo\NEXT_STEPS.md` fully and run its "Before touching anything" block. If the test count or top commit differs from what it says, find out why first.
3. Read the memory files named in `~/.claude/projects/C--Users-compl/memory/MEMORY.md` under "Biofluids venture / biorheo".

Then execute NEXT_STEPS.md "Open work, in order": finish any item marked in progress, then take the next one. For each item: pre-register the gate before running it, never loosen a bar after seeing the result, print every number that appears in prose from a script, add or update tests that pin the result, run `python -m pytest -q -m "not slow"` green, then commit from `C:\Users\compl` with `git add Desktop/Projects/biorheo && git -c core.safecrlf=false commit -m "..."`. A failed gate is a result: record it as FAILED with the numbers and move on. Run the matching reviewer (physics-referee, science-code-reviewer, units-auditor, numerics-verifier) on anything that will be believed, and act on its findings. When an item turns out not to be runnable or its premise is wrong, close it in NEXT_STEPS.md with the reason instead of working around it.

Handoff protocol. When roughly 85 % of the session budget is used, or earlier if a natural stopping point arrives, stop starting new items and do all of the following, then commit:
- `NEXT_STEPS.md`: mark done/closed items with dates and results, rewrite the next open item so a model with zero context can start it, update the expected test count and top commit in "Before touching anything".
- Status card `Biofluids Venture.md`: §1 state line and test count, §2 new findings with their numbers, §3 next actions, §5 commands, §9 history with commit hashes.
- Memory: update `biorheo-phase0-state.md` and any finding worth a new memory file; keep `MEMORY.md` pointers current.
- `~/.claude/lessons.md`: one line per corrected mistake this session, with the why and the date.
- If something is half done, say exactly what is on disk uncommitted and what the next command is.

Rules that hold throughout: no parameter or constant from memory (DOI-resolved or tagged UNVERIFIED); records under `registry/records` are generated, change the script; every number in a `.md` under `data/` or `case_studies/` is script-printed; no `si_fit`; only `Status/` in the Obsidian vault is writable; do not edit `.tex`/regex/Windows-path files with sed or perl; write patch scripts with the Write tool, not bash heredocs. Keep chat replies short; put detail in files.
