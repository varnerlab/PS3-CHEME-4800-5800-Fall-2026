# PS3 cold read, second pass (October 2, 2026)

Scope: the student-facing README (all four sections), `responses.md`, `RUBRIC.md`,
and the `supply_ranges` docstring and TODO in `src/Compute.jl`. The first review
([../ps3-cold-read/rounds.md](../ps3-cold-read/rounds.md)) predates Part 3, the
five written questions, and the fourth function, so none of these had been read
cold before. The originals are preserved in [before/](before/).

Each round used fresh readers who saw only the reader file: everything a student
would have read before the audited part, then the audited blocks with IDs. The
persona was a senior undergraduate or first-year graduate student in chemical
engineering who has seen FBA in lecture. FVA, γ, and the assignment's own terms
were deliberately left out of the readers' assumed knowledge. The responses reader
also saw the header rows of the four generated CSV tables.

## Flag counts (severity-3 / severity-2)

| Target | R1 | R2 | R3 | R4 | R5 |
|---|---:|---:|---:|---:|---:|
| README intro, logistics, model | 0 / 12 | 0 / 8 | 0 / 8 | 0 / 16 | 0 / 10 |
| README Parts 1–2 | 0 / 7 | 0 / 8 | 0 / 7 | 0 / 3 | 0 / 3 |
| README Part 3, FVA definition | 0 / 11 | 0 / 6 | 0 / 4 | 0 / 7 | 0 / 3 |
| README Part 3, supply_ranges to submission | 0 / 1 | 0 / 6 | 0 / 4 | 0 / 6 | 0 / 6 |
| responses.md | 1 / 13 | 0 / 10 | 0 / 8 | 0 / 4 | 0 / 5 |
| RUBRIC.md | 0 / 10 | 0 / 7 | 0 / 3 | 0 / 3 | — |
| supply_ranges docstring and TODO | 0 / 2 | 0 / 2 | 0 / 2 | 0 / 2 | 0 / 1 |
| **Total** | **1 / 56** | **0 / 47** | **0 / 36** | **0 / 41** | **0 / 28** |

The rubric was not changed after round 4, so it was not re-read in round 5. Raw
reports are in `r<N>/reports/`; reader files and block manifests are in `r<N>/`.
The round-4 rise in the intro came from detail added in rounds 2–3; round 4
removed it again.

## Main repairs

- **Severity 3 (round 1).** Question 2 asked students to compare production with
  the translation capacity to 0.001 µM/h, but no file or terminal line gave that
  value. `src/Reporting.jl` now prints `Translation capacity: 11.175812 μM/h of deGFP`,
  and the README and question 2(a) point to it. This is the only code change.
- **README Part 3** was reorganized: the scope of FVA is stated before the
  definition, the fixed-flux rule comes before the function, `supply_ranges` is
  specified with its fields as a list, the method is three numbered steps (z* in
  mM/h, the `vcat` row −cᵀv ≤ −γz*, ±1 objectives with the sign flip), and the
  no-tolerance warning gives its reason.
- **responses.md** questions now have lettered parts and name the exact file and
  columns. Question 1 has five parts and question 3 four. Question 4(b) no longer
  names glutamine and threonine, so it does not reveal the answer to 3(a).
- **README Parts 1–2**: the empty `A` (0×n) and one-row `A` / `b = [U]` shapes are
  explicit, `solve_lp` is said to maximize, the driver's files and checks are
  described, and the tuple fields of `production_curve` are spelled out.
- **README intro**: the budget is defined as a limit on the sum of supply rates,
  `load_model` usage and the 0–30 mM/h supply bounds are stated, and the archive
  name has an example.
- **RUBRIC.md** clarity only: the first review bullet was split, "preserve inputs"
  and "supplied model" were reworded, and the clean-copy consequence was stated.
  No scoring rule changed.
- The `supply_ranges` docstring (starter and reference) gained a `# Method`
  section matching the README steps, and TODO 4 was brought in line with it.
  The reference prompts in `solution/responses.md` were kept identical to the
  student prompts.

Facts added to the prose were checked against code and data: supply bounds of
0–30 mM/h (`data/conditions.tsv`), the export and degradation blocks,
`upper[translation] = translation_limit` (`src/Model.jl`), one deGFP per
initiation (`translation_deGFP` in `data/reactions.tsv`), SHA-256 fingerprints
and output regeneration only when all suites pass (`check_submission.jl`), and
unique equal-allocation supplies at U = 1 (`supply-ranges.csv`). The glutamate
deduction in question 4(b) was confirmed by re-solving with glutamate pinned at
its maximum: glutamine fell to its minimum, threonine to zero, and production
stayed at 7.4432 µM/h.

## Author decisions (applied after round 5)

1. Grading reruns the checker in a clean copy, confirms the tests ran, then reviews
   the written answers. The rubric now says the plots and tables are regenerated
   rather than taken from the submitted `outputs/`, and the README's checker
   paragraph points to this.
2. Review items count only at 72/72, and all tests passing with unanswered
   questions scores 3. The rubric now states this under the score table; the
   redundant "pending instructor review" paragraph was removed.
3. Lecture references (L6a, L6b) were removed from the README and the data notes.
4. Anyone who submits something by the initial deadline is eligible to revise,
   even if it does not run. The README's revision bullet now says so.

A round-6 reader checked the revised rubric only (0 / 6). Three flags came from
the new wording and were fixed: the rerun is now described before the score
table, `responses.md` is named as the source of the written answers, and an
unverified reason for the `MANIFEST.txt` fingerprints was removed from the
README. The author confirmed that students should still zip `outputs/` and
`MANIFEST.txt` even though grading regenerates the outputs. Still open: the
standard for "comments that explain nonobvious choices" is not defined.

## Remaining round-5 severity-2 flags

Most ask for optional detail. These would be cheap to fix but were not applied,
because every applied fix needs a fresh re-read:

- Intro: "at most 1/20" can still read as forcing equal rates; FVA is not tied to
  a strategy and budget; "Use these conditions for both strategies" reads as
  conflicting with strategies that change bounds.
- Part 3: why γ < 1 is useful is no longer stated next to γ; "Steps 2 and 3 below"
  is a forward reference; the six supply-range column names are given as a
  prefix/suffix pattern; the `push!`/`vcat` sentence (README and docstring) still
  needs a reread; the `ArgumentError` sentence is split by "as in Part 1"; the
  checker paragraph packs four facts into one sentence.
- responses.md: 1(e)'s note that equal-allocation values are unique sits uneasily
  with the README's "either strategy can have multiple" optima; 2(c) asks three
  things; 3(a) asks for two of three groups; 4(b) chains four steps; 5(a)
  compares two thresholds with different tolerances.

## Structure and validation

Links, code blocks, display equations, and headings are unchanged from
[before/](before/). The only numbers added are 30 (supply bound), 256 (SHA-256),
and 10 (10⁻⁵ in responses.md). Prose word counts: README 2178 → 2821, RUBRIC
312 → 381, responses.md 457 → 664. Pandoc renders all three files; the README's
two math warnings come from the unchanged displays and also occur in the original.
`instructor/validate.py` passes every case after the edits, and the student ZIP was
rebuilt.

## Tightening pass (after the author decisions)

The README, `responses.md`, rubric, and `supply_ranges` docstring were tightened
to remove repetition. A round-7 read of the tightened text (0 severity-3; 41
severity-2: intro 9, Parts 1–2 5, Part 3 definition 7, Part 3 implementation 10,
responses 4, rubric 4, docstring 2) showed some cuts had gone too far: signposts
("described below", "the steps below show how", how z* is found), the named
driver script, and several sentences packed too densely. Those were restored or
split, question 1(e) became 1(e)–(f), and 4(b) was loosened. The restorations
mostly bring back wording that earlier rounds had passed; they were not re-read
separately.

Final prose word counts against the pre-review originals: README 2178 → 2406
(down from 2821 before tightening), RUBRIC 312 → 400 (the grading-procedure
paragraph is new), responses.md 457 → 646 (lettered parts). Links, code blocks,
displays, and headings are unchanged; the validator passes and the student ZIP
was rebuilt.

## Rubric simplified (final)

At the author's direction, RUBRIC.md was replaced with the five-level scheme used
in grading: 0 nothing runs; 1 some tests pass but not a majority (1–36); 2 a
majority but not all (37–71); 3 all 72 pass but questions unanswered; 4 all pass
and questions answered. The review checklist (docstrings, comments, hard-coding,
units) was removed, which also closes the open "nonobvious choices" item.
