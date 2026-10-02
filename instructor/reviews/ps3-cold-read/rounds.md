# PS3 cold-read and code-documentation review

Reviewed the assignment README, rubric, and all 15 active Julia/Python source
files, including the reference implementation and instructor utilities. Used
the cold-read skill for prose and the Cornell Julia style skill for source
documentation. Fresh readers received only the assigned prose and its preceding
context; they did not receive implementations, earlier reports, or editor reasoning.

The final fresh-reader passes have no remaining clarity flags. The rubric needed
no edits. All 22 named functions have docstrings, and the three student TODOs now
describe the required steps, return shapes, units, and error handling. Executable
syntax is unchanged, and the reference solution passes all 48 scored tests.

## Fresh-reader rounds

| Round | Target | Severity 3 | Severity 2 | Severity 1 | Raw report |
|---|---|---:|---:|---:|---|
| 1 | README | 0 | 2 | 0 | [Report](r1/readme-report.txt) |
| 1 | Rubric, after reading the README | 0 | 0 | 0 | [Report](r1/rubric-report.txt) |
| 2 | Revised README | 0 | 1 | 0 | [Report](r2/readme-report.txt) |
| 2 | Student-facing docstrings and TODOs, after reading the README | 0 | 0 | 0 | [Report](r2/contracts-report.txt) |
| 3 | Revised README | 0 | 0 | 1 | [Report](r3/readme-report.txt) |
| 4 | Final README | 0 | 0 | 0 | [Report](r4/readme-report.txt) |

The counts are raw reader flags, not separate votes on every earlier issue.
There were four specific flags across the rounds. Two concerned different
named-tuple return descriptions; two concerned different ambiguities in the
student editing instructions. All were repaired in place. None was dismissed
as optional detail.

## Repairs and their evidence

| Flag | Decision and repair | Evidence |
|---|---|---|
| Round 1, README R1, severity 2 | Name the editable files and keep other supplied files protected. Preserve the paragraph's job of directing student work. | Students implement Compute.jl, may load helpers through Include.jl, and answer responses.md; the rubric already specifies those files. |
| Round 1, README R2, severity 2 | State explicitly that allocation_constraints returns a named tuple with the four required fields. | The supplied solver example uses limits.lower, limits.upper, limits.A, and limits.b. |
| Round 2, README B011, severity 2 | Separate implementing functions from answering prompts within the same sentence, so responses.md cannot be read as the destination for code. | The starter functions are in src/Compute.jl; responses.md contains the three writing prompts. |
| Round 3, README B038, severity 1 | Use Julia's named-tuple notation, (; budget, equal, optimized), in the production-curve return description. | The driver and tests access each returned row by these field names. |

The source review added missing Python function docstrings, expanded abbreviated
Julia helper docstrings, documented script usage and module roles, and explained
the supplied kinetic constants and unit conversions in comments. The starter
and reference function contracts are synchronized. The TODOs remain instructions;
all three starter placeholder errors remain in the student package.

## Structure, length, and behavior

| File | Prose paragraphs before → after | Prose words before → after |
|---|---:|---:|
| README.md | 30 → 30 | 1434 → 1433 |
| RUBRIC.md | 4 → 4 | 276 → 276 |

Word counts include prose in lists and tables and exclude headings, fenced code,
displayed/inline mathematics, and link destinations. No page limit was specified.
The user explicitly requested docstring additions and clearer TODOs, so source
documentation was allowed to grow.

The original files are preserved in [before/](before/). A Markdown block inventory
was used because the skill's LaTeX structure guard does not parse this format.
Block kinds and order, paragraph boundaries, headings, list-item structure,
tables, code examples, equations, numbers, citations, and link destinations are
preserved. Only README blocks B011, B029, and B038 changed. Their original roles
remain editing instructions, allocation-function contract, and production-curve
return contract. No claim, numerical result, deadline, or model assumption changed.
The rubric is byte-for-byte unchanged.

- [Structure and length results](structure-results.json)
- [Original paragraph inventory: README](README.md.paragraph-inventory.json)
- [Original paragraph inventory: rubric](RUBRIC.md.paragraph-inventory.json)
- [Skill invariant-check output](prose-invariants.txt): no flags
- [Julia syntax comparison](julia-invariants.txt): executable syntax unchanged;
  all 17 named functions have attached docstrings, including the three reference functions
- [Python syntax comparison](python-invariants.json): executable syntax unchanged;
  all five named functions and four modules have docstrings

Anonymous test functions retain their descriptive case labels. Julia scripts
without named functions have usage and side-effect comments. No functions were
added, removed, renamed, or refactored.

## Validation and rendering

The [reference checker](reference-checker.log) passes both 24-test suites and
generates the plot and tables. The checker also handles the expected failure
cases: [starter](starter-checker.log), [partial implementation](partial-checker.log),
[syntax error](syntax-error-checker.log), [missing responses](missing-responses-checker.log),
and [reporting error](reporting-error-checker.log). These were regression checks
of the existing assignment; no new parameter study or scientific audit was added.

Pandoc rendered the [README](README.html) and [rubric](RUBRIC.html) without warnings.
The [HTML inventory](render-check.json) confirms unchanged block and MathML element
counts. Six added inline-code spans identify editable filenames and individual
return fields. A visual browser preview could not be completed: the browser's
security policy blocks local-file URLs, and the sandbox blocks a local HTTP
listener. No workaround was attempted after the browser-policy rejection.
There is no PDF or page budget for this assignment.

The [student ZIP](../../artifacts/PS3-CHEME-4800-5800-Fall-2026-student.zip) was rebuilt
with the reviewed sources. It contains 24 student files and excludes this review,
the reference solution, and generated answers.
