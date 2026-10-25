This release contains the completed PS3 reference solution for CHEME 4800/5800,
Fall 2026: **Feeding a cell-free protein factory**.

Download **PS3-CHEME-4800-5800-Fall-2026-solution.zip** from the release assets
and extract it. This prepared ZIP has the completed code in `src/`, model
answers in `responses.md`, the original starter in `starter/`, and the
assignment, data, tests, and grading rubric. GitHub's automatic source archives
do not include the solution; use the named solution ZIP.

With Julia `1.12.7`, open a terminal in the extracted folder and run:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

All 72 tests should pass. The test scripts and supplied code were rewritten to
be easier to read; they run the same checks and give the same results as the
versions in the assignment release.

Use the reference to understand and debug your own work; you may not copy it.
The assignment's independent-work and revision policies still apply.
