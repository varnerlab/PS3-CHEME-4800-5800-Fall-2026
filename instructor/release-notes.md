# CHEME 4800/5800 Problem Set 3

**Feeding a cell-free protein factory**

Download **Source code (zip)** under **Assets** and extract the archive.
In VS Code, open the folder containing `Project.toml`, `README.md`, and
`check_submission.jl`. Begin with `README.md`.

Complete the functions in `src/Compute.jl` and the questions in `responses.md`
for your track. `TRACK.txt` starts as `standard`; change it to `advanced` only
when you start the optional Part 3. Install the packages with
`julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'`, then run
`julia --startup-file=no --project=. check_submission.jl` to check your work.
The starter functions are unfinished, so the tests initially fail.

The release includes the starter code, the deGFP model data, the supplied LP
solver, test scripts, the plotting driver, the submission checker, and the
grading rubric. The reference solution is not included. Use Julia **1.12.7**.

**Release:** Saturday, October 3, 2026.

**Due:** Saturday, October 17, 2026, at **11:59 PM ET** on Canvas.
