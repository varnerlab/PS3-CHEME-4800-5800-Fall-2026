# PS3 reference solution

This directory holds the reference solution for PS3: Feeding a cell-free
protein factory.

- [src/Compute.jl](src/Compute.jl) completes the four functions:
  `protein_objective`, `allocation_constraints`, and `production_curve` for the
  Standard track, and `supply_ranges` for the Advanced track.
- [responses.md](responses.md) has model answers to the three written questions.

Use the reference solution to understand mistakes and debug your own work, but
do not copy it. Copying the reference solution results in a score of `0` for
this assignment.

## Running the reference solution

The easiest way is to download `PS3-CHEME-4800-5800-Fall-2026-solution.zip`
from the repository's solution release. That ZIP already has this `Compute.jl`
in `src/` and the model answers in `responses.md`; see its README.

To run it in your own PS3 folder instead, copy this `solution` directory into
that folder, next to `src/`. Then run this command from your PS3 folder:

```bash
julia --startup-file=no --project=. check_submission.jl --solution
```

The `--solution` flag loads `solution/src/Compute.jl` and `solution/responses.md`
in place of your files and always checks the Advanced track. All 72 tests should
pass, and the status should read `REFERENCE SOLUTION PASSED`. The plots and
tables go to `solution/outputs/`; your own code, responses, and outputs are not
changed.

| Budget (mM/h) | Equal allocation (µM/h) | Optimized allocation (µM/h) |
|---:|---:|---:|
| 0.0 | 2.661 | 2.661 |
| 1.0 | 6.004 | 7.443 |
| 2.0 | 8.927 | 11.176 |
| 4.0 | 11.176 | 11.176 |
