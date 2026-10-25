# PS3 reference solution: Feeding a cell-free protein factory

This archive contains the completed reference solution for CHEME 4800/5800,
Fall 2026. The solution is already in [src/](src/), and [TRACK.txt](TRACK.txt)
is set to `advanced`, so the checker runs all three parts without moving or
editing any files.

Use Julia `1.12.7`. From the folder containing this README, run:

```sh
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

The first command installs the packages and takes several minutes. The checker
should report that all 72 tests pass and that the status is `READY TO PACKAGE`.
It writes the plots and tables to `outputs/` and a local `MANIFEST.txt`; it
does not upload anything.

The archive includes:

- [src/Compute.jl](src/Compute.jl): the completed, commented reference code for
  `protein_objective`, `allocation_constraints`, `production_curve`, and
  `supply_ranges`. The other files in `src/` are the supplied code.
- [responses.md](responses.md): model answers to the three written questions.
- [starter/Compute.jl](starter/Compute.jl): the original starter, for comparison.
- [ASSIGNMENT.md](ASSIGNMENT.md): the original problem statement. Its
  directions to complete the code apply to your own submission.
- [data/](data/README.md): the deGFP model data.
- [RUBRIC.md](RUBRIC.md): grading rules and completion requirements.
- [testme_part_1.jl](testme_part_1.jl), [testme_part_2.jl](testme_part_2.jl),
  and [testme_part_3.jl](testme_part_3.jl): the public tests.

To write only the plots and tables, run
`julia --startup-file=no --project=. runproduction.jl`.

| Budget (mM/h) | Equal allocation (µM/h) | Optimized allocation (µM/h) |
|---:|---:|---:|
| 0.0 | 2.661 | 2.661 |
| 1.0 | 6.004 | 7.443 |
| 2.0 | 8.927 | 11.176 |
| 4.0 | 11.176 | 11.176 |

Use this reference to understand mistakes and debug your own work. You may
study it, but you may not copy it. Copying the reference solution results in a
score of `0` for this assignment.
