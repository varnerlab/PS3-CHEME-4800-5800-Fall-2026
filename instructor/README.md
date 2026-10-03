# PS3 instructor notes

The student task is fixed: one protein (deGFP), equal versus optimized allocation
of a total amino-acid supply budget and a production-versus-budget plot.
Parts 1–2 require three functions and two written responses, graded out of four
points using 48 tests. Optional Advanced work adds flux variability analysis
(FVA) of the 20 supply fluxes: `supply_ranges`, 24 Part 3 tests, and question 3.
Complete Advanced work earns one Magic Point after a regular score of 4 and
teaching-team review; unfinished Advanced work does not reduce the regular score.
The package uses one source file and one response file for both tracks;
`TRACK.txt`, shipped as `standard`, selects which checks and reports run. The
supplied LP solver uses JuMP and GLPK. It follows L6a's urea-cycle interface:
`build(MyPrimalFluxBalanceAnalysisCalculationModel, data)` and `solve(calculation)`
with the same five model fields and dictionary result keys. The standalone
package supplies the type and builder locally. The `A` and `b` keywords add
PS3's inequalities; nonoptimal solves return the status with `nothing` values
so student functions can identify the failing budget and strategy. Part 3 cites Mahadevan and Schilling
(2003), Metabolic Engineering 5(4):264–276,
[DOI: 10.1016/j.ymben.2003.09.002](https://doi.org/10.1016/j.ymben.2003.09.002).

## Reviewing the assignment

- Read [the problem statement](../README.md) and [rubric](../RUBRIC.md).
- Inspect [the reference plot](reference/production-versus-budget.png) and its
  [numerical table](reference/production-versus-budget.csv).
- Inspect the Part 3 [supply ranges](reference/supply-ranges.png),
  [range table](reference/supply-ranges.csv),
  [range-width plot](reference/range-width-versus-budget.png), and
  [range-width table](reference/range-width-versus-budget.csv).
- The local, Git-ignored [reference solution](../solution/README.md) contains
  completed code and sample written responses.
- [The student ZIP](artifacts/PS3-CHEME-4800-5800-Fall-2026-student.zip)
  contains starter code and data, with no instructor files, outputs, or solutions.

Release is Saturday, October 3, 2026. The initial submission deadline is
11:59 PM ET on Saturday, October 17, 2026. No commit, push, tag, GitHub release, or Canvas upload is
performed by these scripts.

## Reproducing the checks

To run the reference solution directly from the assignment root:

```bash
julia --startup-file=no --project=. check_submission.jl --solution
```

The reference solution always uses the Advanced track, whatever `TRACK.txt`
says. This runs all 72 tests, checks `solution/responses.md`, and writes plots and
tables to `solution/outputs/` and a manifest to `solution/MANIFEST.txt`.
The student files and outputs stay unchanged. The same flag works with
`testme_part_1.jl`, `testme_part_2.jl`, `testme_part_3.jl`, and `runproduction.jl`;
`runproduction.jl --solution` writes all reference reports.
The default remains the student implementation; a missing local solution causes
an error rather than falling back to student code.

From the assignment root, install the Julia dependencies and run:

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
python3 instructor/validate.py
python3 instructor/build_student_package.py
```

The validator first checks the supplied type, builder, and solver interface with
19 assertions in [test_solver.jl](test_solver.jl), separate from student scores.
It then copies the student files into temporary directories and substitutes
the reference implementation only in those copies. It checks all 72 tests, the
plot driver, and the submission checker's handling of unfinished code, partial
implementations, source syntax errors, missing responses, and reporting errors.
It also checks that a submission with only Parts 1 and 2 and questions 1–2
complete passes the checker on the Standard track. On the Advanced track,
Part 3 test, writing, and output failures are reported separately. On the
Standard track, the driver writes only production outputs. An invalid
`TRACK.txt` stops the checker and driver with a clear message.
Logs are in `instructor/validation-output/`; reference plots and tables are copied
to `instructor/reference/`. Both directories are ignored by Git and excluded from
the student package.

An independent check reconstructs the original LP from the publication files
and compares both strategies at every plotted budget using SciPy/HiGHS. It also
recomputes the FVA endpoints at U=1 mM/h and the total range widths and
fixed-supply counts at every plotted budget:

```bash
python3 instructor/crosscheck.py /path/to/Sequence-Specific-FBA-CFPS-Publication-Code
```

That check requires NumPy and SciPy and the source revision recorded in
[data/provenance.json](../data/provenance.json). Students need only Julia.
The source repository is read only. To regenerate the teaching data from that
pinned checkout, run `python3 instructor/import_model.py /path/to/checkout`.

## Modeling choices

The stoichiometric matrix and original kinetic formulas are preserved. Glucose
supply is reduced to 1 mM/h so the amino-acid budget has a visible effect before
the fixed translation-capacity bound becomes active. The original Case 1 uses a
30 mM/h glucose limit, which makes the proposed budget experiment less informative.

Equal allocation means equal *upper bounds*, not prescribed consumption. The
optimized strategy keeps the original individual limits and adds one shared
inequality. It can leave some amino-acid supply fluxes at zero while the network
synthesizes those amino acids internally. Different optimal uptake vectors may
exist; grading must not require one particular vector.

This is a teaching adaptation of a static model, not a new experimental
calibration. The original deGFP kinetic transcript-length parameter is 683,
whereas the transcription reaction consumes 678 nucleotides. The source also
lists a protein-length parameter of 229, while its translation reaction consumes
225 amino acids. The supplied kinetic formula uses the transcript-length
parameter; the unused protein-length parameter is not introduced into PS3.
These original conventions are retained, so no student task infers length from
the coefficients or modifies the sequence. The data notes disclose the kinetic
transcript-length distinction.

FVA ranges are unique even when optimal flux vectors are not, so Part 3 tests
compare endpoints with reference values. The protein requirement is
`c'v >= fraction*z` with no tolerance. GLPK solves every Part 3 scenario with
the exact requirement. Relaxing it by 1e-9 mM/h creates false ranges as wide as
1.6e-5 mM/h (optimized allocation at U=0.4), so the README forbids a tolerance and defines a fixed flux as one whose
width is at most 1e-5 mM/h. At the exact requirement, every computed width is
either below 1e-12 or above 1e-4 mM/h. FVA is restricted to the supply fluxes:
at U=1 with optimized allocation, 66 of the 265 columns span at least 29 mM/h because separate forward and
reverse columns can carry offsetting flux.

At U=1 with optimized allocation, 17 of 20 supply fluxes are fixed. Glutamate,
glutamine, and threonine share one 0.05472 mM/h portion of the budget; the
allocation table's threonine value is one vertex of that set. The supply ranges
open at U=2.0 (optimized) and U=2.9 (equal), the sampled budgets at which each
strategy reaches translation capacity. At 99% of maximum output, no supply flux
is fixed. The optional third written response asks students to identify one
variable supply and explain why one optimal flux vector is not unique.

The shared budget is molar supply per time and volume, with equal weights per
amino acid. It does not represent economic cost or initial reagent inventory.
The LP predicts a feasible steady-state production rate; final batch titer needs
assumptions about how long that rate persists.

## Preparing a release

Review the dates in the main README, the plot, and student workload, then run
the validation and build scripts. The builder refuses a
release package with dates still marked to be announced. It uses an explicit
allowlist and checks that the four starter errors remain in `src/Compute.jl`.

The reference code, sample answers, and generated reference outputs are kept
outside the student archive. Preserve the original attribution in
[THIRD_PARTY_NOTICES.txt](../THIRD_PARTY_NOTICES.txt). The source and data can be
released under the included permission notices when publication is authorized.
