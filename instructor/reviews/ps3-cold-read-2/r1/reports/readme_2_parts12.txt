README.md:B021 | UNDEFINED | 1 | "Define $\mathcal{I}$ as the set of amino-acid supply columns" | The text never ties the math symbols to the code. It does not say that I is `model.amino_acid_uptake`, that p is `model.protein`, or that U is the `budget` argument. B030's "supplied column indices" does not name the fields either. | I = model.amino_acid_uptake, p = model.protein, U = budget, ℓ/u = model.lower/model.upper.
README.md:B021 | REREAD | 1 | "hold any additional supply inequality" | "Any" with the singular "inequality" left me unsure whether A always has rows, sometimes has none, or can have several. | A may be empty (equal) or have one row (optimized).
README.md:B023 | INCONSISTENT | 1 | "u_j^{\mathrm{equal}}=\min(u_j,U/K)" | In B020, **u** is the bound vector the LP actually uses. Here u_j means the model's original upper bound, so one symbol stands for two things. | u_j = the model's original upper bound before the budget is applied.
README.md:B024 | FORWARD | 1 | "so return an `A` with zero rows and an empty `b`" | "Return" from what? No function has been introduced yet; allocation_constraints first appears in B029. | allocation_constraints(...; strategy=:equal) returns A with no rows and an empty b.
README.md:B024 | UNDEFINED | 2 | "return an `A` with zero rows and an empty `b`" | The shape and type of the empty A are not given: 0×(number of reactions) or 0×0? Matrix{Float64}? I can't tell what solve_lp or the tests will accept. | A = zeros(0, number_of_reactions), b = Float64[].
README.md:B027 | REREAD | 2 | "One row of `A` has a coefficient of one at each" | "One row of A" reads like one of several rows, which is hard to square with B025's "add one shared supply constraint". It is also not stated that A must be a 1×n matrix rather than a vector, or that b is a one-element vector rather than the scalar U. | A is exactly one 1×n row; b = [U].
README.md:B029 | REREAD | 1 | "so that solving one scenario cannot change the next" | "Scenario" is undefined. The risk is mutating model.upper while building the constraints, not while solving, so the stated reason confused me. | Don't modify model.lower/upper in place; return copies so later (budget, strategy) calls start from the original bounds.
README.md:B029 | OVERLOAD | 1 | "Copy the bound arrays so that solving one scenario cannot change the next." | One list item packs the signature, return fields, allowed strategies, copying rule and error rules into four sentences. | Split into: returns, strategies, copy rule, error rule.
README.md:B030 | FORWARD | 1 | "Your functions must also work on the small test networks" | This is the first mention of small test networks; nothing earlier said the tests use models other than the supplied one. | The test scripts build toy models with the same named-tuple fields.
README.md:B030 | PURPOSE | 1 | "at least one supply column" | The reason for this guarantee is not given. I worked out afterwards that it means K ≥ 1, so U/K never divides by zero. | You need not handle K = 0.
README.md:B030 | REREAD | 1 | "Keep all lower bounds and unrelated upper bounds unchanged." | "Unrelated" depends on the strategy: under equal it means non-supply bounds, and under optimized it means all bounds. | Only the supply upper bounds change, and only under :equal.
README.md:B033 | ORDER | 1 | "Explain these numbers to yourself before testing the larger network." | The hand check comes after the command to run the Part 1 tests. "Larger network" is ambiguous: the full model, or the test networks just called "small"? | Do the hand check before running testme_part_1.jl; "larger network" means the full 146×265 model.
README.md:B036 | UNDEFINED | 2 | "result = solve_lp(model.S, limits.lower, limits.upper, c;" | Nothing says whether solve_lp maximizes or minimizes c'v. Many LP solvers minimize by default, so I would wonder whether to negate c. | solve_lp maximizes c'v, matching B020.
README.md:B037 | ORDER | 1 | "In that case the solver returns `nothing` for `flux` and `objective`." | The fact that motivates the check arrives after the instruction to throw. | State first that non-optimal solves return nothing, then say to check and throw.
README.md:B037 | UNDEFINED | 1 | "throw an `ErrorException` naming the budget, strategy, and status" | The message format and the type of `status` are not given, and I can't tell whether the tests check the message text. | error("... budget=$budget strategy=$strategy status=$(result.status)").
README.md:B037 | INCONSISTENT | 1 | "Check `result.optimal` before reading the fluxes." | production_curve needs the objective value, not the fluxes, so "fluxes" left me unsure which field to read. | Check optimal before reading flux or objective; production = objective (or flux[p]) × 1000.
README.md:B038 | REREAD | 1 | "one named tuple `(; budget, equal, optimized)` per input budget" | The `(; a, b)` shorthand is not basic Julia for many students. It is also not stated that `equal`/`optimized` hold the production rates. | Each element is (budget=…, equal=µM/h rate, optimized=µM/h rate).
README.md:B038 | UNDEFINED | 1 | "An empty input returns an empty vector. Negative or nonfinite budgets raise" | The element type of the empty vector is not given, nor whether all budgets must be validated before any solve or only when each is reached. | Validate all budgets up front; the element type is unspecified.
README.md:B041 | REFERENT | 1 | "The supplied driver evaluates budgets" | runproduction.jl has never been called "the driver"; the connection is left to inference. | Driver = runproduction.jl.
README.md:B041 | FORWARD | 2 | "Until `supply_ranges` is complete, it skips the Part 3 files" | supply_ranges has never been introduced; I don't know what it is or that it is a function I must write. | supply_ranges is the Part 3 function in src/Compute.jl (presumably the fourth of the four).
README.md:B041 | REFERENT | 2 | "it skips the Part 3 files and says so:" | The colon after the skipping sentence makes the list below look like the skipped Part 3 files, not the files that get written. I also can't tell which output, if any, belongs to Part 3 (amino-acid-allocation.csv?). | The list is the Part 2 outputs; Part 3 files are separate and not listed here.
README.md:B042 | UNDEFINED | 2 | "with independent feasibility checks reported in the terminal" | "Independent feasibility checks" is never explained: independent of what, what is checked, and what should I look for? | The driver re-checks S v ≈ 0, the bounds and the budget for the reported fluxes and prints pass/fail.
README.md:B042 | PURPOSE | 1 | "and the fixed translation-capacity limit" | The limit is on translation initiation in mM/h (B014). It is not explained how it compares to protein production in µM/h, or why it appears on the production plot. | A horizontal line in µM/h showing the maximum production that translation capacity allows.
README.md:B043 | REFERENT | 2 | "use it with the tables to answer questions 1 and 2" | "The tables" was never introduced; B042 lists CSV files, not tables. | The two CSV files in outputs/.
README.md:B043 | PURPOSE | 1 | "Tests check feasibility and production, so they accept alternative optimal allocations." | production_curve returns only production rates, so I can't tell which allocation the tests could be checking for feasibility. | The tests and the driver check the supply fluxes behind each solve, not just the returned numbers.

README.md:B019: Introduces the flux vector v and p as the index of the protein-output column.
README.md:B020: The LP: maximize protein output subject to steady state, flux bounds and an extra inequality Av ≤ b.
README.md:B021: ℓ/u are the flux bounds, A/b are the extra supply constraint, I is the set of supply columns, K is their count and U is the total budget in mM/h.
README.md:B022: Under equal allocation each amino acid gets at most one equal share of the budget.
README.md:B023: The equal-allocation upper bound is the smaller of the original bound and U/K.
README.md:B024: Unused shares are not shared out, and equal allocation needs no extra rows in A or b (though A's shape is unclear).
README.md:B025: Optimized allocation keeps the original bounds and adds one shared budget constraint.
README.md:B026: The sum of all amino-acid supply fluxes is at most U.
README.md:B027: How to build the optimized row of A and its b entry; metabolic synthesis does not count against the budget.
README.md:B028: Two Part 1 functions go in src/Compute.jl.
README.md:B029: Specifications of protein_objective and allocation_constraints, including copying and error rules.
README.md:B030: Use the model's index fields rather than hard-coded positions, handle the varied test networks, and change only the supply upper bounds.
README.md:B031: Introduces the Part 1 test command.
README.md:B033: A two-amino-acid toy example with expected answers of 1.5 (equal) and 2 (optimized).
README.md:B035: production_curve solves both strategies for each budget using the Part 1 functions and solve_lp.
README.md:B037: The fields of the solve_lp result; throw an ErrorException on a non-optimal solve.
README.md:B038: The return format of production_curve, units, ordering, empty-input behaviour and error rules.
README.md:B039: Introduces the Part 2 commands.
README.md:B041: UNCLEAR - It gives the driver's budget grid, but refers to an unintroduced supply_ranges and ambiguous "Part 3 files", and the colon makes the following list look like the skipped files.
README.md:B042: Lists the plot and CSV outputs; "independent feasibility checks" is unexplained.
README.md:B043: Use the plot and outputs to answer Q1–Q2; optimized allocations are not unique and tests allow for that; Part 3 quantifies the spread.

TOTALS: severity-3 = 0; severity-2 = 7; severity-1 = 18
