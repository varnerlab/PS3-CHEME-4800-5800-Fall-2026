R1 | INCONSISTENT | 2 | "Keep the data, supplied code, and test files unchanged." | This appears to prohibit editing Compute.jl and Include.jl, although this same paragraph requires implementing supplied functions and loading helpers from Include.jl. The intended exceptions are not stated. | Best guess: edit the three functions in Compute.jl and add helper-loading statements to Include.jl as needed; leave all other supplied code unchanged.
R2 | FORWARD | 2 | "returns `(lower, upper, A, b)`" | This looks like an ordinary Julia tuple return specification, but the Part 2 example accesses limits.lower and limits.A, which require named fields. I would discover the required return structure only when reaching that example. | Best guess: return a named tuple with fields lower, upper, A, and b.

Paragraph roles:
B002: Explain what a cell-free protein synthesis system needs and how it obtains amino acids.
B003: State the prediction task, introduce the two supply strategies, and identify the final plot.
B004: Connect the assignment to prior coursework and identify the model's sources.
B008: Specify Julia and the working directory for setup and testing.
B010: Explain the setup commands and why the initial tests fail.
B011: Identify the student work, helper-file location, and protected materials.
B013: Describe the model's size, processes, loader, and returned fields.
B015: Set the common operating conditions for comparing strategies.
B016: Explain supply-flux signs, protein output, and steady-state balances.
B017: Define flux units, output conversion, and the interpretation of a supply budget.
B019: Introduce the flux vector and protein-output index for the optimization problem.
B021: Define bounds, additional supply constraints, supply-column indices, and the total budget.
B022: Introduce the equal-share upper-bound rule.
B024: Explain unused equal shares and why no extra total-budget inequality is needed.
B025: Introduce the shared constraint while retaining original individual bounds.
B027: Explain how to encode optimized allocation and distinguish external supply from synthesis.
B028: Introduce the two Part 1 implementation tasks.
B030: Require index-based, general implementations that preserve unrelated bounds.
B031: Introduce the command for Part 1 checks.
B033: Provide a small numerical example for checking understanding of the strategies.
B035: Introduce production_curve and demonstrate use of the supplied solver.
B037: Explain solver-result fields and required handling of unsuccessful solves.
B038: Specify production_curve's output format, units, ordering, and input edge cases.
B039: Introduce the Part 2 test and production-driver commands.
B041: Specify the driver's budget grid and introduce its output files.
B043: Connect generated results to written responses and explain acceptable alternative optima.
B045: Introduce the final complete-check command.
B047: Explain the checker's tests, manifest, warnings, generated outputs, and instructor review.
B048: Specify archive contents, naming, and the student's Canvas upload responsibility.
B049: Explain what to submit if tests still fail at the deadline.
TOTALS: severity-3 = 0; severity-2 = 2; severity-1 = 0
