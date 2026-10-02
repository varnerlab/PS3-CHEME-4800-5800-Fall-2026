ALREADY READ: Assignment README
[README.md:B001 heading]
# Problem Set 3 (PS3): Feeding a cell-free protein factory

[README.md:B002 paragraph]
A cell-free protein synthesis system uses cellular machinery in a reaction
mixture to make a protein. That machinery needs amino acids, nucleotides, and
energy. The metabolic reactions in the mixture can make amino acids from other
available compounds, and amino acids can also be supplied directly.

[README.md:B003 paragraph]
You will use linear programming to predict the production rate of deGFP, a green
fluorescent protein, under two amino-acid supply strategies. With **equal
allocation**, every amino acid receives the same supply limit. With **optimized
allocation**, the solver distributes a shared supply budget among the amino acids.
Your final result is a plot of protein production against the total supply budget.

[README.md:B004 paragraph]
This assignment extends the flux balance analysis in L6b. The supplied model is
adapted from the model described by **Vilkhovoy et al. (2018), “Sequence Specific
Modeling of E. coli Cell-Free Protein Synthesis,” ACS Synthetic Biology, 7(8),
1844–1857** ([paper](https://pubmed.ncbi.nlm.nih.gov/29944340/),
[DOI: 10.1021/acssynbio.7b00465](https://doi.org/10.1021/acssynbio.7b00465)).
The code and model data are available in the Varner laboratory's
[sequence-specific cell-free protein synthesis repository](https://github.com/varnerlab/Sequence-Specific-FBA-CFPS-Publication-Code).
The [data notes](data/README.md) identify the source and the PS3 operating conditions.

[README.md:B005 heading]
## Logistics

[README.md:B006 list]
- **Release:** Saturday, October 3, 2026. **Initial submission deadline:**
  **11:59 PM ET on Saturday, October 17, 2026**.
- **Submission:** Upload `CHEME-4800-5800-PS3-<your netid>.zip` to Canvas.
  Replace the complete placeholder with your NetID.
- **Infinite-revision policy:** Submit something by the initial deadline to
  participate. Eligible students may then revise and resubmit until the end of
  the semester; the highest score is retained. No initial submission means a
  score of `0` and no access to later revisions.
- **Reference solution:** Released after the initial deadline. You may consult
  it to understand mistakes and debug your work, but may not copy it. Copying
  the reference solution results in a score of `0`.
- **Group and AI policy:** Submit independent work. You may discuss ideas with
  classmates but may not share code or solutions directly. You may use Julia
  documentation, AI tools, and internet resources. You are responsible for
  understanding and explaining your implementation.
- **Grading:** Maximum `4` points; see [RUBRIC.md](RUBRIC.md).

[README.md:B007 heading]
## Getting started

[README.md:B008 paragraph]
Use Julia `1.12.7`. The assignment contains its own data and source files.
Open the extracted assignment folder in VS Code. Run these commands from the
folder containing this README and [Project.toml](Project.toml):

[README.md:B009 code]
```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B010 paragraph]
The first command installs the required packages and can take several minutes.
The second runs the assignment tests. The three functions in
[src/Compute.jl](src/Compute.jl) initially raise errors because you must implement
them. Tests will fail until that work is complete.

[README.md:B011 paragraph]
Complete the functions and prompts in [responses.md](responses.md).
Keep helper files in `src/`; edit [Include.jl](Include.jl) to load them.
You may edit `src/Compute.jl`, `Include.jl`, and `responses.md` and add helpers.
Keep all other supplied files unchanged.

[README.md:B012 heading]
## The model and its operating conditions

[README.md:B013 paragraph]
The model contains 146 species and 265 reactions. It includes metabolism,
transcription, translation, amino-acid supply, and protein output. The supplied
[load_model function](src/Model.jl) reads the stoichiometric matrix and sets the
operating conditions. Its returned named tuple contains:

[README.md:B014 table]
| Field | Meaning |
|---|---|
| `S` | Stoichiometric matrix; rows are species and columns are reactions. |
| `reactions`, `metabolites` | Identifiers in matrix column and row order. |
| `lower`, `upper` | One lower and upper flux bound per reaction, in mM/h. |
| `amino_acids` | Names of the 20 supplied amino acids. |
| `amino_acid_uptake` | Column indices of their supply reactions. |
| `protein` | Column index of the deGFP output reaction. |
| `translation` | Column index of translation initiation. |
| `translation_limit` | Fixed translation-capacity bound, in mM/h. |

[README.md:B015 paragraph]
The glucose supply limit is **1 mM/h** and the oxygen supply limit is
**100 mM/h**. Amino-acid synthesis is enabled; the model's designated
amino-acid degradation reactions and amino-acid removal reactions are blocked.
Transcription and translation bounds are calculated from fixed expression
parameters before solving. Use the supplied conditions for both strategies.

[README.md:B016 paragraph]
Supply reactions are written as `[] --> species` and have **nonnegative**
fluxes. Each amino-acid supply flux therefore contributes positively to the
budget. The deGFP output reaction removes protein from the balanced network and
records its production rate. Metabolic intermediates satisfy steady-state
balances while protein leaves through that output reaction.

[README.md:B017 paragraph]
All model fluxes are in **mM/h**, millimoles per liter per hour. Report protein
production in **µM/h** by multiplying its model flux by 1000. The supply budget
is a rate limit; it does not specify an initial concentration or a monetary cost.
This static model predicts feasible production rates under fixed conditions.

[README.md:B018 heading]
## Part 1: Construct the objective and supply constraints

[README.md:B019 paragraph]
Let the vector $\mathbf{v}$ contain the reaction fluxes, and let $p$ be the
protein-output column. The production problem is given by:

[README.md:B020 math]
$$
\begin{aligned}
\underset{\mathbf{v}}{\operatorname{maximize}}\quad & v_p\\
\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\\
& \boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\\
& \mathbf{A}\mathbf{v}\leq\mathbf{b}.
\end{aligned}
$$

[README.md:B021 paragraph]
The vectors $\boldsymbol{\ell}$ and $\mathbf{u}$ specify the flux bounds.
The matrix $\mathbf{A}$ and vector $\mathbf{b}$ hold any additional supply
inequality. Define $\mathcal{I}$ as the set of amino-acid supply columns, let
$K=|\mathcal{I}|$, and let $U\geq0$ be the total supply budget in mM/h.

[README.md:B022 paragraph]
**Equal allocation.** Each amino acid may use at most one equal share of the
budget. Its upper bound is given by:

[README.md:B023 math]
$$
u_j^{\mathrm{equal}}=\min(u_j,U/K),\qquad j\in\mathcal{I}.
$$

[README.md:B024 paragraph]
The supply fluxes can be smaller than their limits. Any unused share remains
unavailable to the other amino acids. These individual limits already enforce
the total budget, so return an `A` with zero rows and an empty `b`.

[README.md:B025 paragraph]
**Optimized allocation.** Retain the original individual bounds and add one
shared supply constraint, given by:

[README.md:B026 math]
$$
\sum_{j\in\mathcal{I}}v_j\leq U.
$$

[README.md:B027 paragraph]
The solver may distribute the available supply unequally. One row of `A` has
a coefficient of one at each amino-acid supply column and zero elsewhere;
the corresponding entry of `b` is `U`. Metabolic synthesis of amino acids remains
available under both strategies and does not count as external supply.

[README.md:B028 paragraph]
Implement two functions in [src/Compute.jl](src/Compute.jl):

[README.md:B029 list]
1. **protein_objective(model)** returns a `Vector{Float64}` with one coefficient
   per reaction. Its protein-output coefficient is one; the others are zero.
2. **allocation_constraints(model, budget; strategy=:equal)** returns a named tuple
   with fields `lower`, `upper`, `A`, and `b`. Support `:equal` and `:optimized`.
   Copy the bound arrays so that solving one scenario cannot change the next.
   Reject a negative or nonfinite budget or an unknown strategy with `ArgumentError`.

[README.md:B030 paragraph]
Use the model's supplied column indices. Your functions must also work on the
small test networks, which have different numbers and orders of reactions and
amino-acid supplies. All test models have nonnegative supply fluxes and at least
one supply column. Keep all lower bounds and unrelated upper bounds unchanged.

[README.md:B031 paragraph]
Run the Part 1 checks with:

[README.md:B032 code]
```bash
julia --startup-file=no --project=. testme_part_1.jl
```

[README.md:B033 paragraph]
For a hand check, consider two amino acids with supply fluxes $v_1$ and $v_2$.
One unit of protein output $v_p$ consumes two units of the first amino acid and
one unit of the second. The balances are $v_1=2v_p$ and $v_2=v_p$.
With a total supply budget of 6 and all other bounds loose, equal allocation
allows at most **1.5** units of protein output; optimized allocation allows
**2**. Explain these numbers to yourself before testing the larger network.

[README.md:B034 heading]
## Part 2: Compute and interpret the production curves

[README.md:B035 paragraph]
Implement **production_curve(model, budgets)** in
[src/Compute.jl](src/Compute.jl). For each budget, construct both sets of
constraints and solve both LPs. Use your Part 1 functions and the supplied
[solve_lp function](src/Solver.jl). A single solve has this form:

[README.md:B036 code]
```julia
limits = allocation_constraints(model, budget; strategy=:equal);
c = protein_objective(model);
result = solve_lp(model.S, limits.lower, limits.upper, c;
    A=limits.A, b=limits.b);
```

[README.md:B037 paragraph]
The result contains `status`, `optimal`, `flux`, and `objective`. Check
`result.optimal` before reading the fluxes. If a solve is not optimal, throw an
`ErrorException` naming the budget, strategy, and status. In that case the solver
returns `nothing` for `flux` and `objective`.

[README.md:B038 paragraph]
Return a vector containing one named tuple `(budget, equal, optimized)` per
input budget. The budget remains in mM/h; both production rates must be in µM/h.
Preserve the input order, including repeated budgets. An empty input returns an
empty vector. Negative or nonfinite budgets raise `ArgumentError`.

[README.md:B039 paragraph]
After completing the function, run:

[README.md:B040 code]
```bash
julia --startup-file=no --project=. testme_part_2.jl
julia --startup-file=no --project=. runproduction.jl
```

[README.md:B041 paragraph]
The supplied driver evaluates budgets from **0 to 4 mM/h in steps of 0.1** and
writes these files:

[README.md:B042 list]
- `outputs/production-versus-budget.png` and `.svg`: the two production curves
  and the fixed translation-capacity limit.
- `outputs/production-versus-budget.csv`: the numerical values used in the plot.
- `outputs/amino-acid-allocation.csv`: actual supply fluxes for both strategies
  at a total budget of **1 mM/h**, with independent feasibility checks reported
  in the terminal.

[README.md:B043 paragraph]
Open the plot and use it with the tables to answer [responses.md](responses.md).
The optimized allocation can have multiple equally productive flux distributions.
Tests check feasibility and production, so they accept alternative optimal allocations.

[README.md:B044 heading]
## Checking and submitting your work

[README.md:B045 paragraph]
Run the complete check after finishing the code and responses:

[README.md:B046 code]
```bash
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B047 paragraph]
The checker runs **24 tests per part, 48 total**, records their results and source
file fingerprints in `MANIFEST.txt`, and warns about unanswered discussion prompts.
If both test suites pass, it generates the plot and data tables from your code.
The instructor still reviews the implementation and written explanations.

[README.md:B048 paragraph]
Zip the entire assignment folder, including your source, data, responses,
`MANIFEST.txt`, and generated files in `outputs/`. Name the archive
`CHEME-4800-5800-PS3-<your netid>.zip` and upload it to Canvas. The checker runs
locally; you must upload the ZIP yourself.

[README.md:B049 paragraph]
If some tests still fail at the deadline, submit your current work to receive
partial credit and participate in the infinite-revision policy.

AUDIT THIS PART: Student-facing function contracts and TODOs


[src/Compute.jl:D01]
protein_objective(model) -> Vector{Float64}

Build the objective for maximizing protein output in mM/h.

# Arguments
- `model`: a named tuple with a stoichiometric matrix `S` and the integer
  column index `protein` of the protein-output reaction.

# Returns
- A new `Vector{Float64}` with one entry per column of `model.S`. The entry at
  `model.protein` is one; every other entry is zero. Leave the model unchanged.

[src/Compute.jl:D02]
allocation_constraints(model, budget::Real; strategy::Symbol=:equal) -> NamedTuple

Construct the bounds and supply constraints for one allocation strategy.

# Arguments
- `model`: a named tuple with `S`, bound vectors `lower` and `upper`, and the
  distinct supply-column indices `amino_acid_uptake`. Supply fluxes are
  nonnegative; there is at least one supply column.
- `budget::Real`: a finite, nonnegative total amino-acid supply limit in mM/h.
- `strategy::Symbol`: `:equal` (default) or `:optimized`.

# Returns
A named tuple with fields `lower`, `upper`, `A`, and `b`, suitable for `solve_lp`.
Let K be the number of supply columns and n the number of columns in `model.S`.
Copy both bound vectors; leave the model unchanged.

For `:equal`, set each supply upper bound to `min(original_upper, budget/K)`.
Return `A=zeros(0, n)` and `b=Float64[]`; the individual limits enforce the budget.
For `:optimized`, keep the model bounds and add one inequality limiting the
sum of the K supply fluxes to `budget`. Return `A` of size (1,n), with ones in
the supply columns and zeros elsewhere, and a length-one vector `b=[budget]`.
Keep all lower bounds and all other upper bounds unchanged.

# Errors
Throw `ArgumentError` for a negative or nonfinite budget (`Inf`, `-Inf`, or
`NaN`), or a strategy other than `:equal` or `:optimized`.

[src/Compute.jl:D03]
production_curve(model, budgets::AbstractVector{<:Real}) -> Vector{NamedTuple}

Compute protein production under both allocation strategies at each budget.

# Arguments
- `model`: a model accepted by `protein_objective` and `allocation_constraints`.
- `budgets`: a vector of finite, nonnegative supply limits in mM/h. Its entries
  may be unsorted or repeated; an empty vector is allowed.

# Returns
A vector of named tuples with fields `(budget, equal, optimized)`, one per input
entry in the original order, including duplicates. The budget is in mM/h;
`equal` and `optimized` are protein production rates in micromolar/h (multiply
the protein flux in mM/h by 1000). An empty input returns an empty vector.
Leave the model and budget vector unchanged.

Use `protein_objective`, `allocation_constraints`, and the supplied `solve_lp`
to solve both strategies independently at each budget.

# Errors
Throw `ArgumentError` for any negative or nonfinite budget. If either solve is
not optimal, throw `ErrorException` with the budget, strategy, and solver status
in its message. Check `result.optimal` before reading a flux or objective value.

[src/Compute.jl:TODO1]
    # TODO 1: Build and return the objective vector.
    # Use the number of columns in model.S for its length and model.protein
    # for the nonzero entry; the tests also use smaller, reordered networks.


[src/Compute.jl:TODO2]
    # TODO 2: Build and return the named tuple (lower, upper, A, b).
    # Check the budget and strategy, then copy the model's bound vectors.
    # For :equal, change only the supply upper bounds using the actual number
    # of supply columns. For :optimized, build the single shared-budget row.
    # Return arrays with the dimensions specified in the docstring.


[src/Compute.jl:TODO3]
    # TODO 3: Return one production row per input budget.
    # Validate the budgets and build the protein objective with your first function.
    # At each budget, construct and solve :equal and :optimized separately.
    # Pass the returned lower, upper, A, and b arrays to solve_lp, then check
    # its status before reading the protein flux. Convert each rate to μM/h.
    # Collect named tuples with fields budget, equal, and optimized in input order.


[src/Solver.jl:D01]
solve_lp(S, lower, upper, c; A=zeros(0, length(c)), b=Float64[]) -> NamedTuple

Maximize `dot(c, v)` subject to `S*v == 0`, the flux bounds, and `A*v <= b`.

# Arguments
- `S`: an m-by-n stoichiometric matrix, with species in rows and reactions in columns.
- `lower`, `upper`: length-n vectors defining `lower[j] <= v[j] <= upper[j]`.
- `c`: a length-n objective-coefficient vector in the same reaction order.
- `A`, `b`: a k-by-n matrix and length-k vector specifying extra inequalities.
  The defaults add no inequalities. All matrix and vector entries must be finite.

# Returns
A named tuple with fields `status` (solver termination status), `optimal`
(whether the status is `OPTIMAL`), `flux` (reaction-flux vector), and `objective`
(the scalar value of `dot(c, flux)`). If the status is not `OPTIMAL`, both `flux`
and `objective` are `nothing`. Leave all inputs unchanged.

# Errors
Throw `DimensionMismatch` if dimensions disagree, or `ArgumentError` for
nonfinite entries or a lower bound greater than its upper bound.

[src/Solver.jl:D02]
check_solution(model, result, constraints; atol=1e-7) -> NamedTuple

Check the returned fluxes against balances, bounds, and allocation inequalities.

# Arguments
- `model`: the model whose `S` defines the species balances.
- `result`: an optimal result from `solve_lp`, with a reaction-ordered `flux` vector.
- `constraints`: the named tuple `lower`, `upper`, `A`, `b` used in that solve.
- `atol::Real`: finite, nonnegative absolute tolerance for each balance, bound,
  and inequality check (mM/h for the PS3 constraints). Dimensions must agree.

# Returns
A named tuple containing `balance_ok`, `bounds_ok`, `allocation_ok`, their
conjunction `valid`, and `maximum_residual` (the largest absolute entry of
`model.S * result.flux`). These checks establish feasibility; optimality is
reported by the solver. Leave all inputs unchanged.

# Errors
Throw `ArgumentError` for a nonoptimal result or an invalid tolerance.

[src/Model.jl:D01]
reaction_index(model, id::AbstractString) -> Int

Find the reaction column identified by `id`.

# Arguments
- `model`: a named tuple whose `reactions` vector lists identifiers in matrix
  column order.
- `id::AbstractString`: the reaction identifier to locate.

# Returns and errors
Return its one-based column index as an `Int`. Throw `ArgumentError` if the
identifier is absent. Leave the model unchanged.

[src/Model.jl:D02]
load_model(; glucose_limit=1.0) -> NamedTuple

Load the supplied deGFP model with the PS3 operating conditions.

# Arguments
- `glucose_limit::Real`: a finite, nonnegative glucose supply limit in mM/h;
  the default is 1.0. Read model files from the data directory beside src/.

# Returns
A named tuple containing newly loaded model arrays and fixed expression rates:
- `S`: stoichiometric matrix, with species in rows and reactions in columns.
- `reactions`, `metabolites`: identifiers in column and row order.
- `lower`, `upper`: reaction-flux bounds in mM/h.
- `amino_acids`, `amino_acid_uptake`: names and matching supply-column indices.
- `protein`, `translation`: column indices for protein output and translation
  initiation, respectively.
- `transcription_rate`, `translation_limit`: fixed expression rate and capacity
  in mM/h. Separate supply and removal reactions have nonnegative fluxes.

# Errors
Throw `ArgumentError` for an invalid glucose limit and `ErrorException` if the
matrix dimensions disagree with the identifiers. File-reading and parsing errors
propagate to the caller. The data files are read without being changed.

[src/Reporting.jl:D01]
Supplied plotting and table output for the two PS3 allocation strategies.

[src/Reporting.jl:D02]
write_outputs(assignment, directory) -> Vector{NamedTuple}

Generate the production curves and the amino-acid uptake table.

# Arguments
- `assignment`: the loaded `CellFreeProduction` module, including the student's
  completed objective, allocation, and production-curve functions.
- `directory`: output-directory path; create it if needed.

# Returns and files
Return the production-curve rows for budgets 0:0.1:4 mM/h. Write or overwrite
`production-versus-budget.csv`, `.png`, and `.svg`, plus
`amino-acid-allocation.csv` for U=1 mM/h. Protein rates are in μM/h; supply
rates are in mM/h. Print the U=1 production rates, total uptake, and balance
residuals. Reset the plot theme to the default light appearance.

# Errors
Throw `ErrorException` for a curve with the wrong row count or budget order,
nonfinite production rates, or a nonoptimal or infeasible U=1 allocation.
Errors from student functions, file writing, or plotting propagate to the caller.

[Include.jl:D01]
CellFreeProduction

Load the supplied model and solver utilities and the three student functions.
Include this file once per Julia session to make the exported names available.
Restart Julia after editing source files, or run the scripts in a new process.