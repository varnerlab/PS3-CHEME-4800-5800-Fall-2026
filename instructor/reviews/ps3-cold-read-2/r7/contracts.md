ALREADY READ: the complete README.md, then src/Compute.jl up to the audited part

[README.md:B001 heading]
# Problem Set 3 (PS3): Feeding a cell-free protein factory

[README.md:B002 paragraph]
A cell-free protein synthesis system uses cellular machinery in a reaction
mixture to make a protein. That machinery needs amino acids, nucleotides, and
energy. Metabolic reactions in the mixture can make amino acids from other
compounds, and amino acids can also be supplied directly.

[README.md:B003 paragraph]
You will use linear programming to find the maximum production rate of deGFP, a
green fluorescent protein, when the 20 amino-acid supply rates share one budget,
a limit in mM/h on their sum. With **equal allocation**, each amino acid may use
at most 1/20 of the budget. With **optimized allocation**, the solver divides the
budget. A supplied script plots the maximum rate against the budget for both
strategies.

[README.md:B004 paragraph]
The same maximum can often be reached with different supply rates. Flux
variability analysis finds the smallest and largest value of each supply rate
among all solutions that reach it.

[README.md:B005 paragraph]
The model is adapted from **Vilkhovoy et al. (2018), “Sequence Specific
Modeling of E. coli Cell-Free Protein Synthesis,” ACS Synthetic Biology, 7(8),
1844–1857** ([paper](https://pubmed.ncbi.nlm.nih.gov/29944340/),
[DOI: 10.1021/acssynbio.7b00465](https://doi.org/10.1021/acssynbio.7b00465)).
The publication's code and data are in the Varner laboratory's
[sequence-specific cell-free protein synthesis repository](https://github.com/varnerlab/Sequence-Specific-FBA-CFPS-Publication-Code);
you do not need them. The [data notes](data/README.md) describe the source and
the PS3 operating conditions.

[README.md:B006 heading]
## Logistics

[README.md:B007 list]
- **Release:** Saturday, October 3, 2026. **Initial submission deadline:**
  **11:59 PM ET on Saturday, October 17, 2026**.
- **Submission:** Upload a ZIP of the assignment folder to Canvas, named
  `CHEME-4800-5800-PS3-<your netid>.zip` with your NetID in place of
  `<your netid>` (for example, `CHEME-4800-5800-PS3-abc123.zip`). The last
  section lists what it must contain.
- **Infinite-revision policy:** Submit something by the initial deadline; it does
  not need to run or pass any tests. You may then revise and resubmit until the
  end of the semester, and the highest score is retained. No initial submission
  means a score of `0` and no access to later revisions.
- **Reference solution:** Released after the initial deadline. You may consult
  it to understand mistakes and debug your work, but copying it results in a
  score of `0`.
- **Group and AI policy:** Submit independent work. You may discuss ideas with
  classmates but may not share code or solutions directly. You may use Julia
  documentation, AI tools, and internet resources. You are responsible for
  understanding and explaining your implementation.
- **Grading:** Maximum `4` points; see [RUBRIC.md](RUBRIC.md).

[README.md:B008 heading]
## Getting started

[README.md:B009 paragraph]
Use Julia `1.12.7`; the folder contains all the data and source files you need.
Open it in VS Code and run these commands from the folder containing this README
and [Project.toml](Project.toml):

[README.md:B010 code]
```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B011 paragraph]
The first command installs the required packages and can take several minutes.
The second runs the tests, which fail until you implement the four functions in
[src/Compute.jl](src/Compute.jl) described in Parts 1–3.

[README.md:B012 paragraph]
Implement those functions and answer the five questions in
[responses.md](responses.md). You may add helper files in `src/` and load them
from [Include.jl](Include.jl). Keep all other supplied files unchanged.

[README.md:B013 heading]
## The model and its operating conditions

[README.md:B014 paragraph]
The model has 146 species and 265 reactions covering metabolism, transcription,
translation, amino-acid supply, and protein output. The supplied
[load_model function](src/Model.jl), called as `model = load_model()`, returns a
named tuple with these fields:

[README.md:B015 table]
| Field | Meaning |
|---|---|
| `S` | Stoichiometric matrix; rows are species and columns are reactions. |
| `reactions`, `metabolites` | Identifiers in column and row order. |
| `lower`, `upper` | Lower and upper flux bound of each reaction, in mM/h. |
| `amino_acids` | Names of the 20 supplied amino acids. |
| `amino_acid_uptake` | Column indices of their supply reactions. |
| `protein` | Column index of the deGFP output reaction. |
| `translation` | Column index of translation initiation; each initiation makes one deGFP. |
| `translation_limit` | Translation capacity in mM/h, already set as `upper[translation]`. |

[README.md:B016 paragraph]
`load_model` already applies the operating conditions: glucose supply of at most
**1 mM/h**, oxygen supply of at most **100 mM/h**, and amino-acid supply between
0 and 30 mM/h. Amino-acid synthesis is enabled; reactions that degrade amino
acids or export them from the mixture are blocked. Transcription and translation
bounds come from fixed gene-expression parameters. Both strategies start from
these conditions; Part 1 describes how each one limits amino-acid supply.

[README.md:B017 paragraph]
Supply reactions are written as `[] --> species` and have **nonnegative**
fluxes, so each adds to the total that the budget limits. The deGFP output
reaction removes the protein that translation makes; its flux is the production
rate. Every species, including deGFP, is at steady state.

[README.md:B018 paragraph]
Model fluxes are in **mM/h**, millimoles per liter per hour. Protein production
is reported in **µM/h** in function results, plots, and answers: multiply the
protein-output flux by 1000. The supply budget is a rate limit, not an initial
concentration or a monetary cost. This static model predicts feasible production
rates under fixed conditions.

[README.md:B019 heading]
## Part 1: Construct the objective and supply constraints

[README.md:B020 paragraph]
Let $\mathbf{v}$ contain the reaction fluxes and let $p$ be the protein-output
column, `model.protein`. The production problem is:

[README.md:B021 display]
$$
\begin{aligned}
\underset{\mathbf{v}}{\operatorname{maximize}}\quad & v_p\\
\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\\
& \boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\\
& \mathbf{A}\mathbf{v}\leq\mathbf{b}.
\end{aligned}
$$

[README.md:B022 paragraph]
The bounds $\boldsymbol{\ell}$ and $\mathbf{u}$ start from the model's `lower` and
`upper`, and $\mathbf{A}$ and $\mathbf{b}$ hold any extra supply inequality. Let
$\mathcal{I}$ be the set of amino-acid supply columns, `model.amino_acid_uptake`,
let $K=|\mathcal{I}|$ be their number, and let the scalar $U\geq0$ be the total
supply budget in mM/h, the `budget` argument of your functions.

[README.md:B023 paragraph]
**Equal allocation.** Each supply flux's upper bound becomes the smaller of its
original bound $u_j$ and an equal share of the budget:

[README.md:B024 display]
$$
u_j^{\mathrm{equal}}=\min(u_j,U/K),\qquad j\in\mathcal{I}.
$$

[README.md:B025 paragraph]
Supply fluxes can stay below their limits, and an unused share is not available
to other amino acids. These limits already enforce the budget, so equal
allocation adds no inequality: `A = zeros(0, n)`, where `n` is the number of
reactions, and `b = Float64[]`.

[README.md:B026 paragraph]
**Optimized allocation.** Keep the original bounds and add one shared supply
constraint:

[README.md:B027 display]
$$
\sum_{j\in\mathcal{I}}v_j\leq U.
$$

[README.md:B028 paragraph]
The solver may divide the supply unequally. `A` is a 1×n matrix with ones at the
amino-acid supply columns and zeros elsewhere, and `b = [U]`. Under both
strategies, amino acids made by metabolism do not count as external supply.

[README.md:B029 paragraph]
Implement two functions in [src/Compute.jl](src/Compute.jl):

[README.md:B030 list]
1. **protein_objective(model)** returns a `Vector{Float64}` with one entry per
   reaction: one at the protein-output column and zero elsewhere.
2. **allocation_constraints(model, budget; strategy=:equal)** returns a named tuple
   with fields `lower`, `upper`, `A`, and `b` for `:equal` or `:optimized`.
   Return copies of the bound arrays so that `model.lower` and `model.upper` stay
   unchanged. Reject a negative or nonfinite budget or an unknown strategy with
   `ArgumentError`.

[README.md:B031 paragraph]
Use the model's column indices: the test scripts also build small networks with
the same fields but different numbers and orders of reactions and supplies. All
test networks have nonnegative supply fluxes and at least one supply column. Only
the supply upper bounds change, and only under `:equal`.

[README.md:B032 paragraph]
Run the Part 1 checks with:

[README.md:B033 code]
```bash
julia --startup-file=no --project=. testme_part_1.jl
```

[README.md:B034 paragraph]
For a hand check, take two amino acids with supply fluxes $v_1$ and $v_2$. One
unit of protein output $v_p$ consumes two units of the first and one of the
second, so $v_1=2v_p$ and $v_2=v_p$. With a total supply budget of 6 and all other
bounds loose, equal allocation allows at most **1.5** units of protein output and
optimized allocation allows **2**. Make sure you can explain these numbers.

[README.md:B035 heading]
## Part 2: Compute and interpret the production curves

[README.md:B036 paragraph]
Implement **production_curve(model, budgets)** in
[src/Compute.jl](src/Compute.jl). For each budget, build both strategies'
constraints with your Part 1 functions and solve each with the supplied
[solve_lp function](src/Solver.jl), which maximizes $\mathbf{c}^{\top}\mathbf{v}$
subject to the balances, bounds, and $\mathbf{A}\mathbf{v}\leq\mathbf{b}$. A single
solve has this form:

[README.md:B037 code]
```julia
limits = allocation_constraints(model, budget; strategy=:equal);
c = protein_objective(model);
result = solve_lp(model.S, limits.lower, limits.upper, c;
    A=limits.A, b=limits.b);
```

[README.md:B038 paragraph]
The result contains `status`, `optimal`, `flux`, and `objective`. When a solve is
not optimal, `flux` and `objective` are `nothing`, so check `result.optimal` first
and throw an `ErrorException` naming the budget, strategy, and status.

[README.md:B039 paragraph]
Return one named tuple per input budget, in input order and including repeats,
with fields `budget` (mM/h) and `equal` and `optimized` (maximum production rates
in µM/h). An empty input returns an empty vector; a negative or nonfinite budget
raises `ArgumentError`.

[README.md:B040 paragraph]
After completing the function, run:

[README.md:B041 code]
```bash
julia --startup-file=no --project=. testme_part_2.jl
julia --startup-file=no --project=. runproduction.jl
```

[README.md:B042 paragraph]
The supplied driver, `runproduction.jl`, evaluates budgets from **0 to 4 mM/h in
steps of 0.1** and writes these files:

[README.md:B043 list]
- `outputs/production-versus-budget.png` and `.svg`: both production curves and
  the translation capacity, the highest rate translation allows. The terminal
  prints this capacity in µM/h.
- `outputs/production-versus-budget.csv`: the plotted values.
- `outputs/amino-acid-allocation.csv`: the supply fluxes of one optimal solution
  per strategy at a budget of **1 mM/h**. The driver re-solves these with your
  Part 1 functions, checks the balances, bounds, and budget, and prints each
  production rate and total supply.

[README.md:B044 paragraph]
Until `supply_ranges` (Part 3) is complete, the driver skips the Part 3 files and
says so. Use the plot and CSV tables to answer questions 1 and 2 in
[responses.md](responses.md).

[README.md:B045 paragraph]
The optimized column of your `amino-acid-allocation.csv` may differ from a
classmate's: optimal solutions are often not unique, and the tests accept any
that satisfies the constraints and gives the optimal production rate. Part 3
measures how much they can differ.

[README.md:B046 heading]
## Part 3: Which supply fluxes does the optimum determine?

[README.md:B047 paragraph]
The optimal production rate of a linear program is unique, but the flux vector
that achieves it often is not. **Flux variability analysis (FVA)** finds the
smallest and largest value of each flux among all flux vectors that reach the
optimum or, more generally, at least a fraction $\gamma$ of it. FVA, including
this fractional form, was introduced by **Mahadevan and Schilling (2003),
“The effects of alternate optimal solutions in constraint-based genome-scale
metabolic models,” Metabolic Engineering, 5(4), 264–276**
([paper](https://pubmed.ncbi.nlm.nih.gov/14642354/),
[DOI: 10.1016/j.ymben.2003.09.002](https://doi.org/10.1016/j.ymben.2003.09.002)).

[README.md:B048 paragraph]
In PS3, $\gamma$ is the `fraction` argument of `supply_ranges`, with
$0\leq\gamma\leq1$. The default, $\gamma=1$, keeps only optimal flux vectors;
$\gamma=0.99$ shows how far the supplies can shift if production may fall 1%.

[README.md:B049 paragraph]
PS3 applies FVA only to the amino-acid supply fluxes. Many internal reactions
have separate forward and reverse columns; running both at the same rate changes
no balance, so their ranges mostly reflect their bounds. Supply reactions have no
reverse partner.

[README.md:B050 paragraph]
Let $\mathbf{c}$ be the protein objective and $z^\star$ the maximum
protein-output flux, in mM/h rather than µM/h, for one strategy at budget $U$.
For each supply column $j\in\mathcal{I}$, the range endpoints are:

[README.md:B051 display]
$$
\begin{aligned}
v_j^{\min}=\min_{\mathbf{v}}\ v_j,\qquad v_j^{\max}=\max_{\mathbf{v}}\ v_j
\quad\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\quad
\boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\quad \mathbf{A}\mathbf{v}\leq\mathbf{b},\\
& \mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star.
\end{aligned}
$$

[README.md:B052 paragraph]
The bounds, $\mathbf{A}$, and $\mathbf{b}$ are those that `allocation_constraints`
returns for the same strategy and budget.

[README.md:B053 paragraph]
A supply flux is **fixed** when $v_j^{\max}-v_j^{\min}\leq10^{-5}$ mM/h; at
$\gamma=1$, the optimum then determines its value. The driver and questions 3–5
use this rule; `supply_ranges` returns only the endpoints.

[README.md:B054 paragraph]
Implement **supply_ranges(model, budget; strategy=:equal, fraction=1.0)** in
[src/Compute.jl](src/Compute.jl), where `fraction` is $\gamma$. Return one named
tuple per supply column, in the order of `model.amino_acid_uptake`, with fields:

[README.md:B055 list]
- `column`: the column index in `S`;
- `minimum` and `maximum`: $v_j^{\min}$ and $v_j^{\max}$, in mM/h.

[README.md:B056 paragraph]
Because `solve_lp` only maximizes and accepts only
$\mathbf{A}\mathbf{v}\leq\mathbf{b}$, compute the ranges in three steps:

[README.md:B057 list]
1. Build the strategy's limits and $\mathbf{c}$ with your Part 1 functions and
   solve once. Take $z^\star$ as this solve's `objective`, in mM/h.
2. Write $\mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star$ as
   $-\mathbf{c}^{\top}\mathbf{v}\leq-\gamma z^\star$ and append it as one row:
   `vcat(A, -c')` and `vcat(b, -fraction * z)`, where `z` holds $z^\star$. Use
   `vcat`, not `push!`, which cannot add a row to `A` and fails if an integer
   budget made your `b` a `Vector{Int}`.
3. For each supply column $j$, solve twice with the step 1 bounds and the new `A`
   and `b`. For $v_j^{\max}$, maximize an objective that is one at column $j$ and
   zero elsewhere. For $v_j^{\min}$, use minus one at column $j$ and negate that
   solve's `objective`.

[README.md:B058 paragraph]
Use $\gamma z^\star$ exactly. Some FVA codes lower it by a small tolerance to
guard against round-off, but near the optimum a tiny loss of production can
allow a large shift in supply: lowering $\gamma z^\star$ by only $10^{-9}$ mM/h
can turn a fixed supply flux into one whose range is wider than $10^{-5}$ mM/h.
With the supplied solver, every step 3 solve is feasible at exactly
$\gamma z^\star$.

[README.md:B059 paragraph]
Throw `ArgumentError` for a negative or nonfinite budget, an unknown strategy, or
a `fraction` that is nonfinite or outside $[0,1]$. If any solve is not optimal,
throw an `ErrorException` naming the budget, strategy, and status.

[README.md:B060 paragraph]
For a hand check, return to the two amino acids, but now one unit of $v_p$
consumes one unit of each, and a conversion flux $v_c\geq0$ turns the second into
the first, so $v_1+v_c=v_p$ and $v_2=v_c+v_p$. With a budget of 6 and all other
bounds loose, both strategies reach $v_p=3$. Equal allocation fixes both supplies
at **3**. Optimized allocation fixes only their sum at 6, so $v_1$ ranges over
**[0, 3]** and $v_2$ over **[3, 6]**. The optimum sets the total supply but not
its split.

[README.md:B061 paragraph]
Run the Part 3 checks and then the driver again:

[README.md:B062 code]
```bash
julia --startup-file=no --project=. testme_part_3.jl
julia --startup-file=no --project=. runproduction.jl
```

[README.md:B063 paragraph]
The driver adds these files to `outputs/`:

[README.md:B064 list]
- `supply-ranges.png`, `.svg`, and `.csv`: supply ranges at **1 mM/h** for equal
  allocation, optimized allocation, and optimized allocation with $\gamma=0.99$.
  The table's columns are `amino_acid`, `equal_min_mM_per_h`,
  `equal_max_mM_per_h`, `optimized_min_mM_per_h`, `optimized_max_mM_per_h`,
  `optimized_99pct_min_mM_per_h`, and `optimized_99pct_max_mM_per_h`. The
  terminal reports how many supply fluxes are fixed in each case.
- `range-width-versus-budget.png`, `.svg`, and `.csv`: one row per budget from 0
  to 4 mM/h, both strategies at $\gamma=1$. `equal_total_width_mM_per_h` and
  `optimized_total_width_mM_per_h` sum the 20 range widths
  $v_j^{\max}-v_j^{\min}$; `equal_fixed_supplies` and `optimized_fixed_supplies`
  count the fixed supply fluxes.

[README.md:B065 paragraph]
Use these files to answer questions 3–5 in [responses.md](responses.md). Unlike
optimal flux vectors, range endpoints are unique, so the tests compare them with
reference values.

[README.md:B066 heading]
## Checking and submitting your work

[README.md:B067 paragraph]
Run the complete check after finishing the code and responses:

[README.md:B068 code]
```bash
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B069 paragraph]
The checker runs **24 tests per part, 72 total**, records the results and a
SHA-256 fingerprint of each source file in `MANIFEST.txt`, and warns about
unanswered questions in `responses.md`. If all three test suites pass, it also
regenerates `outputs/` from your code; otherwise, run the driver after your final
edits and then run the checker last. For grading, the instructor reruns the
checker on your code and, if all tests pass, reviews the implementation and
written answers; see [RUBRIC.md](RUBRIC.md).

[README.md:B070 paragraph]
Zip the entire assignment folder, including your source, data, responses,
`MANIFEST.txt`, and `outputs/`. Name the archive
`CHEME-4800-5800-PS3-<your netid>.zip` and upload it to Canvas yourself; the
checker runs only locally.

[README.md:B071 paragraph]
If some tests still fail at the deadline, submit your current work to receive
partial credit and participate in the infinite-revision policy.

[src/Compute.jl:H01 file comment]
# Complete these four functions. Keep their arguments and return fields.
# Replace each placeholder error(...) with your implementation and return value.
# Use the supplied solve_lp(...) function; keep the model and its arrays unchanged.

[src/Compute.jl:D01 docstring]
"""
    protein_objective(model) -> Vector{Float64}

Build the objective for maximizing protein output in mM/h.

# Arguments
- `model`: a named tuple with a stoichiometric matrix `S` and the integer
  column index `protein` of the protein-output reaction.

# Returns
- A new `Vector{Float64}` with one entry per column of `model.S`. The entry at
  `model.protein` is one; every other entry is zero. Leave the model unchanged.
"""
function protein_objective(model)

[src/Compute.jl:TODO1 TODO comment]
    # TODO 1: Build and return the objective vector.
    # Use the number of columns in model.S for its length and model.protein
    # for the nonzero entry; the tests also use smaller, reordered networks.
    error("Complete protein_objective in src/Compute.jl");

[src/Compute.jl:D02 docstring]
"""
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
"""
function allocation_constraints(model, budget::Real; strategy::Symbol=:equal)

[src/Compute.jl:TODO2 TODO comment]
    # TODO 2: Build and return the named tuple (lower, upper, A, b).
    # Check the budget and strategy, then copy the model's bound vectors.
    # For :equal, change only the supply upper bounds using the actual number
    # of supply columns. For :optimized, build the single shared-budget row.
    # Return arrays with the dimensions specified in the docstring.
    error("Complete allocation_constraints in src/Compute.jl");

[src/Compute.jl:D03 docstring]
"""
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
"""
function production_curve(model, budgets::AbstractVector{<:Real})

[src/Compute.jl:TODO3 TODO comment]
    # TODO 3: Return one production row per input budget.
    # Validate the budgets and build the protein objective with your first function.
    # At each budget, construct and solve :equal and :optimized separately.
    # Pass the returned lower, upper, A, and b arrays to solve_lp, then check
    # its status before reading the protein flux. Convert each rate to μM/h.
    # Collect named tuples with fields budget, equal, and optimized in input order.
    error("Complete production_curve in src/Compute.jl");

AUDIT THIS PART: src/Compute.jl:D04 to src/Compute.jl:TODO4

[src/Compute.jl:D04 docstring]
"""
    supply_ranges(model, budget::Real; strategy::Symbol=:equal, fraction::Real=1.0) -> Vector{NamedTuple}

Compute the range of every amino-acid supply flux for one budget and strategy.

# Arguments
- `model`: a model accepted by `protein_objective` and `allocation_constraints`.
- `budget::Real`: a finite, nonnegative total amino-acid supply limit in mM/h.
- `strategy::Symbol`: `:equal` (default) or `:optimized`.
- `fraction::Real`: the required fraction γ of the maximum protein output, with
  `0 <= fraction <= 1`. The default 1.0 keeps only optimal flux vectors.

# Returns
One named tuple per supply column, in the order of `model.amino_acid_uptake`,
with fields `column` (the column index in `model.S`), `minimum`, and `maximum`.
The endpoints are the smallest and largest values of that supply flux, in mM/h,
over all flux vectors that satisfy the strategy's constraints and reach at least
`fraction` times the maximum protein output. Leave the model unchanged.

# Method
1. Build the limits with `allocation_constraints` and `c = protein_objective(model)`,
   then maximize protein output with `solve_lp`; `z` is its `objective`, in mM/h.
2. Build `vcat(A, -c')` and `vcat(b, -fraction*z)` from the limits, which require
   protein output of at least `fraction*z`. Use `vcat`, not `push!`, which cannot
   add a row to `A` and fails if an integer budget made `b` a `Vector{Int}`. Do
   not relax the requirement by a tolerance.
3. For each supply column `j`, solve twice with the strategy's bounds and the new
   arrays. For the maximum, maximize an objective that is one at `j` and zero
   elsewhere; for the minimum, use minus one at `j` and negate the `objective`.

# Errors
Throw `ArgumentError` for an invalid budget or strategy, or for a `fraction` that
is nonfinite or outside [0, 1]. If any solve is not optimal, throw `ErrorException`
with the budget, strategy, and solver status in its message.
"""
function supply_ranges(model, budget::Real; strategy::Symbol=:equal, fraction::Real=1.0)

[src/Compute.jl:TODO4 TODO comment]
    # TODO 4: Return one (column, minimum, maximum) named tuple per supply column.
    # Check the fraction, build the strategy's limits, and solve once for the
    # maximum protein output z in mM/h. With vcat (not push!), append one row
    # requiring protein output of at least fraction*z. For each supply column j,
    # maximize v[j] for the maximum; maximize -v[j] and negate result.objective
    # for the minimum. Check result.optimal after every solve, before reading it.
    # Collect named tuples (column=j, minimum=lo, maximum=hi) in mM/h.
    error("Complete supply_ranges in src/Compute.jl");
