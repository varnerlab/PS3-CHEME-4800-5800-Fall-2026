ALREADY READ: the complete README.md, then src/Compute.jl up to the audited part

[README.md:B001 heading]
# Problem Set 3 (PS3): Feeding a cell-free protein factory

[README.md:B002 paragraph]
A cell-free protein synthesis system uses cellular machinery in a reaction
mixture to make a protein. That machinery needs amino acids, nucleotides, and
energy. The metabolic reactions in the mixture can make amino acids from other
available compounds, and amino acids can also be supplied directly.

[README.md:B003 paragraph]
You will use linear programming to find the maximum production rate of deGFP, a
green fluorescent protein, when the amino-acid supply rates share one total
budget: a limit, in mM/h, on the sum of the 20 supply rates. Two strategies divide
that budget. With **equal allocation**, each amino acid may use at most an equal
share. With **optimized allocation**, the solver decides how to divide the budget
among the amino acids. A supplied driver plots the maximum production rate against
the budget for both strategies. Different sets of supply rates can reach the same
maximum, so you will then use flux variability analysis to find, for each strategy
and budget, the smallest and largest value each supply rate can take at that
maximum.

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
- **Submission:** Zip the assignment folder, name the archive
  `CHEME-4800-5800-PS3-<your netid>.zip` with `<your netid>`, brackets included,
  replaced by your NetID, and upload it to Canvas. The last section of this README
  lists what the archive must contain.
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
The second runs the assignment tests. The four functions in
[src/Compute.jl](src/Compute.jl) initially raise errors because you must implement
them. Tests will fail until that work is complete.

[README.md:B011 paragraph]
Implement the four functions in `src/Compute.jl` and answer the five questions in
[responses.md](responses.md). You may add helper files in `src/`; if you do, load
them from [Include.jl](Include.jl). Keep all other supplied files unchanged.

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
| `translation` | Column index of translation initiation, whose capacity caps protein production. |
| `translation_limit` | That capacity, in mM/h; `upper[translation]` already equals it. |

[README.md:B015 paragraph]
The `load_model` function already sets these operating conditions in `lower` and
`upper`. The glucose supply limit is **1 mM/h** and the oxygen supply limit is
**100 mM/h**. Each amino-acid supply flux starts with bounds of 0 and 30 mM/h.
Amino-acid synthesis is enabled, while the reactions that degrade amino acids or
export them from the mixture are blocked. The transcription and translation bounds
are computed from fixed gene-expression parameters in the supplied data. Use these
conditions for both strategies. Your functions change only the amino-acid supply
limits, through their upper bounds or one added inequality, as Part 1 describes.

[README.md:B016 paragraph]
Supply reactions are written as `[] --> species` and have **nonnegative**
fluxes, so each amino-acid supply flux adds a nonnegative amount to the total that
the budget limits. The deGFP output reaction removes the protein that translation
makes, so its flux is the production rate. Every species, including deGFP,
satisfies a steady-state balance.

[README.md:B017 paragraph]
All model fluxes are in **mM/h**, millimoles per liter per hour. Report protein
production in **µM/h** by multiplying its model flux by 1000. The supply budget
is a rate limit; it does not specify an initial concentration or a monetary cost.
This static model predicts feasible production rates under fixed conditions.

[README.md:B018 heading]
## Part 1: Construct the objective and supply constraints

[README.md:B019 paragraph]
Let the vector $\mathbf{v}$ contain the reaction fluxes, and let $p$ be the
protein-output column, `model.protein`. The production problem is given by:

[README.md:B020 display]
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
inequality. Define $\mathcal{I}$ as the set of amino-acid supply columns,
`model.amino_acid_uptake`, and let $K=|\mathcal{I}|$ be their number. Let
$U\geq0$ be the total supply budget in mM/h, the `budget` argument of your functions.

[README.md:B022 paragraph]
**Equal allocation.** Each amino acid may use at most one equal share of the
budget. The upper bound on its supply flux becomes the smaller of the model's
original bound $u_j$ and that share:

[README.md:B023 display]
$$
u_j^{\mathrm{equal}}=\min(u_j,U/K),\qquad j\in\mathcal{I}.
$$

[README.md:B024 paragraph]
The supply fluxes can be smaller than their limits. Any unused share remains
unavailable to the other amino acids. These individual limits already enforce
the total budget, so equal allocation adds no inequality: `A = zeros(0, n)`,
where `n` is the number of reactions, and `b = Float64[]`.

[README.md:B025 paragraph]
**Optimized allocation.** Retain the original individual bounds and add one
shared supply constraint, given by:

[README.md:B026 display]
$$
\sum_{j\in\mathcal{I}}v_j\leq U.
$$

[README.md:B027 paragraph]
The solver may distribute the available supply unequally. `A` is a single row,
a 1×n matrix, with a coefficient of one at each amino-acid supply column and zero
elsewhere, and `b = [U]`. Metabolic synthesis of amino acids remains
available under both strategies and does not count as external supply.

[README.md:B028 paragraph]
Implement two functions in [src/Compute.jl](src/Compute.jl):

[README.md:B029 list]
1. **protein_objective(model)** returns a `Vector{Float64}` with one coefficient
   per reaction. Its protein-output coefficient is one; the others are zero.
2. **allocation_constraints(model, budget; strategy=:equal)** returns a named tuple
   with fields `lower`, `upper`, `A`, and `b`. Support `:equal` and `:optimized`.
   Return copies of the bound arrays, so that building one set of constraints
   never changes the model used for the next.
   Reject a negative or nonfinite budget or an unknown strategy with `ArgumentError`.

[README.md:B030 paragraph]
Use the model's supplied column indices. Your functions must also work on the
small test networks that the test scripts build. These have the same fields as the
model from `load_model` but different numbers and orders of reactions and
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
**2**. Explain these numbers to yourself before moving to the full model in Part 2.

[README.md:B034 heading]
## Part 2: Compute and interpret the production curves

[README.md:B035 paragraph]
Implement **production_curve(model, budgets)** in
[src/Compute.jl](src/Compute.jl). For each budget, construct both sets of
constraints and solve both LPs. Use your Part 1 functions and the supplied
[solve_lp function](src/Solver.jl), which maximizes $\mathbf{c}^{\top}\mathbf{v}$
subject to the balances, bounds, and $\mathbf{A}\mathbf{v}\leq\mathbf{b}$. A single
solve has this form:

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
Return a vector containing one named tuple `(; budget, equal, optimized)` per
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
The supplied driver, `runproduction.jl`, evaluates budgets from **0 to 4 mM/h in
steps of 0.1** and writes these files:

[README.md:B042 list]
- `outputs/production-versus-budget.png` and `.svg`: the two production curves
  and the translation capacity as a protein production rate in µM/h; each
  translation initiation makes one deGFP. The terminal also prints this capacity.
- `outputs/production-versus-budget.csv`: the numerical values used in the plot.
- `outputs/amino-acid-allocation.csv`: actual supply fluxes for both strategies
  at a total budget of **1 mM/h**. The driver first checks that each solution
  satisfies the balances, bounds, and budget, then prints its production rate
  and total supply.

[README.md:B043 paragraph]
Until you complete `supply_ranges`, the fourth function, described in Part 3, the
driver skips the Part 3 files and says so. Open the plot and use it with the CSV
tables to answer questions 1 and 2 in [responses.md](responses.md).

[README.md:B044 paragraph]
Your allocation table may differ from a classmate's. The optimized allocation can
have multiple equally productive flux distributions, and the tests accept any that
satisfy the constraints and give the optimal production rate. Part 3 measures how
much those distributions differ.

[README.md:B045 heading]
## Part 3: Which supply fluxes does the optimum determine?

[README.md:B046 paragraph]
The optimal production rate of a linear program is unique, but the flux vector
that achieves it often is not. L6a introduced **flux variability analysis (FVA)**,
which finds the smallest and largest value of each flux among all flux vectors
that reach the optimum. FVA can also include flux vectors whose protein output is
at least a fraction $\gamma$ of the optimum. In PS3, $\gamma$ is an argument of
the function you write, with $0\leq\gamma\leq1$; $\gamma=1$ keeps only optimal
flux vectors. Most of PS3 uses $\gamma=1$, and one driver case uses
$\gamma=0.99$ to show what a 1% loss of production allows. FVA, including this
fractional form, was introduced by **Mahadevan and Schilling (2003),
“The effects of alternate optimal solutions in constraint-based genome-scale
metabolic models,” Metabolic Engineering, 5(4), 264–276**
([paper](https://pubmed.ncbi.nlm.nih.gov/14642354/),
[DOI: 10.1016/j.ymben.2003.09.002](https://doi.org/10.1016/j.ymben.2003.09.002)).

[README.md:B047 paragraph]
PS3 applies FVA only to the amino-acid supply fluxes, 20 in the full model. Many
internal reactions are written as separate forward and reverse columns. Running
both columns of such a pair at the same rate changes no balance, so each column
can rise to its upper bound whatever the optimum is, and its range says little.

[README.md:B048 paragraph]
Let $\mathbf{c}$ be the protein objective and $z^\star$ the maximum protein-output
flux, in mM/h, found by solving the production LP for one allocation strategy at
budget $U$. For each supply column $j\in\mathcal{I}$, the range endpoints are
given by:

[README.md:B049 display]
$$
\begin{aligned}
v_j^{\min}=\min_{\mathbf{v}}\ v_j,\qquad v_j^{\max}=\max_{\mathbf{v}}\ v_j
\quad\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\quad
\boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\quad \mathbf{A}\mathbf{v}\leq\mathbf{b},\\
& \mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star.
\end{aligned}
$$

[README.md:B050 paragraph]
The bounds, $\mathbf{A}$, and $\mathbf{b}$ are those that `allocation_constraints`
returns for the chosen strategy.

[README.md:B051 paragraph]
The driver and questions 3–5 call a supply flux **fixed** when
$v_j^{\max}-v_j^{\min}\leq10^{-5}$ mM/h; at $\gamma=1$, the optimum then
determines its value. This is a reporting rule only; `supply_ranges`, described
next, returns the range endpoints and no fixed flag.

[README.md:B052 paragraph]
Implement **supply_ranges(model, budget; strategy=:equal, fraction=1.0)** in
[src/Compute.jl](src/Compute.jl), where `fraction` is $\gamma$. It returns a
vector with one named tuple per supply column, in the order of
`model.amino_acid_uptake`. Each tuple has three fields:

[README.md:B053 list]
- `column`: the supply flux's column index in `S`;
- `minimum` and `maximum`: $v_j^{\min}$ and $v_j^{\max}$, in mM/h.

[README.md:B054 paragraph]
Compute the ranges in three steps:

[README.md:B055 list]
1. Build the strategy's limits and the protein objective $\mathbf{c}$ with your
   Part 1 functions, and solve once with `solve_lp`. Take $z^\star$ as this solve's
   `objective`, in mM/h, not the µM/h rate from `production_curve`.
2. Add the requirement $\mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star$ as one more
   row of $\mathbf{A}\mathbf{v}\leq\mathbf{b}$: build new arrays with `vcat`,
   appending $-\mathbf{c}^{\top}$ to `A` and $-\gamma z^\star$ to `b`. Avoid
   `push!`. With `:optimized` and an integer budget, your Part 1 `b` may be a
   `Vector{Int}`, which cannot hold $-\gamma z^\star$; `vcat` returns a new
   `Float64` vector instead.
3. For each supply column $j$, solve twice with the step 1 bounds and the new
   `A` and `b`. For $v_j^{\max}$, maximize with an objective that is one at column
   $j$ and zero elsewhere. Because `solve_lp` only maximizes, find $v_j^{\min}$ by
   maximizing with an objective that is minus one at column $j$ and negating that
   solve's `objective`.

[README.md:B056 paragraph]
Some FVA codes lower $\gamma z^\star$ by a small tolerance to avoid infeasible
solves. Do not: lowering it by only $10^{-9}$ mM/h can turn a fixed supply flux into
one whose range is wider than $10^{-5}$ mM/h. With the supplied solver, every solve
in step 3 is feasible at exactly $\gamma z^\star$.

[README.md:B057 paragraph]
Reject an invalid budget or strategy, as in Part 1, or a `fraction` that is
nonfinite or outside $[0,1]$, with `ArgumentError`. If the solve in step 1 or any
solve in step 3 is not optimal, throw an `ErrorException` naming the budget,
strategy, and status.

[README.md:B058 paragraph]
For a hand check, consider the two amino acids again, but now one unit of $v_p$
consumes one unit of each species, and a conversion flux $v_c\geq0$ turns the
second species into the first. The balances are $v_1+v_c=v_p$ and
$v_2=v_c+v_p$. With a budget of 6 and all other bounds loose, both strategies
reach $v_p=3$. Equal allocation fixes both supplies at **3**. Optimized
allocation fixes only their sum at 6, so $v_1$ ranges over **[0, 3]** and $v_2$
over **[3, 6]**. The optimum chooses the total supply but not how it is split.

[README.md:B059 paragraph]
Run the Part 3 checks and then the driver again:

[README.md:B060 code]
```bash
julia --startup-file=no --project=. testme_part_3.jl
julia --startup-file=no --project=. runproduction.jl
```

[README.md:B061 paragraph]
The driver adds these files to `outputs/`:

[README.md:B062 list]
- `supply-ranges.png`, `.svg`, and `.csv`: supply ranges at a budget of **1 mM/h**
  for equal allocation, optimized allocation, and optimized allocation with
  $\gamma=0.99$. The table has an `amino_acid` column and a minimum and maximum
  column for each case, such as `optimized_min_mM_per_h` and
  `optimized_99pct_max_mM_per_h`. The terminal reports how many supply fluxes
  are fixed in each case.
- `range-width-versus-budget.png`, `.svg`, and `.csv`: one row per budget from 0
  to 4 mM/h, for both strategies at $\gamma=1$. The columns
  `equal_total_width_mM_per_h` and `optimized_total_width_mM_per_h` hold the sum
  of the 20 range widths $v_j^{\max}-v_j^{\min}$. The columns
  `equal_fixed_supplies` and `optimized_fixed_supplies` count the fixed supply
  fluxes.

[README.md:B063 paragraph]
Use these files to answer questions 3–5 in [responses.md](responses.md).
Unlike the optimal flux vectors, the range endpoints are unique, so the tests
compare them with reference values.

[README.md:B064 heading]
## Checking and submitting your work

[README.md:B065 paragraph]
Run the complete check after finishing the code and responses:

[README.md:B066 code]
```bash
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B067 paragraph]
The checker runs **24 tests per part, 72 total**, records their results and source
file fingerprints in `MANIFEST.txt`, and warns about unanswered discussion prompts.
If all three test suites pass, it also regenerates the plots and data tables in
`outputs/` from your code; otherwise it leaves `outputs/` as your last
`runproduction.jl` run wrote it.
The instructor still reviews the implementation and written explanations.

[README.md:B068 paragraph]
Zip the entire assignment folder, including your source, data, responses,
`MANIFEST.txt`, and generated files in `outputs/`. Name the archive
`CHEME-4800-5800-PS3-<your netid>.zip` and upload it to Canvas. The checker runs
locally; you must upload the ZIP yourself.

[README.md:B069 paragraph]
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

Compute the flux variability of every amino-acid supply flux at one budget.

# Arguments
- `model`: a model accepted by `protein_objective` and `allocation_constraints`.
- `budget::Real`: a finite, nonnegative total amino-acid supply limit in mM/h.
- `strategy::Symbol`: `:equal` (default) or `:optimized`.
- `fraction::Real`: the required fraction γ of the maximum protein output, with
  `0 <= fraction <= 1`. The default 1.0 keeps only optimal flux distributions.

# Returns
A vector with one named tuple, with fields `column`, `minimum`, and `maximum`, per
supply column, in the order of `model.amino_acid_uptake`. `column` is the column
index in `model.S`.
`minimum` and `maximum` are the smallest and largest values of that supply flux,
in mM/h, over all flux vectors that satisfy the strategy's constraints and reach
at least `fraction` times the maximum protein output. Leave the model unchanged.

First maximize protein output under the strategy's constraints; `z` is that
solve's `objective`, in mM/h. Then build new arrays with `vcat`, appending the row
`-c'` to `A` and the value `-fraction*z` to `b`, where `c` is
`protein_objective(model)`; the new row requires protein output of at least
`fraction*z`. Do not relax this requirement by a tolerance. Avoid `push!`, because
the strategy's `b` can hold integers. With the extra row, solve twice for each
supply column `j`. For the maximum, maximize with an objective that is one at `j`
and zero elsewhere. For the minimum, maximize with an objective that is minus one
at `j` and zero elsewhere, then negate that solve's `objective` value.

# Errors
Throw `ArgumentError` for an invalid budget or strategy, or for a `fraction` that
is nonfinite or outside [0, 1]. If any solve is not optimal, throw `ErrorException`
with the budget, strategy, and solver status in its message.
"""
function supply_ranges(model, budget::Real; strategy::Symbol=:equal, fraction::Real=1.0)

[src/Compute.jl:TODO4 TODO comment]
    # TODO 4: Return one flux-variability row per supply column.
    # Check the fraction, build the strategy's limits, and solve once for the
    # maximum protein output z in mM/h. Use vcat to build new A and b with one
    # more row requiring protein output of at least fraction*z; push! can fail
    # because Part 1's b can hold integers. For each supply column j, maximize
    # v[j] for the maximum; maximize -v[j] and negate result.objective for the
    # minimum (solve_lp only maximizes). Check result.optimal after every solve.
    # Collect named tuples (column=j, minimum=lo, maximum=hi) in mM/h.
    error("Complete supply_ranges in src/Compute.jl");
