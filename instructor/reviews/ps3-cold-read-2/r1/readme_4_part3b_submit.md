ALREADY READ: README.md up to the audited part

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
You will plot protein production against the total supply budget. You will then
use flux variability analysis to find which supply fluxes the optimum determines
and which it leaves open.

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
The second runs the assignment tests. The four functions in
[src/Compute.jl](src/Compute.jl) initially raise errors because you must implement
them. Tests will fail until that work is complete.

[README.md:B011 paragraph]
Implement the functions; answer the prompts in [responses.md](responses.md).
Keep helper files in `src/`; edit [Include.jl](Include.jl) to load them.
Edit `src/Compute.jl`, `Include.jl`, and `responses.md` and add helpers.
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
inequality. Define $\mathcal{I}$ as the set of amino-acid supply columns, let
$K=|\mathcal{I}|$, and let $U\geq0$ be the total supply budget in mM/h.

[README.md:B022 paragraph]
**Equal allocation.** Each amino acid may use at most one equal share of the
budget. Its upper bound is given by:

[README.md:B023 display]
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

[README.md:B026 display]
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
The supplied driver evaluates budgets from **0 to 4 mM/h in steps of 0.1** and
writes these files. Until `supply_ranges` is complete, it skips the Part 3 files
and says so:

[README.md:B042 list]
- `outputs/production-versus-budget.png` and `.svg`: the two production curves
  and the fixed translation-capacity limit.
- `outputs/production-versus-budget.csv`: the numerical values used in the plot.
- `outputs/amino-acid-allocation.csv`: actual supply fluxes for both strategies
  at a total budget of **1 mM/h**, with independent feasibility checks reported
  in the terminal.

[README.md:B043 paragraph]
Open the plot and use it with the tables to answer questions 1 and 2 in
[responses.md](responses.md). The optimized allocation can have multiple equally
productive flux distributions. Tests check feasibility and production, so they
accept alternative optimal allocations. Part 3 measures how much those
distributions differ.

[README.md:B044 heading]
## Part 3: Which supply fluxes does the optimum determine?

[README.md:B045 paragraph]
The optimal production rate of a linear program is unique, but the flux vector
that achieves it often is not. L6a introduced **flux variability analysis (FVA)**,
which finds the smallest and largest value of each flux among all flux vectors
that reach the optimum. FVA was introduced by **Mahadevan and Schilling (2003),
“The effects of alternate optimal solutions in constraint-based genome-scale
metabolic models,” Metabolic Engineering, 5(4), 264–276**
([paper](https://pubmed.ncbi.nlm.nih.gov/14642354/),
[DOI: 10.1016/j.ymben.2003.09.002](https://doi.org/10.1016/j.ymben.2003.09.002)).
Their method also covers suboptimal objective values, which PS3 represents with
a required fraction $\gamma$ of the optimum.

[README.md:B046 paragraph]
Let $\mathbf{c}$ be the protein objective and $z^\star$ the maximum production
rate under one allocation strategy. For each supply column $j\in\mathcal{I}$, the
range endpoints are given by:

[README.md:B047 display]
$$
\begin{aligned}
v_j^{\min}=\min_{\mathbf{v}}\ v_j,\qquad v_j^{\max}=\max_{\mathbf{v}}\ v_j
\quad\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\quad
\boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\quad \mathbf{A}\mathbf{v}\leq\mathbf{b},\\
& \mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star.
\end{aligned}
$$

[README.md:B048 paragraph]
The bounds, $\mathbf{A}$, and $\mathbf{b}$ are those of the chosen strategy.
With $\gamma=1$, the ranges cover only optimal flux vectors; with $\gamma<1$,
they also cover flux vectors that produce at least the fraction $\gamma$ of the optimum.
The last constraint is one more row of $\mathbf{A}\mathbf{v}\leq\mathbf{b}$:
append $-\mathbf{c}^{\top}$ to `A` and $-\gamma z^\star$ to `b`. Build the new
arrays with `vcat`; your Part 1 `b` may hold integers when the budget is an
integer. Because `solve_lp` only maximizes, find $v_j^{\min}$ by maximizing
$-v_j$ and changing the sign of the result. Take $z^\star$ as the `objective`
of the first solve and use $\gamma z^\star$ exactly. Do not lower it by a
tolerance: in this model, lowering it by only $10^{-9}$ mM/h can widen a range
that the optimum fixes by more than $10^{-5}$ mM/h. With the supplied solver,
every range solve in this assignment is feasible without a tolerance.

[README.md:B049 paragraph]
The driver and questions 3–5 call a supply flux **fixed** when
$v_j^{\max}-v_j^{\min}\leq10^{-5}$ mM/h. This is a reporting rule only;
`supply_ranges` returns the three fields below and no fixed flag.

[README.md:B050 paragraph]
PS3 applies FVA only to the 20 supply fluxes. Many internal reactions are
written as separate forward and reverse columns. Such a pair can carry equal
and opposite flux without changing any balance, so its range mostly reflects
its bounds rather than the optimum.

AUDIT THIS PART: README.md:B051 to README.md:B063

[README.md:B051 paragraph]
Implement **supply_ranges(model, budget; strategy=:equal, fraction=1.0)** in
[src/Compute.jl](src/Compute.jl), where `fraction` is $\gamma$. Solve once for
$z^\star$ with your Part 1 functions, then solve two LPs per supply column.
Return a vector containing one named tuple `(; column, minimum, maximum)` per
supply column, in the order of `model.amino_acid_uptake`. The `column` field is
the reaction-column index; the range endpoints are in mM/h. Reject an invalid
budget or strategy, or a `fraction` that is nonfinite or outside $[0,1]$, with
`ArgumentError`. If any solve is not optimal, throw an `ErrorException` naming
the budget, strategy, and status.

[README.md:B052 paragraph]
For a hand check, consider two supplies again, but now one unit of $v_p$
consumes one unit of each species, and a conversion flux $v_c\geq0$ turns the
second species into the first. The balances are $v_1+v_c=v_p$ and
$v_2=v_c+v_p$. With a budget of 6 and all other bounds loose, both strategies
reach $v_p=3$. Equal allocation fixes both supplies at **3**. Optimized
allocation fixes only their sum at 6, so $v_1$ ranges over **[0, 3]** and $v_2$
over **[3, 6]**. The optimum chooses the total supply but not how it is split.

[README.md:B053 paragraph]
Run the Part 3 checks and then the driver again:

[README.md:B054 code]
```bash
julia --startup-file=no --project=. testme_part_3.jl
julia --startup-file=no --project=. runproduction.jl
```

[README.md:B055 paragraph]
The driver adds these files to `outputs/`:

[README.md:B056 list]
- `supply-ranges.png`, `.svg`, and `.csv`: supply ranges at a budget of **1 mM/h**
  for equal allocation, optimized allocation, and optimized allocation with
  $\gamma=0.99$. The table has an `amino_acid` column and a minimum and maximum
  column for each case, such as `optimized_min_mM_per_h` and
  `optimized_99pct_max_mM_per_h`. The terminal reports how many supply fluxes
  are fixed in each case.
- `range-width-versus-budget.png`, `.svg`, and `.csv`: for each budget and
  strategy with $\gamma=1$, the sum of the 20 range widths
  (`equal_total_width_mM_per_h`, `optimized_total_width_mM_per_h`) and the number
  of fixed supply fluxes (`equal_fixed_supplies`, `optimized_fixed_supplies`).

[README.md:B057 paragraph]
Use these files to answer questions 3–5 in [responses.md](responses.md).
Unlike the optimal flux vectors, the range endpoints are unique, so the tests
compare them with reference values.

[README.md:B058 heading]
## Checking and submitting your work

[README.md:B059 paragraph]
Run the complete check after finishing the code and responses:

[README.md:B060 code]
```bash
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B061 paragraph]
The checker runs **24 tests per part, 72 total**, records their results and source
file fingerprints in `MANIFEST.txt`, and warns about unanswered discussion prompts.
If all three test suites pass, it generates the plots and data tables from your code.
The instructor still reviews the implementation and written explanations.

[README.md:B062 paragraph]
Zip the entire assignment folder, including your source, data, responses,
`MANIFEST.txt`, and generated files in `outputs/`. Name the archive
`CHEME-4800-5800-PS3-<your netid>.zip` and upload it to Canvas. The checker runs
locally; you must upload the ZIP yourself.

[README.md:B063 paragraph]
If some tests still fail at the deadline, submit your current work to receive
partial credit and participate in the infinite-revision policy.
