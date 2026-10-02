ALREADY READ: the complete README.md

[README.md:B001 heading]
# Problem Set 3 (PS3): Feeding a cell-free protein factory

[README.md:B002 paragraph]
A cell-free protein synthesis system uses cellular machinery in a reaction
mixture to make a protein. That machinery needs amino acids, nucleotides, and
energy. The metabolic reactions in the mixture can make amino acids from other
available compounds, and amino acids can also be supplied directly.

[README.md:B003 paragraph]
You will use linear programming to find the maximum production rate of deGFP, a
green fluorescent protein, when the supply rates of the 20 amino acids share one
total budget, a limit in mM/h on their sum. With **equal allocation**, each amino
acid may use at most 1/20 of the budget. With **optimized allocation**, the solver
decides how to divide the budget among the amino acids. Your functions compute the
maximum production rate at each budget, and a supplied script plots it for both
strategies.

[README.md:B004 paragraph]
The same maximum can often be reached with different supply rates. You will use
flux variability analysis to find the smallest and largest value of each supply
rate among all solutions that reach the maximum.

[README.md:B005 paragraph]
This assignment extends the flux balance analysis in L6b. The supplied model is
adapted from the model described by **Vilkhovoy et al. (2018), “Sequence Specific
Modeling of E. coli Cell-Free Protein Synthesis,” ACS Synthetic Biology, 7(8),
1844–1857** ([paper](https://pubmed.ncbi.nlm.nih.gov/29944340/),
[DOI: 10.1021/acssynbio.7b00465](https://doi.org/10.1021/acssynbio.7b00465)).
The publication's code and model data are available in the Varner laboratory's
[sequence-specific cell-free protein synthesis repository](https://github.com/varnerlab/Sequence-Specific-FBA-CFPS-Publication-Code);
you do not need them for PS3.
The [data notes](data/README.md) identify the source and the PS3 operating conditions.

[README.md:B006 heading]
## Logistics

[README.md:B007 list]
- **Release:** Saturday, October 3, 2026. **Initial submission deadline:**
  **11:59 PM ET on Saturday, October 17, 2026**.
- **Submission:** Zip the assignment folder, name the archive
  `CHEME-4800-5800-PS3-<your netid>.zip` with your NetID in place of
  `<your netid>` (for example, `CHEME-4800-5800-PS3-abc123.zip`), and upload it to
  Canvas. The last section of this README
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

[README.md:B008 heading]
## Getting started

[README.md:B009 paragraph]
Use Julia `1.12.7`. The assignment contains its own data and source files.
Open the extracted assignment folder in VS Code. Run these commands from the
folder containing this README and [Project.toml](Project.toml):

[README.md:B010 code]
```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B011 paragraph]
The first command installs the required packages and can take several minutes.
The second runs the assignment tests. The four functions in
[src/Compute.jl](src/Compute.jl), described in Parts 1–3, initially raise errors because you must implement
them. Tests will fail until that work is complete.

[README.md:B012 paragraph]
Implement the four functions in `src/Compute.jl` and answer the five questions in
[responses.md](responses.md). You may add helper files in `src/`; if you do, load
them from [Include.jl](Include.jl). Keep all other supplied files unchanged.

[README.md:B013 heading]
## The model and its operating conditions

[README.md:B014 paragraph]
The model contains 146 species and 265 reactions. It includes metabolism,
transcription, translation, amino-acid supply, and protein output. The supplied
[load_model function](src/Model.jl), called as `model = load_model()`, reads the
stoichiometric matrix and sets the operating conditions. Its returned named tuple
contains:

[README.md:B015 table]
| Field | Meaning |
|---|---|
| `S` | Stoichiometric matrix; rows are species and columns are reactions. |
| `reactions`, `metabolites` | Identifiers in matrix column and row order. |
| `lower`, `upper` | One lower and upper flux bound per reaction, in mM/h. |
| `amino_acids` | Names of the 20 supplied amino acids. |
| `amino_acid_uptake` | Column indices of their supply reactions. |
| `protein` | Column index of the deGFP output reaction. |
| `translation` | Column index of translation initiation. Each initiation makes one deGFP, so this step's capacity caps protein production. |
| `translation_limit` | That capacity, in mM/h; `upper[translation]` already equals it. |

[README.md:B016 paragraph]
The `load_model` function already applies the operating conditions to `lower`
and `upper`: glucose supply of at most **1 mM/h**, oxygen supply of at most
**100 mM/h**, and amino-acid supply fluxes between 0 and 30 mM/h. Amino-acid
synthesis is enabled, while the reactions that degrade amino acids or export them
from the mixture are blocked. The transcription and translation bounds are
computed from fixed gene-expression parameters. Use these conditions for both
strategies; Part 1 describes how each strategy limits the amino-acid supply.

[README.md:B017 paragraph]
Supply reactions are written as `[] --> species` and have **nonnegative**
fluxes, so each amino-acid supply flux adds a nonnegative amount to the total that
the budget limits. The deGFP output reaction removes the protein that translation
makes, so its flux is the production rate. Every species, including deGFP,
satisfies a steady-state balance.

[README.md:B018 paragraph]
All model fluxes are in **mM/h**, millimoles per liter per hour. Protein
production is reported in **µM/h**, in function results, plots, and answers:
multiply the protein-output flux by 1000. The supply budget
is a rate limit; it does not specify an initial concentration or a monetary cost.
This static model predicts feasible production rates under fixed conditions.

[README.md:B019 heading]
## Part 1: Construct the objective and supply constraints

[README.md:B020 paragraph]
Let the vector $\mathbf{v}$ contain the reaction fluxes, and let $p$ be the
protein-output column, `model.protein`. The production problem is given by:

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
The vectors $\boldsymbol{\ell}$ and $\mathbf{u}$ specify the flux bounds. They
start from the model's `lower` and `upper`; equal allocation lowers some entries
of $\mathbf{u}$.
The matrix $\mathbf{A}$ and vector $\mathbf{b}$ hold any additional supply
inequality. Define $\mathcal{I}$ as the set of amino-acid supply columns,
`model.amino_acid_uptake`, and let $K=|\mathcal{I}|$ be their number. Let the
scalar $U\geq0$ be the total supply budget in mM/h, the `budget` argument of your
functions.

[README.md:B023 paragraph]
**Equal allocation.** Each amino acid may use at most one equal share of the
budget. The upper bound on its supply flux becomes the smaller of the model's
original bound $u_j$ and that share:

[README.md:B024 display]
$$
u_j^{\mathrm{equal}}=\min(u_j,U/K),\qquad j\in\mathcal{I}.
$$

[README.md:B025 paragraph]
The supply fluxes can be smaller than their limits. Any unused share remains
unavailable to the other amino acids. These individual limits already enforce
the total budget, so equal allocation adds no inequality: `A = zeros(0, n)`,
where `n` is the number of reactions, and `b = Float64[]`.

[README.md:B026 paragraph]
**Optimized allocation.** Retain the original individual bounds and add one
shared supply constraint, given by:

[README.md:B027 display]
$$
\sum_{j\in\mathcal{I}}v_j\leq U.
$$

[README.md:B028 paragraph]
The solver may distribute the available supply unequally. `A` is a single row,
a 1×n matrix, with a coefficient of one at each amino-acid supply column and zero
elsewhere, and `b = [U]`. Metabolic synthesis of amino acids remains
available under both strategies and does not count as external supply.

[README.md:B029 paragraph]
Implement two functions in [src/Compute.jl](src/Compute.jl):

[README.md:B030 list]
1. **protein_objective(model)** returns a `Vector{Float64}` with one coefficient
   per reaction. Its protein-output coefficient is one; the others are zero.
2. **allocation_constraints(model, budget; strategy=:equal)** returns a named tuple
   with fields `lower`, `upper`, `A`, and `b`. Support `:equal` and `:optimized`.
   Return copies of the bound arrays and leave `model.lower` and `model.upper`
   unchanged, so later calls start from the original bounds.
   Reject a negative or nonfinite budget or an unknown strategy with `ArgumentError`.

[README.md:B031 paragraph]
Use the model's supplied column indices. Your functions must also work on the
small test networks that the test scripts build. These have the same fields as the
model from `load_model` but different numbers and orders of reactions and
amino-acid supplies. All test models have nonnegative supply fluxes and at least
one supply column. Keep all lower bounds and unrelated upper bounds unchanged.

[README.md:B032 paragraph]
Run the Part 1 checks with:

[README.md:B033 code]
```bash
julia --startup-file=no --project=. testme_part_1.jl
```

[README.md:B034 paragraph]
For a hand check, consider two amino acids with supply fluxes $v_1$ and $v_2$.
One unit of protein output $v_p$ consumes two units of the first amino acid and
one unit of the second. The balances are $v_1=2v_p$ and $v_2=v_p$.
With a total supply budget of 6 and all other bounds loose, equal allocation
allows at most **1.5** units of protein output; optimized allocation allows
**2**. Explain these numbers to yourself before moving to the full model in Part 2.

[README.md:B035 heading]
## Part 2: Compute and interpret the production curves

[README.md:B036 paragraph]
Implement **production_curve(model, budgets)** in
[src/Compute.jl](src/Compute.jl). For each budget, construct both sets of
constraints and solve both LPs. Use your Part 1 functions and the supplied
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
The result contains `status`, `optimal`, `flux`, and `objective`. Check
`result.optimal` before reading the fluxes. If a solve is not optimal, throw an
`ErrorException` naming the budget, strategy, and status. In that case the solver
returns `nothing` for `flux` and `objective`.

[README.md:B039 paragraph]
Return a vector containing one named tuple per input budget, with fields
`budget`, `equal`, and `optimized`: the budget in mM/h and the maximum production
rate under each strategy in µM/h.
Preserve the input order, including repeated budgets. An empty input returns an
empty vector. Negative or nonfinite budgets raise `ArgumentError`.

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
- `outputs/production-versus-budget.png` and `.svg`: the two production curves
  and the translation capacity, the highest production rate that translation
  allows. The terminal prints this capacity in µM/h.
- `outputs/production-versus-budget.csv`: the numerical values used in the plot.
- `outputs/amino-acid-allocation.csv`: the supply fluxes of one optimal solution
  for each strategy at a total budget of **1 mM/h**. Because `production_curve`
  returns only rates, the driver solves these two problems again with your Part 1
  functions. It checks each solution against the balances, bounds, and budget,
  and prints its production rate and total supply.

[README.md:B044 paragraph]
Until you complete `supply_ranges`, the fourth function, described in Part 3, the
driver skips the Part 3 files (the supply-range plots and tables) and says so. Open the plot and use it with the CSV
tables to answer questions 1 and 2 in [responses.md](responses.md).

[README.md:B045 paragraph]
Your `amino-acid-allocation.csv` may differ from a classmate's. Either strategy
can have multiple equally productive flux distributions, and small differences in
how you build the same constraints can lead the solver to a different one. The
tests accept any solution that satisfies the constraints and gives the optimal
production rate. Part 3 measures how much those distributions differ.

[README.md:B046 heading]
## Part 3: Which supply fluxes does the optimum determine?

[README.md:B047 paragraph]
The optimal production rate of a linear program is unique, but the flux vector
that achieves it often is not. L6a introduced **flux variability analysis (FVA)**,
which finds the smallest and largest value of each flux among all flux vectors
that reach the optimum. FVA can also be run over a larger set: every flux vector
whose protein output is at least a fraction $\gamma$ of the optimum. FVA, including this fractional form,
was introduced by **Mahadevan and Schilling (2003),
“The effects of alternate optimal solutions in constraint-based genome-scale
metabolic models,” Metabolic Engineering, 5(4), 264–276**
([paper](https://pubmed.ncbi.nlm.nih.gov/14642354/),
[DOI: 10.1016/j.ymben.2003.09.002](https://doi.org/10.1016/j.ymben.2003.09.002)).

[README.md:B048 paragraph]
In PS3, $\gamma$ is the `fraction` keyword argument of `supply_ranges`, described
below. It must satisfy $0\leq\gamma\leq1$, and its default, $\gamma=1$, keeps only
optimal flux vectors.

[README.md:B049 paragraph]
PS3 applies FVA only to the amino-acid supply fluxes, 20 in the full model. Many
internal reactions are written as separate forward and reverse columns. Running
both columns of such a pair at the same rate forms a cycle that changes no
balance, so each column can rise to its upper bound no matter how much protein is
made. Their ranges would mostly reflect those bounds. Supply reactions have no
reverse partner, so their ranges are informative.

[README.md:B050 paragraph]
Let $\mathbf{c}$ be the protein objective. Let $z^\star$ be the maximum
protein-output flux for one allocation strategy at budget $U$, found by solving
the production LP. It is in mM/h, not the µM/h rate that `production_curve`
returns. For each supply column $j\in\mathcal{I}$, the range endpoints are
given by:

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
returns for the same strategy and budget. Steps 2 and 3 below show how to pose
these problems for `solve_lp`, which only maximizes and accepts only
$\mathbf{A}\mathbf{v}\leq\mathbf{b}$.

[README.md:B053 paragraph]
The driver and questions 3–5 call a supply flux **fixed** when
$v_j^{\max}-v_j^{\min}\leq10^{-5}$ mM/h, whatever $\gamma$ the ranges were
computed with. At $\gamma=1$, a fixed flux is one whose value the optimum
determines. This is a reporting rule only; `supply_ranges`, described
next, returns the range endpoints and no fixed flag.

[README.md:B054 paragraph]
Implement **supply_ranges(model, budget; strategy=:equal, fraction=1.0)** in
[src/Compute.jl](src/Compute.jl), where `fraction` is $\gamma$. It returns a
vector with one named tuple per supply column, in the order of
`model.amino_acid_uptake`. Each tuple has three fields:

[README.md:B055 list]
- `column`: the supply flux's column index in `S`;
- `minimum` and `maximum`: $v_j^{\min}$ and $v_j^{\max}$, in mM/h.

[README.md:B056 paragraph]
Compute the ranges in three steps:

[README.md:B057 list]
1. Build the strategy's limits and the protein objective $\mathbf{c}$ with your
   Part 1 functions, and solve once with `solve_lp`. Take $z^\star$ as this solve's
   `objective`, in mM/h, not the µM/h rate from `production_curve`.
2. Write the requirement $\mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star$ as
   $-\mathbf{c}^{\top}\mathbf{v}\leq-\gamma z^\star$ and add it as one more row:
   `vcat(A, -c')` and `vcat(b, -fraction * z)`, where `z` holds $z^\star$. Use
   `vcat`, which returns new arrays. `push!` cannot add a row to `A`, and it fails
   when an integer budget has made your Part 1 `b` a `Vector{Int}`; Part 1 does
   not need to change.
3. For each supply column $j$, solve twice with the strategy's `lower` and
   `upper` from step 1 and the new `A` and `b`. For $v_j^{\max}$, maximize with an objective that is one at column
   $j$ and zero elsewhere. Because `solve_lp` only maximizes, find $v_j^{\min}$ by
   maximizing with an objective that is minus one at column $j$ and negating that
   solve's `objective`.

[README.md:B058 paragraph]
Some FVA codes lower $\gamma z^\star$ by a small tolerance, because round-off can
make the exact requirement slightly unreachable. Do not. Near the optimum, giving up a tiny amount of production can allow
a large change in some supply fluxes: lowering $\gamma z^\star$ by only
$10^{-9}$ mM/h can turn a fixed supply flux into one whose range is wider than
$10^{-5}$ mM/h. With the supplied solver, every solve
in step 3 is feasible at exactly $\gamma z^\star$.

[README.md:B059 paragraph]
Reject an invalid budget or strategy, as in Part 1, or a `fraction` that is
nonfinite or outside $[0,1]$, with `ArgumentError`. If the solve in step 1 or any
solve in step 3 is not optimal, throw an `ErrorException` naming the budget,
strategy, and status.

[README.md:B060 paragraph]
For a hand check, consider the two amino acids again, but now one unit of $v_p$
consumes one unit of each species, and a conversion flux $v_c\geq0$ turns the
second species into the first. The balances are $v_1+v_c=v_p$ and
$v_2=v_c+v_p$. With a budget of 6 and all other bounds loose, both strategies
reach $v_p=3$. Equal allocation fixes both supplies at **3**. Optimized
allocation fixes only their sum at 6, so $v_1$ ranges over **[0, 3]** and $v_2$
over **[3, 6]**. The optimum chooses the total supply but not how it is split.

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
- `supply-ranges.png`, `.svg`, and `.csv`: supply ranges at a budget of **1 mM/h**
  for equal allocation, optimized allocation, and optimized allocation with
  $\gamma=0.99$, which shows how the ranges widen when production may fall 1%
  below the optimum. The table has an `amino_acid` column and, for each case, a
  minimum and maximum column: `equal_`, `optimized_`, or `optimized_99pct_`
  followed by `min_mM_per_h` or `max_mM_per_h`. The terminal reports how many supply fluxes
  are fixed in each case.
- `range-width-versus-budget.png`, `.svg`, and `.csv`: one row per budget from 0
  to 4 mM/h, for both strategies at $\gamma=1$. The columns
  `equal_total_width_mM_per_h` and `optimized_total_width_mM_per_h` hold the sum
  of the 20 range widths $v_j^{\max}-v_j^{\min}$. The columns
  `equal_fixed_supplies` and `optimized_fixed_supplies` count the fixed supply
  fluxes.

[README.md:B065 paragraph]
Use these files to answer questions 3–5 in [responses.md](responses.md).
Unlike the optimal flux vectors, the range endpoints are unique, so the tests
compare them with reference values.

[README.md:B066 heading]
## Checking and submitting your work

[README.md:B067 paragraph]
Run the complete check after finishing the code and responses:

[README.md:B068 code]
```bash
julia --startup-file=no --project=. check_submission.jl
```

[README.md:B069 paragraph]
The checker runs **24 tests per part, 72 total**. It records their results in
`MANIFEST.txt`, with a SHA-256 fingerprint of each source file so the instructor
can match the results to the submitted code, and it warns about unanswered
questions in `responses.md`.
If all three test suites pass, it also regenerates the plots and data tables in
`outputs/` from your code. Otherwise it leaves `outputs/` as your last
`runproduction.jl` run wrote it; in that case, run the driver after your final
edits and then run the checker last.
The instructor still reviews the implementation and written explanations.

[README.md:B070 paragraph]
Zip the entire assignment folder, including your source, data, responses,
`MANIFEST.txt`, and generated files in `outputs/`. Name the archive
`CHEME-4800-5800-PS3-<your netid>.zip` and upload it to Canvas. The checker runs
locally; you must upload the ZIP yourself.

[README.md:B071 paragraph]
If some tests still fail at the deadline, submit your current work to receive
partial credit and participate in the infinite-revision policy.

FILES THE STUDENT HAS OPEN: after running the driver, the student has these generated
tables in outputs/ (header rows shown; each table has one row per budget or per amino acid):
- amino-acid-allocation.csv: amino_acid,equal_uptake_mM_per_h,optimized_uptake_mM_per_h,equal_limit_mM_per_h
- production-versus-budget.csv: budget_mM_per_h,equal_uM_per_h,optimized_uM_per_h
- range-width-versus-budget.csv: budget_mM_per_h,equal_total_width_mM_per_h,optimized_total_width_mM_per_h,equal_fixed_supplies,optimized_fixed_supplies
- supply-ranges.csv: amino_acid,equal_min_mM_per_h,equal_max_mM_per_h,optimized_min_mM_per_h,optimized_max_mM_per_h,optimized_99pct_min_mM_per_h,optimized_99pct_max_mM_per_h
The student has also viewed the generated plots described in the README.

AUDIT THIS PART: responses.md:B001 to responses.md:B012

[responses.md:B001 heading]
# PS3 written responses

[responses.md:B002 paragraph]
Keep the five numbered questions and replace each `TODO` line with your answer.
Answer every lettered part and give the requested numbers with units. A few
sentences per part is enough.

[responses.md:B003 list]
1. **Why can optimized allocation match or improve protein production?**
   - (a) Compare the feasible sets of the two strategies at the same budget.
   - (b) Explain why equal supply limits do not require equal actual uptake.
   - (c) Explain why each strategy's problem is a linear program, even though the
     equal-allocation bound min(u_j, U/K) takes a minimum.
   - (d) From `production-versus-budget.csv` at a budget of 1.0 mM/h, report both
     production rates and the percentage improvement of optimized over equal
     allocation.
   - (e) Add up each uptake column of `amino-acid-allocation.csv` and compare the
     total supply each strategy uses with the 1.0 mM/h budget. Name the amino
     acids whose `equal_uptake_mM_per_h` is below `equal_limit_mM_per_h`
     (0.05 mM/h) by more than 10⁻⁵ mM/h. These equal-allocation values are the
     same in every optimal solution at this budget.

[responses.md:B004 paragraph]
   TODO: Write your response.

[responses.md:B005 list]
2. **Why do the curves level off?**
   - (a) Using `production-versus-budget.csv` and the translation capacity printed
     by the driver, identify the first sampled budget at which each strategy
     comes within 0.001 µM/h of that capacity.
   - (b) Explain why relaxing a supply constraint cannot reduce the optimal
     production rate.
   - (c) Explain why a larger amino-acid budget eventually stops raising
     production in this model. Which model bound would you relax next to raise
     the plateau, and which constraint might limit production after that? Reason
     from the operating conditions in the README; no new computation is needed.
   - (d) Does this steady-state calculation determine the final protein
     concentration after several hours? Explain.

[responses.md:B006 paragraph]
   TODO: Write your response.

[responses.md:B007 list]
3. **What does the optimum determine?** Use the `optimized_min_mM_per_h` and
   `optimized_max_mM_per_h` columns of `supply-ranges.csv` (budget 1.0 mM/h, γ = 1).
   - (a) A supply is fixed when its range width is at most 10⁻⁵ mM/h. List the
     amino acids that are never supplied externally (both endpoints within
     10⁻⁵ mM/h of zero) and those whose supply can vary.
   - (b) Name any two amino acids whose supply is fixed above the equal share
     U/K = 0.05 mM/h.
   - (c) The optimized uptake in `amino-acid-allocation.csv` is one optimal
     solution. Does each value lie within its range, and which values could a
     different optimal solution change?
   - (d) Explain how internal synthesis allows different supply vectors to give
     the same production rate.

[responses.md:B008 paragraph]
   TODO: Write your response.

[responses.md:B009 list]
4. **Are the ranges independent?** Parts (a) and (b) use the same columns as
   question 3.
   - (a) Add the 20 values of `optimized_max_mM_per_h` and compare the total with
     the budget. Can every supply sit at its maximum at the same time?
   - (b) Subtract the sum of the fixed supplies from the budget. How much can the
     supplies that vary use at most? Suppose glutamate, one of them, is supplied
     at its maximum. Using this budget arithmetic and the minimum of each range,
     show which value, in mM/h, each of the other varying supplies is then
     forced to take.
   - (c) Add the 20 values of `optimized_99pct_max_mM_per_h` (γ = 0.99) and
     compare the total with the budget.
   - (d) Explain why the 20 ranges, taken together, do not tell you which
     combinations of supply values reach the required production, cᵀv ≥ γz*.

[responses.md:B010 paragraph]
   TODO: Write your response.

[responses.md:B011 list]
5. **How far can you trust one flux vector?**
   - (a) Use `range-width-versus-budget.csv` to identify the first sampled budget
     at which each strategy has no fixed supply flux. Compare these budgets with
     your answer to question 2(a) and explain the relationship.
   - (b) Using the `optimized_99pct` columns of `supply-ranges.csv` (optimized
     allocation, 1.0 mM/h, γ = 0.99), report how many supply fluxes remain fixed.
   - (c) Report valine's `optimized_99pct_min_mM_per_h` and how far it falls
     below valine's fixed optimized supply at γ = 1, in mM/h.
   - (d) FVA ranges are sometimes described as flux uncertainty. Explain what
     these ranges represent and what they do not.

[responses.md:B012 paragraph]
   TODO: Write your response.
