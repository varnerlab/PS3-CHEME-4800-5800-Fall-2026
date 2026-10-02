ALREADY READ: the complete README.md

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
Keep the five numbered questions and replace each `TODO` line with your answers
to every lettered part, giving the requested numbers with units. A few sentences
per part is enough.

[responses.md:B003 list]
1. **Why can optimized allocation match or improve protein production?**
   - (a) Compare the feasible sets of the two strategies at the same budget.
   - (b) Explain why equal supply limits do not require equal actual uptake.
   - (c) Explain why each strategy is a linear program, even though the
     equal-allocation bound min(u_j, U/K) takes a minimum.
   - (d) From `production-versus-budget.csv` at 1.0 mM/h, report both production
     rates and the percentage improvement of optimized over equal allocation.
   - (e) Sum each uptake column of `amino-acid-allocation.csv` and compare each
     total with the 1.0 mM/h budget. Name the amino acids whose
     `equal_uptake_mM_per_h` is more than 10⁻⁵ mM/h below `equal_limit_mM_per_h`
     (0.05 mM/h); these values are the same in every optimal solution.

[responses.md:B004 paragraph]
   TODO: Write your response.

[responses.md:B005 list]
2. **Why do the curves level off?**
   - (a) Using `production-versus-budget.csv` and the translation capacity printed
     by the driver, find the first sampled budget at which each strategy comes
     within 0.001 µM/h of the capacity.
   - (b) Explain why relaxing a supply constraint cannot reduce the optimal
     production rate.
   - (c) Explain why a larger budget eventually stops raising production. Which
     model bound would you relax next to raise the plateau, and which constraint
     might then limit production? Reason from the README's operating conditions;
     no new computation is needed.
   - (d) Does this steady-state calculation determine the final protein
     concentration after several hours? Explain.

[responses.md:B006 paragraph]
   TODO: Write your response.

[responses.md:B007 list]
3. **What does the optimum determine?** Use the `optimized_min_mM_per_h` and
   `optimized_max_mM_per_h` columns of `supply-ranges.csv` (1.0 mM/h, γ = 1). A
   supply is fixed when its range width is at most 10⁻⁵ mM/h.
   - (a) List the amino acids never supplied externally (both endpoints within
     10⁻⁵ mM/h of zero) and those whose supply can vary.
   - (b) Name any two amino acids whose supply is fixed above the equal share
     U/K = 0.05 mM/h.
   - (c) The optimized uptake in `amino-acid-allocation.csv` is one optimal
     solution. Does each value lie within its range, and which values could
     another optimal solution change?
   - (d) Explain how internal synthesis allows different supply vectors to give
     the same production rate.

[responses.md:B008 paragraph]
   TODO: Write your response.

[responses.md:B009 list]
4. **Are the ranges independent?** Parts (a) and (b) use the same columns as
   question 3.
   - (a) Sum the 20 values of `optimized_max_mM_per_h` and compare the total with
     the budget. Can every supply sit at its maximum at once?
   - (b) Subtract the sum of the fixed supplies from the budget: how much can the
     varying supplies use at most? If glutamate, one of them, is supplied at its
     maximum, use this arithmetic and each range's minimum to find the value, in
     mM/h, that each other varying supply is forced to take.
   - (c) Sum the 20 values of `optimized_99pct_max_mM_per_h` (γ = 0.99) and
     compare the total with the budget.
   - (d) Explain why the 20 ranges, taken together, do not tell you which
     combinations of supply values reach the required production, cᵀv ≥ γz*.

[responses.md:B010 paragraph]
   TODO: Write your response.

[responses.md:B011 list]
5. **How far can you trust one flux vector?**
   - (a) From `range-width-versus-budget.csv`, find the first sampled budget at
     which each strategy has no fixed supply flux. Compare these budgets with your
     answer to 2(a) and explain the relationship.
   - (b) From the `optimized_99pct` columns of `supply-ranges.csv` (optimized
     allocation, 1.0 mM/h, γ = 0.99), report how many supply fluxes remain fixed.
   - (c) Report valine's `optimized_99pct_min_mM_per_h` and how far it falls
     below valine's fixed optimized supply at γ = 1, in mM/h.
   - (d) FVA ranges are sometimes described as flux uncertainty. Explain what
     these ranges represent and what they do not.

[responses.md:B012 paragraph]
   TODO: Write your response.
