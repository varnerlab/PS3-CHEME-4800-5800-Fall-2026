# Problem Set 3 (PS3): Feeding a cell-free protein factory

A cell-free protein synthesis system uses cellular machinery in a reaction
mixture to make a protein. That machinery needs amino acids, nucleotides, and
energy. Metabolic reactions in the mixture can make amino acids from other
compounds, and amino acids can also be supplied directly.

You will use linear programming to find the maximum production rate of deGFP, a
green fluorescent protein. The 20 amino-acid supply rates share one budget, a
limit in mM/h on their sum. With **equal allocation**, each amino acid's supply is
capped at 1/20 of the budget and may be less. With **optimized allocation**, the
solver divides the budget freely. The supplied script `runproduction.jl` plots
the maximum rate against the budget for both strategies.

The same maximum can often be reached with different supply rates. You will then
use flux variability analysis to find the smallest and largest value each supply
rate can take among all solutions that reach the maximum.

The model is adapted from **Vilkhovoy et al. (2018), “Sequence Specific
Modeling of E. coli Cell-Free Protein Synthesis,” ACS Synthetic Biology, 7(8),
1844–1857** ([paper](https://pubmed.ncbi.nlm.nih.gov/29944340/),
[DOI: 10.1021/acssynbio.7b00465](https://doi.org/10.1021/acssynbio.7b00465)).
The publication's code and data are in the Varner laboratory's
[sequence-specific cell-free protein synthesis repository](https://github.com/varnerlab/Sequence-Specific-FBA-CFPS-Publication-Code);
you do not need them. The [data notes](data/README.md) describe the source and
the PS3 operating conditions.

## Logistics

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

## Getting started

Use Julia `1.12.7`; the folder contains all the data and source files you need.
Open it in VS Code and run these commands from the folder containing this README
and [Project.toml](Project.toml):

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

The first command installs the required packages and can take several minutes.
The second runs the tests, which fail until you implement the four functions in
[src/Compute.jl](src/Compute.jl) described in Parts 1–3.

Implement those functions and answer the five questions in
[responses.md](responses.md). You may add helper files in `src/` and load them
from [Include.jl](Include.jl). Keep all other supplied files unchanged.

## The model and its operating conditions

The model has 146 species and 265 reactions covering metabolism, transcription,
translation, amino-acid supply, and protein output. The supplied
[load_model function](src/Model.jl), called as `model = load_model()`, returns a
named tuple with these fields:

| Field | Meaning |
|---|---|
| `S` | Stoichiometric matrix; rows are species and columns are reactions. |
| `reactions`, `metabolites` | Identifiers in column and row order. |
| `lower`, `upper` | Lower and upper flux bound of each reaction, in mM/h. |
| `amino_acids` | Names of the 20 supplied amino acids. |
| `amino_acid_uptake` | Column indices of their supply reactions. |
| `protein` | Column index of the deGFP output reaction. |
| `translation` | Column index of translation initiation; in this model, each initiation yields one deGFP. |
| `translation_limit` | Translation capacity in mM/h, already set as `upper[translation]`. |

`load_model` already applies the operating conditions: glucose supply of at most
**1 mM/h**, oxygen supply of at most **100 mM/h**, and each amino-acid supply
between 0 and 30 mM/h. Amino-acid synthesis is enabled; reactions that degrade amino
acids or export them from the mixture are blocked. Transcription and translation
bounds come from fixed gene-expression parameters. Both strategies start from
these conditions; Part 1 describes how each one limits amino-acid supply.

Amino-acid supply reactions are written as `[] --> species` and have
**nonnegative** fluxes, so the budget can simply cap their sum. The deGFP output
reaction removes the protein that translation makes; its flux is the production
rate. Every species, including deGFP, is at steady state.

Model fluxes are in **mM/h**, millimoles per liter per hour. Protein production
rates are reported in **µM/h** in plots, answers, and the Part 2 function results:
multiply the protein-output flux by 1000. The supply budget is a rate limit, not an initial
concentration or a monetary cost. This static model predicts feasible production
rates under fixed conditions.

## Part 1: Construct the objective and supply constraints

Let $\mathbf{v}$ contain the reaction fluxes and let $p$ be the protein-output
column, `model.protein`. The production problem is:

$$
\begin{aligned}
\underset{\mathbf{v}}{\operatorname{maximize}}\quad & v_p\\
\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\\
& \boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\\
& \mathbf{A}\mathbf{v}\leq\mathbf{b}.
\end{aligned}
$$

The bounds $\boldsymbol{\ell}$ and $\mathbf{u}$ start from the model's `lower` and
`upper`, and $\mathbf{A}$ and $\mathbf{b}$ hold any extra supply inequality. Let
$\mathcal{I}$ be the set of amino-acid supply columns, `model.amino_acid_uptake`,
and $K=|\mathcal{I}|$ their number. The scalar $U\geq0$ is the total supply budget
in mM/h, the `budget` argument of your functions.

**Equal allocation.** Each supply flux's upper bound becomes the smaller of its
original bound $u_j$ and an equal share of the budget:

$$
u_j^{\mathrm{equal}}=\min(u_j,U/K),\qquad j\in\mathcal{I}.
$$

Supply fluxes can stay below their limits, and an unused share is not available
to other amino acids. These limits already enforce the budget, so equal
allocation adds no inequality: `A = zeros(0, n)`, where `n` is the number of
reactions, and `b = Float64[]`.

**Optimized allocation.** Keep the original bounds and add one shared supply
constraint:

$$
\sum_{j\in\mathcal{I}}v_j\leq U.
$$

The solver may divide the supply unequally. `A` is a 1×n matrix with ones at the
amino-acid supply columns and zeros elsewhere, and `b = [U]`. Under both
strategies, amino acids made by metabolism do not count as external supply.

Implement two functions in [src/Compute.jl](src/Compute.jl):

1. **protein_objective(model)** returns the objective vector $\mathbf{c}$, with
   $\mathbf{c}^{\top}\mathbf{v}=v_p$: a `Vector{Float64}` with 1 at the
   protein-output column and 0 elsewhere.
2. **allocation_constraints(model, budget; strategy=:equal)** returns a named tuple
   with fields `lower`, `upper`, `A`, and `b` for `:equal` or `:optimized`.
   Return copies of the bound arrays so that `model.lower` and `model.upper` stay
   unchanged. Reject a negative or nonfinite budget or an unknown strategy with
   `ArgumentError`.

Use the model's column indices: the test scripts also build small networks with
the same fields but different numbers and orders of reactions and supplies. All
test networks have nonnegative supply fluxes and at least one supply column. In
the returned bounds, only the supply upper bounds differ from the model's, and
only under `:equal`.

Run the Part 1 checks with:

```bash
julia --startup-file=no --project=. testme_part_1.jl
```

For a hand check, take two amino acids with supply fluxes $v_1$ and $v_2$. One
unit of protein output $v_p$ consumes two units of the first and one of the
second, so $v_1=2v_p$ and $v_2=v_p$. With a total supply budget of 6 and all other
bounds loose, equal allocation allows at most **1.5** units of protein output and
optimized allocation allows **2**. Make sure you can explain these numbers.

## Part 2: Compute and interpret the production curves

Implement **production_curve(model, budgets)** in
[src/Compute.jl](src/Compute.jl). For each budget, build both strategies'
constraints with your Part 1 functions and solve each with the supplied
[solve_lp function](src/Solver.jl), which maximizes $\mathbf{c}^{\top}\mathbf{v}$
subject to the balances, bounds, and $\mathbf{A}\mathbf{v}\leq\mathbf{b}$. A single
solve has this form:

```julia
limits = allocation_constraints(model, budget; strategy=:equal);
c = protein_objective(model);
result = solve_lp(model.S, limits.lower, limits.upper, c;
    A=limits.A, b=limits.b);
```

The result contains `status`, `optimal`, `flux`, and `objective`. When a solve is
not optimal, `flux` and `objective` are `nothing`, so check `result.optimal` first
and throw an `ErrorException` naming the budget, strategy, and status.

Return one named tuple per input budget, in input order and including repeats,
with fields `budget` (mM/h) and `equal` and `optimized` (maximum production rates
in µM/h). An empty input returns an empty vector; a negative or nonfinite budget
raises `ArgumentError`.

After completing the function, run:

```bash
julia --startup-file=no --project=. testme_part_2.jl
julia --startup-file=no --project=. runproduction.jl
```

The supplied driver, `runproduction.jl`, evaluates budgets from **0 to 4 mM/h in
steps of 0.1** and writes these files:

- `outputs/production-versus-budget.png` and `.svg`: both production curves and
  the translation capacity, the highest rate translation allows. The terminal
  prints this capacity in µM/h.
- `outputs/production-versus-budget.csv`: the plotted values.
- `outputs/amino-acid-allocation.csv`: the supply fluxes of one optimal solution
  per strategy at a budget of **1 mM/h**. Because `production_curve` returns only
  rates, the driver solves both strategies again with your Part 1 functions. It
  checks each solution against the balances, bounds, and budget and prints its
  production rate and total supply.

Until `supply_ranges` (Part 3) is complete, the driver skips the Part 3 files and
says so. Use the plot and CSV tables to answer questions 1 and 2 in
[responses.md](responses.md).

Your `amino-acid-allocation.csv` may differ from a classmate's: optimal solutions
are often not unique, and the tests accept any that satisfies the constraints and
gives the optimal production rate. Part 3 measures how much they can differ.

## Part 3: Which supply fluxes does the optimum determine?

The optimal production rate of a linear program is unique, but the flux vector
that achieves it often is not. **Flux variability analysis (FVA)** finds the
smallest and largest value of each flux among all flux vectors that reach the
optimum or, more generally, at least a fraction $\gamma$ of it. FVA, including
this fractional form, was introduced by **Mahadevan and Schilling (2003),
“The effects of alternate optimal solutions in constraint-based genome-scale
metabolic models,” Metabolic Engineering, 5(4), 264–276**
([paper](https://pubmed.ncbi.nlm.nih.gov/14642354/),
[DOI: 10.1016/j.ymben.2003.09.002](https://doi.org/10.1016/j.ymben.2003.09.002)).

In PS3, $\gamma$ is the `fraction` argument of `supply_ranges`, described below, with
$0\leq\gamma\leq1$. The default, $\gamma=1$, keeps only optimal flux vectors;
$\gamma=0.99$ shows how far the supplies can shift if production may fall 1%.

PS3 applies FVA only to the amino-acid supply fluxes. Many internal reactions
have separate forward and reverse columns. Running both columns of a pair at the
same rate forms a cycle that changes no balance, so each can rise to its upper
bound and its range says little. Supply reactions have no reverse partner.

Let $\mathbf{c}$ be the protein objective and $z^\star$ the maximum
protein-output flux for one strategy at budget $U$, found by solving the
production problem once. It is in mM/h, not µM/h.
For each supply column $j\in\mathcal{I}$, the range endpoints are:

$$
\begin{aligned}
v_j^{\min}=\min_{\mathbf{v}}\ v_j,\qquad v_j^{\max}=\max_{\mathbf{v}}\ v_j
\quad\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\quad
\boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\quad \mathbf{A}\mathbf{v}\leq\mathbf{b},\\
& \mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star.
\end{aligned}
$$

The bounds, $\mathbf{A}$, and $\mathbf{b}$ are those that `allocation_constraints`
returns for the same strategy and budget. The steps below show how to pose these
problems for `solve_lp`.

A supply flux is **fixed** when $v_j^{\max}-v_j^{\min}\leq10^{-5}$ mM/h. At
$\gamma=1$, this means the optimum determines its value. The driver and questions
3–5 use this rule; `supply_ranges`, described next, returns only the endpoints.

Implement **supply_ranges(model, budget; strategy=:equal, fraction=1.0)** in
[src/Compute.jl](src/Compute.jl), where `fraction` is $\gamma$. Return one named
tuple per supply column, in the order of `model.amino_acid_uptake`, with fields:

- `column`: the column index in `S`;
- `minimum` and `maximum`: $v_j^{\min}$ and $v_j^{\max}$, in mM/h.

Because `solve_lp` only maximizes and accepts only
$\mathbf{A}\mathbf{v}\leq\mathbf{b}$, compute the ranges in three steps:

1. Build the strategy's limits and $\mathbf{c}$ with your Part 1 functions and
   solve once. Take $z^\star$ as this solve's `objective`, in mM/h.
2. Write $\mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star$ as
   $-\mathbf{c}^{\top}\mathbf{v}\leq-\gamma z^\star$ and append it as one row:
   `vcat(A, -c')` and `vcat(b, -fraction * z)`, where `z` holds $z^\star$. `vcat`
   also converts a Part 1 `b` that holds integers to `Float64`, so Part 1 needs
   no change.
3. For each supply column $j$, solve twice with the step 1 bounds and the new `A`
   and `b`. For $v_j^{\max}$, maximize an objective vector that is one at column
   $j$ and zero elsewhere. For $v_j^{\min}$, maximize one that is minus one at
   column $j$, and negate that solve's `objective`.

Use $\gamma z^\star$ exactly; do not lower it by a small tolerance, as some FVA
codes do. Near the optimum, a tiny loss of production can allow a large shift in
supply: lowering $\gamma z^\star$ by only $10^{-9}$ mM/h can make a fixed supply
flux's range wider than $10^{-5}$ mM/h. With the supplied solver, every step 3
solve is feasible at exactly $\gamma z^\star$.

Throw `ArgumentError` for a negative or nonfinite budget, an unknown strategy, or
a `fraction` that is nonfinite or outside $[0,1]$. If any solve is not optimal,
throw an `ErrorException` naming the budget, strategy, and status.

For a hand check, return to the two amino acids, but now one unit of $v_p$
consumes one unit of each, and a conversion flux $v_c\geq0$ turns the second into
the first, so $v_1+v_c=v_p$ and $v_2=v_c+v_p$. Adding these gives
$v_1+v_2=2v_p$, so with a budget of 6 and all other bounds loose, both strategies
reach $v_p=3$. Equal allocation fixes both supplies
at **3**. Optimized allocation fixes only their sum at 6, so $v_1$ ranges over
**[0, 3]** and $v_2$ over **[3, 6]**. The optimum sets the total supply but not
its split.

Run the Part 3 checks and then the driver again:

```bash
julia --startup-file=no --project=. testme_part_3.jl
julia --startup-file=no --project=. runproduction.jl
```

The driver adds these files to `outputs/`:

- `supply-ranges.png`, `.svg`, and `.csv`: supply ranges at **1 mM/h** for equal
  allocation, optimized allocation, and optimized allocation with $\gamma=0.99$.
  The CSV's columns are `amino_acid`, `equal_min_mM_per_h`,
  `equal_max_mM_per_h`, `optimized_min_mM_per_h`, `optimized_max_mM_per_h`,
  `optimized_99pct_min_mM_per_h`, and `optimized_99pct_max_mM_per_h`. The
  terminal reports how many supply fluxes are fixed in each case.
- `range-width-versus-budget.png`, `.svg`, and `.csv`: one row per budget from 0
  to 4 mM/h, both strategies at $\gamma=1$. `equal_total_width_mM_per_h` and
  `optimized_total_width_mM_per_h` sum the 20 range widths
  $v_j^{\max}-v_j^{\min}$; `equal_fixed_supplies` and `optimized_fixed_supplies`
  count the fixed supply fluxes.

Use these files to answer questions 3–5 in [responses.md](responses.md). Unlike
optimal flux vectors, range endpoints are unique, so the tests compare them with
reference values.

## Checking and submitting your work

Run the complete check after finishing the code and responses:

```bash
julia --startup-file=no --project=. check_submission.jl
```

The checker runs **24 tests per part, 72 total**. It records the results, with a
SHA-256 fingerprint of each source file, in `MANIFEST.txt` and warns about
unanswered questions in `responses.md`. If all three test suites pass, it also
regenerates `outputs/` from your code. If not, run the driver yourself after your
final edits so `outputs/` is current, then run the checker last. For grading, the
instructor reruns the checker on your code; [RUBRIC.md](RUBRIC.md) explains how
the score follows from the result.

Zip the entire assignment folder, including your source, data, responses,
`MANIFEST.txt`, and `outputs/`. Name the archive
`CHEME-4800-5800-PS3-<your netid>.zip` and upload it to Canvas yourself; the
checker runs only locally.

If some tests still fail at the deadline, submit your current work to receive
partial credit and participate in the infinite-revision policy.
