ALREADY READ: README.md up to the audited part

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

AUDIT THIS PART: README.md:B019 to README.md:B045

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
