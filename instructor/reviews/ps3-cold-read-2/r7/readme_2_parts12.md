ALREADY READ: README.md up to the audited part

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

AUDIT THIS PART: README.md:B019 to README.md:B045

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
