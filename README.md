# Problem Set 3 (PS3): Feeding a cell-free protein factory

A cell-free protein synthesis system makes protein using cellular machinery in a
reaction mixture, outside living cells. The machinery needs amino acids, which
the mixture can make from other compounds or receive as a direct supply.

You will use linear programming to find the maximum production rate of deGFP, a
green fluorescent protein, in an _E. coli_ cell-free system. The 20 amino-acid
supply rates share one **budget**: a limit, in mM/h, on their sum. You will
compare two ways to use that budget:

- **Equal allocation:** each amino acid may use at most 1/20 of the budget.
- **Optimized allocation:** the solver divides the budget among the amino acids
  to maximize production.

**Due: 11:59 PM ET, Saturday, October 17, 2026.** Released October 3, 2026.

## Choose your track

[TRACK.txt](TRACK.txt) holds your track, `standard` or `advanced`. Leave it as
`standard` until you start Part 3. The checker and `runproduction.jl` read this
file to decide what to run.

| Track | Code in [src/Compute.jl](src/Compute.jl) | Questions in [responses.md](responses.md) | Credit |
|:--|:--|:--|:--|
| **Standard** | Parts 1–2: `protein_objective`, `allocation_constraints`, `production_curve` | 1 and 2 | Up to 4 points |
| **Advanced** | Parts 1–3: the three Standard functions plus `supply_ranges` | 1, 2, and 3 | Up to 4 points plus one Magic Point |

On either track, Parts 1–2 and questions 1–2 set your score out of 4, so
choosing Advanced never lowers it. The Magic Point requires a score of 4,
all Part 3 tests passing, a satisfactory answer to question 3, and
teaching-team review. You may change tracks in a later revision. See
[RUBRIC.md](RUBRIC.md) for the grading scale.

On the Standard track, leave `supply_ranges` and the question 3 `TODO`
unchanged.

## Setup

Use Julia `1.12.7`. Open this folder in VS Code and run these commands from it:

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

The first command installs the packages and takes several minutes. The second
runs the 48 Standard tests, which fail until you finish Parts 1–2.

Edit only [TRACK.txt](TRACK.txt), [src/Compute.jl](src/Compute.jl),
[responses.md](responses.md), and, if you add helper files in `src/`,
[Include.jl](Include.jl).

## The model

The model has 146 species and 265 reactions covering metabolism, transcription,
translation, amino-acid supply, and protein output. Calling
`model = load_model()` ([src/Model.jl](src/Model.jl)) returns a named tuple:

| Field | Meaning |
|---|---|
| `S` | Stoichiometric matrix; rows are species and columns are reactions. |
| `reactions`, `metabolites` | Identifiers in column and row order. |
| `lower`, `upper` | Lower and upper flux bound of each reaction, in mM/h. |
| `amino_acids` | Names of the 20 supplied amino acids. |
| `amino_acid_uptake` | Column indices of their supply reactions. |
| `protein` | Column index of the deGFP output reaction. |
| `translation` | Column index of translation initiation; each initiation yields one deGFP. |
| `translation_limit` | Translation capacity in mM/h, already set as `upper[translation]`. |

The `load_model(...)` function already applies the operating conditions: glucose supply of at most
1 mM/h, oxygen supply of at most 100 mM/h, and each amino-acid supply between 0
and 30 mM/h. The mixture can make amino acids but cannot degrade or export them.

All fluxes are in **mM/h**. Amino-acid supply reactions (`[] --> species`) have
nonnegative fluxes, so the budget can simply cap their sum. The flux of the
deGFP output reaction is the production rate; report it in **µM/h** by
multiplying by 1000.

## Part 1 (Standard): Objective and supply constraints

Let $\mathbf{v}$ be the reaction fluxes and $p$ the protein-output column
(`model.protein`). The production problem is the linear program

$$
\begin{aligned}
\underset{\mathbf{v}}{\text{maximize}}\quad & v_p\\
\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\\
& \boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\\
& \mathbf{A}\mathbf{v}\leq\mathbf{b}.
\end{aligned}
$$

Here $\boldsymbol{\ell}$ and $\mathbf{u}$ are `model.lower` and `model.upper`,
and $\mathbf{A}\mathbf{v}\leq\mathbf{b}$ holds any extra supply constraint. Let
$\mathcal{I}$ be the supply columns (`model.amino_acid_uptake`), $K=|\mathcal{I}|$
their number, and $U\geq0$ the budget in mM/h. The budget counts only amino
acids supplied from outside the mixture, not those it makes.

**Equal allocation** caps each supply at its share of the budget and keeps any
tighter original bound:

$$
u_j^{\mathrm{equal}}=\min(u_j,U/K),\qquad j\in\mathcal{I}.
$$

There is no extra constraint, so return `A = zeros(0, n)` (zero rows, one column
per flux) and `b = Float64[]`.

**Optimized allocation** keeps the original bounds and lets the amino acids share
the budget:

$$
\sum_{j\in\mathcal{I}}v_j\leq U.
$$

Return this as a one-row `A` with `1` in each supply column and `0` elsewhere,
and `b = [U]`.

Implement in [src/Compute.jl](src/Compute.jl):

1. **`protein_objective(model)`** returns $\mathbf{c}$, with `1` at the
   protein-output column and `0` elsewhere, so that $\mathbf{c}^{\top}\mathbf{v}=v_p$.
2. **`allocation_constraints(model, budget; strategy=:equal)`** returns the named
   tuple `(lower, upper, A, b)` for `:equal` or `:optimized`.

Each docstring lists the required input checks and errors. Test with:

```bash
julia --startup-file=no --project=. testme_part_1.jl
```

## Part 2 (Standard): Production curves

Implement **`production_curve(model, budgets)`**. For each budget, solve once with
each strategy, using your Part 1 functions for the objective and constraints.
Return one entry per budget, in the given order:

```julia
(budget=U, equal=equal_rate, optimized=optimized_rate) # rates in µM/h
```

Follow the [build](src/Factory.jl) and [solve](src/Solver.jl) pattern from the
L6a urea-cycle example. This code solves one strategy at one budget:

```julia
# Needed only outside src/Compute.jl; Include.jl already imports JuMP there.
import JuMP;
strategy = :equal; # repeat with :optimized
limits = allocation_constraints(model, budget; strategy=strategy);
c = protein_objective(model);
calculation = build(MyPrimalFluxBalanceAnalysisCalculationModel, (
    S=model.S,
    fluxbounds=hcat(limits.lower, limits.upper), # lower and upper bound columns
    objective=c,
    species=model.metabolites,
    reactions=model.reactions
));
result = solve(calculation; A=limits.A, b=limits.b);
status = result["termination_status"];
status == JuMP.MOI.OPTIMAL || error("budget=$budget, strategy=$strategy, status=$status");
flux = result["argmax"];
rate = 1000.0 * flux[model.protein]; # mM/h to μM/h
```

Pass `A=limits.A, b=limits.b` for both strategies; under optimized allocation,
they enforce the shared budget. Always check the status before reading
`"argmax"` or `"objective_value"`, because a nonoptimal solve returns `nothing`.
The docstring lists the remaining requirements. Then run:

```bash
julia --startup-file=no --project=. testme_part_2.jl
julia --startup-file=no --project=. runproduction.jl
```

`runproduction.jl` evaluates budgets from 0 to 4 mM/h in steps of 0.1 and writes
to `outputs/`:

- `production-versus-budget.png`, `.svg`, and `.csv`: both production curves. The
  plot also marks the translation capacity, which the terminal prints in µM/h.
- `amino-acid-allocation.csv`: the supply fluxes of one optimal solution per
  strategy at a budget of 1 mM/h.

Use these files to answer questions 1 and 2. Your allocation CSV may differ from a
classmate's: optimal flux vectors are often not unique, and the tests accept any
optimal one.

**Standard track:** Your code is complete; go to
[Check and submit](#check-and-submit). **Advanced track:** Continue to Part 3,
which measures how much optimal supply fluxes can differ.

## Part 3 (Advanced): Which supply fluxes does the optimum determine?

**Advanced track only.** The maximum production rate is unique, but the flux
vector that achieves it often is not. **Flux variability analysis (FVA)** finds
the smallest and largest value of each flux over all flux vectors that reach at
least a fraction $\gamma$ of the optimum. FVA was introduced by **Mahadevan and Schilling (2003), “The effects of
alternate optimal solutions in constraint-based genome-scale metabolic models,”
Metabolic Engineering, 5(4), 264–276**
([paper](https://pubmed.ncbi.nlm.nih.gov/14642354/),
[DOI: 10.1016/j.ymben.2003.09.002](https://doi.org/10.1016/j.ymben.2003.09.002)).

Implement **`supply_ranges(model, budget; strategy=:equal, fraction=1.0)`**, where
`fraction` is $\gamma$. Return one `(column, minimum, maximum)` named tuple per
supply column, in the order of `model.amino_acid_uptake`, with endpoints in mM/h.

Because `solve` only maximizes and accepts extra constraints only as
$\mathbf{A}\mathbf{v}\leq\mathbf{b}$, use three steps:

1. Build the strategy's limits and $\mathbf{c}$ with your Part 1 functions, solve
   once as in Part 2, and set `z = result["objective_value"]`. This is the
   maximum protein flux $z^\star$, in mM/h.
2. Write the production requirement $\mathbf{c}^{\top}\mathbf{v}\geq\gamma z^\star$
   as $-\mathbf{c}^{\top}\mathbf{v}\leq-\gamma z^\star$ and append it as one row:
   `A = vcat(limits.A, -c')` and `b = vcat(limits.b, -fraction * z)`.
3. For each supply column $j$, solve twice with the step 1 bounds and the new
   `A` and `b`, setting `calculation.objective` before each `solve`:
   - **Maximum:** maximize an objective that is `1` at column $j$ and `0` elsewhere.
   - **Minimum:** maximize one that is `-1` at column $j$, then negate
     `result["objective_value"]`.

Set [TRACK.txt](TRACK.txt) to `advanced`, then run:

```bash
julia --startup-file=no --project=. testme_part_3.jl
julia --startup-file=no --project=. runproduction.jl
```

The driver adds supply-range plots and tables to `outputs/`. Use
`supply-ranges.csv` to answer question 3. Unlike optimal flux vectors, range
endpoints are unique, so the tests compare them with reference values.

## Check and submit

Run the checker. It reads [TRACK.txt](TRACK.txt) and runs the checks for your
track:

```bash
julia --startup-file=no --project=. check_submission.jl
```

The Standard track runs the 48 Standard tests and checks that questions 1–2 are
answered. The Advanced track adds the 24 Part 3 tests and question 3. The
checker always writes `MANIFEST.txt`. It clears old reports from `outputs/` and
regenerates them only when the tests pass. It runs only on your computer; it
does not submit anything.

Zip the whole folder, including `TRACK.txt`, `MANIFEST.txt`, and `outputs/`, as
`CHEME-4800-5800-PS3-<netid>.zip` (for example, `CHEME-4800-5800-PS3-abc123.zip`)
and upload it to Canvas. If some tests still fail at the deadline, submit anyway
for partial credit.

## Policies

- **Revisions:** Submit something by the deadline, even if it does not run. You
  may then revise and resubmit until the end of the semester; your highest score
  counts. No initial submission means a score of `0` and no revisions.
- **Reference solution:** Released after the deadline. Use it to understand and
  debug your work; copying it results in a score of `0`.
- **Collaboration and AI:** Discuss ideas with classmates, but do not share code
  or solutions. You may use Julia documentation, AI tools, and internet
  resources. You must be able to explain your code.

## Source

The model is adapted from **Vilkhovoy et al. (2018), “Sequence Specific Modeling
of _E. coli_ Cell-Free Protein Synthesis,” ACS Synthetic Biology, 7(8),
1844–1857** ([paper](https://pubmed.ncbi.nlm.nih.gov/29944340/),
[DOI: 10.1021/acssynbio.7b00465](https://doi.org/10.1021/acssynbio.7b00465)).
The [data notes](data/README.md) link the source code and describe the PS3
operating conditions.
