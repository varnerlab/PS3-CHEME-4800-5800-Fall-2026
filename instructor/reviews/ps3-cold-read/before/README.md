# Problem Set 3 (PS3): Feeding a cell-free protein factory

A cell-free protein synthesis system uses cellular machinery in a reaction
mixture to make a protein. That machinery needs amino acids, nucleotides, and
energy. The metabolic reactions in the mixture can make amino acids from other
available compounds, and amino acids can also be supplied directly.

You will use linear programming to predict the production rate of deGFP, a green
fluorescent protein, under two amino-acid supply strategies. With **equal
allocation**, every amino acid receives the same supply limit. With **optimized
allocation**, the solver distributes a shared supply budget among the amino acids.
Your final result is a plot of protein production against the total supply budget.

This assignment extends the flux balance analysis in L6b. The supplied model is
adapted from the model described by **Vilkhovoy et al. (2018), “Sequence Specific
Modeling of E. coli Cell-Free Protein Synthesis,” ACS Synthetic Biology, 7(8),
1844–1857** ([paper](https://pubmed.ncbi.nlm.nih.gov/29944340/),
[DOI: 10.1021/acssynbio.7b00465](https://doi.org/10.1021/acssynbio.7b00465)).
The code and model data are available in the Varner laboratory's
[sequence-specific cell-free protein synthesis repository](https://github.com/varnerlab/Sequence-Specific-FBA-CFPS-Publication-Code).
The [data notes](data/README.md) identify the source and the PS3 operating conditions.

## Logistics

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

## Getting started

Use Julia `1.12.7`. The assignment contains its own data and source files.
Open the extracted assignment folder in VS Code. Run these commands from the
folder containing this README and [Project.toml](Project.toml):

```bash
julia --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate()'
julia --startup-file=no --project=. check_submission.jl
```

The first command installs the required packages and can take several minutes.
The second runs the assignment tests. The three functions in
[src/Compute.jl](src/Compute.jl) initially raise errors because you must implement
them. Tests will fail until that work is complete.

Complete those functions and the three prompts in [responses.md](responses.md).
Keep any helper files in `src/` and load them from [Include.jl](Include.jl).
The model loader, LP solver, and plotting code are supplied. Keep the data,
supplied code, and test files unchanged.

## The model and its operating conditions

The model contains 146 species and 265 reactions. It includes metabolism,
transcription, translation, amino-acid supply, and protein output. The supplied
[load_model function](src/Model.jl) reads the stoichiometric matrix and sets the
operating conditions. Its returned named tuple contains:

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

The glucose supply limit is **1 mM/h** and the oxygen supply limit is
**100 mM/h**. Amino-acid synthesis is enabled; the model's designated
amino-acid degradation reactions and amino-acid removal reactions are blocked.
Transcription and translation bounds are calculated from fixed expression
parameters before solving. Use the supplied conditions for both strategies.

Supply reactions are written as `[] --> species` and have **nonnegative**
fluxes. Each amino-acid supply flux therefore contributes positively to the
budget. The deGFP output reaction removes protein from the balanced network and
records its production rate. Metabolic intermediates satisfy steady-state
balances while protein leaves through that output reaction.

All model fluxes are in **mM/h**, millimoles per liter per hour. Report protein
production in **µM/h** by multiplying its model flux by 1000. The supply budget
is a rate limit; it does not specify an initial concentration or a monetary cost.
This static model predicts feasible production rates under fixed conditions.

## Part 1: Construct the objective and supply constraints

Let the vector $\mathbf{v}$ contain the reaction fluxes, and let $p$ be the
protein-output column. The production problem is given by:

$$
\begin{aligned}
\underset{\mathbf{v}}{\operatorname{maximize}}\quad & v_p\\
\text{subject to}\quad & \mathbf{S}\mathbf{v}=\mathbf{0},\\
& \boldsymbol{\ell}\leq\mathbf{v}\leq\mathbf{u},\\
& \mathbf{A}\mathbf{v}\leq\mathbf{b}.
\end{aligned}
$$

The vectors $\boldsymbol{\ell}$ and $\mathbf{u}$ specify the flux bounds.
The matrix $\mathbf{A}$ and vector $\mathbf{b}$ hold any additional supply
inequality. Define $\mathcal{I}$ as the set of amino-acid supply columns, let
$K=|\mathcal{I}|$, and let $U\geq0$ be the total supply budget in mM/h.

**Equal allocation.** Each amino acid may use at most one equal share of the
budget. Its upper bound is given by:

$$
u_j^{\mathrm{equal}}=\min(u_j,U/K),\qquad j\in\mathcal{I}.
$$

The supply fluxes can be smaller than their limits. Any unused share remains
unavailable to the other amino acids. These individual limits already enforce
the total budget, so return an `A` with zero rows and an empty `b`.

**Optimized allocation.** Retain the original individual bounds and add one
shared supply constraint, given by:

$$
\sum_{j\in\mathcal{I}}v_j\leq U.
$$

The solver may distribute the available supply unequally. One row of `A` has
a coefficient of one at each amino-acid supply column and zero elsewhere;
the corresponding entry of `b` is `U`. Metabolic synthesis of amino acids remains
available under both strategies and does not count as external supply.

Implement two functions in [src/Compute.jl](src/Compute.jl):

1. **protein_objective(model)** returns a `Vector{Float64}` with one coefficient
   per reaction. Its protein-output coefficient is one; the others are zero.
2. **allocation_constraints(model, budget; strategy=:equal)** returns
   `(lower, upper, A, b)` using the rules above. Support `:equal` and `:optimized`.
   Copy the bound arrays so that solving one scenario cannot change the next.
   Reject a negative or nonfinite budget or an unknown strategy with `ArgumentError`.

Use the model's supplied column indices. Your functions must also work on the
small test networks, which have different numbers and orders of reactions and
amino-acid supplies. All test models have nonnegative supply fluxes and at least
one supply column. Keep all lower bounds and unrelated upper bounds unchanged.

Run the Part 1 checks with:

```bash
julia --startup-file=no --project=. testme_part_1.jl
```

For a hand check, consider two amino acids with supply fluxes $v_1$ and $v_2$.
One unit of protein output $v_p$ consumes two units of the first amino acid and
one unit of the second. The balances are $v_1=2v_p$ and $v_2=v_p$.
With a total supply budget of 6 and all other bounds loose, equal allocation
allows at most **1.5** units of protein output; optimized allocation allows
**2**. Explain these numbers to yourself before testing the larger network.

## Part 2: Compute and interpret the production curves

Implement **production_curve(model, budgets)** in
[src/Compute.jl](src/Compute.jl). For each budget, construct both sets of
constraints and solve both LPs. Use your Part 1 functions and the supplied
[solve_lp function](src/Solver.jl). A single solve has this form:

```julia
limits = allocation_constraints(model, budget; strategy=:equal);
c = protein_objective(model);
result = solve_lp(model.S, limits.lower, limits.upper, c;
    A=limits.A, b=limits.b);
```

The result contains `status`, `optimal`, `flux`, and `objective`. Check
`result.optimal` before reading the fluxes. If a solve is not optimal, throw an
`ErrorException` naming the budget, strategy, and status. In that case the solver
returns `nothing` for `flux` and `objective`.

Return a vector containing one named tuple `(budget, equal, optimized)` per
input budget. The budget remains in mM/h; both production rates must be in µM/h.
Preserve the input order, including repeated budgets. An empty input returns an
empty vector. Negative or nonfinite budgets raise `ArgumentError`.

After completing the function, run:

```bash
julia --startup-file=no --project=. testme_part_2.jl
julia --startup-file=no --project=. runproduction.jl
```

The supplied driver evaluates budgets from **0 to 4 mM/h in steps of 0.1** and
writes these files:

- `outputs/production-versus-budget.png` and `.svg`: the two production curves
  and the fixed translation-capacity limit.
- `outputs/production-versus-budget.csv`: the numerical values used in the plot.
- `outputs/amino-acid-allocation.csv`: actual supply fluxes for both strategies
  at a total budget of **1 mM/h**, with independent feasibility checks reported
  in the terminal.

Open the plot and use it with the tables to answer [responses.md](responses.md).
The optimized allocation can have multiple equally productive flux distributions.
Tests check feasibility and production, so they accept alternative optimal allocations.

## Checking and submitting your work

Run the complete check after finishing the code and responses:

```bash
julia --startup-file=no --project=. check_submission.jl
```

The checker runs **24 tests per part, 48 total**, records their results and source
file fingerprints in `MANIFEST.txt`, and warns about unanswered discussion prompts.
If both test suites pass, it generates the plot and data tables from your code.
The instructor still reviews the implementation and written explanations.

Zip the entire assignment folder, including your source, data, responses,
`MANIFEST.txt`, and generated files in `outputs/`. Name the archive
`CHEME-4800-5800-PS3-<your netid>.zip` and upload it to Canvas. The checker runs
locally; you must upload the ZIP yourself.

If some tests still fail at the deadline, submit your current work to receive
partial credit and participate in the infinite-revision policy.
