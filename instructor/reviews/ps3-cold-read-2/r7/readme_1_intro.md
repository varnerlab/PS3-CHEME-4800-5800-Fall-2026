AUDIT THIS PART: README.md:B001 to README.md:B018

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
