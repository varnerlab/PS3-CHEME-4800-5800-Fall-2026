# Load the assignment from its own directory, independent of the working directory.
"""
    CellFreeProduction

Load the supplied model and solver utilities and the four student functions.
`load_model()` returns the baseline data. Use `build` to put a scenario's arrays
in a `MyPrimalFluxBalanceAnalysisCalculationModel`, then pass it to `solve`.
The student functions and the reference implementation use the same utilities.
Pass `--solution` after a script's filename to load `solution/src/Compute.jl`
instead of the student file. The local solution must exist; no files are copied
or replaced. Without that flag, always load `src/Compute.jl`.
Include this file once per Julia session to make the exported names available.
Restart Julia after editing source files, or run the scripts in a new process.
JuMP and GLPK are imported inside this module; a separate script that names
`JuMP.MOI.OPTIMAL` also needs its own `import JuMP`.
"""
module CellFreeProduction

using DelimitedFiles # supplied stoichiometric matrix and reaction tables
using LinearAlgebra  # matrix and vector operations
import JuMP          # supplied linear-programming interface
import GLPK          # supplied LP solver

# Load the calculation type before the builder and solver that dispatch on it -
include(joinpath(@__DIR__, "src", "Types.jl"));
include(joinpath(@__DIR__, "src", "Factory.jl"));
include(joinpath(@__DIR__, "src", "Model.jl"));
include(joinpath(@__DIR__, "src", "Solver.jl"));
# If you add helper files, include them here in dependency order before Compute.jl.
# For example, after creating src/MyHelpers.jl:
# include(joinpath(@__DIR__, "src", "MyHelpers.jl"));
# Select exactly one implementation; never fall back when the solution is absent -
const COMPUTE_PATH = "--solution" in ARGS ?
    joinpath(@__DIR__, "solution", "src", "Compute.jl") :
    joinpath(@__DIR__, "src", "Compute.jl");
isfile(COMPUTE_PATH) || error("Implementation file is missing: $COMPUTE_PATH. --solution requires the local instructor solution.");
include(COMPUTE_PATH);

export AbstractFluxCalculationModel, MyPrimalFluxBalanceAnalysisCalculationModel, build;
export load_model, reaction_index, solve, check_solution;
export protein_objective, allocation_constraints, production_curve, supply_ranges;

end

using .CellFreeProduction; # make the assignment functions available to the calling script
