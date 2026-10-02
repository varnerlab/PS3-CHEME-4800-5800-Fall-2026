# Load the assignment from its own directory, independent of the working directory.
"""
    CellFreeProduction

Load the supplied model and solver utilities and the four student functions.
Include this file once per Julia session to make the exported names available.
Restart Julia after editing source files, or run the scripts in a new process.
"""
module CellFreeProduction

using DelimitedFiles # supplied stoichiometric matrix and reaction tables
using LinearAlgebra  # matrix and vector operations
import JuMP          # supplied linear-programming interface
import GLPK          # supplied LP solver

include(joinpath(@__DIR__, "src", "Model.jl"));
include(joinpath(@__DIR__, "src", "Solver.jl"));
# If you add helper files, include them here in dependency order before Compute.jl.
# For example, after creating src/MyHelpers.jl:
# include(joinpath(@__DIR__, "src", "MyHelpers.jl"));
include(joinpath(@__DIR__, "src", "Compute.jl"));

export load_model, reaction_index, solve_lp, check_solution;
export protein_objective, allocation_constraints, production_curve, supply_ranges;

end

using .CellFreeProduction; # make the assignment functions available to the calling script
