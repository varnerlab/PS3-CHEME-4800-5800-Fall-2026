# Load the assignment from its own directory, independent of the working directory.
module CellFreeProduction

using DelimitedFiles # supplied stoichiometric matrix and reaction tables
using LinearAlgebra  # matrix and vector operations
import JuMP          # supplied linear-programming interface
import GLPK          # supplied LP solver

include(joinpath(@__DIR__, "src", "Model.jl"));
include(joinpath(@__DIR__, "src", "Solver.jl"));
# Load any student helper files from src/ here, before Compute.jl.
include(joinpath(@__DIR__, "src", "Compute.jl"));

export load_model, reaction_index, solve_lp, check_solution;
export protein_objective, allocation_constraints, production_curve;

end

using .CellFreeProduction;
