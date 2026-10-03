# ----------------------------------------------------------------------------
# Supplied type | Store the arrays for one flux balance calculation.
# ----------------------------------------------------------------------------

"""Abstract base type for flux calculation models, as in the L6a urea example."""
abstract type AbstractFluxCalculationModel end

"""
    MyPrimalFluxBalanceAnalysisCalculationModel <: AbstractFluxCalculationModel

Store a flux balance calculation using the same fields as the L6a urea example.
Use `build` to populate the zero-argument constructor before calling `solve`.
The empty constructor leaves every field undefined; it is not a solvable model.
The solver maximizes `objective' * flux` subject to steady-state balances and
flux bounds, with optional extra inequalities supplied to `solve`.

# Fields

Let m be the number of balanced species and n the number of reaction fluxes.
- `S::Matrix{Float64}`: m-by-n stoichiometric matrix, species in rows.
- `fluxbounds::Matrix{Float64}`: n-by-2 matrix, lower bounds in column 1 and
  upper bounds in column 2.
- `objective::Vector{Float64}`: n objective coefficients in reaction order.
- `species::Vector{String}`: m identifiers in matrix row order.
- `reactions::Vector{String}`: n identifiers in matrix column order.

Fluxes and bounds use mM/h in PS3. Positive flux follows the written reaction.
The positive protein-output flux is selected with objective coefficient +1.
The type stores calculation arrays only; the named tuple from `load_model`
also holds the amino-acid indices and expression parameters used to set them.

# Mutation

Fields can be replaced to change a calculation, as in the course example.
`build` shares the supplied arrays, so use copied bounds and a new objective
when editing a scenario. `solve` reads the arrays without changing them.
"""
mutable struct MyPrimalFluxBalanceAnalysisCalculationModel <: AbstractFluxCalculationModel

    # Numerical data -
    S::Array{Float64,2};
    fluxbounds::Array{Float64,2}; # one reaction per row; lower then upper bound
    objective::Array{Float64,1};

    # Matrix row and column labels -
    species::Array{String,1};
    reactions::Array{String,1};

    # Empty constructor; build fills all five fields -
    MyPrimalFluxBalanceAnalysisCalculationModel() = new();
end
