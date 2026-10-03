# ----------------------------------------------------------------------------
# Supplied builder | Populate the calculation type used in the urea example.
# ----------------------------------------------------------------------------

"""
    build(::Type{MyPrimalFluxBalanceAnalysisCalculationModel}, data::NamedTuple)

Populate a calculation model as in the L6a urea example.

# Arguments

- `modeltype`: `MyPrimalFluxBalanceAnalysisCalculationModel`.
- `data`: named tuple with `S::Matrix{Float64}` (m-by-n),
  `fluxbounds::Matrix{Float64}` (n-by-2), `objective::Vector{Float64}` (length n),
  `species::Vector{String}` (length m), and `reactions::Vector{String}` (length n).
  Here m counts species and n counts reactions. Species follow matrix rows;
  reactions follow matrix columns, bound rows, and objective entries.

# Returns and mutation

Return a `MyPrimalFluxBalanceAnalysisCalculationModel`. Arrays of the declared
field types are attached without copying. Editing their entries through the
calculation also changes the supplied arrays. Replacing a field, such as
`calculation.objective = c`, attaches that array instead.

For PS3, form `fluxbounds=hcat(limits.lower, limits.upper)` from the copied
Part 1 bounds, use a fresh objective vector, and set `species=model.metabolites`.
The network matrix and labels can be shared because the solver only reads them.
Pass the extra inequalities separately to `solve(calculation; A=limits.A, b=limits.b)`.

# Errors

Missing fields or values that cannot be converted to the field types raise
Julia errors. Construction does not check dimensions, finite coefficients, or
bound order; `solve` checks those before optimization.
"""
function build(
    modeltype::Type{MyPrimalFluxBalanceAnalysisCalculationModel},
    data::NamedTuple,
)

    # Initialize an empty calculation model -
    model = modeltype();

    # Attach the matrix, bounds, and objective in reaction order -
    model.S = data.S;
    model.fluxbounds = data.fluxbounds;
    model.objective = data.objective;

    # Attach labels in the same order as the matrix rows and columns -
    model.species = data.species;
    model.reactions = data.reactions;

    # Return the populated calculation -
    return model;
end
