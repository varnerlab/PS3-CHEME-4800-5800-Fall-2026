# Complete these three functions. Keep their arguments and return fields.

"""
    protein_objective(model) -> Vector{Float64}

Build an objective vector with one coefficient per column of `model.S`.
Set the coefficient at `model.protein` to one and every other coefficient to zero.
Maximizing this objective maximizes protein output in mM/h. Leave the model unchanged.
"""
function protein_objective(model)
    # TODO: Construct the protein-production objective.
    error("Complete protein_objective in src/Compute.jl");
end

"""
    allocation_constraints(model, budget::Real; strategy::Symbol=:equal) -> NamedTuple

Construct `(lower, upper, A, b)` for a finite, nonnegative amino-acid supply
`budget` in mM/h. Copy the model's bounds; leave the model unchanged.
`model.amino_acid_uptake` lists the nonnegative supply-flux columns.
Let K be the length of that list and n the number of columns in `model.S`.

For `:equal`, set each supply upper bound to `min(original_upper, budget/K)`.
Return `A=zeros(0, n)` and `b=Float64[]`; the individual limits enforce the budget.
For `:optimized`, retain the model bounds and add one inequality limiting the
sum of the K supply fluxes to `budget`. Return A of size (1,n) and b of length 1.
Keep all lower bounds and all other upper bounds unchanged.

Throw `ArgumentError` for a negative or nonfinite budget, or an unknown strategy.
Models supplied to this function have at least one amino-acid supply column.
"""
function allocation_constraints(model, budget::Real; strategy::Symbol=:equal)
    # TODO: Validate the arguments and construct the bounds and inequalities.
    error("Complete allocation_constraints in src/Compute.jl");
end

"""
    production_curve(model, budgets::AbstractVector{<:Real}) -> Vector{NamedTuple}

For each supply budget in mM/h, solve both allocation strategies using
`protein_objective`, `allocation_constraints`, and the supplied `solve_lp`.
Return one named tuple `(budget, equal, optimized)` per input entry. `budget`
is in mM/h; `equal` and `optimized` are protein production rates in micromolar/h
(multiply the protein flux in mM/h by 1000). Preserve input order and duplicates.
Return an empty vector for an empty input. Leave the model and budgets unchanged.

Throw `ArgumentError` for any negative or nonfinite budget. If either LP is
not optimal, throw `ErrorException` with the budget, strategy, and solver status
in its message. Check `result.optimal` before using a flux or objective value.
"""
function production_curve(model, budgets::AbstractVector{<:Real})
    # TODO: Solve the two LPs at each budget and collect the production rates.
    error("Complete production_curve in src/Compute.jl");
end
