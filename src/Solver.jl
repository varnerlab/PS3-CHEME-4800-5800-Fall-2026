# ----------------------------------------------------------------------------
# Supplied solver | Solve one calculation and check its constraints.
# ----------------------------------------------------------------------------

"""
    solve(model::MyPrimalFluxBalanceAnalysisCalculationModel;
        A=zeros(0, length(model.objective)), b=Float64[]) -> Dict{String,Any}

Maximize `dot(model.objective, v)` subject to `model.S*v == 0`, the flux
bounds, and `A*v <= b`.

# Arguments

- `model`: calculation model built as in the L6a urea example. Its `S` is m-by-n,
  `fluxbounds` is n-by-2 (lower, upper), and `objective` is length n. The
  `species` and `reactions` identifiers follow the matrix row and column order.
- `A`, `b`: a k-by-n matrix and length-k vector specifying extra inequalities.
  Each row of `A` has one coefficient per reaction. The defaults have zero rows
  and add no inequalities; use these keywords for the shared supply budget or
  the protein-output requirement. All numerical entries must be finite.
  Bounds and `b` use mM/h for PS3's dimensionless constraint coefficients.

# Returns

A dictionary with the course example's keys: `"termination_status"`, `"argmax"`
(reaction-flux vector in mM/h), and `"objective_value"` (the scalar objective).
For the protein objective, the objective value is the production rate in mM/h.
`"argmax"` holds a vector attaining the maximum, not a reaction index. Different
optimal vectors can have the same objective value.

Check `result["termination_status"] == JuMP.MOI.OPTIMAL` before reading values.
Unlike the urea example's assertion on failure, a nonoptimal solve returns
`nothing` for both values so callers can report the budget and strategy.
Leave the model and all input arrays unchanged.

# Errors

Throw `DimensionMismatch` if dimensions disagree, or `ArgumentError` for
nonfinite entries or a lower bound greater than its upper bound.
Other model-construction and solver errors propagate to the caller. A nonoptimal
termination status is returned in the dictionary rather than thrown as an error.
"""
function solve(model::MyPrimalFluxBalanceAnalysisCalculationModel;
    A=zeros(0, length(model.objective)),
    b=Float64[],
)

    S, c = model.S, model.objective; # shared inputs; this function only reads them

    # 1. Check dimensions and bounds -
    n = size(S, 2);

    size(model.fluxbounds) == (n, 2) ||
        throw(DimensionMismatch("fluxbounds must have one row per reaction and two columns"));

    length(model.species) == size(S, 1) && length(model.reactions) == n ||
        throw(DimensionMismatch("species and reaction identifiers must match S"));

    lower, upper = model.fluxbounds[:, 1], model.fluxbounds[:, 2];

    length(lower) == length(upper) == length(c) == n ||
        throw(DimensionMismatch("one bound and objective coefficient is required per reaction"));

    size(A) == (length(b), n) ||
        throw(DimensionMismatch("A and b dimensions do not agree"));

    all(all(isfinite, x) for x in (S, lower, upper, c, A, b)) ||
        throw(ArgumentError("LP inputs must be finite"));

    all(lower .<= upper) ||
        throw(ArgumentError("lower bounds must not exceed upper bounds"));

    # 2. Build the continuous linear program -
    lp = JuMP.Model(GLPK.Optimizer);
    JuMP.set_silent(lp);
    JuMP.@variable(lp, lower[j] <= v[j=1:n] <= upper[j]);
    JuMP.@constraint(lp, S * v .== 0); # no net accumulation of any balanced species
    JuMP.@constraint(lp, A * v .<= b); # zero rows add no constraints
    JuMP.@objective(lp, Max, sum(c[j] * v[j] for j in 1:n));

    # 3. Solve and inspect the termination status -
    JuMP.optimize!(lp);

    # Return fluxes only when the solver reports an optimal solution -
    status = JuMP.termination_status(lp);
    status == JuMP.MOI.OPTIMAL ||
        return Dict{String,Any}(
            "termination_status" => status,
            "argmax" => nothing,
            "objective_value" => nothing,
        );

    # 4. Return the course example's result fields -
    return Dict{String,Any}(
        "termination_status" => status,
        "argmax" => JuMP.value.(v),
        "objective_value" => JuMP.objective_value(lp),
    );
end


"""
    check_solution(model, result, constraints; atol=1e-7) -> NamedTuple

Check the returned fluxes against balances, bounds, and allocation inequalities.

# Arguments

- `model`: loaded data or a calculation model whose `S` defines the balances.
- `result`: an optimal result from `solve`, with reaction-ordered `result["argmax"]`.
- `constraints`: the named tuple `lower`, `upper`, `A`, `b` used in that solve.
  Pass the scenario's bounds, which may be tighter than the original bounds.
- `atol::Real`: finite, nonnegative absolute tolerance for each balance, bound,
  and inequality check (mM/h for the PS3 constraints). Dimensions must agree.

# Returns

A named tuple containing `balance_ok`, `bounds_ok`, `allocation_ok`, their
conjunction `valid`, and `maximum_residual` (the largest absolute entry of
`model.S * result["argmax"]`). These checks establish feasibility; optimality is
reported by the solver. Leave all inputs unchanged.
The residual measures balance error only, not the largest bound or budget error.

# Errors

Throw `ArgumentError` for a nonoptimal result or an invalid tolerance.
"""
function check_solution(model, result, constraints; atol::Real=1e-7)

    # Validate the result and verification tolerance -
    result["termination_status"] == JuMP.MOI.OPTIMAL ||
        throw(ArgumentError("check the solver status before reading fluxes"));

    isfinite(atol) && atol >= 0 ||
        throw(ArgumentError("atol must be finite and nonnegative"));

    v = result["argmax"];

    # Check the solved scenario; the tolerance is for verification, not optimization -
    maximum_residual = maximum(abs, model.S * v; init=0.0); # mM/h imbalance
    balance_ok = all(isfinite, v) && maximum_residual <= atol;

    bounds_ok = all(constraints.lower .- atol .<= v) &&
        all(v .<= constraints.upper .+ atol);

    allocation_ok = all(constraints.A * v .<= constraints.b .+ atol);

    # Return each check separately so a failed constraint can be identified -
    return (
        valid=balance_ok && bounds_ok && allocation_ok,
        balance_ok=balance_ok,
        bounds_ok=bounds_ok,
        allocation_ok=allocation_ok,
        maximum_residual=maximum_residual,
    );
end
