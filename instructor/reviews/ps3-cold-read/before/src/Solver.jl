"""
    solve_lp(S, lower, upper, c; A=zeros(0, length(c)), b=Float64[]) -> NamedTuple

Maximize `dot(c, v)` subject to `S*v == 0`, `lower <= v <= upper`, and `A*v <= b`.
All inputs are unchanged. `A` has one row per additional inequality. Bounds and
coefficients must be finite, dimensions must agree, and lower bounds cannot
exceed upper bounds. Invalid values raise `ArgumentError`; invalid dimensions
raise `DimensionMismatch`.

Return `(status, optimal, flux, objective)`. If the status is not `OPTIMAL`,
`optimal` is false and both `flux` and `objective` are `nothing`.
"""
function solve_lp(S, lower, upper, c; A=zeros(0, length(c)), b=Float64[])
    n = size(S, 2);
    length(lower) == length(upper) == length(c) == n ||
        throw(DimensionMismatch("one bound and objective coefficient is required per reaction"));
    size(A) == (length(b), n) || throw(DimensionMismatch("A and b dimensions do not agree"));
    all(all(isfinite, x) for x in (S, lower, upper, c, A, b)) ||
        throw(ArgumentError("LP inputs must be finite"));
    all(lower .<= upper) || throw(ArgumentError("lower bounds must not exceed upper bounds"));
    lp = JuMP.Model(GLPK.Optimizer);
    JuMP.set_silent(lp);
    JuMP.@variable(lp, lower[j] <= v[j=1:n] <= upper[j]);
    JuMP.@constraint(lp, S * v .== 0);
    JuMP.@constraint(lp, A * v .<= b);
    JuMP.@objective(lp, Max, sum(c[j] * v[j] for j in 1:n));
    JuMP.optimize!(lp);
    status = JuMP.termination_status(lp);
    status == JuMP.MOI.OPTIMAL ||
        return (status=status, optimal=false, flux=nothing, objective=nothing);
    return (status=status, optimal=true, flux=JuMP.value.(v), objective=JuMP.objective_value(lp));
end

"""
    check_solution(model, result, constraints; atol=1e-7) -> NamedTuple

Check balances, bounds, and additional inequalities independently of the solver.
`constraints` contains `lower`, `upper`, `A`, and `b`. Return `valid`, `balance_ok`,
`bounds_ok`, `allocation_ok`, and `maximum_residual`. This checks feasibility;
the solver status reports optimality. Throw `ArgumentError` for a nonoptimal
result or an invalid tolerance. Inputs are unchanged.
"""
function check_solution(model, result, constraints; atol::Real=1e-7)
    result.optimal || throw(ArgumentError("check the solver status before reading fluxes"));
    isfinite(atol) && atol >= 0 || throw(ArgumentError("atol must be finite and nonnegative"));
    v = result.flux;
    maximum_residual = maximum(abs, model.S * v; init=0.0);
    balance_ok = all(isfinite, v) && maximum_residual <= atol;
    bounds_ok = all(constraints.lower .- atol .<= v) && all(v .<= constraints.upper .+ atol);
    allocation_ok = all(constraints.A * v .<= constraints.b .+ atol);
    return (valid=balance_ok && bounds_ok && allocation_ok, balance_ok=balance_ok,
        bounds_ok=bounds_ok, allocation_ok=allocation_ok, maximum_residual=maximum_residual);
end
