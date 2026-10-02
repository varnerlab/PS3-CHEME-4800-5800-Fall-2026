"""
    solve_lp(S, lower, upper, c; A=zeros(0, length(c)), b=Float64[]) -> NamedTuple

Maximize `dot(c, v)` subject to `S*v == 0`, the flux bounds, and `A*v <= b`.

# Arguments
- `S`: an m-by-n stoichiometric matrix, with species in rows and reactions in columns.
- `lower`, `upper`: length-n vectors defining `lower[j] <= v[j] <= upper[j]`.
- `c`: a length-n objective-coefficient vector in the same reaction order.
- `A`, `b`: a k-by-n matrix and length-k vector specifying extra inequalities.
  The defaults add no inequalities. All matrix and vector entries must be finite.

# Returns
A named tuple with fields `status` (solver termination status), `optimal`
(whether the status is `OPTIMAL`), `flux` (reaction-flux vector), and `objective`
(the scalar value of `dot(c, flux)`). If the status is not `OPTIMAL`, both `flux`
and `objective` are `nothing`. Leave all inputs unchanged.

# Errors
Throw `DimensionMismatch` if dimensions disagree, or `ArgumentError` for
nonfinite entries or a lower bound greater than its upper bound.
"""
function solve_lp(S, lower, upper, c; A=zeros(0, length(c)), b=Float64[])
    # Check that every coefficient and bound refers to a valid reaction column -
    n = size(S, 2);
    length(lower) == length(upper) == length(c) == n ||
        throw(DimensionMismatch("one bound and objective coefficient is required per reaction"));
    size(A) == (length(b), n) || throw(DimensionMismatch("A and b dimensions do not agree"));
    all(all(isfinite, x) for x in (S, lower, upper, c, A, b)) ||
        throw(ArgumentError("LP inputs must be finite"));
    all(lower .<= upper) || throw(ArgumentError("lower bounds must not exceed upper bounds"));
    # Build the continuous linear program -
    lp = JuMP.Model(GLPK.Optimizer);
    JuMP.set_silent(lp);
    JuMP.@variable(lp, lower[j] <= v[j=1:n] <= upper[j]);
    JuMP.@constraint(lp, S * v .== 0);
    JuMP.@constraint(lp, A * v .<= b);
    JuMP.@objective(lp, Max, sum(c[j] * v[j] for j in 1:n));
    JuMP.optimize!(lp);
    # Return fluxes only when the solver reports an optimal solution -
    status = JuMP.termination_status(lp);
    status == JuMP.MOI.OPTIMAL ||
        return (status=status, optimal=false, flux=nothing, objective=nothing);
    return (status=status, optimal=true, flux=JuMP.value.(v), objective=JuMP.objective_value(lp));
end

"""
    check_solution(model, result, constraints; atol=1e-7) -> NamedTuple

Check the returned fluxes against balances, bounds, and allocation inequalities.

# Arguments
- `model`: the model whose `S` defines the species balances.
- `result`: an optimal result from `solve_lp`, with a reaction-ordered `flux` vector.
- `constraints`: the named tuple `lower`, `upper`, `A`, `b` used in that solve.
- `atol::Real`: finite, nonnegative absolute tolerance for each balance, bound,
  and inequality check (mM/h for the PS3 constraints). Dimensions must agree.

# Returns
A named tuple containing `balance_ok`, `bounds_ok`, `allocation_ok`, their
conjunction `valid`, and `maximum_residual` (the largest absolute entry of
`model.S * result.flux`). These checks establish feasibility; optimality is
reported by the solver. Leave all inputs unchanged.

# Errors
Throw `ArgumentError` for a nonoptimal result or an invalid tolerance.
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
