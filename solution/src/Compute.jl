# PS3 reference solution | Completed versions of the four functions in src/Compute.jl.
# Parts 1 and 2 (Standard): protein_objective, allocation_constraints, production_curve.
# Part 3 (Advanced): supply_ranges. The docstrings are the same as in the starter.

# ----------------------------------------------------------------------------
# Part 1a: Select protein production
# ----------------------------------------------------------------------------

"""
    protein_objective(model) -> Vector{Float64}

Build the objective for maximizing protein output in mM/h.

# Arguments

- `model`: a named tuple with a stoichiometric matrix `S` and the integer
  column index `protein` of the protein-output reaction. The index must refer
  to a column of `S`; the caller supplies consistent model data.

# Returns

- A new `Vector{Float64}` with one entry per column of `model.S`. The entry at
  `model.protein` is one; every other entry is zero. Leave the model unchanged.

The coefficients are dimensionless, so `dot(c, flux)` is the protein rate in
mM/h. This output reaction carries positive flux; the urea example uses -1
because its uptake-positive exchange carries a negative flux during export.
"""
function protein_objective(model)

    # Start with a zero coefficient for every reaction column -
    n = size(model.S, 2);
    c = zeros(Float64, n);

    # Select the protein-output flux, so c' * v is the protein rate in mM/h -
    c[model.protein] = 1.0;

    return c;
end


# ----------------------------------------------------------------------------
# Part 1b: Set the amino-acid supply limits
# ----------------------------------------------------------------------------

"""
    allocation_constraints(model, budget::Real; strategy::Symbol=:equal) -> NamedTuple

Construct the bounds and supply constraints for one allocation strategy.

# Arguments

- `model`: a named tuple with `S`, bound vectors `lower` and `upper`, and the
  distinct supply-column indices `amino_acid_uptake`. Supply fluxes are
  nonnegative; there is at least one supply column.
- `budget::Real`: a finite, nonnegative total amino-acid supply limit in mM/h.
- `strategy::Symbol`: `:equal` (default) or `:optimized`.

# Returns

A named tuple with fields `lower`, `upper`, `A`, and `b`. Use
`hcat(lower, upper)` for `calculation.fluxbounds` and pass `A`, `b` to `solve`.
Let K be the number of supply columns and n the number of columns in `model.S`.
Copy both bound vectors; leave the model unchanged. The returned bounds are
length n, and their order must match the columns of `S`.

**Equal allocation (`:equal`).** Set each supply upper bound to
`min(original_upper, budget/K)`.
Return `A=zeros(0, n)` and `b=Float64[]`; the individual limits enforce the budget.

**Optimized allocation (`:optimized`).** Keep the model bounds and add one
inequality limiting the sum of the K supply fluxes to `budget`. Return `A` of size (1,n), with ones in
the supply columns and zeros elsewhere, and a length-one vector `b=[budget]`.

Keep all lower bounds and all other upper bounds unchanged. Both strategies
limit supply rates rather than prescribing them: an upper bound may be unused.
This function constructs constraints only; it does not solve or test feasibility.

# Errors

Throw `ArgumentError` for a negative or nonfinite budget (`Inf`, `-Inf`, or
`NaN`), or a strategy other than `:equal` or `:optimized`.
"""
function allocation_constraints(model, budget::Real; strategy::Symbol=:equal)

    # Check the inputs before building anything -
    isfinite(budget) && budget >= 0 ||
        throw(ArgumentError("budget must be finite and nonnegative"));

    strategy in (:equal, :optimized) ||
        throw(ArgumentError("strategy must be :equal or :optimized"));

    # Copy the original bounds, so this scenario never changes the model -
    lower = copy(model.lower);
    upper = copy(model.upper);

    # Read the sizes from the model; the tests also use small, reordered networks -
    n = size(model.S, 2); # number of reactions
    uptake = model.amino_acid_uptake; # columns of the supply reactions
    K = length(uptake); # number of supplied amino acids

    if strategy == :equal

        # Equal allocation: cap each supply at its share of the budget -
        share = budget / K; # mM/h per amino acid
        for j in uptake
            upper[j] = min(upper[j], share); # keep an original bound that is tighter
        end

        # There is no extra inequality: A has zero rows and one column per reaction -
        A = zeros(0, n);
        b = Float64[];

    else

        # Optimized allocation: one inequality, sum of the supply fluxes <= budget -
        A = zeros(1, n);
        for j in uptake
            A[1, j] = 1.0;
        end
        b = [Float64(budget)];
    end

    return (lower=lower, upper=upper, A=A, b=b);
end


# ----------------------------------------------------------------------------
# Part 2: Compare production across supply budgets
# ----------------------------------------------------------------------------

"""
    production_curve(model, budgets::AbstractVector{<:Real}) -> Vector{NamedTuple}

Compute protein production under both allocation strategies at each budget.

# Arguments

- `model`: data accepted by `protein_objective` and `allocation_constraints`,
  including `metabolites` and `reactions` identifiers for the calculation model.
- `budgets`: a vector of finite, nonnegative supply limits in mM/h. Its entries
  may be unsorted or repeated; an empty vector is allowed.

# Returns

A vector of named tuples with fields `(budget, equal, optimized)`, one per input
entry in the original order, including duplicates. The budget is in mM/h;
`equal` and `optimized` are protein production rates in micromolar/h (multiply
the protein flux in mM/h by 1000). An empty input returns an empty vector.
Leave the model and budget vector unchanged.

# Method

Set `c = protein_objective(model)` and obtain `limits` from
`allocation_constraints` for each scenario. Following the urea example, use
`build(MyPrimalFluxBalanceAnalysisCalculationModel, data)` with these fields:

- `S = model.S`
- `fluxbounds = hcat(limits.lower, limits.upper)`
- `objective = c`
- `species = model.metabolites`
- `reactions = model.reactions`

Call `solve(calculation; A=limits.A, b=limits.b)` and check the solver status.
Then read the protein entry of `result["argmax"]` and convert it to μM/h.

Start each strategy from the original bounds, so an equal-allocation cap
cannot restrict the optimized solve.

# Errors

Throw `ArgumentError` for any negative or nonfinite budget. If either solve is
not optimal, throw `ErrorException` with the budget, strategy, and solver status
in its message. Check `result["termination_status"] == JuMP.MOI.OPTIMAL` before
reading a flux or objective value.
"""
function production_curve(model, budgets::AbstractVector{<:Real})

    # Check every budget before solving anything -
    for budget in budgets
        isfinite(budget) && budget >= 0 ||
            throw(ArgumentError("all budgets must be finite and nonnegative"));
    end

    # The protein objective is the same for every scenario -
    c = protein_objective(model);
    rows = NamedTuple[];

    # Solve both strategies at each budget, in the input order -
    for budget in budgets

        rates = Dict{Symbol,Float64}(); # production rate of each strategy, in μM/h

        for strategy in (:equal, :optimized)

            # Each strategy starts from the original model bounds -
            limits = allocation_constraints(model, budget; strategy=strategy);

            # Build the calculation model for this scenario -
            calculation = build(MyPrimalFluxBalanceAnalysisCalculationModel, (
                S = model.S,
                fluxbounds = hcat(limits.lower, limits.upper),
                objective = c,
                species = model.metabolites,
                reactions = model.reactions,
            ));

            # Solve, then check the status before reading any flux -
            result = solve(calculation; A=limits.A, b=limits.b);
            status = result["termination_status"];
            status == JuMP.MOI.OPTIMAL ||
                error("budget=$budget, strategy=$strategy, status=$status");

            # Read the protein flux and convert it from mM/h to μM/h -
            flux = result["argmax"];
            rates[strategy] = 1000.0 * flux[model.protein];
        end

        # Add one row per budget, so a repeated budget gives a repeated row -
        push!(rows, (budget=Float64(budget), equal=rates[:equal], optimized=rates[:optimized]));
    end

    return rows;
end


# ----------------------------------------------------------------------------
# Advanced (Part 3, optional): Find the range of each supply flux
# ----------------------------------------------------------------------------

"""
    supply_ranges(model, budget::Real;
        strategy::Symbol=:equal, fraction::Real=1.0) -> Vector{NamedTuple}

Compute the range of every amino-acid supply flux for one budget and strategy.

# Arguments

- `model`: data accepted by `protein_objective` and `allocation_constraints`,
  including `metabolites` and `reactions` identifiers for the calculation model.
- `budget::Real`: a finite, nonnegative total amino-acid supply limit in mM/h.
- `strategy::Symbol`: `:equal` (default) or `:optimized`.
- `fraction::Real`: the required fraction γ of the maximum protein output, with
  `0 <= fraction <= 1`. The default 1.0 keeps only optimal flux vectors.

# Returns

One named tuple per supply column, in the order of `model.amino_acid_uptake`,
with fields `column` (the column index in `model.S`), `minimum`, and `maximum`.
The endpoints are the smallest and largest values of that supply flux, in mM/h,
over all flux vectors that satisfy the strategy's constraints and reach at least
`fraction` times the maximum protein output. Leave the model unchanged.

Each endpoint comes from its own optimization; the complete set of minima or
maxima need not occur in one feasible flux vector. These ranges describe freedom
within the model, not measurement uncertainty. A zero fraction removes the
positive production requirement; balances and allocation constraints still apply.

# Method

1. Build the limits with `allocation_constraints` and `c = protein_objective(model)`,
   then build a calculation model and call `solve`. Set `z` to
   `result["objective_value"]`, in mM/h.

2. Build `vcat(A, -c')` and `vcat(b, -fraction*z)` from the limits' `A` and `b`.
   The new row requires protein output of at least `fraction*z`; do not relax it
   by a tolerance. `vcat` also converts a `b` that holds integers to `Float64`.

3. For each supply column `j`, solve twice with the strategy's bounds and the new
   arrays. Set the calculation's objective for each solve, either by replacing
   its `objective` field or by building a new calculation model. For the maximum,
   use a vector that is one at `j` and zero elsewhere. For the minimum, use minus
   one at `j` and negate `result["objective_value"]`: the solver always maximizes.
   Keep the original protein vector `c` when forming the production requirement;
   the objectives selecting individual supplies serve a different purpose.

# Errors

Throw `ArgumentError` for an invalid budget or strategy, or for a `fraction` that
is nonfinite or outside [0, 1]. If any solve is not optimal, throw `ErrorException`
with the budget, strategy, and solver status in its message. Check
`result["termination_status"] == JuMP.MOI.OPTIMAL` after every solve before
reading its objective value.
"""
function supply_ranges(model, budget::Real; strategy::Symbol=:equal, fraction::Real=1.0)

    # Check the fraction; allocation_constraints checks the budget and strategy -
    isfinite(fraction) && 0 <= fraction <= 1 ||
        throw(ArgumentError("fraction must be finite and between 0 and 1"));

    # Step 1: solve once for the maximum protein flux z, in mM/h -
    limits = allocation_constraints(model, budget; strategy=strategy);
    c = protein_objective(model);

    calculation = build(MyPrimalFluxBalanceAnalysisCalculationModel, (
        S = model.S,
        fluxbounds = hcat(limits.lower, limits.upper),
        objective = c,
        species = model.metabolites,
        reactions = model.reactions,
    ));

    result = solve(calculation; A=limits.A, b=limits.b);
    status = result["termination_status"];
    status == JuMP.MOI.OPTIMAL ||
        error("budget=$budget, strategy=$strategy, status=$status");
    z = result["objective_value"];

    # Step 2: require c' * v >= fraction * z, written as -c' * v <= -fraction * z -
    A = vcat(limits.A, -c');
    b = vcat(limits.b, -fraction * z); # no tolerance; one would widen the ranges

    # Step 3: find the smallest and largest value of each supply flux -
    rows = NamedTuple[];
    for j in model.amino_acid_uptake

        # This vector selects supply flux j -
        e = zeros(length(c));
        e[j] = 1.0;

        # Minimum: the solver only maximizes, so maximize -v[j] and negate the result -
        calculation.objective = -e;
        result = solve(calculation; A=A, b=b);
        status = result["termination_status"];
        status == JuMP.MOI.OPTIMAL ||
            error("budget=$budget, strategy=$strategy, status=$status");
        lowest = -result["objective_value"];

        # Maximum: maximize v[j] -
        calculation.objective = e;
        result = solve(calculation; A=A, b=b);
        status = result["termination_status"];
        status == JuMP.MOI.OPTIMAL ||
            error("budget=$budget, strategy=$strategy, status=$status");
        highest = result["objective_value"];

        push!(rows, (column=j, minimum=lowest, maximum=highest));
    end

    return rows;
end
