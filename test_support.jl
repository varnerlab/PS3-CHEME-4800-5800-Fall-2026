# Independent checks and small networks used by the scored tests.
"""Supplied test fixtures and feasibility checks for the PS3 student functions."""
module PS3Checks

using ..CellFreeProduction; # exercise the implementation loaded by the current test suite
export toy_model, swap_model, throws_type, curve_close, ranges_close, valid_real_solution;

"""
    toy_model(; protein_first=false, ceiling=100.0) -> NamedTuple

Build a two-species, three-reaction test model whose balances require
`supply1 = 2*protein` and `supply2 = protein`. `ceiling` is the protein-flux
upper bound; tests supply finite, nonnegative values. Each supply is capped at 100.

Return fresh arrays in `S`, `lower`, and `upper`, plus `amino_acid_uptake` and
`protein` column indices. With `protein_first=false`, the columns are supply1,
supply2, protein; with `true`, they are protein, supply1, supply2.
"""
function toy_model(; protein_first=false, ceiling=100.0)
    S = [1.0 0.0 -2.0; 0.0 1.0 -1.0];
    order = protein_first ? [3, 1, 2] : [1, 2, 3];
    return (S=S[:, order], lower=zeros(3), upper=[100.0, 100.0, ceiling][order],
        amino_acid_uptake=protein_first ? [2, 3] : [1, 2],
        protein=protein_first ? 1 : 3);
end

"""
    swap_model(; protein_first=false, ceiling=100.0) -> NamedTuple

Build a two-species, four-reaction test model with interchangeable supplies.
supply1 adds species A, supply2 adds species B, a conversion reaction turns one B
into one A, and protein output consumes one A and one B. The balances are
`supply1 + convert = protein` and `supply2 = convert + protein`. `ceiling` is the
protein-flux upper bound; tests supply finite, nonnegative values. The other
three reactions are capped at 100.

Return fresh arrays in `S`, `lower`, and `upper`, plus `amino_acid_uptake` and
`protein` column indices. With `protein_first=false`, the columns are supply1,
supply2, convert, protein; with `true`, they are protein, supply1, supply2, convert.
"""
function swap_model(; protein_first=false, ceiling=100.0)
    S = [1.0 0.0 1.0 -1.0; 0.0 1.0 -1.0 -1.0];
    order = protein_first ? [4, 1, 2, 3] : [1, 2, 3, 4];
    return (S=S[:, order], lower=zeros(4), upper=[100.0, 100.0, 100.0, ceiling][order],
        amino_acid_uptake=protein_first ? [2, 3] : [1, 2],
        protein=protein_first ? 1 : 4);
end

"""
    throws_type(f, T; contains=String[]) -> Bool

Call the zero-argument function `f`. Return true only if it throws an exception
of type `T` whose printed message contains every string in `contains`. Return
false if no exception is thrown, the type differs, or a required string is absent.
Catch the exception for the test; any side effects of calling `f` still occur.
"""
function throws_type(f, T; contains=String[])
    try
        f();
    catch e
        return e isa T && all(s -> occursin(s, sprint(showerror, e)), contains);
    end
    return false;
end

"""
    curve_close(row, budget, equal, optimized; atol=1e-3) -> Bool

Compare a returned named tuple `row` with the expected budget (mM/h) and two
protein production rates (μM/h). Require the fields `budget`, `equal`, and
`optimized` in that order. Use absolute tolerance 1e-10 for the budget and
`atol` for each rate, with zero relative tolerance. Return whether every check
passes. Inputs are unchanged; malformed rows can raise field-access errors.
"""
function curve_close(row, budget, equal, optimized; atol=1e-3)
    return keys(row) == (:budget, :equal, :optimized) &&
        isapprox(row.budget, budget; atol=1e-10, rtol=0) &&
        isapprox(row.equal, equal; atol=atol, rtol=0) &&
        isapprox(row.optimized, optimized; atol=atol, rtol=0);
end

"""
    ranges_close(rows, columns, expected; atol=1e-8) -> Bool

Compare the rows returned by `supply_ranges` with expected supply columns and
`(minimum, maximum)` pairs in mM/h. Require one row per expected column, in the
same order, with the fields `column`, `minimum`, and `maximum` in that order.
Use absolute tolerance `atol` and zero relative tolerance. Return whether every
check passes. Inputs are unchanged; malformed rows can raise field-access errors.
"""
function ranges_close(rows, columns, expected; atol=1e-8)
    return length(rows) == length(columns) == length(expected) &&
        all(keys(row) == (:column, :minimum, :maximum) for row in rows) &&
        [row.column for row in rows] == columns &&
        all(isapprox(row.minimum, lo; atol=atol, rtol=0) && isapprox(row.maximum, hi; atol=atol, rtol=0)
            for (row, (lo, hi)) in zip(rows, expected));
end

"""
    valid_real_solution(model, budget, strategy) -> Bool

Use the student's objective and allocation functions to solve the supplied
20-amino-acid deGFP model at `budget` mM/h with `:equal` or `:optimized` allocation.
Return false for a nonoptimal solve. Otherwise check its flux vector directly
against `model.S`, the original model bounds, and the stated supply-budget rules,
using an absolute tolerance of 1e-7 mM/h. The supplied test inputs are valid;
exceptions from the student functions or solver propagate to the test framework.
"""
function valid_real_solution(model, budget, strategy)
    limits = allocation_constraints(model, budget; strategy=strategy);
    result = solve_lp(model.S, limits.lower, limits.upper, protein_objective(model);
        A=limits.A, b=limits.b);
    result.optimal || return false;
    v = result.flux;
    supplies = v[model.amino_acid_uptake];
    common = maximum(abs, model.S * v) <= 1e-7 &&
        all(model.lower .- 1e-7 .<= v) && all(v .<= model.upper .+ 1e-7) &&
        all(supplies .>= -1e-7) && sum(supplies) <= budget + 1e-7;
    return common && (strategy == :optimized || all(supplies .<= budget / 20 + 1e-7));
end

end
