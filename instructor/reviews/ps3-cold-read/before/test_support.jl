# Independent checks and small networks used by the scored tests.
module PS3Checks

using ..CellFreeProduction;
export toy_model, throws_type, curve_close, valid_real_solution;

"""Make a two-supply model with balances supply1=2*protein and supply2=protein."""
function toy_model(; protein_first=false, ceiling=100.0)
    S = [1.0 0.0 -2.0; 0.0 1.0 -1.0];
    order = protein_first ? [3, 1, 2] : [1, 2, 3];
    return (S=S[:, order], lower=zeros(3), upper=[100.0, 100.0, ceiling][order],
        amino_acid_uptake=protein_first ? [2, 3] : [1, 2],
        protein=protein_first ? 1 : 3);
end

"""Return whether f raises the requested exception; optionally check its message."""
function throws_type(f, T; contains=String[])
    try
        f();
    catch e
        return e isa T && all(s -> occursin(s, sprint(showerror, e)), contains);
    end
    return false;
end

"""Check a production row's budget and rates against independent expected values."""
function curve_close(row, budget, equal, optimized; atol=1e-3)
    return keys(row) == (:budget, :equal, :optimized) &&
        isapprox(row.budget, budget; atol=1e-10, rtol=0) &&
        isapprox(row.equal, equal; atol=atol, rtol=0) &&
        isapprox(row.optimized, optimized; atol=atol, rtol=0);
end

"""Check an actual solution directly against balances, base bounds, and budget rules."""
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
