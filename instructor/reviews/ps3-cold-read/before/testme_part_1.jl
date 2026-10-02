using Test;
Base.include(@__MODULE__, joinpath(@__DIR__, "Include.jl"));
Base.include(@__MODULE__, joinpath(@__DIR__, "test_support.jl"));
using .PS3Checks;

# Each case contributes one scored test, including cases whose implementation raises an error.
model = load_model();
checks = [
    "objective dimensions and type" => () -> begin
        c = protein_objective(model); c isa Vector{Float64} && length(c) == 265
    end,
    "objective selects only protein" => () -> begin
        c = protein_objective(model); c[model.protein] == 1 && count(!iszero, c) == 1
    end,
    "objective follows a changed column order" => () -> protein_objective(toy_model(protein_first=true)) == [1, 0, 0],
    "objective preserves the input" => () -> begin
        m = deepcopy(model); before = deepcopy(m); protein_objective(m); isequal(m, before)
    end,
    "equal bounds are independent copies" => () -> begin
        m = deepcopy(model); a = allocation_constraints(m, 1.0); a.lower[1] = -99; a.upper[1] = -99;
        m.lower == model.lower && m.upper == model.upper
    end,
    "equal supply limits" => () -> all(allocation_constraints(model, 1.0).upper[model.amino_acid_uptake] .== 0.05),
    "equal preserves lower bounds" => () -> allocation_constraints(model, 1.0).lower == model.lower,
    "equal preserves other upper bounds" => () -> begin
        j = setdiff(1:265, model.amino_acid_uptake);
        allocation_constraints(model, 1.0).upper[j] == model.upper[j]
    end,
    "equal needs no extra inequality" => () -> begin
        a = allocation_constraints(model, 1.0); size(a.A) == (0, 265) && isempty(a.b)
    end,
    "equal hand calculation on two supplies" => () -> begin
        m = toy_model(); a = allocation_constraints(m, 6.0);
        r = solve_lp(m.S, a.lower, a.upper, protein_objective(m); A=a.A, b=a.b);
        r.optimal && isapprox(r.objective, 1.5; atol=1e-8)
    end,
    "zero equal budget blocks supply" => () -> all(iszero, allocation_constraints(model, 0).upper[model.amino_acid_uptake]),
    "equal retains tighter original bounds" => () -> allocation_constraints(model, 10000).upper == model.upper,
    "optimized retains original bounds" => () -> begin
        a = allocation_constraints(model, 1; strategy=:optimized); a.lower == model.lower && a.upper == model.upper
    end,
    "optimized inequality selects precisely the supply columns" => () -> begin
        a = allocation_constraints(model, 1; strategy=:optimized);
        size(a.A) == (1, 265) && all(a.A[1, model.amino_acid_uptake] .== 1) &&
            count(!iszero, a.A) == 20
    end,
    "optimized right-hand side uses the supplied budget" => () -> allocation_constraints(model, 2.75; strategy=:optimized).b == [2.75],
    "optimized hand calculation on two supplies" => () -> begin
        m = toy_model(); a = allocation_constraints(m, 6.0; strategy=:optimized);
        r = solve_lp(m.S, a.lower, a.upper, protein_objective(m); A=a.A, b=a.b);
        r.optimal && isapprox(r.objective, 2.0; atol=1e-8)
    end,
    "optimized bounds are independent copies" => () -> begin
        m = deepcopy(model); a = allocation_constraints(m, 1; strategy=:optimized);
        a.lower[1] = -99; a.upper[1] = -99; m.lower == model.lower && m.upper == model.upper
    end,
    "strategies preserve the complete model" => () -> begin
        m = deepcopy(model); before = deepcopy(m);
        allocation_constraints(m, 0.5); allocation_constraints(m, 2; strategy=:optimized); isequal(m, before)
    end,
    "negative budget rejected" => () -> throws_type(() -> allocation_constraints(model, -1), ArgumentError),
    "infinite budget rejected" => () -> throws_type(() -> allocation_constraints(model, Inf; strategy=:optimized), ArgumentError),
    "NaN budget rejected" => () -> throws_type(() -> allocation_constraints(model, NaN), ArgumentError),
    "unknown strategy rejected" => () -> throws_type(() -> allocation_constraints(model, 1; strategy=:other), ArgumentError),
    "repeated equal scenarios start from original bounds" => () -> begin
        m = deepcopy(model); allocation_constraints(m, 0);
        all(allocation_constraints(m, 2).upper[m.amino_acid_uptake] .== 0.1)
    end,
    "optimized zero budget is a zero shared allowance" => () -> begin
        a = allocation_constraints(toy_model(protein_first=true), 0; strategy=:optimized);
        a.A == [0.0 1.0 1.0] && a.b == [0.0]
    end,
];

@testset "PS3 Part 1 (24 tests)" begin
    for (label, check) in checks
        @testset "$label" begin
            @test check();
        end
    end
end;
