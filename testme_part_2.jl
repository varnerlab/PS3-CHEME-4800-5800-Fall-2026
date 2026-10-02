# Run: julia --startup-file=no --project=. testme_part_2.jl
# Each labeled case checks a production-curve requirement or solution feasibility.
using Test; # record one scored assertion per labeled case
Base.include(@__MODULE__, joinpath(@__DIR__, "Include.jl"));
Base.include(@__MODULE__, joinpath(@__DIR__, "test_support.jl"));
using .PS3Checks; # small networks and independent result checks

model = load_model();
checks = [
    "toy production rates and unit conversion" => () -> curve_close(only(production_curve(toy_model(), [6.0])), 6, 1500, 2000),
    "input budget order retained" => () -> begin
        rows = production_curve(toy_model(), [6.0, 0.0, 3.0]);
        length(rows) == 3 && all(curve_close(row, b, 250b, 1000b/3) for (row, b) in zip(rows, [6, 0, 3]))
    end,
    "duplicate budgets retained" => () -> begin
        rows = production_curve(toy_model(), [3.0, 3.0]); length(rows) == 2 && rows[1] == rows[2]
    end,
    "empty budget vector" => () -> begin
        rows = production_curve(model, Float64[]); rows isa AbstractVector && isempty(rows)
    end,
    "model and budget vector preserved" => () -> begin
        m = deepcopy(model); before = deepcopy(m); budgets = [1.0, 0.0];
        production_curve(m, budgets); isequal(m, before) && budgets == [1.0, 0.0]
    end,
    "deGFP production without external amino acids" => () -> curve_close(only(production_curve(model, [0])), 0, 2.6611, 2.6611),
    "equal production at 0.5" => () -> isapprox(only(production_curve(model, [0.5])).equal, 4.4197; atol=1e-3),
    "optimized production at 0.5" => () -> isapprox(only(production_curve(model, [0.5])).optimized, 5.3399; atol=1e-3),
    "equal production at 1.0" => () -> isapprox(only(production_curve(model, [1])).equal, 6.0035; atol=1e-3),
    "optimized production at 1.0" => () -> isapprox(only(production_curve(model, [1])).optimized, 7.4432; atol=1e-3),
    "equal production at 2.0" => () -> isapprox(only(production_curve(model, [2])).equal, 8.9265; atol=1e-3),
    "optimized production at 2.0" => () -> isapprox(only(production_curve(model, [2])).optimized, 11.1758; atol=1e-3),
    "both strategies reach the expression limit" => () -> curve_close(only(production_curve(model, [4])), 4, 11.1758, 11.1758),
    "production is nondecreasing" => () -> begin
        rows = production_curve(model, collect(0:0.5:4));
        all(diff([r.equal for r in rows]) .>= -1e-5) && all(diff([r.optimized for r in rows]) .>= -1e-5)
    end,
    "optimized production is at least equal production" => () -> all(r.optimized >= r.equal - 1e-5 for r in production_curve(model, collect(0:0.5:4))),
    "tighter protein bound is respected" => () -> curve_close(only(production_curve(toy_model(ceiling=0.4), [6])), 6, 400, 400),
    "curve follows a changed protein column" => () -> curve_close(only(production_curve(toy_model(protein_first=true), [6])), 6, 1500, 2000),
    "different glucose condition changes production" => () -> curve_close(only(production_curve(load_model(glucose_limit=0.5), [1])), 1, 3.3952, 4.7493),
    "negative budget rejected" => () -> throws_type(() -> production_curve(model, [-1.0]), ArgumentError),
    "infinite budget rejected" => () -> throws_type(() -> production_curve(model, [Inf]), ArgumentError),
    "NaN budget rejected" => () -> throws_type(() -> production_curve(model, [NaN]), ArgumentError),
    "nonoptimal status is handled before reading fluxes" => () -> begin
        m = toy_model(); m.lower[m.protein] = 1;
        throws_type(() -> production_curve(m, [0]), ErrorException; contains=["0", "equal", "INFEASIBLE"])
    end,
    "invalid later budget rejected" => () -> throws_type(() -> production_curve(model, [1.0, -0.1]), ArgumentError),
    "real allocations satisfy balances, bounds, and budgets" => () -> all(
        valid_real_solution(model, b, strategy) for b in (0.0, 0.5, 1.0, 2.0, 4.0) for strategy in (:equal, :optimized)),
];

@testset "PS3 Part 2 (24 tests)" begin
    for (label, check) in checks
        @testset "$label" begin
            @test check();
        end
    end
end;
