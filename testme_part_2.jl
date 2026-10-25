# Run: julia --startup-file=no --project=. testme_part_2.jl
# Add --solution to test solution/src/Compute.jl without replacing student code.
# Each labeled case checks a production-curve requirement or solution feasibility.


# ----------------------------------------------------------------------------
# Setup | Load the selected implementation and test helpers
# ----------------------------------------------------------------------------

using Test;  # record one scored assertion per labeled case

Base.include(@__MODULE__, joinpath(@__DIR__, "Include.jl"));
Base.include(@__MODULE__, joinpath(@__DIR__, "test_support.jl"));

using .PS3Checks; # small networks and independent result checks

model = load_model();


# ----------------------------------------------------------------------------
# Run the 24 scored checks | One labeled case, one @test each
# ----------------------------------------------------------------------------
# A case whose code raises an error counts as one failed test.

@testset "PS3 Part 2 (24 tests)" begin

    # Output layout, units, and input preservation -
    @testset "toy production rates and unit conversion" begin
        rows = production_curve(toy_model(), [6.0]);
        @test curve_close(only(rows), 6, 1500, 2000);
    end

    @testset "input budget order retained" begin
        rows = production_curve(toy_model(), [6.0, 0.0, 3.0]);

        # Each row must match its own budget; stop checking at the first mismatch -
        all_close = length(rows) == 3;
        for (row, b) in zip(rows, [6, 0, 3])
            all_close = all_close && curve_close(row, b, 250b, 1000b/3);
        end
        @test all_close;
    end

    @testset "duplicate budgets retained" begin
        rows = production_curve(toy_model(), [3.0, 3.0]);
        @test length(rows) == 2 && rows[1] == rows[2];
    end

    @testset "empty budget vector" begin
        rows = production_curve(model, Float64[]);
        @test rows isa AbstractVector && isempty(rows);
    end

    @testset "model and budget vector preserved" begin
        m = deepcopy(model);
        before = deepcopy(m);
        budgets = [1.0, 0.0];
        production_curve(m, budgets);
        @test isequal(m, before) && budgets == [1.0, 0.0];
    end

    # Production rates for the supplied deGFP model -
    @testset "deGFP production without external amino acids" begin
        rows = production_curve(model, [0]);
        @test curve_close(only(rows), 0, 2.6611, 2.6611);
    end

    @testset "equal production at 0.5" begin
        rows = production_curve(model, [0.5]);
        @test isapprox(only(rows).equal, 4.4197; atol=1e-3);
    end

    @testset "optimized production at 0.5" begin
        rows = production_curve(model, [0.5]);
        @test isapprox(only(rows).optimized, 5.3399; atol=1e-3);
    end

    @testset "equal production at 1.0" begin
        rows = production_curve(model, [1]);
        @test isapprox(only(rows).equal, 6.0035; atol=1e-3);
    end

    @testset "optimized production at 1.0" begin
        rows = production_curve(model, [1]);
        @test isapprox(only(rows).optimized, 7.4432; atol=1e-3);
    end

    @testset "equal production at 2.0" begin
        rows = production_curve(model, [2]);
        @test isapprox(only(rows).equal, 8.9265; atol=1e-3);
    end

    @testset "optimized production at 2.0" begin
        rows = production_curve(model, [2]);
        @test isapprox(only(rows).optimized, 11.1758; atol=1e-3);
    end

    @testset "both strategies reach the expression limit" begin
        rows = production_curve(model, [4]);
        @test curve_close(only(rows), 4, 11.1758, 11.1758);
    end

    # Curve behavior and changed model conditions -
    @testset "production is nondecreasing" begin
        rows = production_curve(model, collect(0:0.5:4));

        # Collect each strategy's rates, then compare neighboring budgets -
        equal = Float64[];
        optimized = Float64[];
        for r in rows
            push!(equal, r.equal);
            push!(optimized, r.optimized);
        end
        @test all(diff(equal) .>= -1e-5) && all(diff(optimized) .>= -1e-5);
    end

    @testset "optimized production is at least equal production" begin
        rows = production_curve(model, collect(0:0.5:4));

        # Stop checking at the first budget where optimized falls below equal -
        at_least_equal = true;
        for r in rows
            at_least_equal = at_least_equal && r.optimized >= r.equal - 1e-5;
        end
        @test at_least_equal;
    end

    @testset "tighter protein bound is respected" begin
        rows = production_curve(toy_model(ceiling=0.4), [6]);
        @test curve_close(only(rows), 6, 400, 400);
    end

    @testset "curve follows a changed protein column" begin
        rows = production_curve(toy_model(protein_first=true), [6]);
        @test curve_close(only(rows), 6, 1500, 2000);
    end

    @testset "different glucose condition changes production" begin
        rows = production_curve(load_model(glucose_limit=0.5), [1]);
        @test curve_close(only(rows), 1, 3.3952, 4.7493);
    end

    # Invalid inputs and solver failures -
    @testset "negative budget rejected" begin
        @test_throws ArgumentError production_curve(model, [-1.0]);
    end

    @testset "infinite budget rejected" begin
        @test_throws ArgumentError production_curve(model, [Inf]);
    end

    @testset "NaN budget rejected" begin
        @test_throws ArgumentError production_curve(model, [NaN]);
    end

    @testset "nonoptimal status is handled before reading fluxes" begin
        m = toy_model();
        m.lower[m.protein] = 1;

        # Either strategy may be solved first; the message names the one that failed -
        message = "";
        try
            production_curve(m, [0]);
        catch caught
            if caught isa ErrorException
                message = sprint(showerror, caught);
            end
        end
        @test occursin("0", message) && occursin(r"equal|optimized", message) &&
            occursin("INFEASIBLE", message);
    end

    @testset "invalid later budget rejected" begin
        @test_throws ArgumentError production_curve(model, [1.0, -0.1]);
    end

    # Independent feasibility checks -
    @testset "real allocations satisfy balances, bounds, and budgets" begin

        # Stop checking at the first scenario that fails -
        all_valid = true;
        for b in (0.0, 0.5, 1.0, 2.0, 4.0)
            for strategy in (:equal, :optimized)
                all_valid = all_valid && valid_real_solution(model, b, strategy);
            end
        end
        @test all_valid;
    end
end;
