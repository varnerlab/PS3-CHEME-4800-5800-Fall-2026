# Run: julia --startup-file=no --project=. testme_part_1.jl
# Add --solution to test solution/src/Compute.jl without replacing student code.
# Each labeled case checks one objective or allocation requirement.


# ----------------------------------------------------------------------------
# Setup | Load the selected implementation and test helpers
# ----------------------------------------------------------------------------

import JuMP; # solver termination status
using Test;  # record one scored assertion per labeled case

Base.include(@__MODULE__, joinpath(@__DIR__, "Include.jl"));
Base.include(@__MODULE__, joinpath(@__DIR__, "test_support.jl"));

using .PS3Checks; # small networks and independent result checks

model = load_model();


# ----------------------------------------------------------------------------
# Run the 24 scored checks | One labeled case, one @test each
# ----------------------------------------------------------------------------
# A case whose code raises an error counts as one failed test.

@testset "PS3 Part 1 (24 tests)" begin

    # Protein objective -
    @testset "objective dimensions and type" begin
        c = protein_objective(model);
        @test c isa Vector{Float64} && length(c) == 265;
    end

    @testset "objective selects only protein" begin
        c = protein_objective(model);
        @test c[model.protein] == 1 && count(c .!= 0) == 1;
    end

    @testset "objective follows a changed column order" begin
        @test protein_objective(toy_model(protein_first=true)) == [1, 0, 0];
    end

    @testset "objective preserves the input" begin
        m = deepcopy(model);
        before = deepcopy(m);
        protein_objective(m);
        @test isequal(m, before);
    end

    # Equal-allocation bounds -
    @testset "equal bounds are independent copies" begin
        m = deepcopy(model);
        a = allocation_constraints(m, 1.0);
        a.lower[1] = -99;
        a.upper[1] = -99;
        @test m.lower == model.lower && m.upper == model.upper;
    end

    @testset "equal supply limits" begin
        a = allocation_constraints(model, 1.0);
        @test all(a.upper[model.amino_acid_uptake] .== 0.05);
    end

    @testset "equal preserves lower bounds" begin
        @test allocation_constraints(model, 1.0).lower == model.lower;
    end

    @testset "equal preserves other upper bounds" begin
        j = setdiff(1:265, model.amino_acid_uptake);
        @test allocation_constraints(model, 1.0).upper[j] == model.upper[j];
    end

    @testset "equal needs no extra inequality" begin
        a = allocation_constraints(model, 1.0);
        @test size(a.A) == (0, 265) && isempty(a.b);
    end

    @testset "equal hand calculation on two supplies" begin
        m = toy_model();
        a = allocation_constraints(m, 6.0);

        calculation = build(MyPrimalFluxBalanceAnalysisCalculationModel, (
            S = m.S,
            fluxbounds = hcat(a.lower, a.upper),
            objective = protein_objective(m),
            species = m.metabolites,
            reactions = m.reactions,
        ));

        r = solve(calculation; A=a.A, b=a.b);
        @test r["termination_status"] == JuMP.MOI.OPTIMAL &&
            isapprox(r["objective_value"], 1.5; atol=1e-8);
    end

    @testset "zero equal budget blocks supply" begin
        a = allocation_constraints(model, 0);
        @test all(iszero, a.upper[model.amino_acid_uptake]);
    end

    @testset "equal retains tighter original bounds" begin
        @test allocation_constraints(model, 10000).upper == model.upper;
    end

    # Optimized-allocation constraint -
    @testset "optimized retains original bounds" begin
        a = allocation_constraints(model, 1; strategy=:optimized);
        @test a.lower == model.lower && a.upper == model.upper;
    end

    @testset "optimized inequality selects precisely the supply columns" begin
        a = allocation_constraints(model, 1; strategy=:optimized);
        @test size(a.A) == (1, 265) && all(a.A[1, model.amino_acid_uptake] .== 1) &&
            count(a.A .!= 0) == 20;
    end

    @testset "optimized right-hand side uses the supplied budget" begin
        @test allocation_constraints(model, 2.75; strategy=:optimized).b == [2.75];
    end

    @testset "optimized hand calculation on two supplies" begin
        m = toy_model();
        a = allocation_constraints(m, 6.0; strategy=:optimized);

        calculation = build(MyPrimalFluxBalanceAnalysisCalculationModel, (
            S = m.S,
            fluxbounds = hcat(a.lower, a.upper),
            objective = protein_objective(m),
            species = m.metabolites,
            reactions = m.reactions,
        ));

        r = solve(calculation; A=a.A, b=a.b);
        @test r["termination_status"] == JuMP.MOI.OPTIMAL &&
            isapprox(r["objective_value"], 2.0; atol=1e-8);
    end

    @testset "optimized bounds are independent copies" begin
        m = deepcopy(model);
        a = allocation_constraints(m, 1; strategy=:optimized);
        a.lower[1] = -99;
        a.upper[1] = -99;
        @test m.lower == model.lower && m.upper == model.upper;
    end

    # Input preservation and invalid arguments -
    @testset "strategies preserve the complete model" begin
        m = deepcopy(model);
        before = deepcopy(m);
        allocation_constraints(m, 0.5);
        allocation_constraints(m, 2; strategy=:optimized);
        @test isequal(m, before);
    end

    @testset "negative budget rejected" begin
        @test_throws ArgumentError allocation_constraints(model, -1);
    end

    @testset "infinite budget rejected" begin
        @test_throws ArgumentError allocation_constraints(model, Inf; strategy=:optimized);
    end

    @testset "NaN budget rejected" begin
        @test_throws ArgumentError allocation_constraints(model, NaN);
    end

    @testset "unknown strategy rejected" begin
        @test_throws ArgumentError allocation_constraints(model, 1; strategy=:other);
    end

    # Repeated calls and a zero shared budget -
    @testset "repeated equal scenarios start from original bounds" begin
        m = deepcopy(model);
        allocation_constraints(m, 0);
        a = allocation_constraints(m, 2);
        @test all(a.upper[m.amino_acid_uptake] .== 0.1);
    end

    @testset "optimized zero budget is a zero shared allowance" begin
        a = allocation_constraints(toy_model(protein_first=true), 0; strategy=:optimized);
        @test a.A == [0.0 1.0 1.0] && a.b == [0.0];
    end
end;
