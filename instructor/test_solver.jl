# Supplied solver checks, separate from the 72 scored student tests.
using Test; # interface and feasibility assertions
import JuMP; # termination-status constants
include(joinpath(@__DIR__, "..", "Include.jl"));

@testset "Course-style FBA interface" begin
    # Uptake-positive exchanges make export negative, as in the urea example -
    data = (S=[1.0 1.0], fluxbounds=[0.0 3.0; -10.0 0.0],
        objective=[0.0, -1.0], species=["A"], reactions=["uptake", "export"]);
    calculation = build(MyPrimalFluxBalanceAnalysisCalculationModel, data);
    @test calculation isa AbstractFluxCalculationModel;
    @test calculation.S === data.S && calculation.objective === data.objective;
    before = deepcopy(data);
    result = solve(calculation);
    @test Set(keys(result)) == Set(["argmax", "objective_value", "termination_status"]);
    @test result["termination_status"] == JuMP.MOI.OPTIMAL;
    @test result["argmax"] ≈ [3.0, -3.0];
    @test result["objective_value"] ≈ 3.0;
    @test all(isequal(getfield(calculation, key), getfield(before, key)) for key in keys(data));

    # Extra inequalities can limit production or make a scenario infeasible -
    A, b = [1.0 0.0], [2.0];
    constrained = solve(calculation; A=A, b=b);
    @test constrained["objective_value"] ≈ 2.0;
    @test A == [1.0 0.0] && b == [2.0];
    limits = (lower=data.fluxbounds[:, 1], upper=data.fluxbounds[:, 2], A=A, b=b);
    @test check_solution(calculation, constrained, limits).valid;
    infeasible = solve(calculation; A=A, b=[-1.0]);
    @test infeasible["termination_status"] == JuMP.MOI.INFEASIBLE;
    @test isnothing(infeasible["argmax"]) && isnothing(infeasible["objective_value"]);
    @test_throws ArgumentError check_solution(calculation, infeasible, limits);
    @test_throws DimensionMismatch solve(calculation; A=zeros(1, 3), b=[1.0]);
    @test_throws ArgumentError solve(calculation; A=A, b=[NaN]);

    # Field updates select a new objective, just as in the course notebook -
    calculation.objective = [0.0, 1.0];
    @test solve(calculation)["objective_value"] ≈ 0.0;
    calculation.fluxbounds = zeros(2, 1);
    @test_throws DimensionMismatch solve(calculation);
    calculation.fluxbounds = [4.0 3.0; -10.0 0.0];
    @test_throws ArgumentError solve(calculation);
    calculation.fluxbounds = copy(data.fluxbounds);
    calculation.reactions = ["missing column"];
    @test_throws DimensionMismatch solve(calculation);
end
