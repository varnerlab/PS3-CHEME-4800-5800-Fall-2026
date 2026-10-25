# Run: julia --startup-file=no --project=. testme_part_3.jl
# Add --solution to test solution/src/Compute.jl without replacing student code.
# Each labeled case checks a flux-variability requirement or a range on the deGFP model.


# ----------------------------------------------------------------------------
# Setup | Load the selected implementation and test helpers
# ----------------------------------------------------------------------------

import JuMP; # solver termination status
using Test;  # record one scored assertion per labeled case

Base.include(@__MODULE__, joinpath(@__DIR__, "Include.jl"));
Base.include(@__MODULE__, joinpath(@__DIR__, "test_support.jl"));

using .PS3Checks; # small networks and independent result checks

model = load_model();


"""
    range_of(model, rows, name) -> Tuple{Float64,Float64}

Return the `(minimum, maximum)` pair, in mM/h, for the amino acid called `name`.
`rows` comes from `supply_ranges` and follows the order of `model.amino_acids`.
"""
function range_of(model, rows, name)

    for k in 1:length(model.amino_acids)
        if model.amino_acids[k] == name
            return (rows[k].minimum, rows[k].maximum);
        end
    end

    error("unknown amino acid: $name");
end


"""
    close_to(pair, lo, hi) -> Bool

Return true when both endpoints of `pair` are within 1e-5 mM/h of `lo` and `hi`.
"""
function close_to(pair, lo, hi)
    return isapprox(pair[1], lo; atol=1e-5, rtol=0) &&
        isapprox(pair[2], hi; atol=1e-5, rtol=0);
end


"""
    fixed(row) -> Bool

Return true when a supply range is no wider than 1e-5 mM/h.
"""
function fixed(row)
    return row.maximum - row.minimum <= 1e-5;
end


# ----------------------------------------------------------------------------
# Run the 24 scored checks | One labeled case, one @test each
# ----------------------------------------------------------------------------
# A case whose code raises an error counts as one failed test.

@testset "PS3 Part 3 (24 tests)" begin

    # Hand-checkable ranges on the small networks -
    @testset "toy optimized ranges split one shared total" begin
        rows = supply_ranges(swap_model(), 6; strategy=:optimized);
        @test ranges_close(rows, [1, 2], [(0, 3), (3, 6)]);
    end

    @testset "toy equal ranges are fixed" begin
        rows = supply_ranges(swap_model(), 6; strategy=:equal);
        @test ranges_close(rows, [1, 2], [(3, 3), (3, 3)]);
    end

    @testset "default strategy is equal" begin
        rows = supply_ranges(swap_model(), 6);
        @test ranges_close(rows, [1, 2], [(3, 3), (3, 3)]);
    end

    @testset "toy optimized near-optimal ranges" begin
        rows = supply_ranges(swap_model(), 6; strategy=:optimized, fraction=0.5);
        @test ranges_close(rows, [1, 2], [(0, 3), (1.5, 6)]);
    end

    @testset "toy equal near-optimal ranges" begin
        rows = supply_ranges(swap_model(), 6; fraction=0.5);
        @test ranges_close(rows, [1, 2], [(0, 3), (1.5, 3)]);
    end

    @testset "zero fraction gives feasibility ranges" begin
        rows = supply_ranges(swap_model(), 6; strategy=:optimized, fraction=0);
        @test ranges_close(rows, [1, 2], [(0, 3), (0, 6)]);
    end

    @testset "ranges follow a changed column order" begin
        rows = supply_ranges(swap_model(protein_first=true), 6; strategy=:optimized);
        @test ranges_close(rows, [2, 3], [(0, 3), (3, 6)]);
    end

    @testset "tighter protein bound is respected" begin
        rows = supply_ranges(swap_model(ceiling=1), 6; strategy=:optimized);
        @test ranges_close(rows, [1, 2], [(0, 1), (1, 2)]);
    end

    # Input preservation and repeated calls -
    @testset "model is preserved" begin
        m = deepcopy(model);
        before = deepcopy(m);
        supply_ranges(m, 1.0; strategy=:optimized);
        supply_ranges(m, 1.0; fraction=0.99);
        @test isequal(m, before);
    end

    @testset "repeated calls agree" begin
        m = swap_model();
        initial = supply_ranges(m, 6; strategy=:optimized);
        supply_ranges(m, 0; fraction=0.5);
        @test initial == supply_ranges(m, 6; strategy=:optimized);
    end

    # Invalid inputs and solver failures -
    @testset "negative budget rejected" begin
        @test_throws ArgumentError supply_ranges(model, -1.0);
    end

    @testset "infinite budget rejected" begin
        @test_throws ArgumentError supply_ranges(model, Inf; strategy=:optimized);
    end

    @testset "unknown strategy rejected" begin
        @test_throws ArgumentError supply_ranges(model, 1.0; strategy=:other);
    end

    @testset "invalid fractions rejected" begin

        # Every invalid fraction must raise an ArgumentError; stop at the first that does not -
        rejected = true;
        for f in (1.5, -0.1, NaN, Inf)
            try
                supply_ranges(swap_model(), 6; fraction=f);
                rejected = false; # no error was raised
            catch caught
                rejected = caught isa ArgumentError;
            end
            if !rejected
                break;
            end
        end
        @test rejected;
    end

    @testset "nonoptimal status is handled before reading fluxes" begin
        m = swap_model();
        m.lower[m.protein] = 1;

        # The error message must name the budget, strategy, and status -
        message = "";
        try
            supply_ranges(m, 0; strategy=:optimized);
        catch caught
            if caught isa ErrorException
                message = sprint(showerror, caught);
            end
        end
        @test occursin("0", message) && occursin("optimized", message) &&
            occursin("INFEASIBLE", message);
    end

    # Reference ranges for the deGFP model -
    @testset "deGFP equal supplies are fixed at 1.0" begin
        rows = supply_ranges(model, 1.0);
        @test all(fixed, rows) &&
            close_to(range_of(model, rows, "Cysteine"), 0.012007, 0.012007) &&
            close_to(range_of(model, rows, "Leucine"), 0.05, 0.05) &&
            close_to(range_of(model, rows, "Tryptophan"), 0.006004, 0.006004);
    end

    @testset "deGFP optimized fixed supplies at 1.0" begin
        rows = supply_ranges(model, 1.0; strategy=:optimized);
        @test close_to(range_of(model, rows, "Leucine"), 0.141421, 0.141421) &&
            close_to(range_of(model, rows, "Lysine"), 0.133977, 0.133977) &&
            close_to(range_of(model, rows, "Arginine"), 0.044659, 0.044659);
    end

    @testset "deGFP supplies that are never used at 1.0" begin
        rows = supply_ranges(model, 1.0; strategy=:optimized);

        # Stop checking at the first amino acid whose range is not zero -
        unused = true;
        for a in ("Alanine", "Asparagine", "Aspartate", "Cysteine", "Glycine", "Serine")
            unused = unused && close_to(range_of(model, rows, a), 0, 0);
        end
        @test unused;
    end

    @testset "deGFP interchangeable supplies at 1.0" begin
        rows = supply_ranges(model, 1.0; strategy=:optimized);
        @test close_to(range_of(model, rows, "Glutamate"), 0, 0.054715) &&
            close_to(range_of(model, rows, "Glutamine"), 0.059545, 0.114261) &&
            close_to(range_of(model, rows, "Threonine"), 0, 0.054715) &&
            count(fixed, rows) == 17;
    end

    @testset "deGFP ranges open at translation capacity" begin
        rows = supply_ranges(model, 3.0; strategy=:optimized);

        # Every minimum must be zero; stop checking at the first that is not -
        all_zero = true;
        for r in rows
            all_zero = all_zero && isapprox(r.minimum, 0; atol=1e-5);
        end
        @test all_zero && close_to(range_of(model, rows, "Glutamate"), 0, 3.0);
    end

    @testset "deGFP near-optimal ranges at 1.0" begin
        rows = supply_ranges(model, 1.0; strategy=:optimized, fraction=0.99);
        @test !any(fixed, rows) && close_to(range_of(model, rows, "Valine"), 0.041459, 0.126534);
    end

    @testset "deGFP equal ranges open at translation capacity" begin
        rows = supply_ranges(model, 3.0);
        @test close_to(range_of(model, rows, "Leucine"), 0.112226, 0.15);
    end

    @testset "different glucose condition changes ranges" begin
        rows = supply_ranges(load_model(glucose_limit=0.5), 1.0; strategy=:optimized);
        @test close_to(range_of(model, rows, "Glutamate"), 0, 0.396840);
    end

    # Consistency between an optimal solution and its ranges -
    @testset "ranges contain the reported optimal allocations" begin

        # For each strategy, every supply flux of one optimal solution must lie
        # inside its range; stop checking at the first that does not -
        contained = true;
        for strategy in (:equal, :optimized)
            limits = allocation_constraints(model, 1.0; strategy=strategy);

            calculation = build(MyPrimalFluxBalanceAnalysisCalculationModel, (
                S = model.S,
                fluxbounds = hcat(limits.lower, limits.upper),
                objective = protein_objective(model),
                species = model.metabolites,
                reactions = model.reactions,
            ));

            result = solve(calculation; A=limits.A, b=limits.b);
            rows = supply_ranges(model, 1.0; strategy=strategy);

            contained = result["termination_status"] == JuMP.MOI.OPTIMAL;
            for r in rows
                contained = contained &&
                    r.minimum - 1e-6 <= result["argmax"][r.column] <= r.maximum + 1e-6;
            end

            if !contained
                break;
            end
        end
        @test contained;
    end
end;
