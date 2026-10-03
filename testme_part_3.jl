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

# Return the (minimum, maximum) pair for one named amino acid -
range_of(rows, name) = (
    r = rows[findfirst(==(name), model.amino_acids)];
    (r.minimum, r.maximum)
);

# Compare both endpoints in mM/h -
close_to(pair, lo, hi) =
    isapprox(pair[1], lo; atol=1e-5, rtol=0) &&
    isapprox(pair[2], hi; atol=1e-5, rtol=0);

fixed(row) = row.maximum - row.minimum <= 1e-5; # width at most 1e-5 mM/h

# ----------------------------------------------------------------------------
# Test cases | One scored assertion per labeled case
# ----------------------------------------------------------------------------

checks = [
    # Hand-checkable ranges on the small networks -
    "toy optimized ranges split one shared total" => () ->
        ranges_close(
            supply_ranges(swap_model(), 6; strategy=:optimized),
            [1, 2],
            [(0, 3), (3, 6)],
        ),

    "toy equal ranges are fixed" => () ->
        ranges_close(
            supply_ranges(swap_model(), 6; strategy=:equal),
            [1, 2],
            [(3, 3), (3, 3)],
        ),

    "default strategy is equal" => () ->
        ranges_close(supply_ranges(swap_model(), 6), [1, 2], [(3, 3), (3, 3)]),

    "toy optimized near-optimal ranges" => () ->
        ranges_close(
            supply_ranges(swap_model(), 6; strategy=:optimized, fraction=0.5),
            [1, 2],
            [(0, 3), (1.5, 6)],
        ),

    "toy equal near-optimal ranges" => () ->
        ranges_close(
            supply_ranges(swap_model(), 6; fraction=0.5),
            [1, 2],
            [(0, 3), (1.5, 3)],
        ),

    "zero fraction gives feasibility ranges" => () ->
        ranges_close(
            supply_ranges(swap_model(), 6; strategy=:optimized, fraction=0),
            [1, 2],
            [(0, 3), (0, 6)],
        ),

    "ranges follow a changed column order" => () ->
        ranges_close(
            supply_ranges(swap_model(protein_first=true), 6; strategy=:optimized),
            [2, 3],
            [(0, 3), (3, 6)],
        ),

    "tighter protein bound is respected" => () ->
        ranges_close(
            supply_ranges(swap_model(ceiling=1), 6; strategy=:optimized),
            [1, 2],
            [(0, 1), (1, 2)],
        ),

    # Input preservation and repeated calls -
    "model is preserved" => () -> begin
        m = deepcopy(model);
        before = deepcopy(m);
        supply_ranges(m, 1.0; strategy=:optimized);
        supply_ranges(m, 1.0; fraction=0.99);

        isequal(m, before)
    end,

    "repeated calls agree" => () -> begin
        m = swap_model();
        initial = supply_ranges(m, 6; strategy=:optimized);
        supply_ranges(m, 0; fraction=0.5);

        initial == supply_ranges(m, 6; strategy=:optimized)
    end,

    # Invalid inputs and solver failures -
    "negative budget rejected" => () ->
        throws_type(() -> supply_ranges(model, -1.0), ArgumentError),

    "infinite budget rejected" => () ->
        throws_type(() -> supply_ranges(model, Inf; strategy=:optimized), ArgumentError),

    "unknown strategy rejected" => () ->
        throws_type(() -> supply_ranges(model, 1.0; strategy=:other), ArgumentError),

    "invalid fractions rejected" => () ->
        all(
            throws_type(() -> supply_ranges(swap_model(), 6; fraction=f), ArgumentError)
            for f in (1.5, -0.1, NaN, Inf)
        ),

    "nonoptimal status is handled before reading fluxes" => () -> begin
        m = swap_model();
        m.lower[m.protein] = 1;

        throws_type(
            () -> supply_ranges(m, 0; strategy=:optimized),
            ErrorException;
            contains=["0", "optimized", "INFEASIBLE"],
        )
    end,

    # Reference ranges for the deGFP model -
    "deGFP equal supplies are fixed at 1.0" => () -> begin
        rows = supply_ranges(model, 1.0);

        all(fixed, rows) && close_to(range_of(rows, "Cysteine"), 0.012007, 0.012007) &&
            close_to(range_of(rows, "Leucine"), 0.05, 0.05) &&
            close_to(range_of(rows, "Tryptophan"), 0.006004, 0.006004)
    end,

    "deGFP optimized fixed supplies at 1.0" => () -> begin
        rows = supply_ranges(model, 1.0; strategy=:optimized);

        close_to(range_of(rows, "Leucine"), 0.141421, 0.141421) &&
            close_to(range_of(rows, "Lysine"), 0.133977, 0.133977) &&
            close_to(range_of(rows, "Arginine"), 0.044659, 0.044659)
    end,

    "deGFP supplies that are never used at 1.0" => () -> begin
        rows = supply_ranges(model, 1.0; strategy=:optimized);

        all(
            close_to(range_of(rows, a), 0, 0)
            for a in ("Alanine", "Asparagine", "Aspartate", "Cysteine", "Glycine", "Serine")
        )
    end,

    "deGFP interchangeable supplies at 1.0" => () -> begin
        rows = supply_ranges(model, 1.0; strategy=:optimized);

        close_to(range_of(rows, "Glutamate"), 0, 0.054715) &&
            close_to(range_of(rows, "Glutamine"), 0.059545, 0.114261) &&
            close_to(range_of(rows, "Threonine"), 0, 0.054715) &&
            count(fixed, rows) == 17
    end,

    "deGFP ranges open at translation capacity" => () -> begin
        rows = supply_ranges(model, 3.0; strategy=:optimized);

        all(isapprox(r.minimum, 0; atol=1e-5) for r in rows) &&
            close_to(range_of(rows, "Glutamate"), 0, 3.0)
    end,

    "deGFP near-optimal ranges at 1.0" => () -> begin
        rows = supply_ranges(model, 1.0; strategy=:optimized, fraction=0.99);

        !any(fixed, rows) && close_to(range_of(rows, "Valine"), 0.041459, 0.126534)
    end,

    "deGFP equal ranges open at translation capacity" => () ->
        close_to(range_of(supply_ranges(model, 3.0), "Leucine"), 0.112226, 0.15),

    "different glucose condition changes ranges" => () ->
        close_to(
            range_of(
                supply_ranges(load_model(glucose_limit=0.5), 1.0; strategy=:optimized),
                "Glutamate",
            ),
            0,
            0.396840,
        ),

    # Consistency between an optimal solution and its ranges -
    "ranges contain the reported optimal allocations" => () ->
        all(
            begin
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

                result["termination_status"] == JuMP.MOI.OPTIMAL &&
                    all(
                        r.minimum - 1e-6 <= result["argmax"][r.column] <= r.maximum + 1e-6
                        for r in rows
                    )
            end
            for strategy in (:equal, :optimized)
        ),
];


# ----------------------------------------------------------------------------
# Run the 24 scored checks
# ----------------------------------------------------------------------------

@testset "PS3 Part 3 (24 tests)" begin
    for (label, check) in checks
        @testset "$label" begin
            @test check();
        end
    end
end;
