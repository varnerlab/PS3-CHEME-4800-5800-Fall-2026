# ----------------------------------------------------------------------------
# Supplied reports | Write production and supply-range plots and tables.
# ----------------------------------------------------------------------------

"""Supplied plotting and table output for the PS3 allocation strategies and supply ranges."""
module ProductionReports

import JuMP  # solver termination status
using Plots  # production-versus-budget figure and portable PNG/SVG export
using Printf # numerical tables with explicit units

# ----------------------------------------------------------------------------
# Part 2 reports: production curves and allocations
# ----------------------------------------------------------------------------

"""
    write_outputs(assignment, directory) -> Vector{NamedTuple}

Generate the production curves and the amino-acid uptake table.

# Arguments

- `assignment`: the loaded `CellFreeProduction` module, including the student's
  completed objective, allocation, and production-curve functions.
- `directory`: output-directory path; create it if needed.

# Returns and files

Return the production-curve rows for budgets 0:0.1:4 mM/h. Write or overwrite
`production-versus-budget.csv`, `.png`, and `.svg`, plus
`amino-acid-allocation.csv` for U=1 mM/h. Protein rates are in μM/h; supply
rates are in mM/h. Print the translation capacity in μM/h and the U=1 production
rates, total uptake, and balance residuals. The uptake table records one optimal
flux vector per strategy; other vectors may achieve the same production rate.
Reset the plot theme to the default light appearance.

# Errors

Throw `ErrorException` for a curve with the wrong row count or budget order,
nonfinite production rates, or a nonoptimal or infeasible U=1 allocation.
Errors from student functions, file writing, or plotting propagate to the caller.
"""
function write_outputs(assignment, directory)

    # Compute the production rates and check the returned table layout -
    mkpath(directory);
    model = assignment.load_model();
    budgets = collect(0.0:0.1:4.0);
    rows = assignment.production_curve(model, budgets);

    # Check the returned curve before writing files -
    length(rows) == length(budgets) ||
        error("production_curve must return one row per budget");
    all(isapprox(row.budget, b; atol=1e-10) for (row, b) in zip(rows, budgets)) ||
        error("production_curve changed the budget order");
    all(isfinite(row.equal) && isfinite(row.optimized) for row in rows) ||
        error("production_curve returned a nonfinite rate");

    # Write the production table -
    open(joinpath(directory, "production-versus-budget.csv"), "w") do io
        println(io, "budget_mM_per_h,equal_uM_per_h,optimized_uM_per_h");
        for row in rows
            @printf(io, "%.6f,%.9f,%.9f\n", row.budget, row.equal, row.optimized);
        end
    end

    # Draw both strategies with distinct colors and line styles -
    theme(:default);
    p = plot(
        budgets,
        [r.optimized for r in rows];
        label="Optimized allocation",
        color="#D55E00",
        linewidth=3,
        xlabel="Total amino-acid supply budget (mM/h)",
        ylabel="deGFP production rate (μM/h)",
        title="Feeding a cell-free protein factory",
        size=(1000, 650),
        dpi=160,
        xlims=(0, 4),
        ylims=(0, 12.5),
        legend=:bottomright,
        legendfontsize=11,
        guidefontsize=12,
        tickfontsize=10,
        titlefontsize=15,
        gridalpha=0.18,
        framestyle=:box,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=5Plots.mm,
    );

    # Place the dashed curve above the solid curve so both colors show on the plateau.
    plot!(
        p,
        budgets,
        [r.equal for r in rows];
        label="Equal allocation",
        color="#0072B2",
        linewidth=3,
        linestyle=:dash,
    );
    hline!(
        p,
        [1000model.translation_limit];
        label="Translation capacity",
        color="#555555",
        linewidth=1.5,
        linestyle=:dot,
    );

    savefig(p, joinpath(directory, "production-versus-budget.png"));

    savefig(p, joinpath(directory, "production-versus-budget.svg"));
    @printf("Translation capacity: %.6f μM/h of deGFP\n", 1000model.translation_limit);

    # Inspect actual uptake at one shared budget; allocations need not be unique -
    results = Dict();

    for strategy in (:equal, :optimized)
        a = assignment.allocation_constraints(model, 1.0; strategy=strategy);

        # Use the supplied course type; hcat places each reaction's bounds in one row -
        calculation = assignment.build(
            assignment.MyPrimalFluxBalanceAnalysisCalculationModel,
            (
                S = model.S,
                fluxbounds = hcat(a.lower, a.upper),
                objective = assignment.protein_objective(model),
                species = model.metabolites,
                reactions = model.reactions,
            ),
        );

        result = assignment.solve(calculation; A=a.A, b=a.b);

        result["termination_status"] == JuMP.MOI.OPTIMAL ||
            error("allocation at U=1, $strategy: $(result["termination_status"])");

        checks = assignment.check_solution(model, result, a);
        checks.valid ||
            error("allocation at U=1, $strategy failed the feasibility check");

        results[strategy] = result;
        @printf(
            "%s at U=1: %.6f μM/h; total amino-acid uptake %.6f mM/h; maximum balance residual %.3e\n",
            String(strategy),
            1000result["argmax"][model.protein],
            sum(result["argmax"][model.amino_acid_uptake]),
            checks.maximum_residual,
        );
    end

    # Write one optimal allocation for each strategy -
    open(joinpath(directory, "amino-acid-allocation.csv"), "w") do io
        println(io, "amino_acid,equal_uptake_mM_per_h,optimized_uptake_mM_per_h,equal_limit_mM_per_h");
        for (name, j) in zip(model.amino_acids, model.amino_acid_uptake)
            @printf(
                io,
                "%s,%.9f,%.9f,%.9f\n",
                name,
                results[:equal]["argmax"][j],
                results[:optimized]["argmax"][j],
                1/20,
            );
        end
    end

    println("Wrote production plot and data tables to ", abspath(directory));
    return rows;
end


# ----------------------------------------------------------------------------
# Part 3 reports: supply ranges and total range widths
# ----------------------------------------------------------------------------

"""
    write_range_outputs(assignment, directory) -> Vector

Generate the Part 3 supply-range table and plots from the student's `supply_ranges`.

# Arguments

- `assignment`: the loaded `CellFreeProduction` module, including the student's
  completed objective, allocation, and supply-range functions.
- `directory`: output-directory path; create it if needed.

# Returns and files

Return the three U=1 mM/h range vectors: equal, optimized, and optimized with
fraction 0.99. Write or overwrite `supply-ranges.csv`, `.png`, and `.svg` for
those scenarios, and `range-width-versus-budget.csv`, `.png`, and `.svg` for
budgets 0:0.1:4 mM/h with fraction 1. Ranges and widths are in mM/h. Print the
number of fixed supply fluxes in each U=1 scenario. Reset the plot theme to the
default light appearance. Each interval contains the feasible values of one
supply flux under the production requirement; its endpoints need not occur
alongside the endpoints of the other supplies in a single solution.

# Errors

Throw `ErrorException` for ranges with the wrong columns, nonfinite values, or a
minimum above its maximum. Errors from student functions, file writing, or
plotting propagate to the caller.
"""
function write_range_outputs(assignment, directory)

    # Initialize the model and reporting budgets -
    mkpath(directory);
    model = assignment.load_model();
    budgets = collect(0.0:0.1:4.0);
    theme(:default);

    # Report the flux variability of every supply flux at U=1 mM/h -
    share = 1.0 / length(model.amino_acid_uptake);
    labels = ["equal", "optimized", "optimized_99pct"];
    ranges = [
        checked_ranges(assignment, model, 1.0, :equal, 1.0),
        checked_ranges(assignment, model, 1.0, :optimized, 1.0),
        checked_ranges(assignment, model, 1.0, :optimized, 0.99),
    ];

    open(joinpath(directory, "supply-ranges.csv"), "w") do io
        println(
            io,
            "amino_acid,",
            join(
                ["$(l)_$(e)_mM_per_h" for l in labels for e in ("min", "max")],
                ",",
            ),
        );
        for (k, name) in enumerate(model.amino_acids)
            print(io, name);
            for r in ranges
                @printf(io, ",%.9f,%.9f", tidy(r[k].minimum), tidy(r[k].maximum));
            end
            println(io);
        end
    end

    # Report how many supplies are fixed in each scenario -
    for (label, r) in zip(labels, ranges)
        @printf(
            "supply ranges at U=1, %s: %d of %d supply fluxes fixed (width <= %.0e mM/h)\n",
            label,
            count(isfixed, r),
            length(r),
            FIXED_WIDTH,
        );
    end

    # Draw each range as a horizontal interval; the 99% band sits behind the optimum -
    y = collect(1:length(model.amino_acids));
    q = plot(;
        xlabel="Supply flux at a 1 mM/h budget (mM/h)",
        ylabel="",
        title="Amino-acid supply ranges at a 1 mM/h budget",
        yticks=(y, model.amino_acids),
        yflip=true,
        ylims=(0.3, length(y) + 0.7),
        xlims=(0, 1.08 * maximum(r.maximum for r in ranges[3])),
        size=(1000, 820),
        dpi=160,
        legend=:topright,
        legendfontsize=10,
        guidefontsize=12,
        tickfontsize=10,
        titlefontsize=15,
        gridalpha=0.18,
        framestyle=:box,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=5Plots.mm,
    );
    plot!(
        q,
        interval_path(ranges[3], y)...;
        label="Optimized, 99% of maximum output",
        color="#D55E00",
        alpha=0.25,
        linewidth=10,
    );
    plot!(
        q,
        interval_path(ranges[2], y)...;
        label="Optimized, maximum output",
        color="#D55E00",
        linewidth=3,
        markershape=:circle,
        markersize=5,
        markerstrokewidth=0,
    );
    plot!(
        q,
        interval_path(ranges[1], y .- 0.3)...;
        label="Equal, maximum output",
        color="#0072B2",
        linewidth=3,
        markershape=:diamond,
        markersize=6,
        markerstrokewidth=0,
    );
    vline!(q, [share]; label="Equal share", color="#555555", linewidth=1.5, linestyle=:dot);

    savefig(q, joinpath(directory, "supply-ranges.png"));

    savefig(q, joinpath(directory, "supply-ranges.svg"));

    # Sum individual widths to describe variation; this is not a joint feasible volume -
    widths = Dict(s => Float64[] for s in (:equal, :optimized));
    counts = Dict(s => Int[] for s in (:equal, :optimized));

    for budget in budgets, strategy in (:equal, :optimized)
        r = checked_ranges(assignment, model, budget, strategy, 1.0);
        push!(widths[strategy], sum(x.maximum - x.minimum for x in r));
        push!(counts[strategy], count(isfixed, r));
    end

    # Write the total widths and fixed-supply counts -
    open(joinpath(directory, "range-width-versus-budget.csv"), "w") do io
        println(io, "budget_mM_per_h,equal_total_width_mM_per_h,optimized_total_width_mM_per_h,equal_fixed_supplies,optimized_fixed_supplies");
        for (k, budget) in enumerate(budgets)
            @printf(
                io,
                "%.6f,%.9f,%.9f,%d,%d\n",
                budget,
                tidy(widths[:equal][k]),
                tidy(widths[:optimized][k]),
                counts[:equal][k],
                counts[:optimized][k],
            );
        end
    end

    # Draw the total width across supply budgets -
    w = plot(
        budgets,
        widths[:optimized];
        label="Optimized allocation",
        color="#D55E00",
        linewidth=3,
        xlabel="Total amino-acid supply budget (mM/h)",
        ylabel="Sum of supply-range widths (mM/h)",
        title="Total supply-range width at maximum protein output",
        size=(1000, 650),
        dpi=160,
        xlims=(0, 4),
        ylims=(0, :auto),
        legend=:topleft,
        legendfontsize=11,
        guidefontsize=12,
        tickfontsize=10,
        titlefontsize=15,
        gridalpha=0.18,
        framestyle=:box,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
        top_margin=5Plots.mm,
    );
    plot!(
        w,
        budgets,
        widths[:equal];
        label="Equal allocation",
        color="#0072B2",
        linewidth=3,
        linestyle=:dash,
    );

    savefig(w, joinpath(directory, "range-width-versus-budget.png"));

    savefig(w, joinpath(directory, "range-width-versus-budget.svg"));
    println("Wrote supply-range plots and data tables to ", abspath(directory));
    return ranges;
end

# A supply flux is reported as fixed when its range is no wider than this, in mM/h.
const FIXED_WIDTH = 1e-5;

# ----------------------------------------------------------------------------
# Supporting functions: range checks and plotting coordinates
# ----------------------------------------------------------------------------

"""
    isfixed(row) -> Bool

Return whether a supply-range row's width, `row.maximum - row.minimum`, is at
most `FIXED_WIDTH` mM/h. Use a row already validated by `checked_ranges`;
this function does not check finiteness or endpoint order. The row is unchanged.
"""
isfixed(row) = row.maximum - row.minimum <= FIXED_WIDTH;

"""
    tidy(x) -> Float64

Return `x`, or zero when `abs(x) < 5e-10`, so that solver round-off does not
print as `-0.000000000` in the nine-decimal tables. The input is unchanged.
This is display rounding only; it does not relax the optimization constraints.
"""
tidy(x) = abs(x) < 5e-10 ? 0.0 : Float64(x);

"""
    interval_path(rows, y) -> Tuple{Vector{Float64},Vector{Float64}}

Build one plotting path that draws every supply range as a horizontal segment.
`rows` holds `(column, minimum, maximum)` named tuples and `y` the matching
vertical positions. Segments are separated by `NaN` so they do not connect.
The caller supplies equal-length collections; `zip` would otherwise stop at
the shorter one. Return new x and y vectors; the inputs are unchanged.
"""
function interval_path(rows, y)

    xs, ys = Float64[], Float64[];

    for (r, h) in zip(rows, y)
        append!(xs, (r.minimum, r.maximum, NaN));
        append!(ys, (h, h, NaN));
    end

    return xs, ys;
end


"""
    checked_ranges(assignment, model, budget, strategy, fraction) -> Vector

Call the student's `supply_ranges` for one scenario and check the returned layout.
Return the rows unchanged. Throw `ErrorException` unless there is one row per supply
column in `model.amino_acid_uptake` order, every value is finite, and no minimum
exceeds its maximum by more than 1e-7 mM/h. Errors from `supply_ranges` propagate.
This checks the output's layout and finite, ordered endpoints; it does not
independently establish that the student's endpoints are feasible or optimal.
"""
function checked_ranges(assignment, model, budget, strategy, fraction)

    r = assignment.supply_ranges(model, budget; strategy=strategy, fraction=fraction);
    scenario = "budget=$budget, strategy=$strategy, fraction=$fraction";

    # Check row count and reaction-column order -
    length(r) == length(model.amino_acid_uptake) &&
        [x.column for x in r] == model.amino_acid_uptake ||
        error("supply_ranges must return one row per supply column in model order ($scenario)");

    # Check each pair of endpoints -
    all(isfinite(x.minimum) && isfinite(x.maximum) && x.minimum <= x.maximum + 1e-7 for x in r) ||
        error("supply_ranges returned a nonfinite or reversed range ($scenario)");

    return r;
end

end
