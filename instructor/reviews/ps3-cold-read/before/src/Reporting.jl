module ProductionReports

using Plots  # production-versus-budget figure and portable PNG/SVG export
using Printf # numerical tables with explicit units

"""
    write_outputs(assignment, directory) -> Vector{NamedTuple}

Run the supplied PS3 experiment with the functions in module `assignment`.
Save the curve, rates, and uptake allocations in `directory`, creating it if
needed. Independently check the allocation solutions at U=1 mM/h. Return the
production-curve rows. A nonoptimal or invalid allocation raises an error.
"""
function write_outputs(assignment, directory)
    mkpath(directory);
    model = assignment.load_model();
    budgets = collect(0.0:0.1:4.0);
    rows = assignment.production_curve(model, budgets);
    length(rows) == length(budgets) || error("production_curve must return one row per budget");
    all(isapprox(row.budget, b; atol=1e-10) for (row, b) in zip(rows, budgets)) ||
        error("production_curve changed the budget order");
    all(isfinite(row.equal) && isfinite(row.optimized) for row in rows) ||
        error("production_curve returned a nonfinite rate");
    open(joinpath(directory, "production-versus-budget.csv"), "w") do io
        println(io, "budget_mM_per_h,equal_uM_per_h,optimized_uM_per_h");
        for row in rows
            @printf(io, "%.6f,%.9f,%.9f\n", row.budget, row.equal, row.optimized);
        end
    end

    # Draw both strategies with distinct colors and line styles -
    theme(:default);
    p = plot(budgets, [r.optimized for r in rows]; label="Optimized allocation",
        color="#D55E00", linewidth=3,
        xlabel="Total amino-acid supply budget (mM/h)",
        ylabel="deGFP production rate (μM/h)",
        title="Feeding a cell-free protein factory",
        size=(1000, 650), dpi=160, xlims=(0, 4), ylims=(0, 12.5),
        legend=:bottomright, legendfontsize=11, guidefontsize=12, tickfontsize=10,
        titlefontsize=15, gridalpha=0.18, framestyle=:box,
        left_margin=7Plots.mm, bottom_margin=7Plots.mm, top_margin=5Plots.mm);
    # Draw the dashed curve last so both colors remain visible on the shared plateau.
    plot!(p, budgets, [r.equal for r in rows]; label="Equal allocation",
        color="#0072B2", linewidth=3, linestyle=:dash);
    hline!(p, [1000model.translation_limit]; label="Translation capacity",
        color="#555555", linewidth=1.5, linestyle=:dot);
    savefig(p, joinpath(directory, "production-versus-budget.png"));
    savefig(p, joinpath(directory, "production-versus-budget.svg"));

    # Inspect actual uptake at one shared budget; allocations need not be unique -
    results = Dict();
    for strategy in (:equal, :optimized)
        a = assignment.allocation_constraints(model, 1.0; strategy=strategy);
        result = assignment.solve_lp(model.S, a.lower, a.upper, assignment.protein_objective(model); A=a.A, b=a.b);
        result.optimal || error("allocation at U=1, $strategy: $(result.status)");
        checks = assignment.check_solution(model, result, a);
        checks.valid || error("allocation at U=1, $strategy failed the feasibility check");
        results[strategy] = result;
        @printf("%s at U=1: %.6f μM/h; total amino-acid uptake %.6f mM/h; maximum balance residual %.3e\n",
            String(strategy), 1000result.flux[model.protein], sum(result.flux[model.amino_acid_uptake]), checks.maximum_residual);
    end
    open(joinpath(directory, "amino-acid-allocation.csv"), "w") do io
        println(io, "amino_acid,equal_uptake_mM_per_h,optimized_uptake_mM_per_h,equal_limit_mM_per_h");
        for (name, j) in zip(model.amino_acids, model.amino_acid_uptake)
            @printf(io, "%s,%.9f,%.9f,%.9f\n", name, results[:equal].flux[j], results[:optimized].flux[j], 1/20);
        end
    end
    println("Wrote plot and data tables to ", abspath(directory));
    return rows;
end

end
