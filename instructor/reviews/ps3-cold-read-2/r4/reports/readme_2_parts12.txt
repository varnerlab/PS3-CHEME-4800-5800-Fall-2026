README.md:B019 | FORWARD | 1 | "let $p$ be the protein-output column, `model.protein`" | The name `model` is used for the first time here. Earlier text only spoke of "its returned named tuple" from load_model. | `model` is the named tuple that `load_model` returns.
README.md:B021 | REREAD | 2 | "Let $U\geq0$ be the total supply budget in mM/h" | Capital U (scalar budget) sits next to bold u (upper-bound vector) and u_j (its entries). In B023 both appear in min(u_j, U/K), so I had to check the case of each letter. | U is the single budget number; u_j is the per-reaction original upper bound.
README.md:B021 | FORWARD | 1 | "the `budget` argument of your functions" | No function has been named yet, so "your functions" has no referent until B029/B035. | `allocation_constraints(model, budget; ...)`, and indirectly `production_curve`.
README.md:B021 | INCONSISTENT | 1 | "`model.amino_acid_uptake`" | The field name says "uptake" but every sentence says "supply", so the reader must keep two words for one thing. | The uptake columns are the supply reactions.
README.md:B021 | REREAD | 1 | "The matrix $\mathbf{A}$ and vector $\mathbf{b}$ hold any additional supply inequality." | "any ... inequality" (singular) leaves it open whether A has zero, one, or several rows. This only becomes clear in B024 and B027. | Zero rows for :equal, exactly one row for :optimized.
README.md:B023 | ORDER | 1 | "u_j^{\mathrm{equal}}=\min(u_j,U/K)" | In the real model u_j = 30 always exceeds U/K, so the min looks pointless. Its reason (test networks with different bounds) only arrives in B030. | The min protects against small original bounds in test networks.
README.md:B024 | PURPOSE | 1 | "The supply fluxes can be smaller than their limits." | This sentence states the obvious, and its point only shows in the next sentence. | It sets up the claim that unused share is not redistributed.
README.md:B024 | REREAD | 1 | "Any unused share remains unavailable to the other amino acids." | "Remains unavailable" reads as if the share had been unavailable all along, which made me pause. | An amino acid's unused share is not transferred to the others.
README.md:B027 | ORDER | 1 | "Metabolic synthesis of amino acids remains available under both strategies" | A statement about both strategies is buried in the optimized-allocation paragraph. It also introduces "external supply" without definition and leaves open why the budget limits production at all if synthesis is on. | Synthesis reactions are not in I and do not count toward U, under either strategy.
README.md:B029 | ORDER | 1 | "returns a `Vector{Float64}` with one coefficient per reaction" | Part 1 writes the objective only as v_p, never as c^T v. Why the function returns a coefficient vector only becomes clear in B035. | This is c with c_p = 1 and 0 elsewhere, so c^T v = v_p.
README.md:B030 | ORDER | 1 | "Use the model's supplied column indices." | The instruction is vague ("as opposed to what?") and comes before its reason, the test networks with different orders. | Do not hard-code 265/20 or specific indices; use model.protein and model.amino_acid_uptake.
README.md:B030 | PURPOSE | 1 | "All test models have nonnegative supply fluxes and at least one supply column." | The text does not say why I need to know this or what it lets me skip. | K ≥ 1, so U/K is defined, and negative supply lower bounds need no handling.
README.md:B030 | REREAD | 1 | "Keep all lower bounds and unrelated upper bounds unchanged." | This repeats a B029 requirement in a different place, and "unrelated" is vague. | Only supply-column upper bounds change, and only under :equal.
README.md:B033 | ORDER | 1 | "For a hand check, consider two amino acids with supply fluxes" | The hand check comes after the "run the Part 1 checks" command. It would help more before implementing or testing. Its v_1, v_2 indices also clash with the j∈I column-index notation. | Use it to sanity-check the two formulas before coding.
README.md:B033 | INCONSISTENT | 1 | "The balances are $v_1=2v_p$ and $v_2=v_p$." | B027 says synthesis remains available, but this toy silently has no amino-acid synthesis, and it never says so. | The toy network has no synthesis route, so supply is the only source.
README.md:B037 | REREAD | 1 | "Check `result.optimal` before reading the fluxes." | production_curve needs the rate, not the flux vector. It is unclear whether to use result.objective or result.flux[model.protein]. | Either works, since the objective equals v_p in mM/h.
README.md:B037 | ORDER | 1 | "In that case the solver returns `nothing` for `flux` and `objective`." | The reason for the optimality check comes after the instruction to throw. | Check first, because non-optimal results carry no numbers.
README.md:B038 | REREAD | 1 | "the budget in mM/h and the maximum production rate under each strategy" | Three fields are matched to two descriptions. The ×1000 conversion from the solver's mM/h objective is not restated here; it appears only back in B017. | equal and optimized = 1000 × the LP objective, in µM/h.
README.md:B038 | REREAD | 1 | "Negative or nonfinite budgets raise `ArgumentError`." | It is unclear whether one bad entry in the vector fails the whole call, and whether to check before solving. It is also unclear whether letting allocation_constraints' error propagate is enough. | Any bad entry makes the call throw ArgumentError.
README.md:B042 | OVERLOAD | 2 | "the supply fluxes chosen by the driver's own solves, built with your Part 1 functions" | One bullet packs in file content, who solves, which functions, which budget, validation checks, and two terminal prints. "Driver's own solves" makes me wonder how these differ from production_curve. | production_curve returns only rates, so the driver re-solves at U = 1 with allocation_constraints and protein_objective to get the flux vectors.
README.md:B042 | INCONSISTENT | 1 | "the translation capacity, the highest production rate that translation allows" | B014 described this as the capacity of translation *initiation*, `translation_limit`, in mM/h. Here it is a protein production rate in µM/h, so the reader must infer a 1:1 link and the ×1000. | Plotted line = 1000 × translation_limit.
README.md:B042 | PURPOSE | 1 | "for both strategies at a total budget of **1 mM/h**" | No reason is given for choosing 1 mM/h over any other budget. | It is presumably the budget that questions 1–2 ask about.
README.md:B043 | FORWARD | 2 | "the driver skips the Part 3 files and says so" | "Part 3 files" have never been listed, and supply_ranges is named for the first time. | Output files built from supply_ranges (the FVA min/max supply ranges).
README.md:B043 | ORDER | 1 | "Open the plot and use it with the CSV tables to answer questions 1 and 2" | The driver-skip notice and the instruction to answer questions share one paragraph. The questions are not shown, so I do not know what to look for in the plot. | Go to responses.md for Q1–Q2 and use the plot plus both CSVs.
README.md:B044 | REFERENT | 1 | "The tests accept any solution that satisfies the constraints" | It is unclear which tests are meant, and which "solution": the allocation CSV is written by the driver, not returned by a tested function. It is also unclear whether the CSV itself is graded. | The tests check feasibility and the optimal rate, not exact flux values.

README.md:B019: Names the flux vector v and the protein-output index p = model.protein.
README.md:B020: States the production LP: maximize v_p subject to Sv = 0, bounds ℓ ≤ v ≤ u, and Av ≤ b.
README.md:B021: Says ℓ and u start from model bounds and A, b carry extra supply inequalities. Defines I (supply columns), K = |I| and budget U.
README.md:B022: Equal allocation caps each supply flux at the smaller of its original bound and an equal share.
README.md:B023: Gives the formula u_j_equal = min(u_j, U/K) for supply columns.
README.md:B024: Unused equal shares are not redistributed, and equal allocation returns an empty A (0×n) and an empty b.
README.md:B025: Optimized allocation keeps the original bounds and adds one shared constraint.
README.md:B026: The shared constraint is that the supply fluxes sum to at most U.
README.md:B027: A is a 1×n row of ones at supply columns and b = [U]. Synthesis is not counted as supply.
README.md:B028: Introduces the two Part 1 functions to implement.
README.md:B029: Specifies protein_objective (unit vector at protein) and allocation_constraints (return fields, strategies, copy semantics, ArgumentError cases).
README.md:B030: Functions must use the model's indices and work on small test networks, and must leave lower bounds and non-supply upper bounds unchanged.
README.md:B031: Introduces the Part 1 test command.
README.md:B033: A two-amino-acid toy example with budget 6 gives 1.5 (equal) vs 2 (optimized).
README.md:B035: production_curve must solve both strategies per budget using the Part 1 functions and solve_lp.
README.md:B037: Describes solve_lp result fields; a non-optimal result must raise an ErrorException with budget, strategy and status.
README.md:B038: Return format is one (budget, equal, optimized) tuple per input in µM/h, keeping order and repeats. Empty input gives an empty vector; bad budgets raise ArgumentError.
README.md:B039: Introduces the Part 2 test and driver commands.
README.md:B041: The driver sweeps budgets from 0 to 4 mM/h in steps of 0.1 and writes output files.
README.md:B042: Lists the outputs: plot with translation capacity, curve CSV, and the allocation CSV at 1 mM/h with feasibility checks and printed totals.
README.md:B043: The driver skips Part 3 outputs until supply_ranges exists. Use the plot and CSVs for Q1–Q2.
README.md:B044: Allocation fluxes can differ between students because optima are not unique. Tests accept any feasible optimal solution, and Part 3 quantifies the spread.

TOTALS: severity-3 = 0; severity-2 = 3; severity-1 = 22
