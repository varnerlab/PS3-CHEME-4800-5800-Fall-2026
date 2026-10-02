# PS3 written responses

Keep the five numbered questions and replace each `TODO` line with your answers
to every lettered part, giving the requested numbers with units. A few sentences
per part is enough.

1. **Why can optimized allocation match or improve protein production?**
   - (a) Compare the feasible sets of the two strategies at the same budget.
   - (b) Explain why equal supply limits do not require equal actual uptake.
   - (c) Explain why each strategy is a linear program, even though the
     equal-allocation bound min(u_j, U/K) takes a minimum.
   - (d) From `production-versus-budget.csv` at 1.0 mM/h, report both production
     rates and the percentage improvement of optimized over equal allocation.
   - (e) Sum each uptake column of `amino-acid-allocation.csv` and compare each
     total with the 1.0 mM/h budget.
   - (f) Name the amino acids whose `equal_uptake_mM_per_h` is more than
     10⁻⁵ mM/h below their limit of 0.05 mM/h. The equal-allocation column is the
     same in every optimal solution, so this list does not depend on the solver.

   TODO: Write your response.

2. **Why do the curves level off?**
   - (a) Using `production-versus-budget.csv` and the translation capacity printed
     by the driver, find the first sampled budget at which each strategy comes
     within 0.001 µM/h of the capacity.
   - (b) Explain why relaxing a supply constraint cannot reduce the optimal
     production rate.
   - (c) Explain why a larger budget eventually stops raising production. Which
     model bound would you relax next to raise the plateau, and which constraint
     might then limit production? Reason from the README's operating conditions;
     no new computation is needed.
   - (d) Does this steady-state calculation determine the final protein
     concentration after several hours? Explain.

   TODO: Write your response.

3. **What does the optimum determine?** Use the `optimized_min_mM_per_h` and
   `optimized_max_mM_per_h` columns of `supply-ranges.csv` (1.0 mM/h, γ = 1). A
   supply is fixed when its range width is at most 10⁻⁵ mM/h.
   - (a) List the amino acids never supplied externally (both endpoints within
     10⁻⁵ mM/h of zero) and those whose supply can vary.
   - (b) Name any two amino acids whose supply is fixed above the equal share
     U/K = 0.05 mM/h.
   - (c) The optimized uptake in `amino-acid-allocation.csv` is one optimal
     solution. Does each value lie within its range, and which values could
     another optimal solution change?
   - (d) Explain how internal synthesis allows different supply vectors to give
     the same production rate.

   TODO: Write your response.

4. **Are the ranges independent?** Parts (a) and (b) use the same columns as
   question 3.
   - (a) Sum the 20 values of `optimized_max_mM_per_h` and compare the total with
     the budget. Can every supply sit at its maximum at once?
   - (b) Subtract the sum of the fixed supplies from the budget to find how much
     the varying supplies can use at most. Now suppose glutamate, one of them, is
     supplied at its maximum. Compare what is left with the other varying
     supplies' minimums, and give the value, in mM/h, that each is forced to take.
   - (c) Sum the 20 values of `optimized_99pct_max_mM_per_h` (γ = 0.99) and
     compare the total with the budget.
   - (d) Explain why the 20 ranges, taken together, do not tell you which
     combinations of supply values reach the required production, cᵀv ≥ γz*.

   TODO: Write your response.

5. **How far can you trust one flux vector?**
   - (a) From `range-width-versus-budget.csv`, find the first sampled budget at
     which each strategy has no fixed supply flux. Compare these budgets with your
     answer to 2(a) and explain the relationship.
   - (b) From the `optimized_99pct` columns of `supply-ranges.csv` (optimized
     allocation, 1.0 mM/h, γ = 0.99), report how many supply fluxes remain fixed.
   - (c) Report valine's `optimized_99pct_min_mM_per_h` and how far it falls
     below valine's fixed optimized supply at γ = 1, in mM/h.
   - (d) FVA ranges are sometimes described as flux uncertainty. Explain what
     these ranges represent and what they do not.

   TODO: Write your response.
