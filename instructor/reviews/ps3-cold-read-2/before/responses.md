# PS3 written responses

Keep the five numbered questions and replace each placeholder with your answer.
A short paragraph per question, supported by the requested numbers, is sufficient.

1. **Why can optimized allocation match or improve protein production?** Compare
   the feasible sets of the two strategies at the same budget. Explain why equal
   supply limits do not require equal actual uptake and why the problem remains linear.
   At 1.0 mM/h, report both production rates and the percentage improvement over
   equal allocation. Use the allocation table to compare the total supply each
   strategy actually uses, and name the amino acids that leave part of their
   equal share unused.

   TODO: Write your response.

2. **Why do the curves level off?** Identify the first sampled budget where each
   strategy reaches the translation-capacity limit within 0.001 µM/h. Explain why
   relaxing a supply constraint cannot reduce the optimal production rate, why more
   amino acids eventually stop helping here, and which model bound you would
   investigate next. Does this steady-state calculation determine the final protein
   concentration after several hours? Explain.

   TODO: Write your response.

3. **What does the optimum determine?** Use the supply-range table at 1.0 mM/h.
   For optimized allocation at maximum output, name the amino acids that are never
   supplied externally, two amino acids whose supply is fixed above the equal share
   of U/K = 0.05 mM/h, and the amino acids whose supply can vary. Compare these ranges with the single
   solution in the allocation table. Explain how internal synthesis allows
   different supply vectors to give the same production rate.

   TODO: Write your response.

4. **Are the ranges independent?** Use the optimized columns of the supply-range
   table at 1.0 mM/h. Add the 20 maximum supplies and compare the total with the
   budget. Can every supply sit at its maximum at the same time? Use the fixed
   supplies to find how much budget remains for the amino acids whose supply can
   vary. If glutamate is supplied at its maximum, what does that imply for
   glutamine and threonine? Repeat the sum for the `optimized_99pct` columns
   (γ = 0.99). Explain why a set of FVA ranges is not the set of feasible supply
   vectors.

   TODO: Write your response.

5. **How far can you trust one flux vector?** Use the range-width table to identify
   the first sampled budget at which each strategy has no fixed supply flux, and
   relate these budgets to your answer to question 2. Using the `optimized_99pct`
   columns of the supply-range table (optimized allocation, 1.0 mM/h, γ = 0.99),
   report how many supply fluxes remain fixed. Report valine's minimum supply and
   its decrease, in mM/h, from its value at maximum output.
   FVA ranges are sometimes described as flux uncertainty. Explain what these
   ranges represent and what they do not.

   TODO: Write your response.
