# PS3 reference responses

Model answers to the three written questions. The numbers come from the
reference outputs; yours should agree to the digits shown.

1. **Does sharing the amino-acid budget improve protein production?**
   At a total budget of 1.0 mM/h, report the production rate for each strategy
   from `production-versus-budget.csv`. Under equal allocation, each amino acid
   can use at most 0.05 mM/h; under optimized allocation, the 20 amino acids share
   the total budget. Does either strategy produce more protein at this budget?
   Explain your result in terms of how the two strategies limit amino-acid supply.

   At a budget of 1.0 mM/h, equal allocation produces 6.00 µM/h of deGFP and
   optimized allocation produces 7.44 µM/h, about 24% more. deGFP does not use
   the 20 amino acids equally: it contains 19 leucines but only one tryptophan.
   Equal allocation still caps every amino acid at 0.05 mM/h, so arginine,
   cysteine, methionine, and tryptophan leave part of their share unused, while
   leucine and lysine are held below what the protein needs. Optimized
   allocation uses the whole 1.0 mM/h where it helps: for example, it supplies
   0.141 mM/h of leucine and none of six amino acids, such as alanine and
   glycine, that the mixture makes itself. Any supply that meets the equal caps
   also meets the shared budget, so optimized allocation can never produce less.

2. **Why does protein production stop increasing as the budget grows?**
   Use `production-versus-budget.png` to estimate the production rate at which
   the curves level off. Compare this rate with the translation capacity shown
   on the plot. Explain why supplying more amino acids no longer increases
   production once this capacity is reached.

   Both curves level off at about 11.18 µM/h, which is the translation capacity
   marked on the plot. Optimized allocation reaches it at a budget of 2.0 mM/h
   and equal allocation at about 2.9 mM/h. Each translation initiation makes one
   deGFP, and the model caps the initiation rate at 11.18 µM/h, so protein
   output cannot go higher no matter how many amino acids are supplied. Past
   that budget, extra amino acids go unused. Only a higher translation capacity,
   for example more ribosomes, would raise the plateau.

3. **Advanced (optional): Can different amino-acid supplies give the same maximum production rate?**
   Use `supply-ranges.csv` at a budget of 1.0 mM/h. The
   `optimized_min_mM_per_h` and `optimized_max_mM_per_h` columns give the lowest
   and highest supply rates that still allow maximum protein production. Find
   one amino acid whose supply can vary (maximum minus minimum greater than
   10⁻⁵ mM/h), and report its range. What does this example tell you about
   whether the supply rates returned by a single optimization are the only
   way to achieve maximum production?

   With optimized allocation at 1.0 mM/h, glutamate supply can be anywhere from
   0 to 0.0547 mM/h while production stays at its maximum of 7.44 µM/h.
   Threonine has the same range, and glutamine can vary from 0.0595 to
   0.1143 mM/h. These three supplies share one 0.0547 mM/h part of the budget:
   if glutamate takes all of it, threonine drops to 0 and glutamine to
   0.0595 mM/h, and the rest of the network adjusts to keep production the same.
   So a single optimization returns only one of many supply plans that reach the
   maximum rate, and another correct program can report a different glutamate
   supply.
