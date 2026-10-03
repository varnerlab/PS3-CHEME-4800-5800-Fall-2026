# PS3 written responses

Answer questions 1 and 2 for the Standard track. Question 3 is for the optional
Advanced track (one Magic Point); leave its `TODO` if you are not attempting it.
Replace each `TODO` you answer with one short paragraph (about 3–5 sentences),
keeping the blank line after the question. Use your results from `outputs/` to
support your explanations, and include units with numerical values.

1. **Does sharing the amino-acid budget improve protein production?**
   At a total budget of 1.0 mM/h, report the production rate for each strategy
   from `production-versus-budget.csv`. Under equal allocation, each amino acid
   can use at most 0.05 mM/h; under optimized allocation, the 20 amino acids share
   the total budget. Does either strategy produce more protein at this budget?
   Explain your result in terms of how the two strategies limit amino-acid supply.

   TODO: Write your response.

2. **Why does protein production stop increasing as the budget grows?**
   Use `production-versus-budget.png` to estimate the production rate at which
   the curves level off. Compare this rate with the translation capacity shown
   on the plot. Explain why supplying more amino acids no longer increases
   production once this capacity is reached.

   TODO: Write your response.

3. **Advanced (optional): Can different amino-acid supplies give the same maximum production rate?**
   Use `supply-ranges.csv` at a budget of 1.0 mM/h. The
   `optimized_min_mM_per_h` and `optimized_max_mM_per_h` columns give the lowest
   and highest supply rates that still allow maximum protein production. Find
   one amino acid whose supply can vary (maximum minus minimum greater than
   10⁻⁵ mM/h), and report its range. What does this example tell you about
   whether the supply rates returned by a single optimization are the only
   way to achieve maximum production?

   TODO: Write your response.
