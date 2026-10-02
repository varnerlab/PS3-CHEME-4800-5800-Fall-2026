"""
    reaction_index(model, id::AbstractString) -> Int

Find the column of `model.S` named `id` in `model.reactions`.
Throw `ArgumentError` when the identifier is absent. The model is unchanged.
"""
function reaction_index(model, id::AbstractString)
    j = findfirst(==(id), model.reactions);
    isnothing(j) && throw(ArgumentError("unknown reaction: $id"));
    return j;
end

"""
    load_model(; glucose_limit=1.0) -> NamedTuple

Load the supplied deGFP model with the PS3 operating conditions. `glucose_limit`
is a finite, nonnegative glucose supply limit in mM/h. Invalid limits raise
`ArgumentError`. Files are located relative to this source file.

The returned fields are `S` (species by reactions), `reactions`, `metabolites`,
`lower`, `upper`, `amino_acids`, `amino_acid_uptake` (column indices),
`protein` (protein-output column), `translation` (translation-initiation column),
`transcription_rate`, and `translation_limit`. All fluxes are in mM/h.
Separate supply and removal reactions have nonnegative fluxes.
"""
function load_model(; glucose_limit::Real=1.0)
    isfinite(glucose_limit) && glucose_limit >= 0 ||
        throw(ArgumentError("glucose_limit must be finite and nonnegative"));
    directory = joinpath(@__DIR__, "..", "data");
    S = readdlm(joinpath(directory, "stoichiometry.tsv"), '\t', Float64);
    rows = [split(line, '\t') for line in readlines(joinpath(directory, "reactions.tsv"))[2:end]];
    reactions = String[row[1] for row in rows];
    lower = [parse(Float64, row[2]) for row in rows];
    upper = [parse(Float64, row[3]) for row in rows];
    metabolites = readlines(joinpath(directory, "metabolites.txt"));
    size(S) == (length(metabolites), length(reactions)) || error("model dimensions do not agree");
    aa_rows = [split(line, '\t') for line in readlines(joinpath(directory, "amino_acids.tsv"))[2:end]];
    amino_acids = String[row[1] for row in aa_rows];
    lookup = Dict(id => j for (j, id) in enumerate(reactions));
    amino_acid_uptake = [lookup[row[2]] for row in aa_rows];

    # Apply the original Case 1 restrictions, with a smaller glucose supply -
    # Reaction names in the conditions table replace legacy numeric indices.
    for row in [split(line, '\t') for line in readlines(joinpath(directory, "conditions.tsv"))[2:end]]
        j = lookup[row[1]];
        lower[j] = parse(Float64, row[2]);
        upper[j] = parse(Float64, row[3]);
    end
    upper[lookup["M_glc_D_c_exchange"]] = Float64(glucose_limit);

    # Calculate fixed expression bounds from the original deGFP parameters -
    # These constants are evaluated before optimization; the optimization is linear.
    transcript_length = 683.0; # kinetic length parameter retained from the source
    transcription_rate = (25.0 / transcript_length) * (75.0 / 1e6) *
        (5.0 / (3.5 + 5.0)) * 3600.0 * (10.0 / 11.0); # mM/h
    mrna = transcription_rate / 5.2; # mM, steady-state transcript concentration
    translation_limit = (10.0 * 3.0 * 2.0 / transcript_length) * 3600.0 *
        0.0016 * mrna / (0.045 + mrna); # mM/h
    for id in ("transcriptional_initiation_deGFP", "mRNA_degradation_deGFP")
        lower[lookup[id]] = transcription_rate;
        upper[lookup[id]] = transcription_rate;
    end
    translation = lookup["translation_initiation_deGFP"];
    upper[translation] = translation_limit;
    return (; S, reactions, metabolites, lower, upper, amino_acids,
        amino_acid_uptake, protein=lookup["PROTEIN_export_deGFP"], translation,
        transcription_rate, translation_limit);
end
