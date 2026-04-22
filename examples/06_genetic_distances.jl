# Example: Calculate genetic distances of one sample to all other samples.
# Requires previous execution of example 05_add_sample.jl.

using CSV, EigenstratFormat, DataFrames

# ADJUST basedir TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/")

# Filename for the results.
resultfile = "zzz_genetic_distances.csv"

# Database that was created in the previous example.
genofile = joinpath(basedir, "zzz_database.geno")
indfile = joinpath(basedir, "zzz_database.ind")

# Annotations file of the AADR database.
annofile = joinpath(basedir, "v62.0_HO_public.anno")

# Column indices of important columns.
anno_id_col = 1
anno_age_col = 10
anno_country_col = 16

# We consider only samples with a minimum coverage.
# The chosen 0.5 in this example is just arbitrary.
# You can adjust this to your needs.
min_coverage = 0.5

# Read minimal information about individuals.
ind = read_eigenstrat_ind(indfile)

# Read full information about individuals from the anno file
# and create an index for faster access.
anno = read_eigenstrat_anno(annofile)
anno_idx = Dict{String, Int64}()
for (i, id) in enumerate(anno[:, anno_id_col])
    anno_idx[id] = i
end

# Read genotypes.
geno = read_eigenstrat_geno(genofile)
rows, cols = size(geno)

# Define a nice DataFrame for our samples.
result = DataFrame(
    index = Int64[],
    distance = Float64[],
    coverage = Float64[],
    id = String[],
    age = Int64[], 
    country = String[],
    population = String[]
)

# Compare the last sample to all other individuals in the database.
for i in 1:cols
    age = 0
    anno_id = ""
    country = ""
    population = ind.Status[i]

    # Extract ID and age from annotations.
    ind_id = ind.ID[i]
    if ind_id in keys(anno_idx)
        idx = anno_idx[ind_id]
        age = anno[idx, anno_age_col]
        country = anno[idx, anno_country_col]
    end

    # Add samples to result.
    c = coverage(geno[:, i])
    # Filter samples.
    if c >= min_coverage && age >= 0
        d = distance(geno[:, cols], geno[:, i])
        push!(result, [i, d, c, ind_id, age, country, population])
    end
end
sort!(result, :distance)
CSV.write(resultfile, result)

# Now you should have a nice table that can easily be opened with any
# spreadsheet program.

