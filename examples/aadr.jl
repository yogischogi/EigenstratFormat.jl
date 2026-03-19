# Example: Extract data from the AADR database.
#
# This example shows how to extract data of modern day individuals
# from the AADR database. It should be easy to adjust it to your
# own needs.
#
# Information about the database:
# Mallick, S., Micco, A., Mah, M. et al.
# The Allen Ancient DNA Resource (AADR) a curated compendium of ancient human genomes.
# Sci Data 11, 182 (2024). https://doi.org/10.1038/s41597-024-03031-7
#
# The database is very large. You need to download it from:
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
# Files used in this example:
# v62.0_HO_public.anno
# v62.0_HO_public.geno
# v62.0_HO_public.ind
# v62.0_HO_public.snp

using CSV, DataFrames, EigenstratFormat

# Input files.
# ADJUST basedir TO THE PATH ON YOUR COMPUTER.
const basedir = normpath("/home/dirk/Geno/AADR/database/")
const annofile = joinpath(basedir, "v62.0_HO_public.anno")
const genofile = joinpath(basedir, "v62.0_HO_public.geno")
const indfile = joinpath(basedir, "v62.0_HO_public.ind")
const snpfile = joinpath(basedir, "v62.0_HO_public.snp")

# Output files that will contain modern individuals.
const indfileout = joinpath(basedir, "modern.ind")
const snpfileout = joinpath(basedir, "modern.snp")
const annofileout = joinpath(basedir, "modern.anno")
const genofileout = joinpath(basedir, "modern.geno")

# Let's define some helper functions.

# Extract all modern inndividuals from .ind file.
function extract_modern_ind(indices::Vector{<:Integer}, infile, outfile)
    inds = read_eigenstrat_ind(infile)
    result = inds[indices, :]
    write_eigenstrat_ind(outfile, result)
end

# Extract all modern inndividuals from .anno file.
function extract_modern_anno(indices::Vector{<:Integer}, infile, outfile)
    annotations = read_eigenstrat_anno(infile)
    result = annotations[indices, :]
    write_eigenstrat_ind(outfile, result)
end

# Extract genotypes of modern individuals.
function extract_modern_geno(indices::Vector{Int64}, nsnp, nind, genofile, indfile, snpfile, outfile)
    geno = read_eigenstrat_geno(genofile, nsnp, nind; ind_idx = indices)
    indhash = hash_ids(indfile)
    snphash = hash_ids(snpfile)
    write_eigenstrat_geno(outfile, geno; ind_hash = indhash, snp_hash = snphash)
end

# Create a new database containing modern individuals.

# Determine which samples are modern by checking the age.
annotations = read_eigenstrat_anno(annofile)
# Create a vector that contains true or false for each entry which's age == 0.
modern = annotations[!, :("Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]")] .== 0

# Determine which samples are valid (not marked as "Ignore")
individuals = read_eigenstrat_ind(indfile)
valid = .~startswith.(individuals[!, :Status], "Ignore")
outlier = contains.(individuals[!, :Status], "outlier")

# Create a vector of indices we can use to access the database.
all_idxs = [i for i = 1:nrow(individuals)]
hits = modern .& valid .& (.~outlier)
# Filter index vector for all true entries in the hits vector.
idxs = filter(i -> hits[i], all_idxs)

# Use indices to create now .ind and .anno files.
extract_modern_ind(idxs, indfile, indfileout)
extract_modern_anno(idxs, annofile, annofileout)

# SNP file remains untouched.
# Copy it for consistency and overwrite an old version.
cp(snpfile, snpfileout; force = true)

nsnp = countlines(snpfile)
nind = countlines(indfile)
extract_modern_geno(idxs, nsnp, nind, genofile, indfile, snpfile, genofileout)

# That's it! You should now have 4 new files:
# modern.ind
# modern.snp
# modern.anno
# modern.geno


