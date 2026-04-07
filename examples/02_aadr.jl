# Example: Extract data from the AADR database.
#
# This example shows how to extract data of modern day individuals (HGDP)
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

# ADJUST basedir TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/")

# Filenames of the database.
annofile = joinpath(basedir, "v62.0_HO_public.anno")
genofile = joinpath(basedir, "v62.0_HO_public.geno")
indfile = joinpath(basedir, "v62.0_HO_public.ind")
snpfile = joinpath(basedir, "v62.0_HO_public.snp")

# Output files for individuals from the HGDP project.
indfileout = joinpath(basedir, "HGDP.ind")
snpfileout = joinpath(basedir, "HGDP.snp")
annofileout = joinpath(basedir, "HGDP.anno")
genofileout = joinpath(basedir, "HGDP.geno")


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
function extract_modern_geno(indices::Vector{Int64}, genofile, indfile, snpfile, outfile)
    geno = read_eigenstrat_geno(genofile; ind_idxs = indices)
    indhash = hash_ids(indfile)
    snphash = hash_ids(snpfile)
    write_eigenstrat_geno(outfile, geno; ind_hash = indhash, snp_hash = snphash)
end

# Determine which samples belong to the HGDP project.
individuals = read_eigenstrat_ind(indfile)
is_valid = startswith.(individuals[!, :ID], "HGDP")

# Create a vector of indices we can use to access the database.
all_idxs = [i for i = 1:nrow(individuals)]
# Filter index vector for all true entries in the hits vector.
idxs = filter(i -> is_valid[i], all_idxs)

# Shrink database by selecting a subset of indices.
# This is useful if the resulting database gets
# too large or computations take too much time.
idxs = [idxs[i] for i in 10:10:length(idxs)]

# Use indices to create now .ind and .anno files.
extract_modern_ind(idxs, indfile, indfileout)
extract_modern_anno(idxs, annofile, annofileout)

# SNP file remains untouched.
# Copy it for consistency and overwrite an old version.
cp(snpfile, snpfileout; force = true)

nsnp = countlines(snpfile)
nind = countlines(indfile)
extract_modern_geno(idxs, genofile, indfile, snpfile, genofileout)

# That's it! You should now have 4 new files:
# HGDP.ind
# HGDP.snp
# HGDP.anno
# HGDP.geno

