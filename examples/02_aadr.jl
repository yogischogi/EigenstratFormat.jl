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
# v66.p1_HO.aadr.PUB.anno
# v66.p1_HO.aadr.patch.PUB.geno
# v66.p1_HO.aadr.patch.PUB.ind
# v66.p1_HO.aadr.patch.PUB.snp

using CSV, DataFrames, EigenstratFormat

# ADJUST basedir TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/v66.p1/")

# Filenames of the database.
annofile = joinpath(basedir, "v66.p1_HO.aadr.PUB.anno")
genofile = joinpath(basedir, "v66.p1_HO.aadr.patch.PUB.geno")
indfile = joinpath(basedir, "v66.p1_HO.aadr.patch.PUB.ind")
snpfile = joinpath(basedir, "v66.p1_HO.aadr.patch.PUB.snp")

# Output files for individuals from the HGDP project.
indfileout = joinpath(basedir, "HGDP.ind")
snpfileout = joinpath(basedir, "HGDP.snp")
annofileout = joinpath(basedir, "HGDP.anno")
genofileout = joinpath(basedir, "HGDP.geno")

# Determine which samples belong to the HGDP project.
individuals = read_eigenstrat_ind(indfile)

# Create a vector of indices we can use to access the database.
all_idxs = [i for i = 1:nrow(individuals)]
# Filter index vector for all true entries in the hits vector.
is_valid = startswith.(individuals[!, :ID], "HGDP")
idxs = filter(i -> is_valid[i], all_idxs)

# Shrink database by selecting a subset of indices.
# This is useful if the resulting database gets
# too large or computations take too much time.
idxs = [idxs[i] for i = 10:10:length(idxs)]
write_eigenstrat_ind(indfileout, individuals[idxs, :])

# SNP file remains untouched.
# Copy it for consistency and overwrite an old version.
cp(snpfile, snpfileout; force = true)

# Annofile may also remain untouched.
# Samples in the annofile may not be in the same order
# as in the .ind file.
cp(annofile, annofileout; force = true)

# Select genotypes and write them to new .geno file.
geno = read_eigenstrat_geno(genofile; ind_idxs = idxs)
indhash = hash_ids(indfileout)
snphash = hash_ids(snpfileout)
write_eigenstrat_geno(genofileout, geno; ind_hash = indhash, snp_hash = snphash)

# That's it! You should now have 4 new files:
# HGDP.ind
# HGDP.snp
# HGDP.anno
# HGDP.geno
