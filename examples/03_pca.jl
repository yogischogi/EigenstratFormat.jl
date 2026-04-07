# Example: Make a simple PCA (Principal Component Analysis).

using CairoMakie, CSV, EigenstratFormat, DataFrames, Statistics

# Database that was created in the 02_aadr.jl example.
# ADJUST basedir TO THE PATH ON YOUR COMPUTER.
basedir = normpath("/home/dirk/Geno/AADR/database/")

indfile = joinpath(basedir, "HGDP.ind")
snpfile = joinpath(basedir, "HGDP.snp")
annofile = joinpath(basedir, "HGDP.anno")
genofile = joinpath(basedir, "HGDP.geno")
coordinatesfile = joinpath(basedir, "HGDP_coordinates.csv")

# Load database.
individuals = read_eigenstrat_ind(indfile) 
genotypes = read_eigenstrat_geno(genofile)

# Remove invariant markers.
genotypes = remove_invariant!(genotypes)

# Impute missing allele values by taking the mean.
genotypes = impute_missing(genotypes)

# Calculate PCA model.
m = pca!(genotypes)

# Use m to get PCA coordinates for each sample.
coordinates = pca_coordinates(m, genotypes, individuals[:, :ID])

# Write coordinates to file.
# Not needed here but useful for the next example.
CSV.write(coordinatesfile, coordinates)

# Plot first two PCA coordinates using CairoMakie.
scatter(coordinates[:, :PC1], coordinates[:, :PC2])

