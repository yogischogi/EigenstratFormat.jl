# Example: Calculating population distances.
# Requires a file containing genetic distances that was
# produced in example 06_genetic_distances.jl.

using CSV, EigenstratFormat, DataFrames, Statistics

# ADJUST basedir TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/")

# Output file for genetic distances from one individual to populations.
popdistancefile = "zzz_population_distances.csv"
popselectionfile = "zzz_population_selection.csv"

# Input file for genetic distances from one individual to all others.
# See example 06_genetic_distances.jl.
distancefile = "zzz_genetic_distances.csv"
distances = DataFrame(CSV.File(distancefile))

# Create populations by filtering according to individual age and origin.
selected = subset(distances, :age => x -> x .>= 150 .&& x .<= 1000)
populations = population_idxs(Vector(selected.country))

# Create a DataFrame including the population,
# the population distance to an individual and the
# number of samples in the population.
population_distances =
    DataFrame(population = String[], distance = Float64[], nsamples = Int64[])

# Calculate mean population distances and save them to file.
for (country, idxs) in populations
    d = mean(selected.distance[idxs])
    nsamples = length(idxs)
    push!(population_distances, [country, d, nsamples])
end
sort!(population_distances, :distance)
CSV.write(popdistancefile, population_distances)
CSV.write(popselectionfile, selected)
