# Example: Nice PCA plot using Makie.jl.

using CairoMakie, CSV, DataFrames, EigenstratFormat

# ADJUST basedir TO THE PATH ON YOUR COMPUTER.
basedir = normpath("/home/dirk/Geno/AADR/database/")

indfile = joinpath(basedir, "HGDP.ind")

# This file was created in example 03_pca.jl
coordinatesfile = joinpath(basedir, "HGDP_coordinates.csv")


# Determine populations and associated individuals.
individuals = read_eigenstrat_ind(indfile)
pop_indices = population_idxs(individuals.Status)

# Remove small populations to get a nice plot.
for (name, indices) in pop_indices
    if length(indices) < 3
        delete!(pop_indices, name)
    end
end

# Load PCA coordinates we created in example 03_pca.jl.
coordinates = DataFrame(CSV.File(coordinatesfile))

# Plot populations.
f = Figure()
ax = Axis(
    f[1, 1],
    title = "PCA coordinates of modern humans",
    xlabel = "PC1",
    ylabel = "PC2",
)
m = getmarkers(length(pop_indices))
pop_names = String[]
for (i, pop_name) in enumerate(keys(pop_indices))
    push!(pop_names, pop_name)
    # Get PCA coordinates for population.
    pc1 = Float64[]
    pc2 = Float64[]
    for popi in pop_indices[pop_name]
        push!(pc1, coordinates[popi, :PC1])
        push!(pc2, coordinates[popi, :PC2])
    end
    # Plot coordinates.    
    scatter!(ax, pc1, pc2, label = pop_name, marker = m.markers[i], color = m.colors[i])
end
# Add legend to plot.
marker_elements = [
    MarkerElement(marker = m.markers[i], color = m.colors[i]) for i = 1:length(pop_indices)
]
names = [pop_names[i] for i = 1:length(pop_indices)]
Legend(f[1, 2], marker_elements, names)
f
# Now we have plot that is similar to plots that have been around
# at least since the publication of the 1,000 genomes project.
# http://massgenomics.org/2012/11/human-genetic-variation-1000-genomes.html
# https://www.nature.com/articles/nature11632
