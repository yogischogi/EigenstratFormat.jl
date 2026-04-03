# Example: Nice PCA plot using Makie.jl.

using CairoMakie, CSV, DataFrames, EigenstratFormat

# ADJUST basedir TO THE PATH ON YOUR COMPUTER.
const basedir = normpath("/home/dirk/Geno/AADR/database/")
const indfile = joinpath(basedir, "HGDP.ind")
# This file was created in example 03_pca.jl
const coordinatesfile = joinpath(basedir, "HGDP_coordinates.csv")

"""
    _getmarkers(n::Integer; smile = false)

Create a set of n markers for a plot. If smile == true
the last marker is a smiley.

Return a named tuple (markers, colors) that contains an array
of markers and an array of colors for the plot.

The number of markers is limited to 17 (18 including the smiley)
and 9 different colors.
"""
function _getmarkers(n::Integer; smile = false)
    markers = [:circle, :rect, :diamond, :hexagon, :cross,
        :xcross, :utriangle, :dtriangle, :ltriangle, :rtriangle,
        :pentagon, :star4, :star5, :star6, :star8,
        :vline, :hline]
    smiley = '😄'

    colors = [:grey0, :green1, :blue, :chocolate, :cyan3,
        :red, :fuchsia, :orange, :darkred]

    # Calculate markers and colors.
    plotmarkers = []
    plotcolors = []
    for i in 1:n
        j = i % length(markers) + 1
        push!(plotmarkers, markers[j])
        k = i % length(colors) + 1
        push!(plotcolors, colors[k])
    end

    # Replace last element with a smiley.
    if smile == true
        pop!(plotmarkers)
        push!(plotmarkers, smiley)        
    end
    return (markers = plotmarkers, colors = plotcolors)
end

"""
    population_indices(eigenstrat_inds::DataFrame)

Return a Dictionary that contains a Vector of indices for each population.

Population name => [indices of individuals]
"""
function population_indices(eigenstrat_inds::DataFrame)
    # A Dictionary that contains population names and a list of sample indices.
    result = Dict{String, Vector{Integer}}()

    # Create empy Dictionary of population names.
    for p in eigenstrat_inds.Status
        push!(result, p => String[])
    end

    # Fill population entries with sample indices.
    # The population name is in the :Status column.
    for i in 1:nrow(eigenstrat_inds)
        push!(result[eigenstrat_inds[i, :Status]], i)
    end
    return result
end

# Determine populations and associated individuals.
individuals = read_eigenstrat_ind(indfile)
pop_indices = population_indices(individuals)

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
ax = Axis(f[1, 1],
    title = "PCA coordinates of modern humans",
    xlabel = "PC1",
    ylabel = "PC2",
)
m = _getmarkers(length(pop_indices))
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
marker_elements = [MarkerElement(marker = m.markers[i], color = m.colors[i]) for i = 1:length(pop_indices)]
names = [pop_names[i] for i = 1:length(pop_indices)]
Legend(f[1, 2], marker_elements, names)
f
# Now we have plot that is similar to plots that have been around
# at least since the publication of the 1,000 genomes project.
# http://massgenomics.org/2012/11/human-genetic-variation-1000-genomes.html
# https://www.nature.com/articles/nature11632


