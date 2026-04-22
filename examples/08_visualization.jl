# Example: Plotting population distances.

using CSV, DataFrames, CairoMakie, EigenstratFormat

# Genetic distances from one individual to all others.
# See example 07_population_distances.jl.
popselectionfile = "zzz_population_selection.csv"

# Distances from one individual to populations.
# See example 06_genetic_distances.jl.
popdistancefile = "zzz_population_distances.csv"

distances = DataFrame(CSV.File(popselectionfile))
popdistances = DataFrame(CSV.File(popdistancefile))


# Create a graph with genetic distances on the x-axis
# and populations on the y-axis.
f = Figure()
ax = Axis(f[1, 2],
    title = "Genetic distances for a single individual",
    subtitle = "Populations averages = X, individuals = O,  1000 - 1800 CE",
    xlabel = "Distance",
    yreversed = true,
    limits = (nothing, nothing, 0, 21)
)
hideydecorations!(ax, ticks=false)

# Display firt 20 populations.

# Display individuals for each population.
pop_idxs = population_idxs(Vector(distances.country))
for i in 1:20
    population = popdistances.population[i]
    idxs = pop_idxs[population]
    xs = distances.distance[idxs]
    ys = [i for j in 1:length(xs)]
    scatter!(ax, xs, ys; markersize = 20)
end

# Average values
for i in 1:20
    scatter!(ax, popdistances.distance[i], i; marker = 'x', color = :black, markersize = 20)
end

# Population labels
poplabels = Axis(f[1, 1];
    limits = (0, nothing, 0, 21),
    yreversed = true
)
hidespines!(poplabels)
hideydecorations!(poplabels, grid = false)
hidexdecorations!(poplabels, grid = false)
#hideydecorations!(poplabels, ticks=false)
#hidexdecorations!(poplabels, ticks=false)
for i = 1:20
    text!(poplabels, 0, i + 0.5; text = popdistances.population[i])
end

# Adjust column sizes and show on screen.
colsize!(f.layout, 1, Aspect(1, 0.3))
f





