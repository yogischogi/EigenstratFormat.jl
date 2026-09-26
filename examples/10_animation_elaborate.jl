# Example: Create an animation from a table of samples with genetic distances.
# Use different colors for sample ages.
# Requires: Example 06_genetic_distances.jl
using CSV, DataFrames, EigenstratFormat
using CairoMakie, GeoMakie

# ADJUST basedir and annofile TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/v66.p1/")
annofile = joinpath(basedir, "v66.p1_1240K.aadr.PUB.anno")
distancesfile = "zzz_genetic_distances.csv"
outfile = "zzz_relatives.mp4"

annos = read_eigenstrat_anno(annofile)
distances = DataFrame(CSV.File(distancesfile))

# Add latitude and longitude.
id_col = 1
lat_col = 18
long_col = 19
sample_coordinates = DataFrame(id = annos[:, id_col], lat = annos[:, lat_col], long = annos[:, long_col])
# Make sure that lat and long exist.
coordinates = DataFrame(id = String[], lat = Float64[], long = Float64[])
for row in eachrow(sample_coordinates)
    try
        lat = parse(Float64, row.lat)
        long = parse(Float64, row.long)
        push!(coordinates, [row.id, lat, long])
    catch
        # Do nothing because coordinates are not always listed.
    end
end
# Add coordinates to distances.
distances = leftjoin(distances, coordinates, on = :id)
# Remove entries where lat or long are missing.
distances = subset(distances, :lat => x -> (!ismissing).(x) , :long => x -> (!ismissing).(x) )

# Filter for age of samples.
# We use only historical samples here.
# Modern samples often miss geographic coordinates.
distances = subset(distances, :age => x -> (x .> 0) .& (x .<= 5000))
sort!(distances, :distance)

# Animation parameters.
min = distances.distance[1]
max = distances.distance[nrow(distances)]
steps = 120
stepwidth = (max - min) / steps
iterator = min:stepwidth:max

# Create a figure that should be animated.
fig = Figure()
ga = GeoAxis(
    fig[1, 1];
    dest = "+proj=wintri", # the CRS in which you want to plot
    limits = (-20, 40, 35, 70), # Europe
    #limits = (-20, 120, 20, 70), # Europe and Asia
    #limits = (-180, 179, -60, 80), # Earth
)

# Add earth map.
img = rotr90(GeoMakie.earth())
meshimage!(ga, -180..180, -90..90, img; npoints = 500)
lines!(ga, GeoMakie.coastlines())

# Create markers for the plot and the legend.
markers = [:circle, :rect, :diamond, :hexagon, :cross]
colors = [:red, :orange, :green, :blue, :black]
labels = [
    "last 1000 years",
    "1000 - 2000 years",
    "2000 - 3000 years",
    "3000 - 4000 years",
    "4000 - 5000 years",
]
marker_elements = [
    MarkerElement(marker = markers[i], color = colors[i]) for i = 1:length(labels)
]
Legend(fig[1, 2], marker_elements, labels, "Sample age")

# Create an animation by modifying the figure parameters frame by frame.
record(fig, outfile, iterator; framerate = 2) do step
    println(step)
    ga.title = "Relatives, genetic distance <= $(trunc(step; digits = 1))"

    samples = subset(distances, :distance => d -> d .<= step)
    samples1000 = subset(samples, :age => a -> a .<= 1000) 
    samples2000 = subset(samples, :age => a -> (a .<= 2000) .& (a .> 1000))
    samples3000 = subset(samples, :age => a -> (a .<= 3000) .& (a .> 2000))
    samples4000 = subset(samples, :age => a -> (a .<= 4000) .& (a .> 3000))
    samples5000 = subset(samples, :age => a -> (a .<= 5000) .& (a .> 4000))

    scatter!(ga, samples1000[!, :long], samples1000[!, :lat]; marker = markers[1], color = colors[1])
    scatter!(ga, samples2000[!, :long], samples2000[!, :lat]; marker = markers[2], color = colors[2])
    scatter!(ga, samples3000[!, :long], samples3000[!, :lat]; marker = markers[3], color = colors[3])
    scatter!(ga, samples4000[!, :long], samples4000[!, :lat]; marker = markers[4], color = colors[4])
    scatter!(ga, samples5000[!, :long], samples5000[!, :lat]; marker = markers[5], color = colors[5])
end


