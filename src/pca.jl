# Functions for PCA plots.

"""
    getmarkers(n::Integer; smile = false)

Create a set of n markers for a plot. If smile == true
the last marker is a smiley.

Return a named tuple (markers, colors) that contains an array
of markers and an array of colors for the plot.

The number of markers is limited to 17 (18 including the smiley)
and 9 different colors.
"""
function getmarkers(n::Integer; smile = false)
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
    pca!(geno::Matrix{<:Real}; contains_reference = true, ncomponents = 25, scaling = "genetic drift")

Perform PCA analysis on Eigenstrat data using SVD (Single value decomposition).

This method does some computations in place to save memory.
This changes the genomatrix.

You can also chose if missing values in the genomatrix should be imputed.
if `impute = true` the genomatrix is not changed.

Return a `MulitvariateStats.PCA` that can be used to compute the
PCA components of a genomatrix by calling `pca_components`.

`geno` is a genomatrix that was read by `read_eigenstrat_geno`.
Usually the matrix contains the number of reference alleles for each SNP
and individual. `geno` must be a valid matrix with no missing values
or invariant markers. If you are not sure call `remove_invariant!(geno)` and
`impute_missing(geno)` before this function.

`ncomponents` is the maximum number of Eigenvalues for the PCA. Default
is 25.

`scaling` scales the number of alleles acccording to different methods.
Standard is `genetic drift` like recommended by Nick Patterson in
[Population Structure and Eigenanalysis](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020190).
"""
function pca!(geno::Matrix{Float64}; ncomponents = 25, scaling = "genetic drift")
    if scaling != "z-score" && scaling != "genetic drift"
        throw("EigenstratFormat.pca: scaling must be z-score or genetic drift")
    end
    nrow, ncol = size(geno)
    adjusted::Matrix{Float64} = geno

    # Mean number of derived alleles and standard deviations.
    m = zeros(Float64, nrow)
    for i in 1:nrow
        m[i] = mean(adjusted[i, :])
    end

    if scaling == "z-score"
        # Standardize: Center and scale (z-score).
        for i in 1:nrow
            s = std(adjusted[i, :])
            adjusted[i, :] = (adjusted[i, :] .- m[i]) / s
        end
    else 
        # Center and adjust for genetic drift.
        p = m ./ 2  # Allele frequemcy (see Patterson 2006).
        for i in 1:nrow
            adjusted[i, :] = (adjusted[i, :] .- m[i]) / sqrt( p[i] * (1 - p[i]) )
        end
    end

    # Use MultivariateStats.jl for PCA.
    return fit(PCA, adjusted; maxoutdim = ncomponents)
end

"""
    pca_coordinates(M::MultivariateStats.PCA, geno::Matrix{<:Real}, sample_ids::Vector{<:AbstractString})

Return PCA coordinates for all samples of a given genomatrix and `PCA M`.

The result is a Dataframe which contains a sample ID followed by the PCA coordinates
in each row.

`ID PC1 PC2 ... `
"""
function pca_coordinates(M::MultivariateStats.PCA, geno::Matrix{<:Real}, sample_ids::Vector{<:AbstractString})
    coordinates = predict(M, geno)
    coordinates = coordinates'
    _, ncol = size(coordinates)
    data = hcat(sample_ids, coordinates)
    names = vcat("ID", ["PC$i" for i = 1:ncol])
    return DataFrame(data, names)
end












