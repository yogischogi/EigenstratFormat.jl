"""
    remove_invariant!(geno::Matrix{<:Real})

Remove invariant markers from a genomatrix.
Return a view to the cahanged geno matrix that contains only valid markers.

Because geno matrices can get very big the original matrix
is changed in place and is no longer valid.

`geno` is the matrix. Each row represents one single marker for multiple samples.
"""
function remove_invariant!(geno::Matrix{<:Real})
    # Remove invariant markers by copying non-invariant markers in place
    # into the old matrix.
    missing_value = 3
    nrow, ncol = size(geno)
    count = 0
    for i in 1:nrow
        values = filter(x -> x != missing_value, geno[i, :])
        for v in values
            if v != values[1]
                count += 1
                geno[count, :] = geno[i, :]
                break
            end
        end
    end
    return geno[1:count, :]
end

"""
    impute_missing(genomatrix::Matrix{<:Real})

Impute missing values by replacing them with mean values.

Return a Matrix{Float64}.

It is assumed that all SNPs are biallelic.
This function can really eat up lots of memory because the computation
is not done in place and the result is a Float64 matrix.
Maybe we should stay with UInt8 and use rounded values for missing alleles?

`geno` contains the number of reference (or derived) alleles for each sample and SNP.
Each row represents one SNP. Each sample is represented by a column.
"""
function impute_missing(geno::Matrix{<:Real})
    missing_value = 3
    nrow, ncol = size(geno)
    result = Matrix{Float64}(undef, nrow, ncol) 
    for i in 1:nrow
        # Compute mean for existing values.
        m = mean(filter(x -> x != missing_value, geno[i, :]))
        for j in 1:ncol
            if geno[i, j] == missing_value
                result[i, j] = m
            else
                result[i, j] = geno[i, j]
            end
        end
    end
    return result
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

`contains_reference` describes if the genomatrix contains the number
of reference alleles or derived alleles. Standard is `true` like documented
in David Reichs's [Input File Formats](https://reich.hms.harvard.edu/software/InputFileFormats).
This is also how Nick Patterson's
[smartpca](https://github.com/DReichLab/EIG/tree/master/src/eigensrc)
program interprets the file format.

However the R package
[smartsnp](https://github.com/ChristianHuber/smartsnp)
as described in 
[smartsnp, an r package for fast multivariate analyses of big genomic data](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13684)
assumes that the genomatrix contains the number of derived alleles.
To stay compatible with the R package choose `contains_reference = false`.
"""
function pca!(geno::Matrix{Float64}; contains_reference = true, ncomponents = 25, scaling = "genetic drift")
    if scaling != "z-score" && scaling != "genetic drift"
        throw("EigenstratFormat.pca: scaling must be z-score or genetic drift")
    end
    nrow, ncol = size(geno)
    adjusted::Matrix{Float64} = geno

    # If geno contains number of reference alleles
    # calculate number of derived alleles, assuming biallelic markers.
    # According to the Eigensoft manual the genomatrix contains the number
    # of reference alleles.
    # However, some researches seem to interpret this differently.
    # For compatibillity reasons we allow both options here.
    if contains_reference == true
        for i in 1:nrow, j in 1:ncol
            # Calculate number of derived alleles in place to save memory.
            adjusted[i, j] = 2 - adjusted[i, j]
        end
    end

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












