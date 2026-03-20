"""
    impute_missing(genomatrix::Matrix){UInt8}

Impute missing values by replacing them with mean values.

Return a Matrix{Float64}.

It is assumed that all SNPs are biallelic.

`genomatrix` contains the number of reference (or derived) alleles for each sample and SNP.
Each row represents one SNP. Each sample is represented by a column.
"""
function impute_missing(genomatrix::Matrix{UInt8})
    missing_value = 3
    nrow, ncol = size(genomatrix)
    result = Matrix{Float64}(undef, nrow, ncol)
    
    for i in 1:nrow
        # Compute mean for existing values.
        m = mean(filter(x -> x != missing_value, genomatrix[i, :]))
        for j in 1:ncol
            if genomatrix[i, j] == missing_value
                result[i, j] = m
            else
                result[i, j] = genomatrix[i, j]
            end
        end
    end
    return result
end

"""
    pca(genomatrix, indframe; ncomponents = 25, scaling = "genetic drift", geno_contains_reference = true)

Perform PCA analysis on Eigenstrat data using SVD (Single value decomposition).
Missing SNP values are imputed as mean values.

Return a DataFrame where the first column contains individual IDs and the following
the PCA coordinates, so that each individual is represented by a single row.

`genomatrix` is a genomatrix that was read by `read_eigenstrat_geno`.
Usually the matrix contains the number of reference alleles for each SNP
and individual.

`indframe` is a DataFrame that was read by 'read_eigenstrat_ind`.
It contains data for the individuals in the Eigenstrat database.

`ncomponents` is the maximum number of Eigenvalues for the PCA. Default
is 25.

`scaling` scales the number of alleles acccording to different methods.
Standard is `genetic drift` like recommended by Nick Patterson in
[Population Structure and Eigenanalysis](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020190).

`geno_contains_reference` describes if the genomatrix contains the number
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
To stay compatible with the R package choose `geno_contains_reference = false`.
"""
function pca(genomatrix::Matrix{UInt8}, indframe; ncomponents = 25, scaling = "genetic drift", geno_contains_reference = true)
    if scaling != "z-score" && scaling != "genetic drift"
        throw("EigenstratFormat.pca: scaling must be z-score or genetic drift")
    end

    nrow, ncol = size(genomatrix)
    adjusted = impute_missing(genomatrix)    

    # If genomatrix contains number of reference alleles
    # calculate number of derived alleles, assuming biallelic markers.
    # According to the Eigensoft manual the genomatrix contains the number
    # of reference alleles.
    # However, some researches seem to interpret this differently.
    # For compatibillity reasons we allow both options here.
    if geno_contains_reference == true
        for i in 1:nrow, j in 1:ncol
            # Calculate number of derived alleles in place to save memory.
            adjusted[i, j] = 2 - adjusted[i, j]
        end
    end

    # Remove invariant markers by copying non-invariant markers in place
    # into the old matrix.
    count = 0
    for i in 1:nrow
        first = adjusted[i, 1]    
        for j in 1:ncol
            if adjusted[i, j] != first
                count += 1
                adjusted[count, :] = adjusted[i, :]
                break
            end
        end
    end
    adjusted = adjusted[1:count, :]
    nrow, ncol = size(adjusted)

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
    M = fit(PCA, adjusted; maxoutdim = ncomponents)
    coordinates = predict(M, adjusted)
    co_cols = size(coordinates', 2)
    colnames = vcat("ID", ["C$i" for i = 1:co_cols])
    result_data = hcat(indframe[!, :ID], coordinates')
    return DataFrame(result_data, colnames)
end


