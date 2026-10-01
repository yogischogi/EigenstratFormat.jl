# Utility functions to work with Eigenstrat databases.

"""
    remove_invariant!(geno::Matrix{<:Real})

Remove invariant markers from a genomatrix.
Return a view to the changed geno matrix that contains only valid markers.

Because geno matrices can get very big the original matrix
is changed in place and is no longer valid.

`geno` is the matrix. Each row represents one single marker for multiple samples.
"""
function remove_invariant!(geno::Matrix{<:Real})
    # Remove invariant markers by copying non-invariant markers in place
    # into the old matrix.
    nrow, ncol = size(geno)
    count = 0
    for i = 1:nrow
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
    impute_missing(geno::Matrix{UInt8}; ind_idxs = Int64[])

Copy the geno matrix and call `impute_missing!(...)`.
"""
function impute_missing(geno::Matrix{UInt8}; ind_idxs = Int64[])
    return impute_missing!(copy(geno); ind_idxs = ind_idxs)
end

"""
    impute_missing!(geno::Matrix{UInt8}; ind_idxs = Int64[])

Impute missing values by replacing them with mean values.

`geno` contains the number of reference (or derived) alleles for each sample and SNP.
Each row represents one SNP. Each sample is represented by a column.

`ind_idxs` specifies the indices of the individual samples that should
be used for the computation. This allows population-wise imputations.
if `ind_idx` is empty all samples are used.

`ind_idxs` must contain a minumum number of valid samples (currently 5) for
population-wise imputations. The this is not the case all samples are used
for the imputation.
"""
function impute_missing!(geno::Matrix{UInt8}; ind_idxs = Int64[])
    nrow, ncol = size(geno)

    # If indices are not specified use whole genomatrix.
    if length(ind_idxs) == 0
        ind_idxs = [i for i = 1:ncol]
    end

    # Compute mean values, 1 row represents 1 SNP.
    for i = 1:nrow
        # Calculate mean value using the specified values.
        m = _mean(geno[i, ind_idxs])
        if m == missing_value
            # Include all samples to compute the average.
            m = _mean(geno[i, :])
        end
        # If m is still missing set it to 0 (no mutations).
        if m == missing_value
            m = 0
        end

        # Fill up the geno matrix for the specified indices.
        for j in ind_idxs
            if geno[i, j] == missing_value
                geno[i, j] = m
            end
        end
    end
    return geno
end

"""
    population_idxs(population_names::Vector{<:String})

Return a Dictionary that contains a Vector of indices for each population.

Population name => [indices of individuals]

`population_names` represents a vector that contains a population
name for each sample in the .ind file and must be in the same order
to address individuals in the genomatrix properly.

Often the `Status` field of the .ind file contains population names.
"""
function population_idxs(population_names::Vector{<:Union{AbstractString,Missing}})
    # A Dictionary that contains population names and a list of sample indices.
    result = Dict{String,Vector{Integer}}()

    # Create a Dictionary with population names.
    for p in population_names
        if !ismissing(p)
            push!(result, p => String[])
        end
    end

    # Fill population entries with sample indices.
    for i = 1:length(population_names)
        if !ismissing(population_names[i])
            push!(result[population_names[i]], i)
        end
    end
    return result
end

"""
    distance(genotype1::Vector, genotype2::Vector; metric = "geometric")

Approximate the distance between two genotypes. The default metric
is the Euklidian/geometric metric.

Because the genotypes often contain missing values only valid
entries are used for the calculation and the rest is approximated.

It appears that because of the quadratic terms in the geometric distance
the Manhattan distance often yields better results when it comes to
approximation.

`metric` can be "geometric" or "manhattan".

To compare genotypes that were retrieved by haploid and pseudo-haploid calls
it is a good idea to make them compatible by calling pseudohaploid(genotype)
before calculating the distance.
"""
function distance(genotype1::Vector, genotype2::Vector; metric = "geometric")
    if metric == "geometric"
        return _geometric_distance(genotype1, genotype2)
    elseif metric == "manhattan"
        return _manhattan_distance(genotype1, genotype2)
    else
        throw("distance() only supports metrics 'geometric' and 'manhattan'.")
    end
end

"""
    _manhattan_distance(genotype1::Vector, genotype2::Vector)

Approximate the Manhattan distance between two genotypes.

Because the genotypes often contain missing values only valid
entries are used for the calculation and the rest is approximated.
"""
function _manhattan_distance(genotype1::Vector, genotype2::Vector)
    distance = 0
    comparisons = 0
    for i = 1:length(genotype1)
        # Make sure that a >= b.
        a = missing_value
        b = missing_value
        if genotype1[i] >= genotype2[i]
            a = genotype1[i]
            b = genotype2[i]
        else
            a = genotype2[i]
            b = genotype1[i]
        end
        # Calculate distance
        if a != missing_value && b != missing_value
            comparisons += 1
            distance += a - b
        end
    end
    approx = distance / comparisons * length(genotype1)
    return approx
end

"""
    _geometric_distance(genotype1::Vector, genotype2::Vector)

Approximate the geometric distance between two genotypes using
the Euklidian metric.

Because the genotypes often contain missing values only valid
entries are used for the calculation and the rest is approximated.
"""
function _geometric_distance(genotype1::Vector, genotype2::Vector)
    distance = 0
    comparisons = 0
    for i = 1:length(genotype1)
        if genotype1[i] != missing_value && genotype2[i] != missing_value
            comparisons += 1
            # This line works with unsigned Integers UInt8.
            distance += (genotype1[i] - genotype2[i])^2
        end
    end
    approx = distance / comparisons * length(genotype1)
    return sqrt(approx)
end

"""
    coverage(genotype::Vector)

Return the fraction of valid entries in the genotype.

Ancient DNA is often significantly degraded. So it is worth
to check the coverage. This may differ from the coverage value
in the metadata because the genotypes are often filtered before
use.
"""
function coverage(genotype::Vector)
    valids = 0
    for g in genotype
        if g != missing_value
            valids += 1
        end
    end
    return valids / length(genotype)
end

"""
    _mode(values::Vector{UInt8}, should_flip::Bool)

Calculate the mode/modal value of a list of values.
The missing_value is not used for the calculation.
"""
function _mode(values::Vector{UInt8})
    counts = zeros(Integer, 4)
    # Count values (0, 1, 2, 3).
    for v in values
        counts[v+1] += 1
    end
    # Only missing entries.
    if counts[4] == length(values)
        return missing_value
    end
    # Determine maximum.
    max = 0
    imax = 4  # 4 = index of missing value.
    for i = 1:3
        if counts[i] > max
            max = counts[i]
            imax = i
        elseif counts[i] == max
            # At least two values occur equally often.
            imax = 4
        end
    end
    return UInt8(imax - 1)
end

"""
    _mean(values::Vector{UInt8})

Calculate the mean value of a list of values.

This method returns only UInt8 integer values (0, 1, 2, 3).
The missing_value is not used for the calculation.
"""
function _mean(values::Vector{UInt8})
    count = 0
    total = 0
    for v in values
        if v != missing_value
            total += v
            count += 1
        end
    end
    if count == 0
        return missing_value
    end
    a = floor(total / count + 0.5)
    return UInt8(a)
end

"""
    mode_genotype(geno::Matrix{UInt8})

Return the modal genotpye.
"""
function mode_genotype(geno::Matrix{UInt8})
    nrow, _ = size(geno)
    result = zeros(UInt8, nrow)
    for i = 1:nrow
        result[i] = _mode(geno[i, :])
    end
    return result
end

"""
    mean_genotype(geno::Matrix{UInt8})

Return the mean genotpye.
"""
function mean_genotype(geno::Matrix{UInt8})
    nrow, _ = size(geno)
    result = zeros(UInt8, nrow)
    for i = 1:nrow
        result[i] = _mean(geno[i, :])
    end
    return result
end

"""
    pseudohaploid(genotype::Vector{UInt8})

Return a genotye that simulates pseudo-haploid calls
on a diploid genotype.

Note that this method introduces statistical noise. So the
results are not fully reproducible.

Gentotypes that were allready retrieved by pseude-haploid calls
remain the untouched.
"""
function pseudohaploid(genotype::Vector{UInt8})
    result = zeros(UInt8, length(genotype))
    for i = 1:length(genotype)
        if genotype[i] == 1
            # Introduce uncertainty to simulate pseudo haploid calling.
            # In the AADR databse samples which were called by pseudo-haploid
            # always produce 0 or 2.
            result[i] = rand([0, 2])
        else
            result[i] = genotype[i]
        end
    end
    return result
end
