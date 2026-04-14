"Methods for the Eigenstrat format that is commonly used in genetics."
module EigenstratFormat

using CSV, DelimitedFiles, DataFrames, MultivariateStats, Statistics

# File operations
export read_eigenstrat_geno, write_eigenstrat_geno
export read_eigenstrat_snp, write_eigenstrat_snp
export read_eigenstrat_ind, write_eigenstrat_ind
export read_eigenstrat_anno
export read_vendor_data, write_23andMe
export add_individual, hash_ids

# Computations
export distance, getmarkers, impute_missing, impute_missing!, pca!, pca_coordinates
export population_idxs, remove_invariant!

include("files.jl")
include("utils.jl")
include("pca.jl")

end # module EigenstratFormat
