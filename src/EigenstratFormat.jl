"Methods for the Eigenstrat format that is commonly used in genetics."
module EigenstratFormat

using CSV, DelimitedFiles, DataFrames

export read_eigenstrat_geno, write_eigenstrat_geno
export read_eigenstrat_snp, write_eigenstrat_snp
export read_eigenstrat_ind, write_eigenstrat_ind
export read_snp_file, write_23andMe
export add_individual, hash_ids

include("files.jl")

end # module EigenstratFormat
