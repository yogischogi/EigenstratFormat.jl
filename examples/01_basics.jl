# Examples for basic IO operations using the EigenstratFormat package.
#
# Eigenstrat databases can get very large. Thus it is neccessary to
# install a database separately from this package.
#
# Some smaller datasets can be found at David Reich's laboratory:
# https://reich.hms.harvard.edu/datasets
#
# A much more complete database is the
# Allen Ancient DNA Resource (AADR), a curated compendium of ancient human genomes
# https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW
#
#
# For testing purposes we use genotypes from the Lower Rhine-Meuse paper
# which can be found at:
# https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/Lower%20Rhine-Meuse%20paper-20260212T134708Z-1-001.zip
#
# 1. Download the database. You should now have a file named
#    Lower Rhine-Meuse paper-20260212T134708Z-1-001.zip
# 2. Unzip the above file. You should now have a directory named
#    Lower Rhine-Meuse paper including the database files.
# 3. Adjust the following file paths.

# Genotypes
const genopath = normpath("Lower Rhine-Meuse paper/LowerRhine.geno")
# Information about individuals
const indpath = normpath("Lower Rhine-Meuse paper/LowerRhine.ind")
# SNP information
const snppath = normpath("Lower Rhine-Meuse paper/LowerRhine.snp")


using DataFrames, EigenstratFormat

"""
    read_eigenstrat_database(genofile, indfile, snpfile)

Read database files an print some information to the screen.
"""
function read_eigenstrat_database(genofile, indfile, snpfile)
    # Read database files.
    snps = read_eigenstrat_snp(snpfile)
    individuals = read_eigenstrat_ind(indfile) 
    genotypes = read_eigenstrat_geno(genofile, nrow(snps), nrow(individuals))
    
    # Show some content on the screen.
    println("SNPs")
    print(first(snps, 5))
    println()
    println("Individuals")
    print(first(individuals, 5))
    println()
    println("Genotypes $(typeof(genotypes))")
    print(first(genotypes, 5))
end

"""
    reduce_eigenstrat_database(genofile, indfile, snpfile)

Reduce the size of the database.

Eigenstrat databases can get really huge and may not fit into
the computer's memory. Sometimes it is useful to work on a smaller
subset of individuls that are easier to handle with respect to
memory and computing time.
"""
function reduce_eigenstrat_database(genofile, indfile, snpfile)
    # Read database files.
    # Select each fourth individual when reading the genofile.
    snps = read_eigenstrat_snp(snpfile)
    individuals = read_eigenstrat_ind(indfile)
    indices = 4:4:nrow(individuals)
    genotypes = read_eigenstrat_geno(genofile, nrow(snps), nrow(individuals); ind_idx = [i for i in indices])

    # Generate temporary filenames for the reduced database.
    # These files are deleted automatically when the Julia interpreter exits.
    tempsnp = tempname("."; suffix = ".snp")
    tempind = tempname("."; suffix = ".ind")
    tempgeno = tempname("."; suffix = ".geno")

    # Write selected individuals to new data base.

    # SNPs are still the same.
    write_eigenstrat_snp(tempsnp, snps)

    # Write only selected individuals that are defined by indices.
    few_individuals = individuals[indices, 1:end]
    write_eigenstrat_ind(tempind, few_individuals)

    # Write geno file.
    ihash = hash_ids(tempind)
    shash = hash_ids(tempsnp)
    write_eigenstrat_geno(tempgeno, genotypes; ind_hash = ihash, snp_hash = shash)
end


# Using our methods.
read_eigenstrat_database(genopath, indpath, snppath)
reduce_eigenstrat_database(genopath, indpath, snppath)













