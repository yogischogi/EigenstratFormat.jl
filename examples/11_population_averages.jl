# Example: Calculate mean averages for populations and add them to the database.
using EigenstratFormat

# Database that was created in the 02_aadr.jl example.
# ADJUST basedir TO THE PATH ON YOUR COMPUTER.
basedir = normpath("/home/dirk/Geno/AADR/database/v66.p1/")

indfile = joinpath(basedir, "zzz_database.ind")
snpfile = joinpath(basedir, "zzz_database.snp")
genofile = joinpath(basedir, "zzz_database.geno")

indfileout = joinpath(basedir, "zzz_database_mean.ind")
snpfileout = joinpath(basedir, "zzz_database_mean.snp")
genofileout = joinpath(basedir, "zzz_database_mean.geno")

# Load database.
individuals = read_eigenstrat_ind(indfile)

# Get populations.
populations = population_idxs(individuals.Status)
min_size = 3
for (name, idxs) in populations
    if length(idxs) <= min_size
        delete!(populations, name)
    end
end

# Compute average genotype for each population.
genotypes = read_eigenstrat_geno(genofile)
i = 1
for (name, idxs) in populations
    println("Population: $name")
    avg = mean_genotype(genotypes[:, idxs])
    println("coverage: $(coverage(avg))")
    println()
    global genotypes = hcat(genotypes, avg)
    # Name for population, includes number of individuals.
    avg_name = "m$(i)_$(length(idxs))" 
    push!(individuals, [avg_name, "U", name])
    global i += 1
end

# Write new database containing average genotypes.
# SNPs remain untouched.
snps = read_eigenstrat_snp(snpfile)
write_eigenstrat_snp(snpfileout, snps)
write_eigenstrat_ind(indfileout, individuals)
snp_hash = hash_ids(snpfileout)
ind_hash = hash_ids(indfileout)
write_eigenstrat_geno(genofileout, genotypes; ind_hash = ind_hash, snp_hash = snp_hash)



