# Experimental: Add modal genotypes to database.
using EigenstratFormat

# Database that was created in the 02_aadr.jl example.
# ADJUST basedir TO THE PATH ON YOUR COMPUTER.
basedir = normpath("/home/dirk/Geno/AADR/database/v66.p1/")

indfile = joinpath(basedir, "zzz_database.ind")
snpfile = joinpath(basedir, "zzz_database.snp")
genofile = joinpath(basedir, "zzz_database.geno")

outindfile = joinpath(basedir, "zzz_modal_database.ind")
outsnpfile = joinpath(basedir, "zzz_modal_database.snp")
outgenofile = joinpath(basedir, "zzz_modal_database.geno")

# Load database.
individuals = read_eigenstrat_ind(indfile)

# Get populations.
populations = population_idxs(individuals.Status)
min_size = 20
for (name, idxs) in populations
    if length(idxs) <= min_size
        delete!(populations, name)
    end
end

# Compute modal genotype for each population.
genotypes = read_eigenstrat_geno(genofile)
i = 1
for (name, idxs) in populations
    println("Population: $name")
    modal = mode(genotypes[:, idxs])
    println("coverage: $(coverage(modal))")
    println()
    global genotypes = hcat(genotypes, modal)
    modal_name = "m$i" 
    push!(individuals, [modal_name, "U", name])
    global i += 1
end

# Write new database containing modal genotypes.
# SNPs remain untouched.
snps = read_eigenstrat_snp(snpfile)
write_eigenstrat_snp(outsnpfile, snps)
write_eigenstrat_ind(outindfile, individuals)
snp_hash = hash_ids(outsnpfile)
ind_hash = hash_ids(outindfile)
write_eigenstrat_geno(outgenofile, genotypes; ind_hash = ind_hash, snp_hash = snp_hash)










