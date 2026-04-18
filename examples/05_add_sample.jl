# Example: Add a new sample/individual to an existing database.
# The AADR database must exist on your computer.
# Requires autosomal results from Family Tree DNA or MyHeritage.
# Other vendors should be possible but not tested with this example.

using  EigenstratFormat

# ADJUST basedir TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/")

# ADJUST vendorfile to your autosomal results!
vendorfile = "my-heritage.csv"

# Prefixes for database names.
database_in = "v62.0_1240k_public"
database_out = "zzz_database"

# Full AADR database.
infileprefix = joinpath(basedir, database_in)
outfileprefix = joinpath(basedir, database_out)

add_individual(infileprefix, outfileprefix, vendorfile, "Me"; status = "Me")

