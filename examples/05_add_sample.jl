# Example: Add a new sample/individual to an existing database.
# The AADR database must exist on your computer.
# Requires autosomal results from Family Tree DNA or MyHeritage.
# Other vendors should be possible but not tested with this example.

using EigenstratFormat

# ADJUST basedir TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/v66.p1/")

# ADJUST vendorfile to your autosomal results!
vendorfile = "family-finder.csv"
#vendorfile = "my-heritage.csv"

# Prefixes for database names.
database_in = "v66.p1_1240K.aadr.patch.PUB"
database_out = "zzz_database"

# Full AADR database.
infileprefix = joinpath(basedir, database_in)
outfileprefix = joinpath(basedir, database_out)

add_individual(infileprefix, outfileprefix, vendorfile, "Me"; status = "Me")
