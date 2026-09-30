# Example: Find your closest relatives for different time periods.
# Requires: 06_genetic_distances.jl

using CSV, DataFrames, EigenstratFormat

# ADJUST basedir and annofile TO THE PATH ON YOUR COMPUTER!
basedir = normpath("/home/dirk/Geno/AADR/database/v66.p1/")
distancesfile = "zzz_genetic_distances.csv"


# Define parameters for a set of time periods.
# Here we divide the last 5000 years by intervals of 500 years.
start_year = -3000
interval_length = 500
final_year = 2000
intervals = (final_year - start_year) / interval_length

# Return the interval number according to the sample age.
interval(sample_age) = floor(Integer, (sample_age - start_year) / interval_length) + 1

# Read list of samples and distances.
samples = DataFrame(CSV.File(distancesfile))

# Create an array that holds a DataFrame of samples for each time period.
relatives_in_time = [similar(samples, 0) for _ in 1:intervals]

# Fill intervals with ancient samples.
for sample in eachrow(samples)
    year = 1950 - sample.age
    if !ismissing(year)
        i = interval(year)
        if i > 0 && i <= intervals
            push!(relatives_in_time[i], sample)
        end
    end
end

# Calculate genetic distances for each table and print them to screen.
for (i, table) in enumerate(relatives_in_time)
    sort!(table, :distance)
    start = start_year + (i - 1) * interval_length
    final = start + interval_length
    println()
    println("Closest matches from year: $start to: $final")
    println(first(table, 10))
end





