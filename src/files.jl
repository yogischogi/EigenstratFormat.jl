# Functions to read and write Eigenstrat data.

"""
    _hashit(sequence::String)

Calculate the hashsum of a String.

This is basically the same method as in Nick Patterson's original
[C code](https://github.com/DReichLab/EIG/blob/master/src/admutils.c).

However, the original C version uses 32 bit integer values and integer
overflows which may result in undefined behavior on some machines. 
This method uses 64 bit integers to circumvent this problem.

This method will probably fail if the sequence contains non-ASCII
Unicode characters. I do not know if this is defined in any way.
"""
function _hashit(sequence::AbstractString)
    hash::Int64 = 0
    for c in sequence
        hash *= 23
        hash = hash % (2^32) 
        hash += Int64(c)
    end
    return hash
end

"""
    _hashsum(sequences::Vector{<:AbstractString})

Calculate the hashsum of a Vector of Strings using the _hashit() method.

This is the same method as the `hasharr` function in
Nick Patterson's [C code](https://github.com/DReichLab/EIG/blob/master/src/admutils.c).
"""
function _hashsum(sequences::Vector{<:AbstractString})
    hash::Int64 = 0
    for s in sequences
        thash = _hashit(s)
        hash *= 17
        hash = hash % (2^32) 
        hash = xor(hash, thash)
    end
    return hash
end

"""
    hash_ids(filename::AbstractString)

Create the hashsum of .snp and .ind files.
SNP and Individual files contain an ID in the first column.
This method uses those IDs to calculate a hashsum.
The sums are needed to store genotype data in packed format.
"""
function hash_ids(filename::AbstractString)
    data = readdlm(filename, String)
    hash = _hashsum(data[:, 1])
    return hash
end

"""
    _bitpair(byte::UInt8, pos::Int64)

Extract 2 bits from a byte. The bit positions start at 0.
Allowed are values 0, 1, 2, 3.
"""
function _bitpair(byte::UInt8, pos::Integer)
    if pos == 0
        (byte & 0b11000000) >> 6
    elseif pos == 1
        (byte & 0b00110000) >> 4
    elseif pos == 2
        (byte & 0b00001100) >> 2
    elseif pos == 3
        (byte & 0b00000011)
    end
end

"""
    _alleles(bytes::Vector{UInt8}, idx::Int64)

Return the number of variant alleles for a specified individual.

`bytes`: Row of bytes that encode an SNP for all individuals.

`idx`: Index of the individual.
"""
function _alleles(bytes::Vector{UInt8}, idx::Int64)
    # Find byte that encodes the specified SNP.
    pos = Int64(ceil(idx / 4))
    byte = bytes[pos]

    # Extract bitpair that contains individual's data.
    # bitpair index starts at 0.
    bitpair_no = (idx - 1) % 4
    bits = _bitpair(byte, bitpair_no)
    return UInt8(bits)
end

"""
    read_eigenstrat_geno(
        genofile::AbstractString,
        nsnp::Int64,
        nind::Int64;
        ind_idx::Vector{Int64} = [i for i = 1:nind]
    )

Read a genofile in PackedAncestryMap format. The file must be unzipped.

`genofile`: filename

`nsnp`: number of SNPs listed in .snp file.

`nind`: number of individuals listed in .ind file.

`ind_idx`: Indices of individuals that should be read from the file.

XXX Check for comment lines in .snp and .ind files.

File description: File header starts with GENO or TGENO (transposed GENO).
So far files in the AADR archive seem to be GENO. So this method
does not support the transposed TGENO format.

The text format contains one line per genotype:

SNP_ID  Sample_ID   Number_of_variant_alleles

The packed format:

Each SNP entry has 2 bits: 0, 1, 2, 3=missing, that denote the number
of variant alleles as described at [David Reich's laboratory](https://reich.hms.harvard.edu/software/InputFileFormats).
"""
function read_eigenstrat_geno(
    genofile::AbstractString,
    nsnp::Int64,
    nind::Int64;
    ind_idx::Vector{Int64} = [i for i = 1:nind]
)
    result = zeros(UInt8, nsnp, length(ind_idx))

    # 1 SNP value for 4 individuals is encoded as 1 byte.
    bytes_per_line = Int64(ceil(nind / 4))
    # Header size must be at least 48 bytes.
    if bytes_per_line < 48
        bytes_per_line = 48
    end

    # Read file line by line and extract SNP alleles for all individuals.
    open(genofile) do io
        buffer = Array{UInt8}(undef, bytes_per_line)
        # Skip first line because it is a header.
        read!(io, buffer)
        for snp = 1:nsnp
            read!(io, buffer)
            # Extract alleles for all individuals.
            for i in 1:length(ind_idx)
                result[snp, i] = _alleles(buffer, ind_idx[i])
            end
        end
    end
    return result
end

"""
    write_eigenstrat_geno(
        genofile::AbstractString,
        genomatrix::Matrix{UInt8};
        ind_hash::Int64 = 0,
        snp_hash::Int64 = 0
    )

Write a genotype matrix to file in PackedAncestryMap format.

`ind_hash`: Hashsum of .ind file.

`snp_hash`: Hashsum of .snp file.
"""
function write_eigenstrat_geno(
    genofile::AbstractString,
    genomatrix::Matrix{UInt8};
    ind_hash::Int64 = 0,
    snp_hash::Int64 = 0
)
    (rows, cols) = size(genomatrix)

    # 4 columns are encoded as 1 byte.
    bytes_per_line = Int64(ceil(cols / 4))
    # Header size must be at least 48 bytes.
    minlen = 48
    fillbytes = zeros(UInt8, 0)
    if bytes_per_line < minlen
        fillbytes = zeros(UInt8, minlen - bytes_per_line)
        bytes_per_line = minlen
    end

    # Create header.
    ihash = string(ind_hash, base = 16)
    shash = string(snp_hash, base = 16)
    buffer = zeros(UInt8, bytes_per_line)
    io = IOBuffer(buffer, write = true)
    header = "GENO $cols $rows $ihash $shash"
    # This should never happen, but the file would still be
    # useable if the hash sums are not checked.
    if length(header) > bytes_per_line
        header = header[1: bytes_per_line]
    end
    write(io, header)

    # Write to file.
    open(genofile, create = true, write = true) do file
        write(file, buffer)
        for row in 1:rows
            byte::UInt8 = 0
            for col in 1:cols
                value = genomatrix[row, col]
                # Encode in byte.
                bitpair_pos = (col - 1) % 4
                shift = 6 - 2 * bitpair_pos
                byte = byte | (value << shift)

                # Write byte to file.
                if  bitpair_pos == 3 || col == cols
                    write(file, byte)
                    byte = 0
                end
            end
            write(file, fillbytes)
        end
        flush(file)    
    end
end

"""
    read_snp_file(filename::AbstractString)

Read file with autosomal results from FTDNA Family Finder, My Heritage
or LivingDNA. Should also work with 24andMe files but not tested.

Return a `DataFrame` containing the columns:

rsid  chromosome  position  genotype
"""
function read_snp_file(filename::AbstractString)
    snptable = CSV.read(filename, DataFrame; header = ["rsid", "chromosome", "position", "genotype"],
        types = [String, String, String, String], comment = "#")
    # Drop header line (MyHeritage, 23andMe).
    if snptable.rsid[1] == "RSID" || snptable.rsid[1] == "rsid"
        delete!(snptable, 1)
    end
    return snptable
end

"""
    read_eigenstrat_snp(snpfile::AbstractString)

Read Eigenstrat .snp file. The SNP file contains information about
each SNP.

Return a `DataFrame` containing the columns:

chromosome, rsid, cM, position, allele1, allele2
"""
function read_eigenstrat_snp(snpfile::AbstractString)
    local eigenstrat_snps::DataFrame
    if endswith(snpfile, ".bim")
        eigenstrat_snps = CSV.read(snpfile, DataFrame; header = ["chromosome", "rsid", "cM", "position", "allele1", "allele2"])
    elseif endswith(snpfile, ".snp")
        eigenstrat_snps = CSV.read(snpfile, DataFrame;
            header = ["rsid", "chromosome", "cM", "position", "allele1", "allele2"],
            delim = ' ', ignorerepeated = true)
    else
        throw("read_eigenstrat: Wrong file format! File must end in .bim or .snp.")
    end
    return eigenstrat_snps
end

"""
    write_eigenstrat_snp(filename::AbstractString, snps::DataFrame)

Write .snp file in Eigenstrat format. The SNP file contains information
about each SNP. The SNPs are provided as a `DataFrame` in parameter `snps`.

The DataFrame must consist of the following columns:

chromosome, rsid, cM, position, allele1, allele2
"""
function write_eigenstrat_snp(filename::AbstractString, snps::AbstractDataFrame)
    CSV.write(filename, snps; writeheader = false, delim = ' ')
end

"""
    read_eigenstrat_ind(indfile::AbstractString)

Read individuals from Eigenstrat .ind file.
The IND flle contains information about each individual in the
database.

Return a `DataFrame` consisting of the columns:

ID, Gender, Status

where

Gender: M (male), F (Female) or U (unknown).

Status: Case, Control or population label.
"""
function read_eigenstrat_ind(indfile::AbstractString)
    inds = CSV.read(indfile, DataFrame; header = ["ID", "Gender", "Status"],
        delim = ' ', ignorerepeated = true)
    return inds
end

"""
    write_eigenstrat_ind(filename::AbstractString, inds::AbstractDataFrame)

Write information about each individual to an .ind file.

The parameter `inds` contains information about each individual.

The `DataFrame` must have the columns

ID, Gender, Status
"""
function write_eigenstrat_ind(filename::AbstractString, inds::AbstractDataFrame)
    CSV.write(filename, inds; writeheader = false, delim = ' ')
end

"""
    _encode(
        genotype::Tuple{Char, Char},
        byte::UInt8,
        position::Int64,
        reference_snp::Tuple{Char, Char}
    )

Encode the given genotype in a byte at the given bitpair position (0..3).

Each SNP is characterized by a Tuple of two alleles. The `reference_snp`
contains the reference allele and the derived allele.
"""
function _encode(
    genotype::Tuple{Char, Char},
    byte::UInt8,
    position::Int64,
    reference_snp::Tuple{Char, Char}
)
    # Calculate number of reference alleles.
    n = 0x3 # 3 = no data/invalid
    if (genotype[1] == reference_snp[1]) && (genotype[2] == reference_snp[1])
        n = 0x2
    elseif (genotype[1] == reference_snp[1]) || (genotype[2] == reference_snp[1])
        n = 0x1
    elseif (genotype[1] != reference_snp[1]) && (genotype[2] != reference_snp[1])
        n = 0x0
    end

    # Check for invalid and triallelic markers.
    if (genotype[1] != reference_snp[1]) && (genotype[2] != reference_snp[1]) &&
       (genotype[1] != reference_snp[2]) && (genotype[2] != reference_snp[2])
        n = 0x3
    end

    # Encode value in byte.
    result = byte
    if position == 0
        result = byte | (n << 6)
    elseif position == 1
        result = byte | (n << 4)
    elseif position == 2
        result = byte | (n << 2)
    elseif position == 3
        result = byte | n
    end
    return result
end

"""
    add_individual(
        inprefix::AbstractString,
        outprefix::AbstractString,
        ind_snp_file::AbstractString,
        id::AbstractString;
        gender = "U",
        status = "Control"
    )

Add an individual to a database in Eigenstrat format. 

The SNPs in the database remain untouched. If the individual displayes SNPs
that are not listed in the database or multiallelic ones those SNPs are removed.

`inprefix`: Prefix of the input database.

`outprefix`: Prefix of the output database.

`ind_snp_file`: File containing SNP results for the individual. This should
    work with files from Family Tree DNA Family Finder, MyHeritage, LivingDNA
    and 23andMe.

`id`: ID of the individual. For living persons I recommed the name.

`gender`: U, F or M (Unknown, Female or Male)

`status`: Control, Case or a population label.
"""
function add_individual(
    inprefix::AbstractString,
    outprefix::AbstractString,
    ind_snp_file::AbstractString,
    id::AbstractString;
    gender = "U",
    status = "Control"
)
    indsuffix = ".ind"
    snpsuffix = ".snp"
    genosuffix = ".geno"
    snps = read_eigenstrat_snp(inprefix * snpsuffix)
    nsnp = nrow(snps)
    inds = read_eigenstrat_ind(inprefix * indsuffix)
    nind = nrow(inds)

    # Write .ind file.
    push!(inds, [id, gender, status])
    write_eigenstrat_ind(outprefix * indsuffix, inds)
    ind_hash = hash_ids(outprefix * indsuffix)

    # .snp file remains untouched.
    write_eigenstrat_snp(outprefix * snpsuffix, snps)
    snp_hash = hash_ids(outprefix * snpsuffix)

    # Geno file.

    # Put individual's SNPs into dictionary.
    ind_snps = read_snp_file(ind_snp_file)
    # Put individual's SNPs into a Dictionary.
    ind_dict = Dict{String, String}()
    for snp in eachrow(ind_snps)
        ind_dict[snp.rsid] = snp.genotype
    end

    # 1 SNP value for 4 individuals is encoded as 1 byte.
    bytes_per_line = Int64(ceil(nind / 4))
    # Header size must be at least 48 bytes.
    if bytes_per_line < 48
        bytes_per_line = 48
    end

    # The output must contain enough space for 1 extra individual.
    out_bytes_per_line = Int64(ceil((nind + 1) / 4))
    # Adjust for minimum header size.
    if out_bytes_per_line < bytes_per_line
        out_bytes_per_line = bytes_per_line
    end

    # Read file line by line, add SNP for individual and write to outfile.
    open(inprefix * genosuffix) do infile
        open(outprefix * genosuffix, create = true, write = true) do outfile
            inbuffer = Array{UInt8}(undef, bytes_per_line)
            # Skip first line because it is a header.
            read!(infile, inbuffer)
            # Write new header.
            ihash = string(ind_hash, base = 16)
            shash = string(snp_hash, base = 16)
            outbuf = zeros(UInt8, out_bytes_per_line)
            cols = nind + 1
            header = Array{UInt8}("GENO $cols $nsnp $ihash $shash")
            outbuf[1: length(header)] = header
            write(outfile, outbuf)
            flush(outfile)
            outbuf .= 0
            
            # Write SNPs.
            for i = 1:nsnp
                # Copy input row to output row.
                read!(infile, inbuffer)
                outbuf[1:bytes_per_line] = inbuffer

                # Add individual's SNP value.
                pos = Int64(ceil((nind + 1) / 4))
                byte = outbuf[pos]
                bitpair_no = nind % 4
                rsid = snps.rsid[i]
                reference = (snps.allele1[i][1], snps.allele2[i][1])
                if haskey(ind_dict, rsid)
                    alleles = ind_dict[rsid]
                    genotype = (alleles[1], alleles[2])
                    byte = _encode(genotype, byte, bitpair_no, reference)
                else
                    # missing genotype
                    genotype = ('-', '-')
                    byte = _encode(genotype, byte, bitpair_no, reference)
                end
                outbuf[pos]  = byte
                # Write to file.
                write(outfile, outbuf)
                flush(outfile)
                outbuf .= 0
            end
        end
    end
end

"""
    write_23andMe(filename::AbstractString, snptable)

Write a table of SNPs in 23andMe file format.

This method is included to satisfy users who use
[Plink](https://www.cog-genomics.org/plink2/).
Plink supports 23andMe files.

The snp table must satisfy the
[Tables.jl interface](https://github.com/JuliaData/Tables.jl).

The table must contain the columns:

rsid, chromosome, position, genotype

However exact spelling is not mandatory.
"""
function write_23andMe(filename::AbstractString, snptable)
    CSV.write(filename, snptable; delim = "\t", header = ["#rsid", "chromosome", "position", "genotype"])
end




