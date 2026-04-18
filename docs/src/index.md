# EigenstratFormat.jl
```@meta
CurrentModule = EigenstratFormat
```

## Introduction

EigenstratFormat.jl is a library for the [Julia](https://julialang.org/)
programming language to work with the Eigenstrat data format.

The Eigenstrat format is often used in genetics. Its definition
can be found at [David Reich's laboratory](https://reich.hms.harvard.edu/software/InputFileFormats)
and in more detail in the [Eigensoft package](https://github.com/DReichLab/EIG).

## Getting started

You can install EigenstratFormat from the Julia interpreter.

```julia
julia> ]
pkg> add EigenstratFormat
pkg> Press BACKSPACE
julia> using EigenstratFormat
```
After that you may want to take a look at the
[code examples](https://github.com/yogischogi/EigenstratFormat.jl/tree/main/examples).

- `01_files.jl` demonstrates basic file operations.
- `02_aadr.jl` shows how to extract samples from the AADR database.
- `03_pca.jl` performs a PCA and creates a simple PCA plot.
- `04_pca_nice.jl` creates a nice looking PCA plot with colored markers and a legend.
- `05_add_sample.jl` adds a sample to the AADR database using a file from a big DNA
   testing company.
- `06_genetic_distances.jl` calculates the genetic distances from one sample to
   all other samples in the database.


## Human DNA samples

Some sources for human DNA samples in Eigenstrat format are:

- [DNA samples](https://reich.hms.harvard.edu/datasets) from single papers at David Reich's laboratory.
- [Allen Ancient DNA Resource (AADR)](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW) 

Most of these samples are ancient, some are from modern pupulations.

## Functions
```@autodocs
Modules = [EigenstratFormat]
```

## Index
```@index
```

