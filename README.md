## Summary

This is a combination of two scripts (one shell, one R) and an auxiliary data file that can generate 
the genomic reference data required for calculating conditional and conunctional false discovery rates 
based on pleiotropic GWAS results. The result files are in R's compressed binary `.RData` format and
intended for use with the R package `cfdr.pleio`.

NOTE: this is only provided to allow replication and modification of the reference data. If you just want to
run `cfdr.pleio`, you can use the the preprocessed, ready-to-use version of the data available from
[https://zenodo.org/record/5750318/files/genref.zip](https://zenodo.org/record/5750318/files/genref.zip)

## Installation

Prepare a project directory and clone the repository:
```
mkdir cfdr_pleio_refdata
cd cfdr_pleio_refdata
git clone https://github.com/alexploner/genref_cfdr.pleio.git .
```

`make_genref.sh` will run parts of the analysis in parallel: by default, it will use the maximum number
of available computing units (as determined by `nproc`) minus one; if you want to override this, you can set the 
variable `n_jobs` in the script manually. 

## Usage

If you are happy with the settings, change into the project directory and run the main script:
```
./make_genref.sh
```

This will download ca. 15.5 GB of raw data and put it through a multi-step processing pipeline, so this 
will not be fast. 

## Results

The preprocessed reference data is stored in the directory that variable `REFDAT_BINARY_DIR` in the main
script `make_genref.sh` points at (default: `bindata` in the project directory). Unless you want to re-run 
the data generation process, this is all you need to keep for e.g. working with `pleio.cfdr`. 


## Requirements

* A working internect connection, for downloading 15.5 GB of raw data
* Ca. 50 GB of free hard disk space (the final reference data are only ca. 3 GB, but the scripts generate large intermediate files) 
* Enough memory (this _may_ work with 16 GB, but has only been tested with 32 GB or more)
* plink v1.9  (available from https://www.cog-genomics.org/plink2/)
* A recent version of `R`, with packages `data.table` and `R.utils`installed (available from  https://cran.r-project.org/)
* A bash environment with 

    * GNU coreutils and utilities `parallel`, `wget`, `zgrep` and `unzip` 
    * `awk`
    * `pigz (replace with `gzip` in `make_genref.sh` if not available)


## Credit

This code is based on a fork (in February 2021) of the pre-processing code given in 

https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt

which is part  current repository for the original matlab implementation
of the conditional & conjunctional FDR at 

https://github.com/precimed/pleiofdr

For more backgroun on theory and application of conditional & conjuncational FDR in
a pleiotropy setting, see [Andreassen et al (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/23375658/)


