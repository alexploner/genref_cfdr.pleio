#!/usr/bin/bash

# BACKGROUND
# This script downloads genetic data from the 1000 genome project and pre-processes
# it to be used as reference data for the calculatung conditional and conjunctional
# FDR for a pleiotropy-informed GWAS as described in Andreassen et al (2013)
# https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/23375658/
#
# This script generates (1) generic reference data as compressed text files that
# can be used with other implementations (as described below), (2) specific 
# reference data as .rda files for use with the R package pleio.cfdr available 
# from https://github.com/alexploner/pleio.cfdr
#
# The code below is based on pre-processing code given in 
# https://precimed.s3-eu-west-1.amazonaws.com/pleiofdr/about.txt
# which is linked to from the current repository for the original matlab implementation
# of the conditional & cojunctional FDR at https://github.com/precimed/pleiofdr
#
# DATA
# This code is based on version v5a of the 1000 genome project, which was 
# withdrawn in March 2021. As the new version (current: v5b) does not contain
# rs-identifiers, and in order to keep the original reference data replicable,
# this code works with an archived copy of the v5a version.
#
# The original data was downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502
# on 2021-02-09. The archived version of the data is available from 
# https://zenodo.org/record/5750318/files/genref_rawdata.zip, see below. 
#
# USAGE
# Put this shell script and the two related files make_binary_ref.R and 
# magic_intergenic.txt into a directory and run the shell script, e.g.
#
# 	git clone https://github.com/alexploner/genref_pleio.cfdr
# 	cd genref_pleio.cfdr
# 	./make_genref.sh
#
# RESULTS
# Directory ${REFDAT_BINARY_DIR} (see below) and its content can be directly
# used as reference data by cfdr.pleio
#
# REQUIREMENTS
# * plink v1.9  https://www.cog-genomics.org/plink2/ 
# * A recent version of R with packages data.table and R.utils https://cran.r-project.org/
# * parallel, wget, zgrep, unzip
# * pigz (parallel version of gzip, use gzip below if your system does not have it)
# * ca. 50 GB of hard disk space: a total download of ca. 16 GB of data, 
#   generating another 20 GB of data for a total of ca. 35GB of data after 
#   a complete run - however, this will generate some large intermediate files
# * 16 GB RAM *should* be enough, though this has been tested only with 32GB
# * More cores are better (see setup)
#
# FIXME
# * nasty homecooked temp files
#
# LICENSE 
# GPL-3 (as for the original pleiofdr software)
#
# AUTHOR
# For this script: Alexander.Ploner@ki.se (2021-12-1)
# For the original code: Oleksandr Frei https://github.com/ofrei


#---- SETUP ----------

# Saner error handling, also for pipelines
# Also, mirror commands to stdout (slow run)
set -euxo pipefail

# Set number of jobs / threads used with parallel and plink
# See man parallel for details
export n_jobs=`(echo \`nproc\` - 1) | bc` 

# Set directories
export RAW_DIR="./raw"
export CHR_DIR="./chrdata"
export R2_DIR="./r2data"
export LOG_DIR="./log"

# This is the final target directory: its contents can be directly used by R package pleio.fdr
export REFDAT_BINARY_DIR="./bindata"


#------ DOWNLOAD RAW DATA ---------

# Set up directory if required
[ ! -d ${RAW_DIR} ]  && mkdir ${RAW_DIR}


# Download
wget --directory-prefix=${RAW_DIR} https://zenodo.org/record/5750318/files/genref_rawdata.zip
unzip ${RAW_DIR}/genref_rawdata.zip -d ${RAW_DIR}/
rm ${RAW_DIR}/genref_rawdata.zip


#------- PREPROCESS ---------------

# Set up directory if required
[ ! -d ${CHR_DIR} ]  && mkdir ${CHR_DIR}

# Extract the duplicated variants for all 22 chromosomes
seq 22 | parallel --jobs ${n_jobs} "zgrep -v '^#' ${RAW_DIR}/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"\
" | cut -f3 | sort | uniq -d > ${CHR_DIR}/chr{}.dups"

# Get European sample IDs (set family ID = individual ID): all-founder population of unrelated individuals
# This will produce 670 samples, however this will reduce to 503 samples after applying additional filters 
# in the "plink --vcf" step below.
awk 'BEGIN{OFS="\t"; FS="\t"} {if ($7 ~ "IBS|TSI|GBR|CEU|FIN") print($2,$2);}' ${RAW_DIR}/integrated_call_samples_v3.20200731.ALL.ped > samples.eur

# Extract sex for samples
tail -n+2 ${RAW_DIR}/integrated_call_samples_v3.20200731.ALL.ped | awk 'BEGIN{OFS="\t"; FS="\t";} {print($2,$2,$5)}' > samples.sex

# Convert the vcf files for 22 chromosomes to plink format: taking only biallelic reliable variants 
# for European individuals and setting sex from file created in the previous step
seq 22 | parallel --jobs ${n_jobs} \
"plink --vcf ${RAW_DIR}/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz "\
"--biallelic-only strict --out ${CHR_DIR}/chr{} --make-bed --mind 0.1 --geno 0.1 --hwe 1.E-20 midp --maf 0.01 "\
"--keep samples.eur --exclude ${CHR_DIR}/chr{}.dups --update-sex samples.sex"


#-------- AUGMENT VARIANT DESCRIPTION ----------

# Calculate MAFs 
seq 22 | parallel --jobs ${n_jobs} "plink --bfile ${CHR_DIR}/chr{} --freq --out ${CHR_DIR}/chr{}.maf"

# Keep identifier and MAF only, drop the header for the shorter file
seq 22 | parallel --jobs ${n_jobs} "awk 'BEGIN{OFS=\"\t\"} {print(\$2, \$5)}' ${CHR_DIR}/chr{}.maf.frq | tail -n+2 > ${CHR_DIR}/chr{}.maf"

# Join the MAFs to the general variant information (per chromosome)
# We first check that the snp names are identical, then we paste
seq 22 | parallel --jobs ${n_jobs} "cmp --silent <(cut -f2 ${CHR_DIR}/chr{}.bim) <(cut -f1 ${CHR_DIR}/chr{}.maf) "\
"&& paste <(cut -f1,2,4,5,6 ${CHR_DIR}/chr{}.bim) <(cut -f2 ${CHR_DIR}/chr{}.maf) > ${CHR_DIR}/chr{}.ref"\
"|| echo Mismatch for chr{}"

# Create a reference file  containing all SNPs, sorted by chromosome, position & name
echo -e "CHR\tSNP\tBP\tA1\tA2\tMAF" > _tmp_ && cat ${CHR_DIR}/chr*.ref >> _tmp_
sort  -k 1,1n -k 3,3n -k 2,2 _tmp_ > _tmp2_

# Magic file: intergenic flag, see comments in the original reference how this 
# flag was generated
cmp <(cut -f2 <(tail -n+2 _tmp2_ )) <(cut -f1 <(tail -n+2 <(pigz -cd magic_intergenic.txt.gz))) && \
paste _tmp2_ <(cut -f2 <(pigz -cd magic_intergenic.txt.gz)) > all_chr_9524.ref || \
echo key mismatch
rm _tmp_ _tmp2_


#-------- CALCULATE LD R2 -----------------

# Create directory if required
[ ! -d ${R2_DIR} ]  && mkdir ${R2_DIR}

# Calculate r2 for 22 chromosomes: r2 coefficients are calculated based on genotypes {0,1,2}, 
# (for r2 values based on maximum likelihood haplotypes 'dprime' option should be included)
#
# Note that here we loop instead of invoking parallel - not sure why, but parallel only finishes a 
# few small chromosomes and terminates - some kind of race condition? Anyway, this is rather 
# computational intense, and top suggests that plink parallelizes rather efficiently within-task.
#
# Note: this will generate large intermediate files ca. 10-25GB 
#
# Note: using pigz instead of gzip here, to speed up things; use gzip if your system does not have it

for ((i=1;i<23;i++)); do 
	plink --bfile ${CHR_DIR}/chr${i} --threads ${n_jobs} --r2 yes-really \
	--ld-window 1000000 --ld-window-kb 20000 --ld-window-r2 0.05  \
	--out ${R2_DIR}/chr${i}.r2 
	awk '{print($3,$6,$7)}' ${R2_DIR}/chr${i}.r2.ld | pigz -c > ${R2_DIR}/chr${i}.r2.ldshort.gz
	rm ${R2_DIR}/chr${i}.r2.ld
done


#------------- EXPORT LD FILES TO R FORMAT -----------------------

Rscript --file make_binary_ref.R


#---------- CLEAN UP -------------

# Make directory if required
[ ! -d ${LOG_DIR} ]  && mkdir ${LOG_DIR}

# Explicint log files
mv ${CHR_DIR}/*.log ${LOG_DIR}/
mv ${R2_DIR}/*.log ${LOG_DIR}/

# Not really log files, but intermediate
mv samples.eur ${LOG_DIR}/
mv samples.sex ${LOG_DIR}/
mv all_chr_9524.ref ${LOG_DIR}/

