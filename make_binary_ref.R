#' make_binary_ref.R
#' 
#' An R script called as part of the shell script make_gztext_ref.sh which
#' takes genetic reference data generated as gzipped text files and turns
#' them into .RData files 
#' 
#' This script is expected to be run via Rscript, see use of Sys.getenv()
#' below
#' 
#' Alexander.Ploner@ki.se 2021-02-12


#' Setup
#' =====

#' Get the name of the chromosome directory from an environment variable
R2DIR <- Sys.getenv("R2_DIR")
if (nchar(R2DIR) == 0) {
  R2DIR <- "./r2data"
  warning("No LD data directory specified via $R2_DIR; falling back to default ", R2DIR)
}

#' Get the name of the target directory from an environment variable
OUTDIR <- Sys.getenv("REFDAT_BINARY_DIR")
if (nchar(OUTDIR) == 0) {
  OUTDIR <- "./bindata"
  warning("No target directory specified via $REFDAT_BINARY_DIR; falling back to default ", OUTDIR)
}

#' By default, do not overwrite exisiting files
FORCE  <- FALSE     ## Remove any files found in OUTDIR? 

#' Check that local default output directory exists, and is empty; creates the
#' directory if necessary, throws error if the directory is not empty and 
#' FORCE = FALSE, otherwise deletes any existing files
if ( !dir.exists(OUTDIR) ) {
  message("Output directory '", OUTDIR, "' does not exist, creating it...")
  dir.create(OUTDIR)
} else {
  if ( length(dir(OUTDIR)) > 0 ) {
    if (FORCE) {
      warning("Removing files found in output directory '", OUTDIR, "'")
    } else {
      stop("Files found in output directory '", OUTDIR, "' - remove, re-think, or set FORCE=TRUE")
    }
  }
}

#' Load required packages; stop if not available
libs <- c("data.table",   ## fast read of text files via fread()
          "R.utils",      ## gzip support for fread()
          "Matrix"        ## sparse data structures
          )
sapply(libs, function(x) stopifnot( require(x, character.only = TRUE) ) )


#' Per-variant reference data
#' ==========================

#' Read the per-variant reference data file
ref1   <- fread(file = "./all_chr_9524.ref")
#' assert: sorted by chromosome, bp, SNP

#' Save this for future work 
fn = file.path(OUTDIR, "all_chr_perVariant.rds")
saveRDS(ref1, file = fn)


#' Split by chromosome for easier processing; get the order by chromosome right,
#' to facilitate looping. Note that we use 22 chromosomes baked in
chrvec   <- as.character(ref1$CHR)
ref1   <- split(ref1, chrvec)
ref1   <- ref1[order(as.numeric(names(ref1)))]

#' LD per-chromosome reference data
#' ================================

#' Name for chromosomes, properly zero-padded
chr_incl <- sort(as.numeric(unique(chrvec)))
chr_namu <- paste0( "chr", (formatC(chr_incl, flag = "0", width = 2)) )

cat( "Chromosome: ")
for (i in 1:22)  {
  cat(i, "... ") ## simple counter, fread shows a progress bar for larger files

  #' Read the zipped text data
  fn  <- file.path(R2DIR, paste0("chr", i, ".r2.ldshort.gz"))
  tmp <- fread(file = fn)
  
  #' Build a sparse matrix
  ndx1 <- match(tmp$SNP_A, ref1[[i]]$SNP)
  stopifnot(all(!is.na(ndx1)))
  ndx2 <- match(tmp$SNP_B, ref1[[i]]$SNP)
  stopifnot(all(!is.na(ndx2)))
  
  #' Drop the tmp object (no longer required) to ease the memory burden
  rm(tmp)
  
  #' Instead of relaying on the symmetric flag for storing the upper triangle
  #' only, and instead of the LD + t(LD), we just add the flipped row- and
  #' column indices to the pair definitions
  LD <- sparseMatrix(i = c(ndx1, ndx2), j = c(ndx2,ndx1), ## no value: logical pattern matrix!
                     dims = rep(nrow(ref1[[i]]), 2),
                     ## Trying to be clever: only save SNP names once, as 
                     ## rownames==colnames; choose colnames for faster 
                     ## by-col access? I think so
                     dimnames = list(NULL, ref1[[i]]$SNP), 
                     ## Do not make symmetric, try to speed it up
                     symmetric = FALSE )

  ## compression factor
  ## as.numeric(object.size(tmp)/object.size(LD))
  
  #' Drop the indices (no longer required) to ease the memory burden
  rm(ndx1, ndx2)
  
  #' Save as binary / compressed data
  fn = file.path(OUTDIR, paste0(chr_namu[i], "_LDpairs.rds"))
  saveRDS(LD, file = fn)
}
