#'*BAD SCRIPT*
#'*BC the internal cores can failed*

library(data.table)
library(fread)


foldername = "/n/groups/patel/shakson_ukb/UK_Biobank/BScripts/Module5/deCODE/slurm/"
jobid = 28695110

phrase <- "Done with providing MR analysis of UKB exposures, deCODE proteomics, FinnGen Diseases"

files <- sprintf("%sModule5_%s_%d.out", foldername, jobid, 1:2000)

has_phrase <- vapply(files, function(f) {
  if (!file.exists(f)) return(NA)            # missing file
  any(grepl(phrase, readLines(f, warn = FALSE), fixed = TRUE))
}, logical(1))

# indices finished (have the phrase)
which(has_phrase %in% TRUE)
length(which(has_phrase %in% TRUE))

# missing files
which(is.na(has_phrase))

# not finished (file exists but phrase not found)
which(has_phrase %in% FALSE)
