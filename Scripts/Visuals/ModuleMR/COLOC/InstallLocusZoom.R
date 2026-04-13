if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}

BiocManager::install(
  c("AnnotationFilter", "ensembldb", "GenomicRanges", "rtracklayer"),
  ask = FALSE, update = FALSE
)

install.packages("locuszoomr", repos="https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos="https://cloud.r-project.org")
}

BiocManager::install("EnsDb.Hsapiens.v86", ask = FALSE, update = FALSE)
library(EnsDb.Hsapiens.v86)

