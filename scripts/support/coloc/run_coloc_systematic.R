#!/usr/bin/env Rscript
# ============================================================================
# support/coloc/run_coloc_systematic.R
# ----------------------------------------------------------------------------
# Systematic colocalization (coloc.abf) over the cis-pQTL edges that the tier
# table flagged `coloc_status == "pending"` (lane==pQTL & edge_class==cis &
# Tier1). Two flavours of locus, both anchored on the protein's cis lead SNP
# (±window):
#   Pcis_to_D : pQTL (quant) x FinnGen disease (cc, s from manifest)
#   Pcis_to_E : pQTL (quant) x exposure REGENIE GWAS (quant)
#
# Engine logic (harmonize + coloc.abf) ported verbatim from the working legacy
# ModuleMR/COLOC/runColoc.R. Adds a REGENIE reader (UKB pQTL + exposure GWAS
# share the REGENIE step2 schema) alongside the deCODE/FinnGen readers.
#
# Output (canonical):
#   support/coloc/per_locus/<arm>__<prot>__<target>_coloc_summary.tsv
#   support/coloc/coloc_results.tsv     one row per locus: ...,PP.H3,PP.H4,status
#
# Usage:
#   module load gcc/14.2.0 R/4.4.2
#   HEAP_PATHS_FILE=.../workflow/00_paths.R \
#     Rscript scripts/support/coloc/run_coloc_systematic.R [--only <protID>] [--window_kb 500]
# ============================================================================
suppressPackageStartupMessages({ library(data.table); library(coloc) })
local({
  cand <- c(Sys.getenv("HEAP_PATHS_FILE", unset = ""),
            "/n/groups/patel/shakson_ukb/HEAP/workflow/00_paths.R")
  hit <- cand[nzchar(cand) & file.exists(cand)][1]
  if (is.na(hit)) stop("Cannot locate workflow/00_paths.R (set HEAP_PATHS_FILE).")
  source(hit)
})
have_arrow <- requireNamespace("arrow", quietly = TRUE)
have_dplyr <- requireNamespace("dplyr", quietly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  w <- which(args == flag); if (!length(w)) return(default); args[w + 1] }
ONLY      <- get_arg("--only", NA_character_)
WINDOW_KB <- suppressWarnings(as.numeric(get_arg("--window_kb", "500")))

IGLOO_HEAP <- heap_project_output()                       # .../UKB/HEAP/output
OUT_DIR    <- heap_project_output("support", "coloc")
PL_DIR     <- file.path(OUT_DIR, "per_locus")
dir.create(PL_DIR, recursive = TRUE, showWarnings = FALSE)

FINNGEN_DIR  <- igloo_path("FinnGen", "SummaryStats")
FINNGEN_MAN  <- fread(igloo_path("FinnGen", "finngen_R12_manifest.tsv"))
SOMA_MAP     <- fread(igloo_path("DECODE", "pQTLmetadata", "somascan_protein_map.tsv"))
DECODE_PQTL  <- igloo_path("DECODE", "pQTL", "final_somascan_smp")
UKB_PQTL     <- igloo_path("UKB", "pQTL")
EXPO_DIR     <- heap_gwas("regenie_step2")
CIS_UKB      <- heap_project_output("mr", "protein_inst")
CIS_DEC      <- heap_project_output("mr", "protein_inst_decode")

# ---------------------------------------------------------------- readers ----
set_if <- function(DT, from, to) { if (from %in% names(DT) && !(to %in% names(DT))) setnames(DT, from, to); DT }
to_maf <- function(x){ v<-suppressWarnings(as.numeric(x)); ifelse(is.finite(v), pmin(v,1-v), NA_real_) }
is_palin <- function(a1,a2){ a1<-toupper(a1);a2<-toupper(a2); (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C") }
infer_N <- function(x){ x<-suppressWarnings(as.numeric(x)); x<-x[is.finite(x)&x>0]; if(!length(x)) return(NA_real_); ux<-unique(x); ux[which.max(tabulate(match(x,ux)))] }

# Region-filtered TEXT read: stream the file through awk keeping only the header +
# rows on chr0 within [lo,hi], selecting wanted cols. Keeps memory at ~window size
# (MBs) instead of the genome-wide table (GBs) that OOM'd the 50G cgroup.
read_text_region <- function(path, chrcol, poscol, chr0, lo, hi, want) {
  hdr <- names(fread(path, nrows = 0))
  ci <- match(chrcol, hdr); pi <- match(poscol, hdr)
  if (is.na(ci) || is.na(pi)) stop("chr/pos col missing in ", basename(path))
  rdr <- if (grepl("\\.gz$", path)) "zcat" else "cat"
  # FS = one-or-more space/tab so this works for BOTH tab-delimited (FinnGen/deCODE)
  # and space-delimited (REGENIE exposure GWAS) files; $0 is printed unchanged so
  # fread auto-detects the original separator.
  awk <- sprintf("%s %s | awk -F'[ \\t]+' 'NR==1 || ((($%d==\"%s\")||($%d==\"chr%s\")) && ($%d+0)>=%d && ($%d+0)<=%d)'",
                 rdr, shQuote(path), ci, chr0, ci, chr0, pi, lo, pi, hi)
  fread(cmd = awk, select = intersect(want, hdr), showProgress = FALSE)
}

read_regenie <- function(path, chr0, lo, hi) {   # UKB pQTL parquet OR exposure REGENIE text, WINDOWED
  want <- c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","N","BETA","SE","LOG10P","rsid")
  if (grepl("\\.parquet$", path)) {
    if (!have_arrow) stop("arrow needed for parquet: ", path)
    sc <- names(arrow::open_dataset(path)$schema)
    ds <- arrow::open_dataset(path)
    DT <- if (have_dplyr)
            as.data.table(dplyr::collect(dplyr::filter(ds, GENPOS >= lo & GENPOS <= hi)))
          else as.data.table(arrow::read_parquet(path, col_select = intersect(want, sc)))
    DT <- DT[, intersect(want, names(DT)), with = FALSE]
    if ("CHROM" %in% names(DT)) DT <- DT[as.character(CHROM)==as.character(chr0) | as.character(CHROM)==paste0("chr",chr0)]
    DT <- DT[GENPOS >= lo & GENPOS <= hi]
  } else {
    DT <- read_text_region(path, "CHROM", "GENPOS", chr0, lo, hi, want)
  }
  set_if(DT, "rsid", "SNP"); if (!"SNP" %in% names(DT)) set_if(DT, "ID", "SNP")
  set_if(DT, "CHROM", "chr"); set_if(DT, "GENPOS", "pos")
  set_if(DT, "BETA", "beta"); set_if(DT, "SE", "se")
  set_if(DT, "ALLELE1", "effect_allele"); set_if(DT, "ALLELE0", "other_allele")  # REGENIE: ALLELE1 = tested
  set_if(DT, "A1FREQ", "eaf"); set_if(DT, "N", "N")
  if ("LOG10P" %in% names(DT) && !"pval" %in% names(DT)) DT[, pval := 10^(-as.numeric(LOG10P))]
  DT[, chr := gsub("^chr","",as.character(chr))]; DT[, pos := as.integer(pos)]
  DT
}
read_decode <- function(path, chr0, lo, hi) {
  want <- c("rsids","Name","Chrom","Pos","Beta","SE","Pval","effectAllele","otherAllele","ImpMAF","N")
  DT <- read_text_region(path, "Chrom", "Pos", chr0, lo, hi, want)
  if ("rsids" %in% names(DT)) { DT[, SNP := tstrsplit(rsids, ",", fixed=TRUE, keep=1)]; DT[SNP %in% c("",".",NA), SNP := NA_character_] }
  if (!"SNP" %in% names(DT) || all(is.na(DT$SNP))) set_if(DT, "Name", "SNP")
  set_if(DT,"Chrom","chr"); if ("chr"%in%names(DT)) DT[,chr:=gsub("^chr","",chr)]
  set_if(DT,"Pos","pos"); set_if(DT,"Beta","beta"); set_if(DT,"SE","se"); set_if(DT,"Pval","pval")
  set_if(DT,"effectAllele","effect_allele"); set_if(DT,"otherAllele","other_allele")
  set_if(DT,"ImpMAF","eaf"); set_if(DT,"N","N")
  DT[, chr := as.character(chr)]; DT[, pos := as.integer(pos)]
  DT
}
read_finngen <- function(path, chr0, lo, hi) {
  want <- c("#chrom","pos","ref","alt","rsids","rsid","beta","sebeta","pval","af_alt","n_total")
  DT <- read_text_region(path, "#chrom", "pos", chr0, lo, hi, want)
  if (!"rsid" %in% names(DT)) { if ("rsids"%in%names(DT)) DT[,rsid:=tstrsplit(rsids,",",fixed=TRUE,keep=1)] }
  set_if(DT,"rsid","SNP"); if (!"SNP"%in%names(DT) && "#chrom"%in%names(DT)) DT[,SNP:=paste0(get("#chrom"),":",pos,"_",ref,"_",alt)]
  set_if(DT,"sebeta","se"); set_if(DT,"alt","effect_allele"); set_if(DT,"ref","other_allele")
  set_if(DT,"af_alt","eaf"); set_if(DT,"#chrom","chr")
  DT[, chr := gsub("^chr","",as.character(chr))]; DT[, pos := as.integer(pos)]
  DT
}

dedup_snp <- function(DT){ DT<-as.data.table(DT); DT[, .pr := if("pval"%in%names(DT)) suppressWarnings(as.numeric(pval)) else NA_real_]
  setorder(DT, SNP, .pr, na.last=TRUE); out<-DT[, .SD[1], by=SNP]; out[, .pr:=NULL]; out }
harmonize_two <- function(d1,d2){
  d1<-dedup_snp(d1); d2<-dedup_snp(d2)
  m<-merge(d1,d2,by="SNP",suffixes=c(".1",".2"),allow.cartesian=FALSE); if(!nrow(m)) return(NULL)
  for (c in c("effect_allele.1","other_allele.1","effect_allele.2","other_allele.2")) m[[c]]<-toupper(m[[c]])
  same<-m$effect_allele.1==m$effect_allele.2 & m$other_allele.1==m$other_allele.2
  swap<-m$effect_allele.1==m$other_allele.2  & m$other_allele.1==m$effect_allele.2
  keep<-same|swap; m<-m[keep]; swap<-swap[keep]; if(!nrow(m)) return(NULL)
  if (any(swap)) m[swap, beta.2 := -as.numeric(beta.2)]
  pal<-is_palin(m$effect_allele.1,m$other_allele.1)
  if ("eaf.1"%in%names(m)) { e1<-suppressWarnings(as.numeric(m$eaf.1)); m<-m[!(pal & is.finite(e1) & e1>0.42 & e1<0.58)] } else m<-m[!pal]
  m
}
mk_ds <- function(m, which, type, s=NULL, N=NULL){
  beta<-suppressWarnings(as.numeric(m[[paste0("beta.",which)]])); se<-suppressWarnings(as.numeric(m[[paste0("se.",which)]]))
  eaf_col<-paste0("eaf.",which); maf<-if(eaf_col%in%names(m)) to_maf(m[[eaf_col]]) else rep(NA_real_,nrow(m))
  ok <- is.finite(beta)&is.finite(se)&se>0
  ds<-list(snp=m$SNP[ok], beta=beta[ok], varbeta=(se[ok])^2, MAF=maf[ok], type=type)
  if (!is.null(N)) ds$N<-N
  if (type=="cc") ds$s<-s
  ds
}

# ------------------------------------------------------------ locus table ----
read_pending <- function(armtag, f) {
  if (!file.exists(f)) return(NULL)
  d <- fread(f)
  d <- d[coloc_status=="pending" & edge_class=="cis"]
  if (!nrow(d)) return(NULL)
  data.table(arm=armtag, protID=d$src_id, target=d$tgt_id, edge_dir=d$edge_dir)
}
loci <- rbindlist(list(
  read_pending("UKB",    file.path(IGLOO_HEAP,"mr_edges","summary","mr_sensitivity_long.tsv")),
  read_pending("DECODE", file.path(IGLOO_HEAP,"mr_edges","summary","DECODE","mr_sensitivity_long.tsv"))
), use.names=TRUE)
loci <- unique(loci)
if (!is.na(ONLY)) loci <- loci[protID==ONLY]
cat("Pending cis loci:", nrow(loci), "(", loci[,.N,arm][,paste(arm,N,collapse=", ")], ")\n")

# resolve pQTL file + cis lead SNP per (arm,protID)
ukb_pq_file <- function(g){ f<-list.files(UKB_PQTL, pattern=paste0("^",g,"_.*\\.parquet$"), full.names=TRUE); if(length(f)) f[1] else NA_character_ }
dec_pq_file <- function(g){ r<-SOMA_MAP[EntrezGeneSymbol==g | hgnc_symbol==g]; if(nrow(r)) file.path(DECODE_PQTL, r$file[1]) else NA_character_ }
lead_info <- function(arm,g){ d<-if(arm=="DECODE") CIS_DEC else CIS_UKB; f<-file.path(d, paste0("protein_",g,"_cis_clumped.tsv"))
  if(!file.exists(f)) return(NULL); x<-fread(f); if(!nrow(x) || !"chr.exposure"%in%names(x)) return(NULL)
  i<-which.min(suppressWarnings(as.numeric(x$pval.exposure))); if(!length(i)) return(NULL)
  list(snp=x$SNP[i], chr=gsub("^chr","",as.character(x$chr.exposure[i])), pos=suppressWarnings(as.integer(x$pos.exposure[i]))) }

fg_pheno <- function(t) sub("^finngen_R12_","",t)
fg_file  <- function(t) file.path(FINNGEN_DIR, paste0(t,".gz"))
expo_file <- function(t){ c1<-file.path(EXPO_DIR,t,paste0(t,".regenie")); c2<-file.path(EXPO_DIR,paste0(t,".regenie"))
  if(file.exists(c1)) c1 else if(file.exists(c2)) c2 else { g<-Sys.glob(file.path(EXPO_DIR,t,"*.regenie")); if(length(g)) g[1] else NA_character_ } }

# ------------------------------------------------------------------ run ------
# Region-filtered reads (cis ±window only, via the lead's chr/pos from the cis
# cache) => ~MB per read, no genome-wide tables, no caches. Memory-trivial.
loci[, out_kind := fifelse(edge_dir=="Pcis_to_D","disease","exposure")]
results <- list()
setorder(loci, out_kind, target, protID)
# Optional contiguous chunking for on-node parallelism (preserves target grouping
# so each outcome file is still read once within a chunk). Parallel chunks write
# only per-locus files; aggregate per_locus/* afterwards.
NCH <- suppressWarnings(as.integer(get_arg("--nchunks","1")))
CH  <- suppressWarnings(as.integer(get_arg("--chunk","1")))
if (is.finite(NCH) && NCH > 1) {
  n <- nrow(loci); sz <- ceiling(n/NCH); lo <- (CH-1)*sz+1; hi <- min(CH*sz, n)
  if (lo > n) { cat("chunk", CH, "of", NCH, "-> 0 loci\n"); quit(save="no") }
  loci <- loci[lo:hi]; cat("chunk", CH, "of", NCH, "->", nrow(loci), "loci (rows", lo, "-", hi, ")\n")
}
for (i in seq_len(nrow(loci))) {
  L <- loci[i]; tag <- sprintf("%s__%s__%s", L$arm, L$protID, L$target)
  pl_f <- file.path(PL_DIR, paste0(tag, "_coloc_summary.tsv"))
  if (file.exists(pl_f)) {                                  # idempotent: reuse a prior run's result
    pr <- tryCatch(fread(pl_f), error=function(e) NULL)
    if (!is.null(pr) && nrow(pr)) {
      g <- function(c) if (c %in% names(pr)) pr[[c]][1] else NA
      results[[i]] <- data.table(arm=L$arm,protID=L$protID,target=L$target,edge_dir=L$edge_dir,
                                 lead_snp=g("lead_snp"), nsnps=g("nsnps"), PP.H3=g("PP.H3"),
                                 PP.H4=g("PP.H4"), status="OK_cached")
      cat(sprintf("[%3d/%3d] %-42s CACHED  PP.H4=%s\n", i, nrow(loci), tag,
                  ifelse("PP.H4"%in%names(pr), sprintf("%.3f", as.numeric(pr$PP.H4[1])), "NA")))
      next
    }
  }
  st <- "OK"; ph3<-ph4<-NA_real_; nsnp<-0L; lead<-NA_character_
  res <- tryCatch({
    li <- lead_info(L$arm, L$protID)
    if (is.null(li) || is.na(li$chr) || is.na(li$pos)) { st<-"missing_lead"; stop("no lead") }
    lead <- li$snp; chr0 <- li$chr; pos0 <- li$pos
    lo <- max(1L, pos0 - WINDOW_KB*1000L); hi <- pos0 + WINDOW_KB*1000L
    pqf <- if (L$arm=="DECODE") dec_pq_file(L$protID) else ukb_pq_file(L$protID)
    if (is.na(pqf) || !file.exists(pqf)) { st<-"missing_pqtl"; stop("no pqtl") }
    pqr <- if (L$arm=="DECODE") read_decode(pqf, chr0, lo, hi) else read_regenie(pqf, chr0, lo, hi)
    if (is.null(pqr) || !nrow(pqr)) { st<-"pqtl_empty"; stop("pqtl empty") }
    if (L$out_kind=="disease") { of<-fg_file(L$target); if(!file.exists(of)){st<-"missing_gwas";stop("no fg")}; otr<-read_finngen(of, chr0, lo, hi) }
    else { of<-expo_file(L$target); if(is.na(of)){st<-"missing_expo_gwas";stop("no expo")}; otr<-read_regenie(of, chr0, lo, hi) }
    if (is.null(otr) || !nrow(otr)) { st<-"outcome_empty"; stop("outcome empty") }
    m<-harmonize_two(pqr, otr); if (is.null(m) || nrow(m)<50) { st<-"too_few_snps"; stop("few snps") }
    nsnp<-nrow(m)
    Npq<-infer_N(if("N.1"%in%names(m)) m$N.1 else NA); if(!is.finite(Npq)) Npq<-infer_N(pqr$N)
    d1<-mk_ds(m,"1","quant",N=if(is.finite(Npq))Npq else NULL)
    if (L$out_kind=="disease") {
      ph<-fg_pheno(L$target); mr<-FINNGEN_MAN[phenocode==ph]
      if(!nrow(mr)){ st<-"no_manifest"; stop("no manifest") }
      nca<-as.numeric(mr$num_cases[1]); nco<-as.numeric(mr$num_controls[1]); Ntot<-nca+nco; s<-nca/Ntot
      d2<-mk_ds(m,"2","cc",s=s,N=Ntot)
    } else {
      Nx<-infer_N(if("N.2"%in%names(m)) m$N.2 else NA); if(!is.finite(Nx)) Nx<-infer_N(otr$N)
      d2<-mk_ds(m,"2","quant",N=if(is.finite(Nx))Nx else NULL)
    }
    cc<-coloc::coloc.abf(dataset1=d1, dataset2=d2)
    sm<-cc$summary; ph3<-unname(sm["PP.H3.abf"]); ph4<-unname(sm["PP.H4.abf"])
    fwrite(data.table(arm=L$arm,protID=L$protID,target=L$target,edge_dir=L$edge_dir,lead_snp=lead,
                      chr=chr0,pos=pos0,nsnps=nsnp,PP.H3=ph3,PP.H4=ph4,status="OK"),
           file.path(PL_DIR, paste0(tag,"_coloc_summary.tsv")), sep="\t")
    TRUE
  }, error=function(e) { FALSE })
  results[[i]] <- data.table(arm=L$arm,protID=L$protID,target=L$target,edge_dir=L$edge_dir,
                             lead_snp=lead,nsnps=nsnp,PP.H3=ph3,PP.H4=ph4,status=st)
  cat(sprintf("[%3d/%3d] %-42s %-12s nsnp=%-5s PP.H4=%-7s %s\n", i, nrow(loci), tag, L$edge_dir,
              nsnp, ifelse(is.na(ph4),"NA",sprintf("%.3f",ph4)), st))
}
R <- rbindlist(results)
# Only the single-process run writes the canonical results table; parallel chunks
# leave per-locus files for a separate aggregation pass (avoids overwrites).
if (!is.finite(NCH) || NCH <= 1) fwrite(R, file.path(OUT_DIR, "coloc_results.tsv"), sep="\t") else fwrite(R, file.path(OUT_DIR, paste0("coloc_results_chunk", CH, ".tsv")), sep="\t")
cat("\n=== SUMMARY ===\n")
cat("loci attempted:", nrow(R), "\n")
print(R[, .N, by=status][order(-N)])
cat("colocalized (PP.H4>=0.8):", R[is.finite(PP.H4) & PP.H4>=0.8, .N], "\n")
cat("not colocalized (PP.H4<0.8):", R[is.finite(PP.H4) & PP.H4<0.8, .N], "\n")
cat("wrote:", file.path(OUT_DIR, "coloc_results.tsv"), "\n")
