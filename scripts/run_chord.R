# scripts/run_chord.R
suppressPackageStartupMessages({
  library(tools); library(utils); library(methods)
})

args <- list(
  vcf_dir   = snakemake@input[["vcf_dir"]],
  outdir    = snakemake@params[["outdir"]],
  sv_caller = snakemake@params[["sv_caller"]],
  recursive = as.logical(snakemake@params[["recursive"]]),
  ref_ucsc  = snakemake@params[["ref_ucsc"]]
)

log_msg <- function(...) { cat("[CHORD]", sprintf(...), "\n") }

ensure_pkg <- function(pkg, github=NULL, bioc=NULL){
  if (!requireNamespace(pkg, quietly=TRUE)){
    if (!is.null(github)){
      if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org")
      remotes::install_github(github, upgrade="never", dependencies=TRUE)
    } else if (!is.null(bioc)){
      if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org")
      BiocManager::install(bioc, update=FALSE, ask=FALSE)
    } else {
      install.packages(pkg, repos="https://cloud.r-project.org")
    }
  }
}

ensure_pkg("data.table"); ensure_pkg("stringr")
ensure_pkg("CHORD", github="UMCUGenetics/CHORD")
ensure_pkg("mutSigExtractor", github="UMCUGenetics/mutSigExtractor")

suppressPackageStartupMessages({
  library(data.table); library(stringr); library(CHORD); library(mutSigExtractor)
})

if (args$ref_ucsc == "hg38") {
  ensure_pkg("BSgenome.Hsapiens.UCSC.hg38", bioc="BSgenome.Hsapiens.UCSC.hg38")
  library(BSgenome.Hsapiens.UCSC.hg38); ref.genome <- BSgenome.Hsapiens.UCSC.hg38
} else if (args$ref_ucsc == "hg19") {
  ensure_pkg("BSgenome.Hsapiens.UCSC.hg19", bioc="BSgenome.Hsapiens.UCSC.hg19")
  library(BSgenome.Hsapiens.UCSC.hg19); ref.genome <- BSgenome.Hsapiens.UCSC.hg19
} else { stop("ref_ucsc must be 'hg38' or 'hg19'") }

vcf_dir <- args$vcf_dir; outdir <- args$outdir
dir.create(file.path(outdir, "contexts"), recursive=TRUE, showWarnings=FALSE)
contexts_dir <- file.path(outdir, "contexts")

snv_indel_files <- list.files(vcf_dir, pattern="snv_indel\.vcf(\.gz)?$", full.names=TRUE, recursive=args$recursive)
sv_files        <- list.files(vcf_dir, pattern="sv\.vcf(\.gz)?$",        full.names=TRUE, recursive=args$recursive)
if (length(snv_indel_files) == 0 || length(sv_files) == 0) stop("No VCFs found in: ", vcf_dir)

samples <- basename(snv_indel_files)
sample_names <- sapply(strsplit(samples, "_"), `[`, 1)
dt <- data.table(sample=sample_names, snv_indel=snv_indel_files)
dt[, sv := vapply(sample, function(s){
  cand <- sv_files[grepl(paste0("^", s, "_"), basename(sv_files))]
  if (length(cand)==0) NA_character_ else cand[1]
}, character(1))]
dt <- dt[!is.na(sv)]
if (nrow(dt)==0) stop("No matching SNV/indel and SV VCF pairs found.")

log_msg("Found %d sample(s) with both SNV/indel and SV VCFs.", nrow(dt))

for (i in seq_len(nrow(dt))) {
  snv <- dt$snv_indel[i]; sv <- dt$sv[i]; sm <- dt$sample[i]
  out_path <- file.path(contexts_dir, paste0(sm, "_contexts.txt"))
  if (!file.exists(out_path)) {
    log_msg("Extracting contexts for sample: %s", sm)
    tryCatch({
      extractSigsChord(
        vcf.snv      = snv,
        vcf.sv       = sv,
        sv.caller    = args$sv_caller,
        sample.name  = sm,
        output.path  = out_path,
        ref.genome   = ref.genome,
        verbose      = FALSE
      )
    }, error=function(e){
      message(sprintf("Error processing sample %s: %s", sm, e$message))
    })
  }
}

context_files <- list.files(contexts_dir, full.names=TRUE, pattern="_contexts\.txt$")
if (length(context_files) == 0) stop("No contexts files produced.")
l_contexts <- lapply(context_files, function(p) data.table::fread(p))
merged_contexts <- data.table::rbindlist(l_contexts, use.names=TRUE, fill=TRUE)

data.table::fwrite(merged_contexts, file.path(outdir, "merged_contexts.txt"), sep="\t", quote=FALSE)

chord_output <- CHORD::chordPredict(merged_contexts, do.bootstrap=TRUE, verbose=FALSE)
data.table::fwrite(chord_output, file.path(outdir, "chord_pred.txt"), sep="\t", quote=FALSE)

log_msg("Finished. Outputs written to: %s", outdir)
