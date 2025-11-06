#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
  library(locuszoomr)
  library(EnsDb.Hsapiens.v75)   # brings in ensembldb, etc.
  library(rtracklayer)
})

# ----------------- Paths -----------------
input_dir <- "input_data"

gwas_file <- file.path(
  input_dir,
  "29892016-GCST006085-EFO_0001663-build37.f.tsv.gz"
)

ld_candidates <- c(
  file.path(input_dir, "B3GAT1_LD_r2_clean.tsv"),
  file.path(input_dir, "B3GAT1_LD_r2.tsv")
)
ld_file <- ld_candidates[file.exists(ld_candidates)][1]

recomb_bw <- file.path(
  input_dir,
  "recomb_hg19_hapMapRelease24CombinedRecombMap.bw"
)

if (!file.exists(gwas_file)) {
  stop("GWAS file not found: ", gwas_file)
}
if (is.na(ld_file)) {
  stop("No LD file found. Expected one of:\n  ",
       paste(ld_candidates, collapse = "\n  "))
}
if (!file.exists(recomb_bw)) {
  stop("Recombination .bw file not found: ", recomb_bw)
}

cat("Using GWAS file: ", gwas_file, "\n")
cat("Using LD file:   ", ld_file, "\n")
cat("Using recomb BW: ", recomb_bw, "\n")

# ----------------- 1) Read GWAS -----------------
gwas <- fread(gwas_file)
gwas <- as.data.frame(gwas)

needed_cols <- c(
  "variant_id", "chromosome", "base_pair_location", "p_value"
)
missing <- setdiff(needed_cols, colnames(gwas))
if (length(missing) > 0L) {
  stop("GWAS is missing required columns: ",
       paste(missing, collapse = ", "),
       "\nAvailable: ",
       paste(colnames(gwas), collapse = ", "))
}

cat("GWAS columns (first 12):\n")
print(colnames(gwas)[1:min(12, ncol(gwas))])

# ----------------- 2) Read LD and normalize -----------------
ld_df <- fread(ld_file)
ld_names <- colnames(ld_df)

# Map typical LDassoc columns to canonical names
if ("RS_Number" %in% ld_names && !("variant_id" %in% ld_names)) {
  setnames(ld_df, "RS_Number", "variant_id")
}
if ("SNP" %in% ld_names && !("variant_id" %in% ld_names)) {
  setnames(ld_df, "SNP", "variant_id")
}
if ("R2" %in% ld_names && !("r2" %in% ld_names)) {
  setnames(ld_df, "R2", "r2")
}

if (!all(c("variant_id", "r2") %in% colnames(ld_df))) {
  stop(
    "LD file must have 'variant_id' and 'r2' columns (or mappable names).\n",
    "Current LD columns: ",
    paste(colnames(ld_df), collapse = ", ")
  )
}

ld_df[, variant_id := as.character(variant_id)]
ld_clean <- unique(ld_df[, .(variant_id, r2)], by = "variant_id")

cat("LD rows:", nrow(ld_clean), "\n")

# ----------------- 3) Merge GWAS + LD -----------------
gwas_dt <- as.data.table(gwas)
gwas_dt[, variant_id := as.character(variant_id)]

gwas_ld <- merge(
  gwas_dt,
  ld_clean,
  by = "variant_id",
  all.x = TRUE
)

cat("Variants with non-NA r2:",
    sum(!is.na(gwas_ld$r2)), "\n")

# ----------------- 4) Define B3GAT1 locus & index SNP -----------------
region_chr <- 11L
center_bp  <- 134281332L
flank_bp   <- 150000L

start_bp <- center_bp - flank_bp
end_bp   <- center_bp + flank_bp

cat("Locus window: chr", region_chr, ":",
    start_bp, "-", end_bp, "\n", sep = "")

# Subset to locus window for CSV and for index SNP search
region_sub <- gwas_ld[
  chromosome == region_chr &
    base_pair_location >= start_bp &
    base_pair_location <= end_bp
]

cat("Variants in locus window:", nrow(region_sub), "\n")

if (nrow(region_sub) == 0L) {
  stop("No variants in B3GAT1 window; check coordinates / build.")
}

# Try exact coordinate first
sentinel_row <- region_sub[
  base_pair_location == center_bp
]

if (nrow(sentinel_row) == 0L) {
  message("No SNP exactly at center_bp; using min p-value in window as index SNP.")
  sentinel_row <- region_sub[which.min(p_value)]
}

index_id <- sentinel_row$variant_id[1]
cat("Index SNP variant_id:", index_id, "\n")

# ----------------- 5) Build locuszoomr locus object (for plotting) -----------------
loc <- locus(
  data       = gwas_ld,
  ens_db     = "EnsDb.Hsapiens.v75",
  chrom      = "chromosome",
  pos        = "base_pair_location",
  p          = "p_value",
  labs       = "variant_id",
  index_snp  = index_id,
  flank      = flank_bp,
  LD         = "r2",       # use local LD from LDassoc
  std_filter = TRUE
)

cat("Basic locus summary:\n")
print(summary(loc))

# ----------------- 6) Add recombination from local .bw -----------------
cat("Importing recombination track...\n")
recomb.hg19 <- import.bw(recomb_bw)
loc <- link_recomb(loc, recomb = recomb.hg19)

# ----------------- 7) Write output.csv from region_sub -----------------
output_df <- data.table(
  variant_id   = region_sub$variant_id,
  chromosome   = as.integer(region_sub$chromosome),
  position_bp  = as.integer(region_sub$base_pair_location),
  p_value      = as.numeric(region_sub$p_value),
  neg_log10_p  = -log10(as.numeric(region_sub$p_value)),
  r2           = as.numeric(region_sub$r2),
  is_index_snp = region_sub$variant_id == index_id
)

fwrite(output_df, file = "output.csv")
cat("Wrote output.csv with", nrow(output_df), "rows.\n")

# ----------------- 8) Write output.png -----------------
png("output.png", width = 900, height = 700)
locus_plot(
  loc,
  ylab          = "-log10(p) Prostate cancer",
  main          = "B3GAT1 locus – Prostate cancer GWAS (GCST006085)",
  labels        = "index",
  recomb_offset = 0.1
)
dev.off()

cat("Wrote output.png\n")
