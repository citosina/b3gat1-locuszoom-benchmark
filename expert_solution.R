#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
  library(locuszoomr)
  library(ensembldb)
  library(EnsDb.Hsapiens.v75)
  library(rtracklayer)
})

# ---------- Paths ----------
input_dir  <- "input_data"
gwas_file  <- file.path(input_dir, "29892016-GCST006085-EFO_0001663-build37.f.tsv.gz")
ld_file    <- file.path(input_dir, "B3GAT1_LD_r2_clean.tsv")
recomb_bw  <- file.path(input_dir, "recomb_hg19_hapMapRelease24CombinedRecombMap.bw")

# ---------- Region parameters ----------
region_chr <- 11L
center_bp  <- 134281332L
flank_bp   <- 150000L
start_bp   <- center_bp - flank_bp
end_bp     <- center_bp + flank_bp
index_snp_preferred <- "rs78760579"

# ---------- Sanity: files exist ----------
for (f in c(gwas_file, ld_file, recomb_bw)) {
  if (!file.exists(f)) {
    stop("Required file not found: ", f)
  }
}

# ---------- 1. Load GWAS ----------
cat("Reading GWAS from:", gwas_file, "\n")
gwas <- fread(gwas_file)

if (!("variant_id" %in% names(gwas))) {
  stop("GWAS file must contain a 'variant_id' column. Found: ", paste(names(gwas), collapse = ", "))
}
if (!all(c("chromosome", "base_pair_location", "p_value") %in% names(gwas))) {
  stop("GWAS file must contain 'chromosome', 'base_pair_location', and 'p_value' columns.")
}

gwas[, variant_id := as.character(variant_id)]

# ---------- 2. Load LD file and normalize columns ----------
cat("Reading LD from:", ld_file, "\n")
ld <- fread(ld_file)

# Try to standardize to variant_id + r2
ld_names <- names(ld)

# Map common LDassoc formats to our standard names
if ("RS_Number" %in% ld_names && !("variant_id" %in% ld_names)) {
  setnames(ld, "RS_Number", "variant_id")
}
if ("SNP" %in% ld_names && !("variant_id" %in% ld_names)) {
  setnames(ld, "SNP", "variant_id")
}
if ("R2" %in% ld_names && !("r2" %in% ld_names)) {
  setnames(ld, "R2", "r2")
}

if (!all(c("variant_id", "r2") %in% names(ld))) {
  stop(
    "LD file must have 'variant_id' and 'r2' columns (or something mappable, e.g. 'RS_Number' + 'R2'). ",
    "Current columns: ", paste(names(ld), collapse = ", ")
  )
}

ld[, variant_id := as.character(variant_id)]

# Keep only what we need
ld <- ld[, .(variant_id, r2)]

# ---------- 3. Merge GWAS + LD ----------
cat("Merging GWAS with LD by variant_id...\n")
gwas_ld <- merge(
  gwas,
  ld,
  by = "variant_id",
  all.x = TRUE
)

cat("Variants with non-NA r2:", sum(!is.na(gwas_ld$r2)), "\n")

# ---------- 4. Subset to B3GAT1 locus ----------
region <- gwas_ld[
  chromosome == region_chr &
    base_pair_location >= start_bp &
    base_pair_location <= end_bp
]

cat("Variants in region:", nrow(region), "\n")

if (nrow(region) == 0L) {
  stop("No variants found in the specified B3GAT1 locus window.")
}

# ---------- 5. Choose index SNP ----------
idx_row <- region[variant_id == index_snp_preferred]

if (nrow(idx_row) == 0L) {
  message("Preferred index SNP rs78760579 not found in region; using min p-value as index SNP.")
  idx_row <- region[which.min(p_value)]
}

index_id <- idx_row$variant_id[1]
cat("Using index SNP:", index_id, "\n")

# ---------- 6. Build locuszoomr locus object ----------
loc <- locus(
  data       = region,
  ens_db     = "EnsDb.Hsapiens.v75",
  chrom      = "chromosome",
  pos        = "base_pair_location",
  p          = "p_value",
  labs       = "variant_id",
  index_snp  = index_id,
  flank      = flank_bp,
  LD         = "r2",
  std_filter = FALSE
)

# ---------- 7. Add recombination rate from local bigWig ----------
cat("Importing recombination track from:", recomb_bw, "\n")
recomb_hg19 <- import.bw(recomb_bw)
loc <- link_recomb(loc, recomb = recomb_hg19)

# ---------- 8. Write output.csv ----------
loc_data <- loc$data

output_df <- data.table(
  variant_id   = loc_data$variant_id,
  chromosome   = as.integer(loc_data$chromosome),
  position_bp  = as.integer(loc_data$base_pair_location),
  p_value      = as.numeric(loc_data$p_value),
  neg_log10_p  = -log10(as.numeric(loc_data$p_value)),
  r2           = as.numeric(loc_data$r2),
  is_index_snp = loc_data$variant_id == index_id
)

fwrite(output_df, file = "output.csv")
cat("Wrote output.csv with", nrow(output_df), "rows.\n")

# ---------- 9. Write output.png ----------
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
