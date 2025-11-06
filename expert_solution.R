#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(locuszoomr)
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

# ---------- Load data ----------
if (!file.exists(gwas_file)) {
  stop("GWAS file not found: ", gwas_file)
}
if (!file.exists(ld_file)) {
  stop("LD file not found: ", ld_file)
}
if (!file.exists(recomb_bw)) {
  stop("Recombination map file not found: ", recomb_bw)
}

gwas <- fread(gwas_file)
ld   <- fread(ld_file)

gwas[, variant_id := as.character(variant_id)]
ld[, variant_id   := as.character(variant_id)]

# Merge LD into GWAS by variant_id
gwas_ld <- merge(
  gwas,
  ld,
  by = "variant_id",
  all.x = TRUE
)

# ---------- Subset to locus ----------
region <- gwas_ld[
  chromosome == region_chr &
    base_pair_location >= start_bp &
    base_pair_location <= end_bp
]

if (nrow(region) == 0L) {
  stop("No variants found in the specified B3GAT1 locus window.")
}

# ---------- Choose index SNP ----------
idx_row <- region[variant_id == index_snp_preferred]

if (nrow(idx_row) == 0L) {
  message("Preferred index SNP rs78760579 not found; using min p-value in region.")
  idx_row <- region[which.min(p_value)]
}

index_id <- idx_row$variant_id[1]
cat("Using index SNP:", index_id, "\n")

# ---------- Build locuszoomr locus object ----------
loc <- locus(
  data       = region,
  ens_db     = "EnsDb.Hsapiens.v75",    # hg19 / GRCh37
  chrom      = "chromosome",
  pos        = "base_pair_location",
  p          = "p_value",
  labs       = "variant_id",
  index_snp  = index_id,
  flank      = flank_bp,
  LD         = "r2",                    # use local LD column
  std_filter = FALSE                    # already region-filtered
)

# ---------- Add recombination rate from local bigWig ----------
recomb_hg19 <- import.bw(recomb_bw)
loc <- link_recomb(loc, recomb = recomb_hg19)

# ---------- Write output.csv ----------
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

# ---------- Write output.png ----------
png("output.png", width = 900, height = 700)
locus_plot(
  loc,
  ylab          = "-log10(p) Prostate cancer",
  main          = "B3GAT1 locus – Prostate cancer GWAS (GCST006085)",
  labels        = "index",
  recomb_offset = 0.1
)
dev.off()

cat("Wrote output.csv and output.png\n")
