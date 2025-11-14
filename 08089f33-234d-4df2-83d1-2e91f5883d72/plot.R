#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

csv_file <- "extracted_solution.csv"

if (!file.exists(csv_file)) {
  stop("Could not find ", csv_file, " in working directory.")
}

dat <- fread(csv_file)

required_cols <- c("chromosome", "position_bp", "minus_log10_p", "r2", "is_index")
missing <- setdiff(required_cols, colnames(dat))
if (length(missing) > 0L) {
  stop("Missing columns in extracted_solution.csv: ",
       paste(missing, collapse = ", "))
}

dat[, pos_Mb := position_bp / 1e6]

p <- ggplot(dat, aes(x = pos_Mb, y = minus_log10_p)) +
  geom_point(aes(color = r2),
             size  = 1.6,
             alpha = 0.9,
             na.rm = FALSE) +
  geom_point(
    data = dat[is_index == TRUE],
    shape = 21,
    size  = 3.2,
    stroke = 0.7,
    colour = "black",
    fill   = "yellow"
  ) +
  scale_color_viridis_c(
    option = "plasma",
    na.value = "grey80",
    name = expression(r^2)
  ) +
  labs(
    x     = "Position on chr11 (Mb)",
    y     = expression(-log[10](p)),
    title = "B3GAT1 locus – Prostate cancer GWAS (GCST006085)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

outfile <- "extracted_solution.jpg"
ggsave(outfile, plot = p, width = 9, height = 7, units = "in", dpi = 300)

cat("Wrote", outfile, "\n")
