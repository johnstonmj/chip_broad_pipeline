# Load libraries
library(Hmisc)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(corrplot)

# Input data
value_file=snakemake@input[["counts"]]
chip_read_table = read_tsv(value_file) %>%
  unite(bin, c(`#'chr'`,`'start'`,`'end'`)) %>%
  select(-bin) %>%
  na.omit()

# Calculate statistical correlation values
stat_vals = rcorr(as.matrix(chip_read_table), type = "pearson")

# Plot and save
cairo_pdf(filename = snakemake@output[["pdf"]],
          height = 16,
          width = 16,
          family = "sans",
          fallback_resolution = 600)
corrplot(stat_vals$r, method="number", type="upper", order="hclust")
dev.off()

# Export pearson values as TSV
write_tsv(as_tibble(stat_vals$r) %>%
            add_column(id = rownames(stat_vals$r), .before = 1),
          col_names = TRUE, path=snakemake@output[["tsv"]])

# # Export pearson values as TSV
# write_tsv(as_tibble(stat_vals$r, rownames = NA), 
#           col_names = TRUE, path=snakemake@output[["tsv"]])