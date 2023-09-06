library(readr)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))
library(purrr)

input_metadata <- "Downloads/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz"
output_table <- "reference_proc/hgdp_1kg.popdata.tsv.gz"
output_smap <- "reference_proc/hgdp_1kg.smap"

input_metadata <- snakemake@input[["metadata"]]
input_samplist <- snakemake@input[["samplist"]]
output_table <- snakemake@output[["table"]]
output_smap <- snakemake@output[["smap"]]


parse_gnomad <- Vectorize(function(x) {
    if (is.na(x)) {
        return(NA)
    } else {
        return(jsonlite::fromJSON(x))
    }
}, SIMPLIFY = F, USE.NAMES = F)

samplist <- input_samplist |>
  read_lines()

poptab <- input_metadata |>
  read_tsv() |>
  select(s, gnomad_population_inference, high_quality) |>
  mutate(popinf = map_vec(gnomad_population_inference, parse_gnomad)) |>
  filter(!is.na(popinf) & !is.na(high_quality), high_quality) |>
  select(-high_quality, -gnomad_population_inference) |>
  unnest_wider(popinf) |>
  rename(pc = pca_scores) |>
  unnest_wider(pc, names_sep = "") |>
  filter(s %in% samplist) |>
  write_tsv(output_table)

poptab |>
  select("#Sample" = "s", "Panel" = "pop") |>
  write_tsv(output_smap)

