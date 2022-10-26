library(tidyverse)
library(glue)

awk_filter_fb <- function(fb, start, end, ancestry = "AFR") {
  "awk 'NR == 2 || ($2 >= {start} && $2 <= {end})' {fb}" %>%
    glue() %>%
    system(intern = TRUE) %>%
    I() %>%
    read_tsv(col_types = cols(.default = "d")) %>%
    select(chr = chromosome, pos =  "physical position", ends_with(ancestry))
}

awk_filter_msp <- function(msp, start, end) {
  "awk 'NR == 2 || ($2 >= {start} && $3 <= {end})' {msp}" %>%
    glue() %>%
    system(intern = TRUE) %>%
    I() %>%
    read_tsv(col_types = cols(.default = "d")) %>%
    select(chr = "#chm", spos, epos, matches("\\.[10]$"))
}

pivot_msp <- function(msp) {
  pvt <- function(df, suffix) {
    df %>%
      select(chr, spos, epos, ends_with(suffix)) %>%
      pivot_longer(
        cols = ends_with(suffix),
        names_to = "sample",
        values_to = "ancestry") %>%
      mutate(sample = str_replace(sample, paste0("\\", suffix), "")) %>%
      rename_with(~ paste0(.x, suffix), "ancestry")
  }

  pvt(msp, ".0") %>%
    left_join(pvt(msp, ".1"), by = c("chr", "spos", "epos", "sample"))
}

get_hom_afr <- . %>%
  group_by(sample) %>%
  summarize(hom_afr_all = all(hom_afr)) %>%
  filter(hom_afr_all) %>%
  pull(sample)

start <- 43384997
end <- 44913630

## Regeneron

regeneron_gsa_afr <- "output/merged/BioMe_array-regeneron_GSA_chr17.fb" %>%
  awk_filter_fb(start, end)

regeneron_gsa_msp <- "output/merged/BioMe_array-regeneron_GSA_chr17.msp" %>%
  awk_filter_msp(start, end)

regeneron_gsa_hc <- regeneron_gsa_msp %>%
  pivot_msp() %>%
  mutate(hom_afr = (ancestry.0 == 3 & ancestry.1 == 3),
         het_afr = ((ancestry.0 == 3 | ancestry.1 == 3) & !hom_afr),
         sample_single = str_sub(sample, end = (nchar(sample) - 1) / 2),
         sample_reconstruct = paste(sample_single, sample_single, sep = "_"))

stopifnot(nrow(filter(regeneron_gsa_hc, sample_reconstruct != sample)) == 0)

regeneron_gsa_hc <- regeneron_gsa_hc %>%
  select(-sample, -sample_reconstruct, sample = sample_single)

regeneron_gsa_hapsum <- regeneron_gsa_hc %>%
  group_by(sample) %>%
  summarise(hom_afr = sum(hom_afr), het_afr = sum(het_afr))

regeneron_gsa_hapsum %>%
  count(hom_afr, het_afr) %>%
  write_excel_csv("output/merged/rg_counts.csv")

regeneron_gsa_msp %>%
  select(chr, spos, epos) %>%
  write_excel_csv("output/merged/rg_hardcall_ranges.csv")

regeneron_gsa_AFR %>%
  select(chr, pos) %>%
  write_excel_csv("output/merged/rg_prob_positions.csv")

regeneron_gsa_hc %>%
  get_hom_afr %>%
  write_lines("output/merged/rg_afrhom_samples.txt")

## Sema4

sema4_gda_afr <- "output/merged/BioMe_array-sema4_GDA_chr17.fb" %>%
  awk_filter_fb(start, end)

sema4_gda_msp <- "output/merged/BioMe_array-sema4_GDA_chr17.msp" %>%
  awk_filter_msp(start, end)

sema4_gda_hc <- sema4_gda_msp %>%
  pivot_msp() %>%
  mutate(hom_afr = (ancestry.0 == 3 & ancestry.1 == 3),
         het_afr = ((ancestry.0 == 3 | ancestry.1 == 3) & !hom_afr),
         sample_single = str_sub(sample, end = (nchar(sample) - 1) / 2),
         sample_reconstruct = paste(sample_single, sample_single, sep = "_"))

stopifnot(nrow(filter(sema4_gda_hc, sample_reconstruct != sample)) == 0)

sema4_gda_hc <- sema4_gda_hc %>%
  select(-sample, -sample_reconstruct, sample = sample_single)

sema4_gda_hapsum <- sema4_gda_hc %>%
  group_by(sample) %>%
  summarise(hom_afr = sum(hom_afr), het_afr = sum(het_afr))

sema4_gda_hapsum %>%
  count(hom_afr, het_afr) %>%
  write_excel_csv("output/merged/s4_counts.csv")

sema4_gda_msp %>%
  select(chr, spos, epos) %>%
  write_excel_csv("output/merged/s4_hardcall_ranges.csv")

sema4_gda_AFR %>%
  select(chr, pos) %>%
  write_excel_csv("output/merged/s4_prob_positions.csv")

sema4_gda_hc %>%
  get_hom_afr %>%
  write_lines("output/merged/s4_afrhom_samples.txt")
