#!/usr/bin/env Rscript

library(dplyr)
library(vroom)
library(purrr)
library(magrittr)

cohorts <- c(regeneron = "output/merged/BioMe_array-regeneron_GSA_chr17",
             sema4 = "output/merged/BioMe_array-sema4_GDA_chr17")

add_extension <- function(x, nm, ext) setNames(paste0(x, ext), nm)

joins <- function(x, by_cols) {
  out_tab <- x[[1]]
  for (i in seq_along(x)) {
    out_tab <- full_join(out_tab, x[[i]], by = by_cols)
  }
  return(out_tab)
}


fbs <- cohorts %>%
  imap(add_extension, ".fb") %>%
  map(vroom, comment = "#ref", col_types = cols(.default = "d"), na = ".")

msps <- cohorts %>%
  imap(add_extension, ".msp") %>%
  map(vroom, comment = "#Subpop", na = ".",
      col_types = cols(.default = "i", sgpos = "d", egpos = "d"))

msp <- msps %>%
  map(~ rename(.x, chr = "#chm", n = "n snps")) %>%
  joins(c("chr", "spos", "epos", "sgpos", "egpos", "n"))
  
fb <- fbs %>%
  map(~ rename(.x, chr = "chromosome",
      pos = "physical position", gpos = "genetic_position")) %>%
  joins(c("chrom", "pos", "gpos", "genetic_marker_index"))

names_fb <- fbs %>%
  map(~ .x[5:length(.x)]) %>%
  map(names)

names_msp <- msps %>%
  map(~ .x[7:length(.x)]) %>%
  map(names)

rm("msps", "fbs")

save.image("ancestry.rda")