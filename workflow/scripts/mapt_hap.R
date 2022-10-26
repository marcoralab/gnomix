#!/usr/bin/env Rscript
library(stringr)
library(readr)
library(tidyr)
library(tibble)
library(purrr)
suppressPackageStartupMessages(library(dplyr))
library("writexl")

calc_hap <- function(rs1052553 = NA, rs8070723 = NA) {
  count_allele <- . %>%
    str_extract_all("[0-2]") %>%
    map_dbl(~ sum(as.numeric(.x)))
  if (all(is.na(rs1052553)) && all(is.na(rs8070723))) {
    stop("rs1052553 and rs8070723 not genotyped")
  } else if (all(is.na(rs1052553))) {
    warning("rs1052553 not genotyped. Using only rs8070723")
  } else if (all(is.na(rs8070723))) {
    warning("rs8070723 not genotyped. Using only rs1052553")
  }

  tibble(rs1 = rs1052553, rs8 = rs8070723) %>%
    mutate_all(count_allele) %>%
    mutate_all(~ case_when(
      .x == 0 ~ 11,
      .x == 1 ~ 12,
      .x == 2 ~ 22,
      TRUE ~ NA_real_)) %>%
    mutate(allele = case_when(
      is.na(rs1) ~ rs8,
      is.na(rs8) ~ rs1,
      rs1 == rs8 ~ rs1,
      TRUE ~ NA_real_)) %>%
    pull(allele)
}

add_hap <- function(genotypes, alleles, colname) {
  # Add haplotype column to genotype df using calc_hap

  # Default value for missing columns is NA
  default <- alleles %>%
    select("variant") %>%
    mutate(def = NA_character_) %>%
    tibble::deframe()

  # calculate haplotype using all genotype columns
  calc <- . %>%
    select(-VCF_ID) %>%
    as.list %>%
    rlang::exec(calc_hap, !!!.)

  genotypes %>%
    #Add missing columns
    add_column(!!!default[setdiff(names(default), names(.))]) %>%
    mutate(hap = calc(.)) %>%
    rename(!!colname := hap)
}

get_hap <- function(vcf, alleles, hapname) {
  cpra <- alleles %>%
    select("variant", "cpra") %>%
    tibble::deframe() %>%
    as.list()

  rsids <- alleles %>% pull(variant)

  pos_bcftools <- paste(alleles$chrom_pos, collapse = ",")
  query_format <- "%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE=%GT, ]\n"
  bcftools <- sprintf("bcftools query -r '%s' -f '%s' %s",
                      pos_bcftools, query_format, vcf)

  system(command = bcftools, intern = TRUE) %>%
    str_replace(", $", "") %>%
    I() %>%
    read_tsv(col_names = c("chromsome", "position", "ref", "alt", "genotype"),
             col_types = "ciccc") %>%
    separate_rows(genotype, sep = ", ") %>%
    separate(genotype, sep = "=", c("VCF_ID", "genotype")) %>%
    #combine chromsome, position, ref, and alt into a new ID
    unite("CPRA", chromsome, position, ref, alt, sep = ":", remove = TRUE) %>%
    pivot_wider(names_from = CPRA, values_from = genotype) %>%
    rename(any_of(unlist(cpra))) %>%
    add_hap(alleles, hapname) %>%
    select(-all_of(rsids))
}

build <- "hg19"

positions <- list(
    hg19 = c(rs1052553 = "17:44073889", rs8070723 = "17:44081064"),
    hg38 = c(rs1052553 = "chr17:45996523", rs8070723 = "chr17:46003698")) %>%
  .[[build]]

alleles <- tribble(~variant, ~ref, ~alt,
                   "rs1052553", "A", "G",
                   "rs8070723", "A", "G") %>%
  mutate(chrom_pos = positions[variant]) %>%
  unite(cpra, chrom_pos, ref, alt, sep = ":", remove = FALSE)

mapt_process <- function(vcfname, cohortname) {
  message(paste("Calculating", cohortname))
  genos <- vcfname %>%
    get_hap(alleles, "MAPT") %>%
    mutate(mapt2dose = stringr::str_count(as.character(MAPT), "2"),
           sample = str_sub(VCF_ID, end = (nchar(VCF_ID) - 1) / 2),
           sample_reconstruct = paste(sample, sample, sep = "_"),
           cohort = cohortname)

  stopifnot(nrow(filter(genos, sample_reconstruct != VCF_ID)) == 0)

  genos %>%
    select(-VCF_ID, -sample_reconstruct)
}

mapt_s4 <- "output/BioMe_array-sema4_GDA_chrall.vcf.gz" %>%
  mapt_process("Sema4") %>%
  mutate(afr = sample %in% read_lines("output/merged/s4_afrhom_samples.txt"))

mapt_rg <- "output/BioMe_array-regeneron_GSA_chrall.vcf.gz" %>%
  mapt_process("Regeneron") %>%
  mutate(afr = sample %in% read_lines("output/merged/rg_afrhom_samples.txt"))

bind_rows(mapt_s4, mapt_rg) %>%
  filter(afr & mapt2dose == 2) %>%
  list(Samples = bind_rows(mapt_s4, mapt_rg),
       Passes = .) %>%
  write_xlsx("output/merged/h2hom_afr_samples.xlsx")

get_hap_rs8070723 <- function(vcf, alleles, hapname) {
  cpra <- alleles %>%
    select("variant", "cpra") %>%
    tibble::deframe() %>%
    as.list()

  rsids <- alleles %>% pull(variant)

  pos_bcftools <- paste(alleles$chrom_pos, collapse = ",")
  query_format <- "%CHROM\t%POS\t%REF\t%ALT\t[%SAMPLE=%GT, ]\n"
  bcftools <- sprintf("bcftools query -r '%s' -f '%s' %s",
                      pos_bcftools, query_format, vcf)

  system(command = bcftools, intern = TRUE) %>%
    str_replace(", $", "") %>%
    I() %>%
    read_tsv(col_names = c("chromsome", "position", "ref", "alt", "genotype"),
             col_types = "ciccc") %>%
    separate_rows(genotype, sep = ", ") %>%
    separate(genotype, sep = "=", c("VCF_ID", "genotype")) %>%
    #combine chromsome, position, ref, and alt into a new ID
    unite("CPRA", chromsome, position, ref, alt, sep = ":", remove = TRUE) %>%
    pivot_wider(names_from = CPRA, values_from = genotype) %>%
    rename(any_of(unlist(cpra))) %>%
    mutate(rs1052553 = NA_character_) %>%
    add_hap(alleles, hapname) %>%
    select(-all_of(rsids))
}

"output/BioMe_array-sema4_GDA_chrall.vcf.gz" %>%
  get_hap_rs8070723(alleles, "MAPT") %>%
  mutate(mapt2dose = stringr::str_count(as.character(MAPT), "2"),
         sample = str_sub(VCF_ID, end = (nchar(VCF_ID) - 1) / 2),
         sample_reconstruct = paste(sample, sample, sep = "_"),
         afr = sample %in% read_lines("output/merged/s4_afrhom_samples.txt")) %>%
  count(MAPT, afr)
