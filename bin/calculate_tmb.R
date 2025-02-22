#!/usr/bin/env Rscript

#' title: "Estimating Tumor Mutational Burden (TMB)"
#' author: "Kelsey Zhu"
#' date: "September 24, 2024"

library(dplyr)
library(optparse)

option_list = list(
  make_option(c("-m", "--maf"), type="character", default=NULL,
              help="maf file"),
  make_option(c("-s", "--hotspot_maf"), type="character", default=NULL,
              help="path to the hotspot variant file"),
  make_option(c("-o", "--out"), type="character", default="tmb.tsv",
              help="output file name [default= %default]"),
  make_option(c("-d", "--depth"), type="integer", default=250,
              help="read depth in tumor"),
  make_option(c("-c", "--count"), type="integer", default=5,
              help="tumor alt count"),
  make_option(c("-v", "--vaf"), type="double", default=0.05,
              help="variant allele frequency in tumor"),
  make_option(c("-s", "--callable"), type="integer", default=1408279,
              help="TMB regions covered")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

calculate_tmb <- function(maf_file, hotspot_maf, out_tsv, tumor_coverage, alt_count, tumor_vaf, callable) {
  hotspot_df <- read.delim(hotspot_maf, header = TRUE, sep = "\t",
                           comment.char = "#", stringsAsFactors = FALSE)
  hotspot_df <- hotspot_df %>% mutate(var_id=paste(Chromosome, Start_Position,
                                            End_Position,Reference_Allele,Tumor_Seq_Allele2,sep="-"))
  tmb_maf <- read.delim(maf_file, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  tmb_df <- tmb_maf %>% 
                    select(Hugo_Symbol, Center, Chromosome, Start_Position, End_Position,
                               Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2,
                               Strand, Variant_Classification, Variant_Type,Tumor_Sample_Barcode,
                               HGVSc, HGVSp, t_depth, t_ref_count, t_alt_count,
                               n_depth, n_ref_count, n_alt_count, Consequence,
                               CLIN_SIG, Existing_variation, IMPACT, AF) %>%
                    filter( grepl('missense|synonymous|frameshift|inframe|stop_gained', Consequence)
                          & t_depth >= tumor_coverage
                          & t_alt_count >= alt_count
                          & AF <= 0.001) %>%
                    mutate(var_id=paste(Chromosome, Start_Position, End_Position,
                            Reference_Allele,Tumor_Seq_Allele2,sep="-"),
                            t_vaf=round(t_alt_count/t_depth,2),
                            n_vaf=round(n_alt_count/n_depth,2),
                            cosmic_count=str_count(Existing_variation, pattern = "COSM")) %>%
                    filter( t_vaf >= tumor_vaf, cosmic_count < 3,
                            !var_id %in% hotspot_df$var_id) %>%
                    arrange(Hugo_Symbol, var_id, Center) 
  # remove duplicates
  tmb_df <- tmb_df[!duplicated(tmb_df$var_id),]
  # calcualte TMB
  tmb <- data.frame(unclass(table(tmb_df$Tumor_Sample_Barcode)))
  t_value <- apply(tmb,1, function(x, size=callable) round(x/round(size/10^6,2),2)  )
  write.table(data.frame(t_value), file=out_tsv, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
}

# Estimate Tumor Mutational Burden
calculate_tmb(opt$maf, opt$hotspot_maf, opt$out, opt$depth, opt$count,
              opt$vaf, opt$callable)