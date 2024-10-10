#load libraries

library(tidyverse)

#library(xlsx)

####CUSTOM FUNCTION####

#function to read snpEff data and merge snps and indels of the same sample

read_snpEff = function(sample_name){
  
  header = c("CHROM","POS","REF","ALT","GT","AD","ADF","ADR",
             "GQ","FREQ","PVAL","EFF")
  
  snps = read.table(paste0(i,"_snps_snpEff_annotation.tsv"),
                    sep = "\t", 
                    header = FALSE)
  colnames(snps) = header
  snps$VC = "SNV"
  snps$SAMPLE_ID = i
  
  indels = read.table(paste0(i,"_indels_snpEff_annotation.tsv"),
                      sep = "\t", 
                      header = FALSE)
  colnames(indels) = header
  indels$VC = "INDEL"
  indels$SAMPLE_ID = i
  
  assign(gsub("-","_",paste0("sample",i)), 
         rbind(snps, indels), 
         .GlobalEnv)
  
}

####MAIN SCRIPT####

#Read files

sample_id = unique(gsub("_indels_snpEff_annotation.tsv|_snps_snpEff_annotation.tsv",
                        "",
                        list.files(pattern = ".tsv")))

for(i in sample_id){
  read_snpEff(i)
}
rm(i)

#unify sample files

sample_id = paste0("sample",gsub("-","_",sample_id))

main_dataset = eval(parse(
  text = paste0("rbind(",paste0(sample_id, collapse = ","),")")))

eval(parse(
  text = paste0("rm(",paste0(sample_id, collapse = ","),")")))

rm(sample_id)

#Adding gene name to the main_dataset

snpeff_annotation = main_dataset %>%
  select(EFF) %>%
  rowwise() %>%
  mutate(EFF = gsub("\\(", "|", EFF)) %>%
  mutate(EFF = gsub("\\)","", EFF)) %>%
  separate(., 
           col = EFF, 
           into = paste0("V", rep(1:10)), 
           sep = "\\|")

main_dataset = cbind(main_dataset, snpeff_annotation)
