#This code allows calculating log2 fold-changes in protein abundance
#for single gene deletion mutants of E. coli

#Load required packages
library(tidyverse)
library(vsn)


#Load sample mapping file
#(contains information of what is in each TMT channel in each sample)
sample_mapping <- read_delim('sample_mapping.csv',delim = ',')


#Make sure that protein files are in MS_data folder
#Data can be downloaded from PRIDE (accession numbers: PXD023945 and PXD016589)
#Only files ending with '_proteins.txt' are necessary
MS_files <- unique(sample_mapping$protein_file)


#Import protein files and perform vsn normalization per sample
prot_melt <- bind_rows(lapply(MS_files, function(MS_file){

  #Get sample number
  sample_number <- unique(filter(sample_mapping, protein_file == MS_file) %>% pull(MS_sample_number))
  
  #Only files from temperatures 1 and 2 from PXD016589 are necessary
  if(sample_number %in% unique(sample_mapping$MS_sample_number)) {
    
    #Import protein file, remove contaminants and reverse database hits and filter for proteins with at least 2 unique peptides
    prot_tab <- read_delim(file.path('MS_data',MS_file), delim = "\t") %>%
      filter(!grepl("##", protein_id),
             qupm > 1) %>%
      group_by(gene_name) %>%
      filter(qupm == max(qupm), top3 == max(top3)) %>%
      ungroup
    
    #Get which channels were used for each experiment
    TMT_labels <- sample_mapping %>%
      filter(MS_sample_number == sample_number) %>%
      pull(TMT_label)
    
    #Perform vsn normalization
    signal_sum_mat <- as.matrix(prot_tab %>%
                                  dplyr::select(paste0("signal_sum_",TMT_labels)))
    
    vsn_fit <- vsn2(signal_sum_mat)
    vsn_norm_mat <- as.data.frame(predict(vsn_fit, signal_sum_mat))
    vsn_norm_mat$gene_name <- prot_tab$gene_name
    
    vsn_norm_mat %>%
      as_tibble() %>%
      gather(key,value,-gene_name) %>%
      mutate(MS_sample_number = sample_number)
  }
  
}))

#Merging data from all experiments and annotating each data point
prot_merged <- prot_melt %>%
  mutate(key = sub('signal_sum_','',key)) %>%
  left_join(sample_mapping,
            by = c('key' = 'TMT_label','MS_sample_number' = 'MS_sample_number')) %>%
  dplyr::select(-key) %>%
  filter(!is.na(KO_mutant)) %>%
  group_by(KO_mutant,gene_name,sample_type) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  group_by(gene_name,MS_sample_number) %>%
  mutate(median_gene = median(value)) %>%
  ungroup() %>%
  mutate(rep_logFC = value - median_gene) %>%
  group_by(KO_mutant,gene_name,sample_type) %>%
  mutate(mean_logFC = mean(rep_logFC)) %>%
  ungroup() %>%
  distinct() %>%
  filter(KO_mutant != gene_name)