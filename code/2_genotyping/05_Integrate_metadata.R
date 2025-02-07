# Make metadata file for TARGET-seq+ multiomic analysis

# This integrates the FACS indexing and genotyping results to make a metadata table that can be used in single cell genomics analysis packages.

## Load the relevant libraries ####
library(tidyverse)


## Read in data ####

# Read in the FACS indexing data
facs_index_data <- read_tsv("data/FACS_indexing/FACS_index_data.tsv")


# Read in the metadata files which includes blank wells
metadata <- read_tsv("data/genotyping/genotyping_metadata.txt")


# Paths to genotyping data files
  
genotyping_filepaths <- list.files(path = "data/genotyping/", 
                                     pattern = "genotypes.*",
                                     full.names = T)

# Read in genotyping files for all samples
genotyping_data <- genotyping_filepaths %>% 
  map_df(~read_tsv(.))


# Make genotype groups for cross sample analysis
genotyping_data <- genotyping_data %>% 
  select(Cell, Sample, Sample_type, Plate, Clone) %>% 
  mutate(Genotype = case_when(
    str_detect(Clone, "DNMT3A") ~ "DNMT3A",
    Sample_type == "CH" & Clone == "WT" ~ "WT-CH",
    Sample_type == "Control" ~ "WT-control",
    TRUE ~ Clone
  ))


## Join all the metadata ####

metadata <- metadata %>% 
  left_join(facs_index_data) %>% 
  left_join(genotyping_data)

metadata %>% write_tsv('data/rna/integrated_metadata.txt')


