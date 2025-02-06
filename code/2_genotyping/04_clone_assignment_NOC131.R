# Clone assignment from TARGET-seq single cell genotyping data

# For integrating single cell genotypes across multiple mutations amplified within the same cells to call clone identities.

# Sample NOC131

#### Load the relevant libraries ####
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(forcats)

#### Set key variables ####

# Paths to genotype summary files
mutation_1_path <- "data/genotyping/sample_NOC131/genotyping_summary_NOC131_DNMT3A_pQ606X.tsv"
mutation_2_path <- "data/genotyping/sample_NOC131/genotyping_summary_NOC131_DNMT3A_pI780T.tsv"

# Name of control sample
control <- "NOC153"

# Sample being analysed
sample <- "NOC131"

## Read in files #####

# Define a function that reads file and takes columns needed
read_genotyping_file <- function(file_path) {
  file <- read_tsv(file_path) %>% 
    select(Cell, Sample, Sample_type, Plate, mutation, genotype)
  mutation_name <- file$mutation[1]
  file <- file %>% rename(!!mutation_name := genotype) %>% 
    select(-mutation)
  return(file)
}


# Read in the genotyping data
mutation_1 <- read_genotyping_file(mutation_1_path)
mutation_2 <- read_genotyping_file(mutation_2_path)

# Join them together in a single dataframe
genotypes <- left_join(mutation_1, mutation_2)


## Visualise the mutation calls #### 

# We make a raster plot to visualise the mutation calls in each single cell. This gives us an idea of the clonal structure.

# First remove cells from the control sample and sort the cells by the genotype of each mutation
# Then make Cell ID a factor with this order - this makes the raster neater
genotyping_plot <- genotypes %>% 
  filter(Sample ==sample) %>% 
  mutate(DNMT3A_pQ606X = factor(DNMT3A_pQ606X, levels = c("Mut", "WT", "ADO", "Undetermined","Undetected")),
         DNMT3A_pI780T = factor(DNMT3A_pI780T, levels = c("Mut", "WT", "ADO", "Undetermined","Undetected")),
         ) %>% 
  arrange(DNMT3A_pQ606X, DNMT3A_pI780T) %>% 
  mutate(Cell=factor(Cell, levels=unique(Cell)))

# Then make a long format table
genotyping_plot_long <- genotyping_plot %>% 
  pivot_longer(c(DNMT3A_pQ606X, DNMT3A_pI780T), names_to = "Mutation", values_to  = "genotype") %>% 
  # Call ADO cells undetermined
  mutate(genotype = if_else(genotype =="ADO", "Undetermined",genotype)) %>% 
  mutate(genotype = factor(genotype, levels = c("Mut", "WT", "ADO", "Undetermined","Undetected")))

genotyping_plot_long %>% 
  mutate(Mutation = str_replace(Mutation, "_", " ")) %>% 
  ggplot(aes(x =  Cell, y = Mutation, fill=genotype))+
  geom_tile()+
  scale_fill_manual(values = c("Undetermined" = "black","Undetected" = "#f0f0f0", "WT" = "#d9d9d9", "Mut" = "#e41a1c"))+ 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "left", axis.title.x = element_text(hjust = 0))+
  scale_y_discrete(position = "right")+scale_x_discrete(position = "top")+
  labs(title = sample, x = paste("Single cells", sprintf('\u2192')), fill = NULL)

# In this plot, each column represents a single cell.
# From this, we can see that the DNMT3A_pQ606X and DNMT3A_pI780T mutations are mutually exclusive, so form a branched clonal structure.

# The most likely clonal structure can be identified by running infSCITE.


## Assign the clonal identity ####

# We then assign a consensus genotype based on the clonal structure identified from the pattern of mutation co-occurence from the raster and infSCITE

# In this case, the two mutations occur in independent clones
genotypes <- genotypes %>% 
  mutate(Clone = case_when(
    # Cells with either mutation are DNMT3A mut
    DNMT3A_pQ606X %in% c("WT","Undetected", "ADO","Undetermined") & DNMT3A_pI780T == "Mut" ~ "DNMT3A-I780T",
    DNMT3A_pQ606X == "Mut" & DNMT3A_pI780T %in% c("WT","Undetected","Undetermined") ~ "DNMT3A-Q606X",
    # Cells with indeterminate genotype for either mutation are Undetermined
    DNMT3A_pQ606X %in% c("Undetected", "ADO","Undetermined")   & DNMT3A_pI780T %in% c("WT","Undetected","Undetermined") ~ "Undetermined",
    DNMT3A_pQ606X  %in% c("WT","Undetected", "ADO","Undetermined") & DNMT3A_pI780T %in% c("Undetermined","Undetected")  ~ "Undetermined",
    # All other cells are WT
    TRUE ~ "WT"))

genotypes %>% 
  group_by(Sample, Clone) %>% count



## Plot the mutation calls with the assigned clonal identity #### 

# We make a raster plot to visualise the mutation calls and clonal identity in each single cell.

# First remove cells from the control sample and sort the cells by the genotype of each mutation
# Then make Cell ID a factor with this order - this makes the raster neater
genotyping_plot <- genotypes %>% 
  filter(Sample ==sample) %>% 
  mutate(DNMT3A_pQ606X = factor(DNMT3A_pQ606X, levels = c("Mut", "WT", "ADO", "Undetermined","Undetected")),
         DNMT3A_pI780T = factor(DNMT3A_pI780T, levels = c("Mut", "WT", "ADO", "Undetermined","Undetected")),
         Clone = factor(Clone, levels = c("Undetermined", "WT", "DNMT3A-I780T", "DNMT3A-Q606X"))) %>% 
  arrange(Clone, DNMT3A_pQ606X, DNMT3A_pI780T) %>% 
  mutate(Cell=factor(Cell, levels=unique(Cell)))


# Then tidy the table
genotyping_plot_long <- genotyping_plot %>% 
  pivot_longer(c(DNMT3A_pQ606X, DNMT3A_pI780T), names_to = "Mutation", values_to  = "genotype") %>% 
  # Call ADO cells undetermined
  mutate(genotype = if_else(genotype =="ADO", "Undetermined",genotype)) %>% 
  mutate(genotype = factor(genotype, levels = c("Mut", "WT", "ADO", "Undetermined","Undetected")))

cells <- genotyping_plot %>% filter(Sample ==sample) %>% count()

p1 <- genotyping_plot_long %>% 
  mutate(Mutation = str_replace(Mutation, "_", " ")) %>% 
  ggplot(aes(x =  Cell, y = Mutation, fill=genotype))+
  geom_tile()+
  scale_fill_manual(values = c("Undetermined" = "black","Undetected" = "#f0f0f0", "WT" = "#d9d9d9", "Mut" = "#e41a1c"))+ 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "left",
        axis.title.x = element_text(hjust = 0))+
  scale_y_discrete(position = "right")+scale_x_discrete(position = "top")+
  labs(title = sample, x = paste("Single cells", sprintf('\u2192')), fill = NULL)
p1

p2 <- genotyping_plot_long %>% 
  ggplot(aes(x =  Cell, y = Sample, fill= Clone))+
  geom_tile()+
  scale_fill_manual(values = c("Undetermined" = "black", "WT" = "#d9d9d9", "DNMT3A-Q606X" = "#4eb3d3", "DNMT3A-I780T" = "#2c7fb8"))+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), 
        legend.position = "left",
        axis.line = element_blank())+
  scale_y_discrete(position = "right")+
  labs(title = NULL, x = paste("n =", cells, "cells"), y = "Clone", fill = NULL)
p2

pcol <- plot_grid(p1, p2, 
                  align = "v", axis = "b", 
                  nrow = 2,
                  rel_heights = c(2, 1.2))
pcol

## Write clonal assignments to a file for use in multiomic analysis ####

genotypes %>% write_tsv(paste("data/genotyping/genotypes_", sample,".txt", sep = ""))

