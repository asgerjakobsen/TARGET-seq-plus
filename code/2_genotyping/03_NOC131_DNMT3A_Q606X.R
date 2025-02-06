# Analysis of TARGET-seq single cell genotyping data

# For calling single cell genotypes at a single mutation locus, integrating data from gDNA and mRNA/cDNA amplicons.

# Sample NOC131 DNMT3A_Q606X
# In this case, a heterozygous germline SNP within the gDNA amplicon is used for allelic dropout (ADO) determination.

#### Load the relevant libraries ####
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

#### Set key variables ####
# This makes it easier to use code for other files

# Paths to genotype call files
gDNA_filepath <- "data/genotyping/sample_NOC131/gDNA_DNMT3A_pQ606X.txt"
mRNA_filepath <- "data/genotyping/sample_NOC131/mRNA_DNMT3A_pQ606X.txt"
gDNA_SNP_filepath <- "data/genotyping/sample_NOC131/gDNA_DNMT3A_rs2289093.txt"

# Path to metadata file
metadata_filepath <- "data/genotyping/genotyping_metadata.txt"

# Name of control sample
control <- "NOC153"

# Test sample being analysed
sample <- "NOC131"

# Mutation base change
# This sets which columns in the file to use as reference and variant alleles
ref <- quo(G)
mut <- quo(A)
SNPref <- quo(T)
SNPalt <- quo(G)


## READ IN FILES #####

# Read in the csv file for each sequencing run:
genotyping_gDNA <- read_tsv(gDNA_filepath)
genotyping_mRNA <- read_tsv(mRNA_filepath)

# Read in SNP files
genotyping_SNP <- read_tsv(gDNA_SNP_filepath)

# Read in the metadata file
metadata <- read_tsv(metadata_filepath)


# Separate the id column to get the cell ID and tidy up
genotyping_gDNA <- genotyping_gDNA %>% 
  separate(id, c("id", "mutation"), sep = "[.]") %>% 
  separate(id, c("Cell", NA), sep = "[_]") %>% 
  rename(chromosome = bp, position = ref)
str(genotyping_gDNA)
genotyping_SNP <- genotyping_SNP %>% 
  separate(id, c("id", "mutation"), sep = "[.]") %>% 
  separate(id, c("Cell", NA), sep = "[_]") %>% 
  rename(chromosome = bp, position = ref)
str(genotyping_SNP)

genotyping_mRNA<- genotyping_mRNA %>% 
  separate(id, c("id", "mutation"), sep = "[.]") %>% 
  separate(id, c("Cell", NA), sep = "[_]") %>% 
  rename(chromosome = bp, position = ref)
str(genotyping_mRNA)


# Mutation name
mutation_name <- genotyping_gDNA$mutation[1]
SNP_name <- genotyping_SNP$mutation[1]


### JOIN METADATA TO VARIANT CALLS ####

# Join the library metadata to the genotyping data]
genotyping_gDNA <- genotyping_gDNA %>% 
  full_join(metadata, by = "Cell")
genotyping_SNP <- genotyping_SNP %>% 
  full_join(metadata, by = "Cell")
genotyping_mRNA <- genotyping_mRNA %>% 
  full_join(metadata, by = "Cell") 


# List of plates for the relevant sample
sample_plates <- genotyping_gDNA %>% 
  filter(Sample == sample)  %>% 
  pull(Plate) %>% 
  unique()

# Keep only the plates relevant to this mutation
genotyping_gDNA <- genotyping_gDNA %>% 
  filter(Plate %in% sample_plates)
genotyping_SNP <- genotyping_SNP %>% 
  filter(Plate %in% sample_plates)
genotyping_mRNA <- genotyping_mRNA %>% 
  filter(Plate %in% sample_plates)


# Where no reads were detected and no fastq file exists, the sample will not be in the genotyping file
# Fill in the read counts with zeros for these, to allow them to be included, with no coverage

genotyping_gDNA <- genotyping_gDNA %>% 
  replace_na(list(A = 0, G = 0, C = 0, T = 0, del = 0, ins = 0))
genotyping_SNP <- genotyping_SNP %>% 
  replace_na(list(A = 0, G = 0, C = 0, T = 0, del = 0, ins = 0))
genotyping_mRNA <- genotyping_mRNA %>% 
  replace_na(list(A = 0, G = 0, C = 0, T = 0, del = 0, ins = 0))

# Calculate the coverage and mutant VAF 
genotyping_gDNA <- genotyping_gDNA %>% 
  mutate(coverage = A + G + C + T + del + ins) %>% 
  mutate(VAF = !!mut / coverage)
genotyping_SNP <- genotyping_SNP %>% 
  mutate(coverage = A + G + C + T + del + ins) %>% 
  mutate(VAF = !!SNPalt / coverage)
genotyping_mRNA <- genotyping_mRNA %>% 
  mutate(coverage = A + G + C + T + del + ins) %>% 
  mutate(VAF = !!mut / coverage)



## FILTERING #####


# Calculate the coverage in blank wells to set the threshold for detection
genotyping_gDNA_blanks <- genotyping_gDNA %>% 
  filter(Sample == "Blank")

gDNA_blank_cov <- genotyping_gDNA_blanks %>% 
  group_by(Plate) %>% 
  summarise(mean_coverage = mean(coverage),
            max_coverage = max(coverage),
            min_coverage = min(coverage),
  ) %>% arrange(desc(max_coverage))
gDNA_blank_cov


genotyping_SNP_blanks <- genotyping_SNP %>% 
  filter(Sample == "Blank")
SNP_blank_cov <- genotyping_SNP_blanks %>% 
  group_by(Plate) %>% 
  summarise(mean_coverage = mean(coverage),
            max_coverage = max(coverage),
            min_coverage = min(coverage),
  ) %>% 
  arrange(desc(max_coverage))
SNP_blank_cov

genotyping_mRNA_blanks <- genotyping_mRNA %>% 
  filter(Sample == "Blank")

mRNA_blank_cov <- genotyping_mRNA_blanks %>% 
  group_by(Plate) %>% 
  summarise(mean_coverage = mean(coverage),
            max_coverage = max(coverage),
            min_coverage = min(coverage),
  ) %>% 
  arrange(desc(max_coverage))
mRNA_blank_cov

# Take a look at the coverage for non-blank samples
genotyping_gDNA %>% 
  filter(Sample != "Blank") %>% 
  select(Cell, A, G, C, T, coverage, VAF) %>% 
  arrange(desc(coverage))

genotyping_SNP %>% 
  filter(Sample != "Blank") %>% 
  select(Cell, A, G, C, T, coverage, VAF) %>% 
  arrange(desc(coverage))

genotyping_mRNA %>% 
  filter(Sample != "Blank") %>%   
  select(Cell, A, G, C, T, coverage, VAF) %>% 
  arrange(desc(coverage))

# Plot the coverage
# gDNA
ggplot(genotyping_gDNA %>% mutate(Sample=if_else(Sample == "NOC153", "Control", Sample)),
       aes(x= Plate, y=(coverage+1), fill = Sample, colour = Sample)) +
  geom_violin(scale = "width", alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.9), alpha = 0.2)+
  labs(title = paste(mutation_name, "gDNA amplicon"), y = "Coverage")+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# gDNA SNP
ggplot(genotyping_SNP %>% mutate(Sample=if_else(Sample == "NOC153", "Control", Sample)),
       aes(x= Plate, y=(coverage+1), fill = Sample, colour = Sample)) +
  geom_violin(scale = "width", alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.9), alpha = 0.2)+
  labs(title = paste(SNP_name, "gDNA amplicon"), y = "Coverage")+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# mRNA
ggplot(genotyping_mRNA %>% mutate(Sample=if_else(Sample == "NOC153", "Control", Sample)), 
       aes(x= Plate, y=(coverage+1), fill = Sample, colour = Sample)) +
  geom_violin(scale = "width", alpha = 0.4) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.9), alpha = 0.2)+
  labs(title = paste(mutation_name, "mRNA amplicon"), y = "Coverage")+
  scale_y_log10()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




##### Join the mutation data to SNP calls ####

genotyping_gDNA <- genotyping_SNP %>% 
  select(Cell, mutation, !!SNPref, !!SNPalt, coverage, VAF) %>% 
  rename(SNP = mutation,
         SNP.ref = !!SNPref, 
         SNP.alt = !!SNPalt, 
         SNP.coverage = coverage, 
         SNP.VAF = VAF) %>% 
  full_join(genotyping_gDNA)

# Plot the coverage of SNP vs mutation - this should be well correlated
ggplot(genotyping_gDNA %>% mutate(Sample=if_else(Sample == "NOC153", "Control", Sample)),
       aes(x= coverage, y=(SNP.coverage), colour = Sample)) +
  geom_point()+
  labs(title = paste(mutation_name, "gDNA amplicon"), x = "Mutation coverage", y = "SNP coverage")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



##### Filter the dataset to remove failed libraries ####

# We set the minimum coverage for gDNA and mRNA amplicons
gDNA_min_coverage <- 50
mRNA_min_coverage <- 30

# Filter out cells with coverage below the threshold
genotyping_gDNA_filtered <- genotyping_gDNA %>% 
  filter(coverage >= gDNA_min_coverage)
genotyping_mRNA_filtered <- genotyping_mRNA %>% 
  filter(coverage >= mRNA_min_coverage)


###### Check how many libraries pass QC ####

qc_results <- genotyping_gDNA %>% 
  group_by(Sample) %>% 
  summarise(Total = n(),
            Passed = sum(coverage >= gDNA_min_coverage),
            Failed = sum(coverage < gDNA_min_coverage),
            Pass_rate = Passed/Total,
            Mean_coverage = mean(coverage)) %>% 
  mutate(Amplicon = "gDNA")

qc_results <- genotyping_mRNA %>% 
  group_by(Sample) %>% 
  summarise(Total = n(),
            Passed = sum(coverage >= mRNA_min_coverage),
            Failed = sum(coverage < mRNA_min_coverage),
            Pass_rate = Passed/Total,
            Mean_coverage = mean(coverage)) %>% 
  mutate(Amplicon = "mRNA") %>% 
  bind_rows(qc_results) %>% 
  filter(Sample != "Blank") %>% 
  mutate(Mutation = mutation_name)
qc_results

ggplot(qc_results, aes(Amplicon, Pass_rate, fill = Sample))+
  geom_col(position = "dodge")+
  labs(title = mutation_name, y = "Cells passing QC %")+
  scale_fill_brewer()




## CALLING MUTANT CELLS #####

# Plot mutant VAF

ggplot(genotyping_gDNA_filtered, aes(x = Sample, y = VAF, fill = Sample)) +
  geom_point(position = "jitter")+
  labs(title = paste(mutation_name, "gDNA amplicon"), y = "VAF")

ggplot(genotyping_mRNA_filtered, aes(x = Sample, y = VAF, fill = Sample)) +
  geom_point(position = "jitter") +
  labs(title = paste(mutation_name, "mRNA amplicon"), y = "VAF")


# SNP
ggplot(genotyping_gDNA_filtered, aes(x = Sample, y = SNP.VAF, fill = Sample)) +
  geom_point(position = "jitter")+
  labs(title = paste(SNP_name, "gDNA amplicon"), y = "VAF")

# Examine SNP VAF vs gDNA mutation VAF to determine phasing
# Here the SNP Alt allele is out of phase with the mutation
ggplot(genotyping_gDNA_filtered, aes(x = SNP.VAF, y = VAF, colour = Sample)) +
  geom_point(alpha = 0.3)+
  scale_color_brewer(palette = "Set1")+
  labs(title = paste(mutation_name, "gDNA vs ", SNP_name), x = "Mutation VAF", y = "SNP VAF")



# To define the appropriate VAF thresholds for mutation calling, we use the VAF distribution from WT control cells to determine the noise at each locus.

# We define two thresholds for calling mutations:
# 1. The lower VAF threshold for calling a cell WT is defined as: mean(VAF) + 3 * SD(VAF)
# 2. The upper VAF threshold for calling a cell mutant is defined as: mean(VAF) + 3 * SD(VAF) + 0.01
# Cells between these thresholds are called "Undetermined"

# For the SNP we define two thresholds:
# 1. The VAF threshold for detection of the SNP Alt allele (SNP_HomRef_threshold) is defined as: mean(VAF) + 3 * SD(VAF) + 0.01
# 2. The VAF threshold for calling a cell homozygous alternate (Hom Alt) for the SNP is defined as the inverse of the Alt threshold: 1 - SNP_HomRef_threshold
# Cells with scVAF above this threshold are called Hom Alt (i.e. ADO of the SNP Ref allele had occurred).
# Cells with VAF between the two thresholds are called heterozygous (Het; i.e. biallelic detection).

# Calculate the mean and SD of the WT control VAF distribution and use these to calculate the VAF thresholds
gDNA_control_VAF <- genotyping_gDNA_filtered %>% 
  filter(Sample == control) %>% 
  summarise(mean_VAF = mean(VAF), SD_VAF = sd(VAF),
            max_VAF = max(VAF),
            min_VAF = min(VAF),
            WT_threshold = mean_VAF + 3*SD_VAF,
            MUT_threshold = WT_threshold + 0.01,
            # Stats for the SNP
            mean_SNP_VAF = mean(SNP.VAF),
            SD_VAF_SNP = sd(SNP.VAF),
            max_SNP_VAF = max(SNP.VAF),
            SNP_HomRef_threshold = mean_SNP_VAF+ 3*SD_VAF_SNP + 0.01,
            SNP_HomAlt_threshold = 1 - SNP_HomRef_threshold
  )
gDNA_control_VAF

mRNA_control_VAF <- genotyping_mRNA_filtered %>% 
  filter(Sample == control) %>% 
  summarise(mean_VAF = mean(VAF), SD_VAF = sd(VAF),
            max_VAF = max(VAF),
            min_VAF = min(VAF),
            WT_threshold = mean_VAF + 3*SD_VAF,
            MUT_threshold = WT_threshold + 0.01)
mRNA_control_VAF


# Furthermore, we required a minimum number of 10 mutant reads for a cell to be called mutant. 
gDNA_min_mut_reads <- 10
mRNA_min_mut_reads <- 10


# Plot these thresholds on the VAF distribution
ggplot(genotyping_gDNA_filtered, aes(x = fct_recode(Sample, "Control" ="NOC153"), y = VAF)) +
  geom_jitter(height = 0, size = 0.75)+
  geom_hline(aes(yintercept = gDNA_control_VAF$WT_threshold), linetype = "longdash", color = "#E69F00") +
  geom_hline(aes(yintercept = gDNA_control_VAF$MUT_threshold), linetype = "longdash", color = '#56B4E9') +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = paste(str_replace(mutation_name, "_", " "), "gDNA"), y = "VAF", x = 'Sample')

ggplot(genotyping_gDNA_filtered, aes(x = fct_recode(Sample, "Control" ="NOC153"), y = SNP.VAF)) +
  geom_jitter(height = 0, size = 0.75)+
  geom_hline(aes(yintercept = gDNA_control_VAF$SNP_HomRef_threshold), linetype = "longdash", color = "#E69F00") +
  geom_hline(aes(yintercept = gDNA_control_VAF$SNP_HomAlt_threshold), linetype = "longdash", color = '#56B4E9') +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = paste(str_replace(SNP_name, "_", " "), "gDNA"), y = "VAF", x = 'Sample')

ggplot(genotyping_mRNA_filtered, aes(x = fct_recode(Sample, "Control" ="NOC153"), y = VAF)) +
  geom_jitter(height = 0, size = 0.75)+
  geom_hline(aes(yintercept = mRNA_control_VAF$WT_threshold), linetype = "longdash", color = "#E69F00") +
  geom_hline(aes(yintercept = mRNA_control_VAF$MUT_threshold), linetype = "longdash", color = '#56B4E9') +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(title = paste(str_replace(mutation_name, "_", " "), "mRNA"), y = "VAF", x = 'Sample')


##### Assign genotypes based on VAF thresholds ####

# We then use the SNP genotype calls and the phasing information to establish whether the allele harboring the mutation was detected:
# - When the mutation VAF is above the upper threshold, the cell is called ‘mutant’ and the SNP analysis is not required.

# For cells that appear WT for the mutation, we can use the SNP genotype to determine if they are truly WT or if there was allelic dropout of the mutant allele.
# The criteria depend on whether the SNP Alt allele is in-phase or out-of-phase with the mutation.

# When the SNP Alt allele and the mutation are in-phase:
# - Cells that appear WT for the mutation, and Hom_Ref for the SNP are called ‘ADO’
# - Cells that are Het or Hom_Alt for the SNP are called ‘WT’ for the mutation (because the allele containing the CH mutation was detected).
# When the SNP Alt allele and the mutation are out-of-phase (as in this case):
# - Cells that appear WT for the mutation, and Hom_Alt for the SNP are called 'ADO'
# - Cells that are Het or Hom_Ref for the SNP are called ‘WT’ for the mutation (because the allele containing the CH mutation was detected).

# Call genotypes in gDNA amplicons
# First we assign a SNP genotype to determine whether there was biallelic detection or ADO of either allele
genotyping_gDNA_filtered <- genotyping_gDNA_filtered %>% 
  mutate(SNP.genotype = case_when(
    SNP.VAF >= gDNA_control_VAF$SNP_HomRef_threshold & SNP.VAF <= gDNA_control_VAF$SNP_HomAlt_threshold ~ "Het",
    SNP.VAF < gDNA_control_VAF$SNP_HomRef_threshold ~ "Hom_Ref",
    SNP.VAF > gDNA_control_VAF$SNP_HomAlt_threshold ~ "Hom_Alt")) %>% 
  # Then we assign genotypes for the mutation, while accounting for ADO
  # If SNP is Hom Ref, allelic dropout of allele carrying mutation has occurred so these cells are Undetected
  mutate(genotype = case_when(
    # If mutation VAF is above MUT threshold and has ≥ 10 MUT reads, cell is MUT
    VAF >= gDNA_control_VAF$MUT_threshold & !!mut >= gDNA_min_mut_reads ~ "Mut", 
    # If a cell has allelic dropout of the mutant allele (Hom_Alt for the SNP), we do not know the genotype, so it is called ADO
    Sample == sample & SNP.genotype == "Hom_Alt" & (VAF < gDNA_control_VAF$MUT_threshold | !!mut < gDNA_min_mut_reads) ~ "ADO",
    # All other cells are WT
    TRUE ~ "WT"))

# Call genotypes in mRNA amplicons
genotyping_mRNA_filtered <- genotyping_mRNA_filtered %>% 
  mutate(genotype = case_when(
    # If VAF is below WT threshold, cell is WT
    VAF < mRNA_control_VAF$WT_threshold ~ "WT",  
    # If VAF is above MUT threshold and has ≥ 10 MUT reads, cell is MUT
    VAF >= mRNA_control_VAF$MUT_threshold & !!mut >= mRNA_min_mut_reads ~ "Mut", 
    # If VAF is above MUT threshold but has < 10 MUT reads, cell is Undetermined
    VAF >= mRNA_control_VAF$MUT_threshold & !!mut < mRNA_min_mut_reads  ~ "Undetermined", 
    # If VAF is between MUT and WT thresholds, it is Undetermined
    VAF < mRNA_control_VAF$MUT_threshold & VAF >= mRNA_control_VAF$WT_threshold ~ "Undetermined")) 



# Plot thresholds on a log scale of VAF and mutant read count to check the distribution of test and control cells
ggplot(genotyping_gDNA_filtered, aes(x = !!mut+0.1, y = VAF+0.0001, colour = Sample)) +
  geom_rect(aes(xmin = gDNA_min_mut_reads, xmax =Inf,  ymin = gDNA_control_VAF$MUT_threshold, ymax = Inf),
            fill = "#f0f0f0", colour = NA)+geom_point(position = "jitter")+
  geom_vline(aes(xintercept = gDNA_min_mut_reads),  linetype = "longdash")+
  geom_hline(aes(yintercept = gDNA_control_VAF$MUT_threshold),  linetype = "longdash") +
  geom_hline(aes(yintercept = gDNA_control_VAF$WT_threshold),  linetype = "longdash") +
  labs(title = paste(mutation_name, "gDNA threshold"), y = "VAF", x = "Mutant reads")+
  scale_x_log10()+scale_y_log10()
ggplot(genotyping_gDNA_filtered %>% filter(Sample !=control), 
       aes(x = !!mut+0.1, y = VAF+0.0001, colour = genotype)) +
  geom_rect(aes(xmin = gDNA_min_mut_reads, xmax =Inf,  ymin = gDNA_control_VAF$MUT_threshold, ymax = Inf),
            fill = "#f0f0f0", colour = NA)+geom_point(position = "jitter")+
  geom_vline(aes(xintercept = gDNA_min_mut_reads),  linetype = "longdash")+
  geom_hline(aes(yintercept = gDNA_control_VAF$MUT_threshold),  linetype = "longdash") +
  geom_hline(aes(yintercept = gDNA_control_VAF$WT_threshold),  linetype = "longdash") +
  labs(title = paste(mutation_name, "gDNA threshold"), y = "VAF", x = "Mutant reads")+
  scale_x_log10()+scale_y_log10()

# Plot the SNP genotypes against the SNP VAF and mutation VAF. We invert the SNP VAF for easier visualisation
ggplot(genotyping_gDNA_filtered %>% filter(Sample !=control), 
       aes(x = 1-SNP.VAF+0.0001, y = VAF+0.0001, colour = SNP.genotype)) +
  geom_point(alpha = 0.8)+
  geom_vline(aes(xintercept = gDNA_control_VAF$MUT_threshold),  linetype = "longdash")+
  geom_vline(aes(xintercept = gDNA_control_VAF$WT_threshold),  linetype = "longdash")+
  geom_hline(aes(yintercept = gDNA_control_VAF$SNP_HomRef_threshold),  linetype = "longdash") +
  scale_color_brewer(palette = "Set1")+
  labs(title = paste(mutation_name, "gDNA vs ", SNP_name), x = "SNP VAF", y = "Mutation VAF")+
  scale_x_log10()+scale_y_log10()

# Plot the mutation genotypes against the SNP VAF and mutation VAF
ggplot(genotyping_gDNA_filtered %>% filter(Sample !=control), 
       aes(x = 1-SNP.VAF+0.0001, y = VAF+0.0001, colour = genotype)) +
  geom_point(alpha = 0.8)+
  geom_vline(aes(xintercept = gDNA_control_VAF$MUT_threshold),  linetype = "longdash")+
  geom_vline(aes(xintercept = gDNA_control_VAF$WT_threshold),  linetype = "longdash")+
  geom_hline(aes(yintercept = gDNA_control_VAF$SNP_HomRef_threshold),  linetype = "longdash") +
  labs(title = paste(mutation_name, "gDNA vs ", SNP_name), x = "SNP VAF", y = "Mutation VAF")+
  scale_x_log10()+scale_y_log10()

# Plot thresholds on a log scale of VAF and mutant read count to check the distribution of test and control cells
ggplot(genotyping_mRNA_filtered, aes(x = !!mut+0.1, y = VAF+0.0001, colour = Sample)) +
  geom_rect(aes(xmin = mRNA_min_mut_reads, xmax =Inf,  ymin = mRNA_control_VAF$MUT_threshold, ymax = Inf),
            fill = "#f0f0f0", colour = NA)+geom_point(position = "jitter")+
  geom_vline(aes(xintercept = mRNA_min_mut_reads),  linetype = "longdash")+
  geom_hline(aes(yintercept = mRNA_control_VAF$MUT_threshold),  linetype = "longdash") +
  geom_hline(aes(yintercept = mRNA_control_VAF$WT_threshold),  linetype = "longdash") +
  labs(title = paste(mutation_name, "mRNA threshold"), y = "VAF", x = "Mutant reads")+
  scale_x_log10()+scale_y_log10()
ggplot(genotyping_mRNA_filtered %>% filter(Sample !=control), 
       aes(x = !!mut+0.1, y = VAF+0.0001, colour = genotype)) +
  geom_rect(aes(xmin = mRNA_min_mut_reads, xmax =Inf,  ymin = mRNA_control_VAF$MUT_threshold, ymax = Inf),
            fill = "#f0f0f0", colour = NA)+geom_point(position = "jitter")+
  geom_vline(aes(xintercept = mRNA_min_mut_reads),  linetype = "longdash")+
  geom_hline(aes(yintercept = mRNA_control_VAF$MUT_threshold),  linetype = "longdash") +
  geom_hline(aes(yintercept = mRNA_control_VAF$WT_threshold),  linetype = "longdash") +
  labs(title = paste(mutation_name, "mRNA threshold"), y = "VAF", x = "Mutant reads")+
  scale_x_log10()+scale_y_log10()

# Print the number of mutant and WT cells 
genotyping_gDNA_filtered %>% 
  group_by(Sample, genotype) %>% 
  count() %>% 
  spread(genotype, n)

genotyping_mRNA_filtered %>% 
  group_by(Sample, genotype) %>% 
  count() %>% 
  spread(genotype, n)



##### CONSENSUS between gDNA and mRNA amplicons ######

# Join the mutation calling to the full table
genotyping_gDNA <- genotyping_gDNA_filtered %>% 
  select(Cell, genotype, SNP.genotype) %>% 
  right_join(genotyping_gDNA)

genotyping_mRNA <- genotyping_mRNA_filtered %>% 
  select(Cell, genotype) %>% 
  right_join(genotyping_mRNA)

# Fill in the undetected cells as Undetected
genotyping_gDNA <- genotyping_gDNA %>% 
  mutate(genotype = replace_na(genotype,"Undetected"))
genotyping_mRNA <- genotyping_mRNA %>% 
  mutate(genotype = replace_na(genotype,"Undetected"))

# Create a summary table
genotyping_gDNA_summary <- genotyping_gDNA %>% 
  select(Cell, Sample, Sample_type, Plate, mutation, coverage, !!ref, !!mut, VAF, genotype, SNP, SNP.coverage, SNP.ref, SNP.alt, SNP.VAF, SNP.genotype) %>% 
  rename(Ref = !!ref, 
         Mut = !!mut)
genotyping_mRNA_summary <- genotyping_mRNA %>% 
  select(Cell, Sample, Plate, mutation, coverage, !!ref, !!mut, VAF, genotype) %>% 
  rename(Ref = !!ref, 
         Mut = !!mut)

# Join the gDNA and mRNA tables together
genotyping_joined <- genotyping_gDNA_summary %>% 
  full_join(genotyping_mRNA_summary, by = c("Cell", "Sample", "Plate", "mutation"), suffix = c(".gDNA", ".mRNA"))


# Call a consensus genotype based on gDNA and mRNA
# 1. If the mutation was identified in either the gDNA or cDNA amplicon, the cell is called ‘mutant’, unless gDNA is WT (with proven no ADO)
# 2. If both amplicons were WT, the cell is called ‘WT’.
# 3. If the gDNA amplicon was WT but the cDNA amplicon was undetected or undetermined, the cell is called ‘WT’.
# 4. If the cDNA amplicon was WT but the gDNA amplicon was undetected or undetermined/ADO, the cell is called ‘undetermined’,
# due to the high allelic dropout rate of cDNA amplicons.


genotyping_joined <- genotyping_joined %>% 
  mutate(genotype = case_when(
    # If gDNA is MUT cell is MUT
    genotype.gDNA == "Mut" & genotype.mRNA %in% c("Undetected","Undetermined", "WT", "Mut") ~ "Mut",
    # If mRNA is MUT cell is MUT unless gDNA is WT
    genotype.gDNA %in% c("Undetected","Undetermined", "ADO", "Mut") & genotype.mRNA == "Mut" ~ "Mut",
    genotype.gDNA == "WT" & genotype.mRNA == "Mut" ~ "Undetermined",
    # If both are WT, cell is WT
    genotype.gDNA == "WT" & genotype.mRNA == "WT" ~ "WT",
    # IF gDNA is WT and mRNA anything but Mut, cell is WT
    genotype.gDNA == "WT" & genotype.mRNA %in% c("WT", "Undetected", "Undetermined") ~ "WT",
    # If gDNA is undetected or undetermined, we call cell undetermined
    genotype.gDNA == "Undetected" & genotype.mRNA %in% c("Undetected","Undetermined", "WT") ~ "Undetected",
    genotype.gDNA == "ADO" & genotype.mRNA %in% c("WT", "Undetected", "Undetermined") ~ "ADO",
    TRUE ~ "NA"))


## Write genotyping calls to a summary file ####

genotyping_joined %>% write_tsv(paste("data/genotyping/sample_", sample, "/genotyping_summary_", sample, "_", mutation_name, ".tsv", sep = "")) 

