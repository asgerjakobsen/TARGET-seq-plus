# Analysis of TARGET-seq+ FACS indexing data


# This step gathers cell surface immunophenotyping data from FACS index sorting files into a single table. 
# From these data, we know which cell type was sorted into each well and the sample donor ID. 
# This is important for defining which cells came from control and test samples when performing downstream analysis.

# In this example, we import FACS index files exported from the Sony MA900 cell sorter.

## Load the relevant libraries ####
library(tidyverse)

## Import data ####

# Paths to indexing data files
Sort5_index_data_path <- "data/FACS_indexing"

# Index data files
files_Sort5 <- list.files(path = Sort5_index_data_path, pattern = ".*.csv",
                       full.names = T)

# Can append files from other FACS sorts
# all_files <- append(files_Sort5)

all_files <- files_Sort5


# Functions to extract plate and sample Sample name
index_name <- function(column) {
  row <- str_extract(column, "^\\D")
  col <- str_extract(column, "\\d+$")
  index <- if_else(str_length(col) == 1, paste("0", col, row, sep = ""), paste(col, row, sep = ""))
  return(index)
}

extract_sort_number <- function(file_path) {
  sort <- str_extract(file_path, "Sort \\d+?") %>% 
    str_remove("ort ")
  return(sort)
}

extract_plate_number <- function(file_path) {
  plate_no <- str_extract(file_path, "Plate - \\d+") %>%
    str_remove("Plate - ")
  plate_no <- sprintf("%02d", as.integer(plate_no))
  plate <- paste("PL", plate_no, sep = "")
  return(plate)
}

extract_sample_name <- function(file_path) {
  sample <- str_extract(file_path, "Bone Marrow_NOC\\d+") %>% 
    str_remove("Bone Marrow_")
  return(sample)
}


###  Read in files ######

# Define a function that reads file and reformats index columns to obtain sample ID
read_index_file <- function(file_path) {
  file <- read_csv(file_path) %>% 
    mutate(Index = index_name(Index)) %>% 
    mutate(Sort = extract_sort_number(file_path)) %>% 
    mutate(Plate = paste(extract_sort_number(file_path),
                         extract_plate_number(file_path), sep = "")) %>% 
    mutate(Cell = paste(Plate,
                        Index,
                        extract_sample_name(file_path),
                        sep = "-")) %>% 
    mutate(Sample = extract_sample_name(file_path))
  return(file)
}


# Read in each file
index_data_all_plates <- all_files %>% 
  map_df(~read_index_file(.))

head(index_data_all_plates$Cell)

# Rename the cell surface markers
index_data_all_plates <- index_data_all_plates %>% 
  dplyr::rename(CD45RA_BB515 = `CD45RA: BB515-A-Compensated`,
         CD123_PE = `CD123: PE-A-Compensated`,
         CD49f_PEDaz = `CD49f: PE-Dazzle-A-Compensated`,
         Lin_PECy5_7AAD = `Lin Live/Dead: PE-Cy5-A-Compensated`,
         CD90_PECy7 = `CD90: PE-Cy7-A-Compensated`,
         CD38_BV421 = `CD38: Brilliant Violet 421-A-Compensated`,
         CD10_BV605 = `CD10: Brilliant Violet 605-A-Compensated`,
         CD34_APC = `CD34: APC-A-Compensated`,
         CD117_BV785 = `CD117: Brilliant Violet 785-A-Compensated`)

## Define the genotyping control samples ####

# In this experiment, Sample NOC153 was used as the WT genotyping control
# All others were from clonal haematopoiesis (CH) samples
control_samples <- c("NOC153")

index_data_all_plates <- index_data_all_plates %>% 
  mutate(Sample_type = ifelse(Sample %in% control_samples, "Control", "CH"))


## Define immunophenotypic populations based on FACS gates ####

# First we set gates for defining positive and negative populations
# Due to variability in instrument settings, we set these separately for each FACS sort 

index_data_all_plates <- index_data_all_plates %>% 
  mutate(Lin = case_when(
    Sort == "S5" & Lin_PECy5_7AAD <2900 ~ "Lin-", 
    TRUE ~ "Lin+")) %>% 
  mutate(CD34 = case_when(
    Sort == "S5" & CD34_APC <1000 ~ "CD34-", 
    TRUE ~ "CD34+")) %>% 
  mutate(CD38 = case_when(
    Sort == "S5" & CD38_BV421 <900 ~ "CD38-", 
    TRUE ~ "CD38+")) %>% 
  mutate(CD10 = case_when(
    Sort == "S5" & CD38 == "CD38-" & CD10_BV605 <5000 ~ "CD10-", 
    Sort == "S5" & CD38 == "CD38+" & CD10_BV605 <5200 ~ "CD10-", 
    TRUE ~ "CD10+")) %>% 
  mutate(CD45RA = case_when(
    Sort == "S5" & CD45RA_BB515 <1000 ~ "CD45RA-", 
    TRUE ~ "CD45RA+")) %>% 
  mutate(CD90 = case_when(
    Sort == "S5" & CD90_PECy7  <720 ~ "CD90-", 
    TRUE ~ "CD90+")) %>% 
  mutate(CD123 = case_when(
    Sort == "S5" & CD123_PE <300 ~ "CD123-",
    Sort == "S5" & between(CD123_PE, 300, 12000) ~ "CD123mid",
    TRUE ~ "CD123hi")) %>% 
  mutate(CD49f = case_when(
    Sort == "S5" & CD49f_PEDaz < 1000 ~ "CD49f-",
    TRUE ~ "CD49f+")) %>% 
  mutate(CD117 = case_when(
    Sort == "S5" & CD117_BV785 < 1800 ~ "CD117-",
    TRUE ~ "CD117+"))


# Then define immunophenotypic populations based on the protein expression pattern

index_data_all_plates <- index_data_all_plates %>% 
  mutate(Lin_CD34_CD38 = case_when(
    Lin == "Lin+" ~ "Lin+",
    Lin == "Lin-" & CD34 == "CD34-" ~ "Lin-CD34-",
    Lin == "Lin-" & CD34 == "CD34+" & CD38 == "CD38-" ~ "CD34+CD38-",
    Lin == "Lin-" & CD34 == "CD34+" & CD38 == "CD38+" ~ "CD34+CD38+",
    TRUE ~ "NA")) %>% 
  mutate(Immunophenotype = case_when(
    Lin_CD34_CD38 == "Lin+" ~ "Lin+",
    Lin_CD34_CD38 == "Lin-CD34-" ~ "Lin-CD34-",
    Lin_CD34_CD38 == "CD34+CD38-" & CD10 == "CD10+" ~ "MLP",
    Lin_CD34_CD38 == "CD34+CD38-" & CD10 == "CD10-" & CD45RA == "CD45RA+" ~ "LMPP",
    Lin_CD34_CD38 == "CD34+CD38-" & CD10 == "CD10-" & CD45RA == "CD45RA-" & CD90 == "CD90+" ~ "HSC",
    Lin_CD34_CD38 == "CD34+CD38-" & CD10 == "CD10-" & CD45RA == "CD45RA-" & CD90 == "CD90-" ~ "MPP",
    Lin_CD34_CD38 == "CD34+CD38+" & CD10 == "CD10+" ~ "B-NK",
    Lin_CD34_CD38 == "CD34+CD38+" & CD10 == "CD10-" & CD45RA == "CD45RA+" & CD123 == "CD123mid" ~ "GMP",
    Lin_CD34_CD38 == "CD34+CD38+" & CD10 == "CD10-" & CD45RA == "CD45RA-" & CD123 == "CD123mid" ~ "CMP",
    Lin_CD34_CD38 == "CD34+CD38+" & CD10 == "CD10-" & CD45RA == "CD45RA-" & CD123 == "CD123hi" ~ "CMP",
    Lin_CD34_CD38 == "CD34+CD38+" & CD10 == "CD10-" & CD45RA == "CD45RA-" & CD123 == "CD123-" ~ "MEP",
    Lin_CD34_CD38 == "CD34+CD38+" & CD10 == "CD10-" & CD45RA == "CD45RA+" & CD123 == "CD123hi" ~ "pDC",
    Lin_CD34_CD38 == "CD34+CD38+" & CD10 == "CD10-" & CD45RA == "CD45RA+" & CD123 == "CD123-" ~ "CD38+CD10-CD123-"))


# Check for wells containing more than one cell
duplicates <- index_data_all_plates %>% group_by(Index, Plate) %>% dplyr::count() %>% filter(n>1) %>% left_join(index_data_all_plates)

# Add a column so these can be excluded from analysis
duplicates <- duplicates %>% 
  select(Cell, Index, Plate) %>% mutate(Single_cell = "FALSE")

index_data_all_plates <- index_data_all_plates %>% 
  left_join(duplicates) %>% 
  mutate(Single_cell = if_else(is.na(Single_cell), "TRUE", Single_cell))


### Write the FACS data to a file ####

index_data_all_plates %>% 
  select(Cell, Sample, Sample_type, Sort, Plate, Index, 1:16, 23:34) %>% 
  write_tsv("data/FACS_indexing/FACS_index_data.tsv")


