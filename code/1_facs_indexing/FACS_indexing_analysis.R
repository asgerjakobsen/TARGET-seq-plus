# Analysis of TARGET-seq+ FACS indexing data

## Load the relevant libraries ####
library(tidyverse)

## Import data ####

# Paths to indexing data files
Sort5_index_data_path <- "data/FACS_indexing"


# Index data files
files_Sort5 <- list.files(path = Sort5_index_data_path, pattern = ".*.csv",
                       full.names = T)
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
  file <- read_csv(all_files[8]) %>% 
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
    Sort == "S3" & Lin_PECy5_7AAD <3750 ~ "Lin-", 
    Sort == "S4" & Lin_PECy5_7AAD <2950 ~ "Lin-", 
    Sort == "S5" & Lin_PECy5_7AAD <2900 ~ "Lin-", 
    Sort == "S6" & Sample %in% c("NOC153") & Lin_PECy5_7AAD <2600 ~ "Lin-", 
    Sort == "S6" & Sample %in% c("NOC062", "NOC071", "NOC115", "NOC117") & Lin_PECy5_7AAD <2300 ~ "Lin-", 
    Sample %in% c("NOC156") & Lin_PECy5_7AAD <2300 ~ "Lin-", 
    TRUE ~ "Lin+")) %>% 
  mutate(CD34 = case_when(
    Sort == "S3" & CD34_APC <1300 ~ "CD34-", 
    Sort == "S4" & CD34_APC <1000 ~ "CD34-", 
    Sort == "S5" & CD34_APC <1000 ~ "CD34-", 
    Sort == "S6" & CD34_APC <1000 ~ "CD34-", 
    Sample %in% c("NOC156") & CD34_APC <900 ~ "CD34-", 
    TRUE ~ "CD34+")) %>% 
  mutate(CD38 = case_when(
    Sort == "S3" & CD38_BV421 <1600 ~ "CD38-", 
    Sort == "S4" & CD38_BV421 <1000 ~ "CD38-", 
    Sort == "S5" & CD38_BV421 <900 ~ "CD38-", 
    Sort == "S6" & CD38_BV421 <850 ~ "CD38-", 
    Sample %in% c("NOC156") & CD38_BV421 <900 ~ "CD38-", 
    TRUE ~ "CD38+")) %>% 
  mutate(CD10 = case_when(
    Sort == "S3" & CD38 == "CD38-" & CD10_BV605 <5500 ~ "CD10-", 
    Sort == "S3" & CD38 == "CD38+" & CD10_BV605 <13000 ~ "CD10-", 
    Sort == "S4" & CD38 == "CD38-" & CD10_BV605 <4600 ~ "CD10-", 
    Sort == "S4" & CD38 == "CD38+" & CD10_BV605 <5500 ~ "CD10-", 
    Sort == "S5" & CD38 == "CD38-" & CD10_BV605 <5000 ~ "CD10-", 
    Sort == "S5" & CD38 == "CD38+" & CD10_BV605 <5200 ~ "CD10-", 
    Sort == "S6" & CD38 == "CD38-" & CD10_BV605 <4400 ~ "CD10-", 
    Sort == "S6" & CD38 == "CD38+" & CD10_BV605 <4900 ~ "CD10-", 
    Sample == "NOC156" & CD38 == "CD38-" & CD10_BV605 <3000 ~ "CD10-", 
    Sample == "NOC156" & CD38 == "CD38+" & CD10_BV605 <3500 ~ "CD10-", 
    TRUE ~ "CD10+")) %>% 
  mutate(CD45RA = case_when(
    Sort == "S3" & CD45RA_BB515 <1500 ~ "CD45RA-", 
    Sort == "S4" & CD45RA_BB515 <980 ~ "CD45RA-", 
    Sort == "S5" & Sample %in% c("NOC153", "NOC002") & CD45RA_BB515 <1000 ~ "CD45RA-", 
    Sort == "S5" & Sample %in% c("NOC131", "NOC171", "NOC072") & CD45RA_BB515 <1000 ~ "CD45RA-", 
    Sort == "S6" & Sample %in% c("NOC062") & CD45RA_BB515 <1000 ~ "CD45RA-", 
    Sort == "S6" & Sample %in% c("NOC153", "NOC071", "NOC115", "NOC117") & CD45RA_BB515 <1400 ~ "CD45RA-", 
    Sample %in% c("NOC156") & CD45RA_BB515 <1100 ~ "CD45RA-", 
    TRUE ~ "CD45RA+")) %>% 
  mutate(CD90 = case_when(
    Sort == "S3" & CD90_PECy7  <1500 ~ "CD90-", 
    Sort == "S4" & CD90_PECy7  <700 ~ "CD90-", 
    Sort == "S5" & CD90_PECy7  <720 ~ "CD90-", 
    Sort == "S6" & Sample %in% c("NOC062") & CD90_PECy7  <900 ~ "CD90-",
    Sort == "S6" & Sample %in% c("NOC071", "NOC153", "NOC115", "NOC117") & CD90_PECy7  <1100 ~ "CD90-",
    Sample %in% c("NOC156") & CD90_PECy7  <700 ~ "CD90-",
    TRUE ~ "CD90+")) %>% 
  mutate(CD123 = case_when(
    Sort == "S3" & Sample == "NOC108" & CD123_PE <335 ~ "CD123-",
    Sort == "S3" & Sample == "NOC108" & between(CD123_PE, 335, 11000) ~ "CD123mid",
    Sort == "S3" & Sample == "NOC132" & CD123_PE <650 ~ "CD123-",
    Sort == "S3" & Sample == "NOC132" & between(CD123_PE, 650, 15000) ~ "CD123mid",
    Sort == "S3" & Sample %in% c("NOC153","NOC031") & CD123_PE <450 ~ "CD123-",
    Sort == "S3" & Sample %in% c("NOC153","NOC031") & between(CD123_PE, 450, 11000) ~ "CD123mid",
    Sort == "S4" & Sample == "NOC137" & CD123_PE <350 ~ "CD123-",
    Sort == "S4" & Sample == "NOC137" & between(CD123_PE, 350, 10000) ~ "CD123mid",
    Sort == "S4" & Sample %in% c("NOC153","NOC062","NOC062","NOC131","NOC117") & CD123_PE <450 ~ "CD123-",
    Sort == "S4" & Sample %in% c("NOC153","NOC062","NOC131","NOC117") & between(CD123_PE, 450, 11000) ~ "CD123mid",
    Sort == "S5" & Sample %in% c("NOC002", "NOC072") & CD123_PE < 400 ~ "CD123-",
    Sort == "S5" & Sample == "NOC171" & CD123_PE <150 ~ "CD123-",
    Sort == "S5" & Sample %in% c("NOC153","NOC131") & CD123_PE <300 ~ "CD123-",
    Sort == "S5" & Sample %in% c("NOC002", "NOC072") & between(CD123_PE, 400, 12000) ~ "CD123mid",
    Sort == "S5" & Sample == "NOC171" & between(CD123_PE, 150, 2500) ~ "CD123mid",
    Sort == "S5" & Sample %in% c("NOC153","NOC131") & between(CD123_PE, 300, 12000) ~ "CD123mid",
    Sort == "S6" & Sample %in% c("NOC153",  "NOC115") & CD123_PE <400 ~ "CD123-",
    Sort == "S6" & Sample %in% c("NOC117", "NOC071","NOC062") & CD123_PE <450 ~ "CD123-",
    Sample %in% c("NOC156") & between(CD123_PE, 300, 10000) ~ "CD123mid",
    Sample %in% c("NOC156") & CD123_PE <300 ~ "CD123-",
    Sort == "S6" & between(CD123_PE, 400, 9500) ~ "CD123mid",
    TRUE ~ "CD123hi")) %>% 
  mutate(CD49f = case_when(
    Sort == "S3" & CD49f_PEDaz <1400 ~ "CD49f-", #1900
    Sort == "S4" & CD49f_PEDaz <1250 ~ "CD49f-", #1800
    Sort == "S5" & CD49f_PEDaz < 1000 ~ "CD49f-", #1600
    Sort == "S6" & CD49f_PEDaz <1000 ~ "CD49f-", #1700
    Sample == "NOC156"  ~ "NA",
    TRUE ~ "CD49f+")) %>% 
  mutate(CD117 = case_when(
    Sort == "S3" & CD117_BV785 <3400 ~ "CD117-",
    Sort == "S4" & CD117_BV785 <2000 ~ "CD117-",
    Sort == "S5" & CD117_BV785 < 1800 ~ "CD117-",
    Sort == "S6" & CD117_BV785 <1500 ~ "CD117-",
    Sample == "NOC156" & CD117_BV785 <3000 ~ "CD117-",
    TRUE ~ "CD117+"))



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


